#!/usr/bin/env python
"""Plain-TI pooled SIGN analysis for the direct-FE release-direction coverage wave (run_sign.sh output).

For each member with both legs in runs/fe/out:  ddG_switch = dG_bound - dG_free  (per leg = dG(HIE->HIP)).
  ddG > 0  => His+ destabilizes the BOUND complex more than the free binder => RELEASES at low pH (GOOD).
  ddG < 0  => anti-switch: His+ *stabilizes* binding at low pH = a pH-selectivity (#1-axis) liability.

Estimators:
  (A) PAIRED per-seed  ddG_s = dG_bound_s - dG_free_s over seeds present in BOTH legs; mean/SEM; z=mean/SEM.
      The robust estimator — the double difference cancels correlated per-seed noise.
  (B) INDEPENDENT pool  mean_bound - mean_free, err = quad-sum of the two SEMs.

Confirmation bar (calibrated on the one known anti-switch, des_4098_1_b14: paired ddG -3.9, z=-6.3,
all 6 seeds negative): a CONFIRMED anti-switch needs paired ddG reproducibly negative =
  (all or all-but-one seeds negative) AND |z| > 3   (Bonferroni-safe for a ~20-member portfolio).

Usage:  python scripts/fe/analyze_sign.py [--csv submission.csv] [--out runs/fe/out/sign_summary.json]
With --csv, restricts to designs listed in that CSV (matched by id prefix); otherwise analyzes every
member found in runs/fe/out. Output JSON is written under runs/ (gitignored).
"""
import argparse, json, glob, os, re, math, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "runs", "fe", "out")

def discover_members(outdir):
    bound, free = set(), set()
    for f in glob.glob(os.path.join(outdir, "*_bound_s*.json")):
        b = os.path.basename(f)
        bound.add(re.sub(r"_bound_s\d+\.json$", "", b))
    for f in glob.glob(os.path.join(outdir, "*_free_s*.json")):
        b = os.path.basename(f)
        free.add(re.sub(r"_free_s\d+\.json$", "", b))
    return sorted(bound & free)

def load_leg(outdir, member, leg):
    d = {}
    for f in glob.glob(os.path.join(outdir, f"{member}_{leg}_s*.json")):
        m = re.search(r"_s(\d+)\.json$", f)
        if m:
            d[int(m.group(1))] = json.load(open(f))["dG_kcal"]
    return d

def sem(a):
    a = np.asarray(a, float)
    return float(a.std(ddof=1) / math.sqrt(len(a))) if len(a) > 1 else float("nan")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="restrict to designs in this submission CSV (match by id prefix)")
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--out", default=os.path.join(OUTDIR, "sign_summary.json"))
    args = ap.parse_args()

    members = discover_members(args.outdir)
    if args.csv:
        ids = []
        with open(args.csv) as fh:
            for row in csv.DictReader(fh):
                ids.append(next((row[k] for k in ("id", "name", "design_id") if k in row), "").strip())
        idset = {i for i in ids if i}
        members = [m for m in members if m in idset]
        missing = sorted(idset - set(members))
        if missing:
            print(f"# NOTE: {len(missing)} CSV design(s) have no bound+free coverage in {args.outdir}:")
            for m in missing:
                print(f"#   {m}")

    print(f"{'member':<42} {'ddG_paired':>12} {'z':>6} {'signs p(+/-)':>12} {'ddG_pool':>9}  verdict")
    print("-" * 112)
    rows = []
    for mem in members:
        b = load_leg(args.outdir, mem, "bound"); f = load_leg(args.outdir, mem, "free")
        seeds = sorted(set(b) & set(f))
        if not seeds:
            continue
        bo = np.array([b[s] for s in seeds]); fr = np.array([f[s] for s in seeds])
        paired = bo - fr
        ddg_p = float(paired.mean()); e_p = sem(paired)
        z_p = ddg_p / e_p if e_p == e_p and e_p > 0 else float("nan")
        npos = int((paired > 0).sum()); nneg = int((paired < 0).sum())
        ddg_i = float(bo.mean() - fr.mean())
        e_i = math.sqrt(sem(bo) ** 2 + sem(fr) ** 2) if len(bo) > 1 and len(fr) > 1 else float("nan")
        # A confirmed verdict needs enough seeds for a trustworthy z (the campaign standard is 6; a 2-seed
        # SEM is not reliable). Below MINSEED, report the lean but never "confirmed".
        MINSEED = 4
        enough = len(seeds) >= MINSEED
        anti = enough and (nneg >= len(seeds) - 1) and (z_p == z_p) and (z_p < -3.0)
        rel = enough and (npos >= len(seeds) - 1) and (z_p == z_p) and (z_p > 3.0)
        if not enough:
            verdict = "INSUFFICIENT_SEEDS"
        elif anti:
            verdict = "CONFIRMED_ANTI"
        elif rel:
            verdict = "CONFIRMED_REL"
        elif ddg_p < 0:
            verdict = "NEG_LEAN_UNCONFIRMED"
        else:
            verdict = "POS_NEUTRAL_CLEAN"
        tag = {"CONFIRMED_ANTI": "*** CONFIRMED ANTI-SWITCH ***", "CONFIRMED_REL": "confirmed releaser (+)",
               "NEG_LEAN_UNCONFIRMED": "neg lean, NOT confirmed", "POS_NEUTRAL_CLEAN": "pos/neutral, clean",
               "INSUFFICIENT_SEEDS": f"only {len(seeds)} seeds (<{MINSEED}) — not ranked"}[verdict]
        print(f"{mem:<42} {ddg_p:>+8.2f}+-{e_p:>4.1f} {z_p:>+6.1f} {npos:>5}/{nneg:<5} {ddg_i:>+8.2f}  {tag}")
        rows.append(dict(member=mem, n_seeds=len(seeds), seeds=seeds, ddg_paired=ddg_p, sem_paired=e_p,
                         z_paired=z_p, n_pos=npos, n_neg=nneg, ddg_pool=ddg_i, sem_pool=e_i,
                         bound=bo.tolist(), free=fr.tolist(), verdict=verdict))
    print("-" * 112)
    n_anti = sum(r["verdict"] == "CONFIRMED_ANTI" for r in rows)
    n_rel = sum(r["verdict"] == "CONFIRMED_REL" for r in rows)
    print(f"\n{len(rows)} members | CONFIRMED anti-switch: {n_anti} | CONFIRMED releaser: {n_rel}")
    print("calibration (the one known anti-switch, des_4098_1_b14): paired ddG -3.9, z=-6.3, all 6 seeds neg.")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(rows, open(args.out, "w"), indent=2)
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
