#!/usr/bin/env python3
"""Pick a diverse batch of engineered designs to fold-validate with AF2-ig.

Reads results/obj1_switch.csv, gates on a working switch + a Rosetta interface, ranks by the
conserved-sensor ddG_pH (His433-driven preferred), diversifies across source backbones, and
copies the chosen best-variant PDBs into an input dir for scripts/run_af2ig.py.
"""
from __future__ import annotations
import argparse, csv, os, re, shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def backbone_key(design):
    # e.g. egfr_dom3_H_l_seed0_egfr_dom3_H_l_088 -> spec+seed+id; group by spec to diversify
    m = re.match(r"(egfr_dom3_[A-Za-z]+_[sl])_seed(\d+)", design)
    return m.group(1) if m else design


def fnum(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(REPO / "results" / "obj1_switch.csv"))
    ap.add_argument("--outdir", default=str(REPO / "runs" / "obj1_foldin"))
    ap.add_argument("--n", type=int, default=48)
    ap.add_argument("--min-ddgph", type=float, default=0.4)
    ap.add_argument("--max-iddg", type=float, default=-15.0)
    ap.add_argument("--per-spec-cap", type=int, default=10)
    a = ap.parse_args()

    rows = [r for r in csv.DictReader(open(a.csv)) if not r.get("error")]
    cand = [r for r in rows
            if fnum(r.get("ddg_pH")) >= a.min_ddgph and fnum(r.get("interface_ddg"), 1) <= a.max_iddg
            and r.get("best_pdb") and os.path.exists(r["best_pdb"])]
    # rank: His433-driven first (his433 term), then total sensor ddG_pH, then affinity
    cand.sort(key=lambda r: (fnum(r.get("ddg_pH_his433")), fnum(r.get("ddg_pH")), -fnum(r.get("interface_ddg"))),
              reverse=True)

    out = Path(a.outdir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    per = {}
    picked = []
    for r in cand:
        k = backbone_key(r["design"])
        if per.get(k, 0) >= a.per_spec_cap:
            continue
        per[k] = per.get(k, 0) + 1
        dst = out / (r["design"] + "__best.pdb")
        shutil.copyfile(r["best_pdb"], dst)
        picked.append(r)
        if len(picked) >= a.n:
            break
    print(f"selected {len(picked)} designs -> {out}")
    print("per-spec:", per)
    hd = sum(1 for r in picked if fnum(r.get("ddg_pH_his433")) > 0.3)
    print(f"His433-driven in batch: {hd}/{len(picked)}")
    print("ddg_pH range:", round(min(fnum(r['ddg_pH']) for r in picked), 2),
          "..", round(max(fnum(r['ddg_pH']) for r in picked), 2))


if __name__ == "__main__":
    main()
