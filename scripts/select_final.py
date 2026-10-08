#!/usr/bin/env python3
"""Final Track-3 selection (<=20). For each AF2-validated passer, consider BOTH submittable
versions and keep the best one that BINDS and carries a pH switch:
  - ENGINEERED: binder with the installed His433 acid (seq from runs/eng_foldin/<d>__best.pdb,
                binding = af2ig_eng pae, switch = passers_switch ddg_pH). Stronger switch, but the
                mutation sometimes breaks binding.
  - ORIGINAL  : the native BoltzGen binder (seq from runs/boltz_more_foldin/<d>.pdb, binding =
                af2_passers pae, NATURAL switch = passers_switch orig_ddg_pH). No mutation, no risk.
Pick per design the version with the higher ddG_pH among those that bind (pae < --max-pae).

Ranking = competition priority pH-selective > cross-reactive > affinity, His433-driven preferred.
Cross-reactivity from the designed-complex conserved/divergent contacts (analyzed.csv). Diversify by
sequence; write submission CSV + a rich manifest for the methods write-up.
"""
from __future__ import annotations
import argparse, csv, glob, os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
D3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}


def chainA_seq(pdb):
    seq = []
    try:
        for l in open(pdb):
            if l[:4] == "ATOM" and l[21] == "A" and l[12:16].strip() == "CA":
                seq.append(D3.get(l[17:20].strip(), "X"))
    except OSError:
        return ""
    return "".join(seq)


def load_sc(paths):
    out = {}
    for pat in paths:
        for f in glob.glob(pat):
            hdr = None
            for ln in open(f):
                t = ln.split()
                if ln.startswith("SCORE:") and "description" in ln:
                    hdr = t[1:]
                elif ln.startswith("SCORE:") and hdr:
                    d = dict(zip(hdr, t[1:]))
                    b = d["description"].replace("_af2pred", "").replace("__best", "")
                    out[b] = d
    return out


def ff(d, k, dv=99.0):
    try:
        return float(d[k])
    except (KeyError, ValueError, TypeError):
        return dv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-pae", type=float, default=12.0)
    ap.add_argument("--min-ddgph", type=float, default=0.3)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--out", default=str(REPO / "results" / "submission_track3.csv"))
    ap.add_argument("--manifest", default=str(REPO / "results" / "submission_manifest.csv"))
    a = ap.parse_args()

    ph = {r["design"].replace("__best", ""): r for r in csv.DictReader(open(REPO / "results" / "passers_switch.csv")) if not r.get("error")}
    eng_pae = load_sc([str(REPO / "runs" / "af2ig_eng" / "gpu*" / "out.sc")])
    orig_af2 = {r["design"].replace("__best", ""): r for r in csv.DictReader(open(REPO / "results" / "af2_passers.csv"))}
    analyzed = {}
    for r in csv.DictReader(open(REPO / "results" / "designs_analyzed.csv")):
        stem = f"{r['source_spec']}_{r['source_seed']}_{r['design_id']}"
        analyzed[stem] = r
    cc = {r["design"]: r for r in csv.DictReader(open(REPO / "results" / "contact_conservation.csv"))} \
        if (REPO / "results" / "contact_conservation.csv").exists() else {}

    cand = []
    for d, p in ph.items():
        acons = ff(analyzed.get(d, {}), "n_conserved_contacts", 0)
        adiv = ff(analyzed.get(d, {}), "n_divergent_contacts", 9)
        his = ff(p, "ddg_pH_his433")
        iddg = ff(p, "interface_ddg")
        versions = []
        # engineered
        e_pae = ff(eng_pae.get(d, {}), "pae_interaction")
        e_seq = chainA_seq(REPO / "runs" / "eng_foldin" / f"{d}__best.pdb")
        if e_pae < a.max_pae and e_seq:
            versions.append(("engineered", e_seq, e_pae, ff(p, "ddg_pH"), his))
        # original (natural switch)
        o_pae = ff(orig_af2.get(d, {}), "pae_interaction")
        o_seq = chainA_seq(REPO / "runs" / "boltz_more_foldin" / f"{d}.pdb")
        if o_pae < a.max_pae and o_seq:
            # natural switch has no engineered His433 term; treat his433 term as 0 unless orig raises it
            versions.append(("original", o_seq, o_pae, ff(p, "orig_ddg_pH"), 0.0))
        if not versions:
            continue
        # best viable version by ddG_pH
        versions.sort(key=lambda v: v[3], reverse=True)
        kind, seq, pae, ddgph, histerm = versions[0]
        if ddgph < a.min_ddgph:
            continue
        cand.append({
            "design": d, "version": kind, "sequence": seq, "length": len(seq),
            "ddg_pH": round(ddgph, 2), "his433_term": round(histerm, 2),
            "his433_driven": histerm > 0.3 or (kind == "original" and ff(p, "orig_ddg_pH") > 0.3 and his > 0.3),
            "af2_pae": round(pae, 1), "n_conserved": int(acons), "n_divergent": int(adiv),
            "iddg": round(iddg, 1), "switch_tag": p.get("best_tag", ""),
        })

    # rank: pH (His433-driven, ddg_pH) > cross-react (conserved desc, divergent asc) > affinity
    cand.sort(key=lambda c: (
        0 if c["his433_term"] > 0.3 else 1,
        -c["ddg_pH"],
        -c["n_conserved"], c["n_divergent"],
        c["af2_pae"], c["iddg"]))
    seen = set(); final = []
    for c in cand:
        if c["sequence"] in seen:
            continue
        seen.add(c["sequence"]); final.append(c)
        if len(final) >= a.n:
            break

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["name", "sequence", "molecule_class"])
        for i, c in enumerate(final, 1):
            w.writerow([f"egfr_phsw_{i:02d}", c["sequence"], "protein"])
    with open(a.manifest, "w", newline="") as fh:
        cols = ["rank", "name", "design", "version", "length", "ddg_pH", "his433_term",
                "his433_driven", "af2_pae", "n_conserved", "n_divergent", "iddg", "switch_tag", "sequence"]
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for i, c in enumerate(final, 1):
            w.writerow({"rank": i, "name": f"egfr_phsw_{i:02d}", **{k: c.get(k, "") for k in cols if k not in ("rank", "name")}})

    print(f"viable candidates: {len(cand)} | selected {len(final)} -> {a.out}")
    print(f"  His433-driven: {sum(1 for c in final if c['his433_term']>0.3)}/{len(final)} | "
          f"cross-react(cons>=3,div<=1): {sum(1 for c in final if c['n_conserved']>=3 and c['n_divergent']<=1)}/{len(final)}")
    print(f"{'#':>2} {'design':40s} {'ver':4s} {'ddgpH':>5s} {'his433':>6s} {'pae':>5s} {'cons':>4s} {'div':>3s} {'iddg':>6s} L")
    for i, c in enumerate(final, 1):
        print(f"{i:2d} {c['design'][:40]:40s} {c['version'][:4]:4s} {c['ddg_pH']:5.2f} {c['his433_term']:6.2f} "
              f"{c['af2_pae']:5.1f} {c['n_conserved']:4d} {c['n_divergent']:3d} {c['iddg']:6.1f} {c['length']}")


if __name__ == "__main__":
    main()
