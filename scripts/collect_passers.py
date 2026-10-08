#!/usr/bin/env python3
"""Collect AF2-ig passing binders (their AF2-predicted complex PDBs) from one or more af2ig runs,
so the pH-switch engine can be run on the VALIDATED binding poses.

Gates: pae_interaction < --max-pae, binder_aligned_rmsd < --max-rmsd, plddt_binder > --min-plddt.
Copies <design>_af2pred.pdb into --outdir and writes a manifest of AF2 metrics.
"""
from __future__ import annotations
import argparse, csv, glob, os, shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=["runs/af2ig", "runs/af2ig2", "runs/af2ig_more"])
    ap.add_argument("--outdir", default=str(REPO / "runs" / "passers"))
    ap.add_argument("--manifest", default=str(REPO / "results" / "af2_passers.csv"))
    ap.add_argument("--max-pae", type=float, default=12.0)
    ap.add_argument("--max-rmsd", type=float, default=2.5)
    ap.add_argument("--min-plddt", type=float, default=75.0)
    a = ap.parse_args()

    # gather scores + locate af2pred pdbs
    scores = {}
    pdbmap = {}
    for run in a.runs:
        for sc in glob.glob(os.path.join(run, "gpu*", "out.sc")) + glob.glob(os.path.join(run, "out.sc")):
            hdr = None
            for ln in open(sc):
                t = ln.split()
                if ln.startswith("SCORE:") and "description" in ln:
                    hdr = t[1:]
                elif ln.startswith("SCORE:") and hdr:
                    d = dict(zip(hdr, t[1:]))
                    scores[d["description"]] = d
        for p in glob.glob(os.path.join(run, "gpu*", "pdbs", "*_af2pred.pdb")):
            pdbmap[Path(p).stem] = p

    def f(d, k, dv=99.0):
        try:
            return float(d[k])
        except (KeyError, ValueError, TypeError):
            return dv

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for tag, d in scores.items():
        pae = f(d, "pae_interaction"); rmsd = f(d, "binder_aligned_rmsd"); plddt = f(d, "plddt_binder", 0)
        if pae < a.max_pae and rmsd < a.max_rmsd and plddt > a.min_plddt:
            src = pdbmap.get(tag)
            design = tag.replace("_af2pred", "")
            if src and os.path.exists(src):
                shutil.copyfile(src, out / f"{design}.pdb")
                rows.append({"design": design, "pae_interaction": round(pae, 2),
                             "binder_rmsd": round(rmsd, 2), "plddt_binder": round(plddt, 1),
                             "af2pred_pdb": str(out / f"{design}.pdb")})
    rows.sort(key=lambda r: r["pae_interaction"])
    Path(a.manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(a.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["design", "pae_interaction", "binder_rmsd", "plddt_binder", "af2pred_pdb"])
        w.writeheader()
        w.writerows(rows)
    print(f"collected {len(rows)} passers (pae<{a.max_pae}, rmsd<{a.max_rmsd}, plddt>{a.min_plddt}) -> {out}")
    strict = [r for r in rows if r["pae_interaction"] < 10 and r["binder_rmsd"] < 2 and r["plddt_binder"] > 80]
    print(f"  of which STRICT (pae<10,rmsd<2,plddt>80): {len(strict)}")
    for r in rows[:25]:
        print(f"  {r['design'][:46]:46s} pae={r['pae_interaction']:5.1f} rmsd={r['binder_rmsd']:4.1f} plddt={r['plddt_binder']:4.0f}")


if __name__ == "__main__":
    main()
