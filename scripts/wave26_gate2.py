#!/usr/bin/env python3
"""WAVE 26 gate 2: on af2ig BINDERS, measure how BURIED each His...Asp460 salt bridge is (bridge_burial) and
rank. The buried ones (bridge_sasa_complex < 10 A^2, burial_fraction >= 0.85, bridge_dist < 3.3) are the
candidates that then go to APBS (strong?) + MD-ECC (robust & non-constitutive?). Finds the binder His nearest a
target ASP (the locked motif, = Asp460) per design."""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
import sys
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
import bridge_burial as BB


def load_binders(af2ig_dir):
    keep = {}
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "out.sc")):
        L = open(f).read().splitlines()
        if not L:
            continue
        h = L[0].split()
        try:
            ip, ipl, ir, idc = (h.index("pae_interaction"), h.index("plddt_binder"),
                                h.index("binder_aligned_rmsd"), h.index("description"))
        except ValueError:
            continue
        for ln in L[1:]:
            p = ln.split()
            if len(p) <= max(ip, ipl, ir, idc):
                continue
            try:
                pae, plddt, rms = float(p[ip]), float(p[ipl]), float(p[ir])
            except ValueError:
                continue
            if pae < 10 and plddt > 80 and rms < 2:
                keep[p[idc]] = round(pae, 2)
    return keep


def find_pdb(af2ig_dir, desc):
    hits = glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", desc + "*.pdb"))
    return hits[0] if hits else None


def binder_his_near_asp(pdb):
    m = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    A = [r for r in m["A"] if r.id[0] == " "]; B = [r for r in m["B"] if r.id[0] == " "]
    asps = [r for r in B if r.resname in ("ASP", "GLU")]
    best = None
    for h in A:
        if not h.resname.startswith("HI"):
            continue
        hn = [h[a].coord for a in ("ND1", "NE2") if a in h]
        for ar in asps:
            for o in ("OD1", "OD2", "OE1", "OE2"):
                if o in ar:
                    d = min(np.linalg.norm(ar[o].coord - x) for x in hn)
                    if best is None or d < best[0]:
                        best = (d, h.id[1], ar.id[1])
    return best  # (dist, his_resid, acid_resid)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2ig", default="runs/wave26_buried/af2ig")
    ap.add_argument("--out", default="runs/wave26_buried/gate2.csv")
    a = ap.parse_args()
    binders = load_binders(a.af2ig)
    print(f"af2ig binders: {len(binders)}")
    rows = []
    for desc, pae in sorted(binders.items(), key=lambda kv: kv[1]):
        pdb = find_pdb(a.af2ig, desc)
        if not pdb:
            continue
        nh = binder_his_near_asp(pdb)
        if nh is None or nh[0] > 4.0:
            continue
        try:
            r = BB.run(pdb, "A", "B", nh[1], nh[2], his_atoms=("NE2",))  # inter-ion burial (asymmetric design)
        except Exception as e:
            print(f"  burial FAIL {desc[:40]}: {repr(e)[:60]}"); continue
        r.update(desc=desc, pae=pae, pdb=os.path.abspath(pdb))
        rows.append(r)
        tag = "BURIED" if (r["bridge_sasa_complex"] < 10 and (r["burial_fraction"] or 0) >= 0.85 and r["bridge_dist"] < 3.3) else ""
        print(f"  {desc[:38]:38s} pae={pae:5.2f} His{r['his']}-Asp{r['acid']} dist={r['bridge_dist']} "
              f"bridge_SASA={r['bridge_sasa_complex']:5.1f} burial={r['burial_fraction']} {tag}")
    buried = [r for r in rows if r["bridge_sasa_complex"] < 10 and (r["burial_fraction"] or 0) >= 0.85 and r["bridge_dist"] < 3.3]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["desc", "pae", "his", "acid", "acid_resname", "bridge_dist",
                                           "bridge_sasa_complex", "bridge_sasa_free", "burial_fraction",
                                           "buried", "surface", "pdb"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in w.fieldnames})
    json.dump(buried, open(a.out.replace(".csv", "_buried.json"), "w"), indent=1)
    print(f"\nBURIED bridges (SASA<10 & burial>=0.85 & dist<3.3): {len(buried)} / {len(rows)}")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
