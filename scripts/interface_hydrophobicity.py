#!/usr/bin/env python3
"""Interface hydrophobicity of each submitted binder = fraction of the binder's INTERFACE residues
(any heavy atom within 5 A of the target) that are hydrophobic (AILMFVWY). This is the metric the
GuideFlip paper (Scheres 2026) and binder-design practice flag: direct AF2 optimization buries
hydrophobics at the interface to inflate pae, which predicts poorly for wet-lab binding/solubility.
A polar/charged interface (which a His433 Asp/Glu salt-bridge switch naturally favours) is preferred.
Usage: interface_hydrophobicity.py [manifest.csv]"""
import csv, glob, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
HYD = {"ALA", "ILE", "LEU", "MET", "PHE", "VAL", "TRP", "TYR"}
PASSER_DIRS = ["rfd_passers_pd", "rfd_passers_pd2", "rfd_passers_pd3", "rfd_passers_w1fa", "rfd_passers_w2fa"]
EXTRA = ["obj1_foldin2/{d}__best.pdb", "passers/{d}.pdb", "filtered/{d}.pdb", "cand_final2_canon/{d}.pdb"]


def find_complex(d):
    for pd in PASSER_DIRS:
        p = REPO/"runs"/pd/f"{d}.pdb"
        if p.exists(): return p
    for t in EXTRA:
        p = REPO/"runs"/t.format(d=d)
        if p.exists(): return p
    return None


def atoms_by_chain(pdb):
    ch = {}
    for l in open(pdb):
        if l[:4] == "ATOM":
            c = l[21]; res = (l[22:27], l[17:20].strip())
            try: xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError: continue
            ch.setdefault(c, []).append((res, xyz))
    return ch


def interface_hyd(pdb):
    ch = atoms_by_chain(pdb)
    if len(ch) < 2: return None
    chains = sorted(ch, key=lambda c: len({r for r, _ in ch[c]}))
    binder, target = chains[0], chains[-1]
    tgt = [xyz for _, xyz in ch[target]]
    # binder residues with any heavy atom within 5A of any target heavy atom
    iface = {}
    for res, (x, y, z) in ch[binder]:
        if res in iface and iface[res]: continue
        hit = any((x-tx)**2 + (y-ty)**2 + (z-tz)**2 <= 25.0 for tx, ty, tz in tgt)
        iface[res] = iface.get(res, False) or hit
    inter = [res for res, h in iface.items() if h]
    if not inter: return None
    nhyd = sum(1 for (_, name) in inter if name in HYD)
    return nhyd/len(inter), len(inter)


def main():
    mani = sys.argv[1] if len(sys.argv) > 1 else str(REPO/"results"/"submission_manifest_v2.csv")
    rows = list(csv.DictReader(open(mani)))
    print(f"{'#':>2} {'source':10s} {'ifaceHyd':>8s} {'nIface':>6s} {'design'}")
    vals = []
    for r in rows:
        d = r["design"]; p = find_complex(d)
        if not p: print(f"{r.get('rank','?'):>2} {r['source'][:10]:10s} {'NOCOMPLEX':>8s}      - {d[:40]}"); continue
        res = interface_hyd(p)
        if not res: print(f"{r.get('rank','?'):>2} {r['source'][:10]:10s} {'?':>8s}      - {d[:40]}"); continue
        fr, n = res; vals.append(fr)
        flag = " <== hydrophobic interface" if fr > 0.45 else ""
        print(f"{r.get('rank','?'):>2} {r['source'][:10]:10s} {fr:8.3f} {n:6d} {d[:40]}{flag}")
    if vals:
        import statistics
        print(f"\ninterface hydrophobic fraction: mean={statistics.mean(vals):.3f} "
              f"median={statistics.median(vals):.3f} max={max(vals):.3f}")
        print(f"  >0.45 (hydrophobic-interface risk): {sum(v>0.45 for v in vals)}/{len(vals)}")
        print(f"  natural PPIs ~0.3-0.5; de novo AF2-optimized binders often skew high")


if __name__ == "__main__":
    main()
