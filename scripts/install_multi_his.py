#!/usr/bin/env python3
"""Install His at a LIST of binder (chain A) positions on a complex PDB (keep bb N/CA/C/O + CB, rename
to HIS so af2ig/MPNN rebuilds the ring), preserving all other atoms incl. any existing His. For building
COOPERATIVE multi-His pH-release switches (Baker: every working switch had >=2 His-cation contacts).
Usage: install_multi_his.py <complex.pdb> <out.pdb> <posA,posB[,posC]>  (positions are chain-A resnums)
"""
import sys
KEEP = {"N", "CA", "C", "O", "CB"}
pdb, out, poss = sys.argv[1], sys.argv[2], {int(p) for p in sys.argv[3].split(",")}
lines = []
for ln in open(pdb):
    if ln[:4] != "ATOM":
        lines.append(ln); continue
    ch = ln[21]; rn = int(ln[22:26]); atom = ln[12:16].strip()
    if ch == "A" and rn in poss:
        if atom not in KEEP:
            continue                       # drop sidechain beyond CB
        ln = ln[:17] + "HIS" + ln[20:]     # rename residue -> HIS
    lines.append(ln)
open(out, "w").write("".join(lines))
print(f"installed HIS at chain-A {sorted(poss)} -> {out}")
