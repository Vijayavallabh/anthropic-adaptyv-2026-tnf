#!/usr/bin/env python3
"""MD pre-step for the two-protomer groove: md_release_tnf.py picks binder/target by CHAIN SIZE (2 largest),
so a 3-chain complex (binder A + target B + target C) mis-identifies a protomer as the binder. Merge the two
target protomers into ONE chain so detect_chains() correctly sees (target=merged 304-res, binder=A), and
nearest_lys searches BOTH protomers for the conserved anchor (K65/K98/K112). Correct because the MD restrains
the whole target rigid. Relabels C->B and renumbers C residues +OFFSET to avoid resnum clash.
Usage: merge_protomers.py <in A+B+C .pdb> <out A+mergedB .pdb> [binder=A] [keep=B,C] [offset=200]
"""
import sys
inp, outp = sys.argv[1], sys.argv[2]
BINDER = sys.argv[3] if len(sys.argv) > 3 else "A"
KEEP = (sys.argv[4] if len(sys.argv) > 4 else "B,C").split(",")
OFFSET = int(sys.argv[5]) if len(sys.argv) > 5 else 200   # added to 2nd+ protomer resnums to avoid clash
out = []
for ln in open(inp):
    if ln.startswith("TER"):
        out.append("TER\n"); continue
    if not ln.startswith(("ATOM", "HETATM")) or len(ln) < 22:
        continue
    ch = ln[21]
    if ch == BINDER:
        out.append(ln)  # binder unchanged (chain A)
    elif ch in KEEP:
        idx = KEEP.index(ch)
        if idx == 0:
            out.append(ln)  # first protomer -> keep as-is (chain B)
        else:
            # relabel to the first KEEP chain and renumber +idx*OFFSET
            try:
                rn = int(ln[22:26]) + idx * OFFSET
            except ValueError:
                out.append(ln); continue
            out.append(f"{ln[:21]}{KEEP[0]}{rn:4d}{ln[26:]}")
with open(outp, "w") as fh:
    for ln in out:
        if ln.startswith("TER"):
            fh.write("TER\n")
        else:
            fh.write(ln)
    fh.write("END\n")
# report
import collections
chains = collections.Counter(l[21] for l in out if l.startswith(("ATOM", "HETATM")))
print(f"merged {inp} -> {outp}: chains { {k: v for k, v in chains.items()} } (binder {BINDER} + merged {KEEP[0]})")
