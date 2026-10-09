#!/usr/bin/env python3
"""Thread N His point-mutations (bb+CB kept, sidechain rebuilt by af2ig) onto a reference
binder-target complex PDB. For building buried-His NETWORK candidates (multi-His) from
individually-validated single-His positions. Usage: apply_his_mutations.py <ref.pdb> <pos1,pos2,..> <out.pdb> [binder_chain=A]"""
import sys
ref, poss, out = sys.argv[1], sys.argv[2], sys.argv[3]
bchain = sys.argv[4] if len(sys.argv) > 4 else "A"
mut = set(int(x) for x in poss.split(","))
keep = {"N", "CA", "C", "O", "CB"}
o = []
for ln in open(ref):
    if ln[:4] == "ATOM" and ln[21] == bchain and int(ln[22:26]) in mut:
        if ln[12:16].strip() in keep:
            o.append(ln[:17] + "HIS" + ln[20:])   # rename residue -> HIS, keep bb+CB only
        # drop other sidechain atoms (af2ig rebuilds from sequence)
    else:
        o.append(ln)
open(out, "w").write("".join(o))
print(f"wrote {out}: His at {sorted(mut)} on chain {bchain}")
