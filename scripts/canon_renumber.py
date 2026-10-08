#!/usr/bin/env python3
"""Renumber a binder-target complex to the canonical form obj1_switch expects:
binder = chain A (residues 1..Lb), target = chain B (residues 1..Lt), so His433 = target's
100th residue = B/100. Binder = the SMALLER chain, target = the LARGER (by residue count).
Safe no-op if the structure is already canonical. Use on AF2-ig af2pred outputs before obj1_switch
(af2ig may renumber the target to continuous numbering, which would break the B/100 His433 locator).

Usage:
  canon_renumber.py in.pdb out.pdb
  canon_renumber.py --indir DIR --outdir DIR [--glob '*af2pred.pdb']
"""
import argparse, glob, os
from collections import OrderedDict


def chain_order_and_sizes(path):
    order = []
    nres = {}
    seen = set()
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            ch = ln[21]
            key = (ch, ln[22:27])
            if ch not in nres:
                order.append(ch); nres[ch] = 0
            if key not in seen:
                seen.add(key); nres[ch] += 1
    return order, nres


def canon(in_pdb, out_pdb):
    order, nres = chain_order_and_sizes(in_pdb)
    prot = [c for c in order]
    if len(prot) < 2:
        # nothing to do; copy through
        txt = open(in_pdb).read()
        open(out_pdb, "w").write(txt)
        return None
    target = max(prot, key=lambda c: nres[c])
    binder = min(prot, key=lambda c: nres[c])
    groups = {"binder": [], "target": []}
    for ln in open(in_pdb):
        if ln.startswith(("ATOM", "HETATM")):
            ch = ln[21]
            if ch == binder:
                groups["binder"].append(ln)
            elif ch == target:
                groups["target"].append(ln)
    out = []
    serial = 1
    lb = None
    for grp, newch in (("binder", "A"), ("target", "B")):
        resmap = {}
        nr = 0
        for ln in groups[grp]:
            old = ln[22:27]
            if old not in resmap:
                nr += 1; resmap[old] = nr
            out.append("%s%5d %s%s%4d %s\n" % (
                ln[:6], serial, ln[12:21], newch, resmap[old], ln[27:].rstrip("\n")))
            serial += 1
        out.append("TER\n")
        if grp == "binder":
            lb = nr
    out.append("END\n")
    open(out_pdb, "w").write("".join(out))
    return lb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("args", nargs="*")
    ap.add_argument("--indir")
    ap.add_argument("--outdir")
    ap.add_argument("--glob", default="*.pdb")
    a = ap.parse_args()
    if a.indir:
        os.makedirs(a.outdir, exist_ok=True)
        n = 0
        for p in sorted(glob.glob(os.path.join(a.indir, a.glob))):
            canon(p, os.path.join(a.outdir, os.path.basename(p)))
            n += 1
        print(f"canonicalized {n} -> {a.outdir}")
    elif len(a.args) == 2:
        lb = canon(a.args[0], a.args[1])
        print(f"binder_len={lb} -> {a.args[1]}")
    else:
        ap.error("give in.pdb out.pdb, or --indir --outdir")


if __name__ == "__main__":
    main()
