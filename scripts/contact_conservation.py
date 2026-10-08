#!/usr/bin/env python3
"""Artifact-free cross-reactivity proxy: on each engineered binder's AF2-ig structure, find the
target residues it contacts (<5 A) and classify them as conserved (human=mouse) or divergent.
A binder whose interface lands on conserved anchors (and avoids divergent cetuximab contacts)
will cross-react to mouse EGFR DIII.

af2ig renumbers chains continuously; target ordinal = chainB_pdbnum - (min chainB pdbnum - 1).
Conserved anchors (ordinals): 75,99,100,103,131,156,163 (His433=100, Glu496=163).
Divergent (avoid):            44,109,134,158,159,164.
"""
from __future__ import annotations
import argparse, csv, glob, math, os
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
CONSERVED = {75, 99, 100, 103, 131, 156, 163}
DIVERGENT = {44, 109, 134, 158, 159, 164}
HIS433 = 100


def analyze(pdb):
    chains = {}
    atoms = {"A": [], "B": []}
    tbynum = {}
    for l in open(pdb):
        if l[:4] != "ATOM":
            continue
        ch = l[21]
        try:
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
            num = int(l[22:26])
        except ValueError:
            continue
        if l[12:16].strip() == "CA":
            chains.setdefault(ch, 0)
            chains[ch] += 1
        atoms.setdefault(ch, []).append((num, x, y, z))
    if "A" not in atoms or "B" not in atoms:
        return None
    # binder = fewer CA, target = more; here A=binder, B=target by construction
    binder, target = "A", "B"
    if chains.get("A", 0) > chains.get("B", 0):
        binder, target = "B", "A"
    b_xyz = np.array([[x, y, z] for _, x, y, z in atoms[binder]])
    offset = min(n for n, *_ in atoms[target]) - 1
    tnums = {}
    for n, x, y, z in atoms[target]:
        tnums.setdefault(n, []).append((x, y, z))
    contacts = set()
    for n, xs in tnums.items():
        ta = np.array(xs)
        d = np.sqrt(((ta[:, None, :] - b_xyz[None, :, :]) ** 2).sum(-1)).min()
        if d < 5.0:
            contacts.add(n - offset)          # ordinal
    return {
        "n_conserved": len(contacts & CONSERVED),
        "n_divergent": len(contacts & DIVERGENT),
        "his433_contact": HIS433 in contacts,
        "conserved_hits": ",".join(str(o) for o in sorted(contacts & CONSERVED)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default=str(REPO / "runs" / "af2ig_eng" / "gpu*" / "pdbs" / "*_af2pred.pdb"))
    ap.add_argument("--out", default=str(REPO / "results" / "contact_conservation.csv"))
    a = ap.parse_args()
    rows = []
    for p in sorted(glob.glob(a.glob)):
        design = Path(p).stem.replace("_af2pred", "").replace("__best", "")
        r = analyze(p)
        if r:
            r["design"] = design
            rows.append(r)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["design", "n_conserved", "n_divergent", "his433_contact", "conserved_hits"])
        w.writeheader()
        w.writerows(rows)
    xr = [r for r in rows if r["n_conserved"] >= 3 and r["n_divergent"] <= 1]
    print(f"analyzed {len(rows)} | cross-reactive (conserved>=3 & divergent<=1): {len(xr)}")
    for r in sorted(rows, key=lambda r: (-r["n_conserved"], r["n_divergent"])):
        print(f"  {r['design'][:46]:46s} cons={r['n_conserved']} div={r['n_divergent']} his433={r['his433_contact']} [{r['conserved_hits']}]")


if __name__ == "__main__":
    main()
