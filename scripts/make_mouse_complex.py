#!/usr/bin/env python3
"""Build binder+MOUSE-DIII complexes for cross-reactivity fold-validation.

Threads the mouse Domain-III sequence (1:1 aligned to the human crop; data/targets/
mouse_domainIII_map.json) onto chain B of each human complex PDB, keeping backbone coords and
dropping side chains (af2ig rebuilds them). Chain A (binder) is kept. Human/mouse DIII are 87%
identical with ALL epitope anchors conserved (His433 included), so a correctly cross-reactive
binder should fold-validate against both with similar pae_interaction.

Usage: python scripts/make_mouse_complex.py --indir <dir of human complexes> --outdir <dir>
"""
from __future__ import annotations
import argparse, glob, json, os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AA1to3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
          "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
          "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
BB = {"N", "CA", "C", "O"}


def mouse_seq():
    m = json.load(open(REPO / "data" / "targets" / "mouse_domainIII_map.json"))
    return m["mouse_seq"]


def thread_mouse(src, dst, mouse):
    out = []
    for l in open(src):
        if l[:4] == "ATOM" or l[:6] == "HETATM":
            ch = l[21]
            if ch == "B":                        # target
                atom = l[12:16].strip()
                if atom not in BB:
                    continue                     # drop side chain; af2ig rebuilds
                if mouse is None:                # human control: keep original residue, bb only
                    out.append(l); continue
                try:
                    num = int(l[22:26])
                except ValueError:
                    out.append(l); continue
                if 1 <= num <= len(mouse):
                    out.append(l[:17] + AA1to3[mouse[num - 1]] + l[20:])
                else:
                    out.append(l)
            else:                                # binder chain A: keep backbone only (af2ig rebuilds)
                atom = l[12:16].strip()
                if atom in BB:
                    out.append(l)
        elif l[:3] == "TER":
            out.append(l)
    out.append("END\n")
    with open(dst, "w") as fh:
        fh.writelines(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--target", choices=["mouse", "human"], default="mouse",
                    help="human = control: same backbone-only pipeline but keep the human target")
    a = ap.parse_args()
    mouse = mouse_seq() if a.target == "mouse" else None
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    n = 0
    for pdb in sorted(glob.glob(os.path.join(a.indir, "*.pdb"))):
        dst = out / (Path(pdb).stem + "_mouse.pdb")
        thread_mouse(pdb, dst, mouse)
        n += 1
    print(f"wrote {n} mouse complexes -> {out}")


if __name__ == "__main__":
    main()
