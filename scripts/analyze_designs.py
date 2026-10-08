#!/usr/bin/env python3
"""Cheap structural analysis + pre-filter of generated binder-target complexes (no folding, no PyRosetta).

Parses each BoltzGen inverse-folded complex CIF and computes signals for all three objectives
straight from the generated complex geometry, so we can rank thousands of designs and only
re-fold / pH-score the promising few hundred:

  Obj1 (pH switch)     : binder Asp/Glu carboxyl O within reach of target His433 imidazole N
                         (a protonated-His salt-bridge partner) + binder His count near interface
  Obj2 (cross-react)   : # conserved-anchor contacts (target ordinals 75,99,100,103,131,156,163)
                         vs # divergent contacts (44,109,134,158,159,164); His433=100, Glu496=163
  Obj3 (affinity proxy): interface residue count + binder composition (GLU/LYS)
  novelty              : sequence identity to ubiquitin (BoltzGen's short-length collapse mode)

Target chain in BoltzGen output is renumbered 1..191, so target pdbnum == crop ordinal.
Run: uv run python scripts/analyze_designs.py --runroot runs/gen1 [--runroot ...] --out results/designs_analyzed.csv
"""
from __future__ import annotations

import argparse
import glob
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]

CONSERVED = {75, 99, 100, 103, 131, 156, 163}   # target ordinals of conserved anchors
DIVERGENT = {44, 109, 134, 158, 159, 164}       # divergent cetuximab contacts to avoid
HIS433 = 100
GLU496 = 163
UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
D3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
ACID_O = {"OD1", "OD2", "OE1", "OE2"}     # Asp/Glu carboxyl oxygens
HIS_N = {"ND1", "NE2"}                       # His imidazole nitrogens


def parse_cif_atoms(path):
    header, rows, inloop = [], [], False
    for line in open(path):
        ts = line.strip()
        if ts.startswith("_atom_site."):
            header.append(ts.split(".", 1)[1]); inloop = True; continue
        if inloop:
            if ts[:4] == "ATOM" or ts[:6] == "HETATM":
                f = ts.split()
                if len(f) >= len(header):
                    rows.append(f)
            elif ts == "" or ts[:1] == "#" or ts.startswith("loop_") or ts[:1] == "_":
                if rows:
                    break
    if not rows:
        return None
    idx = {name: i for i, name in enumerate(header)}

    def col(*names):
        for n in names:
            if n in idx:
                return idx[n]
        return None
    ci = col("label_asym_id", "auth_asym_id")
    ni = col("auth_seq_id", "label_seq_id")
    ri = col("label_comp_id", "auth_comp_id")
    ai = col("label_atom_id", "auth_atom_id")
    xi, yi, zi = idx.get("Cartn_x"), idx.get("Cartn_y"), idx.get("Cartn_z")
    if None in (ci, ni, ri, ai, xi, yi, zi):
        return None
    atoms = []
    for f in rows:
        try:
            atoms.append((f[ci], int(f[ni]), f[ri], f[ai].strip('"'),
                          float(f[xi]), float(f[yi]), float(f[zi])))
        except (ValueError, IndexError):
            continue
    return atoms


def analyze(path):
    atoms = parse_cif_atoms(path)
    if not atoms:
        return None
    # chains by residue count
    chain_res = {}
    for ch, num, resn, atn, x, y, z in atoms:
        chain_res.setdefault(ch, set()).add(num)
    if len(chain_res) < 2:
        return None
    binder = min(chain_res, key=lambda c: len(chain_res[c]))
    target = max(chain_res, key=lambda c: len(chain_res[c]))

    b_xyz = np.array([[x, y, z] for ch, n, rn, an, x, y, z in atoms if ch == binder])
    # target residue -> atoms
    t_res = {}
    for ch, n, rn, an, x, y, z in atoms:
        if ch == target:
            t_res.setdefault((n, rn), []).append((an, x, y, z))
    # binder residue seq (CA order) + composition
    b_ca = sorted([(n, rn) for ch, n, rn, an, x, y, z in atoms if ch == binder and an == "CA"])
    seq = "".join(D3.get(rn, "X") for _, rn in b_ca)
    L = len(seq) or 1

    # interface: target residues within 5A of any binder atom
    contacts = []
    for (n, rn), al in t_res.items():
        ta = np.array([[x, y, z] for _, x, y, z in al])
        d = np.sqrt(((ta[:, None, :] - b_xyz[None, :, :]) ** 2).sum(-1)).min()
        if d < 5.0:
            contacts.append(n)
    cset = set(contacts)

    # Obj1: binder Asp/Glu carboxyl O nearest to His433 imidazole N
    his_atoms = [(an, x, y, z) for (n, rn), al in t_res.items() if n == HIS433 and rn == "HIS"
                 for (an, x, y, z) in al if an in HIS_N]
    b_acid = np.array([[x, y, z] for ch, n, rn, an, x, y, z in atoms
                       if ch == binder and rn in ("ASP", "GLU") and an in ACID_O])
    if his_atoms and len(b_acid):
        hn = np.array([[x, y, z] for _, x, y, z in his_atoms])
        acid_his_mindist = float(np.sqrt(((b_acid[:, None, :] - hn[None, :, :]) ** 2).sum(-1)).min())
    else:
        acid_his_mindist = math.inf

    def frac(aa):
        return round(seq.count(aa) / L, 4)

    ident_ubq = round(sum(1 for i in range(min(len(seq), len(UBQ))) if seq[i] == UBQ[i]) / min(len(seq), len(UBQ)), 3) if seq else 0.0

    return {
        "design_id": Path(path).stem,
        "cif": str(path),
        "binder_chain": binder, "target_chain": target,
        "binder_len": len(seq), "binder_seq": seq,
        "n_interface_target": len(contacts),
        "n_conserved_contacts": len(cset & CONSERVED),
        "conserved_hit": ",".join(str(o) for o in sorted(cset & CONSERVED)),
        "his433_contact": HIS433 in cset,
        "glu496_contact": GLU496 in cset,
        "n_divergent_contacts": len(cset & DIVERGENT),
        "acid_near_his433_ang": round(acid_his_mindist, 2) if acid_his_mindist != math.inf else "",
        "ph_partner": acid_his_mindist < 5.0,
        "frac_H": frac("H"), "frac_E": frac("E"), "frac_D": frac("D"),
        "frac_K": frac("K"), "frac_R": frac("R"), "frac_C": frac("C"),
        "ident_to_ubiquitin": ident_ubq,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runroot", action="append", required=True, help="run root(s); repeatable")
    ap.add_argument("--out", default=str(REPO / "results" / "designs_analyzed.csv"))
    a = ap.parse_args()

    cifs = []
    for root in a.runroot:
        found = glob.glob(os.path.join(root, "**", "intermediate_designs_inverse_folded", "*.cif"), recursive=True)
        cifs += [c for c in found if "refold" not in os.path.basename(c) and "/refold" not in c]
    cifs = sorted(set(cifs))
    print(f"analyzing {len(cifs)} complex CIFs ...")

    rows, fail = [], 0
    for i, c in enumerate(cifs):
        try:
            r = analyze(c)
        except Exception:
            r = None
        if r is None:
            fail += 1
            continue
        # provenance from path .../<spec>/<seed>/intermediate_designs_inverse_folded/<id>.cif
        parts = Path(c).parts
        r["source_spec"] = parts[-4] if len(parts) >= 4 else ""
        r["source_seed"] = parts[-3] if len(parts) >= 3 else ""
        rows.append(r)
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(cifs)}")

    if not rows:
        raise SystemExit("no designs analyzed")
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)

    print(f"\n{len(df)} designs ({fail} failed) -> {a.out}")
    print(f"  contact His433: {df['his433_contact'].sum()}  ({100*df['his433_contact'].mean():.0f}%)")
    print(f"  pH-partner (binder acid <5A of His433): {df['ph_partner'].sum()}")
    print(f"  novel (id_ubq<0.5): {(df['ident_to_ubiquitin']<0.5).sum()}")
    print(f"  conserved>=3 & divergent<=1: {((df['n_conserved_contacts']>=3)&(df['n_divergent_contacts']<=1)).sum()}")
    print(f"  mean conserved_contacts={df['n_conserved_contacts'].mean():.2f} divergent={df['n_divergent_contacts'].mean():.2f}")


if __name__ == "__main__":
    main()
