#!/usr/bin/env python3
"""WAVE 25 downstream stage 1: de-novo His433-face CIFs -> staged PDB (binder=A, target=B) + a per-residue
HIS-enrichment bias for MPNN, placing His at binder positions near pKa-RAISING NEUTRAL acceptors.

Goal = a robust multi-His collective Wyman switch (Hill n>1): several binder His whose pKa shifts UP on binding
via charge-assisted H-bonds to NEUTRAL acceptors (so they stack like salt bridges but, having a neutral partner,
are far more ECC-robust). We bias His only at interface positions within reach of the conserved neutral-acceptor
patch -- Ser464 (ord 131, OG), Gln432 (ord 99, OE1 + backbone O), His433 (ord 100, patch) -- and Phe412 (ord 103,
a cation-pi bonus). Emits bias_AA_per_residue_multi {pdb: {"A<pos>": {"H": +BIAS}}} + 0-Cys handled globally.
QC: keep designs with >=2 biasable interface positions (multi-His must be possible). No fixed motif (de novo).
"""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
from moonshot_prep_mpnn import load_cif, write_pdb

TARGET_LEN = 191
HIS_BIAS = 2.0
IFACE_CUT = 8.0     # binder CB within this of a neutral-acceptor atom -> bias His there
MIN_SITES = 2       # need >=2 biasable positions for a multi-His switch
# neutral-acceptor atoms on the conserved His433 face (target ordinals = mature - 309):
ACCEPTORS = {131: ["OG", "O"],          # Ser464 hydroxyl + backbone O
             99: ["OE1", "NE2", "O"],   # Gln432 amide + backbone O
             100: ["ND1", "NE2", "O"],  # His433 patch + backbone O
             103: ["O"]}                 # Phe412 backbone O (ring handled as cation-pi bonus hotspot)


def acceptor_coords(tgt_res):
    coords = []
    for r in tgt_res:
        if r.id[1] in ACCEPTORS:
            for a in ACCEPTORS[r.id[1]]:
                if a in r:
                    coords.append(r[a].coord)
    return np.array(coords) if coords else np.zeros((0, 3))


def analyze(model):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    acc = acceptor_coords(chains[tgt])
    if len(acc) == 0:
        return None
    sites = []
    for r in chains[bnd]:
        cb = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
        if cb is None:
            continue
        d = float(np.linalg.norm(acc - cb.coord, axis=1).min())
        if d <= IFACE_CUT:
            sites.append((r.id[1], round(d, 2)))
    if len(sites) < MIN_SITES:
        return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), sites=sites, chains=chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/wave25_multihis/gen")
    ap.add_argument("--out", default="runs/wave25_multihis/staged")
    ap.add_argument("--max-sites", type=int, default=5, help="cap His-biased positions per design (nearest first)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cifs = sorted(glob.glob(os.path.join(a.gen, "gpu*", "*.cif.gz")))
    print(f"found {len(cifs)} de-novo CIFs")
    pdbs, bias, rows, kept, dropped = [], {}, [], 0, 0
    for cf in cifs:
        name = os.path.basename(os.path.dirname(cf)) + "_" + os.path.basename(cf).replace(".cif.gz", "")
        try:
            info = analyze(load_cif(cf))
        except Exception as e:
            print(f"  parse FAIL {name}: {repr(e)[:70]}"); continue
        if info is None:
            dropped += 1; continue
        out_pdb = os.path.abspath(os.path.join(a.out, f"{name}.pdb"))
        write_pdb(info, out_pdb)
        pdbs.append(out_pdb)
        sites = sorted(info["sites"], key=lambda s: s[1])[:a.max_sites]   # nearest positions
        bias[out_pdb] = {f"A{pos}": {"H": HIS_BIAS} for pos, _ in sites}
        rows.append(dict(name=name, binder_len=info["nbind"], n_sites=len(sites),
                         sites=";".join(f"{p}@{d}" for p, d in sites)))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(bias, open(os.path.join(a.out, "bias_AA_per_residue_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "n_sites", "sites"]); w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (>= {MIN_SITES} interface sites near neutral acceptors); dropped {dropped}")
    if rows:
        import statistics as st
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; "
              f"median His-bias sites/design {st.median(r['n_sites'] for r in rows):.0f}")
    print(f"-> {a.out}/pdbs.json, bias_AA_per_residue_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
