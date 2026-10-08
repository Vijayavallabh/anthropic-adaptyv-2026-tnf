#!/usr/bin/env python3
"""WAVE 23 downstream stage 1: cooperative two-His RFd3 designs -> staged PDB (binder=A, target=B) +
BOTH motif His locked for MPNN.

Each design carries two binder His: His_a salt-bridging conserved Asp460 (target ordinal 127) and His_b
salt-bridging conserved Asp344 (ordinal 35). We keep BOTH when MPNN redesigns the binder, so we find the
binder His nearest each of the two target Asp carboxylates and lock both via fixed_residues_multi
({pdb: "A<a> A<b>"}). QC: drop designs where either bridge is broken (His N .. Asp O > cutoff) or the two
His collapse onto the same residue.
"""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
from moonshot_prep_mpnn import load_cif, write_pdb

TARGET_LEN = 191
HIS = ("HIS", "HIP", "HIE", "HID")
ASP_A_ORD = 127   # Asp460
ASP_B_ORD = 35    # Asp344


def nearest_his_to_asp(chains, bnd, tgt, asp_ord):
    asp = [r for r in chains[tgt] if r.id[1] == asp_ord and r.resname == "ASP"]
    if not asp:
        return None
    best = None
    for hr in chains[bnd]:
        if hr.resname not in HIS:
            continue
        for a in ("ND1", "NE2"):
            if a not in hr:
                continue
            for o in ("OD1", "OD2"):
                if o in asp[0]:
                    d = float(np.linalg.norm(hr[a].coord - asp[0][o].coord))
                    if best is None or d < best[0]:
                        best = (d, hr.id[1])
    return best


def analyze(model, bridge_cut):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    ha = nearest_his_to_asp(chains, bnd, tgt, ASP_A_ORD)
    hb = nearest_his_to_asp(chains, bnd, tgt, ASP_B_ORD)
    if ha is None or hb is None:
        return None
    if ha[1] == hb[1]:                     # same His for both -> not a two-site motif
        return None
    if ha[0] > bridge_cut or hb[0] > bridge_cut:
        return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), his_a=int(ha[1]), his_b=int(hb[1]),
                bridge_a=round(ha[0], 2), bridge_b=round(hb[0], 2), chains=chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/wave23_coop/gen")
    ap.add_argument("--out", default="runs/wave23_coop/staged")
    ap.add_argument("--bridge-cut", type=float, default=4.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cifs = sorted(glob.glob(os.path.join(a.gen, "gpu*", "*.cif.gz")))
    print(f"found {len(cifs)} RFd3 CIFs")
    pdbs, fixed, rows, kept, broken = [], {}, [], 0, 0
    for cf in cifs:
        name = os.path.basename(os.path.dirname(cf)) + "_" + os.path.basename(cf).replace(".cif.gz", "")
        try:
            info = analyze(load_cif(cf), a.bridge_cut)
        except Exception as e:
            print(f"  parse FAIL {name}: {repr(e)[:70]}"); continue
        if info is None:
            broken += 1; continue
        out_pdb = os.path.abspath(os.path.join(a.out, f"{name}.pdb"))
        write_pdb(info, out_pdb)
        pdbs.append(out_pdb)
        fixed[out_pdb] = f"A{info['his_a']} A{info['his_b']}"   # lock BOTH His
        rows.append(dict(name=name, binder_len=info["nbind"], his_a=info["his_a"], his_b=info["his_b"],
                         bridge_a=info["bridge_a"], bridge_b=info["bridge_b"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "his_a", "his_b", "bridge_a", "bridge_b"])
        w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (both bridges <= {a.bridge_cut} A, distinct His); dropped {broken}")
    if rows:
        import statistics as st
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; "
              f"bridge_a median {st.median(r['bridge_a'] for r in rows):.2f}A; "
              f"bridge_b median {st.median(r['bridge_b'] for r in rows):.2f}A")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
