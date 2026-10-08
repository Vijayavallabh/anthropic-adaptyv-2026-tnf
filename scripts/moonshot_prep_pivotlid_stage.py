#!/usr/bin/env python3
"""WAVE 29 PivotLid stage: His...Asp460 staple RFd3 .cif.gz -> staged PDB (binder=A,target=B) + His-lock map.
Lock the binder HIS whose ND1/NE2 is nearest the conserved Asp460 carboxyl (CIF ordinal 127 = mature 436).
Clone of moonshot_prep_dyad_stage.py with the staple = His-N...Asp-O (not a Glu carboxyl dyad)."""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moonshot_prep_mpnn import load_cif, write_pdb

TARGET_LEN = 191
ASP460_ORD = 127
ASP_O = ("OD1", "OD2")
HIS_N = ("ND1", "NE2")


def min_d(res_a, na, res_b, nb):
    a = [res_a[n].coord for n in na if n in res_a]
    b = [res_b[n].coord for n in nb if n in res_b]
    if not a or not b: return None
    return float(min(np.linalg.norm(x - y) for x in a for y in b))


def nearest_his(binder_res, asp):
    best = None
    for hr in binder_res:
        if not hr.resname.startswith("HI"): continue
        d = min_d(hr, HIS_N, asp, ASP_O)
        if d is None: continue
        if best is None or d < best[0]: best = (d, hr.id[1])
    return best


def analyze(model, cc_cut):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders: return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    asp = [r for r in chains[tgt] if r.id[1] == ASP460_ORD and r.resname == "ASP"]
    if not asp: return None
    his = nearest_his(chains[bnd], asp[0])
    if his is None or his[0] > cc_cut: return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), his=int(his[1]), d=round(his[0], 2), chains=chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/wave29_pivotlid/gen")
    ap.add_argument("--out", default="runs/wave29_pivotlid/staged")
    ap.add_argument("--cc-cut", type=float, default=4.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cifs = sorted(glob.glob(os.path.join(a.gen, "gpu*", "*.cif.gz")))
    print(f"found {len(cifs)} RFd3 CIFs")
    pdbs, fixed, rows, kept, dropped = [], {}, [], 0, 0
    for cf in cifs:
        name = os.path.basename(os.path.dirname(cf)) + "_" + os.path.basename(cf).replace(".cif.gz", "")
        try:
            info = analyze(load_cif(cf), a.cc_cut)
        except Exception as e:
            print(f"  parse FAIL {name}: {repr(e)[:60]}"); continue
        if info is None:
            dropped += 1; continue
        out_pdb = os.path.abspath(os.path.join(a.out, f"{name}.pdb"))
        write_pdb(info, out_pdb)
        pdbs.append(out_pdb); fixed[out_pdb] = f"A{info['his']}"
        rows.append(dict(name=name, binder_len=info["nbind"], his=info["his"], d=info["d"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "his", "d"]); w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (His-N..Asp460 within {a.cc_cut} A); dropped {dropped}")
    if rows:
        import statistics as st
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; staple median {st.median(r['d'] for r in rows):.2f}A")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json")


if __name__ == "__main__":
    main()
