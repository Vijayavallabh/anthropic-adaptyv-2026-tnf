#!/usr/bin/env python3
"""WAVE 28 downstream stage 1: carboxylate-dyad RFd3 .cif.gz -> staged PDB (binder=A, target=B) + Glu-lock map.

Each design carries a binder GLU seeded to form a buried dyad with the conserved target Asp460 (mature A436
-> CIF ordinal 436-310+1 = 127). We find the binder GLU whose carboxyl (OE1/OE2) is nearest the Asp460
carboxyl (OD1/OD2), LOCK it for MPNN (fixed_residues_multi {pdb: "A<glu>"}), and QC-drop designs where the
dyad wandered. Inverts moonshot_prep_trp_his433_stage (Glu carboxyl dyad instead of Trp...His cation-pi)."""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
from moonshot_prep_mpnn import load_cif, write_pdb

TARGET_LEN = 191
ASP460_ORD = 127   # Asp460 (P00533 precursor) = mature 436 -> ordinal 436-310+1
ASP_O = ("OD1", "OD2")
GLU_O = ("OE1", "OE2")


def min_oo(res_a, names_a, res_b, names_b):
    oa = [res_a[n].coord for n in names_a if n in res_a]
    ob = [res_b[n].coord for n in names_b if n in res_b]
    if not oa or not ob:
        return None
    return float(min(np.linalg.norm(x - y) for x in oa for y in ob))


def nearest_glu(binder_res, asp):
    best = None
    for gr in binder_res:
        if gr.resname != "GLU":
            continue
        oo = min_oo(gr, GLU_O, asp, ASP_O)
        if oo is None:
            continue
        if best is None or oo < best[0]:
            best = (oo, gr.id[1])
    return best


def analyze(model, cc_cut):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    asp = [r for r in chains[tgt] if r.id[1] == ASP460_ORD and r.resname == "ASP"]
    if not asp:
        return None
    glu = nearest_glu(chains[bnd], asp[0])
    if glu is None or glu[0] > cc_cut:
        return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), glu=int(glu[1]), oo=round(glu[0], 2), chains=chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/wave28_dyad/gen")
    ap.add_argument("--out", default="runs/wave28_dyad/staged")
    ap.add_argument("--cc-cut", type=float, default=3.6, help="max Glu-carboxyl..Asp460-carboxyl O..O to keep")
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
            print(f"  parse FAIL {name}: {repr(e)[:70]}"); continue
        if info is None:
            dropped += 1; continue
        out_pdb = os.path.abspath(os.path.join(a.out, f"{name}.pdb"))
        write_pdb(info, out_pdb)
        pdbs.append(out_pdb)
        fixed[out_pdb] = f"A{info['glu']}"
        rows.append(dict(name=name, binder_len=info["nbind"], glu=info["glu"], oo=info["oo"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "glu", "oo"])
        w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (dyad O..O within {a.cc_cut} A); dropped {dropped}")
    if rows:
        import statistics as st
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; "
              f"dyad O..O median {st.median(r['oo'] for r in rows):.2f}A")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
