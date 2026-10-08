#!/usr/bin/env python3
"""WAVE 27 downstream stage 1: Trp...His433 cation-pi RFd3 .cif.gz -> staged PDB (binder=A, target=B) +
Trp-lock map.

Each design carries a binder TRP seeded to cation-pi-stack the conserved target His433 (the proven pH
sensor; the CATION when protonated at tumor pH). In the generated CIF the target DIII is renumbered to
ordinals 1..191 (ordinal = native_mature - 310 + 1), so His433 (mature 409) -> ordinal 100. We find the
binder TRP whose indole centroid is nearest the His433 ring, LOCK it for MPNN
(fixed_residues_multi {pdb: "A<trp>"}), and QC-drop designs whose Trp wandered off the sensor. Inverts
scripts/moonshot_prep_cationpi_stage.py (aromatic now on the binder, cation on the target). Reuses
load_cif/write_pdb from moonshot_prep_mpnn."""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moonshot_prep_mpnn import load_cif, write_pdb

TARGET_LEN = 191
HIS433_ORD = 100   # His433 (P00533 precursor) = mature 409 -> ordinal 409-310+1
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
INDOLE = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")


def centroid(res, names):
    cs = [res[a].coord for a in names if a in res]
    return np.mean(cs, axis=0) if len(cs) >= 3 else None


def nearest_trp(binder_res, his_centroid):
    best = None
    for tr in binder_res:
        if tr.resname != "TRP":
            continue
        tc = centroid(tr, INDOLE)
        if tc is None:
            continue
        d = float(np.linalg.norm(tc - his_centroid))
        if best is None or d < best[0]:
            best = (d, tr.id[1])
    return best


def analyze(model, cc_cut):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    tgt_res = chains[tgt]
    his433 = [r for r in tgt_res if r.id[1] == HIS433_ORD and r.resname == "HIS"]
    if not his433:
        return None
    chis = centroid(his433[0], HIS_RING)
    if chis is None:
        return None
    trp = nearest_trp(chains[bnd], chis)
    if trp is None or trp[0] > cc_cut:
        return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), trp=int(trp[1]), cc=round(trp[0], 2), chains=chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/wave27_trp_his433/gen")
    ap.add_argument("--out", default="runs/wave27_trp_his433/staged")
    ap.add_argument("--cc-cut", type=float, default=6.5, help="max Trp-indole..His433-ring centroid dist to keep")
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
        fixed[out_pdb] = f"A{info['trp']}"              # lock the His433 cation-pi Trp
        rows.append(dict(name=name, binder_len=info["nbind"], trp=info["trp"], cc=info["cc"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "trp", "cc"])
        w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (His433 Trp within {a.cc_cut} A); dropped {dropped}")
    if rows:
        import statistics as st
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; "
              f"His433 cc median {st.median(r['cc'] for r in rows):.2f}A")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
