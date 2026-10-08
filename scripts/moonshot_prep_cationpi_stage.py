#!/usr/bin/env python3
"""WAVE 24 downstream stage 1: cation-pi RFd3 .cif.gz -> staged PDB (binder=A, target=B) + His-lock map.

Each design carries a binder His seeded to cation-pi-stack the conserved target Phe412. In the generated CIF
the target DIII is renumbered to ordinals 1..191 (ordinal = native_mature - 310 + 1), so Phe412 -> ordinal
103 and Phe380 -> ordinal 71. We find the binder His whose ring centroid is nearest the Phe412 ring, LOCK it
for MPNN (fixed_residues_multi {pdb: "A<his>"}), and also record the nearest binder His to Phe380 (to detect a
serendipitous DUAL later). QC: drop designs whose Phe412 His is beyond ring_cc cutoff (RFd3 occasionally
wanders off the seed). Reuses load_cif/write_pdb from moonshot_prep_mpnn.
"""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
from moonshot_prep_mpnn import load_cif, write_pdb

TARGET_LEN = 191
HIS = ("HIS", "HIP", "HIE", "HID")
PHE412_ORD = 103   # Phe436 precursor / mature 412 -> ordinal
PHE380_ORD = 71    # Phe404 precursor / mature 380 -> ordinal
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
PHE_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")


def centroid(res, names):
    cs = [res[a].coord for a in names if a in res]
    return np.mean(cs, axis=0) if len(cs) >= 3 else None


def nearest_his(binder_res, phe_centroid):
    best = None
    for hr in binder_res:
        if hr.resname not in HIS:
            continue
        hc = centroid(hr, HIS_RING)
        if hc is None:
            continue
        d = float(np.linalg.norm(hc - phe_centroid))
        if best is None or d < best[0]:
            best = (d, hr.id[1])
    return best


def analyze(model, cc_cut):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    tgt_res = chains[tgt]
    phe412 = [r for r in tgt_res if r.id[1] == PHE412_ORD and r.resname == "PHE"]
    phe380 = [r for r in tgt_res if r.id[1] == PHE380_ORD and r.resname == "PHE"]
    if not phe412:
        return None
    c412 = centroid(phe412[0], PHE_RING)
    h412 = nearest_his(chains[bnd], c412)
    if h412 is None or h412[0] > cc_cut:
        return None
    h380 = None
    if phe380:
        c380 = centroid(phe380[0], PHE_RING)
        h380 = nearest_his(chains[bnd], c380)
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), his412=int(h412[1]), cc412=round(h412[0], 2),
                his380=(int(h380[1]) if h380 else None), cc380=(round(h380[0], 2) if h380 else None),
                chains=chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/wave24_cationpi/gen")
    ap.add_argument("--out", default="runs/wave24_cationpi/staged")
    ap.add_argument("--cc-cut", type=float, default=6.5, help="max His-ring..Phe412-ring centroid dist to keep")
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
        fixed[out_pdb] = f"A{info['his412']}"           # lock the Phe412 cation-pi His
        rows.append(dict(name=name, binder_len=info["nbind"], his412=info["his412"], cc412=info["cc412"],
                         his380=info["his380"], cc380=info["cc380"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "his412", "cc412", "his380", "cc380"])
        w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (Phe412 His within {a.cc_cut} A); dropped {dropped}")
    if rows:
        import statistics as st
        near380 = [r for r in rows if r["cc380"] is not None and r["cc380"] <= 6.5]
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; "
              f"Phe412 cc median {st.median(r['cc412'] for r in rows):.2f}A; "
              f"designs with a His also near Phe380 (<=6.5A, potential DUAL): {len(near380)}")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
