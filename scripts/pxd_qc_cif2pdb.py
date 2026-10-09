#!/usr/bin/env python3
"""Convert the 20 uniform re-fold human complexes (Boltz cif) -> PDB with chain A=binder, B=target.
Run with the boltz2 env python (has gemmi). Writes to runs/push/refold_pdb/<design>.pdb."""
import csv, os, glob, gemmi
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
out="runs/push/refold_pdb"; os.makedirs(out,exist_ok=True)
rows=list(csv.DictReader(open("runs/push/final20_v3.csv")))
n=0
for r in rows:
    d=r["design_id"]
    cif=None
    for i in range(4):
        p=f"runs/boltz_rf_h{i}/out/by_seed/{d}/s0/{d}_model_0.cif"
        if os.path.exists(p): cif=p; break
    if not cif: print("MISS",d); continue
    st=gemmi.read_structure(cif)
    st.setup_entities()
    # keep first model; relabel chains A (binder, first) / B (target, second) by order
    m=st[0]
    chains=list(m)
    # Boltz writes binder=A, target=B already; normalize labels
    for idx,ch in enumerate(chains):
        ch.name = "A" if idx==0 else "B"
    st.write_pdb(f"{out}/{d}.pdb")
    n+=1
print(f"wrote {n} complex PDBs -> {out}")
