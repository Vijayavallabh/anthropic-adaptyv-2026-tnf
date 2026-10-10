#!/usr/bin/env python3
"""WAVE-7 (c19 partial-diffusion family): His-bias prep for diffused backbones.
Target B is fixed across all diffused backbones -> take conserved-cation charged-atom coords ONCE from the
seed full-atom target (boltzchk). For each diffused backbone (binder A backbone-only, target B kept), bias
His at binder CA positions within CUT of an interface conserved cation. Emits staged pdbs + bias json.
third_party/LigandMPNN/.venv/bin/python -I scripts/wave7_prep.py"""
import glob, json, os, copy
import numpy as np
from Bio.PDB import PDBParser, PDBIO
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
SEED="runs/hhead/triple_pdb/des_3050_6_c19_boltzchk.pdb"
OUT="runs/wave7"; os.makedirs(f"{OUT}/staged", exist_ok=True)
HIS_BIAS=2.0; CUT=9.0; MIN_SITES=2; MAX_SITES=5
CONS=[1,6,27,39,60,77,85,93,98,107,123]
P=PDBParser(QUIET=True); io=PDBIO()
# seed cations (conserved positions) + seed binder for interface filter
sm=P.get_structure("s",SEED)[0]
tb=sm["B"]; ba=sm["A"]
binder_ca=np.array([r["CA"].coord for r in ba if "CA" in r])
cats=[]  # (bpos-independent) charged-atom coords of CONSERVED interface cations
for r in tb:
    if r.id[1] in CONS and r.resname in("LYS","ARG"):
        chg=[a.coord for a in r if a.name in("NZ","NH1","NH2","NE")]
        for c in chg:
            if binder_ca.size and np.linalg.norm(binder_ca-c,axis=1).min()<13.0:
                cats.append(c)
cats=np.array(cats)
print(f"seed interface conserved-cation charged atoms: {len(cats)}")
bias={}; staged=[]; kept=0
for pdb in sorted(glob.glob("runs/rfd/c19pd_T*/des_*.pdb")):
    if "traj" in pdb: continue
    m=P.get_structure("d",pdb)[0]
    if "A" not in m: continue
    sites=[]
    for r in m["A"]:
        ca=r["CA"].coord if "CA" in r else None
        if ca is None: continue
        d=float(np.linalg.norm(cats-ca,axis=1).min()) if cats.size else 99
        if d<=CUT: sites.append((r.id[1],d))
    if len(sites)<MIN_SITES: continue
    sites.sort(key=lambda x:x[1]); sites=sites[:MAX_SITES]
    tag=pdb.replace("runs/rfd/","").replace("/","_").replace(".pdb","")
    sp=f"{OUT}/staged/{tag}.pdb"; os.symlink(os.path.abspath(pdb), sp) if not os.path.exists(sp) else None
    bias[os.path.abspath(sp)]={f"A{p}":{"H":HIS_BIAS} for p,_ in sites}
    staged.append(os.path.abspath(sp)); kept+=1
json.dump(bias, open(f"{OUT}/bias_AA_per_residue_multi.json","w"))
json.dump(staged, open(f"{OUT}/pdbs.json","w"))
print(f"kept {kept} backbones (>= {MIN_SITES} His-biasable interface sites) -> {OUT}/staged + bias json")
