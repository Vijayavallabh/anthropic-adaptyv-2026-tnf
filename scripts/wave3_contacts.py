#!/usr/bin/env python3
"""WAVE-3 His->cation contact counter. For each scored variant, parse its Boltz HUMAN fold (wv3h*),
identify binder (shorter chain) + target (longer chain), and count how many of the variant's engineered
His (manifest 'his' positions) sit with ND1/NE2 within 6A of a target cation (Lys NZ / Arg NH1/NH2/NE) in
the folded pose. Adds his_cation_contacts + multi_his_dual_strong (dual_strong AND >=2 contacts) to
runs/wave3/scored.json. Run with .venv/bin/python -I."""
import json, os, glob
import numpy as np
from Bio.PDB import MMCIFParser
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
CUT=6.0
scored=json.load(open("runs/wave3/scored.json"))
mani=json.load(open("runs/wave3/manifest.json"))
parser=MMCIFParser(QUIET=True)
def cif(v):
    for t in ("wv3h0","wv3h1","wv3h2","wv3h3"):
        p=f"runs/boltz_{t}/out/by_seed/{v}/s0/{v}_model_0.cif"
        if os.path.exists(p): return p
    return None
def count(v, his):
    c=cif(v)
    if not c or not his: return None
    model=parser.get_structure("x",c)[0]
    chains=[(ch, sum(1 for _ in ch)) for ch in model]
    if len(chains)<2: return None
    chains.sort(key=lambda x:x[1])
    binder=chains[0][0]; target=chains[-1][0]
    catN=[]
    for res in target:
        for a in res:
            if (res.resname=="LYS" and a.name=="NZ") or (res.resname=="ARG" and a.name in("NH1","NH2","NE")):
                catN.append(a.coord)
    if not catN: return 0
    hitn=0
    bres={r.id[1]:r for r in binder}
    for hp in his:
        r=bres.get(hp)
        if r is None or r.resname not in ("HIS","HID","HIE","HIP"): continue
        ha=[a.coord for a in r if a.name in ("ND1","NE2")]
        if ha and min(np.linalg.norm(h-c) for h in ha for c in catN)<CUT: hitn+=1
    return hitn
nmulti=0
for o in scored:
    v=o["variant"]
    his=mani.get(v,{}).get("his",[])
    hc=count(v,his) if his else 0
    o["his_cation_contacts"]=hc
    o["n_his"]=len(his)
    o["multi_his_dual_strong"]=bool(o.get("dual_strong") and (hc or 0)>=2)
    if o["multi_his_dual_strong"]: nmulti+=1
json.dump(scored, open("runs/wave3/scored.json","w"), indent=1)
ms=[o for o in scored if o["multi_his_dual_strong"]]
print(f"\n[wave3 contacts] MULTI-HIS DUAL-STRONG (>=2 His-cation contacts, balanced, keeps binding): {len(ms)}")
print(f"{'variant':24}{'min':>7}{'dHM':>7}{'nHis':>5}{'contacts':>9}")
for o in sorted(ms,key=lambda x:(-x['his_cation_contacts'],-x['min'])):
    print(f"  {o['variant']:24}{o['min']:>7.3f}{o['dHM']:>7.3f}{o['n_his']:>5}{o['his_cation_contacts']:>9}")
