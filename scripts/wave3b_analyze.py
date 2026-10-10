#!/usr/bin/env python3
"""WAVE-3B analyze: for each MPNN interface-adapted design, derive His positions FROM THE SEQUENCE (all H),
measure each His->nearest distinct target cation in the Boltz human fold, compute dual-strong (vs scaffold
_wt baseline) + switch (>=2 distinct conserved cations <8A). Writes runs/wave3b/switch_ranked.json.
Run with .venv/bin/python -I."""
import json, os
import numpy as np
from Bio.PDB import MMCIFParser
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
import csv
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/wave3b/variants.csv"))}
scored={o["variant"]:o for o in json.load(open("runs/wave3b/scored.json"))}
P=MMCIFParser(QUIET=True)
def cif(v):
    for k in range(4):
        p=f"runs/boltz_wv3bh{k}/out/by_seed/{v}/s0/{v}_model_0.cif"
        if os.path.exists(p): return p
    return None
def his_cat(v,his):
    c=cif(v)
    if not c or not his: return None
    model=P.get_structure("x",c)[0]
    chs=sorted(model,key=lambda ch:sum(1 for _ in ch)); binder,target=chs[0],chs[-1]
    cats=[]
    for res in target:
        for a in res:
            if (res.resname=="LYS" and a.name=="NZ") or (res.resname=="ARG" and a.name in("NH1","NH2","NE")):
                cats.append((res.id[1],a.coord))
    bres={r.id[1]:r for r in binder}; perh={}
    for hp in his:
        r=bres.get(hp)
        if r is None or r.resname not in("HIS","HID","HIE","HIP"): continue
        ha=[a.coord for a in r if a.name in("ND1","NE2")]
        if ha and cats:
            d,rid=min((float(np.linalg.norm(h-cc)),rid) for h in ha for rid,cc in cats); perh[hp]=(round(d,2),rid)
    return perh
def base_of(n):
    if n.endswith("_wt"): return n[:-3]
    if "_b" in n: return n.rsplit("_b",1)[0]
    return n
wt={}
for v,o in scored.items():
    if v.endswith("_wt"): wt[v[:-3]]=o["min"]
out=[]
for v,o in scored.items():
    if v.endswith("_wt"): continue
    s=seqs.get(v,""); his=[i+1 for i,ch in enumerate(s) if ch=="H"]
    perh=his_cat(v,his) or {}
    distinct=set(rid for d,rid in perh.values() if d<8.0)
    base=base_of(v); wtmin=wt.get(base)
    keeps=(wtmin is None) or (o["min"]-wtmin>=-0.08)
    dual=o["min"]>=0.55 and o["dHM"]<=0.13 and keeps
    sw=dual and len(distinct)>=2
    out.append(dict(variant=v,min=o["min"],dHM=o["dHM"],h=o["h"],m=o["m"],wt_min=wtmin,
                    nhis=len(his),his=his,ndistinct_cat_h=len(distinct),
                    d_h={k:perh[k][0] for k in perh},dual_strong=bool(dual),switch=bool(sw)))
out.sort(key=lambda z:(-int(z["switch"]),-z["ndistinct_cat_h"],-z["min"]))
json.dump(out, open("runs/wave3b/switch_ranked.json","w"), indent=1)
nsw=sum(1 for z in out if z["switch"])
print(f"{'variant':24}{'min':>7}{'dHM':>6}{'nHis':>5}{'distinctCat':>12}  dists")
for z in out[:40]:
    tag="  SWITCH" if z["switch"] else ("  dual" if z["dual_strong"] else "")
    print(f"{z['variant']:24}{z['min']:>7.3f}{z['dHM']:>6.3f}{z['nhis']:>5}{z['ndistinct_cat_h']:>12}  {z['d_h']}{tag}")
print(f"\n[wave3b] MULTI-HIS SWITCH (interface-adapted, dual-strong + >=2 distinct cations <8A): {nsw}")
