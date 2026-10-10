#!/usr/bin/env python3
"""WAVE-3 contact DIAGNOSTIC + recount. For every dual-strong variant, measure each engineered His's min
distance (ND1/NE2 -> nearest target cation N: Lys NZ / Arg NH1/NH2/NE) in BOTH the human and mouse Boltz
folds, and count contacts at 6/7/8A. Re-flag multi_his_switch = dual_strong AND >=2 His within 8A of a
distinct target cation in the HUMAN fold (electrostatic-repulsion contact, Baker-style). Writes
runs/wave3/switch_ranked.json. .venv/bin/python -I."""
import json, os, glob
import numpy as np
from Bio.PDB import MMCIFParser
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
scored=json.load(open("runs/wave3/scored.json")); mani=json.load(open("runs/wave3/manifest.json"))
P=MMCIFParser(QUIET=True)
def cif(v,sp):
    for k in range(4):
        p=f"runs/boltz_wv3{sp}{k}/out/by_seed/{v}/s0/{v}_model_0.cif"
        if os.path.exists(p): return p
    return None
def analyze(v,his,sp):
    c=cif(v,sp)
    if not c or not his: return None
    model=P.get_structure("x",c)[0]
    chs=sorted(model,key=lambda ch:sum(1 for _ in ch))
    binder,target=chs[0],chs[-1]
    cats=[]   # (resid, coord) sidechain cation N
    for res in target:
        for a in res:
            if (res.resname=="LYS" and a.name=="NZ") or (res.resname=="ARG" and a.name in("NH1","NH2","NE")):
                cats.append((res.id[1],a.coord))
    bres={r.id[1]:r for r in binder}
    perhis={}
    for hp in his:
        r=bres.get(hp)
        if r is None or r.resname not in("HIS","HID","HIE","HIP"): perhis[hp]=None; continue
        ha=[a.coord for a in r if a.name in("ND1","NE2")]
        if not ha or not cats: perhis[hp]=None; continue
        best=min((float(np.linalg.norm(h-cc)),rid) for h in ha for rid,cc in cats)
        perhis[hp]=best   # (dist, cation_resid)
    return perhis
out=[]
for o in scored:
    if not o.get("dual_strong"): continue
    v=o["variant"]; his=mani.get(v,{}).get("his",[])
    ph=analyze(v,his,"h"); pm=analyze(v,his,"m")
    if ph is None: continue
    dists=[d for d in (x[0] if x else None for x in ph.values()) if d is not None]
    cats_h=set(x[1] for x in ph.values() if x and x[0]<8.0)      # distinct cations contacted <8A (human)
    c6=sum(1 for d in dists if d<6); c7=sum(1 for d in dists if d<7); c8=sum(1 for d in dists if d<8)
    o2=dict(variant=v, nhis=len(his), min=o["min"], dHM=o["dHM"], wt_min=o["wt_min"],
            his=his, d_h={k:(round(x[0],2) if x else None) for k,x in ph.items()},
            ndistinct_cat_h=len(cats_h), c6=c6,c7=c7,c8=c8,
            switch=bool(len(cats_h)>=2))
    out.append(o2)
out.sort(key=lambda z:(-z["ndistinct_cat_h"], -z["min"]))
json.dump(out, open("runs/wave3/switch_ranked.json","w"), indent=1)
print(f"{'variant':22}{'nHis':>5}{'min':>7}{'dHM':>6}{'<6':>4}{'<7':>4}{'<8':>4}{'distinctCat':>12}  His->cat dists(human)")
for z in out:
    sw="  SWITCH" if z["switch"] else ""
    print(f"{z['variant']:22}{z['nhis']:>5}{z['min']:>7.3f}{z['dHM']:>6.3f}{z['c6']:>4}{z['c7']:>4}{z['c8']:>4}{z['ndistinct_cat_h']:>12}  {z['d_h']}{sw}")
nsw=sum(1 for z in out if z["switch"])
print(f"\nMULTI-HIS SWITCH (dual-strong AND >=2 distinct target cations <8A in human fold): {nsw}")
