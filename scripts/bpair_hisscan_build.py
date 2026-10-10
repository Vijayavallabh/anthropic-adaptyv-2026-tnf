#!/usr/bin/env python3
"""Build His-switch variants for the top beta-pairing de-novo cross-react winners.
For each design, read its af2pred bound pose, find binder (chain A) interface positions whose CB is
closest to a target (chain B) cation (Lys NZ / Arg CZ/NH), install His at the N closest such positions,
and emit variant sequences (plus wt). Writes a combined CSV and 4 GPU shards under runs/bpair_mature/.
The His-scan then folds each variant human+mouse (Boltz) and scores min-ipSAE to keep cross-react-
preserving switch sites; MD/PROPKA/xtb triple-confirm follows on the best."""
import os, glob, csv
import numpy as np
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
DESIGNS=["des_287_5","des_2140_5","des_962_5","des_84_3","des_2106_2","des_495_2","des_962_6","des_495_1"]
NCLOSE=8   # His at the 8 interface positions closest to a target cation per design
OUT="runs/bpair_mature"; os.makedirs(OUT, exist_ok=True)

seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/bpair/xreact_full.csv"))}
rows=[]
for d in DESIGNS:
    s=seqs.get(d)
    if not s:
        print(f"{d}: NO SEQ, skip"); continue
    cands=glob.glob(f"{R}/runs/bpair/af2ig_full/gpu*/pdbs/{d}_af2pred.pdb")
    if not cands:
        print(f"{d}: NO af2pred, skip"); continue
    f=cands[0]
    bCB={}; tcat=[]
    for l in open(f):
        if l[:4]!="ATOM": continue
        cc=l[21]; res=l[17:20].strip(); at=l[12:16].strip(); rn=int(l[22:26])
        crd=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        if cc=="A" and at=="CB": bCB[rn]=(crd,res)
        if cc=="B" and ((res=="LYS" and at=="NZ") or (res=="ARG" and at in("CZ","NH1","NH2"))): tcat.append(crd)
    if not tcat:
        print(f"{d}: no target cations, skip"); continue
    dmin={}
    for rn,(crd,res) in bCB.items():
        if res in ("GLY","PRO","HIS"): continue
        if rn>len(s): continue
        dmin[rn]=min(np.linalg.norm(crd-tc) for tc in tcat)
    pos=[p for p,_ in sorted(dmin.items(), key=lambda kv: kv[1])[:NCLOSE]]
    pos.sort()
    rows.append((f"{d}_wt", s))
    for p in pos:
        rows.append((f"{d}_H{p}", s[:p-1]+"H"+s[p:]))
    print(f"{d}: {len(s)}mer  His@{pos}  (closest-cation dists {[round(dmin[p],1) for p in pos]})")

# combined + 4 shards
with open(f"{OUT}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in rows]
shards=[[] for _ in range(4)]
for i,r in enumerate(rows): shards[i%4].append(r)
for k,sh in enumerate(shards):
    with open(f"{OUT}/shard{k}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in sh]
print(f"\nTOTAL {len(rows)} variants ({len(DESIGNS)} designs) -> {OUT}/variants.csv + 4 shards (~{len(shards[0])} each)")
