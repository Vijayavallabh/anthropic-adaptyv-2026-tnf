#!/usr/bin/env python3
"""WAVE-2 His-switch scan: select the top BALANCED cross-react de-novo binders from wave-2
(runs/bpair/xreact_full2_results.txt), then install His at each binder interface position whose CB is
<=8A of a target cation in the af2pred bound pose. Prioritizes BALANCED two-species binding (small Δ_hm)
per the stress-test's dual-strong gate. Emits variants + 4 GPU shards under runs/bpair_mature2/."""
import os, glob, csv
import numpy as np
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
NTOP=16; NCLOSE=8; MIN_IPSAE=0.52
OUT="runs/bpair_mature2"; os.makedirs(OUT, exist_ok=True)

# parse wave-2 cross-react results (design h m min)
cand=[]
for l in open("runs/bpair/xreact_full2_results.txt"):
    p=l.split()
    if len(p)!=4 or p[0]=="design": continue
    d,h,m,mn=p[0],float(p[1]),float(p[2]),float(p[3])
    cand.append((d,h,m,mn,abs(h-m)))
# balanced score: reward min-ipSAE, penalize human-mouse gap (Δ_hm) above 0.10 (folder-optimism / mouse weakness)
cand=[c for c in cand if c[3]>=MIN_IPSAE]
cand.sort(key=lambda c: -(c[3] - 0.5*max(0.0, c[4]-0.10)))
sel=cand[:NTOP]
print(f"selected {len(sel)}/{len(cand)} balanced binders (min>={MIN_IPSAE}):")
for d,h,m,mn,dl in sel: print(f"  {d:14} min={mn:.3f} dHM={dl:.3f} (h{h:.3f}/m{m:.3f})")

seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/bpair/xreact_full2.csv"))}
def af2pred(d):
    c=glob.glob(f"{R}/runs/bpair/af2ig_full2/gpu*/pdbs/{d}_af2pred.pdb")+glob.glob(f"{R}/runs/bpair/af2ig_full2/**/{d}_af2pred.pdb",recursive=True)
    return c[0] if c else None
rows=[]
for d,h,m,mn,dl in sel:
    s=seqs.get(d); f=af2pred(d)
    if not s or not f: print(f"  {d}: missing seq/af2pred, skip"); continue
    bCB={}; tcat=[]
    for l in open(f):
        if l[:4]!="ATOM": continue
        cc=l[21]; res=l[17:20].strip(); at=l[12:16].strip(); rn=int(l[22:26])
        crd=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        if cc=="A" and at=="CB": bCB[rn]=(crd,res)
        if cc=="B" and ((res=="LYS" and at=="NZ") or (res=="ARG" and at in("CZ","NH1","NH2"))): tcat.append(crd)
    if not tcat: continue
    dmin={rn:min(np.linalg.norm(crd-tc) for tc in tcat) for rn,(crd,res) in bCB.items()
          if res not in("GLY","PRO","HIS") and rn<=len(s)}
    pos=sorted([p for p,_ in sorted(dmin.items(),key=lambda kv:kv[1])[:NCLOSE] if dmin[p]<8.0])
    if not pos: continue
    rows.append((f"{d}_wt", s))
    for p in pos: rows.append((f"{d}_H{p}", s[:p-1]+"H"+s[p:]))

with open(f"{OUT}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in rows]
shards=[[] for _ in range(4)]
for i,r in enumerate(rows): shards[i%4].append(r)
for k,sh in enumerate(shards):
    with open(f"{OUT}/shard{k}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in sh]
print(f"\nTOTAL {len(rows)} variants ({len(sel)} binders) -> {OUT}/variants.csv + 4 shards")
