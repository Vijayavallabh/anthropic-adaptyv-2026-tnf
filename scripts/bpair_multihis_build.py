#!/usr/bin/env python3
"""WAVE-10: MULTI-His combinatorial switch scan on the NOVEL cross-reactive bpair de-novo folds.
These folds were only ever SINGLE-His scanned (session 9b) or never scanned; a multi-His separated-cation
switch (the Baker SOTA form, like the c19 flagship's 4 His) has never been tried on them. For each design:
read the af2pred bound pose, find binder (chain A) positions whose CB is < DMAX of a target (chain B)
cation, keep one anchor per distinct cation (closest), then enumerate His COMBINATIONS that are mutually
SEPARATED >= SEP Angstrom in binder space (anti-cooperativity -> cooperative Delta n_H) and target DISTINCT
cations. Emit wt + singles + separated pairs + separated triples (+ quad if available), capped per design.
Downstream: Boltz human+mouse fold each variant -> aggregate Wyman linkage + cross-react iptm -> the best
go to Chai + PROPKA+xTB both-fold triple-confirm. Conservation of the repelled cation is validated by the
mouse co-fold (a variant that keeps binding mouse repels a cross-reactive cation).
Usage: third_party/LigandMPNN/.venv/bin/python scripts/bpair_multihis_build.py [design1 design2 ...]
"""
import os, glob, csv, sys, itertools
import numpy as np
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
DMAX=8.0     # binder CB -> target cation: an installable repulsion anchor
SEP=8.0      # min His-His separation in binder space (anti-clustering -> cooperative switch)
MAXVAR=16    # cap variants per design (keep Boltz tractable)
TCHAINS=set(os.environ.get("MH_TARGET_CHAINS","B"))  # target chains to scan for cations (two-protomer groove: "BC")
OUT=os.environ.get("MHOUT","runs/wave10"); os.makedirs(OUT, exist_ok=True)

DEFAULT=["des_962_5","des_84_3","des_2140_5","des_287_5","des_2106_2","des_962_6","des_495_1","des_495_2",
         "des_962_2","des_495_3","des_737_2","des_1447_6","des_2375_4","des_935_4"]
designs=sys.argv[1:] or DEFAULT
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open(os.environ.get("MH_SEQS","runs/bpair/xreact_full.csv")))}

def pose_for(d):
    c=glob.glob(f"{R}/"+os.environ.get("MH_POSEGLOB","runs/bpair/af2ig_full/gpu*/pdbs")+f"/{d}_af2pred.pdb")
    return c[0] if c else None

rows=[("name","sequence")]; manifest=[]
print(f"{'design':12s} {'L':>4} anchors(pos->cat@dist)                          nvar")
for d in designs:
    s=seqs.get(d); pose=pose_for(d)
    if not s or not pose:
        print(f"{d}: missing seq/pose, skip"); continue
    bCB={}; tcat={}
    for l in open(pose):
        if l[:4]!="ATOM": continue
        cc=l[21]; res=l[17:20].strip(); at=l[12:16].strip(); rn=int(l[22:26])
        crd=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        if cc=="A" and at=="CB": bCB[rn]=(crd,res)
        if cc in TCHAINS and ((res=="LYS" and at=="NZ") or (res=="ARG" and at in("CZ","NH1","NH2"))):
            tcat.setdefault((cc,rn),[]).append(crd)   # key by (chain,resnum) so B/C cations don't collide
    if not tcat: print(f"{d}: no target cations, skip"); continue
    # for each binder position, nearest target cation + dist
    pos_cat={}
    for rn,(crd,res) in bCB.items():
        if res in("GLY","PRO","HIS") or rn>len(s) or rn<=1: continue
        best=min(((min(np.linalg.norm(crd-c) for c in cs),cat) for cat,cs in tcat.items()))
        if best[0]<=DMAX: pos_cat[rn]=(best[1],round(best[0],1),crd)   # pos -> (cation, dist, cb)
    if not pos_cat: print(f"{d}: no interface anchors <{DMAX}A, skip"); continue
    # one anchor per distinct cation (closest binder pos), sorted by dist
    by_cat={}
    for pos,(cat,dist,cb) in pos_cat.items():
        if cat not in by_cat or dist<by_cat[cat][1]: by_cat[cat]=(pos,dist,cb)
    anchors=sorted(by_cat.values(), key=lambda x:x[1])   # (pos,dist,cb)
    # enumerate separated combos (>=SEP apart in binder CB space), distinct cations guaranteed by by_cat
    def sep_ok(combo):
        for a,b in itertools.combinations(combo,2):
            if np.linalg.norm(a[2]-b[2])<SEP: return False
        return True
    variants=[]  # (hist_positions, ncat)
    for pos,dist,cb in anchors: variants.append(([pos],1))                       # singles
    for combo in itertools.combinations(anchors,2):
        if sep_ok(combo): variants.append(([c[0] for c in combo],2))             # separated pairs
    for combo in itertools.combinations(anchors,3):
        if sep_ok(combo): variants.append(([c[0] for c in combo],3))             # separated triples
    for combo in itertools.combinations(anchors,4):
        if sep_ok(combo): variants.append(([c[0] for c in combo],4))             # separated quads
    # prioritize higher cardinality (more His = sharper), then tighter mean dist
    dmap={p:d for p,d,_ in anchors}
    variants.sort(key=lambda v:(-v[1], np.mean([dmap[p] for p in v[0]])))
    variants=variants[:MAXVAR]
    rows.append((f"{d}_wt", s))
    for hpos,ncat in variants:
        mut=list(s)
        for p in hpos: mut[p-1]="H"
        name=f"{d}_mh_"+"_".join(map(str,sorted(hpos)))
        rows.append((name,"".join(mut)))
        manifest.append((name,d,sorted(hpos),ncat))
    adesc=" ".join(f"{p}->{by_cat[[c for c in by_cat if by_cat[c][0]==p][0]][0] if False else ''}" for p,_,_ in anchors[:6])
    print(f"{d:12s} {len(s):>4} "+",".join(f"{p}@{dd}" for p,dd,_ in anchors[:7])+f"   {len(variants)}")

# write combined + 8 GPU shards
with open(f"{OUT}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); [w.writerow(r) for r in rows]
data=rows[1:]
NSH=8
for k in range(NSH):
    with open(f"{OUT}/shard{k}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"])
        for r in data[k::NSH]: w.writerow(r)
with open(f"{OUT}/manifest.tsv","w") as fh:
    for name,d,hpos,ncat in manifest: fh.write(f"{name}\t{d}\t{','.join(map(str,hpos))}\t{ncat}\n")
print(f"\nwrote {len(data)} variants ({sum(1 for m in manifest if m[3]>=2)} multi-His) -> {OUT}/variants.csv + {NSH} shards")
