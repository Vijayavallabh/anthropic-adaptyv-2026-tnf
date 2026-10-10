#!/usr/bin/env python3
"""WAVE-3 multi-His builder (Baker 2025 SOTA recipe: >=2-3 SEPARATED His each repelling a DISTINCT
conserved TARGET cation). From the top BALANCED already-binding groove binders (runs/bpair/xreact_full2_results.txt),
install His COMBINATIONS (singles + 2-combos + 3-combos) at binder positions whose CB is <=8A of a
human-mouse CONSERVED target cation (Lys/Arg) in the af2pred bound pose. Conserved cations derived from a
human<->mouse mature-TNF alignment (NOT hardcoded numbers). Emits variants.csv + 4 GPU shards + a manifest
(variant -> base, his_positions, cations) under runs/wave3/. Run with .venv/bin/python -I."""
import os, glob, csv, itertools, json
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
OUT="runs/wave3"; os.makedirs(OUT, exist_ok=True)
NTOP=14          # balanced binders to redesign
MAXCAT=4         # distinct conserved cations per binder (=> up to C(4,2)+C(4,3)=10 multi-His)
CONTACT=8.0      # CB->cation-N contact cutoff (design pose)
CHB_OFFSET=80    # af2pred chain B residue number = human-seq index + 80 (verified: B81=human R1)

# ---- 1. conserved cations from a human<->mouse mature-TNF global alignment (human 1-based index) ----
h=open("runs/buried2/human_tnf.txt").read().strip()
m=open("runs/buried2/mouse_tnf_mature.txt").read().strip()
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load("BLOSUM62")
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode="global"
aln=al.align(h,m)[0]
hm={}  # human_index(1-based) -> aligned mouse char
for (hs,he),(ms,me) in zip(aln.aligned[0],aln.aligned[1]):
    for k in range(he-hs): hm[hs+k+1]=m[ms+k]
conserved_cat={i for i,aa in enumerate(h,1) if aa in "KR" and hm.get(i) in ("K","R")}
print(f"[wave3] human {len(h)} aa, mouse {len(m)} aa; conserved cations (human idx): "
      +",".join(str(i) for i in sorted(conserved_cat)))

# ---- 2. top BALANCED already-binding binders ----
cand=[]
for l in open("runs/bpair/xreact_full2_results.txt"):
    p=l.split()
    if len(p)!=4 or p[0]=="design": continue
    d,hh,mm,mn=p[0],float(p[1]),float(p[2]),float(p[3])
    cand.append((d,hh,mm,mn,abs(hh-mm)))
cand=[c for c in cand if c[3]>=0.55]
cand.sort(key=lambda c:-(c[3]-0.5*max(0.0,c[4]-0.10)))   # balanced score
sel=cand[:NTOP]
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/bpair/xreact_full2.csv"))}
def af2pred(d):
    c=glob.glob(f"{R}/runs/bpair/af2ig_full2/gpu*/pdbs/{d}_af2pred.pdb")
    return c[0] if c else None

rows=[]; mani={}
print(f"\n{'binder':14}{'min':>6}{'dHM':>6}  distinct-cations(binderPos->humanCat@dist)")
for d,hh,mm,mn,dl in sel:
    s=seqs.get(d); f=af2pred(d)
    if not s or not f: print(f"  {d}: missing seq/af2pred"); continue
    aCB={}; catN={}
    for l in open(f):
        if l[:4]!="ATOM": continue
        ch=l[21]; res=l[17:20].strip(); at=l[12:16].strip(); rn=int(l[22:26])
        crd=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        if ch=="A" and at=="CB": aCB[rn]=crd
        if ch=="B":
            hidx=rn-CHB_OFFSET
            if hidx in conserved_cat and ((res=="LYS" and at=="NZ") or (res=="ARG" and at in("NH1","NH2","NE"))):
                catN.setdefault(hidx,[]).append(crd)
    pos2cat={}
    for rn,crd in aCB.items():
        if rn>len(s) or s[rn-1] in "GPH": continue
        best=None
        for hidx,coords in catN.items():
            dmn=min(float(np.linalg.norm(crd-c)) for c in coords)
            if dmn<CONTACT and (best is None or dmn<best[1]): best=(hidx,dmn)
        if best: pos2cat[rn]=best
    cat2pos={}
    for rn,(hidx,dist) in pos2cat.items():
        if hidx not in cat2pos or dist<cat2pos[hidx][1]: cat2pos[hidx]=(rn,dist)
    chosen=sorted(cat2pos.items(), key=lambda kv:kv[1][1])[:MAXCAT]   # closest distinct cations
    cations=[(hidx,rn,round(dist,2)) for hidx,(rn,dist) in chosen]
    posn=[rn for hidx,rn,dist in cations]
    cat_of={rn:hidx for hidx,rn,dist in cations}
    print(f"  {d:12}{mn:>6.3f}{dl:>6.3f}  "+", ".join(f"{rn}->K/R{hidx}@{dist}" for hidx,rn,dist in cations))
    if not posn: continue
    def mk(poslist):
        ss=list(s)
        for p in poslist: ss[p-1]="H"
        return "".join(ss)
    def add(tag,poslist):
        name=f"{d}_H{'-'.join(map(str,poslist))}" if poslist else f"{d}_wt"
        rows.append((name, mk(poslist)))
        mani[name]=dict(base=d, his=list(poslist), cations=[cat_of[p] for p in poslist], min0=mn, dHM0=dl)
    add("wt",[])
    for p in posn: add("s",[p])
    for c in itertools.combinations(sorted(posn),2): add("d",list(c))
    for c in itertools.combinations(sorted(posn),3): add("t",list(c))

# dedup by sequence-name
seen=set(); urows=[]
for n,s in rows:
    if n in seen: continue
    seen.add(n); urows.append((n,s))
with open(f"{OUT}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in urows]
shards=[[] for _ in range(4)]
for i,r in enumerate(urows): shards[i%4].append(r)
for k,sh in enumerate(shards):
    with open(f"{OUT}/shard{k}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in sh]
json.dump(mani, open(f"{OUT}/manifest.json","w"), indent=1)
nmulti=sum(1 for n,_ in urows if n.count("-")>=1)
print(f"\n[wave3] {len(urows)} variants from {len(sel)} binders ({nmulti} multi-His) -> {OUT}/variants.csv + 4 shards + manifest.json")
