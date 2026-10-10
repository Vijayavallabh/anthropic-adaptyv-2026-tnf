#!/usr/bin/env python3
"""WAVE-6 (Campaign H): 3-His SEPARATED-cation flagship builder (full Baker SOTA pH-switch form).
For each target design's validated Chai human complex: pick up to 3 binder positions, each the CB-nearest
to a DISTINCT conserved target cation (within 8A), under the constraint that the chosen His positions are
mutually separated >=8A in binder space (anti-clustering -> cooperative Delta n_H, per
lesson-sharp-ph-switch-separated-his). His-lock them, define the local accommodation shell (binder residues
within 7A of any His CB, excluding the His + termini), emit a LigandMPNN job (core+His+target fixed, shell
redesigned, --omit_AA C). Writes runs/wave6/<name>.pdb + jobs.tsv.
Usage: third_party/LigandMPNN/.venv/bin/python -I scripts/wave6_prep.py <name1> <name2> ...
       (reads each design's Chai human complex; falls back to any chai complex)."""
import os, sys, copy, glob
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
OUT="runs/wave6"; os.makedirs(OUT, exist_ok=True)
P=MMCIFParser(QUIET=True); io=PDBIO()
CONS=[1,6,27,39,60,77,85,93,98,107,123]
# priority order for cations: strong anchors first (Lys6/Lys93), then the rest
PRIO=[93,6,60,77,1,107,98,27,39,85,123]
SEP=8.0   # minimum His-His separation in binder space (anti-clustering)

def find_chai_human(name):
    for d in sorted(glob.glob("runs/chai*_h")+glob.glob("runs/chai*_ch")):
        f=f"{d}/out/fast/{name}/seed_0/pred.model_idx_0.cif"
        if os.path.exists(f): return f
    g=glob.glob(f"runs/chai*/out/fast/{name}/seed_0/pred.model_idx_0.cif")
    return g[0] if g else None

def cb(res):
    v=[a.coord for a in res if a.name=="CB"] or [a.coord for a in res if a.name=="CA"]
    return v[0] if v else None

names=sys.argv[1:]
jobs=[]
print(f"{'design':34s} {'His-lock (pos->cat)':40s} shellN")
for n in names:
    c=find_chai_human(n)
    if not c: print(f"{n}: no chai cif"); continue
    m=P.get_structure("x",c)[0]
    chs=sorted(m,key=lambda ch:sum(1 for _ in ch)); binder,target=chs[0],chs[-1]
    st=Structure.Structure("x"); md=Model.Model(0); st.add(md)
    nb=Chain.Chain("A"); i=1; bmap={}
    for res in binder:
        if res.id[0]!=' ': continue
        r=copy.deepcopy(res); r.id=(' ',i,' '); nb.add(r); bmap[i]=r; i+=1
    nt=Chain.Chain("B"); j=1; cat={}
    for res in target:
        if res.id[0]!=' ': continue
        r=copy.deepcopy(res); r.id=(' ',j,' '); nt.add(r)
        if r.resname in("LYS","ARG") and j in CONS:
            chg=[a.coord for a in r if a.name in("NZ","NH1","NH2","NE","CZ")]
            if chg: cat[j]=chg
        j+=1
    md.add(nb); md.add(nt)
    # greedily assign His positions: for each cation in priority order, take CB-nearest binder pos <8A
    # that is (a) not already used and (b) separated >=SEP from all chosen His positions.
    hislock={}   # cation -> binder pos
    chosen_cb=[]
    for cid in PRIO:
        if cid not in cat or len(hislock)>=3: continue
        ranked=[]
        for p,r in bmap.items():
            if p in hislock.values(): continue
            b=cb(r)
            if b is None: continue
            d=min(np.linalg.norm(b-ch) for ch in cat[cid])
            ranked.append((d,p,b))
        ranked.sort()
        for d,p,b in ranked:
            if d>=8.0: break
            if all(np.linalg.norm(b-cbx)>=SEP for cbx in chosen_cb):
                hislock[cid]=p; chosen_cb.append(b); break
    if len(hislock)<2:
        print(f"{n}: <2 separable anchor His (skip)"); continue
    for cid,pos in hislock.items(): bmap[pos].resname="HIS"
    lockpos=set(hislock.values())
    lockcb=[cb(bmap[p]) for p in lockpos]
    shell=[]
    for p,r in bmap.items():
        if p in lockpos: continue
        if p<=2 or p>=max(bmap)-1: continue
        b=cb(r)
        if b is not None and min(np.linalg.norm(b-lc) for lc in lockcb)<7.0: shell.append(p)
    pdb=f"{OUT}/{n}.pdb"; io.set_structure(st); io.save(pdb)
    redes=" ".join(f"A{p}" for p in sorted(shell))
    jobs.append((n,os.path.abspath(pdb),redes))
    disp={c:hislock[c] for c in hislock}
    print(f"{n:34s} {str(disp):40s} {len(shell)}")
with open(f"{OUT}/jobs.tsv","w") as f:
    for s,p,red in jobs: f.write(f"{s}\t{p}\t{red}\n")
print(f"\nwrote {len(jobs)} jobs -> {OUT}/jobs.tsv")
