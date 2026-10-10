#!/usr/bin/env python3
"""WAVE-5 (F2) prep: extend the proven wave3b/4b multi-His method to DIVERSE backbones (new lineages).
For each diverse v7-floor backbone that contacts BOTH strong anchors (Lys93 + Lys6) in its Chai complex,
His-lock the binder positions nearest Lys93 and Lys6 (+ a 3rd contacted conserved cation if a clean position
exists), then let LigandMPNN redesign the interface shell AROUND the fixed His (core + His + target fixed) so
the fold packs to accommodate the multi-His — exactly what gave the robust c44/b14 on des_4098_1. Writes
runs/wave5/<scaffold>.pdb + jobs.tsv. .venv/bin/python -I."""
import os, copy, csv
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
OUT="runs/wave5"; os.makedirs(OUT, exist_ok=True)
P=MMCIFParser(QUIET=True); io=PDBIO()
CONS=[1,6,27,39,60,77,85,93,98,107,123]; ANCHORS=[93,6]
# candidate diverse backbones (contact both anchors) -> their Chai human complex
CAND=["des_2106_2_H43","des_287_5_H37","des_59_1_H51","b2m0_c70","b1m7_2A57H_K183_r1","TNFpH06",
      "b1m12_A41H_vsR27","b4m2_r1","b1m9_A12H","b1m12_2A67H_K112_r1","b4m2_A57H"]
def cif(n): return f"runs/chai_chai_v7_h/out/fast/{n}/seed_0/pred.model_idx_0.cif"
jobs=[]
print(f"{'scaffold':24s} {'His-lock (pos->cat)':34s} shellN")
for n in CAND:
    c=cif(n)
    if not os.path.exists(c): print(f"{n}: no chai cif"); continue
    m=P.get_structure("x",c)[0]; chs=sorted(m,key=lambda ch:sum(1 for _ in ch)); binder,target=chs[0],chs[-1]
    # renumber binder A 1..L, target B 1..N
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
    # target atoms for shell
    tatoms=[a.coord for res in nt for a in res]
    # His-lock: binder residue (by CB) nearest each anchor cation, within 8A
    hislock={}
    for cid in ANCHORS+[x for x in CONS if x not in ANCHORS]:
        if cid not in cat: continue
        best=None
        for p,r in bmap.items():
            if p in hislock.values(): continue
            cb=[a.coord for a in r if a.name=="CB"] or [a.coord for a in r if a.name=="CA"]
            if not cb: continue
            d=min(np.linalg.norm(cb[0]-ch) for ch in cat[cid])
            if best is None or d<best[1]: best=(p,d)
        if best and best[1]<8.0 and best[0] not in hislock.values():
            hislock[cid]=best[0]
        if len(hislock)>=2: break          # 2 His at the STRONG anchors (Lys93+Lys6), like c44/b14
    if len(hislock)<2: print(f"{n}: <2 anchor His positions (skip)"); continue
    # set His at lock positions (resname only; fold builds rotamer)
    for cid,pos in hislock.items():
        res=bmap[pos]
        res.resname="HIS"
    # shell = LOCAL accommodation around the His-lock (within 7A of a His CB), excluding His-lock + termini
    lockpos=set(hislock.values())
    lockcb=[([a.coord for a in bmap[pp] if a.name=="CB"] or [a.coord for a in bmap[pp] if a.name=="CA"])[0] for pp in lockpos]
    shell=[]
    for p,r in bmap.items():
        if p in lockpos: continue
        if p<=2 or p>=max(bmap)-1: continue
        cb=[a.coord for a in r if a.name in("CB","CA")]
        if cb and min(np.linalg.norm(cb[0]-lc) for lc in lockcb)<7.0: shell.append(p)
    pdb=f"{OUT}/{n}.pdb"; io.set_structure(st); io.save(pdb)
    redes=" ".join(f"A{p}" for p in sorted(shell))
    jobs.append((n,os.path.abspath(pdb),redes))
    print(f"{n:24s} {str({c:hislock[c] for c in hislock}):34s} {len(shell)}")
with open(f"{OUT}/jobs.tsv","w") as f:
    for s,p,red in jobs: f.write(f"{s}\t{p}\t{red}\n")
print(f"\nwrote {len(jobs)} jobs -> {OUT}/jobs.tsv")
