#!/usr/bin/env python3
"""REBUILD: re-mine the pH-switch on backbones that already BIND + cross-react + express + pass novelty,
but whose CURRENT His placement gives a weak or fold-INCONSISTENT switch (e.g. des_818_4 reverses on Chai).
Generalizes bpair_multihis_build to ANY Boltz/Chai complex CIF (here the fresh v12 human complexes). For
each target: read the bound pose, find binder (chain A) positions whose CB is < DMAX of a target (chain B)
cation, keep one anchor per distinct cation, enumerate His COMBINATIONS mutually SEPARATED >= SEP (anti-
cooperativity), and emit wt + separated pairs/triples/quads as NEW sequences. Downstream: Boltz human+mouse
screen -> both-fold Wyman linkage (keep designs NEGATIVE on BOTH folds) + cross-react iptm -> MD + triple-confirm.
Usage: envs/uplifting-esmc/bin/python scripts/rescue_multihis_build.py <complex_glob> <seqs_csv> <out> [names...]
  complex_glob: a python-glob with {name} for the bound CIF, e.g. 'runs/boltz_v12H/out/by_seed/{name}/s0/{name}_model_0.cif'
"""
import os, sys, glob, csv, itertools
import numpy as np
from Bio.PDB import MMCIFParser
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
DMAX=8.0; SEP=8.0; MAXVAR=14
CONS={1,6,27,39,60,77,85,93,98,107,123}   # mature-TNF conserved cations (prefer these for cross-react)

def load(cif):
    m=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
    chs=sorted(m, key=lambda c: sum(1 for r in c if r.id[0]==' '))
    binder,target=chs[0],chs[-1]
    bCB={}
    for r in binder:
        if r.id[0]!=' ': continue
        cb=[a.coord for a in r if a.name=="CB"] or [a.coord for a in r if a.name=="CA"]
        if cb: bCB[r.id[1]]=(cb[0], r.resname)
    tcat={}
    # map target auth_seqid -> mature index via min offset to CONS is unsafe; just record all cations + resid
    for r in target:
        if r.id[0]!=' ': continue
        if r.resname in ("LYS","ARG"):
            chg=[a.coord for a in r if a.name in ("NZ","CZ","NH1","NH2","NE")]
            if chg: tcat[r.id[1]]=chg
    return bCB, tcat

def main():
    cglob, seqs_csv, out = sys.argv[1], sys.argv[2], sys.argv[3]
    names = sys.argv[4:]
    os.makedirs(f"{R}/{out}", exist_ok=True)
    seqs={r['name']:r['sequence'].split(':')[0] for r in csv.DictReader(open(f"{R}/{seqs_csv}"))}
    rows=[("name","sequence")]; manifest=[]
    for n in names:
        cif=glob.glob(f"{R}/"+cglob.format(name=n))
        s=seqs.get(n)
        if not cif or not s: print(f"{n}: missing cif/seq, skip"); continue
        bCB,tcat=load(cif[0])
        if not tcat: print(f"{n}: no target cations"); continue
        pos_cat={}
        for p,(cb,res) in bCB.items():
            if res in ("GLY","PRO") or p>len(s) or p<=1: continue
            best=min(((min(np.linalg.norm(cb-c) for c in cs), cat) for cat,cs in tcat.items()))
            if best[0]<=DMAX: pos_cat[p]=(best[1],best[0],cb)
        if not pos_cat: print(f"{n}: no anchors"); continue
        by_cat={}
        for p,(cat,d,cb) in pos_cat.items():
            if cat not in by_cat or d<by_cat[cat][1]: by_cat[cat]=(p,d,cb)
        anchors=sorted(by_cat.values(), key=lambda x:x[1])
        def sep_ok(combo): return all(np.linalg.norm(a[2]-b[2])>=SEP for a,b in itertools.combinations(combo,2))
        variants=[]
        for a in anchors: variants.append([a[0]])
        for k in (2,3,4):
            for combo in itertools.combinations(anchors,k):
                if sep_ok(combo): variants.append([c[0] for c in combo])
        dmap={p:d for p,d,_ in anchors}
        variants.sort(key=lambda v:(-len(v), np.mean([dmap[p] for p in v])))
        variants=variants[:MAXVAR]
        rows.append((f"{n}_rwt", s))
        for hp in variants:
            mut=list(s)
            for p in hp: mut[p-1]="H"
            nm=f"{n}_rh_"+"_".join(map(str,sorted(hp)))
            rows.append((nm,"".join(mut))); manifest.append((nm,n,sorted(hp)))
        print(f"{n}: {len(s)}mer anchors={[(p,round(d,1)) for p,d,_ in anchors[:6]]} -> {len(variants)} variants")
    with open(f"{R}/{out}/variants.csv","w",newline="") as fh:
        w=csv.writer(fh); [w.writerow(r) for r in rows]
    data=rows[1:]; NSH=8
    for k in range(NSH):
        with open(f"{R}/{out}/shard{k}.csv","w",newline="") as fh:
            w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in data[k::NSH]]
    with open(f"{R}/{out}/manifest.tsv","w") as fh:
        for nm,d,hp in manifest: fh.write(f"{nm}\t{d}\t{','.join(map(str,hp))}\n")
    print(f"\nwrote {len(data)} variants -> {out}/variants.csv + {NSH} shards")

if __name__=="__main__": main()
