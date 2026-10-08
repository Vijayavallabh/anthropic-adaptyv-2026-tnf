#!/usr/bin/env python3
"""Epitope human<->mouse conservation audit (cross-reactivity determinant, Anthropic Aug-2026: x-react
tracks # substituted footprint residues). For each design's human complex: TNF epitope = target residues
within 5A of binder; align that TNF chain to mouse TNF; report conserved vs divergent footprint residues.
Usage: epitope_conservation.py <names_seqs.csv>"""
import sys, os, csv, glob
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner, substitution_matrices
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
MOUSE=open(f"{R}/data/targets/mouse_tnf.seq").read().split()[0]
SEARCH=[f"{R}/runs/boltz_push_h/out/by_seed",f"{R}/runs/boltz_fin_h/out/by_seed",f"{R}/runs/boltz_coop_h/out/by_seed"]
CUT=5.0
aln=PairwiseAligner(); aln.substitution_matrix=substitution_matrices.load("BLOSUM62")
aln.open_gap_score=-11; aln.extend_gap_score=-1; aln.mode="global"
def find_cif(name):
    for base in SEARCH:
        p=f"{base}/{name}/s0/{name}_model_0.cif"
        if os.path.exists(p): return p
    g=glob.glob(f"{R}/runs/boltz_*_h/out/by_seed/{name}/s0/{name}_model_0.cif"); return g[0] if g else None
def chainseq(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
def idf(a,b):
    n=min(len(a),len(b)); return (sum(1 for i in range(n) if a[i]==b[i])/n) if n else 0
def map_to_mouse(hseq):
    """align human TNF chain seq to mouse; return dict human_idx(0-based)->(mouse_res or '-')."""
    a=aln.align(hseq,MOUSE)[0]
    hi=mi=0; m={}
    ha,ma=a[0],a[1]
    for hc,mc in zip(ha,ma):
        if hc!="-":
            m[hi]=mc if mc!="-" else "-"
            hi+=1
        if mc!="-": mi+=1
    return m
for r in csv.DictReader(open(sys.argv[1])):
    name=r.get("name") or r.get("design_id"); seq=r["sequence"]
    cif=find_cif(name)
    if not cif: print(f"{name:20} NO COMPLEX"); continue
    model=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
    chains=[c for c in model if chainseq(c)]
    bch=max(chains,key=lambda c: idf(chainseq(c),seq))
    tchs=[c for c in chains if c.id!=bch.id]
    if not tchs: print(f"{name:20} single chain?"); continue
    tch=max(tchs,key=lambda c: len(chainseq(c)))  # TNF = largest non-binder chain
    tseq=chainseq(tch); tres=[rr for rr in tch if rr.id[0]==" "]
    batoms=[a for rr in bch if rr.id[0]==" " for a in rr if a.element!="H"]
    ns=NeighborSearch(batoms)
    epi=[i for i,rr in enumerate(tres) if any(ns.search(a.coord,CUT) for a in rr if a.element!="H")]
    m=map_to_mouse(tseq)
    cons=div=0; divlist=[]
    for i in epi:
        mc=m.get(i,"-"); hc=tseq[i] if i<len(tseq) else "?"
        if mc==hc: cons+=1
        else: div+=1; divlist.append(f"{hc}{i+1}>{mc}")
    pct=100*cons/len(epi) if epi else 0
    print(f"{name:20} epi={len(epi):2} conserved={cons:2} divergent={div} ({pct:4.0f}% cons) | div: {','.join(divlist) if divlist else 'none'}")
