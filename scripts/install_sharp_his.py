#!/usr/bin/env python3
"""Install a SHARP multi-His switch (Wyman Delta-n_H -> 2-3) per the pH-binder physics
(bioRxiv 2025.09.29.678932): 2-3 His each near a DISTINCT target cation (Lys6/Lys93/Arg on TNF),
SPATIALLY SEPARATED (>SEP A, Cb-Cb) so they titrate additively instead of anti-cooperatively splitting.
Keeps existing His; biases replacement toward polar/small interface positions (not hydrophobic anchors).
For each design it emits a few geometric variants -> fold (Boltz h+m) -> switch_sharpness -> keep sharpest.
Usage: install_sharp_his.py <names_seqs.csv> <out_dir>"""
import sys, os, csv, glob, itertools
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
import numpy as np
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
SEARCH=[f"{R}/runs/boltz_push_h/out/by_seed",f"{R}/runs/boltz_fin_h/out/by_seed",f"{R}/runs/boltz_coop_h/out/by_seed"]
if os.environ.get("BOLTZ_SEARCH"): SEARCH=os.environ["BOLTZ_SEARCH"].split(":")+SEARCH  # extra complex dir(s), colon-separated (e.g. RFd3-winner Boltz folds)
SEP=8.0           # min Cb-Cb separation between switch His (avoid anti-cooperativity)
REACH=(3.0,8.5)   # binder Cb within this of a target cation charge atom -> His can reach it
CAT={"LYS":["NZ"],"ARG":["NH1","NH2","CZ","NE"]}
# residues worth replacing (polar/charged/small; NOT hydrophobic binding anchors)
REPLACEABLE=set("ASTNQDEKRGH".replace("G",""))  # A S T N Q D E K R H (keep existing H, skip G/P/C and hydrophobics)
def find_cif(n):
    for b in SEARCH:
        p=f"{b}/{n}/s0/{n}_model_0.cif"
        if os.path.exists(p): return p
    g=glob.glob(f"{R}/runs/boltz_*_h/out/by_seed/{n}/s0/{n}_model_0.cif"); return g[0] if g else None
def cseq(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
def idf(a,b):
    n=min(len(a),len(b)); return (sum(1 for i in range(n) if a[i]==b[i])/n) if n else 0
def cb(res):
    if "CB" in res: return res["CB"].coord
    if "CA" in res: return res["CA"].coord
    return None
def gen(name,seq,out_rows,man):
    cif=find_cif(name)
    if not cif: print(f"{name}: no complex"); return
    model=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
    chains=[c for c in model if cseq(c)]
    bch=max(chains,key=lambda c: idf(cseq(c),seq))
    if idf(cseq(bch),seq)<0.5: print(f"{name}: no chain match"); return
    tch=max([c for c in chains if c.id!=bch.id],key=lambda c:len(cseq(c)))
    bres=[r for r in bch if r.id[0]==" "]
    # target cation charge atoms
    cations=[]
    for r in tch:
        if r.id[0]!=" " or r.resname not in CAT: continue
        for an in CAT[r.resname]:
            if an in r: cations.append((f"{r.resname}{r.id[1]}", r[an].coord))
    # candidate binder positions: Cb reachable to a cation, replaceable residue, not existing His
    existing_his=[i for i,c in enumerate(seq) if c=="H"]
    cand={}  # pos(0-based) -> (nearest_cation_label, dist)
    for i,res in enumerate(bres):
        if i>=len(seq) or seq[i]=="H": continue
        if seq[i] not in REPLACEABLE: continue
        p=cb(res)
        if p is None: continue
        best=None
        for lbl,cc in cations:
            d=float(np.linalg.norm(p-cc))
            if REACH[0]<=d<=REACH[1] and (best is None or d<best[1]): best=(lbl,d)
        if best: cand[i]=best
    if len(cand)<1: print(f"{name}: no reachable candidate positions"); return
    coords={i:cb(bres[i]) for i in list(cand)+existing_his if i<len(bres) and cb(bres[i]) is not None}
    def sep_ok(ps):
        for a,b in itertools.combinations(ps,2):
            if a in coords and b in coords and np.linalg.norm(coords[a]-coords[b])<SEP: return False
        return True
    # build variant His-sets: existing His + 1-2 new, near DISTINCT cations, separated
    cand_sorted=sorted(cand, key=lambda i:cand[i][1])  # closest-to-cation first
    sets=[]
    # (a) existing + best 1 new (near a distinct cation from existing)
    for i in cand_sorted:
        ps=sorted(set(existing_his)|{i})
        if sep_ok(ps) and ps not in [s[0] for s in sets]:
            sets.append((ps,[i]));
        if len(sets)>=3: break
    # (b) existing + 2 new near two DISTINCT cations, separated from each other + existing
    for i,j in itertools.combinations(cand_sorted,2):
        if cand[i][0]==cand[j][0]: continue   # must target different cations
        ps=sorted(set(existing_his)|{i,j})
        if sep_ok(ps) and ps not in [s[0] for s in sets]:
            sets.append((ps,[i,j]))
        if len(sets)>=6: break
    for k,(ps,new) in enumerate(sets):
        s=list(seq)
        for i in ps: s[i]="H"
        v="".join(s)
        vid=f"{name}_sh{k+1}"
        out_rows.append((vid,v))
        cats=";".join(f"{cand[i][0]}@{cand[i][1]:.1f}" for i in new if i in cand)
        minsep=min((np.linalg.norm(coords[a]-coords[b]) for a,b in itertools.combinations([p for p in ps if p in coords],2)),default=0)
        man.append((vid,name,"+".join(str(p+1) for p in ps),"+".join(str(p+1) for p in new),cats,f"{minsep:.1f}"))
        print(f"  {vid}: His@{[p+1 for p in ps]} new={[p+1 for p in new]} near[{cats}] minsep={minsep:.1f}A nHis={v.count('H')}")
os.makedirs(sys.argv[2],exist_ok=True)
rows=[]; man=[]
for r in csv.DictReader(open(sys.argv[1])):
    name=r.get("name") or r.get("design_id"); seq=r["sequence"]
    print(f"{name} (L={len(seq)}):")
    gen(name,seq,rows,man)
with open(f"{sys.argv[2]}/variants.csv","w") as fh:
    fh.write("name,sequence\n")
    for n,s in rows: fh.write(f"{n},{s}\n")
with open(f"{sys.argv[2]}/manifest.tsv","w") as fh:
    fh.write("variant\tparent\this_positions\tnew_his\tnear_cations\tmin_sep_A\n")
    for row in man: fh.write("\t".join(row)+"\n")
print(f"\nwrote {len(rows)} sharp-His variants -> {sys.argv[2]}/variants.csv")
