#!/usr/bin/env python3
"""Select RFd3-K90 af2ig binders and emit a (name,sequence) CSV ready for switch-install + cross-react.

Reads runs/<tag>_af2ig/out.sc, keeps binders (pae_interaction<PAE, plddt_binder>80, binder_aligned_rmsd<2),
reads each binder's designed sequence from its af2ig INPUT pdb (chain A), dedups near-identical sequences,
and writes winners.csv. These then go: Boltz h+m (cross-react) -> install_sharp_his (His near conserved
Lys90) -> switch_sharpness -> novelty -> expression.

Usage: rfd3_select_winners.py <tag> [pae_cut=14] [max=40]
"""
import sys, os, glob, csv
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import warnings; warnings.filterwarnings("ignore")
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
tag=sys.argv[1] if len(sys.argv)>1 else "rfd3scale"
PAE=float(sys.argv[2]) if len(sys.argv)>2 else 14.0
MAXN=int(sys.argv[3]) if len(sys.argv)>3 else 40
sc=f"{R}/runs/{tag}_af2ig/out.sc"; afin=f"{R}/runs/{tag}_af2in"
rows=[]; hdr=None
for ln in open(sc):
    f=ln.split()
    if hdr is None: hdr={k:i for i,k in enumerate(f)}; continue
    try:
        pae=float(f[hdr["pae_interaction"]]); pl=float(f[hdr["plddt_binder"]]); rm=float(f[hdr["binder_aligned_rmsd"]])
    except: continue
    if pae<PAE and pl>80 and rm<2:
        rows.append((pae,pl,rm,f[hdr["description"]]))
rows.sort()
def binder_seq(desc):
    name=desc.replace("_af2pred","")
    p=f"{afin}/{name}.pdb"
    if not os.path.exists(p):
        g=glob.glob(f"{afin}/{name}*.pdb"); p=g[0] if g else None
    if not p: return None
    m=PDBParser(QUIET=True).get_structure("x",p)[0]
    ch=m["A"] if "A" in [c.id for c in m] else next(iter(m))
    return "".join(seq1(r.resname) for r in ch if r.id[0]==" ")
seen=[]; out=[]
for pae,pl,rm,desc in rows:
    s=binder_seq(desc)
    if not s or len(s)<40: continue
    # dedup: drop if >90% identical to a kept winner
    dup=False
    for ks in seen:
        if len(ks)==len(s):
            idt=sum(1 for a,b in zip(ks,s) if a==b)/len(s)
            if idt>0.90: dup=True; break
    if dup: continue
    seen.append(s); out.append((desc.replace("_af2pred",""),s,round(pae,2),round(pl,1),round(rm,2)))
    if len(out)>=MAXN: break
with open(f"{R}/runs/{tag}_af2ig/winners.csv","w") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","pae_interaction","plddt_binder","rmsd"])
    for r in out: w.writerow(r)
print(f"{len(rows)} binders (pae<{PAE}), {len(out)} unique -> runs/{tag}_af2ig/winners.csv")
for name,s,pae,pl,rm in out[:15]: print(f"  {name:38} pae={pae} plddt={pl} len={len(s)}")
