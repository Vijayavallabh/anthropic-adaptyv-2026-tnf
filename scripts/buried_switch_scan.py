#!/usr/bin/env python3
"""Buried-network release switch (Baker mechanism ii, self-implemented, validatable readout).
Place a His in the binder CORE next to a BURIED binder cation (Lys/Arg): at low pH both protonate
-> buried like-charge repulsion -> local fold destabilization -> release. OFF-INTERFACE, so binding
and cross-reactivity are preserved. Readout = His-cation distance two-state MD (my proven d_repel).
For each binder (chain A) position whose CB is 4.5-8A from a BURIED (iface>8A) binder K/R sidechain-N,
and which is itself buried (>=50 atoms within 10A) + off-interface (iface>8A) + non-terminal + not Gly/Pro,
emit a variant with that position -> HIS (bb+CB kept; af2ig/MPNN rebuild). Skips the anchor itself.
Usage: buried_switch_scan.py <complex.pdb> <outdir> [maxpos]"""
import sys, os, numpy as np
from collections import defaultdict
aa3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
fn=sys.argv[1]; outdir=sys.argv[2]; maxpos=int(sys.argv[3]) if len(sys.argv)>3 else 10
os.makedirs(outdir, exist_ok=True)
atA=defaultdict(dict); resnA={}; allA=[]; B=[]
for ln in open(fn):
    if ln[:4]!="ATOM": continue
    xyz=np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])]); nm=ln[12:16].strip()
    if ln[21]=='A': atA[int(ln[22:26])][nm]=xyz; resnA[int(ln[22:26])]=aa3.get(ln[17:20].strip(),'X'); allA.append(xyz)
    elif ln[21]=='B': B.append(xyz)
allA=np.array(allA); B=np.array(B); rns=sorted(atA)
def cb(r): return atA[r].get('CB',atA[r].get('CA'))
def burial(c): return int((np.linalg.norm(allA-c,axis=1)<10).sum())
def iface(c): return float(np.linalg.norm(B-c,axis=1).min())
term=set(rns[:8])|set(rns[-8:])
# buried binder cations = anchors
anchors=[]
for r in rns:
    if resnA[r] in 'KR':
        Ns=[atA[r][n] for n in ('NZ','NH1','NH2','NE') if n in atA[r]]
        if Ns and iface(cb(r))>8 and burial(cb(r))>60:
            anchors.append((r,Ns))
base=os.path.basename(fn).replace('.pdb','')
written=0; cand=[]
for r in rns:
    if r in term or resnA[r] in 'GP' or r in [a for a,_ in anchors]: continue
    c=cb(r)
    if burial(c)<50 or iface(c)<8: continue  # must be buried core, off-interface
    # nearest buried-cation N
    best=None
    for ar,Ns in anchors:
        d=min(np.linalg.norm(c-n) for n in Ns)
        if 4.5<d<8.0 and (best is None or d<best[1]): best=(ar,d)
    if best: cand.append((r,best[0],round(best[1],1)))
cand.sort(key=lambda x:x[2])
for r,ar,d in cand[:maxpos]:
    out=[]
    for ln in open(fn):
        if ln[:4]=="ATOM" and ln[21]=='A' and int(ln[22:26])==r:
            if ln[12:16].strip() in ('N','CA','C','O','CB'):
                out.append(ln[:17]+'HIS'+ln[20:])
        else: out.append(ln)
    open(f"{outdir}/{base}__A{r}{resnA[r]}2H_vs{resnA[ar]}{ar}.pdb","w").write("".join(out))
    written+=1
    print(f"  A{r}{resnA[r]}->H  (CB {d}A from buried {resnA[ar]}{ar})")
print(f"{base}: {written} buried-His switch variants (anchors: {[(a,resnA[a]) for a,_ in anchors]}) -> {outdir}")
