#!/usr/bin/env python3
"""Emit single-His release-switch variants of a binder:target complex: for each binder (chain A)
residue whose CB is 4.5-9A from a CONSERVED target (chain B) basic sidechain N (K65/K90/K112 ->
renum 60/85/107; also Arg NH/NE), write a variant with that position mutated to HIS (bb+CB kept,
af2ig rebuilds). Skips Gly/Pro. Usage: scan_switch.py <complex.pdb> <outdir> [maxpos]"""
import sys,os,numpy as np
from collections import defaultdict
aa3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
fn=sys.argv[1]; outdir=sys.argv[2]; maxpos=int(sys.argv[3]) if len(sys.argv)>3 else 8
os.makedirs(outdir,exist_ok=True)
atA=defaultdict(list); atB=defaultdict(list); resnA={}; resnB={}
for ln in open(fn):
    if ln[:4]!="ATOM": continue
    ch=ln[21]; rn=int(ln[22:26]); nm=ln[12:16].strip(); xyz=np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    if ch=='A': atA[rn].append((nm,xyz)); resnA[rn]=aa3.get(ln[17:20].strip(),'X')
    elif ch=='B': atB[rn].append((nm,xyz)); resnB[rn]=aa3.get(ln[17:20].strip(),'X')
# conserved basic target residues (human=mouse): K65/K90/K98 (renum 60/85/93) + R32 (renum 27).
# DISTANCE-GATED: only count a cation whose sidechain-N is within 9A of the binder (far conserved
# cations must NOT block near ones). Fallback adds any conserved-ish Lys/Arg near the binder.
cons_K={60,85,93}; cons_R={27}
Aall=np.array([x for rn in atA for _,x in atA[rn]])
basicN=[]
for rn,atoms in atB.items():
    Ns=[(nm,x) for nm,x in atoms if nm in ('NZ','NH1','NH2','NE')]
    if not Ns: continue
    dmin=min(np.sqrt(((x-Aall)**2).sum(1)).min() for _,x in Ns)
    if dmin>9.0: continue
    if resnB[rn]=='K' and rn in cons_K:
        for nm,x in Ns:
            if nm=='NZ': basicN.append((rn,'K'+str(rn+5),x))
    elif resnB[rn]=='R' and rn in cons_R:
        for nm,x in Ns: basicN.append((rn,'R'+str(rn+5),x))
if not basicN:  # fallback: any Lys NZ within 9A (still cross-react-biased; excludes N-term K<10)
    for rn,atoms in atB.items():
        if resnB[rn]=='K' and rn>=10:
            nz=[x for nm,x in atoms if nm=='NZ']
            if nz and np.sqrt(((nz[0]-Aall)**2).sum(1)).min()<9: basicN.append((rn,'K%d'%(rn+5),nz[0]))
cands=[]
for arn,atoms in atA.items():
    if resnA[arn] in ('G','P','H'): continue
    cb=[x for nm,x in atoms if nm=='CB'] or [x for nm,x in atoms if nm=='CA']
    for brn,blabel,bx in basicN:
        d=np.linalg.norm(cb[0]-bx)
        if 4.5<d<9: cands.append((d,arn,resnA[arn],blabel))
cands.sort()
seen=set(); picked=[]
for d,arn,aa,bl in cands:
    if arn in seen: continue
    seen.add(arn); picked.append((arn,aa,bl,d))
    if len(picked)>=maxpos: break
keep_bb={'N','CA','C','O','CB'}
base=os.path.basename(fn).replace('_af2pred','').replace('.pdb','')
for arn,aa,bl,d in picked:
    out=[]
    for ln in open(fn):
        if ln[:4]=="ATOM":
            if ln[21]=='A' and int(ln[22:26])==arn:
                if ln[12:16].strip() not in keep_bb: continue
                ln=ln[:17]+"HIS"+ln[20:]
            out.append(ln)
    open(os.path.join(outdir,f"{base}__A{arn}{aa}2H_vs{bl}.pdb"),'w').write("".join(out)+"END\n")
print(f"{base}: {len(picked)} His variants near conserved basics -> {outdir}")
for arn,aa,bl,d in picked: print(f"  A{arn}{aa}->H  (CB {d:.1f}A from {bl})")
