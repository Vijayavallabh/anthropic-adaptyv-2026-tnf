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
# conserved basic target N atoms (K65/K90/K112 by renum = orig-5; plus any Arg that is conserved-ish)
cons_renum={60,85,107}
basicN=[]
for rn,atoms in atB.items():
    if rn in cons_renum and resnB[rn]=='K':
        for nm,x in atoms:
            if nm=='NZ': basicN.append((rn,'K'+str(rn+5),x))
if not basicN:  # fallback: any Lys NZ near binder
    Aall=np.array([x for rn in atA for _,x in atA[rn]])
    for rn,atoms in atB.items():
        if resnB[rn]=='K':
            nz=[x for nm,x in atoms if nm=='NZ']
            if nz and np.sqrt(((nz[0]-Aall)**2).sum(1)).min()<10: basicN.append((rn,'K%d'%(rn+5),nz[0]))
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
