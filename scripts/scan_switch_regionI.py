#!/usr/bin/env python3
"""Region-I single-His release-switch scan: binder=chain A, target=chains B,C (two protomers).
For each binder residue whose CB is 4.5-9A from a CONSERVED target Lys NZ (K65/K98/K112, which in the
af2ig-renumbered 1..152 frame = 60/93/107, on EITHER protomer B or C), emit a His variant (bb+CB kept,
af2ig rebuilds the ring). Conserved anchors => the His-repulsion switch is cross-react-preserving.
Skips Gly/Pro/His. Usage: scan_switch_regionI.py <complex.pdb> <outdir> [maxpos]"""
import sys, os, numpy as np
from collections import defaultdict
aa3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
     'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
fn=sys.argv[1]; outdir=sys.argv[2]; maxpos=int(sys.argv[3]) if len(sys.argv)>3 else 8
os.makedirs(outdir,exist_ok=True)
CONS_RENUM={60,93,107}   # K65,K98,K112 in the renumbered (orig-5) target frame; all human/mouse conserved
atA=defaultdict(list); resnA={}; tgt=defaultdict(list); resnT={}
for ln in open(fn):
    if ln[:4]!="ATOM": continue
    ch=ln[21]; rn=int(ln[22:26]); nm=ln[12:16].strip()
    xyz=np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    if ch=='A': atA[rn].append((nm,xyz)); resnA[rn]=aa3.get(ln[17:20].strip(),'X')
    elif ch in ('B','C'): tgt[(ch,rn)].append((nm,xyz)); resnT[(ch,rn)]=aa3.get(ln[17:20].strip(),'X')
# conserved basic target N atoms on either protomer
Aall=np.array([x for rn in atA for _,x in atA[rn]]) if atA else np.zeros((1,3))
basicN=[]
for (ch,rn),atoms in tgt.items():
    if rn in CONS_RENUM and resnT[(ch,rn)]=='K':
        for nm,x in atoms:
            if nm=='NZ' and np.sqrt(((x-Aall)**2).sum(1)).min()<12:
                basicN.append((ch,rn,f"{ch}K{rn+5}",x))
cands=[]
for arn,atoms in atA.items():
    if resnA[arn] in ('G','P','H'): continue
    cb=[x for nm,x in atoms if nm=='CB'] or [x for nm,x in atoms if nm=='CA']
    for ch,brn,blabel,bx in basicN:
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
print(f"{base}: {len(picked)} His variants near conserved Lys -> {outdir}")
for arn,aa,bl,d in picked: print(f"  A{arn}{aa}->H  (CB {d:.1f}A from {bl})")
