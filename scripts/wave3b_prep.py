#!/usr/bin/env python3
"""WAVE-3B prep (interface-ADAPTED multi-His; the win point-mutation can't give). For each top balanced
binder, PRE-SET His at its 2-3 closest DISTINCT conserved-cation positions (resname->HIS, keep backbone),
then let LigandMPNN REDESIGN the interface shell AROUND the fixed His (core + His + target stay fixed) so
the interface packs to accommodate the multi-His. Cross-react is still enforced later by folding both
species and gating balanced (same as Track A). Writes runs/wave3b/<scaffold>.pdb + jobs.tsv + manifest.json.
Run with .venv/bin/python -I."""
import os, glob, csv, json
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
OUT="runs/wave3b"; os.makedirs(OUT, exist_ok=True)
NTOP=14; MAXHIS=3; CONTACT=8.0; SHELL=10.0; CHB_OFFSET=80
AA3={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE',
     'L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

h=open("runs/buried2/human_tnf.txt").read().strip(); m=open("runs/buried2/mouse_tnf_mature.txt").read().strip()
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load("BLOSUM62")
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode="global"; aln=al.align(h,m)[0]
hmap={}
for (hs,he),(ms,me) in zip(aln.aligned[0],aln.aligned[1]):
    for k in range(he-hs): hmap[hs+k+1]=m[ms+k]
conserved_cat={i for i,aa in enumerate(h,1) if aa in "KR" and hmap.get(i) in ("K","R")}

cand=[]
for l in open("runs/bpair/xreact_full2_results.txt"):
    p=l.split()
    if len(p)!=4 or p[0]=="design": continue
    cand.append((p[0],float(p[1]),float(p[2]),float(p[3]),abs(float(p[1])-float(p[2]))))
cand=[c for c in cand if c[3]>=0.55]; cand.sort(key=lambda c:-(c[3]-0.5*max(0.0,c[4]-0.10))); sel=cand[:NTOP]
def af2pred(d):
    c=glob.glob(f"{R}/runs/bpair/af2ig_full2/gpu*/pdbs/{d}_af2pred.pdb"); return c[0] if c else None

jobs=[]; mani={}
print(f"{'scaffold':14}{'nHis':>5}  His(pos->cat)   shellN")
for d,hh,mm,mn,dl in sel:
    f=af2pred(d)
    if not f: continue
    lines=open(f).read().splitlines()
    aCB={}; aCA={}; catN={}; bAtoms=[]
    for l in lines:
        if l[:4]!="ATOM": continue
        ch=l[21]; res=l[17:20].strip(); at=l[12:16].strip(); rn=int(l[22:26])
        crd=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        if ch=="A":
            if at=="CB": aCB[rn]=crd
            if at=="CA": aCA[rn]=(crd,res)
        if ch=="B":
            bAtoms.append(crd)
            hidx=rn-CHB_OFFSET
            if hidx in conserved_cat and ((res=="LYS" and at=="NZ") or (res=="ARG" and at in("NH1","NH2","NE"))):
                catN.setdefault(hidx,[]).append(crd)
    bArr=np.array(bAtoms)
    # choose His: closest binder pos per distinct conserved cation (skip Gly/Pro)
    pos2cat={}
    for rn,crd in aCB.items():
        resn=aCA.get(rn,(None,"X"))[1]
        if resn in ("GLY","PRO"): continue
        best=None
        for hidx,coords in catN.items():
            dmn=min(float(np.linalg.norm(crd-c)) for c in coords)
            if dmn<CONTACT and (best is None or dmn<best[1]): best=(hidx,dmn)
        if best: pos2cat[rn]=best
    cat2pos={}
    for rn,(hidx,dist) in pos2cat.items():
        if hidx not in cat2pos or dist<cat2pos[hidx][1]: cat2pos[hidx]=(rn,dist)
    chosen=sorted(cat2pos.items(), key=lambda kv:kv[1][1])[:MAXHIS]
    his=[rn for hidx,(rn,dist) in chosen]; his_cat={rn:hidx for hidx,(rn,dist) in chosen}
    if len(his)<2: print(f"  {d:12} <2 distinct cations, skip"); continue
    # interface shell = binder CA within SHELL of any target atom, minus His, minus Gly/Pro
    shell=[]
    for rn,(crd,res) in aCA.items():
        if rn in his or res in ("GLY","PRO"): continue
        if float(np.min(np.linalg.norm(bArr-crd,axis=1)))<SHELL: shell.append(rn)
    shell=sorted(shell)
    # write scaffold PDB: keep chain A + B; at His positions set resname HIS and keep only backbone+CB
    outp=f"{OUT}/{d}.pdb"; bb={"N","CA","C","O","CB"}
    with open(outp,"w") as fo:
        for l in lines:
            if l[:4] not in ("ATOM","HETA"):
                if l.startswith("TER"): fo.write(l+"\n")
                continue
            ch=l[21]; rn=int(l[22:26]); at=l[12:16].strip()
            if ch=="A" and rn in his:
                if at not in bb: continue
                l=l[:17]+"HIS"+l[20:]
            fo.write(l+"\n")
        fo.write("END\n")
    redes=" ".join(f"A{rn}" for rn in shell)
    jobs.append((d, outp, redes))
    mani[d]=dict(scaffold=d, his=his, his_cat=his_cat, shell=shell, min0=mn, dHM0=dl, pdb=outp)
    print(f"  {d:12}{len(his):>5}  "+",".join(f"{rn}->{his_cat[rn]}" for rn in his)+f"   {len(shell)}")

with open(f"{OUT}/jobs.tsv","w") as fo:
    for d,p,redes in jobs: fo.write(f"{d}\t{p}\t{redes}\n")
json.dump(mani, open(f"{OUT}/manifest.json","w"), indent=1)
print(f"\n[wave3b] {len(jobs)} scaffolds prepped (fixed multi-His + interface shell to redesign) -> {OUT}/jobs.tsv")
