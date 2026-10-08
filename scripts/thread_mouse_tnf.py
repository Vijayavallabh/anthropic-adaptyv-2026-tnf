#!/usr/bin/env python3
"""Build binder+MOUSE-TNF complexes for a matched-method cross-reactivity af2ig co-fold.
Threads the mouse TNF-alpha sequence onto chain B of each human complex PDB: keep binder (chain A)
in full; for chain B keep backbone N/CA/C/O only and rename each residue to its mouse identity
(af2ig rebuilds the mouse side chains). Human/mouse epitope is conserved, so a cross-reactive binder
folds back with pae like the human. Chain-B order index k (1-based) -> human mature position k+5
(1TNF_B starts at human residue 6); mouse identity from a global NW alignment.
Usage: thread_mouse_tnf.py <indir of human complex PDBs> <outdir>
"""
import sys, os, glob, numpy as np
HUMAN="VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
MOUSE="".join(l.strip() for l in open("data/targets/mouse_TNF_P06804.fasta") if not l.startswith(">"))
AA1to3={"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU","G":"GLY","H":"HIS",
        "I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
BB={"N","CA","C","O"}
def nw(a,b,gap=-2):
    n,m=len(a),len(b);H=np.zeros((n+1,m+1))
    for i in range(n+1):H[i][0]=i*gap
    for j in range(m+1):H[0][j]=j*gap
    for i in range(1,n+1):
        for j in range(1,m+1):
            s=1 if a[i-1]==b[j-1] else -1
            H[i][j]=max(H[i-1][j-1]+s,H[i-1][j]+gap,H[i][j-1]+gap)
    i,j=n,m;A=[];B=[]
    while i>0 or j>0:
        if i>0 and j>0 and H[i][j]==H[i-1][j-1]+(1 if a[i-1]==b[j-1] else -1):A.append(a[i-1]);B.append(b[j-1]);i-=1;j-=1
        elif i>0 and H[i][j]==H[i-1][j]+gap:A.append(a[i-1]);B.append('-');i-=1
        else:A.append('-');B.append(b[j-1]);j-=1
    return A[::-1],B[::-1]
# human mature position -> mouse residue (1-letter), '-' if gap
ah,am=nw(HUMAN,MOUSE); hp=0; mouse_at={}
for ca,cm in zip(ah,am):
    if ca!='-': hp+=1; mouse_at[hp]=cm
def thread(src,dst):
    # chain B residues in order -> human pos (k+5) -> mouse residue
    bchain=[]; last=None
    out=[]
    # first pass: collect chain B residue order
    seen=[]
    for ln in open(src):
        if ln[:4]=="ATOM" and ln[21]=='B':
            rn=int(ln[22:26])
            if not seen or seen[-1]!=rn: seen.append(rn)
    order={rn:i+1 for i,rn in enumerate(seen)}   # 1-based order
    nmut=0
    for ln in open(src):
        if ln[:4]!="ATOM": out.append(ln); continue
        ch=ln[21]
        if ch=='A': out.append(ln); continue
        if ch!='B': continue
        rn=int(ln[22:26]); atom=ln[12:16].strip()
        if atom not in BB: continue           # drop side chain
        k=order[rn]; hpos=k+5; maa=mouse_at.get(hpos)
        if maa and maa!='-':
            res3=AA1to3.get(maa, ln[17:20])
            out.append(ln[:17]+f"{res3:>3}"+ln[20:]);
        else:
            out.append(ln)
    open(dst,'w').write("".join(out))
indir,outdir=sys.argv[1],sys.argv[2]; os.makedirs(outdir,exist_ok=True); n=0
for p in sorted(glob.glob(os.path.join(indir,"*.pdb"))):
    thread(p, os.path.join(outdir, os.path.basename(p))); n+=1
print(f"threaded {n} complexes -> mouse TNF (bb + mouse seq) in {outdir}")
