#!/usr/bin/env python3
"""Deterministic interface-strength metrics (the EGFR-post-mortem signals that beat ipSAE):
buried SASA (dSASA), inter-chain heavy-atom contacts, H-bond-capable contacts, interface residues.
Complex PDB: chain A=binder, B=target. Run in .venv (freesasa+numpy)."""
import sys, glob, os, numpy as np, freesasa
freesasa.setVerbosity(freesasa.silent)
def atoms(pdb, chain):
    out=[]
    for ln in open(pdb):
        if ln[:4] in ("ATOM","HETA") and ln[21]==chain:
            el=ln[76:78].strip() or ln[12:16].strip()[0]
            out.append((int(ln[22:26]), ln[12:16].strip(), el,
                        np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])))
    return out
def sasa(pdb, keepchains):
    # write a temp pdb with only keepchains, compute total SASA
    import tempfile
    tmp=tempfile.NamedTemporaryFile('w',suffix='.pdb',delete=False)
    for ln in open(pdb):
        if ln[:4] in ("ATOM","TER") :
            if ln[21] in keepchains or ln[:3]=="TER": tmp.write(ln)
    tmp.write("END\n"); tmp.close()
    try:
        st=freesasa.Structure(tmp.name); r=freesasa.calc(st); v=r.totalArea()
    except Exception: v=float('nan')
    os.unlink(tmp.name); return v
def metrics(pdb):
    A=atoms(pdb,'A'); B=atoms(pdb,'B')
    if not A or not B: return None
    cA=np.array([a[3] for a in A]); cB=np.array([b[3] for b in B])
    D=np.linalg.norm(cA[:,None,:]-cB[None,:,:],axis=2)   # heavy-atom dists
    # residue-pair contacts <5A
    cpairs=set(); hpairs=0; ires=set()
    hbA={'N','O'};
    ii,jj=np.where(D<5.0)
    for i,j in zip(ii,jj):
        cpairs.add((A[i][0],B[j][0])); ires.add(A[i][0])
    ii,jj=np.where(D<3.5)
    for i,j in zip(ii,jj):
        if A[i][2] in hbA and B[j][2] in hbA: hpairs+=1
    sA=sasa(pdb,'A'); sB=sasa(pdb,'B'); sAB=sasa(pdb,'AB')
    dsasa=(sA+sB-sAB) if not any(np.isnan([sA,sB,sAB])) else float('nan')
    return dict(ncontact=len(cpairs), nhbond=hpairs, nires=len(ires), dsasa=dsasa)
print(f"{'design':40}{'dSASA':>8}{'contacts':>9}{'hbond':>6}{'i_res':>6}")
rows=[]
for pdb in sorted(sys.argv[1:]):
    m=metrics(pdb); nm=os.path.basename(pdb).replace('.pdb','')
    if m: rows.append((nm,m));
for nm,m in sorted(rows,key=lambda r:-r[1]['dsasa'] if not np.isnan(r[1]['dsasa']) else 0):
    print(f"{nm[:39]:40}{m['dsasa']:8.0f}{m['ncontact']:9d}{m['nhbond']:6d}{m['nires']:6d}")
