#!/bin/bash
# Select Genie3-regI binders from the funnel's af2ig output for the switch stage. af2ig is a WEAK ranker for
# TNF (Overath) + Boltz ipSAE proved unusable for the Genie3 fold-type (calibration), so this is a LENIENT,
# OVERSAMPLE gate: rank by af2ig pae_interaction, keep the best with a real interface (binder contacts the
# target) and decent fold. Output binders.csv (name, chain-A seq) for switch_screen_genie3_regI.sh.
# The af2pred poses are already A + merged-B(target) — switch-build + MD use them directly (chain B = target).
# Usage: select_genie3_regI.sh [pae_max=16] [topN=40] [min_contacts=4] [tag=genie3_regI]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
PAEMAX="${1:-16}"; TOPN="${2:-40}"; MINC="${3:-4}"; TAG="${4:-genie3_regI}"
SC="$REPO/runs/${TAG}/af2ig/out.sc"; POSEDIR="$REPO/runs/${TAG}/af2ig"
OUT="$REPO/runs/${TAG}/select"; mkdir -p "$OUT"
[ -s "$SC" ] || { echo "no af2ig out.sc yet"; exit 1; }
envs/uplifting-esmc/bin/python - "$SC" "$POSEDIR" "$OUT/binders.csv" "$PAEMAX" "$TOPN" "$MINC" <<'PY'
import sys,glob,os,re,numpy as np
sc,posedir,outcsv,paemax,topn,minc=sys.argv[1],sys.argv[2],sys.argv[3],float(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6])
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def pose(desc):
    g=glob.glob(f"{posedir}/gpu*/pdbs/{desc}_af2pred.pdb")
    return g[0] if g else None
def chainA(pdb):
    seq={}
    for l in open(pdb):
        if l.startswith("ATOM") and l[21]=="A" and l[12:16].strip()=="CA":
            seq[int(l[22:26])]=AA3.get(l[17:20].strip(),"X")
    return "".join(seq[k] for k in sorted(seq))
def contacts(pdb):
    ca={}
    for l in open(pdb):
        if l.startswith("ATOM") and l[12:16].strip()=="CA": ca.setdefault(l[21],[]).append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    A=np.array(ca.get("A",[])); B=np.array(ca.get("B",[]))
    if not len(A) or not len(B): return 0,99.0
    d=np.sqrt(((A[:,None,:]-B[None,:,:])**2).sum(-1)); return int((d<8).any(1).sum()), float(d.min())
rows=[]
for ln in open(sc):
    if not ln.startswith("SCORE:") or "description" in ln: continue
    f=ln.split()
    try: rmsd,pae_int,plddt_b,desc=float(f[1]),float(f[3]),float(f[5]),f[-1]
    except: continue
    if pae_int>paemax or plddt_b<75: continue
    p=pose(desc)
    if not p: continue
    nc,md=contacts(p)
    if nc<minc: continue
    bb=re.sub(r'_\d+$','',desc)   # backbone = design minus trailing _seqnum
    rows.append((pae_int,desc,bb,rmsd,plddt_b,nc,md,p))
# best per backbone, then top-N by pae
best={}
for r in sorted(rows,key=lambda x:x[0]):
    if r[2] not in best: best[r[2]]=r
ranked=sorted(best.values(),key=lambda x:x[0])[:topn]
import csv
w=csv.writer(open(outcsv,"w")); w.writerow(["name","sequence","pae_int","plddt_b","rmsd","contacts","min_d"])
posepaths=open(outcsv.replace(".csv","_poses.txt"),"w")
print(f"{'design':30}{'pae_int':>8}{'plddt':>7}{'rmsd':>6}{'contacts':>9}")
for pae,desc,bb,rmsd,pl,nc,md,p in ranked:
    seq=chainA(p)
    if len(seq)<50: continue
    w.writerow([desc,seq,round(pae,2),round(pl,1),round(rmsd,2),nc,round(md,1)]); posepaths.write(p+"\n")
    print(f"{desc[:30]:30}{pae:>8.2f}{pl:>7.1f}{rmsd:>6.2f}{nc:>9}")
print(f"\nscored-rows passing(pae<{paemax},plddt>75,contacts>={minc})={len(rows)} backbones={len(best)} -> wrote {min(len(ranked),topn)} to {outcsv}")
PY
echo "SELECT_GENIE3_REGI_DONE"
