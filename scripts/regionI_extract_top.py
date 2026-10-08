#!/usr/bin/env python3
"""Extract top-N Region-I binder sequences (chain A of the af2in PDB) for orthogonal Boltz cross-fold.
Dedup near-identical sequences (keep best-ranked per unique seq). Writes a name,sequence CSV.
Usage: regionI_extract_top.py <ranked.csv> <af2in_dir> <out.csv> [N]"""
import sys, csv, os, glob
RANK=sys.argv[1]; AF2IN=sys.argv[2]; OUT=sys.argv[3]; N=int(sys.argv[4]) if len(sys.argv)>4 else 12
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def binderseq(pdb):
    res={}
    for ln in open(pdb):
        if ln.startswith("ATOM") and ln[21]=='A' and ln[12:16].strip()=="CA":
            res[int(ln[22:26])]=AA3.get(ln[17:20].strip(),'X')
    return "".join(res[k] for k in sorted(res))
def find_af2in(desc):
    stem=desc.replace("_af2pred","")
    # desc like gpu2__tnf_regionI_tnf_regionI_8_model_13_4 -> af2in file same stem
    for p in [f"{AF2IN}/{stem}.pdb"]+glob.glob(f"{AF2IN}/*{stem.split('__',1)[-1]}.pdb"):
        if os.path.exists(p): return p
    return None
rows=list(csv.DictReader(open(RANK)))
seen={}; out=[]
for r in rows:
    if len(out)>=N: break
    pdb=find_af2in(r["design"])
    if not pdb: continue
    seq=binderseq(pdb)
    if not seq or 'X' in seq: continue
    if seq in seen: continue
    seen[seq]=r["design"]
    nm=r["design"].replace("_af2pred","").split("__",1)[-1]
    out.append((nm,seq,r.get("epi_contacts","?"),r.get("epi_near","?"),r.get("plddt_binder","?"),r.get("pae_int","?")))
with open(OUT,"w") as fh:
    fh.write("name,sequence\n")
    for nm,seq,ec,en,pl,pae in out: fh.write(f"{nm},{seq}\n")
print(f"extracted {len(out)} unique Region-I binders -> {OUT}")
for nm,seq,ec,en,pl,pae in out: print(f"  {nm:46} len={len(seq):3} epiC={ec} epiNear={en} plddtB={pl} pae={pae}")
