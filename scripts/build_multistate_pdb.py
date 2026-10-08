#!/usr/bin/env python3
"""Build a combined human+mouse complex for MULTI-STATE ProteinMPNN (cross-reactive binder design,
Chamaileon insight). For each Genie3 backbone (binder=chain A CA, human target=chain B), emit a 4-chain
PDB: A=binder, B=human-target, C=binder(copy), D=mouse-target (B coords, residue names threaded to mouse).
Preserves original ATOM-line column formatting; edits only chain ID (col 22) and resName (cols 18-20).
Also emits a tied_positions jsonl linking A_i<->C_i so MPNN designs ONE sequence for both orthologs.
Usage: build_multistate_pdb.py <backbone_glob> <out_pdb_dir> <tied_jsonl>"""
import sys, os, glob, json, re
from Bio.Align import PairwiseAligner
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
AA3={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
human=open(f"{R}/data/targets/human_tnf.seq").read().split()[0]
mouse=open(f"{R}/data/targets/mouse_tnf_mature.seq").read().split()[0]
al=PairwiseAligner(); al.mode='global'; al.open_gap_score=-10; al.extend_gap_score=-0.5; al.match_score=2; al.mismatch_score=-1
aln=al.align(human,mouse)[0]; hmap={}; hi=mi=0
for ch,cm in zip(aln[0],aln[1]):
    if ch!='-' and cm!='-': hmap[hi]=cm; hi+=1; mi+=1
    elif ch!='-': hmap[hi]=ch; hi+=1
    elif cm!='-': mi+=1
def setchain(ln,c): return ln[:21]+c+ln[22:]
def setresn(ln,r3): return ln[:17]+r3.rjust(3)+ln[20:]
outdir=sys.argv[2]; os.makedirs(outdir,exist_ok=True); tied={}
BB=[p for g in sys.argv[1].split() for p in glob.glob(g)]
for bb in BB:
    m=re.search(r"genie3_tnf_(g\d+).*?/([^/]+)\.pdb$", bb)
    name=(m.group(1)+"_"+m.group(2)) if m else os.path.basename(bb).replace(".pdb","")
    A=[ln for ln in open(bb) if ln.startswith("ATOM") and ln[21]=='A']
    B=[ln for ln in open(bb) if ln.startswith("ATOM") and ln[21]=='B']
    if not A or not B: continue
    binder_resis=sorted({int(ln[22:26]) for ln in A})
    out=[]
    for ln in A: out.append(ln)                                  # chain A = binder (human state)
    out.append("TER\n")
    for ln in B: out.append(ln)                                  # chain B = human target
    out.append("TER\n")
    for ln in A: out.append(setchain(ln,'C'))                    # chain C = binder copy (mouse state)
    out.append("TER\n")
    for ln in B:                                                 # chain D = mouse target (threaded names)
        resi=int(ln[22:26]); maa=hmap.get(resi-1)
        ln2=setresn(ln,AA3[maa]) if (maa and maa in AA3) else ln
        out.append(setchain(ln2,'D'))
    out.append("TER\nEND\n")
    open(f"{outdir}/{name}.pdb","w").write("".join(out))
    tied[name]=[{"A":[i],"C":[i]} for i in binder_resis]
json.dump(tied,open(sys.argv[3],"w"))
print(f"wrote {len(tied)} combined human+mouse PDBs -> {outdir}  + tied jsonl")
