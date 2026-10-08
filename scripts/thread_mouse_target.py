#!/usr/bin/env python3
"""Thread the mature MOUSE TNF sequence onto the human 1TNF_B backbone (79% id, same fold) to make a
mouse target for MULTI-STATE ProteinMPNN (cross-reactivity by construction, per Chamaileon insight).
MPNN uses only backbone N/CA/C/O + residue IDENTITY of the fixed target chain, so renaming residues on
the identical backbone is sufficient. Writes data/targets/1TNF_B_mouse.pdb."""
import sys
from Bio.Align import PairwiseAligner
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
AA3={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
AA1={v:k for k,v in AA3.items()}
human=open(f"{R}/data/targets/human_tnf.seq").read().split()[0]        # 1TNF_B chain-B seq (crystal res 6..157)
mouse=open(f"{R}/data/targets/mouse_tnf_mature.seq").read().split()[0]  # mature mouse (res 1..)
# align: find the mouse sub-sequence matching human_tnf.seq (human starts at mature R6; mouse has VRSSSR... too)
al=PairwiseAligner(); al.mode='global'; al.open_gap_score=-10; al.extend_gap_score=-0.5; al.match_score=2; al.mismatch_score=-1
aln=al.align(human,mouse)[0]
# build human_index -> mouse_residue map from the alignment
hmap={}  # human position (0-based in human_tnf.seq) -> mouse AA (or None for gap)
hi=mi=0
a0,a1=aln[0],aln[1]
for c_h,c_m in zip(a0,a1):
    if c_h!='-' and c_m!='-': hmap[hi]=c_m; hi+=1; mi+=1
    elif c_h!='-' and c_m=='-': hmap[hi]=c_h; hi+=1   # no mouse counterpart -> keep human (rare)
    elif c_h=='-' and c_m!='-': mi+=1
ident=sum(1 for i in range(len(human)) if hmap.get(i)==human[i])/len(human)
print(f"human_tnf len={len(human)} mouse len={len(mouse)} identity(threaded)={ident:.1%}")
# 1TNF_B: chain B, residues 6..157 (crystal). human_tnf.seq[i] = residue (6+i). rename to mouse.
out=[]; nren=0; res_order=[]
seen=set()
# first pass: residue order in file to assign human_tnf index by CA order
order=[]
for ln in open(f"{R}/data/targets/1TNF_B.pdb"):
    if ln.startswith("ATOM") and ln[12:16].strip()=="CA" and ln[21]=='B':
        order.append(int(ln[22:26]))
res2idx={rn:i for i,rn in enumerate(sorted(set(order)))}   # i-th residue -> human_tnf.seq[i]
for ln in open(f"{R}/data/targets/1TNF_B.pdb"):
    if ln.startswith("ATOM") and ln[21]=='B':
        rn=int(ln[22:26]); idx=res2idx.get(rn)
        if idx is not None and idx in hmap:
            maa=hmap[idx]; ln=ln[:17]+AA3[maa].ljust(3)+ln[20:]
            if rn not in seen and maa!=human[idx]: nren+=1
            seen.add(rn)
    out.append(ln)
open(f"{R}/data/targets/1TNF_B_mouse.pdb","w").write("".join(out))
print(f"wrote 1TNF_B_mouse.pdb  ({nren} residues renamed vs human; backbone identical)")
