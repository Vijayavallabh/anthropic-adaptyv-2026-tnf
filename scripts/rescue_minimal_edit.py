#!/usr/bin/env python3
"""Minimal-edit expression fix: fewest SURFACE (non-interface,non-His) point substitutions to reach
instab<40 AND GRAVY<0.4, preserving the fold (far less perturbation than MPNN redesign -> keep mouse).
Reads runs/push/rescue/manifest.tsv. Writes runs/push/rescue/minedit.fasta + .csv."""
import csv, re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; OUT=f"{R}/runs/push/rescue"
PANEL=list("DEKRNQST")          # hydrophilic, low-aggregation surface residues
NMAX=8
def penalty(seq):
    pa=ProteinAnalysis(seq); ins=pa.instability_index(); gr=pa.gravy()
    return max(0,ins-40)/10.0 + max(0,gr-0.4)*2.0, ins, gr
def makes_liability(seq,pos,aa):
    s=seq[:pos]+aa+seq[pos+1:]
    # new NG/NS/NT or >=3 consecutive K/R around pos
    lo=max(0,pos-2); hi=min(len(s),pos+3)
    if re.search(r"N[GST]", s[lo:hi]): return True
    if re.search(r"[KR]{3,}", s[lo:hi]): return True
    return False
rows=[]
for r in csv.DictReader(open(f"{OUT}/manifest.tsv"),delimiter="\t"):
    did=r["design_id"]; seq=r["final_seq"]
    fixed=set(int(x[1:])-1 for x in r["fixed_residues"].split() if x)
    free=[i for i in range(len(seq)) if i not in fixed]
    edits=[]
    cur=seq
    for _ in range(NMAX):
        p0,ins0,gr0=penalty(cur)
        if p0<=0: break
        best=None; bestp=p0
        for pos in free:
            if pos in [e[0] for e in edits]: continue
            for aa in PANEL:
                if aa==cur[pos]: continue
                if makes_liability(cur,pos,aa): continue
                s=cur[:pos]+aa+cur[pos+1:]
                p,_,_=penalty(s)
                if p<bestp: bestp=p; best=(pos,aa,s)
        if best is None: break
        edits.append((best[0],best[1])); cur=best[2]
    pf,insf,grf=penalty(cur)
    ok=pf<=0
    nc=cur.count("C")
    rows.append((did,seq,cur,len(edits),round(insf,1),round(grf,2),ok and nc==0))
    print(f"{did}: edits={len(edits)} -> instab {ProteinAnalysis(seq).instability_index():.1f}->{insf:.1f} gravy {ProteinAnalysis(seq).gravy():.2f}->{grf:.2f} {'PASS' if (ok and nc==0) else 'STILL-FAIL'}  changes={edits}")
with open(f"{OUT}/minedit.fasta","w") as fh, open(f"{OUT}/minedit.csv","w") as fc:
    fc.write("name,sequence\n")
    for did,orig,new,ne,ins,gr,ok in rows:
        if ok and ne>0:
            fh.write(f">{did}_me n_edit={ne} instab={ins} gravy={gr}\n{new}\n")
            fc.write(f"{did}_me,{new}\n")
print(f"\nwrote minedit.fasta/.csv ({sum(1 for r in rows if r[6] and r[3]>0)} passing)")
