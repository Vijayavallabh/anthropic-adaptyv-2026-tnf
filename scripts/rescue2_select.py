#!/usr/bin/env python3
"""Select expression-rescued candidates from SolubleMPNN output: per design, keep surface-redesigned
sequences passing the Adaptyv gate (instability<40, GRAVY<0.4, 0-Cys) that RETAIN the switch His, take
the best-expressing few. Writes runs/push/rescue2/candidates.csv (name,sequence) for re-folding.
Usage: rescue2_select.py [top_n=4]"""
import csv, os, sys, glob
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; OUT=f"{R}/runs/push/rescue2"
TOPN=int(sys.argv[1]) if len(sys.argv)>1 else 4
man={r["design_id"]:r for r in csv.DictReader(open(f"{OUT}/manifest.tsv"),delimiter="\t")}
def ex(seq):
    pa=ProteinAnalysis(seq); return pa.instability_index(), pa.gravy()
rows=[]; stats=[]
for did,m in man.items():
    his_pos=set(int(x) for x in m["his_pos"].split(",") if x)  # 1-based switch His
    fa=glob.glob(f"{OUT}/mpnn/{did}/seqs/*.fa")
    if not fa: print(f"{did}: no fa"); continue
    lines=open(fa[0]).read().splitlines()
    cands=[]
    for i in range(0,len(lines)-1,2):
        if not lines[i].startswith(">"): continue
        if i==0: continue  # native
        seq=lines[i+1].split(":")[0].strip()  # chain A binder
        if not seq: continue
        # His must be retained at all switch positions
        if any(p-1>=len(seq) or seq[p-1]!="H" for p in his_pos): continue
        if seq.count("C")!=0: continue
        ins,gr=ex(seq)
        cands.append((seq,round(ins,1),round(gr,2)))
    # passing the gate
    passing=[c for c in cands if c[1]<40 and c[2]<0.4]
    passing.sort(key=lambda c:(c[1],c[2]))  # lowest instability then GRAVY
    # dedup identical seqs
    seen=set(); uniq=[]
    for c in passing:
        if c[0] in seen: continue
        seen.add(c[0]); uniq.append(c)
    kept=uniq[:TOPN]
    best_all=min(cands,key=lambda c:c[1]) if cands else None
    stats.append((did,len(cands),len(passing),best_all[1] if best_all else None))
    for k,(seq,ins,gr) in enumerate(kept):
        rows.append((f"{did}_r{k+1}",seq,ins,gr))
with open(f"{OUT}/candidates.csv","w") as fh:
    fh.write("name,sequence\n")
    for n,s,_,_ in rows: fh.write(f"{n},{s}\n")
with open(f"{OUT}/candidates_expr.tsv","w") as fh:
    fh.write("name\tinstability\tgravy\n")
    for n,s,ins,gr in rows: fh.write(f"{n}\t{ins}\t{gr}\n")
print(f"{'design':22} {'n_cand':>6} {'n_pass':>6} {'best_instab':>11}")
for did,nc,npass,bi in stats:
    print(f"{did:22} {nc:6} {npass:6} {str(bi):>11} {'-> RESCUABLE' if npass>0 else '(no gate-pass; stubborn)'}")
print(f"\nwrote {len(rows)} rescued candidates -> {OUT}/candidates.csv")
