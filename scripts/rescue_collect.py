#!/usr/bin/env python3
"""Accumulate rescue winners: variants passing expression (instab<40,GRAVY<0.4,0-Cys) AND ipSAE_wc>THR.
Reads current rescued.fasta + the just-written ipSAE _10_10.txt in boltz_rescue_{h,m}. Appends to winners.tsv (dedup)."""
import os, csv, re, glob
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; OUT=f"{R}/runs/push/rescue"
THR=float(os.environ.get("WC_THR","0.5"))
HYD=set("FILVWYMAC")
def mrun(s,cs):
    b=c=0
    for ch in s:
        c=c+1 if ch in cs else 0; b=max(b,c)
    return b
def ipsae(sp,name):
    t=f"{R}/runs/boltz_rescue_{sp}/out/by_seed/{name}/s0/{name}_model_0_10_10.txt"
    if not os.path.exists(t): return None
    for ln in open(t):
        p=ln.split()
        if len(p)>5 and p[4]=="max":
            try: return float(p[5])
            except: return None
    return None
# seqs from current rescued.fasta
seqs={}
name=None
for ln in open(f"{OUT}/rescued.fasta"):
    if ln.startswith(">"): name=ln[1:].split()[0]
    elif ln.strip() and name: seqs[name]=ln.strip(); name=None
rows=[]
for nm,seq in seqs.items():
    pa=ProteinAnalysis(seq); ins=pa.instability_index(); gr=pa.gravy(); nc=seq.count("C")
    hard=(ins<40) and (gr<0.4) and (nc==0)
    h=ipsae("h",nm); m=ipsae("m",nm)
    if h is None or m is None: continue
    wc=min(h,m)
    design=nm.rsplit("_r",1)[0]
    rows.append((design,nm,seq,round(ins,1),round(gr,2),round(h,3),round(m,3),round(wc,3),hard and wc>THR))
# load existing winners
wf=f"{OUT}/winners.tsv"; existing={}
if os.path.exists(wf):
    for r in csv.DictReader(open(wf),delimiter="\t"): existing[r["variant"]]=r
# append winners
added=0
for design,nm,seq,ins,gr,h,m,wc,win in rows:
    if win and nm not in existing:
        existing[nm]=dict(design=design,variant=nm,instab=ins,gravy=gr,ipsae_h=h,ipsae_m=m,ipsae_wc=wc,seq=seq); added+=1
cols=["design","variant","instab","gravy","ipsae_h","ipsae_m","ipsae_wc","seq"]
with open(wf,"w") as fh:
    fh.write("\t".join(cols)+"\n")
    for v in sorted(existing.values(),key=lambda r:(r["design"],-float(r["ipsae_wc"]))):
        fh.write("\t".join(str(v[c]) for c in cols)+"\n")
print(f"this round: {sum(1 for r in rows if r[-1])} winners (of {len(rows)}); +{added} new; total winners={len(existing)}")
print(f"{'design':8} {'variant':10} {'instab':>6} {'gravy':>6} {'ipH':>5} {'ipM':>5} {'ipWC':>5} win")
for design,nm,seq,ins,gr,h,m,wc,win in sorted(rows,key=lambda x:(x[0],-x[7])):
    print(f"{design:8} {nm:10} {ins:6.1f} {gr:6.2f} {h:5.2f} {m:5.2f} {wc:5.2f} {'WIN' if win else ''}")
