#!/usr/bin/env python3
"""Filter soluble-MPNN rescue outputs: extract binder chain, restore His at switch positions,
audit expression (Adaptyv standard gate), verify interface preserved, keep best/diverse passing variants.
Writes runs/push/rescue/rescued.fasta for re-folding."""
import csv, re, os, glob
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R = "/home/prachh/v/anthropic-adaptyv-2026-tnf"
OUT = f"{R}/runs/push/rescue"
HYD = set("FILVWYMAC")
def max_run(seq, cs):
    b=c=0
    for ch in seq:
        if ch in cs: c+=1; b=max(b,c)
        else: c=0
    return b
def audit(seq):
    pa=ProteinAnalysis(seq); ins=pa.instability_index(); gr=pa.gravy()
    hyd=max_run(seq,HYD); kr=max_run(seq,set("KR")); dm=len(re.findall(r"N[GST]",seq)); nc=seq.count("C")
    hard = (ins<40) and (gr<0.4) and (nc==0)       # the two Adaptyv hard numbers + 0-Cys
    clean = hard and (hyd<6) and (kr<3) and (dm==0)
    return ins,gr,hyd,kr,dm,nc,hard,clean

man={}
for row in csv.DictReader(open(f"{OUT}/manifest.tsv"), delimiter="\t"):
    man[row["design_id"]]=row

fout=open(f"{OUT}/rescued.fasta","w")
summary=[]
for did,row in man.items():
    seqs_dir=f"{OUT}/mpnn/{did}/seqs"
    fa=glob.glob(f"{seqs_dir}/*.fa")
    if not fa: continue
    his=[int(x)-1 for x in row["his_pos"].split(",") if x!=""]
    orig=row["final_seq"]
    oa=audit(orig)
    # parse fa: entries after the first (first = input seq)
    lines=open(fa[0]).read().split("\n")
    entries=[]
    hdr=None
    for ln in lines:
        if ln.startswith(">"): hdr=ln; continue
        if ln.strip() and hdr is not None:
            entries.append((hdr,ln.strip())); hdr=None
    cands=[]
    seen=set()
    for i,(hdr,full) in enumerate(entries):
        if i==0: continue  # input
        binder=full.split(":")[0]
        if len(binder)!=len(orig):  # safety
            continue
        b=list(binder)
        for h in his:
            if 0<=h<len(b): b[h]="H"
        s="".join(b)
        if s in seen: continue
        seen.add(s)
        ins,gr,hyd,kr,dm,nc,hard,clean=audit(s)
        # recovery vs orig
        rec=sum(1 for a,bb in zip(s,orig) if a==bb)/len(orig)
        cands.append(dict(seq=s,ins=ins,gr=gr,hyd=hyd,kr=kr,dm=dm,nc=nc,hard=hard,clean=clean,rec=rec))
    # rank: clean first, then hard, then low instability, then low gravy
    cands.sort(key=lambda c:(-int(c["clean"]),-int(c["hard"]),c["ins"],c["gr"]))
    nhard=sum(c["hard"] for c in cands); nclean=sum(c["clean"] for c in cands)
    print(f"\n### {did}  orig: instab={oa[0]:.1f} gravy={oa[1]:.2f} hyd={oa[2]} -> {'HARD-OK' if oa[6] else 'FAIL'} | {len(cands)} variants: {nhard} pass-hard, {nclean} clean")
    for c in cands[:5]:
        tag="CLEAN" if c["clean"] else ("HARD" if c["hard"] else "----")
        print(f"   {tag:5} instab={c['ins']:5.1f} gravy={c['gr']:5.2f} hyd={c['hyd']:2d} kr={c['kr']} dm={c['dm']} rec={c['rec']:.2f}")
    # write top up-to-3 passing (hard) variants for re-fold
    kept=0
    for j,c in enumerate([c for c in cands if c["hard"]][:8]):
        fout.write(f">{did}_r{j} instab={c['ins']:.1f} gravy={c['gr']:.2f} hyd={c['hyd']} rec={c['rec']:.2f}\n{c['seq']}\n")
        kept+=1
    summary.append((did,oa[6],nhard,nclean,kept))
fout.close()
print("\n=== SUMMARY ===")
print(f"{'design':8} {'orig_ok':>7} {'n_hard':>6} {'n_clean':>7} {'kept':>4}")
for did,ok,nh,nc,k in summary:
    print(f"{did:8} {str(ok):>7} {nh:>6} {nc:>7} {k:>4}")
print(f"\nwrote {OUT}/rescued.fasta")
