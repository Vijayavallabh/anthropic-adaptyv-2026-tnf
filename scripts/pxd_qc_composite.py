#!/usr/bin/env python3
"""Full composite table for final20_v3 across ALL judging axes, to drive the v4 rebalance.
Axes: pH-switch (#1) | cross-react (#2) | binding (#3) | expression (dominant risk) | novelty gate | lineage diversity.
Reads: final20_v3.csv, uniform re-fold Boltz (human+mouse) ipSAE+iptm, expr_audit flags, ESM-C PLL."""
import csv, os, glob, json, re, subprocess
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
rows=list(csv.DictReader(open("runs/push/final20_v3.csv")))

def boltz(did, species):   # ipSAE_min + iptm from uniform re-fold
    dirs=[f"runs/boltz_rf_{species}{i}/out/by_seed" for i in range(4)]
    for gd in dirs:
        cif=f"{R}/{gd}/{did}/s0/{did}_model_0.cif"
        if not os.path.exists(cif): continue
        pae=os.path.join(os.path.dirname(cif),f"pae_{did}_model_0.npz")
        out=cif.replace(".cif","_10_10.txt")
        if not os.path.exists(out) and os.path.exists(pae):
            subprocess.run([PY,f"{R}/tools/ipsae.py",pae,cif,"10","10"],capture_output=True,cwd=os.path.dirname(cif))
        ipsae=None
        if os.path.exists(out):
            asym=[]
            for ln in open(out):
                f=ln.split()
                if len(f)>6 and f[4]=="asym":
                    try: asym.append(float(f[5]))
                    except: pass
            ipsae=min(asym) if asym else None
        conf=glob.glob(os.path.join(os.path.dirname(cif),f"confidence_{did}_model_0.json"))
        iptm=None
        if conf:
            d=json.load(open(conf[0])); iptm=d.get("iptm")
        return ipsae, iptm
    return None, None

# expr flags
HYD=set("FILVWYMAC")
def maxrun(s,cs):
    b=c=0
    for ch in s:
        if ch in cs: c+=1; b=max(b,c)
        else: c=0
    return b
from Bio.SeqUtils.ProtParam import ProteinAnalysis
def expr(seq):
    pa=ProteinAnalysis(seq); fl=[]
    ins=pa.instability_index(); gr=pa.gravy()
    if ins>=40: fl.append(f"INSTAB{ins:.0f}")
    if gr>=0.4: fl.append(f"GRAVY{gr:.2f}")
    if seq.count("C")%2: fl.append("ODDCYS")
    hr=maxrun(seq,HYD); kr=maxrun(seq,set("KR")); dm=len(re.findall(r"N[GST]",seq))
    if hr>=6: fl.append(f"HYD{hr}")
    if kr>=3: fl.append(f"KR{kr}")
    if dm: fl.append(f"DEAM{dm}")
    return ins,gr,hr,fl

# ESM-C PLL
pll={}
if os.path.exists("runs/push/esmc_pll_v3.tsv"):
    for ln in open("runs/push/esmc_pll_v3.tsv"):
        p=ln.split("\t")
        if len(p)==2 and p[0]!="name":
            try: pll[p[0]]=float(p[1])
            except: pass

def lineage(did):
    for k in ["b1m9","b1m12","b1m7","b4m2","b2m0","TNFpH","tnf_K90"]:
        if did.startswith(k): return "cb1" if k=="TNFpH" else k
    return did.split("_")[0]

recs=[]
for r in rows:
    did=r["design_id"]; seq=r["sequence"]
    hi,ht=boltz(did,"h"); mi,mt=boltz(did,"m")
    ins,gr,hr,fl=expr(seq)
    delta=(hi-mi) if (hi is not None and mi is not None) else None
    recs.append(dict(did=did, lin=lineage(did), fold=r["fold"], nov=r["novelty"],
        fc=r["switch_h_FC"], drep=r["md_d_repel"], caus=r["md_his_causal"],
        hi=hi, mi=mi, delta=delta, ht=ht, mt=mt, ins=round(ins,0), gr=round(gr,2),
        hyd=hr, fl=",".join(fl) or "clean", pll=pll.get(did), role=r["role"][:1]))

# print master table
print(f"{'design':22}{'lin':6}{'nov':4}{'FC':>4}{'drep':>6}{'h_ipS':>7}{'m_ipS':>7}{'Δhm':>6}{'h_tm':>6}{'m_tm':>6}{'ins':>4}{'hyd':>4}{'PLL':>7}  flags")
# rank by a transparent composite: binding gate first, then switch, then cross-react, then expr
def comp(x):
    hi=x['hi'] or 0; mi=x['mi'] or 0
    return (hi>=0.5, hi+mi)  # gate then magnitude, used only for display order
for x in sorted(recs,key=lambda x:-( (x['hi'] or 0)+(x['mi'] or 0) )):
    pll_s=f"{x['pll']:.3f}" if x['pll'] is not None else "-"
    hi=f"{x['hi']:.3f}" if x['hi'] is not None else "-"
    mi=f"{x['mi']:.3f}" if x['mi'] is not None else "-"
    dl=f"{x['delta']:+.2f}" if x['delta'] is not None else "-"
    ht=f"{x['ht']:.2f}" if x['ht'] is not None else "-"
    mt=f"{x['mt']:.2f}" if x['mt'] is not None else "-"
    print(f"{x['did'][:21]:22}{x['lin']:6}{x['nov']:4}{x['fc']:>4}{x['drep']:>6}{hi:>7}{mi:>7}{dl:>6}{ht:>6}{mt:>6}{x['ins']:>4.0f}{x['hyd']:>4}{pll_s:>7}  {x['fl']}")

# lineage counts
from collections import Counter
print("\nlineage counts:", dict(Counter(x['lin'] for x in recs)))
# gate summary
print(f"h_ipSAE>=0.61: {sum(1 for x in recs if (x['hi'] or 0)>=0.61)}/20   "
      f"h>=0.5: {sum(1 for x in recs if (x['hi'] or 0)>=0.5)}/20")
print(f"cross-react (Δhm<=0.15 AND m_ipSAE>=0.45): {sum(1 for x in recs if x['delta'] is not None and x['delta']<=0.15 and (x['mi'] or 0)>=0.45)}/20")
print(f"expr clean: {sum(1 for x in recs if x['fl']=='clean')}/20   DEAM-only: {sum(1 for x in recs if x['fl'].startswith('DEAM') and ',' not in x['fl'])}/20")
