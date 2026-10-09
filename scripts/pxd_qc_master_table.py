#!/usr/bin/env python3
"""Integrated master table for final20_v3 — ALL axes, ipSAE DEMOTED to a weak tiebreaker
(per EGFR post-mortem: confidence score poorly predicts lab binding).
Primary axes: EXPRESSION (instab/SAP/liabilities) | SWITCH (FC/MD) | CROSS-REACT (Δhm) | INTERFACE (dSASA/contacts/Hbond).
Writes runs/push/master_v3.csv; prints a gated tier ranking."""
import csv, os, glob, json, re, subprocess, tempfile, numpy as np, freesasa
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
freesasa.setVerbosity(freesasa.silent)
PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
rows=list(csv.DictReader(open("runs/push/final20_v3.csv")))

def boltz_ipsae(did, sp):
    for i in range(4):
        cif=f"{R}/runs/boltz_rf_{sp}{i}/out/by_seed/{did}/s0/{did}_model_0.cif"
        if not os.path.exists(cif): continue
        out=cif.replace(".cif","_10_10.txt")
        if os.path.exists(out):
            asym=[float(l.split()[5]) for l in open(out) if len(l.split())>6 and l.split()[4]=="asym"]
            return min(asym) if asym else None
    return None

HYD=set("FILVWYMAC")
def maxrun(s,cs):
    b=c=0
    for ch in s:
        c=c+1 if ch in cs else 0; b=max(b,c)
    return b
from Bio.SeqUtils.ProtParam import ProteinAnalysis
def expr(seq):
    pa=ProteinAnalysis(seq); ins=pa.instability_index(); gr=pa.gravy()
    hr=maxrun(seq,HYD); dm=len(re.findall(r"N[GST]",seq))
    # instability driver: count Pro-dipeptides (ProtParam penalizes PP/xP heavily) -> false-positive check
    npro=seq.count("P")
    return round(ins,0),round(gr,2),hr,dm,npro

def atoms(pdb,ch):
    o=[]
    for ln in open(pdb):
        if ln[:4]=="ATOM" and ln[21]==ch:
            el=ln[76:78].strip() or ln[12:16].strip()[0]
            o.append((int(ln[22:26]),el,np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])))
    return o
def sasa(pdb,keep):
    t=tempfile.NamedTemporaryFile('w',suffix='.pdb',delete=False)
    for ln in open(pdb):
        if ln[:4]=="ATOM" and ln[21] in keep: t.write(ln)
    t.write("END\n"); t.close()
    try: v=freesasa.calc(freesasa.Structure(t.name)).totalArea()
    except Exception: v=float('nan')
    os.unlink(t.name); return v
def iface(pdb):
    A=atoms(pdb,'A'); B=atoms(pdb,'B')
    if not A or not B: return None
    cA=np.array([a[2] for a in A]); cB=np.array([b[2] for b in B])
    D=np.linalg.norm(cA[:,None,:]-cB[None,:,:],axis=2)
    cp=set(); ir=set()
    for i,j in zip(*np.where(D<5.0)): cp.add((A[i][0],B[j][0])); ir.add(A[i][0])
    hb=sum(1 for i,j in zip(*np.where(D<3.5)) if A[i][1] in {'N','O'} and B[j][1] in {'N','O'})
    sA,sB,sAB=sasa(pdb,'A'),sasa(pdb,'B'),sasa(pdb,'AB')
    ds=(sA+sB-sAB) if not any(np.isnan([sA,sB,sAB])) else float('nan')
    return len(cp),hb,len(ir),round(ds)
def sap(pdb):
    o=subprocess.run([f"{R}/.venv/bin/python","scripts/sap_score.py",pdb,"A","10"],capture_output=True,text=True).stdout
    m=re.search(r"SAP\+=\s*([\d.]+)",o); return float(m.group(1)) if m else None

def lineage(did):
    for k in ["b1m9","b1m12","b1m7","b4m2","b2m0","TNFpH","tnf_K90"]:
        if did.startswith(k): return "cb1" if k=="TNFpH" else k
    return did.split("_")[0]

recs=[]
for r in rows:
    did=r["design_id"]; seq=r["sequence"]; pdb=f"runs/push/refold_pdb/{did}.pdb"
    hi=boltz_ipsae(did,"h"); mi=boltz_ipsae(did,"m")
    ins,gr,hr,dm,npro=expr(seq)
    nc,hb,ir,ds=iface(pdb) if os.path.exists(pdb) else (None,)*4
    sp=sap(pdb) if os.path.exists(pdb) else None
    try: fc=float(r["switch_h_FC"])
    except: fc=0
    try: drep=float(r["md_d_repel"])
    except: drep=0
    recs.append(dict(did=did,lin=lineage(did),nov=r["nov"] if "nov" in r else r["novelty"],
        fc=fc,drep=drep,caus=r["md_his_causal"],hi=hi,mi=mi,
        delta=(hi-mi) if (hi and mi) else None, ins=ins,gr=gr,hyd=hr,dm=dm,npro=npro,
        sap=sp,ds=ds,nc=nc,hb=hb,seq=seq,role=r["role"]))

# save master csv
with open("runs/push/master_v3.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["design_id","lineage","novelty","switch_FC","md_drep","h_ipSAE","m_ipSAE","delta_hm","instab","gravy","hydrun","deam","npro","SAP","dSASA","contacts","hbond"])
    for x in recs: w.writerow([x['did'],x['lin'],x['nov'],x['fc'],x['drep'],x['hi'],x['mi'],
        f"{x['delta']:.3f}" if x['delta'] is not None else "",x['ins'],x['gr'],x['hyd'],x['dm'],x['npro'],x['sap'],x['ds'],x['nc'],x['hb']])

# gated tiers — expression-first, interface-real, ipSAE minor
def expr_ok(x): return x['ins']<45 and (x['sap'] or 0)<2.0
def iface_ok(x): return (x['ds'] or 0)>=1900 or (x['nc'] or 0)>=60
def switch_strong(x): return x['fc']>=10 or x['drep']>=1.5
print(f"{'design':22}{'lin':6}{'nov':4}{'FC':>4}{'drep':>5}{'Δhm':>6}{'dSASA':>6}{'ct':>4}{'hb':>3}{'ins':>4}{'SAP':>4}{'hyd':>3}{'hi':>5}  verdict")
def score(x):  # display ranking: expression + switch + cross-react + interface (ipSAE weak)
    s=0
    s+= 2.0 if x['ins']<40 else (1.0 if x['ins']<45 else -1.0)
    s+= min(x['fc']/20,3.0) + min(x['drep']/2,1.5)
    s+= (1.5 if (x['delta'] is not None and x['delta']<=0.1) else (0.5 if (x['delta'] is not None and x['delta']<=0.2) else 0))
    s+= min((x['ds'] or 0)/1500,2.0) + min((x['nc'] or 0)/50,1.5)
    s+= 0.3*((x['hi'] or 0)+(x['mi'] or 0))   # ipSAE small weight
    s+= 1.0 if x['nov']=="4/4" else 0
    s-= 1.0 if (x['sap'] or 0)>=2 else 0
    return s
for x in sorted(recs,key=lambda x:-score(x)):
    v=[]
    if not expr_ok(x): v.append(f"EXPR!(ins{x['ins']:.0f})")
    if not iface_ok(x): v.append("weak-iface")
    if switch_strong(x): v.append("SWITCH+")
    if x['delta'] is not None and x['delta']<=0.1: v.append("XR+")
    if x['nov']=="4/4": v.append("nov4")
    dl=f"{x['delta']:+.2f}" if x['delta'] is not None else "-"
    hi=f"{x['hi']:.2f}" if x['hi'] is not None else "-"
    print(f"{x['did'][:21]:22}{x['lin']:6}{x['nov']:4}{x['fc']:>4.0f}{x['drep']:>5.1f}{dl:>6}{x['ds'] or 0:>6.0f}{x['nc'] or 0:>4}{x['hb'] or 0:>3}{x['ins']:>4.0f}{x['sap'] or 0:>4.1f}{x['hyd']:>3}{hi:>5}  {' '.join(v)}")
print("\nwrote runs/push/master_v3.csv")
