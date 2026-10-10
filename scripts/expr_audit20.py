#!/usr/bin/env python3
"""Complete expression audit over the 20 v9 members. Primary axis = ESMFold monomer pLDDT (cached).
Orthogonal axes: SAP+ (structure, Chennamsetty), GRAVY (Kyte-Doolittle), acidic% (D+E), aromatic% (F+W+Y),
netQ (K+R-D-E). Consensus-drop = weak on >=2 axes. Prints a table + any flags."""
import os, glob, csv, numpy as np, freesasa
R0="/home/prachh/v/anthropic-adaptyv-2026-tnf"
freesasa.setVerbosity(freesasa.silent)
HYD={'ALA':0.616,'ARG':0.000,'ASN':0.236,'ASP':0.028,'CYS':0.680,'GLN':0.251,'GLU':0.043,'GLY':0.501,
     'HIS':0.165,'ILE':0.943,'LEU':0.943,'LYS':0.283,'MET':0.738,'PHE':1.000,'PRO':0.711,'SER':0.359,
     'THR':0.450,'TRP':0.878,'TYR':0.880,'VAL':0.825}
MEAN=sum(HYD.values())/20.0
KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,
    'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
# cached ESMFold monomer pLDDT (primary expression axis)
PLDDT={'des_3050_6_c39':0.88,'des_4098_1_b14':0.90,'des_4098_1_c44':0.82,'des_180_6_H9-14-38':0.81,
 'des_818_4_H16-74':0.88,'des_3050_6_H52-54-79':0.79,'des_875_5_H3-57-80':0.86,'des_287_5_H37':0.88,
 'des_2106_2_H43':0.88,'des_59_1_H51':0.90,'TNFpH06':0.85,'TNFpH06_sh3':0.86,'b1m12_A41H_vsR27':0.89,
 'b1m12_2A67H_K112_r1':0.88,'b2m0_c70':0.90,'b1m7_2A57H_K183_r1':0.93,'b4m2_r1':0.94,'b4m2_A57H':0.93,
 'tnf_K90_tnf_K90_2_model_3_3_sh4':0.92,'b1m9_A12H':0.87}
def find_pdb(name):
    for d in ("runs/wave4b/monomer","runs/wave3/monomer","runs/novelty_struct/binders","runs/wave4b/triple_pdb"):
        p=os.path.join(R0,d,name+".pdb")
        if os.path.exists(p): return p
    return None
def sap_plus(pdb,chain='A',R=10.0):
    st=freesasa.Structure(pdb); res=freesasa.calc(st); areas=res.residueAreas()
    cb={}
    for ln in open(pdb):
        if ln[:4]=="ATOM" and ln[12:16].strip() in ("CB","CA"):
            ch=ln[21]; r=int(ln[22:26]); nm=ln[12:16].strip()
            if (r,ch) not in cb or nm=="CB": cb[(r,ch)]=(np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])]),ln[17:20].strip(),ch)
    def rel(r,ch):
        try: return min(areas[ch][str(r)].relativeTotal,1.0)
        except Exception: return 0.5
    tot=0.0
    for (r,ch),(c,rn,_) in cb.items():
        s=sum(rel(r2,ch2)*(HYD.get(rn2,0.4)-MEAN) for (r2,ch2),(c2,rn2,_) in cb.items() if np.linalg.norm(c-c2)<=R)
        if s>0: tot+=s
    return tot
rows=[]
with open(os.path.join(R0,"runs/push/submission_v9.csv")) as f:
    for name,seq,_ in list(csv.reader(f))[1:]:
        L=len(seq); gravy=sum(KD[a] for a in seq)/L
        acidic=100*sum(seq.count(x) for x in "DE")/L; arom=100*sum(seq.count(x) for x in "FWY")/L
        netQ=sum(seq.count(x) for x in "KR")-sum(seq.count(x) for x in "DE")
        pdb=find_pdb(name); sap=sap_plus(pdb) if pdb else None
        pl=PLDDT.get(name)
        rows.append((name,L,pl,sap,gravy,acidic,arom,netQ,pdb is not None))
print(f"{'name':32s} {'L':>3} {'pLDDT':>5} {'SAP+':>6} {'GRAVY':>6} {'acid%':>5} {'arom%':>5} {'netQ':>4}  flags")
drops=[]
saps=[r[3] for r in rows if r[3] is not None]
sap_hi=np.percentile(saps,90) if saps else 999
for name,L,pl,sap,gravy,acidic,arom,netQ,havepdb in rows:
    fl=[]
    if pl is not None and pl<0.75: fl.append("pLDDT")      # primary axis weak
    if sap is not None and sap>sap_hi: fl.append("SAP")     # top-decile aggregation
    if gravy>0.0: fl.append("GRAVY")                        # net hydrophobic
    if acidic<5.0: fl.append("acidic")                      # too few acidic (solubility)
    if arom>12.0: fl.append("arom")                         # aromatic-rich (aggregation)
    drop = len(fl)>=2
    if drop: drops.append(name)
    sp=f"{sap:6.1f}" if sap is not None else "   n/a"
    print(f"{name:32s} {L:3d} {pl if pl else 0:5.2f} {sp} {gravy:6.2f} {acidic:5.1f} {arom:5.1f} {netQ:4d}  {'DROP:' if drop else ''}{','.join(fl)}")
print(f"\nSAP+ 90th pctile (flag threshold) = {sap_hi:.1f}")
print("CONSENSUS-DROPS (weak on >=2 axes):", drops if drops else "NONE — all 20 pass expression consensus")
missing=[r[0] for r in rows if not r[8]]
if missing: print("missing monomer PDB (SAP skipped):", missing)
