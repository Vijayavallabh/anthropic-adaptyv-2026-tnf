#!/usr/bin/env python3
"""Spatial Aggregation Propensity (Chennamsetty 2009, PNAS) — structure-based exposed-hydrophobic-patch
score that GRAVY (a scalar mean) misses. Per residue: sum over neighbors within R of relSASA*(hyd-mean).
SAP_tot+ = sum of positive residue SAP (aggregation-prone exposed hydrophobic surface). Higher = worse.
Usage: sap_score.py <pdb> [chain=A] [R=10]"""
import sys, numpy as np, freesasa
from collections import defaultdict
HYD={'ALA':0.616,'ARG':0.000,'ASN':0.236,'ASP':0.028,'CYS':0.680,'GLN':0.251,'GLU':0.043,'GLY':0.501,
     'HIS':0.165,'ILE':0.943,'LEU':0.943,'LYS':0.283,'MET':0.738,'PHE':1.000,'PRO':0.711,'SER':0.359,
     'THR':0.450,'TRP':0.878,'TYR':0.880,'VAL':0.825}
MEAN=sum(HYD.values())/20.0
pdb=sys.argv[1]; chain=sys.argv[2] if len(sys.argv)>2 else 'A'; R=float(sys.argv[3]) if len(sys.argv)>3 else 10.0
freesasa.setVerbosity(freesasa.silent)
st=freesasa.Structure(pdb); res=freesasa.calc(st); areas=res.residueAreas()
cb={}
for ln in open(pdb):
    if ln[:4]=="ATOM" and ln[21]==chain and ln[12:16].strip() in ("CB","CA"):
        r=int(ln[22:26]); nm=ln[12:16].strip()
        if r not in cb or nm=="CB":
            cb[r]=(np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])]), ln[17:20].strip())
def rel(r):
    try: return min(areas[chain][str(r)].relativeTotal,1.0)
    except Exception: return 0.5
sap={}
for r,(c,rn) in cb.items():
    s=0.0
    for r2,(c2,rn2) in cb.items():
        if np.linalg.norm(c-c2)<=R: s+=rel(r2)*(HYD.get(rn2,0.4)-MEAN)
    sap[r]=s
tot=sum(v for v in sap.values() if v>0)
hot=sorted(sap.items(),key=lambda kv:-kv[1])[:5]
print(f"{pdb.split('/')[-1][:46]:46s} SAP+={tot:6.1f}  hotspots={[(r,round(v,1)) for r,v in hot]}")
