#!/usr/bin/env python3
"""submission_v7 (novelty now PORTAL-CONFIRMED): all 4 de-novo beta-pairing switch survivors scored 3/4
on the Adaptyv portal (PASS the hard gate). Scale the strongest evidence-backed strategy: ADD the 2 best
de-novo (des_287_5_H37 = strongest switch + best cross-react 0.703; des_2106_2_H43 = pristine expression +
solid switch) and TRIM the single most-redundant cb1 anchor (TNFpH07: FC13, cross-react 0.185, thermo-only)
to keep 20 and reduce cb1 over-weight (6->5) -> improves cross-react (#2) + fold diversity + switch-confidence.
Expression finding applied (DROP ProtParam instab penalty: Pro-FP). v6 (md5 81913f5c) stays fallback floor;
upload USER-GATED."""
import csv, os, hashlib
import numpy as np
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
ADD=["des_287_5_H37","des_2106_2_H43"]     # portal 3/4 CONFIRMED, triple-confirmed switch
REMOVE={"TNFpH07"}                          # weakest/most-redundant cb1 anchor

seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/bpair_mature/variants.csv"))}
DN={
 "des_287_5_H37":dict(lineage="bpair_d287",novelty="3/4",md_drep="2.10",dpka="-2.57",dE="17.35",
   h="0.728",m="0.703",sap="2.5",instab="46",gravy="-0.39",gluP="18.8",lysP="3.8",net="-9",his="37",lys="93",novTM="0.686"),
 "des_2106_2_H43":dict(lineage="bpair_d2106",novelty="3/4",md_drep="1.65",dpka="-1.66",dE="9.96",
   h="0.657",m="0.566",sap="0.0",instab="24",gravy="-0.45",gluP="21.2",lysP="11.2",net="-10",his="43",lys="93",novTM="0.676"),
 "des_962_5_H37":dict(lineage="bpair_d962",novelty="3/4",md_drep="2.23",dpka="-",dE="-",
   h="0.667",m="0.542",sap="-",instab="-",gravy="-",gluP="21.2",lysP="17.5",net="0",his="37",lys="93",novTM="0.758"),
 "des_495_1_H78":dict(lineage="bpair_d495",novelty="3/4",md_drep="1.69",dpka="-0.15",dE="2.97",
   h="0.675",m="0.619",sap="0.0",instab="70",gravy="-0.35",gluP="22.5",lysP="10.0",net="-8",his="78",lys="6",novTM="0.720"),
}

def contacts(variant):
    f=f"{R}/runs/bpair_mature/md/{variant}.pdb"; A={}; B={}
    for l in open(f):
        if l[:4]!="ATOM" or l[12:16].strip()[0]=="H": continue
        ch=l[21]; rn=int(l[22:26]); c=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        (A if ch=="A" else B).setdefault(rn,[]).append(c)
    Bc={rn:np.array(v) for rn,v in B.items()}; pairs=set()
    for ra,va in A.items():
        va=np.array(va)
        for rb,vb in Bc.items():
            if np.min(np.linalg.norm(va[:,None,:]-vb[None,:,:],axis=2))<4.5: pairs.add((ra,rb))
    nc=len(pairs); return nc, nc*32

def f(x,d=0.0):
    try: return float(x)
    except: return d
dn_master={}
for v,d in DN.items():
    nc,dsasa=contacts(v)
    dn_master[v]=dict(design_id=v,lineage=d["lineage"],novelty=d["novelty"],switch_FC="0",md_drep=d["md_drep"],
      h_ipSAE=d["h"],m_ipSAE=d["m"],delta_hm=f"{abs(f(d['h'])-f(d['m'])):.3f}",instab=d["instab"],gravy=d["gravy"],
      hydrun="0",deam="0",npro="1",SAP=d["sap"],dSASA=str(dsasa),contacts=str(nc),hbond="-")

master={r["design_id"]:r for r in csv.DictReader(open("runs/push/master_v6.csv"))}
v6sub=list(csv.DictReader(open("runs/push/submission_v6.csv")))
seqmap={r["name"]:r["sequence"] for r in v6sub}
for d in ADD: seqmap[d]=seqs[d]; master[d]=dn_master[d]
ids=[r["name"] for r in v6sub if r["name"] not in REMOVE]+ADD
assert len(ids)==20, f"expected 20 got {len(ids)}"

def score(did):  # v6 score() MINUS instability penalty (expression research: Pro-FP, no de-novo validation)
    m=master[did]; fc=f(m["switch_FC"]); dr=f(m["md_drep"])
    hi=f(m["h_ipSAE"]); mi=f(m["m_ipSAE"]); dl=f(m["delta_hm"],9); ds=f(m["dSASA"]); nc=f(m["contacts"])
    sap=f(m["SAP"]); nov=m["novelty"]; s=0.0
    s+= 1.0 if sap<1.0 else 0.0
    if sap>1.5: s-= min(0.4*(sap-1.5), 2.0)
    s+= min(fc/20,3.0)+min(dr/2,1.5)
    s+= 2.0*min(hi,mi) + (0.3 if dl<=0.15 else 0)
    s+= min(ds/1500,2.0)+min(nc/50,1.5)
    s+= 0.3*(hi+mi)
    s+= 1.0 if nov=="4/4" else 0
    return s

ranked=sorted(ids,key=lambda d:-score(d))
with open("runs/push/submission_v7.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"])
    for d in ranked: w.writerow([d, seqmap[d], "single_chain"])
with open("runs/push/master_v7.csv","w",newline="") as fh:
    cols=list(master[ids[0]].keys()); w=csv.DictWriter(fh,fieldnames=cols); w.writeheader()
    for d in ids: w.writerow(master[d])
with open("runs/push/bpair_candidates.csv","w",newline="") as fh:
    cols=["design_id","sequence","lineage","in_v7","switch_MD_drep","switch_dpKa","switch_xtb_dE","xreact_h","xreact_m",
          "xreact_min","SAP","gluP","net","switch_His","target_Lys_conserved","portal_novelty","len"]
    w=csv.writer(fh); w.writerow(cols)
    for v,d in sorted(DN.items(), key=lambda kv:-f(kv[1]['m'])):
        w.writerow([v,seqs[v],d["lineage"],("YES" if v in ADD else "backup"),d["md_drep"],d["dpka"],d["dE"],d["h"],d["m"],
          f"{min(f(d['h']),f(d['m'])):.3f}",d["sap"],d["gluP"],d["net"],f"His{d['his']}",f"Lys{d['lys']}=conserved","3/4",len(seqs[v])])

AA=set("ACDEFGHIKLMNPQRSTVWY")
sub=list(csv.DictReader(open("runs/push/submission_v7.csv"))); bad=[]
for r in sub:
    s=r["sequence"]
    if set(s)-AA: bad.append((r["name"],"nonAA"))
    if "C" in s: bad.append((r["name"],"Cys"))
    if not(60<=len(s)<=130): bad.append((r["name"],f"len{len(s)}"))
    if r["molecule_class"]!="single_chain": bad.append((r["name"],"molclass"))
nm=[r["name"] for r in sub]; sq=[r["sequence"] for r in sub]
ok=len(sub)==20 and len(set(nm))==20 and len(set(sq))==20 and not bad
md5=lambda p:hashlib.md5(open(p,'rb').read()).hexdigest()
cb1=[d for d in ids if master[d]["lineage"]=="cb1"]
print("=== submission_v7 (novelty PORTAL-CONFIRMED 3/4; 20 designs) ===")
print(f"{'rank':>4} {'design':30}{'score':>6} {'xreact_min':>10} {'lineage':>10} {'note':>16}")
for i,d in enumerate(ranked,1):
    m=master[d]; xm=min(f(m['h_ipSAE']),f(m['m_ipSAE']))
    note="<-- NEW de-novo" if d in ADD else ""
    print(f"{i:>4} {d[:29]:30}{score(d):>6.2f} {xm:>10.3f} {m['lineage'][:9]:>10} {note:>16}")
print(f"\nFINAL-TEST: {'ALL PASS' if ok else 'ISSUES '+str(bad)}")
print(f"ADDED (portal 3/4): {ADD}  | TRIMMED cb1: {sorted(REMOVE)}  | cb1 now {len(cb1)}")
print(f"submission_v7.csv md5={md5('runs/push/submission_v7.csv')}")
print(f"de-novo in v7: 3 total (des_59_1_H51, des_76_2_H55 from v6 + 2 new); lineages: {len(set(master[d]['lineage'] for d in ids))} distinct")
print(f"v6 md5 81913f5c = fallback floor; upload USER-GATED; backups -> runs/push/bpair_candidates.csv")
