#!/usr/bin/env python3
"""final20_v5: swap 2 weakest f8-lineage members (b1m9_1A55H_K93_sh6_r1 rank10, b1m9_1A55H_K93_rx rank20;
weak cross-react 0.36/0.26, unconfirmed switch) for 2 TRIPLE-CONFIRMED cross-react+switch winners
(des_59_1_H51, des_76_2_H55). Switch = MD audit_brutal SURVIVE + PROPKA pKa-drop + qm GFN2-xTB dE>0,
all 3 methods sign-agree. Re-rank with the SAME expression-first score() as v4. Upload USER-GATED."""
import csv, os, hashlib
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)

REMOVE={"b1m9_1A55H_K93_sh6_r1","b1m9_1A55H_K93_rx"}

# new designs: master-schema metrics + final20-schema rows
NEW_SEQ={
 "des_59_1_H51":"LAELVKKVAEHLGIKISDEEIKKMIEEAKARVLALMGREPSQQLLEVFVLHDLLDYLRNEKGFSVEQLNEAVKYWLELVEKE",
 "des_76_2_H55":"IVDIILEALRKMGITLSRERLEAMYAEVRETLAALMATPATEELIFVFLLAAVLHYLREEEGLSVEELNEFVKLVLELKEKA"}
NEW_MASTER={
 "des_59_1_H51":dict(design_id="des_59_1_H51",lineage="f8d1",novelty="3/4",switch_FC="4",md_drep="3.26",
   h_ipSAE="0.765",m_ipSAE="0.521",delta_hm="0.244",instab="36.5",gravy="-0.13",hydrun="0",deam="0",
   npro="1",SAP="2.7",dSASA="2128",contacts="66",hbond="15"),
 "des_76_2_H55":dict(design_id="des_76_2_H55",lineage="f8d2",novelty="4/4",switch_FC="2",md_drep="0.65",
   h_ipSAE="0.780",m_ipSAE="0.633",delta_hm="0.147",instab="44.8",gravy="0.42",hydrun="0",deam="0",
   npro="1",SAP="5.7",dSASA="2245",contacts="66",hbond="9")}
NEW_ROW={
 "des_59_1_H51":dict(design_id="des_59_1_H51",fold="f8d1",novelty="3/4",switch_h_FC="4",switch_h_peak="-",
   switch_m_FC="-",pkai_FC="-",md_his_causal="YES",md_d_repel="3.26",md_dud="-4.11",iptm_h="0.916",
   iptm_m="0.89",instability="36.5",gravy="-0.13",expr="PASS",length="82",n_cys="0",sequence=NEW_SEQ["des_59_1_H51"],
   role="cross-react+switch (partial-diffusion b1m9): triple-confirmed single-His switch His51 vs target K88 (MD d_repel 3.26 brutal-SURVIVE + PROPKA pKa-drop 6.7->5.5 + qm GFN2-xTB dE +14.6, 3-method sign-consensus); x-react min-ipSAE 0.521 (iptm 0.92h/0.89m), 3/4-novel (seq 0% id, struct TM 0.885), expr-clean SAP 2.7 gravy -0.13 0-Cys"),
 "des_76_2_H55":dict(design_id="des_76_2_H55",fold="f8d2",novelty="4/4",switch_h_FC="2",switch_h_peak="-",
   switch_m_FC="-",pkai_FC="-",md_his_causal="YES",md_d_repel="0.65",md_dud="-1.52",iptm_h="0.947",
   iptm_m="0.907",instability="44.8",gravy="0.42",expr="FLAG",length="82",n_cys="0",sequence=NEW_SEQ["des_76_2_H55"],
   role="BEST cross-react + switch (partial-diffusion b1m9): triple-confirmed single-His switch His55 vs target K88 (MD d_repel 0.65 brutal-SURVIVE + PROPKA pKa-drop 6.3->5.0 + qm GFN2-xTB dE +18.2, 3-method sign-consensus); BEST x-react min-ipSAE 0.633 (iptm 0.95h/0.91m), 4/4-novel (seq 0% id, struct TM 0.732), expr-FLAG SAP 5.7 gravy +0.42 (hydrophobic, higher wet-lab risk) 0-Cys")}

# build master_v5
master={r["design_id"]:r for r in csv.DictReader(open("runs/push/master_v3.csv"))}
for did in REMOVE: master.pop(did,None)
master.update(NEW_MASTER)

def f(x,d=0.0):
    try: return float(x)
    except: return d
def score(did):
    m=master[did]; ins=f(m["instab"]); fc=f(m["switch_FC"]); dr=f(m["md_drep"])
    hi=f(m["h_ipSAE"]); mi=f(m["m_ipSAE"]); dl=f(m["delta_hm"],9); ds=f(m["dSASA"]); nc=f(m["contacts"])
    sap=f(m["SAP"]); nov=m["novelty"]; npro=int(f(m["npro"]))
    s=0.0
    # expression: graded SAP risk (continuous) instead of a binary step at 2.0 (which over-penalized
    # a design at 2.7 as much as one at 5.7); clean bonus retained for very low SAP.
    s+= 1.0 if sap<1.0 else 0.0
    if sap>1.5: s-= min(0.4*(sap-1.5), 2.0)
    if ins>=45: s-= (0.2 if npro>=3 else 0.5)
    elif ins>=40: s-= 0.1
    s+= min(fc/20,3.0)+min(dr/2,1.5)
    # cross-react (#2): reward WORST-ortholog binding min(h,m) (the judging criterion), not the h-m gap
    # (the old delta_hm term rewarded weak-but-balanced designs over genuinely cross-reactive ones);
    # keep a small balance bonus.
    s+= 2.0*min(hi,mi) + (0.3 if dl<=0.15 else 0)
    s+= min(ds/1500,2.0)+min(nc/50,1.5)
    s+= 0.3*(hi+mi)
    s+= 1.0 if nov=="4/4" else 0
    return s
def expr_verdict(did):
    m=master[did]; ins=f(m["instab"]); npro=int(f(m["npro"])); sap=f(m["SAP"]); dm=int(f(m["deam"])); hyd=int(f(m["hydrun"]))
    notes=[]
    if ins>=45 and npro>=3: notes.append(f"ins{ins:.0f}=Pro-FP({npro}Pro)")
    elif ins>=45: notes.append(f"ins{ins:.0f}-GENUINE")
    elif ins>=40: notes.append(f"ins{ins:.0f}-borderline")
    if hyd>=6 and sap<2: notes.append(f"HYD{hyd}=buried(SAP{sap:.1f})")
    if dm: notes.append(f"deam{dm}-minor")
    if sap>=2: notes.append(f"SAP{sap:.1f}-AGG")
    return "clean" if not notes else ";".join(notes)

# build final20_v5 rows (v3 schema), swap
v3=list(csv.DictReader(open("runs/push/final20_v3.csv")))
base_cols=list(v3[0].keys())
rows=[r for r in v3 if r["design_id"] not in REMOVE]
for did in ("des_59_1_H51","des_76_2_H55"):
    rows.append({c:NEW_ROW[did].get(c,"-") for c in base_cols})
assert len(rows)==20, f"expected 20, got {len(rows)}"

ranked=sorted(rows,key=lambda r:-score(r["design_id"]))
qc_cols=["submit_rank","comp_score","h_ipSAE","m_ipSAE","delta_hm","dSASA","contacts","hbond","SAP","expr_verdict"]
hdr=base_cols+qc_cols
with open("runs/push/final20_v5.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=hdr); w.writeheader()
    for i,r in enumerate(ranked,1):
        did=r["design_id"]; m=master[did]; r=dict(r)
        r.update(submit_rank=i,comp_score=f"{score(did):.2f}",h_ipSAE=m["h_ipSAE"],m_ipSAE=m["m_ipSAE"],
                 delta_hm=m["delta_hm"],dSASA=m["dSASA"],contacts=m["contacts"],hbond=m["hbond"],
                 SAP=m["SAP"],expr_verdict=expr_verdict(did))
        w.writerow(r)
# submission_v5 (Adaptyv format)
with open("runs/push/submission_v5.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"])
    for i,r in enumerate(ranked,1):
        w.writerow([r["design_id"], NEW_SEQ.get(r["design_id"], r["sequence"]), "single_chain"])
# master_v5
with open("runs/push/master_v5.csv","w",newline="") as fh:
    mc=list(next(iter(master.values())).keys())
    w=csv.DictWriter(fh,fieldnames=mc); w.writeheader()
    for did in master: w.writerow(master[did])

md5=lambda p: hashlib.md5(open(p,'rb').read()).hexdigest()
print("=== final20_v5 ranking ===")
print(f"{'rank':>4} {'design':24}{'score':>6} {'xreact_min':>10} {'drep':>6} {'nov':>4} {'expr':>18}")
for i,r in enumerate(ranked,1):
    m=master[r["design_id"]]; xm=min(f(m["h_ipSAE"]),f(m["m_ipSAE"]))
    print(f"{i:>4} {r['design_id'][:23]:24}{score(r['design_id']):>6.2f} {xm:>10.3f} {m['md_drep']:>6} {m['novelty']:>4} {expr_verdict(r['design_id'])[:18]:>18}")
print(f"\nsubmission_v5.csv md5={md5('runs/push/submission_v5.csv')}")
print(f"new designs ranked at: "+", ".join(f"{i}" for i,r in enumerate(ranked,1) if r['design_id'] in NEW_SEQ))
