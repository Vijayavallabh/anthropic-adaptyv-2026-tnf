#!/usr/bin/env python3
"""submission_v6: apply PORTAL-VERIFIED novelty (ground truth from Adaptyv novelty check on submission_v5)
and REMOVE b1m12_r6_A92S (portal novelty 2/4, BELOW the 3/4 hard gate -> submit-blocker).
Portal corrected several internal 4/4 estimates to 3/4 (only b1m12_A41H_vsR27 is truly 4/4).
Result = 19 portal-confirmed (>=3/4) designs, re-ranked with the SAME expression-first score() as v5.
This is the new guaranteed-valid FLOOR. Upload USER-GATED. v4/v5 kept intact as history."""
import csv, os, hashlib
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)

REMOVE={"b1m12_r6_A92S"}  # portal 2/4 -> DQ
# Adaptyv portal novelty check on submission_v5 (GROUND TRUTH, overrides internal estimates)
PORTAL_NOV={
 "b1m12_r6_A92S":"2/4","TNFpH07_sh4_r4":"3/4","TNFpH06_sh3":"3/4","b1m12_A41H_vsR27":"4/4",
 "b2m0_c70":"3/4","b1m7_A66H":"3/4","b1m9_A12H":"3/4","b4m2_A57H":"3/4","b2m0_c20":"3/4",
 "b4m2_2A87H_K101_r3":"3/4","b1m12_2A67H_K112_r1":"3/4","b1m7_2A57H_K183_r1":"3/4",
 "TNFpH06_sh5":"3/4","TNFpH06_sh2":"3/4","des_59_1_H51":"3/4","b4m2_r1":"3/4",
 "TNFpH06":"3/4","TNFpH07":"3/4","tnf_K90_tnf_K90_2_model_3_3_sh4":"3/4","des_76_2_H55":"3/4"}

f20=list(csv.DictReader(open("runs/push/final20_v5.csv")))
master={r["design_id"]:r for r in csv.DictReader(open("runs/push/master_v5.csv"))}
# apply portal novelty (ground truth) to both tables
for r in f20:
    if r["design_id"] in PORTAL_NOV: r["novelty"]=PORTAL_NOV[r["design_id"]]
for did,m in master.items():
    if did in PORTAL_NOV: m["novelty"]=PORTAL_NOV[did]
# drop the DQ
f20=[r for r in f20 if r["design_id"] not in REMOVE]
for did in REMOVE: master.pop(did,None)
assert len(f20)==19, f"expected 19, got {len(f20)}"

def f(x,d=0.0):
    try: return float(x)
    except: return d
def score(did):  # identical to v5 score(): expression-first, min-based cross-react, graded SAP
    m=master[did]; ins=f(m["instab"]); fc=f(m["switch_FC"]); dr=f(m["md_drep"])
    hi=f(m["h_ipSAE"]); mi=f(m["m_ipSAE"]); dl=f(m["delta_hm"],9); ds=f(m["dSASA"]); nc=f(m["contacts"])
    sap=f(m["SAP"]); nov=m["novelty"]; npro=int(f(m["npro"]))
    s=0.0
    s+= 1.0 if sap<1.0 else 0.0
    if sap>1.5: s-= min(0.4*(sap-1.5), 2.0)
    if ins>=45: s-= (0.2 if npro>=3 else 0.5)
    elif ins>=40: s-= 0.1
    s+= min(fc/20,3.0)+min(dr/2,1.5)
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

ranked=sorted(f20,key=lambda r:-score(r["design_id"]))
# hard-gate assertion: every submitted design >=3/4
for r in ranked:
    assert int(r["novelty"][0])>=3, f"HARD-GATE FAIL: {r['design_id']} novelty {r['novelty']}"

base_cols=[c for c in f20[0].keys() if c not in
           ("submit_rank","comp_score","h_ipSAE","m_ipSAE","delta_hm","dSASA","contacts","hbond","SAP","expr_verdict")]
qc_cols=["submit_rank","comp_score","h_ipSAE","m_ipSAE","delta_hm","dSASA","contacts","hbond","SAP","expr_verdict"]
with open("runs/push/final20_v6.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=base_cols+qc_cols); w.writeheader()
    for i,r in enumerate(ranked,1):
        did=r["design_id"]; m=master[did]; r={k:r.get(k,"-") for k in base_cols}
        r.update(submit_rank=i,comp_score=f"{score(did):.2f}",h_ipSAE=m["h_ipSAE"],m_ipSAE=m["m_ipSAE"],
                 delta_hm=m["delta_hm"],dSASA=m["dSASA"],contacts=m["contacts"],hbond=m["hbond"],
                 SAP=m["SAP"],expr_verdict=expr_verdict(did))
        w.writerow(r)
seqmap={r["design_id"]:r["sequence"] for r in f20}
with open("runs/push/submission_v6.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"])
    for r in ranked: w.writerow([r["design_id"], seqmap[r["design_id"]], "single_chain"])
with open("runs/push/master_v6.csv","w",newline="") as fh:
    mc=list(next(iter(master.values())).keys())
    w=csv.DictWriter(fh,fieldnames=mc); w.writeheader()
    for did in master: w.writerow(master[did])

md5=lambda p: hashlib.md5(open(p,'rb').read()).hexdigest()
print("=== final20_v6 ranking (19 portal-confirmed >=3/4) ===")
print(f"{'rank':>4} {'design':30}{'score':>6} {'xreact_min':>10} {'drep':>6} {'nov':>4} {'expr':>18}")
for i,r in enumerate(ranked,1):
    m=master[r["design_id"]]; xm=min(f(m["h_ipSAE"]),f(m["m_ipSAE"]))
    print(f"{i:>4} {r['design_id'][:29]:30}{score(r['design_id']):>6.2f} {xm:>10.3f} {m['md_drep']:>6} {m['novelty']:>4} {expr_verdict(r['design_id'])[:18]:>18}")
print(f"\nrows={len(ranked)}  all>=3/4: PASS")
print(f"removed: b1m12_r6_A92S (portal 2/4)")
print(f"submission_v6.csv md5={md5('runs/push/submission_v6.csv')}")
print(f"4/4 designs: "+", ".join(r['design_id'] for r in ranked if r['novelty']=='4/4'))
