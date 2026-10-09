#!/usr/bin/env python3
"""final20_v4: SAME 20 validated sequences as v3 (byte-identical), RE-RANKED by the integrated
composite (EXPRESSION-first, switch, cross-react, interface-strength; ipSAE DEMOTED to a weak
tiebreaker per the EGFR post-mortem). Adds orthogonal-metric QC columns. No design changes.
tnf_K90 flagged as the sole genuine weak link (ins54 genuine, not a Pro artifact) -> swap-target
pending a validated generation survivor. b2m0 instability documented as a proline false-positive."""
import csv, os, hashlib
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
v3=list(csv.DictReader(open("runs/push/final20_v3.csv")))
master={r["design_id"]:r for r in csv.DictReader(open("runs/push/master_v3.csv"))}

def f(x,d=0.0):
    try: return float(x)
    except: return d
def score(did):
    m=master[did]; ins=f(m["instab"]); fc=f(m["switch_FC"]); dr=f(m["md_drep"])
    hi=f(m["h_ipSAE"]); mi=f(m["m_ipSAE"]); dl=f(m["delta_hm"],9); ds=f(m["dSASA"]); nc=f(m["contacts"])
    sap=f(m["SAP"]); nov=m["novelty"]; npro=int(f(m["npro"]))
    s=0.0
    # EXPRESSION: SAP (real cell-free aggregation signal) dominates; instability index is a weak soft gate,
    # and a high value driven by >=3 prolines is a documented false-positive (ProtParam penalizes Pro-dipeptides).
    s-= 1.5 if sap>=2.0 else 0.0
    s+= 1.0 if sap<1.0 else 0.0
    if ins>=45: s-= (0.2 if npro>=3 else 0.5)   # Pro-FP softened vs genuine
    elif ins>=40: s-= 0.1
    s+= min(fc/20,3.0)+min(dr/2,1.5)                  # pH-switch (#1 axis)
    s+= 1.5 if dl<=0.1 else (0.5 if dl<=0.2 else 0)   # cross-react (#2)
    s+= min(ds/1500,2.0)+min(nc/50,1.5)              # interface strength (real binding signal)
    s+= 0.3*(hi+mi)                                   # ipSAE: weak tiebreaker only
    s+= 1.0 if nov=="4/4" else 0
    return s

# expression verdict per design (Pro-FP calibrated)
def expr_verdict(did):
    m=master[did]; ins=f(m["instab"]); npro=int(f(m["npro"])); sap=f(m["SAP"]); dm=int(f(m["deam"])); hyd=int(f(m["hydrun"]))
    notes=[]
    if ins>=45 and npro>=3: notes.append(f"ins{ins:.0f}=Pro-FP({npro}Pro)")  # proline false-positive
    elif ins>=45: notes.append(f"ins{ins:.0f}-GENUINE")
    elif ins>=40: notes.append(f"ins{ins:.0f}-borderline")
    if hyd>=6 and sap<2: notes.append(f"HYD{hyd}=buried(SAP{sap:.1f})")  # buried-core false positive
    if dm: notes.append(f"deam{dm}-minor")
    if sap>=2: notes.append(f"SAP{sap:.1f}-AGG")
    return "clean" if not notes else ";".join(notes)

ranked=sorted(v3,key=lambda r:-score(r["design_id"]))

# new schema: v3 cols + QC cols
qc_cols=["submit_rank","comp_score","h_ipSAE","m_ipSAE","delta_hm","dSASA","contacts","hbond","SAP","expr_verdict"]
hdr=list(v3[0].keys())+qc_cols
out="runs/push/final20_v4.csv"
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=hdr); w.writeheader()
    for i,r in enumerate(ranked,1):
        did=r["design_id"]; m=master[did]
        r=dict(r)
        r.update(submit_rank=i, comp_score=f"{score(did):.2f}", h_ipSAE=m["h_ipSAE"], m_ipSAE=m["m_ipSAE"],
                 delta_hm=m["delta_hm"], dSASA=m["dSASA"], contacts=m["contacts"], hbond=m["hbond"],
                 SAP=m["SAP"], expr_verdict=expr_verdict(did))
        w.writerow(r)
md5=lambda p: hashlib.md5(open(p,'rb').read()).hexdigest()
# verify sequences byte-identical to v3 (no design changes)
v3seq={r["design_id"]:r["sequence"] for r in v3}
v4=list(csv.DictReader(open(out)))
assert all(r["sequence"]==v3seq[r["design_id"]] for r in v4), "SEQ CHANGED!"
assert {r["design_id"] for r in v4}=={r["design_id"] for r in v3}, "SET CHANGED!"
print(f"wrote {out}  md5 {md5(out)[:8]}  (v3 {md5('runs/push/final20_v3.csv')[:8]})")
print(f"sequences byte-identical to v3: YES | set identical: YES | 20 designs re-ranked")
print(f"\n{'rank':>4} {'design':22}{'score':>6}{'expr_verdict':30} switch/xreact/iface")
for i,r in enumerate(v4,1):
    m=master[r["design_id"]]
    print(f"{i:>4} {r['design_id'][:21]:22}{r['comp_score']:>6}  {r['expr_verdict'][:28]:30} FC{m['switch_FC']}/drep{m['md_drep']} Δ{m['delta_hm']} dSASA{m['dSASA']}")
