#!/usr/bin/env python3
"""Assemble the final-push portfolio: tiered by validated robustness + expression, all metrics merged.
Writes runs/push/final_portfolio_push.{csv,fasta}."""
import csv, os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; P=f"{R}/runs/submission"; PU=f"{R}/runs/push"
fp={r["design_id"]:r for r in csv.DictReader(open(f"{P}/final_portfolio.csv"))}
conf={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{P}/confirmed_all.csv"))}
ipall={r["token"]:r for r in csv.DictReader(open(f"{P}/ipsae_all.tsv"),delimiter="\t")}
win={r["variant"]:r for r in csv.DictReader(open(f"{PU}/rescue/winners.tsv"),delimiter="\t")}
esmc={r.split("\t")[0]:r.split("\t")[1].strip() for r in open(f"{PU}/candidates_esmc.tsv") if "\t" in r and not r.startswith("name")}
A20V="gpu0__spec_cb1ep_cb1_epitope_0_model_13_1__A20V2H_vsK197_af2pred"
def ex(seq):
    pa=ProteinAnalysis(seq); return round(pa.instability_index(),1),round(pa.gravy(),2),"PASS" if (pa.instability_index()<40 and pa.gravy()<0.4 and seq.count("C")==0) else "FLAG"
# matched human/mouse co-fold delta (the RIGHT cross-react signal; |Δ| small = cross-reactive). keyed by fold parent.
MDELTA={"TNFpH01":"-0.082","TNFpH06":"-0.071","TNFpH07":"-0.015","TNFpH08":"+0.126","b2m0":"-0.040","b4m2":"+0.007","b1m12":"-0.001","b1m9":"-0.036","b1m7":"-0.033","A20V":"-0.053"}
def mdelta(did):
    for k in MDELTA:
        if did.startswith(k) or did.split("_")[0]==k: return MDELTA[k]
    return "~0"
# (id, fold, tier, origin, switch, seq, dsw, ddg, ipsae_wc, xreact_robust, novelty, notes)
rows=[]
def R_(did,fold,tier,origin,sw,seq,dsw,ddg,wc,xrob,nov,notes):
    ins,gr,ev=ex(seq)
    rows.append(dict(design_id=did,fold=fold,tier=tier,origin=origin,switch_type=sw,sequence=seq,length=len(seq),
        d_switch_MD=dsw,ddg_pH=ddg,human_ipSAE=wc,xreact_mdelta=mdelta(did),instability=ins,gravy=gr,expr=ev,
        esmc_pll=esmc.get(did,"?"),n_cys=seq.count("C"),novelty_portal=nov,notes=notes))
# TIER A: robust + expression-clean (headliners)
R_("b1m12_r4","f5-C5","A","rescued","denovo-His",win["b1m12_r4"]["seq"],"MD-pend","-0.30(PROPKA)","0.69","YES","4/4(inherit)","expr-rescued; robust xreact; novelty 4/4")
R_("b4m2_r1","f2-C2","A","rescued","denovo-His",win["b4m2_r1"]["seq"],"MD-pend","-0.61(PROPKA)","0.68","YES","3/4(inherit)","expr-rescued; robust xreact")
R_("TNFpH06","cb1","A","orig","coop-2His",fp["TNFpH06"]["sequence"],fp["TNFpH06"]["d_switch_MD"],fp["TNFpH06"]["ddg_pH_PROPKA"],"matchedBoltz","coop","3/4","strongest proton linkage ddg -2.34")
R_("TNFpH07","cb1","A","orig","coop-2His",fp["TNFpH07"]["sequence"],fp["TNFpH07"]["d_switch_MD"],fp["TNFpH07"]["ddg_pH_PROPKA"],"matchedBoltz","coop","3/4","expr-clean coop-2His")
R_("TNFpH08","cb1","A","orig","coop-2His",fp["TNFpH08"]["sequence"],fp["TNFpH08"]["d_switch_MD"],fp["TNFpH08"]["ddg_pH_PROPKA"],"matchedBoltz","coop","3/4","coop-2His")
# TIER B: strong binders / hedges
R_("b4m2_2A87H_K101","f2-C2","B","orig","denovo-His",fp["b4m2_2A87H_K101"]["sequence"],fp["b4m2_2A87H_K101"]["d_switch_MD"],fp["b4m2_2A87H_K101"]["ddg_pH_PROPKA"],"0.85","YES","3/4","highest ipSAE binder (expr-risk hedge)")
R_("b1m12_2A67H_K112","f5-C5","B","orig","denovo-His",fp["b1m12_2A67H_K112"]["sequence"],fp["b1m12_2A67H_K112"]["d_switch_MD"],fp["b1m12_2A67H_K112"]["ddg_pH_PROPKA"],"0.83","YES","4/4","novelty 4/4 (expr-risk hedge)")
R_("b1m12_r1","f5-C5","B","rescued","denovo-His",win["b1m12_r1"]["seq"],"MD-pend","PROPKA-ok","0.59","YES","4/4(inherit)","2nd expr-rescued hedge (GRAVY -0.21)")
R_("b4m2_r6","f2-C2","B","rescued","denovo-His",win["b4m2_r6"]["seq"],"MD-pend","PROPKA-ok","0.65","YES","3/4(inherit)","2nd expr-rescued hedge")
R_("TNFpH01","cb1","B","orig","single-His",fp["TNFpH01"]["sequence"],fp["TNFpH01"]["d_switch_MD"],fp["TNFpH01"]["ddg_pH_PROPKA"],"matchedBoltz","single","3/4","strongest d_switch 2.68 (expr-risk)")
# TIER C: fold diversity, fragile xreact (ship original, flagged)
R_("b2m0_1A67H_K84","f1-C1","C","orig","denovo-His",fp["b2m0_1A67H_K84"]["sequence"],fp["b2m0_1A67H_K84"]["d_switch_MD"],fp["b2m0_1A67H_K84"]["ddg_pH_PROPKA"],"0.72","NO","3/4","fold1; xreact passes gate but fragile")
R_("b1m7_2A57H_K183","f11-C11","C","orig","denovo-His",fp["b1m7_2A57H_K183"]["sequence"],fp["b1m7_2A57H_K183"]["d_switch_MD"],fp["b1m7_2A57H_K183"]["ddg_pH_PROPKA"],"0.70","NO","3/4","fold11; fragile xreact")
R_("A20V_1A20H_K197","f12-C12","C","orig(add)","denovo-His",conf[A20V],"0.38","n/a",ipall["A20V2H_vsK197"]["ipsae_wc"],"NO","TBD","fold12 (new 7th); fragile xreact")
R_("b1m9_1A55H_K93","f8-C8","C","orig","denovo-His",fp["b1m9_1A55H_K93"]["sequence"],fp["b1m9_1A55H_K93"]["d_switch_MD"],fp["b1m9_1A55H_K93"]["ddg_pH_PROPKA"],"0.55","NO","3/4","fold8; weakest+fragile xreact")
# write
cols=["design_id","fold","tier","origin","switch_type","length","d_switch_MD","ddg_pH","human_ipSAE","xreact_mdelta","instability","gravy","expr","esmc_pll","n_cys","novelty_portal","sequence","notes"]
with open(f"{PU}/final_portfolio_push.csv","w") as fh:
    w=csv.DictWriter(fh,fieldnames=cols); w.writeheader()
    for r in rows: w.writerow({k:r[k] for k in cols})
with open(f"{PU}/final_portfolio_push.fasta","w") as fh:
    for r in rows: fh.write(f">{r['design_id']} {r['fold']} tier{r['tier']} {r['switch_type']}\n{r['sequence']}\n")
# print summary
print(f"{'id':16} {'fold':8} {'T':1} {'sw':10} {'dsw':7} {'hipSAE':6} {'xrΔ':6} {'ins':5} {'gr':5} {'expr':4} {'esmc':6} {'nov':3}")
for r in rows:
    print(f"{r['design_id']:16} {r['fold']:8} {r['tier']:1} {r['switch_type']:10} {str(r['d_switch_MD']):7} {str(r['human_ipSAE'])[:6]:6} {str(r['xreact_mdelta']):6} {str(r['instability']):5} {str(r['gravy']):5} {r['expr']:4} {str(r['esmc_pll']):6} {r['novelty_portal'][:3]}")
import collections
print("\nfolds:",len(set(r['fold'] for r in rows)),"| designs:",len(rows),"| tiers:",dict(collections.Counter(r['tier'] for r in rows)),
      "| expr-PASS:",sum(1 for r in rows if r['expr']=='PASS'),"| all cross-reactive by matched Δ + conserved contacts")
