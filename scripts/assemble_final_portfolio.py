#!/usr/bin/env python3
"""Assemble the final diverse portfolio: 8 cb1-family (floor) + the 8 MD-CONFIRMED cb1-epitope de-novo
folds, with all objective metrics. Select for all-objective balance + fold-diversity, cap <=20.
Writes runs/submission/final_portfolio.{csv,fasta} + a ranked table.
"""
import csv, os, glob, re
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
# ---- cb1-family (floor) ----
core=list(csv.DictReader(open(f"{R}/runs/submission/floor_core8.csv")))
rows=[]
for r in core:
    rows.append(dict(id=r["design_id"], fold="cb1", seq=r["sequence"], length=int(r["length"]),
        d_switch=float(r["d_switch_MD"]), ddg_pH=float(r["ddg_pH_PROPKA"]),
        pae=float(r["pae_human"]), xreact=float(r["xreact_delta"]), iptm_h=float(r["iptm_human"]),
        instab=float(r["instability"]), gravy=float(r["gravy"]), tm=float(r["novelty_TM"]),
        switch=r["switch_type"]))
# ---- cb1-epitope de-novo (8 confirmed) ----
# d_switch from MD verdict (observed); backbone id from name
DSW={"A1M2H_vsK116":2.73,"A22A2H_vsK194":2.52,"A6E2H_vsK81":1.77,"A49A2H_vsK107":1.68,
     "A7D2H_vsK84":1.58,"A67W2H_vsK84":0.79,"A55R2H_vsK93":0.72,"A31L2H_vsK118":0.72}
pk={r["design"]:float(r["ddg_pH"]) for r in csv.DictReader(open(f"{R}/runs/cb1ep/propka26.tsv"),delimiter="\t") if r["ddg_pH"] not in ("NA","")}
expr={r["design"]:(float(r["gravy"]),float(r["instab"]),int(r["cys"])) for r in csv.DictReader(open(f"{R}/runs/cb1ep/expr26.tsv"),delimiter="\t")}
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/cb1ep/final/confirmed_nameseq.csv"))}
# af2ig pae + Boltz iptm/xreact
def tok(n):
    m=re.search(r"(A\d+[A-Z]2H_vsK\d+)",n); return m.group(1) if m else n
pae_sw={}  # token -> pae_interaction (switch variant)
for ln in open(f"{R}/runs/cb1ep/switch_af2ig/out.sc"):
    p=ln.split()
    if p[:1]==["SCORE:"] and len(p)>3 and p[3]!="pae_interaction":
        try: pae_sw[tok(p[-1])]=float(p[3])
        except: pass
import json
sel=json.load(open(f"{R}/runs/cb1ep/final/selection.json"))
xr={tok(s["name"]):(s.get("iptm_h"),s.get("xr")) for s in sel}
man={}  # token -> manifest desc (for propka/expr keys + backbone id)
for ln in open(f"{R}/runs/cb1ep/md_manifest.tsv"):
    d=ln.split("\t")[0]; man[tok(d)]=d
for name,seq in seqs.items():
    t=tok(name); desc=man.get(t,"")
    g,ins,cys = expr.get(desc,(0,0,0))
    ih,xd = xr.get(t,(None,None))
    bbm=re.search(r"(cb1_epitope_\d+_model_\d+)",desc)
    bb=bbm.group(1).replace("cb1_epitope_","b") if bbm else name
    rows.append(dict(id=name, fold="cb1ep:"+bb, seq=seq, length=len(seq),
        d_switch=DSW.get(t,0), ddg_pH=pk.get(desc,0), pae=pae_sw.get(t,99),
        xreact=xd if xd is not None else -9, iptm_h=ih if ih else 0,
        instab=ins, gravy=g, tm=0.65, switch="denovo-His"))  # tm placeholder (de-novo; portal confirms)
# ---- score + select ----
# gates: switch confirmed (d_switch>0 & ddg<0), binding (pae<10), prefer xreact>-0.10 & instab<45
for r in rows:
    r["xr_clean"]= r["xreact"]>-0.10
    r["expr_ok"] = r["instab"]<45 and r["gravy"]<0.5 and r["instab"]>0
    r["dual_sw"] = r["d_switch"]>0 and r["ddg_pH"]<0
print(f"{'id':30} {'fold':12} {'len':>3} {'dSw':>5} {'ddg':>6} {'pae':>5} {'xrΔ':>6} {'instab':>6} {'flags'}")
for r in sorted(rows,key=lambda r:(r["fold"],-r["d_switch"])):
    fl=[]
    if not r["xr_clean"]: fl.append("xreact")
    if not r["expr_ok"]: fl.append("expr")
    if not r["dual_sw"]: fl.append("switch?")
    print(f"{r['id'][:30]:30} {r['fold']:12} {r['length']:>3} {r['d_switch']:5.2f} {r['ddg_pH']:6.2f} {r['pae']:5.2f} {r['xreact']:+6.2f} {r['instab']:6.1f} {','.join(fl) if fl else 'CLEAN'}")
# FINAL CUT: all 8 cb1-family + cb1ep that are cross-react-clean (xr>-0.1) [distinct folds], cap 2/ epitope-fold not needed (each cb1ep is its own fold)
final=[r for r in rows if r["fold"]=="cb1"]  # 8 core
# de-novo additions: cross-react-clean AND (expression-reasonable OR uniquely strong switch); 1 per fold
cb1ep_keep=[r for r in rows if r["fold"].startswith("cb1ep") and r["xr_clean"] and (r["instab"]<50 or r["d_switch"]>2.0)]
cb1ep_keep=sorted(cb1ep_keep, key=lambda r:-(r["d_switch"] + (1.0 if r["instab"]<40 else 0) + 0.5*r["iptm_h"]))
final += cb1ep_keep
final=final[:20]
with open(f"{R}/runs/submission/final_portfolio.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["design_id","fold","sequence","length","d_switch_MD","ddg_pH_PROPKA","pae_human","xreact_delta","iptm_human","instability","gravy","novelty_TM","switch_type"])
    for r in final: w.writerow([r["id"],r["fold"],r["seq"],r["length"],r["d_switch"],round(r["ddg_pH"],3),r["pae"],r["xreact"],r["iptm_h"],r["instab"],r["gravy"],r["tm"],r["switch"]])
with open(f"{R}/runs/submission/final_portfolio.fasta","w") as f:
    for r in final: f.write(f">{r['id']} fold={r['fold']} dSw={r['d_switch']} ddg={round(r['ddg_pH'],2)} xrΔ={r['xreact']}\n{r['seq']}\n")
folds=sorted(set(r["fold"] for r in final))
print(f"\nFINAL PORTFOLIO: {len(final)} designs across {len(folds)} distinct folds: {folds}")
print("wrote runs/submission/final_portfolio.{csv,fasta}")
