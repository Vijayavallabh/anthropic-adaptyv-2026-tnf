#!/usr/bin/env python3
"""Composite rank of MD-CONFIRMED TNF pH-switch binders, weighted by competition objective priority
(pH-selectivity > cross-reactivity > affinity) with expression as a near-gate (EGFR post-mortem:
SolubleMPNN/expression decisive). Merges signal_table + footprint + ESM2 naturalness + ProtParam.
Z-scores each axis within the confirmed set; prints ranked composite + hard flags. No new deps.
"""
import csv, os, math
from Bio.SeqUtils.ProtParam import ProteinAnalysis
ROOT="/home/prachh/v/anthropic-adaptyv-2026-tnf"
def read_tsv(p):
    with open(p) as fh:
        r=csv.DictReader(fh, delimiter="\t"); return list(r)
sig={r["variant"]:r for r in read_tsv(f"{ROOT}/runs/md/signal_table.tsv")}
foot={r["variant"]:r for r in read_tsv(f"{ROOT}/runs/md/footprint.tsv")}
esm={r["id"]:float(r["esm2_ll"]) for r in read_tsv(f"{ROOT}/runs/md/esm_naturalness.tsv")}
seqs={r[0]:r[1] for r in list(csv.reader(open(f"{ROOT}/runs/md/confirmed.csv")))[1:]}
def num(x,d=0.0):
    try: return float(x)
    except: return d
rows=[]
for v,s in sig.items():
    if s.get("md_verdict")!="CONFIRMED": continue
    seq=seqs.get(v);
    if not seq: continue
    pa=ProteinAnalysis(seq); gravy=pa.gravy(); instab=pa.instability_index()
    rows.append(dict(
        variant=v, lineage=s["lineage"], his=s["his"], anchor=s["anchor"], seq=seq,
        d_switch=num(s["d_switch"]), ddg=num(s["ddg_static"]), n_div=num(foot.get(v,{}).get("n_divergent",99)),
        pae=num(s["pae_af2ig"]), plddt=num(s["plddt"]), rmsd=num(s["rmsd"]), tm=num(s["tm_novelty"]),
        esm=esm.get(v,-0.5), gravy=gravy, instab=instab,
        boltz_delta=num(s.get("boltz_delta",0)),
    ))
def z(vals):
    m=sum(vals)/len(vals); sd=(sum((x-m)**2 for x in vals)/len(vals))**0.5 or 1.0
    return [(x-m)/sd for x in vals]
# axes (sign so that higher z = better)
Z={}
Z["d_switch"]=z([r["d_switch"] for r in rows])          # MD release ladder, higher better
Z["ddg"]=z([-r["ddg"] for r in rows])                   # more negative ddg = release, so -ddg higher better
Z["xreact"]=z([-r["n_div"] for r in rows])              # fewer divergent = better cross-react
Z["pae"]=z([-r["pae"] for r in rows])                   # lower pae better
Z["plddt"]=z([r["plddt"] for r in rows])
Z["rmsd"]=z([-r["rmsd"] for r in rows])
Z["instab"]=z([-r["instab"] for r in rows])             # lower instability better
Z["esm"]=z([r["esm"] for r in rows])                    # higher naturalness better
Z["novelty"]=z([-r["tm"] for r in rows])                # lower TM = more novel
for i,r in enumerate(rows):
    pH   = Z["d_switch"][i] + Z["ddg"][i]                # PRIMARY
    xr   = Z["xreact"][i]                                # SECOND
    aff  = (Z["pae"][i] + Z["plddt"][i] + Z["rmsd"][i])/3.0   # THIRD
    expr = (Z["instab"][i] + Z["esm"][i])/2.0            # near-gate
    nov  = Z["novelty"][i]
    r["pH"],r["xr"],r["aff"],r["expr"],r["nov"]=pH,xr,aff,expr,nov
    r["composite"]=3.0*pH + 2.0*xr + 1.0*aff + 1.5*expr + 1.0*nov
    # hard flags (Adaptyv protein-qc expression + cross-react risk)
    fl=[]
    if r["instab"]>=40: fl.append("UNSTABLE")
    if r["gravy"]>=0.4: fl.append("GRAVY")
    if r["n_div"]>=5: fl.append("XREACT")
    if r["tm"]>=0.80: fl.append("NOTNOVEL")
    r["flags"]=",".join(fl) if fl else "clean"
rows.sort(key=lambda r:-r["composite"])
print(f"{'#':>2} {'variant':46} {'lin':4} {'comp':>6} | {'pH':>5} {'xr':>5} {'aff':>5} {'exp':>5} {'nov':>5} | "
      f"{'dsw':>5} {'ddg':>6} {'ndiv':>4} {'pae':>5} {'instab':>6} {'esmLL':>6} {'tm':>5}  flags")
for i,r in enumerate(rows):
    print(f"{i+1:>2} {r['variant'][:46]:46} {r['lineage']:4} {r['composite']:6.2f} | "
          f"{r['pH']:5.2f} {r['xr']:5.2f} {r['aff']:5.2f} {r['expr']:5.2f} {r['nov']:5.2f} | "
          f"{r['d_switch']:5.2f} {r['ddg']:6.3f} {int(r['n_div']):>4} {r['pae']:5.2f} {r['instab']:6.1f} {r['esm']:6.3f} {r['tm']:5.2f}  {r['flags']}")
# lineage summary
print("\nlineage counts (confirmed):", end=" ")
from collections import Counter
print(dict(Counter(r["lineage"] for r in rows)))
with open(f"{ROOT}/runs/md/composite_rank.tsv","w") as fh:
    cols=["variant","lineage","his","anchor","composite","pH","xr","aff","expr","nov",
          "d_switch","ddg","n_div","pae","plddt","rmsd","instab","gravy","esm","tm","boltz_delta","flags","seq"]
    fh.write("\t".join(cols)+"\n")
    for r in rows: fh.write("\t".join(str(r[c]) for c in cols)+"\n")
print("wrote runs/md/composite_rank.tsv")
