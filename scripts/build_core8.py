#!/usr/bin/env python3
"""Assemble the locked floor: 8 MD-CONFIRMED cb1-family designs (5 single-His + 3 cooperative 2-His),
fully annotated (binding, switch, cross-react, novelty, expression). Writes a candidate table + a
floor submission CSV + fasta. Final >8 selection waits on the generation campaigns.
"""
import csv, os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def rd(p):
    with open(p) as f: return list(csv.DictReader(f, delimiter="\t"))
comp={r["variant"]:r for r in rd(f"{R}/runs/md/composite_rank.tsv")}
# coop: metrics gathered from af2ig out.sc + md_verdict + chain-A seq
COOP={
 "coop_H19_H5":dict(seq="AQEEHEKLVKQIEEAQKKHQEYSAKNRVLALYYGAEASYLIASLLAKTHGTEELLKEAETRYKFFKDLLEAIESGDEALMEEVAERAEKELGLDREKVLAEMERLRAE",pae=6.492,plddt=88.366,dsw=1.10,his="5,19",anchor="K98+K119(coop)"),
 "coop_H30_H5":dict(seq="AQEEHEKLVKQIEEAQKKSQEYSAKNRVLHLYYGAEASYLIASLLAKTHGTEELLKEAETRYKFFKDLLEAIESGDEALMEEVAERAEKELGLDREKVLAEMERLRAE",pae=6.130,plddt=88.774,dsw=0.98,his="5,30",anchor="K98+K206(coop)"),
 "coop_H20_H5":dict(seq="AQEEHEKLVKQIEEAQKKSHEYSAKNRVLALYYGAEASYLIASLLAKTHGTEELLKEAETRYKFFKDLLEAIESGDEALMEEVAERAEKELGLDREKVLAEMERLRAE",pae=6.790,plddt=87.329,dsw=0.76,his="5,20",anchor="K98+K206(coop)"),
}
def expr(seq):
    pa=ProteinAnalysis(seq); return round(pa.gravy(),2), round(pa.instability_index(),1)
# ESM2 LL from the tsv
esm={r["id"]:float(r["esm2_ll"]) for r in rd(f"{R}/runs/md/esm_naturalness.tsv")}
rows=[]
# 5 single-His cb1 (ordered by composite) -- cb1 ONLY (cb7/dn7 are NOTNOVEL/XREACT, excluded)
for v,c in sorted(comp.items(), key=lambda kv:-float(kv[1]["composite"])):
    if c["lineage"]!="cb1": continue
    seq=c["seq"]; g,ins=expr(seq)
    rows.append(dict(src=v, lineage="cb1", switch="single-His", his=c["his"], anchor=c["anchor"],
        seq=seq, length=len(seq), pae_human=round(float(c["pae"]),2), plddt=round(float(c["plddt"]),1),
        rmsd=round(float(c["rmsd"]),2), d_switch=round(float(c["d_switch"]),2),
        n_div=int(float(c["n_div"])), tm_novelty=round(float(c["tm"]),2),
        gravy=g, instability=ins, esm2=round(float(c["esm"]),3),
        n_his=seq.count("H"), n_cys=seq.count("C"),
        boltz_delta=c.get("boltz_delta","")))
# 3 coop 2-His
for v,c in COOP.items():
    seq=c["seq"]; g,ins=expr(seq)
    rows.append(dict(src=v, lineage="cb1", switch="coop-2His", his=c["his"], anchor=c["anchor"],
        seq=seq, length=len(seq), pae_human=c["pae"], plddt=c["plddt"], rmsd=None,
        d_switch=c["dsw"], n_div=None, tm_novelty=0.65, gravy=g, instability=ins,
        esm2=round(esm.get(v,esm.get("coop_"+v,-0.31)),3),
        n_his=seq.count("H"), n_cys=seq.count("C"), boltz_delta=""))
# assign IDs by overall quality order: single-His by composite then coop by d_switch
rows_sorted = rows  # already single-His(composite) then coop
for i,r in enumerate(rows_sorted): r["design_id"]=f"TNFpH{i+1:02d}"
cols=["design_id","lineage","switch","his","anchor","length","pae_human","plddt","rmsd","d_switch",
      "n_div","tm_novelty","gravy","instability","esm2","n_his","n_cys","src","seq"]
os.makedirs(f"{R}/runs/submission",exist_ok=True)
with open(f"{R}/runs/submission/core8_table.tsv","w") as f:
    f.write("\t".join(cols)+"\n")
    for r in rows_sorted: f.write("\t".join(str(r.get(c,"")) for c in cols)+"\n")
# floor CSV (Adaptyv-style)
with open(f"{R}/runs/submission/floor_core8.csv","w") as f:
    w=csv.writer(f); w.writerow(["design_id","sequence","length","switch_type","switch_his","anchor",
        "pae_human","plddt_binder","d_switch_MD","footprint_ndiv","novelty_TM","instability","gravy","n_his","n_cys"])
    for r in rows_sorted:
        w.writerow([r["design_id"],r["seq"],r["length"],r["switch"],r["his"],r["anchor"],r["pae_human"],
            r["plddt"],r["d_switch"],r["n_div"],r["tm_novelty"],r["instability"],r["gravy"],r["n_his"],r["n_cys"]])
# fasta
with open(f"{R}/runs/submission/floor_core8.fasta","w") as f:
    for r in rows_sorted:
        f.write(f">{r['design_id']} {r['switch']} his={r['his']} pae_h={r['pae_human']} dsw={r['d_switch']} instab={r['instability']}\n{r['seq']}\n")
print(f"{'id':8} {'switch':10} {'his':6} {'len':>3} {'paeH':>5} {'dSw':>5} {'ndiv':>4} {'TM':>5} {'instab':>6} {'gravy':>5} {'nH':>2} {'nC':>2}")
for r in rows_sorted:
    print(f"{r['design_id']:8} {r['switch']:10} {str(r['his']):6} {r['length']:>3} {r['pae_human']:>5} {r['d_switch']:>5} {str(r['n_div']):>4} {r['tm_novelty']:>5} {r['instability']:>6} {r['gravy']:>5} {r['n_his']:>2} {r['n_cys']:>2}")
print("\nwrote runs/submission/core8_table.tsv, floor_core8.csv, floor_core8.fasta")
