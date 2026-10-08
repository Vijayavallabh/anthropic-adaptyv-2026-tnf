#!/usr/bin/env python3
"""Assemble the FINAL 20-design dominating portfolio: curated for the Challenge-2 judging axes
(pH-selectivity > cross-reactivity > affinity; novelty >=3/4 gate; expression; fold diversity).
Auto-merges every validated metric from the campaign's stable outputs. Writes
runs/push/final20.{csv,fasta} + upload_final20.csv + a ranked summary."""
import csv, os, glob, json
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PU=f"{R}/runs/push"

# ---- curated 20 (id, fold, novelty, role) ----
C=[
 # b1m9 f8 -- the FULLY-VALIDATED lead: sharp switch (PROPKA 22x + pKAI 7x) AND MD His-causal
 #           (d_repel +0.91, Gln control inert) AND strong cross-react binding (0.91/0.82) AND expresses.
 ("b1m9_1A55H_K93_sh4","f8-C8","3/4","LEAD: sharp switch FC 22x, MD His-causal (control inert), strong x-react bind, expr-clean"),
 ("b1m9_1A55H_K93_sh6_r1","f8-C8","3/4","RESCUED sharp b1m9: FC 32x, iptm 0.92/0.83, now expr-PASS"),
 ("b1m9_1A55H_K93","f8-C8","3/4","validated single-His parent (iptm 0.89)"),
 # cb1 -- sharpest THERMODYNAMIC switches (PROPKA+pKAI consensus; what BLI-at-2-pH measures), both species
 ("TNFpH06_sh3","cb1","3/4","CO-LEAD: sharp FC 59x both species, MD His-causal (d_repel+1.18), x-react, expr-clean"),
 ("TNFpH06_sh5","cb1","3/4","sharp switch FC 31x/76x (pKAI 58x, sharpest by pKAI), x-react, expr-clean"),
 ("TNFpH07_sh4_r4","cb1","3/4","RESCUED physics-max switch: FC 69x, now expr-PASS + cross-react (iptm 0.66/0.72)"),
 ("TNFpH06_sh2","cb1","3/4","MD His-causal sharp cb1 (d_repel+1.25), FC 17x/11x both species, best cb1 binding iptm 0.75, expr-clean"),
 ("TNFpH06","cb1","3/4","validated coop-2His original (lineage anchor)"),
 ("TNFpH07","cb1","3/4","validated coop-2His original (2nd anchor)"),
 # b1m12 f5 -- 4/4 novel strong binders (novelty + affinity + x-react; modest switch)
 ("b1m12_r4","f5-C5","4/4","best novel binder iptm 0.94, expr-clean, novelty 4/4"),
 ("b1m12_2A67H_K112_r1","f5-C5","4/4","RESCUED 4/4 strong binder: iptm 0.92/0.92 perfect cross-react, now expr-PASS"),
 ("b1m12_r1","f5-C5","4/4","4/4 expr-clean binder"),
 # b4m2 f2 -- strong cross-reactive binders (epitope 96-97% conserved)
 ("b4m2_r1","f2-C2","3/4","expr-clean binder iptm 0.91, epitope 97% conserved"),
 ("b4m2_2A87H_K101_r3","f2-C2","3/4","RESCUED strongest binder: iptm 0.95/0.95 perfect cross-react, FC 10x, now expr-PASS"),
 ("b4m2_r6","f2-C2","3/4","2nd expr-clean b4m2 hedge (96% conserved)"),
 # b2m0 f1 -- MD-confirmed de-novo switches (His-causal, Gln control inert)
 ("b2m0_c20","f1-C1","3/4","MD His-causal coop switch (d_repel+control inert), f1 fold"),
 ("b2m0_c70","f1-C1","3/4","MD His-causal coop switch, strongest d_repel 3.2A, f1 fold"),
 ("b2m0_1A67H_K84","f1-C1","3/4","f1 single-His original (diversity)"),
 # b1m7 f11
 ("b1m7_2A57H_K183_r1","f11-C11","3/4","RESCUED f11 binder: iptm 0.90/0.92, now expr-PASS"),
 ("b1m7_2A57H_K183_sh5","f11-C11","3/4","f11 multi-His switch variant"),
]

def load_seqs():
    s={}
    for f in ["upload_push.csv","sharp/variants.csv","sharp2/variants.csv","coop/variants.csv","rescue2/candidates.csv"]:
        p=f"{PU}/{f}"
        if os.path.exists(p):
            for r in csv.DictReader(open(p)):
                n=r.get("name") or r.get("design_id"); s.setdefault(n,r["sequence"])
    return s
def load_tsv(path,key="design_id"):
    return {r[key]:r for r in csv.DictReader(open(path),delimiter="\t")} if os.path.exists(path) else {}
seqs=load_seqs()
SH_H={};
for f in ["switch_sharpness.tsv","sharp/sharpness.tsv","sharp2/sharpness.tsv","rescue2/sharpness_h.tsv"]: SH_H.update(load_tsv(f"{PU}/{f}"))
SH_M={}
for f in ["sharp/sharpness_mouse.tsv","sharp2/sharpness_mouse.tsv"]: SH_M.update(load_tsv(f"{PU}/{f}"))
PKAI=load_tsv(f"{PU}/sharp/pkai.tsv")
def iptm(name,species):
    pats=["boltz_*_h","boltz_r2h*"] if species=="h" else ["boltz_*_m","boltz_r2m*"]
    best=None
    for pat in pats:
        for g in glob.glob(f"{R}/runs/{pat}/out/by_seed/{name}/s0/confidence_{name}_model_0.json"):
            try:
                v=json.load(open(g))["iptm"]
                if best is None or v>best: best=v
            except: pass
    return best
def md(name):
    # md_release_tnf writes --out as a directory containing <name>.json; also handle flat (coop)
    cands=[f"{PU}/sharp/md/{name}.json/{name}.json",f"{PU}/sharp/md/{name}.json",f"{PU}/coop/md/{name}.json",
           f"{PU}/rescue2/md/{name}.json/{name}.json",f"{PU}/rescue2/md/{name}.json"]
    for j in cands:
        if os.path.isfile(j) and os.path.getsize(j)>2:
            try:
                d=json.load(open(j)); dr=d.get("d_repel"); du=d.get("dud_repel")
                # His-CAUSAL = repulsion real (d_repel>0.3) AND exceeds the Gln-isostere control by >0.3
                causal = (dr is not None and dr>0.3 and (du is None or dr>du+0.3))
                return causal, dr, du
            except: pass
    return None,None,None
def expr(seq):
    pa=ProteinAnalysis(seq); ins=pa.instability_index(); gr=pa.gravy()
    return round(ins,1),round(gr,2),("PASS" if ins<40 and gr<0.4 and seq.count("C")==0 else "FLAG")

rows=[]
for did,fold,nov,role in C:
    seq=seqs.get(did,"")
    if not seq: print(f"  !! missing seq {did}"); continue
    hfc=SH_H.get(did,{}).get("fold_change","-"); hpk=SH_H.get(did,{}).get("peak_dnH","-")
    mfc=SH_M.get(did,{}).get("fold_change","-")
    pk=PKAI.get(did,{}).get("pkai_fold_change","-")
    ih=iptm(did,"h"); im=iptm(did,"m")
    causal,drep,dud=md(did); ins,gr,ev=expr(seq)
    rows.append(dict(design_id=did,fold=fold,novelty=nov,sequence=seq,length=len(seq),
        switch_h_FC=hfc,switch_h_peak=hpk,switch_m_FC=mfc,pkai_FC=pk,
        md_his_causal=("YES" if causal else ("no" if causal is False else "pend")),md_d_repel=(round(drep,2) if drep is not None else "-"),md_dud=(round(dud,2) if dud is not None else "-"),
        iptm_h=(round(ih,2) if ih else "-"),iptm_m=(round(im,2) if im else "-"),
        instability=ins,gravy=gr,expr=ev,n_cys=seq.count("C"),role=role))
# write
cols=["design_id","fold","novelty","switch_h_FC","switch_h_peak","switch_m_FC","pkai_FC","md_his_causal","md_d_repel","md_dud","iptm_h","iptm_m","instability","gravy","expr","length","n_cys","sequence","role"]
with open(f"{PU}/final20.csv","w") as fh:
    w=csv.DictWriter(fh,fieldnames=cols); w.writeheader()
    for r in rows: w.writerow({k:r[k] for k in cols})
with open(f"{PU}/final20.fasta","w") as fh:
    for r in rows: fh.write(f">{r['design_id']} {r['fold']} nov{r['novelty']}\n{r['sequence']}\n")
with open(f"{PU}/upload_final20.csv","w") as fh:
    fh.write("name,sequence,molecule_class\n")
    for r in rows: fh.write(f"{r['design_id']},{r['sequence']},single_chain\n")
# summary
print(f"{'id':20} {'fold':8} {'nov':3} {'sw_h':>6} {'sw_m':>6} {'pkai':>5} {'MDcausal':>8} {'ih':>4} {'im':>4} {'ins':>5} {'expr':4}")
for r in rows:
    print(f"{r['design_id']:20} {r['fold']:8} {r['novelty']:3} {str(r['switch_h_FC'])+'x':>6} {str(r['switch_m_FC'])+'x':>6} {str(r['pkai_FC'])+'x':>5} {r['md_his_causal']:>8} {str(r['iptm_h']):>4} {str(r['iptm_m']):>4} {str(r['instability']):>5} {r['expr']:4}")
import collections
print(f"\ndesigns: {len(rows)} | folds: {len(set(r['fold'] for r in rows))} {dict(collections.Counter(r['fold'] for r in rows))}")
print(f"expr PASS: {sum(1 for r in rows if r['expr']=='PASS')}/{len(rows)} | novelty 4/4: {sum(1 for r in rows if r['novelty']=='4/4')} | sharp(FC>=19x human): {sum(1 for r in rows if str(r['switch_h_FC']).replace('-','0').isdigit() and int(str(r['switch_h_FC']).replace('-','0'))>=19)}")
