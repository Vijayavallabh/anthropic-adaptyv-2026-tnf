#!/usr/bin/env python3
"""v9 = v8 with 2 wave4b within-lineage UPGRADES (same already-verified backbones/novelty):
  des_3050_6: b9(2-conf) -> c39 (4 confirmed His / 4 distinct cations / min 0.666 / pLDDT 0.88 / SAP 0.0)
  des_4098_1: b1(2-conf,0.787) -> c44 (min 0.861 / dHM 0.024 / pLDDT 0.82 / SAP 0.1)
Both strictly stronger on switch/binding/aggregation. c39,c44 are new seqs on portal-verified-3/4 backbones
(des_3050_6, des_4098_1) -> expected 3/4, USER re-confirms. v7 floor unchanged. Writes submission_v9.csv +
master_v9.csv + new_for_portal_v9.csv. Run .venv/bin/python -I."""
import csv, os, json, hashlib
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
def load(p): return {r.get("name") or r.get("design_id"): r["sequence"] for r in csv.DictReader(open(p))}
wa=load("runs/wave3/variants.csv"); wb=load("runs/wave3b/variants.csv"); wc=load("runs/wave4b/variants.csv")
v7={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/push/submission_v7.csv"))}
def seq(n):
    for src in (wa,wb,wc,v7):
        if n in src: return src[n].split(":")[0].strip()
    raise KeyError(n)
# NEW switch core (variant, lineage, min, dHM, nHis, conf, pLDDT, SAP, portal)
NEW=[
 ("des_3050_6_c39","des_3050_6",0.666,0.063,8,4,0.88,0.0,"PENDING(lin3/4)"),   # UPGRADE (was b9)
 ("des_4098_1_b14","des_4098_1",0.587,0.106,4,3,0.90,0.1,"3/4"),
 ("des_4098_1_c44","des_4098_1",0.861,0.024,6,2,0.82,0.1,"PENDING(lin3/4)"),   # UPGRADE (was b1)
 ("des_180_6_H9-14-38","des_180_6",0.734,0.058,3,2,0.81,0.7,"3/4"),
 ("des_818_4_H16-74","des_818_4",0.713,0.105,2,2,0.88,1.2,"3/4"),
 ("des_3050_6_H52-54-79","des_3050_6",0.631,0.095,3,2,0.79,0.0,"3/4"),
 ("des_875_5_H3-57-80","des_875_5",0.600,0.016,3,1,0.85,0.0,"4/4"),
]
NEW_PICKS=[n for (n,*_) in NEW]
meta={n:(lin,mn,dh,nh,cf,pl,sap,pn) for (n,lin,mn,dh,nh,cf,pl,sap,pn) in NEW}
V7_KEEP=["des_287_5_H37","des_2106_2_H43","des_59_1_H51","TNFpH06","TNFpH06_sh3",
         "b1m12_A41H_vsR27","b1m12_2A67H_K112_r1","b2m0_c70","b1m7_2A57H_K183_r1",
         "b4m2_r1","b4m2_A57H","tnf_K90_tnf_K90_2_model_3_3_sh4","b1m9_A12H"]
assert len(NEW_PICKS)+len(V7_KEEP)==20
rows=[(n,seq(n),"single_chain") for n in NEW_PICKS]+[(n,seq(n),"single_chain") for n in V7_KEEP]
with open("runs/push/submission_v9.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"]); [w.writerow(r) for r in rows]
with open("runs/push/master_v9.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","source","lineage","min_ipSAE","dHM","nHis","conf_His","mono_pLDDT","SAP","portal_novelty"])
    for n in NEW_PICKS:
        lin,mn,dh,nh,cf,pl,sap,pn=meta[n]; src="wave4b" if "_c" in n else ("wave3b" if "_b" in n else "wave3")
        w.writerow([n,src,lin,mn,dh,nh,cf,pl,sap,pn])
    for n in V7_KEEP: w.writerow([n,"v7",n,"","","","","","","portal>=3/4"])
# portal re-check list: only the 2 NEW upgrade seqs (c39,c44); b9/b1 are the portal-verified fallbacks
with open("runs/push/new_for_portal_v9.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","note"])
    w.writerow(["des_3050_6_c39", seq("des_3050_6_c39"), "upgrade over b9; expect 3/4 (des_3050_6 backbone); fallback=des_3050_6_b9"])
    w.writerow(["des_4098_1_c44", seq("des_4098_1_c44"), "upgrade over b1; expect 3/4 (des_4098_1 backbone); fallback=des_4098_1_b1"])
# FINAL-TEST
names=[r[0] for r in rows]; seqs=[r[1] for r in rows]
assert len(names)==20==len(set(names)); assert all(set(s)<=set("ACDEFGHIKLMNPQRSTVWY") for s in seqs)
assert all(60<=len(s)<=130 for s in seqs)
print("FINAL-TEST ok. v9 md5:", hashlib.md5(open("runs/push/submission_v9.csv","rb").read()).hexdigest())
print("swaps vs v8: des_3050_6_b9->c39 (2->4 conf), des_4098_1_b1->c44 (0.787->0.861 min)")
print("NEW switch core conf-His:", {meta[n][4]:0 for n in NEW_PICKS} if False else [(n.split('_H')[0] if '_H' in n else n, meta[n][4]) for n in NEW_PICKS])
print("portal re-check needed: des_3050_6_c39, des_4098_1_c44 (new_for_portal_v9.csv); fallbacks b9/b1")
lins=set(meta[n][0] for n in NEW_PICKS)|set(V7_KEEP)
print("distinct lineages:", len({*[meta[n][0] for n in NEW_PICKS], *[v.split('_H')[0] if '_H' in v else v.rsplit('_',1)[0] for v in V7_KEEP]}))
