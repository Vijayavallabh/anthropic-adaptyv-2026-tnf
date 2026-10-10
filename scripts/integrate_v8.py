#!/usr/bin/env python3
"""Assemble submission_v8 = NEW multi-His triple-confirmed switch core (wave-3, Baker SOTA recipe:
>=2-3 separated His repelling distinct conserved cations; stronger on the #1 pH-selectivity axis than
v7's single-His/thermo designs) + the diverse PORTAL-VERIFIED v7 floor (drop excess cb1 + SAP-flagged +
redundant twins). Cap ~2/lineage. v7 stays the fallback floor; NEW de-novo lineages are USER-portal-novelty
pending. Writes runs/push/submission_v8.csv + master_v8.csv + new_for_portal.csv. Run .venv/bin/python -I."""
import csv, os, json
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)

def load(p):
    return {r.get("name") or r.get("design_id"): r["sequence"] for r in csv.DictReader(open(p))}
wa=load("runs/wave3/variants.csv"); wb=load("runs/wave3b/variants.csv")
v7rows=list(csv.DictReader(open("runs/push/submission_v7.csv")))
v7={r["name"]:r["sequence"] for r in v7rows}

def seq(n):
    for src in (wa,wb,v7):
        if n in src: return src[n].split(":")[0].strip()
    raise KeyError(n)

# --- NEW multi-His switch core (metrics from triple/switch_ranked/plddt/sap) ---
# (variant, lineage, min, dHM, nHis, conf, mono_pLDDT, SAP, novTM, track)
NEW=[
 ("des_4098_1_b14","des_4098_1",0.587,0.106,4,3,0.90,0.1,0.655,"B"),
 ("des_4098_1_b1" ,"des_4098_1",0.787,0.070,4,2,0.81,0.3,0.758,"B"),
 ("des_180_6_H9-14-38","des_180_6",0.734,0.058,3,2,0.81,0.7,0.661,"A"),
 ("des_818_4_H16-74","des_818_4",0.713,0.105,2,2,0.88,1.2,0.635,"A"),
 ("des_3050_6_b9","des_3050_6",0.624,0.073,4,2,0.88,0.3,0.699,"B"),
 ("des_3050_6_H52-54-79","des_3050_6",0.631,0.095,3,2,0.79,0.0,0.691,"A"),
 ("des_875_5_H3-57-80","des_875_5",0.600,0.016,3,1,0.00,0.0,0.588,"A"),  # 1-conf, adds lineage, best balance
 ("des_4098_1_b3","des_4098_1",0.688,0.083,4,2,0.84,0.2,0.748,"B"),      # portal backup for des_4098_1
]
# v8 picks: cap 2/lineage on the new core -> 7 picks (des_4098_1 b14+b1, des_3050_6 b9+H52-54-79,
# des_180_6, des_818_4, des_875_5); des_4098_1_b3 kept only as a portal backup, not in the 20.
NEW_PICKS=["des_4098_1_b14","des_4098_1_b1","des_180_6_H9-14-38","des_818_4_H16-74",
           "des_3050_6_b9","des_3050_6_H52-54-79","des_875_5_H3-57-80"]
NEW_PORTAL=[n for (n,*_) in NEW]   # all 8 for the user's portal novelty check

# --- v7 floor: keep 13 diverse portal-verified; drop 7 (excess cb1 / SAP / twin / decoupled) ---
V7_KEEP=["des_287_5_H37","des_2106_2_H43","des_59_1_H51",      # portal-verified de-novo switches
         "TNFpH06","TNFpH06_sh3",                               # cb1 best 2 (was 5 -> drop 3)
         "b1m12_A41H_vsR27","b1m12_2A67H_K112_r1",
         "b2m0_c70",                                            # drop c20 (97% twin)
         "b1m7_2A57H_K183_r1",                                  # drop b1m7_A66H (switch decoupled)
         "b4m2_r1","b4m2_A57H",                                 # drop b4m2_2A87H (3rd)
         "tnf_K90_tnf_K90_2_model_3_3_sh4","b1m9_A12H"]
V7_DROP=[r["name"] for r in v7rows if r["name"] not in V7_KEEP]  # 7, available for novelty backfill

assert len(NEW_PICKS)+len(V7_KEEP)==20, (len(NEW_PICKS),len(V7_KEEP))

# --- write submission_v8.csv ---
meta={n:(lin,mn,dh,nh,cf,pl,sap,tm,tr) for (n,lin,mn,dh,nh,cf,pl,sap,tm,tr) in NEW}
rows=[]
for n in NEW_PICKS: rows.append((n, seq(n), "single_chain"))
for n in V7_KEEP:   rows.append((n, seq(n), "single_chain"))
with open("runs/push/submission_v8.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"]); [w.writerow(r) for r in rows]

# master table w/ provenance
with open("runs/push/master_v8.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","source","lineage","switch","min_ipSAE","dHM","nHis","conf_His","mono_pLDDT","SAP","novTM_proxy","portal_novelty"])
    for n in NEW_PICKS:
        lin,mn,dh,nh,cf,pl,sap,tm,tr=meta[n]
        w.writerow([n,f"wave3-{tr}",lin,"multiHis",mn,dh,nh,cf,pl,sap,tm,"PENDING-USER"])
    for n in V7_KEEP:
        w.writerow([n,"v7",n.split("_H")[0] if "_H" in n else n.rsplit("_",1)[0],"v7","","","","","","","portal>=3/4(v6)"])

# designs the user must portal-novelty-check (all 8 new, flag the riskiest)
with open("runs/push/new_for_portal.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","novTM_proxy","risk"])
    for (n,lin,mn,dh,nh,cf,pl,sap,tm,tr) in sorted(NEW,key=lambda x:x[8]):
        w.writerow([n, seq(n), tm, "HIGH(near0.8)" if tm>=0.73 else "ok(~3/4)"])

# --- FINAL-TEST ---
import hashlib
names=[r[0] for r in rows]; seqs=[r[1] for r in rows]
assert len(names)==20 and len(set(names))==20, "not 20 unique"
assert all(set(s)<=set("ACDEFGHIKLMNPQRSTVWY") for s in seqs), "bad residue"
assert all(60<=len(s)<=130 for s in seqs), [(n,len(s)) for n,s in zip(names,seqs)]
lin=[meta[n][0] if n in meta else (n.rsplit("_",1)[0]) for n in names]
from collections import Counter
print("FINAL-TEST: 20 unique, residues OK, lengths OK")
print("v8 md5:", hashlib.md5(open("runs/push/submission_v8.csv","rb").read()).hexdigest())
print(f"\nNEW switch core (7, portal-pending): "+", ".join(NEW_PICKS))
print(f"v7 floor kept (13): "+", ".join(V7_KEEP))
print(f"v7 dropped (7, backfill pool): "+", ".join(V7_DROP))
print(f"\nportal-check list -> runs/push/new_for_portal.csv (8 designs incl 1 backup)")
print("distinct lineages in v8:", len(set(lin)), "/ 20 designs")
