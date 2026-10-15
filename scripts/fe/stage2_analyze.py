#!/usr/bin/env python
"""Stage-2 confirmation on FRESH seeds (3,4,5,6) only (split-sample: the 2 selection seeds that
flagged a hit are EXCLUDED from the decision estimate, to avoid winner's-curse regression).
GOLD bar (des_820 standard, multiplicity-safe for ~27 screened): all fresh seeds same sign,
|mean| > 5 kcal, and |z| > 3. Classifies each as CONFIRMED RELEASER / CONFIRMED ANTI-SWITCH /
unresolved. Also prints pooled all-seed value for context (not the decision)."""
import glob, json, math, os
import numpy as np
OUT="/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/fe/out"
cands=[l.split("\t")[0] for l in open("/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/fe/stage2_pool.tsv") if l.strip()]
FRESH={3,4,5,6}
def ddg(d, seeds_filter=None):
    b={int(f.split("_s")[-1].split(".")[0]):json.load(open(f)) for f in glob.glob(f"{OUT}/{d}_bound_s*.json")}
    fr={int(f.split("_s")[-1].split(".")[0]):json.load(open(f)) for f in glob.glob(f"{OUT}/{d}_free_s*.json")}
    seeds=sorted(set(b)&set(fr))
    if seeds_filter is not None: seeds=[s for s in seeds if s in seeds_filter]
    if not seeds: return None
    arr=np.array([b[s]["dG_kcal"]-fr[s]["dG_kcal"] for s in seeds])
    return seeds, arr
print(f"{'design':44s} {'fresh ddG(s3-6)':>16s} {'z':>6s}  {'pooled(all)':>11s}  verdict   per-fresh-seed")
res=[]
for d in cands:
    fr=ddg(d, FRESH); allr=ddg(d)
    if fr is None: print(f"{d:44s}  no fresh seeds yet"); continue
    _, a=fr; mean=float(a.mean()); sem=float(a.std(ddof=1)/math.sqrt(len(a))) if len(a)>1 else float('nan')
    z=mean/sem if sem and not math.isnan(sem) and sem>0 else 0
    pooled=float(allr[1].mean())
    allsame=all(x>0 for x in a) or all(x<0 for x in a)
    gold = (len(a)>=3 and allsame and abs(mean)>5 and abs(z)>3)
    verdict = ("RELEASER!" if mean>0 else "ANTI-SWITCH!") if gold else "unresolved"
    res.append((d,mean,z,pooled,verdict,len(a)))
    print(f"{d:44s} {mean:+9.2f}+/-{sem:4.2f} {z:+6.2f}  {pooled:+11.2f}  {verdict:12s} {np.round(a,1).tolist()}")
conf_rel=[d for d,m,z,p,v,n in res if v=="RELEASER!"]
conf_anti=[d for d,m,z,p,v,n in res if v=="ANTI-SWITCH!"]
print(f"\nCONFIRMED RELEASERS (fresh-seed gold bar): {conf_rel}")
print(f"CONFIRMED ANTI-SWITCHES (fresh-seed gold bar): {conf_anti}")
json.dump({"confirmed_releasers":conf_rel,"confirmed_anti":conf_anti,
           "rows":[{"d":d,"fresh_mean":m,"z":z,"pooled":p,"verdict":v,"n_fresh":n} for d,m,z,p,v,n in res]},
          open(f"{OUT}/stage2_result.json","w"), indent=1, default=float)
print(f"wrote {OUT}/stage2_result.json")
