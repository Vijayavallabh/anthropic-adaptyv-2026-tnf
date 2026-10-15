#!/usr/bin/env python
"""Stage-1 FE-screen analysis = FILTER ONLY (multiplicity-prone, NOT a result).
With 27 noisy candidates, 2-seed "hits" include expected flukes; nothing is confirmed here.
Flag generously for stage-2 (fresh independent seeds + des_820 bar). Reports expected false
positives so the output can't be misread. ddG=dG_bound-dG_free (>0 = low-pH release)."""
import glob, json, math, os
import numpy as np
OUT="/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/fe/out"
designs=sorted(set(os.path.basename(f).split("_bound_s")[0] for f in glob.glob(f"{OUT}/*_bound_s*.json")))
rows=[]
for d in designs:
    b={int(f.split("_s")[-1].split(".")[0]):json.load(open(f)) for f in glob.glob(f"{OUT}/{d}_bound_s*.json")}
    fr={int(f.split("_s")[-1].split(".")[0]):json.load(open(f)) for f in glob.glob(f"{OUT}/{d}_free_s*.json")}
    seeds=sorted(set(b)&set(fr))
    if not seeds: continue
    ddg=np.array([b[s]["dG_kcal"]-fr[s]["dG_kcal"] for s in seeds])
    mean=float(ddg.mean()); sem=float(ddg.std(ddof=1)/math.sqrt(len(ddg))) if len(ddg)>1 else float('nan')
    nh=len(b[seeds[0]].get("switch_his",[]))
    # GOLD (des_820 bar, multiplicity-safe): >=6 seeds, all same positive sign, mean>+5, z>3
    allpos=all(x>0 for x in ddg)
    gold = (len(seeds)>=6 and allpos and mean>5 and not math.isnan(sem) and sem>0 and mean/sem>3)
    # stage-1 FILTER flag (generous; expect flukes): both/all seeds>+1 OR mean>+2.5
    hit = (mean>2.5) or (len(seeds)>=2 and all(x>1 for x in ddg))
    rows.append(dict(d=d,mean=mean,sem=sem,n=len(seeds),nH=nh,allpos=allpos,
                     ddg=[round(x,1) for x in ddg.tolist()],gold=gold,hit=hit))
rows.sort(key=lambda r:-r["mean"])
print(f"{'design':44s} {'ddG':>7s} {'SEM':>5s} {'n':>2s} nH  flag      per-seed")
for r in rows:
    flag="GOLD-CONFIRMED" if r["gold"] else ("stage1-hit" if r["hit"] else "")
    print(f"{r['d']:44s} {r['mean']:+7.2f} {0 if math.isnan(r['sem']) else r['sem']:5.2f} {r['n']:2d} {r['nH']:2d}  {flag:14s} {r['ddg']}")
hits=[r['d'] for r in rows if r["hit"] and not r["gold"] and r["n"]<6]
gold=[r['d'] for r in rows if r["gold"]]
# multiplicity accounting for the 2-seed filter
n2=sum(1 for r in rows if r["n"]==2)
print(f"\n*** STAGE-1 IS A FILTER, NOT A RESULT ***  {n2} candidates at 2 seeds.")
print(f"Per-design chance of a 2-positive-seed fluke under null ~ 0.25 -> ~{0.25*n2:.0f} false 'hits' expected by luck alone.")
print(f"GOLD-CONFIRMED (>=6 seeds, all +, mean>+5, z>3; Bonferroni-safe): {gold}")
print(f"STAGE-1 HITS to re-test with 4 FRESH seeds (split-sample): {hits}")
json.dump({"gold":gold,"stage1_hits":hits,"rows":rows,"n_at_2seed":n2},
          open(f"{OUT}/screen_result.json","w"), indent=1, default=float)
print(f"wrote {OUT}/screen_result.json")
