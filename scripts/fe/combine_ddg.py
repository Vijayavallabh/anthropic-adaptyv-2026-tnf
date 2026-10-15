#!/usr/bin/env python
"""Combine per-leg alchemical dG into ddG_switch per design + rank with the 2-sigma rule.

ddG_switch = dG_bound - dG_free   (kcal/mol; >0 => protonation disfavors binding => releases at low pH)
Larger positive ddG = sharper pH-switch. Pairs seeds by index; SEM combines seed spread + MBAR error.
Ranking: design i strictly > j only if (ddG_i - ddG_j) > 2*sqrt(SEM_i^2 + SEM_j^2); else a tied band.
"""
import glob, json, os, math, collections
import numpy as np
R = "/home/prachh/v/anthropic-adaptyv-2026-tnf"
OUT = f"{R}/runs/fe/out"

def legs(design, leg):
    d = {}
    for f in glob.glob(f"{OUT}/{design}_{leg}_s*.json"):
        s = int(f.split("_s")[-1].split(".")[0]); d[s] = json.load(open(f))
    return d

designs = sorted(set(os.path.basename(f).split("_bound_s")[0]
                     for f in glob.glob(f"{OUT}/*_bound_s*.json")))
rows = []
for dn in designs:
    b, fr = legs(dn, "bound"), legs(dn, "free")
    seeds = sorted(set(b) & set(fr))
    if not seeds: continue
    ddg = [b[s]["dG_kcal"] - fr[s]["dG_kcal"] for s in seeds]
    mbar = [math.hypot(b[s]["err_kcal"], fr[s]["err_kcal"]) for s in seeds]
    mean = float(np.mean(ddg))
    sem_seed = float(np.std(ddg, ddof=1)/math.sqrt(len(ddg))) if len(ddg) > 1 else float("nan")
    sem_mbar = float(math.sqrt(np.mean(np.array(mbar)**2))/math.sqrt(len(seeds)))
    sem = sem_mbar if len(ddg) == 1 else math.hypot(sem_seed, sem_mbar)
    rows.append(dict(design=dn, ddG=mean, sem=sem, nseed=len(seeds),
                     per_seed=[round(x,2) for x in ddg],
                     his=b[seeds[0]].get("switch_his"),
                     trans=round(float(np.mean([b[s].get("mean_self_transition",float('nan')) for s in seeds])),2)))
rows.sort(key=lambda r: -r["ddG"])

print(f"{'design':34s} {'ddG':>7s} {'+/-':>6s} {'nseed':>5s} {'per_seed':>18s} {'His':>16s} mix")
for r in rows:
    print(f"{r['design']:34s} {r['ddG']:+7.2f} {r['sem']:6.2f} {r['nseed']:5d} "
          f"{str(r['per_seed']):>18s} {str(r['his']):>16s} {r['trans']}")

print("\n=== pairwise 2-sigma resolution (row strictly > col) ===")
for i, a in enumerate(rows):
    for b2 in rows[i+1:]:
        gap = a["ddG"] - b2["ddG"]; thr = 2*math.hypot(a["sem"], b2["sem"])
        verdict = "RESOLVED" if gap > thr else "tied"
        print(f"  {a['design'][:24]:24s} > {b2['design'][:24]:24s}  dGap={gap:+.2f}  2sig={thr:.2f}  [{verdict}]")
json.dump(rows, open(f"{OUT}/ddG_ranking.json","w"), indent=2)
print(f"\nwrote {OUT}/ddG_ranking.json  ({len(rows)} designs)")
