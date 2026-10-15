#!/usr/bin/env python
"""Analyze metadynamics-TI ddG results vs the plain-TI baseline.

ddG_switch = <dG_bound> - <dG_free>  (>0 => releases at low pH). Reports, per design:
  - meta per-seed dG for each leg + seed spread (max-min)  vs  baseline plain-TI spread
  - pooled ddG_switch with propagated error, and the z = ddG / SE
  - sign vs the known control truth (des_820 +, b14 -, c44 +)
Validation PASS if: signs correct AND bound/free seed spreads shrink materially vs baseline.
"""
import json, glob, sys, math
import numpy as np

ROOT = "/home/prachh/v/anthropic-adaptyv-2026-tnf"
META = f"{ROOT}/runs/fe/meta"
BASE = f"{ROOT}/runs/fe/out"
TRUTH = {"des_820_3_r3_sa4_r2_8": "+", "des_4098_1_b14": "-", "des_4098_1_c44": "+"}

def legvals(d, pat):
    vs, errs = [], []
    for f in sorted(glob.glob(pat)):
        j = json.load(open(f)); vs.append(j["dG_kcal"]); errs.append(j.get("err_kcal", float("nan")))
    return np.array(vs), np.array(errs)

def summarize(name):
    mb, mbe = legvals(name, f"{META}/{name}_bound_meta_s*.json")
    mf, mfe = legvals(name, f"{META}/{name}_free_meta_s*.json")
    bb, _ = legvals(name, f"{BASE}/{name}_bound_s*.json")
    bf, _ = legvals(name, f"{BASE}/{name}_free_s*.json")
    print(f"\n=== {name}  (truth {TRUTH.get(name,'?')}) ===")
    if len(bb): print(f"  baseline bound n={len(bb)} mean={bb.mean():+7.2f} spread={np.ptp(bb):5.2f}   free n={len(bf)} mean={bf.mean():+7.2f} spread={np.ptp(bf):5.2f}  ->  baseline ddG={bb.mean()-bf.mean():+.2f}")
    if not len(mb) or not len(mf):
        print(f"  META incomplete: bound {len(mb)}  free {len(mf)}"); return None
    sb, sf = np.ptp(mb) if len(mb) > 1 else 0.0, np.ptp(mf) if len(mf) > 1 else 0.0
    ddg = mb.mean() - mf.mean()
    # error: SE of each leg mean across seeds (between-seed), combined
    seb = mb.std(ddof=1) / math.sqrt(len(mb)) if len(mb) > 1 else float(np.nanmean(mbe))
    sef = mf.std(ddof=1) / math.sqrt(len(mf)) if len(mf) > 1 else float(np.nanmean(mfe))
    se = math.sqrt(seb**2 + sef**2)
    z = ddg / se if se > 0 else float("nan")
    print(f"  META     bound n={len(mb)} mean={mb.mean():+7.2f} spread={sb:5.2f} [{', '.join(f'{x:+.1f}' for x in mb)}]")
    print(f"  META     free  n={len(mf)} mean={mf.mean():+7.2f} spread={sf:5.2f} [{', '.join(f'{x:+.1f}' for x in mf)}]")
    print(f"  META ddG_switch = {ddg:+.2f} +/- {se:.2f} kcal/mol  (z={z:+.1f})  sign={'+' if ddg>0 else '-'}  {'OK' if ('+' if ddg>0 else '-')==TRUTH.get(name) else 'WRONG'}")
    base_spread = max(np.ptp(bb) if len(bb)>1 else 0, np.ptp(bf) if len(bf)>1 else 0)
    meta_spread = max(sb, sf)
    print(f"  SPREAD: baseline max-leg {base_spread:.1f}  ->  meta max-leg {meta_spread:.1f}  ({'COLLAPSED' if meta_spread < base_spread*0.6 else 'not collapsed'})")
    return dict(name=name, ddg=ddg, se=se, z=z, sign_ok=('+' if ddg>0 else '-')==TRUTH.get(name),
                meta_spread=meta_spread, base_spread=base_spread)

def main():
    names = sys.argv[1:] or list(TRUTH)
    res = [summarize(n) for n in names]
    res = [r for r in res if r]
    if res:
        ok = sum(r["sign_ok"] for r in res); coll = sum(r["meta_spread"] < r["base_spread"]*0.6 for r in res)
        print(f"\n==== VALIDATION: signs {ok}/{len(res)} correct;  spread collapsed {coll}/{len(res)} ====")
        print("PASS" if ok == len(res) and coll >= max(1, len(res)-1) else "REVIEW")

if __name__ == "__main__":
    main()
