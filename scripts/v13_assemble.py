#!/usr/bin/env python3
"""Assemble v13 from v12 by swapping aggregate-only / MD-refuted members for MD-confirmed, makeable,
novel, distinct-lineage rebuild winners. Applies the strict decision rule (core iff strict>=2 OR
MD-causal; aggregate linkage alone insufficient). Conservative by default: swaps ONLY the forced-out
members (des_962_2_mh MD-refuted + any keep-MD NONRESPOND among TNFpH06/sh3), and prints bonus options.

Reads live MD verdicts, the gauntlet screen, and the novelty proxy. Prints a full decision report,
writes runs/push/submission_v13.csv, and lists bonus swap candidates for review.
Usage: envs/uplifting-esmc/bin/python scripts/v13_assemble.py [--bonus N]
"""
import sys, glob, json, csv, os, collections
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
BONUS=0
if "--bonus" in sys.argv: BONUS=int(sys.argv[sys.argv.index("--bonus")+1])

def md_verdict(jp):
    d=json.load(open(jp)); agg=d.get("agg",{})
    def g(v,m): return agg.get(v,{}).get(m)
    hip,hie,gln=g("HIP","his_lys_mindist"),g("HIE","his_lys_mindist"),g("GLN","his_lys_mindist")
    sem=max([x for x in (g("HIP","sem_his_lys"),g("HIE","sem_his_lys")) if x is not None] or [0.0])
    dsw=(hip-hie) if None not in (hip,hie) else None
    diso=(hip-gln) if None not in (hip,gln) else None
    v="NA"
    if dsw is not None:
        if dsw>2*sem and dsw>0 and (diso is None or diso>0): v="CONFIRMED"
        elif dsw<=2*sem: v="NONRESPOND"
        else: v="MARGINAL"
    return dict(name=d["name"],d_switch=dsw,d_iso=diso,sem=sem,verdict=v)

def load_md(pat):
    out={}
    for f in glob.glob(pat):
        try: r=md_verdict(f); out[r["name"]]=r
        except Exception as e: print(f"  (md parse fail {f}: {e})")
    return out

keepmd=load_md(f"{R}/runs/rebuild/keepmd/*.json")
gauntmd=load_md(f"{R}/runs/rebuild/gaunt/md/*.json")
gauntmd.update(load_md(f"{R}/runs/rebuild/gaunt/md_q3/*.json"))  # queue3 if run
screen={r["name"]:r for r in json.load(open(f"{R}/runs/rebuild/gaunt/screen.json"))}
nov=json.load(open(f"{R}/runs/rebuild/gaunt/novelty_proxy.json"))
novcand=nov["rebuild_candidates"]; novceil=nov["cal_pass_range"][1]
v12=[l.split(",")[0] for l in open(f"{R}/runs/push/submission_v12.csv").read().splitlines()[1:] if l.strip()]
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/rebuild/gauntlet.csv"))}

def lineage(n): return "_".join(n.split("_")[:2])  # des_1060, des_495, des_1559

# ---- swap-OUT baselines (same-pipeline axes: min cross-react, expr tier, novelty proxy, switch) ----
#   per lesson-portfolio-swap-verify-all-axes: a swap-in must not trade a higher-ranked axis down.
SWAPOUT={
 "des_962_2_mh_9_32_35":{"xr":0.846,"expr":"A/sol0.83","nov":0.955,"switch":"MD-NONRESPOND(dSw0.40)/strict0"},
 "TNFpH06":{"xr":0.648,"expr":"A/sol0.84","nov":0.944,"switch":"strict0"},
 "TNFpH06_sh3":{"xr":0.620,"expr":"A/sol0.83","nov":0.944,"switch":"strict0"},
}
# ---- swap-outs (forced) ----
forced=["des_962_2_mh_9_32_35"]   # strict-0 AND MD-NONRESPOND (wave10)
for n in ("TNFpH06","TNFpH06_sh3"):
    r=keepmd.get(n)
    if r and r["verdict"]=="NONRESPOND": forced.append(n)
print("=== keep-MD verdicts (TNFpH06/sh3/des_818_4) ===")
for n in ("TNFpH06","TNFpH06_sh3","des_818_4_c13S"):
    r=keepmd.get(n); print(f"  {n:16s} {r['verdict'] if r else 'PENDING':12s} dSwitch={r['d_switch'] if r else '--'}")

# ---- ranked rebuild winners (CONFIRMED + all gauntlet gates + novelty in-regime) ----
print("\n=== rebuild MD verdicts (gauntlet) ===")
winners=[]
for n,r in sorted(gauntmd.items(), key=lambda kv:-(kv[1]['d_switch'] or -9)):
    sc=screen.get(n,{}); seq=seqs.get(n,"")
    totH=seq.count("H")
    inreg = novcand.get(n,1.0)<=novceil
    gates = (sc.get("hum_link",0)<0 and (sc.get("chai_link") or 0)<0
             and min(sc.get("hum_iptm",0) or 0, sc.get("mouse_iptm",0) or 0)>=0.70
             and (sc.get("esm_plddt",0) or 0)>=0.70 and (sc.get("boltz_mono_plddt",0) or 0)>=0.70
             and (sc.get("cross_tm",0) or 0)>=0.60 and (sc.get("solub",0) or 0)>=0.50)
    ok = (r["verdict"]=="CONFIRMED" and gates and inreg)
    print(f"  {n:28s} {r['verdict']:11s} dSw={r['d_switch']:.2f} H={totH} xr={min(sc.get('hum_iptm',0) or 0,sc.get('mouse_iptm',0) or 0):.2f} nov={novcand.get(n,0):.3f} -> {'WINNER' if ok else 'drop'}")
    if ok:
        # composite: switch strength + cross-react + expression - His penalty
        comp=(r["d_switch"] + min(sc["hum_iptm"],sc["mouse_iptm"]) + sc["esm_plddt"] - 0.05*totH)
        winners.append(dict(name=n,lineage=lineage(n),totH=totH,d_switch=r["d_switch"],
                            xr=min(sc["hum_iptm"],sc["mouse_iptm"]),comp=comp))
winners.sort(key=lambda w:-w["comp"])

# ---- select winners for the forced swaps, lineage-diverse (prefer new lineage, then low His) ----
CAP={"des_1060":2,"des_495":2,"des_1559":1}
used=collections.Counter()
nslots=len(forced)+BONUS
picks=[]
# pass 1: one per distinct lineage (diversity first)
for w in winners:
    if len(picks)>=nslots: break
    if used[w["lineage"]]==0 and used[w["lineage"]]<CAP.get(w["lineage"],1):
        picks.append(w); used[w["lineage"]]+=1
# pass 2: fill remaining up to caps, prefer low His
for w in sorted(winners,key=lambda w:(w["totH"],-w["comp"])):
    if len(picks)>=nslots: break
    if w in picks: continue
    if used[w["lineage"]]<CAP.get(w["lineage"],1):
        picks.append(w); used[w["lineage"]]+=1

# ---- PAIR swap-outs <-> swap-ins greedily by cross-react (highest-xr out gets highest-xr in) ----
#      so each swap-in's #2 (cross-react) >= its paired swap-out's where possible (no #2 traded down).
outs_sorted=sorted(forced, key=lambda n:-SWAPOUT.get(n,{}).get("xr",0))
ins_sorted=sorted(picks, key=lambda w:-w["xr"])
pairs=list(zip(outs_sorted, ins_sorted))
print(f"\n=== SWAP PLAN (conservative; bonus={BONUS}) -- per-axis check (lesson: no higher axis traded down) ===")
print(f"{'OUT':28s}{'xr':>6}{'switch(out)':>22s}   -> {'IN':28s}{'xr':>6}  {'switch(in)':>16s}  axis-verdict")
for o,w in pairs:
    ob=SWAPOUT.get(o,{})
    xr_ok = w["xr"]>=ob.get("xr",0)-0.005
    sw_ok = True  # swap-in is MD-CONFIRMED by construction; swap-out is NONRESPOND/strict0
    verdict = "DOMINATES" if xr_ok else f"#1 up / #2 down ({ob.get('xr'):.2f}->{w['xr']:.2f}); OK lexicographic(#1>#2)"
    print(f"{o:28s}{ob.get('xr',0):>6.2f}{ob.get('switch','')[:22]:>22s}   -> {w['name']:28s}{w['xr']:>6.2f}  {'CONFIRMED dSw%.2f'%w['d_switch']:>16s}  {verdict}")
if len(picks)<len(forced):
    print(f"  !! only {len(picks)} winners for {len(forced)} forced-out; v13 would drop below 20 -- launch queue3 backups or keep some")

# ---- write v13 ----
keep=[n for n in v12 if n not in forced]
newnames=[p["name"] for p in picks]
v13=keep+newnames
if len(v13)!=20:
    print(f"\n!! v13 has {len(v13)} designs (need 20) -- NOT writing, resolve first")
else:
    with open(f"{R}/runs/push/submission_v13.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"])
        allseq={**{r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/push/submission_v12.csv"))}, **seqs}
        for n in v13: w.writerow([n, allseq[n].split(":")[0], "single_chain"])
    print(f"\nwrote runs/push/submission_v13.csv ({len(v13)} designs)")

print(f"\n=== bonus swap candidates (weakest strict-0 v12 binders, if adding more winners) ===")
strict=json.load(open(f"{R}/runs/push/strict_rerank_v12.json"))
s0=[s["design"] for s in strict if (s["robust"]==0) and s["design"] in keep and s["design"] not in ("des_2106_2_mh_37_43_50_75",)]
print(f"  strict-0 binders still in: {s0}")
print(f"  unused confirmed winners: {[w['name'] for w in winners if w['name'] not in newnames]}")
