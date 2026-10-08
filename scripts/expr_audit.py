#!/usr/bin/env python3
"""Adaptyv protein-qc expression/liability audit on the shipped portfolio.
Deterministic flags only (ESM2 PLL / folding dG deferred to the model track)."""
import csv, re, sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis

HYD = set("FILVWYMAC")  # aggregation-prone hydrophobics for run detection

def max_run(seq, charset):
    best = cur = 0
    for c in seq:
        if c in charset: cur += 1; best = max(best, cur)
        else: cur = 0
    return best

def audit(name, seq):
    pa = ProteinAnalysis(seq)
    instab = pa.instability_index()
    gravy = pa.gravy()
    arom = pa.aromaticity()
    ncys = seq.count("C")
    hyd_run = max_run(seq, HYD)
    kr_run = max_run(seq, set("KR"))
    deam = len(re.findall(r"N[GST]", seq))  # NG/NS/NT deamidation/isomer sequons
    # standard gate
    flags = []
    if instab >= 40: flags.append(f"INSTAB({instab:.0f})")
    if gravy >= 0.4: flags.append(f"GRAVY({gravy:.2f})")
    if ncys % 2 == 1: flags.append(f"ODDCYS({ncys})")
    if hyd_run >= 6: flags.append(f"HYDRUN({hyd_run})")
    if kr_run >= 3: flags.append(f"KR{kr_run}")
    if deam >= 1: flags.append(f"DEAM({deam})")
    verdict = "PASS" if not flags else "FLAG"
    # stringent extra
    strflags = []
    if instab >= 30: strflags.append("instab30")
    if gravy >= 0.2: strflags.append("gravy0.2")
    return dict(name=name, L=len(seq), instab=round(instab,1), gravy=round(gravy,2),
                arom=round(arom,3), ncys=ncys, hyd_run=hyd_run, kr_run=kr_run,
                deam=deam, verdict=verdict, flags=",".join(flags) or "clean",
                stringent=",".join(strflags) or "clean")

path = sys.argv[1]
rows = list(csv.DictReader(open(path)))
print(f"{'design':18} {'L':>3} {'instab':>6} {'gravy':>6} {'arom':>5} {'nC':>2} {'hyd':>3} {'kr':>2} {'dm':>2} {'STD':>5}  flags | stringent")
npass = 0
for r in rows:
    a = audit(r["design_id"], r["sequence"])
    npass += a["verdict"] == "PASS"
    print(f"{a['name']:18} {a['L']:>3} {a['instab']:>6} {a['gravy']:>6} {a['arom']:>5} {a['ncys']:>2} "
          f"{a['hyd_run']:>3} {a['kr_run']:>2} {a['deam']:>2} {a['verdict']:>5}  {a['flags']} | {a['stringent']}")
print(f"\nSTANDARD-gate PASS: {npass}/{len(rows)}")
