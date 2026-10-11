#!/usr/bin/env python3
"""Portfolio expression-RISK distribution from the per-design success probabilities in
runs/expr_v11/consensus.json. The critique's point is that each predictor is a probability, not a
proof; so the honest portfolio statement is not "all 20 pass" but the Poisson-binomial distribution
of how many of the 20 actually express (independent Bernoulli(p_i), p_i from expr_consensus.py).

Prints E[yield], P(>=k) for useful k, and the flagship (c19) express probability. Read-only."""
import json, sys
from pathlib import Path
R = Path("/home/prachh/v/anthropic-adaptyv-2026-tnf")
d = json.load(open(R/"runs/expr_v11/consensus.json"))
ps = [x['p_express'] for x in d['designs']]
names = {x['name']: x['p_express'] for x in d['designs']}
# Poisson-binomial PMF via DP
pmf = [1.0]
for p in ps:
    nxt = [0.0]*(len(pmf)+1)
    for k, v in enumerate(pmf):
        nxt[k]   += v*(1-p)
        nxt[k+1] += v*p
    pmf = nxt
n = len(ps); E = sum(ps)
def Pge(k): return sum(pmf[k:])
print(f"portfolio n={n}  E[expressing]={E:.2f}  (literature-anchored p_i; solubility_loaded={d.get('solubility_model_loaded')})")
print("cumulative P(>= k express):")
for k in [8,10,12,14,16,18]:
    if k<=n: print(f"  P(>={k:2d}) = {Pge(k):.3f}")
mode = max(range(n+1), key=lambda k: pmf[k])
print(f"most-likely count (mode) = {mode}   P(exactly {mode}) = {pmf[mode]:.3f}")
flag = names.get('des_3050_6_c19')
print(f"flagship des_3050_6_c19 P(express) = {flag}")
# switch hedge: prob at least one strong multi-His switch design expresses
switch_rich = ['des_3050_6_c19','des_4098_1_b14','des_4098_1_c44','des_4098_1_c44_w6_3','des_818_4_c13S','des_3050_6_H52-54-79','des_180_6_H9-14-38']
sr = [names[n2] for n2 in switch_rich if n2 in names]
p_none = 1.0
for p in sr: p_none *= (1-p)
print(f"P(>=1 of {len(sr)} fold-robust-His switch designs expresses) = {1-p_none:.4f}   (switch-survival hedge)")
