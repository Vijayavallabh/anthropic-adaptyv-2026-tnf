#!/usr/bin/env python3
"""BOLD BET B: Mosaic (Escalante; gradient design in relaxed sequence space, won Adaptyv Nipah/TREM2) de-novo
binder generation against TNF-alpha. Composite Boltz-2 objective; extract SEQUENCES only for INDEPENDENT
validation (our af2ig + Boltz) — never trust Mosaic's own confidence (self-consistency trap).
Run with third_party/mosaic/.venv/bin/python -I (JAX-CUDA). From repo root.
Usage: mosaic_tnf.py <n_designs> <out_csv> [seed0] [target=human|mouse|joint]"""
import sys, os, csv
import numpy as np
import jax
from mosaic.optimizers import simplex_APGM
from mosaic.common import TOKENS
import mosaic.losses.structure_prediction as sp
from mosaic.models.boltz2 import Boltz2
from mosaic.structure_prediction import TargetChain
from mosaic.proteinmpnn.mpnn import ProteinMPNN
from mosaic.losses.protein_mpnn import InverseFoldingSequenceRecovery

R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
N   = int(sys.argv[1]) if len(sys.argv)>1 else 8
OUT = sys.argv[2] if len(sys.argv)>2 else f"{R}/runs/mosaic/designs.csv"
SEED0 = int(sys.argv[3]) if len(sys.argv)>3 else 0
TGT = sys.argv[4] if len(sys.argv)>4 else "human"
os.makedirs(os.path.dirname(OUT), exist_ok=True)
HUMAN = open(f"{R}/runs/buried2/human_tnf.txt").read().strip()
MOUSE = open(f"{R}/runs/buried2/mouse_tnf_mature.txt").read().strip()
BINDER_LEN = 80
print(f"[mosaic] target={TGT} n={N} binder_len={BINDER_LEN}; loading Boltz2 + ProteinMPNN ...", flush=True)
model = Boltz2()
mpnn = ProteinMPNN.from_pretrained()

def mk_loss(target_seq):
    feats, wtr = model.binder_features(binder_length=BINDER_LEN, chains=[TargetChain(target_seq)])
    L = model.build_loss(
        loss=2*sp.BinderTargetContact() + sp.WithinBinderContact()
             + 5.0*InverseFoldingSequenceRecovery(mpnn, temp=jax.numpy.array(0.01)),
        features=feats)
    return L, feats, wtr

lossH, featsH, wtrH = mk_loss(HUMAN)
if TGT in ("joint","mouse"):
    lossM, featsM, wtrM = mk_loss(MOUSE)
# joint = bind BOTH human AND mouse (bispecific cross-react in the objective)
loss = lossH if TGT=="human" else (lossM if TGT=="mouse" else (lambda *a,**k: tuple(x+y for x,y in zip(lossH(*a,**k), lossM(*a,**k)))) )
# simple additive joint: sum the two scalar losses (aux dicts merged loosely)
if TGT=="joint":
    def loss(seq, key):
        vH,aH = lossH(seq,key); vM,aM = lossM(seq,key)
        return vH+vM, {**{f"h_{k}":v for k,v in aH.items()}, **{f"m_{k}":v for k,v in aM.items()}}

def design(seed):
    _, PSSM = simplex_APGM(loss_function=loss, n_steps=75,
        x=jax.nn.softmax(0.5*jax.random.gumbel(key=jax.random.key(seed), shape=(BINDER_LEN,20))),
        stepsize=0.1, momentum=0.9)
    PSSM_sharper, _ = simplex_APGM(loss_function=loss, n_steps=50, x=PSSM, stepsize=0.5, scale=1.5, momentum=0.0)
    seq = "".join(TOKENS[i] for i in np.array(PSSM_sharper.argmax(-1)))
    iptm=None
    try:
        pred = model.predict(PSSM=PSSM_sharper, features=featsH, writer=wtrH, key=jax.random.key(0))
        iptm=float(pred.iptm)
    except Exception as e:
        print(f"  [predict-triage skip] {str(e)[:100]}", flush=True)
    return seq, iptm

rows=[]
for i in range(N):
    try:
        seq, iptm = design(SEED0+i)
        rows.append((f"mos{TGT[0]}_{SEED0+i}", seq, iptm if iptm is not None else ""))
        print(f"[mosaic] design {i}: mosaic_iptm(triage)={iptm} len={len(seq)} {seq}", flush=True)
    except Exception as e:
        print(f"[mosaic] design {i} FAIL: {str(e)[:250]}", flush=True)
with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence","mosaic_iptm_triage"]); [w.writerow(r) for r in rows]
print(f"[mosaic] wrote {len(rows)}/{N} designs -> {OUT}", flush=True)
print("MOSAIC_GEN_DONE", flush=True)
