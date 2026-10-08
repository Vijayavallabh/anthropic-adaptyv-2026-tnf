#!/usr/bin/env python3
"""ESM-C 600M length-normalized pseudo-log-likelihood (judge's #1 binder/non-binder signal).
Usage: esmc_pll.py <csv name,sequence> [gpu]"""
import os, sys, csv
os.environ.setdefault("HF_HOME", "/home/prachh/v/anthropic-adaptyv-2026-tnf/weights/esmc/hf")
os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

dev = "cuda:" + (sys.argv[2] if len(sys.argv) > 2 else "0")
model = ESMC.from_pretrained("esmc_600m").to(dev).eval()

def pll(seq):
    prot = ESMProtein(sequence=seq)
    t = model.encode(prot)
    ids = t.sequence.to(dev)                       # [L+2] incl BOS/EOS
    with torch.no_grad():
        out = model.logits(t, LogitsConfig(sequence=True))
    lg = out.logits.sequence[0].float()            # [L+2, V]
    logp = torch.log_softmax(lg, dim=-1)
    # score real residues (positions 1..L), true token = ids[i]
    idx = torch.arange(1, len(ids) - 1, device=dev)
    tok = ids[idx]
    vals = logp[idx, tok]
    return float(vals.mean())                      # length-normalized PLL

rows = list(csv.DictReader(open(sys.argv[1])))
print("name\tesmc_pll")
for r in rows:
    nm = r.get("name") or r.get("design_id") or r.get("design")
    seq = r["sequence"].split(":")[0]
    try:
        print(f"{nm}\t{pll(seq):.4f}")
    except Exception as e:
        print(f"{nm}\tERR {e}", file=sys.stderr)
