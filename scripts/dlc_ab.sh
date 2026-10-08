#!/bin/bash
# DLC robustness A/B: re-fold the final-20 complexes with af2ig OFF vs ON (DLC, fast mode) and report
# whether the train/eval LayerNorm correction flips any pae<10 gate call or reshuffles the ranking.
# Per the analysis, this is a ROBUSTNESS CHECK (+ pLDDT secondary signal), NOT a re-rank.
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
MANI=${1:-results/submission_manifest_v2.csv}
IN=runs/dlc_ab/in; rm -rf runs/dlc_ab; mkdir -p "$IN"

envs/md/bin/python - "$MANI" "$IN" <<'PY'
import csv, sys, shutil
from pathlib import Path
REPO = Path('.')
MANI, IN = sys.argv[1], sys.argv[2]
DIRS = ['rfd_passers_pd','rfd_passers_pd2','rfd_passers_w1fa','rfd_passers_w2fa']
EXTRA = ['obj1_foldin2/{d}__best.pdb','passers/{d}.pdb','filtered/{d}.pdb','cand_final2_canon/{d}.pdb']
n = 0
for r in csv.DictReader(open(MANI)):
    d = r['design']
    p = next((REPO/'runs'/x/f'{d}.pdb' for x in DIRS if (REPO/'runs'/x/f'{d}.pdb').exists()), None)
    if p is None:
        for t in EXTRA:
            q = REPO/'runs'/t.format(d=d)
            if q.exists(): p = q; break
    if p: shutil.copy(p, f"{IN}/{r['name']}.pdb"); n += 1
    else: print("MISSING complex for", d)
print("gathered", n, "complexes ->", IN)
PY

echo "=== fold OFF (GPUs 0-3) + ON dlc_q=0.9 (GPUs 4-7), fast mode $(date -u +%T) ==="
envs/md/bin/python scripts/run_af2ig.py --pdbdir "$IN" --out runs/dlc_ab/off --gpus 0-3 --mode fast > logs/dlc_ab_off.log 2>&1 &
OFF=$!
AF2IG_DLC_Q=0.9 envs/md/bin/python scripts/run_af2ig.py --pdbdir "$IN" --out runs/dlc_ab/on --gpus 4-7 --mode fast > logs/dlc_ab_on.log 2>&1 &
ON=$!
wait $OFF $ON
echo "=== DLC A/B comparison $(date -u +%T) ==="
envs/md/bin/python - <<'PY'
import glob
def parse(run):
    rows = {}
    for p in glob.glob(f"runs/dlc_ab/{run}/out.sc") + glob.glob(f"runs/dlc_ab/{run}/gpu*/out.sc"):
        hdr = None
        for ln in open(p):
            t = ln.split()
            if ln.startswith('SCORE:') and 'description' in ln: hdr = t[1:]
            elif ln.startswith('SCORE:') and hdr:
                dd = dict(zip(hdr, t[1:])); rows[dd['description'].replace('_af2pred','')] = dd
    return rows
off, on = parse('off'), parse('on')
common = sorted(set(off) & set(on))
flips = 0; churn = []
print(f"{'name':14s} {'pae_off':>7s} {'pae_on':>7s} {'dpae':>6s} {'plddtB_off':>10s} {'plddtB_on':>10s} gate")
for d in common:
    po = float(off[d].get('pae_interaction', 99)); pn = float(on[d].get('pae_interaction', 99))
    lo = float(off[d].get('plddt_binder', 0)); lnn = float(on[d].get('plddt_binder', 0))
    g = 'OK' if (po < 10) == (pn < 10) else 'FLIP'; flips += g == 'FLIP'
    churn.append(pn - po)
    print(f"{d[:14]:14s} {po:7.2f} {pn:7.2f} {pn-po:+6.2f} {lo:10.2f} {lnn:10.2f} {g}")
import statistics
print(f"\nn={len(common)} | pae<10 gate FLIPS: {flips} | mean dpae={statistics.mean(churn):+.3f} "
      f"median={statistics.median(churn):+.3f} | mean dplddt_binder="
      f"{statistics.mean(float(on[d].get('plddt_binder',0))-float(off[d].get('plddt_binder',0)) for d in common):+.2f}")
print("interpretation: 0 flips + small dpae => selection is ROBUST to DLC; keep OFF gate. "
      "Flips => scrutinize those borderline picks.")
PY
echo "DLC_AB_DONE $(date -u +%TZ)"