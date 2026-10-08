#!/bin/bash
# 7D12-seeded downstream: once the 7D12 partial-diffusion wave finishes, run
#   MPNN (redesign binder chain A) -> full-atom-target AF2-ig -> collect passers -> canon + His433
#   obj1 switch -> engineered re-AF2-ig. Produces results/rfd_{af2_passers,obj1_switch}_7d12.csv,
#   runs/rfd_passers_7d12/, runs/rfd_af2ig_eng_7d12/ for select_final_v2.py source "7d12_amp".
# MD + cross-react + Boltz on the competitive 7D12 picks are done after (targeted).
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
echo "waiting for 7D12 wave... $(date -u +%T)"
for i in $(seq 1 240); do
  grep -q 7D12_WAVE_DONE logs/7d12_wave.log 2>/dev/null && { echo "7D12 wave done after ~$((i*15))s"; break; }
  sleep 15
done
N=$(find runs/rfd -path '*7d12_*' -name 'des_*.pdb' 2>/dev/null | grep -v _traj | wc -l)
echo "=== 7D12 downstream: MPNN + full-atom af2ig on $N backbones (binder=chain A) $(date -u +%T) ==="
bash scripts/run_betA_downstream.sh "runs/rfd/7d12_*" "0-7" 2 7d12 A 2>&1 | tail -16
echo "=== 7D12 switch: collect passers (pae<11) + canon + obj1 His433 scan $(date -u +%T) ==="
bash scripts/run_betA_switch.sh 7d12 11 100 2>&1 | tail -16
echo "=== 7D12 engineered re-af2ig (confirm installed switch still binds) $(date -u +%T) ==="
envs/md/bin/python - <<'PY'
import csv, shutil
from pathlib import Path
REPO=Path('.'); out=REPO/'runs'/'eng_foldin_7d12'
if out.exists(): shutil.rmtree(out)
out.mkdir(parents=True)
objdir=REPO/'runs'/'rfd_obj1_7d12'; n=0
obj=REPO/'results'/'rfd_obj1_switch_7d12.csv'
if obj.exists():
    for r in csv.DictReader(open(obj)):
        if r.get('error') or not r.get('binder_seq'): continue
        try:
            if float(r['ddg_pH'])<=0: continue
        except: continue
        best=objdir/f"{r['design']}__best.pdb"
        if not best.exists() and r.get('best_pdb'): best=REPO/r['best_pdb']
        if best.exists(): shutil.copy(best, out/f"{r['design']}.pdb"); n+=1
print("staged",n,"7D12 engineered complexes ->",out)
PY
if [ -d runs/eng_foldin_7d12 ] && [ "$(ls runs/eng_foldin_7d12/*.pdb 2>/dev/null | wc -l)" != "0" ]; then
  envs/md/bin/python scripts/run_af2ig.py --pdbdir runs/eng_foldin_7d12 --out runs/rfd_af2ig_eng_7d12 --gpus 0-2 --mode fast 2>&1 | tail -4
fi
echo "7D12_DOWNSTREAM_DONE $(date -u +%TZ)"
echo "7d12 passers: $(ls runs/rfd_passers_7d12/*.pdb 2>/dev/null | wc -l) | obj1 rows: $(($(wc -l < results/rfd_obj1_switch_7d12.csv 2>/dev/null || echo 1)-1)) | eng af2ig: $(cat runs/rfd_af2ig_eng_7d12/gpu*/out.sc 2>/dev/null | grep -c SCORE: || echo 0)"
