#!/bin/bash
# pd3 amplification downstream: once the pd3 partial-diffusion wave finishes, run
#   MPNN (redesign binder chain A) -> full-atom-target AF2-ig -> collect passers -> canon + His433
#   obj1 switch -> engineered re-AF2-ig (confirm the installed switch still binds; the step that was
#   previously unwired). Produces results/rfd_{af2_passers,obj1_switch}_pd3.csv, runs/rfd_passers_pd3/,
#   and runs/rfd_af2ig_eng_pd3/ -- everything select_final_v2.py needs for source "pd3_amp".
# MD + cross-react + Boltz on the competitive pd3 picks, and the selector source-loop entry, are done
# after this (targeted, only on designs that rank), to avoid validating the whole batch.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
echo "waiting for PD3 wave... $(date -u +%T)"
for i in $(seq 1 240); do
  grep -q PD3_WAVE_DONE logs/pd3_wave.log 2>/dev/null && { echo "pd3 wave done after ~$((i*15))s"; break; }
  sleep 15
done
N=$(find runs/rfd -path '*pd3_*' -name 'des_*.pdb' 2>/dev/null | grep -v _traj | wc -l)
echo "=== pd3 downstream: MPNN + full-atom af2ig on $N backbones (binder=chain A) $(date -u +%T) ==="
bash scripts/run_betA_downstream.sh "runs/rfd/pd3_*" "0-7" 2 pd3 A 2>&1 | tail -16
echo "=== pd3 switch: collect passers (pae<11) + canon + obj1 His433 scan $(date -u +%T) ==="
bash scripts/run_betA_switch.sh pd3 11 100 2>&1 | tail -16
NPASS=$(ls runs/rfd_passers_pd3/*.pdb 2>/dev/null | wc -l)
echo "=== pd3 engineered re-af2ig (confirm installed switch still binds) $(date -u +%T) ==="
# stage engineered best_pdbs (chain A = installed seq) for the pd3 obj1 winners
envs/md/bin/python - <<'PY'
import csv, shutil
from pathlib import Path
REPO=Path('.'); out=REPO/'runs'/'eng_foldin_pd3'
if out.exists(): shutil.rmtree(out)
out.mkdir(parents=True)
objdir=REPO/'runs'/'rfd_obj1_pd3'; n=0
obj=REPO/'results'/'rfd_obj1_switch_pd3.csv'
if obj.exists():
    for r in csv.DictReader(open(obj)):
        if r.get('error') or not r.get('binder_seq'): continue
        try:
            if float(r['ddg_pH'])<=0: continue
        except: continue
        best=objdir/f"{r['design']}__best.pdb"
        if not best.exists() and r.get('best_pdb'): best=REPO/r['best_pdb']
        if best.exists(): shutil.copy(best, out/f"{r['design']}.pdb"); n+=1
print("staged",n,"pd3 engineered complexes ->",out)
PY
if [ -d runs/eng_foldin_pd3 ] && [ "$(ls runs/eng_foldin_pd3/*.pdb 2>/dev/null | wc -l)" != "0" ]; then
  envs/md/bin/python scripts/run_af2ig.py --pdbdir runs/eng_foldin_pd3 --out runs/rfd_af2ig_eng_pd3 --gpus 0-2 --mode fast 2>&1 | tail -4
fi
echo "PD3_DOWNSTREAM_DONE $(date -u +%TZ)"
echo "pd3 passers: $NPASS | obj1 rows: $(($(wc -l < results/rfd_obj1_switch_pd3.csv 2>/dev/null || echo 1)-1)) | eng af2ig: $(cat runs/rfd_af2ig_eng_pd3/gpu*/out.sc 2>/dev/null | grep -c SCORE: || echo 0)"
