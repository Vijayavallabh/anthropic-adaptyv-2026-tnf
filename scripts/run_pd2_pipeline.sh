#!/bin/bash
# Auto-process the diverse-PD wave once it finishes: soluble-MPNN -> full-atom-target af2ig ->
# collect passers -> canon -> obj1 His433 switch. Produces results/rfd_{af2_passers,obj1_switch}_pd2.csv
# and runs/rfd_passers_pd2/ that select_final_v2 reads (source "pd_diverse"). GPU use kicks in only
# after the wave frees the GPUs; obj1 is CPU (PyRosetta), so GPUs are free again for the next wave.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
echo "waiting for PD-diverse wave... $(date -u +%T)"
for i in $(seq 1 600); do
  grep -q PD_DIVERSE_WAVE2_DONE logs/pd_diverse_wave.log 2>/dev/null && { echo "wave done after ~$((i*30))s"; break; }
  sleep 30
done
N=$(find runs/rfd -path '*pd2_*' -name 'des_*.pdb' 2>/dev/null | wc -l)
echo "=== pd2 downstream: MPNN + full-atom af2ig on $N backbones (binder=chain A) $(date -u +%T) ==="
bash scripts/run_betA_downstream.sh "runs/rfd/pd2_*" "0-7" 2 pd2 A 2>&1 | tail -14
echo "=== pd2 switch: collect passers (pae<11) + obj1 His433 scan $(date -u +%T) ==="
bash scripts/run_betA_switch.sh pd2 11 32 2>&1 | tail -14
echo "PD2_PIPELINE_DONE $(date -u +%TZ)"
echo "pd2 af2 passers: $(($(wc -l < results/rfd_af2_passers_pd2.csv 2>/dev/null || echo 1)-1)) | obj1 rows: $(($(wc -l < results/rfd_obj1_switch_pd2.csv 2>/dev/null || echo 1)-1))"