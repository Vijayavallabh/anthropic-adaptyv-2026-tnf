#!/bin/bash
# After pd2 MD frees the GPUs, fold the de novo wave-2 backbones (novel-backbone novelty axis):
# soluble-MPNN (binder=chain B for de novo) -> full-atom-target af2ig -> collect passers -> obj1 switch.
# Produces results/rfd_{af2_passers,obj1_switch}_w2fa.csv + runs/rfd_passers_w2fa/ (select_final_v2
# reads tag w2fa, source rfd_denovo). Standard (DLC-OFF) af2ig — selection uses the calibrated gate.
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
echo "waiting for pd2 MD to finish... $(date -u +%T)"
for i in $(seq 1 240); do
  grep -q PD2_MD_DONE logs/pd2_md.log 2>/dev/null && { echo "pd2 MD done after ~$((i*30))s"; break; }
  sleep 30
done
N=$(ls runs/rfd/w2_*/des_*.pdb 2>/dev/null | grep -vE '_traj' | wc -l)
echo "=== de novo w2 downstream: MPNN + full-atom af2ig on $N backbones (binder=chain B) $(date -u +%T) ==="
bash scripts/run_betA_downstream.sh "runs/rfd/w2_*" "0-7" 2 w2fa B 2>&1 | tail -14
echo "=== w2fa switch: collect passers (pae<11) + obj1 His433 scan $(date -u +%T) ==="
bash scripts/run_betA_switch.sh w2fa 11 32 2>&1 | tail -14
echo "W2FA_PIPELINE_DONE $(date -u +%TZ)"
echo "w2fa af2 passers: $(($(wc -l < results/rfd_af2_passers_w2fa.csv 2>/dev/null || echo 1)-1)) | obj1 rows: $(($(wc -l < results/rfd_obj1_switch_w2fa.csv 2>/dev/null || echo 1)-1))"