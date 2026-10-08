#!/bin/bash
# WAVE 17 win-campaign downstream (user lifted GPU cap to all 8 H100s, 2026-10-03: "win at all costs").
# Waits for the partial-diffusion Wave A (which diversifies the 20 VALIDATED complexes) to free the
# GPUs, then runs the established, known-good mechanical chain on its backbones:
#   soluble-MPNN (design binder chain A) -> full-atom-target af2ig fold-filter
#   -> collect pae<MAXPAE passers -> canonicalize His433->B/100 -> obj1 His433 pH-switch scan.
# Produces results/rfd_af2_passers_w17.csv + results/rfd_obj1_switch_w17.csv + runs/rfd_passers_w17*/.
# Mechanical SCORING only -- the dominate-or-discard SELECTION vs the live 20 is done by hand after,
# on the priority order (pH-switch > cross-react > affinity), never swapping on noise.
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
WAVELOG=${1:-scratchpad/win/waveA_pd.log}
GLOB=${2:-runs/rfd/pd2_egfr_phsw_*}
TAG=${3:-w17}
MAXPAE=${4:-11}
echo "[$TAG] waiting for Wave A (PD_DIVERSE_WAVE2_DONE in $WAVELOG) ... $(date -u +%FT%TZ)"
for i in $(seq 1 960); do   # up to ~8h
  grep -q PD_DIVERSE_WAVE2_DONE "$WAVELOG" 2>/dev/null && { echo "[$TAG] Wave A done after ~$((i*30))s"; break; }
  sleep 30
done
N=$(find runs/rfd -path "*pd2_egfr_phsw_*" -name 'des_*.pdb' 2>/dev/null | grep -v _traj | wc -l)
echo "[$TAG] backbones to process: $N  $(date -u +%FT%TZ)"
[ "$N" -eq 0 ] && { echo "[$TAG] no backbones; abort"; exit 1; }
bash scripts/run_betA_downstream.sh "$GLOB" "0-7" 2 "$TAG" A
bash scripts/run_betA_switch.sh "$TAG" "$MAXPAE" 48
echo "WAVE17_PIPELINE_DONE $(date -u +%FT%TZ)"
np=$(($(wc -l < "results/rfd_af2_passers_$TAG.csv" 2>/dev/null || echo 1)-1))
no=$(($(wc -l < "results/rfd_obj1_switch_$TAG.csv" 2>/dev/null || echo 1)-1))
echo "[$TAG] af2 passers (pae<$MAXPAE): $np | obj1 switch rows: $no"
