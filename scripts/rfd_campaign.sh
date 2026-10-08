#!/bin/bash
# Fan RFdiffusion-1 binder design across GPUs, one process per GPU.
# Usage: rfd_campaign.sh <wave> <ndes_per_gpu> <len_spec> "<gpu list>" [hotspots_override]
#   e.g. rfd_campaign.sh w2 300 60-110 "0 1 2 3 4 5 6 7"
# GPUs other than the last use the extended His433 conserved-face patch; the last uses the tight patch
# (unless hotspots_override is given, which applies to all).  His433 = A409.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WAVE=${1:-w1}; NDES=${2:-300}; LENSPEC=${3:-60-110}; GPUS=${4:-"0 1 2 3 4 5 6 7"}; OVR=${5:-}
EXT="A408,A409,A412,A440,A465"; TIGHT="A408,A409,A412"
mkdir -p "$REPO/logs/rfd"
last=$(echo $GPUS | tr ' ' '\n' | tail -1)
for g in $GPUS; do
  if [ -n "$OVR" ]; then HS="$OVR"; SUB="${WAVE}_g${g}"
  elif [ "$g" = "$last" ]; then HS="$TIGHT"; SUB="${WAVE}_tight_g${g}"
  else HS="$EXT"; SUB="${WAVE}_ext_g${g}"; fi
  LOG="$REPO/logs/rfd/${SUB}.log"
  echo "GPU $g -> runs/rfd/$SUB (hotspots [$HS], $NDES designs, len $LENSPEC)"
  nohup bash "$REPO/scripts/run_rfd.sh" "$SUB" "$LENSPEC" "$HS" "$NDES" "$g" 0 "$REPO/runs/rfd/sched_${SUB}" > "$LOG" 2>&1 &
done
wait
echo "RFD_CAMPAIGN_${WAVE}_DONE"
