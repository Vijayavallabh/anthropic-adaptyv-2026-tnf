#!/bin/bash
# SolubleMPNN-design PXDesign binder backbones (chain A), TNF trimer (B,C,D) as fixed context.
# 0-Cys (expression), negative surface. Reads runs/pxdesign/stage/manifest.csv.
# Usage: pxd_mpnn.sh [nseq=2] [temp=0.15] "gpu list" (e.g. "3 4 5 6 7")
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf
L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/solublempnn_v_48_020.pt"
OUT="$R/runs/pxdesign/mpnn"; mkdir -p "$OUT"
NS=${1:-2}; T=${2:-0.15}; GPUS=(${3:-3 4 5 6 7}); n=${#GPUS[@]}; gi=0; pids=()
while IFS=, read -r did pdb chain length ncontact task; do
  [ "$did" = "design_id" ] && continue
  [ -z "$pdb" ] && continue
  g=${GPUS[$((gi % n))]}; gi=$((gi+1)); mkdir -p "$OUT/$did"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py \
      --model_type soluble_mpnn --checkpoint_soluble_mpnn "$CKPT" \
      --pdb_path "$pdb" --out_folder "$OUT/$did" \
      --chains_to_design "$chain" \
      --bias_AA "C:-100.0" --number_of_batches "$NS" --batch_size 1 \
      --temperature "$T" --omit_AA "C" ) > "$OUT/$did.log" 2>&1 &
  pids+=($!)
  # throttle: keep at most n*2 in flight
  if [ $(( gi % (n*2) )) -eq 0 ]; then wait "${pids[@]}" 2>/dev/null; pids=(); fi
done < "$R/runs/pxdesign/stage/manifest.csv"
wait 2>/dev/null
echo "PXD_MPNN_DONE $(date -u +%FT%TZ) seqs=$(find "$OUT" -name '*.fa' | wc -l)"
