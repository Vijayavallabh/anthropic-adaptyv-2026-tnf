#!/bin/bash
# Expression rescue: soluble-MPNN redesign of chain A SURFACE (interface+His frozen), 0-Cys.
# Reads runs/push/rescue/manifest.tsv (design_id, backbone_pdb, chain, fixed_residues, his_pos, seq).
# Usage: run_rescue_mpnn.sh [n_batches=8] [batch_size=4] [temp=0.25]
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf
L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/solublempnn_v_48_020.pt"
OUT="$R/runs/push/rescue/mpnn"; mkdir -p "$OUT"
NB=${1:-8}; BS=${2:-4}; T=${3:-0.25}
GPUS=(0 1 2 3 4 5 6 7); gi=0; pids=()
while IFS=$'\t' read -r did pdb chain fixed his seq; do
  [ "$did" = "design_id" ] && continue
  [ -z "${fixed// }" ] && { echo "skip $did (no fixed)"; continue; }
  g=${GPUS[$((gi % 8))]}; gi=$((gi+1))
  mkdir -p "$OUT/$did"
  echo "[$did] gpu$g  fixed=$(echo $fixed | wc -w) his=$his"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py \
      --model_type soluble_mpnn --checkpoint_soluble_mpnn "$CKPT" \
      --pdb_path "$pdb" --out_folder "$OUT/$did" \
      --chains_to_design "$chain" --fixed_residues "$fixed" \
      --bias_AA "C:-100.0" --number_of_batches "$NB" --batch_size "$BS" \
      --temperature "$T" ) > "$OUT/$did.log" 2>&1 &
  pids+=($!)
done < "$R/runs/push/rescue/manifest.tsv"
echo "launched ${#pids[@]} soluble-MPNN jobs; waiting..."
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "rescue MPNN done (failed=$fail) $(date -u +%FT%TZ)"
