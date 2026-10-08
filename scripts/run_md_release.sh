#!/bin/bash
# Two-state MD release-switch gate over a manifest (design_id<TAB>his_pos<TAB>anchor<TAB>pdb),
# sharded design-by-design across GPUs. Each md_release_tnf.py runs its own HIP/HIE x replicas.
# Usage: run_md_release.sh <manifest.tsv> <outdir> "<gpu list>" [prod_ps] [replicas]
set -uo pipefail; cd /home/prachh/v/anthropic-adaptyv-2026-tnf
MAN=$1; OUT=$2; GPUS=${3:-"0 1 2 3 4 5 6 7"}; PROD=${4:-3000}; REP=${5:-3}
PY=envs/md2/bin/python   # cuda12 build (CUDA platform works on this H100)
mkdir -p "$OUT"
read -ra GL <<< "$GPUS"; NG=${#GL[@]}
# round-robin assign designs to GPUs
declare -A Q; i=0
while IFS=$'\t' read -r did hp anc pdb; do
  [ -z "$did" ] && continue
  g=${GL[$((i%NG))]}; Q[$g]+="$did|$hp|$pdb"$'\n'; i=$((i+1))
done < "$MAN"
echo "[md] $i designs over GPUs [$GPUS] prod=${PROD}ps rep=$REP $(date -u +%FT%TZ)"
pids=()
for g in "${GL[@]}"; do
  [ -z "${Q[$g]:-}" ] && continue
  ( while IFS='|' read -r did hp pdb; do
      [ -z "$did" ] && continue
      echo "  gpu$g <- $did (his A$hp)"
      CUDA_VISIBLE_DEVICES=$g $PY scripts/md_release_tnf.py --pdb "$pdb" --his-pos "$hp" \
        --out "$OUT" --device 0 --prod-ps "$PROD" --replicas "$REP" \
        > "$OUT/${did}.log" 2>&1
    done <<< "${Q[$g]}" ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[md] done $(date -u +%FT%TZ); json: $(ls "$OUT"/*.json 2>/dev/null|wc -l)"
