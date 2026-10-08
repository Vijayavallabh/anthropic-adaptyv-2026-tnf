#!/bin/bash
# TNF stage: soluble-MPNN sequence design (design binder chain C, 0-Cys) on RFd backbones ->
# full-atom two-protomer af2ig inputs -> AF2-ig binding filter.
# Usage: seq_af2ig_tnf.sh [gpus e.g. 0-7] [nseq_per_backbone] [rfd_glob_subdir]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf
cd "$REPO"
GPUS=${1:-0-7}; NSEQ=${2:-4}; SUB=${3:-g*}
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"  # VANILLA ProteinMPNN: designs interface hydrophobics/aromatics to complement the Trp/Tyr epitope core. soluble_mpnn strips them (42%-charged bundle -> 0 binders)
RFD_OUT="$REPO/runs/rfd_gen"
STAGE="$REPO/runs/seq/staged"; WORK="$REPO/runs/seq/mpnn"; AF2IN="$REPO/runs/seq/af2in"; AF2IG="$REPO/runs/seq/af2ig"
rm -rf "$STAGE" "$WORK" "$AF2IN" "$AF2IG"   # self-clean: stale backbones here would contaminate build_fullatom's glob
mkdir -p "$STAGE" "$WORK"
IFS=- read -r lo hi <<< "$GPUS"; GPULIST=$(seq "$lo" "$hi"); NG=$(echo $GPULIST | wc -w)

# stage every RFd design under a UNIQUE name (des_N repeats across g* dirs otherwise)
n=0
for p in "$RFD_OUT"/$SUB/des_*.pdb; do
  [ -e "$p" ] || continue
  base=$(echo "$p" | sed "s#$RFD_OUT/##; s#/#_#g")   # e.g. g0_L55_des_0.pdb
  ln -sf "$p" "$STAGE/$base"; n=$((n+1))
done
echo "=== [seq] staged $n RFd backbones; soluble-MPNN chain C, 0-Cys, $NSEQ seq/bb, $NG gpus $(date -u +%FT%TZ) ==="
[ "$n" -gt 0 ] || { echo "FATAL: no RFd backbones under $RFD_OUT/$SUB"; exit 1; }
mapfile -t ALL < <(ls "$STAGE"/*.pdb)

# round-robin shard -> per-gpu json list (LigandMPNN list()s the json; a list works)
declare -A LISTS
gi=0
for p in "${ALL[@]}"; do
  g=$(echo $GPULIST | cut -d' ' -f$(( gi % NG + 1 ))); LISTS[$g]+="\"$p\","; gi=$((gi+1))
done
pids=()
for g in $GPULIST; do
  mkdir -p "$WORK/gpu$g"
  printf '[%s]' "${LISTS[$g]%,}" > "$WORK/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py \
      --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$WORK/gpu$g/multi.json" \
      --chains_to_design "C" --bias_AA "C:-100.0" \
      --out_folder "$WORK/gpu$g" \
      --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$g.log" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
NBB=$(find "$WORK" -path '*/backbones/*.pdb' 2>/dev/null | wc -l)
echo "[seq] MPNN done: $NBB designed complexes ($(date -u +%FT%TZ))"
[ "$NBB" -gt 0 ] || { echo "FATAL: no MPNN backbones"; exit 1; }

echo "=== [seq] full-atom two-protomer af2ig inputs $(date -u +%FT%TZ) ==="
uv run python "$REPO/scripts/build_fullatom_af2in_tnf.py" --mpnn "$WORK/gpu*/backbones/*.pdb" \
   --target "$REPO/data/targets/1TNF_AB.pdb" --out "$AF2IN"
NIN=$(find "$AF2IN" -name '*.pdb' | wc -l)
echo "[seq] af2in: $NIN"; [ "$NIN" -gt 0 ] || { echo "FATAL: no af2in inputs"; exit 1; }

echo "=== [seq] AF2-ig binding filter $(date -u +%FT%TZ) ==="
uv run python "$REPO/scripts/run_af2ig.py" --pdbdir "$AF2IN" --out "$AF2IG" --gpus "$GPUS" --mode fast
echo "SEQ_AF2IG_DONE $(date -u +%FT%TZ)"