#!/bin/bash
# RFd3 de-novo binder gen for TNF, sharded across GPUs with distinct seeds. Binder=chain A, target=chain B in output CIFs.
# Usage: gen_rfd3_tnf.sh <spec.json> <outdir> <n_batches> <bsize> "<gpus>"
set -uo pipefail
R="$(cd "$(dirname "$0")/.." && pwd)"
SPEC=$1; OUT=$2; NB=${3:-4}; BS=${4:-16}; GPUS=${5:-"0 1 2 3 4 5 6 7"}
KIT=$R/third_party/uplifting-biomolecular-modeling/rfdiffusion3
mkdir -p "$OUT"
echo "[rfd3] spec=$SPEC out=$OUT n_batches=$NB bsize=$BS gpus=[$GPUS] $(date -u +%FT%TZ)"
pids=()
for i in $GPUS; do
  CUDA_VISIBLE_DEVICES=$i PATH=$R/envs/uplifting-rfd3-kit/bin:$PATH \
    RFDIFFUSION3_STOCK_PYTHON=$R/envs/uplifting-rfd3-stock/bin/python \
    RFD3_CKPT=$R/weights/rfd3/rfd3_latest.ckpt \
    MODEL_OPT_JIT_ROOT=$R/.jit/rfd3_g$i TRITON_HOME=$R/.jit/rfd3_g$i TRITON_CACHE_DIR=$R/.jit/rfd3_g$i/cache \
    bash "$KIT/run.sh" design --config h100 --mode fast \
      inputs="$SPEC" out_dir="$OUT/gpu$i" seed=$(( ${SEED_BASE:-0} + 100*(i+1) )) n_batches="$NB" diffusion_batch_size="$BS" \
      > "$OUT/gpu$i.log" 2>&1 &
  pids+=($!); sleep 3
done
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[rfd3] done $(date -u +%FT%TZ); outputs: $(find "$OUT" -name '*.cif.gz' 2>/dev/null | wc -l)"
