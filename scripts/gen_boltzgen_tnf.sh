#!/bin/bash
set -uo pipefail
R="$(cd "$(dirname "$0")/.." && pwd)"
GPU=$1; MODE=$2; SPEC=$3; OUT=$4; NDES=$5; SEED=$6; BATCH=$7; shift 7
unset VIRTUAL_ENV 2>/dev/null || true
source "$R/envs/uplifting-boltzgen/bin/activate"
[ -f "$R/envs/uplifting.env" ] && { set -a; source "$R/envs/uplifting.env"; set +a; }
export CUDA_VISIBLE_DEVICES=$GPU MODEL_OPT_JIT_ROOT="$R/.jit/bg_g${GPU}" TRITON_HOME="$R/.jit/bg_g${GPU}" TRITON_CACHE_DIR="$R/.jit/bg_g${GPU}/cache"
mkdir -p "$OUT" "$R/.jit/bg_g${GPU}"
cd "$R/third_party/uplifting-biomolecular-modeling/boltzgen" || exit 1
echo "[bg] gpu=$GPU mode=$MODE ndes=$NDES seed=$SEED spec=$SPEC out=$OUT start=$(date -u +%FT%TZ)"
bash run.sh design --config h100 --mode "$MODE" "$SPEC" --output "$OUT" --num_designs "$NDES" --seed "$SEED" --diffusion_batch_size "$BATCH" "$@"
echo "[bg] DONE rc=$? end=$(date -u +%FT%TZ)"
