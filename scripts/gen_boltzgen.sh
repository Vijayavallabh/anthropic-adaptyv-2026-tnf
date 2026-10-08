#!/bin/bash
# Run ONE BoltzGen design job on ONE GPU, with the isolated venv + weights cache wired up.
# Full pipeline by default (design, inverse_folding, design_folding, folding, analysis, filtering);
# pass e.g. `-- --steps design inverse_folding` after the fixed args to shorten it.
#
# usage: gen_boltzgen.sh GPU MODE SPEC OUTDIR NDES SEED BATCH [extra boltzgen args...]
#   MODE = fast | exact | big | off
set -uo pipefail
R="$(cd "$(dirname "$0")/.." && pwd)"
if [ $# -lt 7 ]; then echo "usage: gen_boltzgen.sh GPU MODE SPEC OUTDIR NDES SEED BATCH [extra...]" >&2; exit 2; fi
GPU=$1; MODE=$2; SPEC=$3; OUT=$4; NDES=$5; SEED=$6; BATCH=$7; shift 7
unset VIRTUAL_ENV 2>/dev/null || true
# shellcheck disable=SC1091
source "$R/envs/uplifting-boltzgen/bin/activate"
set -a; source "$R/envs/uplifting.env"; set +a
export CUDA_VISIBLE_DEVICES=$GPU
export MODEL_OPT_JIT_ROOT="$R/.jit"
mkdir -p "$OUT"
cd "$R/third_party/uplifting-biomolecular-modeling/boltzgen" || exit 1
echo "[gen] gpu=$GPU mode=$MODE ndes=$NDES seed=$SEED batch=$BATCH spec=$SPEC"
echo "[gen] out=$OUT start=$(date -u +%FT%TZ)"
bash run.sh design --config h100 --mode "$MODE" "$SPEC" \
  --output "$OUT" --num_designs "$NDES" --seed "$SEED" --diffusion_batch_size "$BATCH" "$@"
rc=$?
echo "[gen] DONE rc=$rc end=$(date -u +%FT%TZ)"
exit $rc
