#!/bin/bash
# Independent MONOMER folder for the expression-consensus gate: Boltz-2 (AF3-style diffusion,
# architecturally independent of ESMFold's ESM2 backbone) folds each binder ALONE, single-sequence
# (msa: empty), giving an orthogonal monomer pLDDT + pTM. Cross-folder TM(ESMFold,Boltz) agreement is
# the Baker/Rocklin field-standard expression signal (two independent folders must agree on the fold).
# Usage: fold_monomer_boltz.sh <input_yaml_dir> <out_dir> <gpu>   (GPUs 3,4,5 ONLY)
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
IN=$1; OUT=$2; GPU=$3
# absolutize BEFORE we cd into the boltz dir (relative paths would resolve from there -> 0 inputs)
case "$IN" in /*) ;; *) IN="$REPO/$IN";; esac
case "$OUT" in /*) ;; *) OUT="$REPO/$OUT";; esac
mkdir -p "$OUT"
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
unset VIRTUAL_ENV 2>/dev/null || true
source "$REPO/envs/uplifting-boltz2/bin/activate"
set -a; source "$REPO/envs/uplifting.env"; set +a
cd "$REPO/third_party/uplifting-biomolecular-modeling/boltz2"
echo "[boltz-mono] gpu=$GPU input=$IN out=$OUT $(date -u +%FT%TZ)"
bash run.sh pred --config h100 --mode fast --input "$IN" --out_dir "$OUT"
echo "BOLTZ_MONO_DONE gpu=$GPU $(date -u +%FT%TZ) | structures: $(find "$OUT" -name '*_model_0.cif' 2>/dev/null | wc -l)"
