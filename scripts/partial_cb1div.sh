#!/bin/bash
# Diversify the 5 MD-CONFIRMED cb1 binders by RFd1 partial diffusion (proven cb1 recipe: noise_scale=0,
# gentle partial_T) to deepen the cb1-family candidate pool with backbone diversity. Binder=chain A,
# target=chain B (human TNF, 152). Output backbones -> runs/partial_cb1div/out/<seed>/pd_*.pdb
# Usage: partial_cb1div.sh [ndes_per_seed] [partial_T] "<gpus>"
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
KIT="$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
unset VIRTUAL_ENV; source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export WEIGHTS="$REPO/weights/rfdiffusion1" RFD_ROOT="$REPO/third_party/RFdiffusion" TMPDIR=/home/prachh/work/tmp
NDES="${1:-24}"; PT="${2:-12}"; read -ra GPUS <<< "${3:-4 5 6 7}"
OUT="$REPO/runs/partial_cb1div/out"; mkdir -p "$OUT"
mapfile -t SEEDS < <(ls "$REPO"/runs/partial/seeds_cb1div/*.pdb)
NG=${#GPUS[@]}
echo "[cb1div] ${#SEEDS[@]} seeds, ndes=$NDES partial_T=$PT, GPUs: ${GPUS[*]} $(date -u +%FT%TZ)"
i=0
for s in "${SEEDS[@]}"; do
  g=${GPUS[$((i % NG))]}
  Lb=$(awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="A" && substr($0,13,4)~/CA/{n++}END{print n}' "$s")
  base=$(basename "$s" .pdb); pre="$OUT/${base}"; mkdir -p "$pre/sched" "$REPO/.jit/triton_g${g}"
  ( cd "$KIT" && CUDA_VISIBLE_DEVICES="$g" OMP_NUM_THREADS=24 MKL_NUM_THREADS=24 \
       MODEL_OPT_JIT_ROOT="$REPO/.jit/Lb${Lb}" TRITON_CACHE_DIR="$REPO/.jit/triton_g${g}/cache" TRITON_HOME="$REPO/.jit/triton_g${g}" \
       bash run.sh design --config h100 --mode fast \
       inference.input_pdb="$s" \
       "contigmap.contigs=[${Lb}-${Lb}/0 B1-152]" \
       diffuser.partial_T="$PT" \
       denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
       inference.num_designs="$NDES" \
       inference.schedule_directory_path="$pre/sched" \
       inference.output_prefix="$pre/pd" ) > "$OUT/log_${base}.log" 2>&1 &
  i=$((i+1))
  if [ $((i % NG)) -eq 0 ]; then wait; fi
done
wait
echo "[cb1div] done $(date -u +%FT%TZ); backbones: $(ls "$OUT"/*/pd_*.pdb 2>/dev/null | wc -l)"
