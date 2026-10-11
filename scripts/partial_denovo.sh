#!/bin/bash
# LONG-SHOT amplification: partial-diffuse the 3 cross-reactive DE-NOVO shapes (des_138_7/4082_8/2044_8)
# -> backbone-diverse cousins of NEW lineages (retains binding mode; high partial yield). 3 seeds x
# {T10,T20} = 6 jobs across the 6 free GPUs (1,3,4,5,6,7). Per-GPU triton isolation. noise_scale_ca=0.
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
KIT="$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
unset VIRTUAL_ENV; source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export WEIGHTS="$REPO/weights/rfdiffusion1" RFD_ROOT="$REPO/third_party/RFdiffusion"
PD_BASE="${PD_BASE:-runs/partial_denovo}"
OUT="$REPO/$PD_BASE/out"; mkdir -p "$OUT"
NDES="${1:-96}"
IFS=' ' read -ra GPUS <<< "${PD_GPUS:-1 3 4 5 6 7}"
NG=${#GPUS[@]}
mapfile -t SEEDS < <(ls "$REPO"/$PD_BASE/seeds/*.pdb)
JOBS=()
for s in "${SEEDS[@]}"; do for T in 10 20; do JOBS+=("$s|$T"); done; done
echo "[dn-partial] ${#SEEDS[@]} seeds x2 T = ${#JOBS[@]} jobs, ndes=$NDES, GPUs: ${GPUS[*]} $(date -u +%FT%TZ)"
i=0
for job in "${JOBS[@]}"; do
  s="${job%|*}"; T="${job#*|}"; g=${GPUS[$((i % NG))]}
  base=$(basename "$s" .pdb)
  Lb=$(awk 'substr($0,1,4)=="ATOM"&&substr($0,22,1)=="A"&&substr($0,13,4)~/CA/{n++}END{print n}' "$s")
  Lt=$(awk 'substr($0,1,4)=="ATOM"&&substr($0,22,1)=="B"&&substr($0,13,4)~/CA/{n++}END{print n}' "$s")
  pre="$OUT/${base}_T${T}"; mkdir -p "$pre/sched" "$REPO/.jit/triton_g${g}"
  ( cd "$KIT" && CUDA_VISIBLE_DEVICES="$g" OMP_NUM_THREADS=24 MKL_NUM_THREADS=24 OPENBLAS_NUM_THREADS=24 \
      MODEL_OPT_JIT_ROOT="$REPO/.jit/Lb${Lb}" TRITON_CACHE_DIR="$REPO/.jit/triton_g${g}/cache" TRITON_HOME="$REPO/.jit/triton_g${g}" \
      bash run.sh design --config h100 --mode fast \
      inference.input_pdb="$s" \
      "contigmap.contigs=[${Lb}-${Lb}/0 B1-${Lt}]" \
      diffuser.partial_T="$T" \
      denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
      inference.num_designs="$NDES" \
      inference.schedule_directory_path="$pre/sched" \
      inference.output_prefix="$pre/pd" ) > "$OUT/log_${base}_T${T}.log" 2>&1 &
  i=$((i+1))
  [ $((i % NG)) -eq 0 ] && wait
done
wait
echo "[dn-partial] DONE $(date -u +%FT%TZ); backbones: $(ls "$OUT"/*/pd_*.pdb 2>/dev/null | wc -l)"
