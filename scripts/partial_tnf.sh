#!/bin/bash
# Partial-diffuse the best de-novo binder seeds (binder=chain A, target=chain B single protomer)
# to refine the interface (Baker 2-stage protocol; my ~35x lesson). Shards seeds across free GPUs
# with per-GPU cache isolation (avoids the Triton/compile deadlock). Usage: partial_tnf.sh [ndes_per_seed] [partial_T]
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
KIT="$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
unset VIRTUAL_ENV; source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export WEIGHTS="$REPO/weights/rfdiffusion1" RFD_ROOT="$REPO/third_party/RFdiffusion" TMPDIR="${TMPDIR:-/tmp}"
NDES="${1:-10}"; PT="${2:-10}"
OUT="$REPO/runs/partial/out"; mkdir -p "$OUT"
mapfile -t SEEDS < <(ls "$REPO"/runs/partial/seeds/*.pdb)
mapfile -t FREE < <(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | awk -F',' '{gsub(/ /,"");u=$2+0; if($1+0<=2 || u<1000) print $1+0}')
NG=${#FREE[@]}
echo "[partial] ${#SEEDS[@]} seeds, ndes=$NDES partial_T=$PT, GPUs: ${FREE[*]}"
i=0
for s in "${SEEDS[@]}"; do
  g=${FREE[$((i % NG))]}
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
  if [ $((i % NG)) -eq 0 ]; then wait; fi   # drain each GPU wave to bound concurrency to NG
done
wait
echo "[partial] done $(date -u +%FT%TZ); backbones: $(ls "$OUT"/*/pd_*.pdb 2>/dev/null | wc -l)"
