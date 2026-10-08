#!/bin/bash
# Parallel RFdiffusion binder-generation sweep for the TNF-alpha A-B receptor-binding groove.
# Shards across free GPUs (reserved 0,1,2 + any idle others), minibinder-first length spread.
# Each shard gets its own schedule dir (RFD1-002) and output prefix. Usage:
#   scripts/run_gen_tnf.sh [designs_per_shard] [hotspot_set]
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
KIT="$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
unset VIRTUAL_ENV
# shellcheck disable=SC1091
source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export WEIGHTS="$REPO/weights/rfdiffusion1" RFD_ROOT="$REPO/third_party/RFdiffusion"
export TMPDIR=/home/prachh/work/tmp   # MODEL_OPT_JIT_ROOT set per-length below (shared root deadlocks parallel shards on its compile lock)

N="${1:-60}"                                   # designs per shard
HOTSPOTS="${2:-[B86,B87,B88,B90,B91]}"          # BEST bindable conserved+EXPOSED patch (proper CB-nbr scan): Tyr87 aromatic + Val91 hydrophobic (de-novo anchors) + Lys90 (switch/cross-react anchor); all conserved h/mouse, CA-span 11.8A. (prev [B114,B115,B141,B65] had buried aromatics -> pae floor 21.7)
TARGET="$REPO/data/targets/1TNF_AB.pdb"
OUT="$REPO/runs/rfd_gen"; mkdir -p "$OUT"
if [ -n "${LENS_OVERRIDE:-}" ]; then read -r -a LENS <<< "$LENS_OVERRIDE"; else LENS=(55 60 65 70 75 80 85 95); fi  # minibinder sweet-spot; override for longer binders on a flat target

# GPUs: honor GPU_OVERRIDE (space/comma list) if set, else 0,1,2 + any idle others (<1 GiB, no foreign job)
if [ -n "${GPU_OVERRIDE:-}" ]; then
  mapfile -t FREE < <(echo "$GPU_OVERRIDE" | tr ', ' '\n' | grep -E '^[0-9]+$')
else
  mapfile -t FREE < <(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
    | awk -F',' '{gsub(/ /,"");u=$2+0; if($1+0<=2 || u<1000) print $1+0}')
fi
echo "[gen] GPUs: ${FREE[*]}  designs/shard=$N  hotspots=$HOTSPOTS"
i=0
for g in "${FREE[@]}"; do
  L="${LENS[$((i % ${#LENS[@]}))]}"
  pre="$OUT/g${g}_L${L}"; mkdir -p "$pre/sched"
  mkdir -p "$REPO/.jit/triton_g${g}"
  ( cd "$KIT" && CUDA_VISIBLE_DEVICES="$g" MODEL_OPT_JIT_ROOT="$REPO/.jit/L${L}" \
       TRITON_CACHE_DIR="$REPO/.jit/triton_g${g}/cache" TRITON_HOME="$REPO/.jit/triton_g${g}" \
       bash run.sh design --config h100 --mode fast \
       inference.input_pdb="$TARGET" \
       "contigmap.contigs=[A6-157/0 B6-157/0 ${L}-${L}]" \
       "ppi.hotspot_res=$HOTSPOTS" \
       denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
       inference.num_designs="$N" \
       inference.schedule_directory_path="$pre/sched" \
       inference.output_prefix="$pre/des" ) > "$OUT/gen_g${g}_L${L}.log" 2>&1 &
  echo "[gen] GPU $g len $L -> $pre (pid $!)"
  i=$((i+1)); sleep 3
done
wait
echo "[gen] sweep done $(date -u +%FT%TZ); pdbs: $(ls "$OUT"/g*/des_*.pdb 2>/dev/null | wc -l)"