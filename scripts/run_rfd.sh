#!/bin/bash
# RFdiffusion-1 binder design driver (uplifting kit, fast mode, one GPU).
# Usage: run_rfd.sh <out_subdir> <len_spec> <hotspots_csv> <num_designs> <gpu_id> [startnum] [sched_dir]
#   len_spec   : e.g. "80-80" (fixed) or "65-110" (range)  -> binder length after the A310-500/0 target contig
#   hotspots   : PDB author numbering, e.g. "A408,A409,A412,A440,A465"  (A409 = His433)
# Writes runs/rfd/<out_subdir>/des_<i>.pdb (+ .trb, run.log).  Absolute paths (run.sh cd's into kit dir).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUTSUB=$1; LENSPEC=$2; HOTS=$3; NDES=$4; GPU=$5; START=${6:-0}; SCHED=${7:-}
source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export RFD_ROOT="$REPO/third_party/RFdiffusion"
export WEIGHTS="$REPO/weights/rfdiffusion1"
export MODEL_OPT_JIT_ROOT="$REPO/.jit"
export CUDA_VISIBLE_DEVICES=$GPU
INP="$REPO/data/targets/6ARU_domainIII.pdb"
OUTDIR="$REPO/runs/rfd/$OUTSUB"; mkdir -p "$OUTDIR"
EXTRA=()
[ -n "$SCHED" ] && { mkdir -p "$SCHED"; SCHED=$(readlink -f "$SCHED"); EXTRA+=("inference.schedule_directory_path=$SCHED"); }
cd "$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
exec bash run.sh design --config h100 --mode fast \
  inference.input_pdb="$INP" \
  "contigmap.contigs=[A310-500/0 $LENSPEC]" \
  "ppi.hotspot_res=[$HOTS]" \
  inference.num_designs="$NDES" \
  inference.design_startnum="$START" \
  "${EXTRA[@]}" \
  inference.output_prefix="$OUTDIR/des"
