#!/bin/bash
# beta-pairing targeted RFdiffusion (Sappington/Baker): condition the binder to form a strand that
# beta-pairs with TNF exposed edge strand B114-120 (Region I, conserved -> backbone-directed x-react).
# Usage: run_bpair_rfd.sh <out_subdir> <num_designs> <gpu> [startnum]
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUTSUB=$1; NDES=$2; GPU=$3; START=${4:-0}
source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export RFD_ROOT="$REPO/third_party/RFdiffusion" WEIGHTS="$REPO/weights/rfdiffusion1"
export MODEL_OPT_JIT_ROOT="$REPO/.jit" CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=${RFD_CPU_THREADS:-24} MKL_NUM_THREADS=${RFD_CPU_THREADS:-24} OPENBLAS_NUM_THREADS=${RFD_CPU_THREADS:-24}
INP="$REPO/data/targets/1TNF_B.pdb"
TENS="$REPO/runs/bpair/tensors"; TLIST="$REPO/runs/bpair/tensors.txt"
OUTDIR="$REPO/runs/bpair/$OUTSUB"; mkdir -p "$OUTDIR"
cd "$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
# --mode off: fold-conditioning runs on upstream scripts/run_inference.py (the kit's fast driver refuses it).
# ckpt_override_path points at the fold-conditioning checkpoint (outside the kit's pins, downloaded manually).
exec bash run.sh design --config h100 --mode off \
  inference.input_pdb="$INP" \
  inference.ckpt_override_path="$WEIGHTS/Complex_Fold_base_ckpt.pt" \
  "contigmap.contigs=[80-80/0 B6-157]" \
  "ppi.hotspot_res=[B114,B115,B118,B119]" \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$TENS" \
  scaffoldguided.scaffold_list="$TLIST" \
  inference.num_designs="$NDES" \
  inference.design_startnum="$START" \
  inference.output_prefix="$OUTDIR/des"
