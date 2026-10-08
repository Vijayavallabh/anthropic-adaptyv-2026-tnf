#!/bin/bash
# Complete the rfdiffusion1 uplifting kit (Route C): the kit's opt wheel is installed, but it
# wraps upstream RFdiffusion, which must be cloned at the pinned commit, pip -e installed into
# the kit venv, and have its weights fetched. setup_uplifting.sh does not do this.
set -uo pipefail
cd "$(dirname "$0")/.."
export TMPDIR=/home/prachh/work/tmp MODEL_OPT_JIT_ROOT="$PWD/.jit"
export PYTHONHASHSEED=0 CFLAGS=-g0 CC=gcc CXX=g++
PIN=86507b6538f51fce57b5a72477165f03999ed7ae
KIT="$PWD/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
VENV="$PWD/envs/uplifting-rfdiffusion1"
RFD="$PWD/third_party/RFdiffusion"
WEIGHTS="$PWD/weights/rfdiffusion1"
mkdir -p "$WEIGHTS"

echo "=== [1/4] clone upstream RFdiffusion @ $PIN ==="
[ -d "$RFD/.git" ] || git clone https://github.com/RosettaCommons/RFdiffusion "$RFD"
git -C "$RFD" checkout "$PIN"

echo "=== [2/4] pip -e SE3Transformer + RFdiffusion into the kit venv ==="
# shellcheck disable=SC1091
source "$VENV/bin/activate"
pip install --no-deps -e "$RFD/env/SE3Transformer" -e "$RFD"

echo "=== [3/4] run.sh install --weights ==="
export RFD_ROOT="$RFD" WEIGHTS="$WEIGHTS"
( cd "$KIT" && bash run.sh install --weights "$WEIGHTS" )

echo "=== [4/4] run.sh check (RFD_ROOT + WEIGHTS set) ==="
( cd "$KIT" && RFD_ROOT="$RFD" WEIGHTS="$WEIGHTS" bash run.sh check --config h100 --mode fast ) \
  && echo "=== rfd1 ACTIVE ===" || echo "=== rfd1 check nonzero (see above) ==="
# record the runtime env vars for the pipeline
{ echo "export RFD_ROOT=$RFD"; echo "export WEIGHTS_RFD=$WEIGHTS"; } >> "$PWD/envs/uplifting.env"