#!/bin/bash
# Bet A downstream, fired once RFdiffusion backbones exist:
#   RFd backbones -> soluble ProteinMPNN (design binder) -> FULL-ATOM-target af2ig inputs -> AF2-ig fold-filter.
# The target MUST be full-atom or af2ig can't template it (plddt_target ~40, 0% yield); build_fullatom_af2in.py
# superposes the real 6ARU DIII full-atom target onto each design's (rigid) target backbone.
# Usage: run_betA_downstream.sh "<rfd dirs glob>" <gpus> <nseq> <tag>
#   e.g. run_betA_downstream.sh 'runs/rfd/w2_*' 1-7 2 w2
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr
cd "$REPO"
GLOB=${1:-runs/rfd/w1_*}
GPUS=${2:-1-7}
NSEQ=${3:-2}
TAG=${4:-w1}
BCHAIN=${5:-B}   # binder chain to redesign: de novo RFd=B, partial diffusion=A
PY="$REPO/envs/md/bin/python"
echo "=== [$TAG] MPNN (design binder chain $BCHAIN, target fixed) $(date -u +%FT%TZ) ==="
"$PY" scripts/rfd_to_af2ig.py --rfd-dirs "$GLOB" --work "runs/rfd_mpnn_$TAG" --gpus "$GPUS" --nseq "$NSEQ" --binder-chain "$BCHAIN" --mpnn-only
echo "=== [$TAG] build FULL-ATOM-target af2ig inputs $(date -u +%FT%TZ) ==="
"$PY" scripts/build_fullatom_af2in.py --mpnn "runs/rfd_mpnn_$TAG/gpu*/backbones/*.pdb" \
    --target data/targets/6ARU_domainIII.pdb --out "runs/rfd_af2in_$TAG"
echo "=== [$TAG] AF2-ig fold-filter $(date -u +%FT%TZ) ==="
"$PY" scripts/run_af2ig.py --pdbdir "runs/rfd_af2in_$TAG" --out "runs/rfd_af2ig_$TAG" --gpus "$GPUS" --mode fast
echo "BETA_DOWNSTREAM_${TAG}_DONE $(date -u +%FT%TZ)"
