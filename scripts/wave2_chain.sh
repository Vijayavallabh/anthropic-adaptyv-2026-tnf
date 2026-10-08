#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")/.."
echo "[w2chain] waiting for RFd3 wave2 gen $(date -u +%T)"
for i in $(seq 1 120); do grep -q '\[rfd3\] done' runs/rfd3/w2_launch.log 2>/dev/null && break; sleep 20; done
NB=$(find runs/rfd3/w2 -name '*.cif.gz' 2>/dev/null|wc -l); echo "[w2chain] $NB RFd3 wave2 designs"
bash scripts/eval_rfd3.sh "runs/rfd3/w2/gpu*/*.cif.gz" w2 "0-7" 8
