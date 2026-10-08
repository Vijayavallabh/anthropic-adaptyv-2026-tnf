#!/bin/bash
set -uo pipefail; cd "$(dirname "$0")/.."
SEED_BASE=20000 bash scripts/gen_rfd3_tnf.sh "$(readlink -f specs/rfd3/tnf_K90.json)" runs/rfd3/w3 16 16 "0 1 2 3 4 5 6 7" > runs/rfd3/w3_gen.log 2>&1
bash scripts/eval_rfd3.sh "runs/rfd3/w3/gpu*/*.cif.gz" w3 "0-7" 8
