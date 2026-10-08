#!/bin/bash
set -uo pipefail; cd /home/prachh/v/anthropic-adaptyv-2026-tnf
SEED_BASE=40000 bash scripts/gen_rfd3_tnf.sh "$(readlink -f specs/rfd3/tnf_K90.json)" runs/rfd3/w4 16 16 "4 5 6 7" > runs/rfd3/w4_gen.log 2>&1
bash scripts/eval_rfd3.sh "runs/rfd3/w4/gpu*/*.cif.gz" w4 "4-7" 8
