#!/bin/bash
# Partial-diffuse the 5 VALIDATED binders (runs/partial/seeds), 80 designs/seed, partial_T=10, then MPNN+af2ig.
# Expected ~11% af2ig pass (validated-binder partial diffusion) = a confident-binder FAMILY.
set -uo pipefail; cd /home/prachh/v/anthropic-adaptyv-2026-tnf
bash scripts/partial_tnf.sh 80 10 > runs/partial_tnf.log 2>&1
echo "[fam] partial gen done; backbones: $(ls runs/partial/out/*/pd_*.pdb 2>/dev/null|wc -l)"
bash scripts/post_partial.sh 0-7 8 > runs/partial_post.log 2>&1
echo "=== PARTIAL-FAMILY PASS RATE ===" 
sed -n '/PARTIAL-DIFFUSION PASS RATE/,$p' runs/partial_post.log 2>/dev/null | grep -viE 'oneDNN|LEVER|deprecat|WARNING:tensor|Instructions|Use fn' | head -20
echo PARTIAL_FAM_DONE
