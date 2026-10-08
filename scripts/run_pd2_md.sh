#!/bin/bash
# Wait for the de novo w2 generation to free the GPUs, then MD the diverse pd2 shortlist so those
# candidates can dual-confirm (His433 PROPKA + MD) and compete in the family-capped selection.
# Additive: run_md_campaign skips designs already scored in md_final2 (the 7 floor picks).
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
echo "waiting for de novo w2 gen to finish... $(date -u +%T)"
for i in $(seq 1 160); do
  grep -q RFD_CAMPAIGN_w2_DONE logs/rfd/w2_campaign.log 2>/dev/null && { echo "de novo gen done after ~$((i*30))s"; break; }
  sleep 30
done
echo "=== MD diverse pd2 shortlist (2 rep x 1000ps, GPUs 0-7) $(date -u +%T) ==="
envs/md/bin/python scripts/run_md_campaign.py --pdb-list runs/md_pd2_paths.txt --out runs/md_final2 \
    --gpus 0-2 --equil-ps 200 --prod-ps 1000 --replicas 2 2>&1 | tail -10
echo "PD2_MD_DONE $(date -u +%TZ)"
echo "md_final2 jsons: $(ls runs/md_final2/*.json 2>/dev/null | grep -v gpu | wc -l)"