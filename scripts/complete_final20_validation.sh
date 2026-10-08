#!/bin/bash
# Complete the orthogonal validation (cross-react + Boltz-2 + MD) for the final-20 PD/RFd picks
# that are missing one or more layers, so every submitted design carries all three. GPU plan:
#   GPU0 = baseline MD (already running)   GPUs 1-7 = cross-react af2ig
#   then  GPU1 = Boltz   ||   GPUs 2-7 = MD (2 rep x 1000ps, matches baseline for comparable d_occ)
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
MD="$REPO/envs/md/bin/python"
mkdir -p logs
awk -F'\t' '{print $2}' runs/md_final2_list.txt > runs/md_final2_paths.txt

echo "=== PHASE A: cross-react (mouse + human af2ig), 7 designs, GPUs 1-7 $(date -u +%T) ==="
bash scripts/run_xreact.sh runs/cand_xr2 xr2 1-7 2>&1 | tail -6
echo "XR2_DONE $(date -u +%TZ)"

echo "=== PHASE B: Boltz-2 (GPU1) in background $(date -u +%T) ==="
nohup bash scripts/run_boltz.sh runs/boltz_final2.csv fin2 1 > logs/boltz_fin2.log 2>&1 &
BPID=$!

echo "=== PHASE C: MD 7 designs, GPUs 2-7 $(date -u +%T) ==="
"$MD" scripts/run_md_campaign.py --pdb-list runs/md_final2_paths.txt --out runs/md_final2 \
    --gpus 2-7 --equil-ps 200 --prod-ps 1000 --replicas 2 2>&1 | tail -10
echo "MD_FIN2_DONE $(date -u +%TZ)"

wait $BPID
echo "=== ALL_COMPLETION_DONE $(date -u +%TZ) ==="
echo "MD jsons:    $(ls runs/md_final2/*.json 2>/dev/null | wc -l)/7"
echo "Boltz jsons: $(find runs/boltz_fin2/out -name 'confidence_*.json' 2>/dev/null | wc -l)/7"
echo "XR mouse sc: $(cat runs/af2ig_xreact_xr2_mouse/out.sc runs/af2ig_xreact_xr2_mouse/gpu*/out.sc 2>/dev/null | grep -c '^SCORE:')  human sc: $(cat runs/af2ig_xreact_xr2_human/out.sc runs/af2ig_xreact_xr2_human/gpu*/out.sc 2>/dev/null | grep -c '^SCORE:')"