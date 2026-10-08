#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")/.."
uv run python scripts/run_af2ig.py --pdbdir runs/seq/af2in_1p --out runs/seq/af2ig_1p --gpus 0-3 --mode fast > runs/af2ig_1p.log 2>&1
echo "=== SINGLE-PROTOMER PASS RATE ===" >> runs/af2ig_1p.log
awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++} END{printf "scored=%d PASS(pae<10,plddt_b>80,rmsd<2)=%d pae<10any=%d pae<15=%d\n",n,p+0,pl+0,p15+0}' runs/seq/af2ig_1p/out.sc >> runs/af2ig_1p.log
echo "--- top 15 by pae_interaction ---" >> runs/af2ig_1p.log
(head -1 runs/seq/af2ig_1p/out.sc; awk 'NR>1 && $4!="nan"' runs/seq/af2ig_1p/out.sc | sort -k4 -n | head -15) >> runs/af2ig_1p.log
echo "AF2IG_1P_DONE" >> runs/af2ig_1p.log
