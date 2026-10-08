#!/bin/bash
# Larger campaign: LONGER binders (90/120/140) on the best epitope (K90/Y87/V91), GPUs 0,1,2.
# gen -> vanilla MPNN (chain C) -> single-protomer af2ig -> pass rate. Harvests any pae<15.
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
mv runs/rfd_gen runs/rfd_gen_e2_$(date +%s) 2>/dev/null || true
echo "[long] gen longer binders $(date -u +%T)"
GPU_OVERRIDE="0 1 2" LENS_OVERRIDE="90 120 140" bash scripts/run_gen_tnf.sh 150 > runs/gen_long.log 2>&1
NB=$(ls runs/rfd_gen/g*/des_*.pdb 2>/dev/null|wc -l); echo "[long] $NB backbones"
uv run python scripts/rawcontacts2.py "runs/rfd_gen/g*/des_*.pdb" 2>&1 | grep -v VIRTUAL_ENV || true
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
STAGE="$REPO/runs/long/staged"; WORK="$REPO/runs/long/mpnn"; AF2IN="$REPO/runs/long/af2in"; AF2IG="$REPO/runs/long/af2ig"
rm -rf "$REPO/runs/long"; mkdir -p "$STAGE" "$WORK"
n=0; for p in runs/rfd_gen/g*/des_*.pdb; do [ -e "$p" ]||continue; b=$(echo "$p"|sed 's#runs/rfd_gen/##;s#/#_#g'); ln -sf "$REPO/$p" "$STAGE/$b"; n=$((n+1)); done
GL="0 1 2"; NG=3; mapfile -t ALL < <(ls "$STAGE"/*.pdb)
declare -A LISTS; gi=0
for p in "${ALL[@]}"; do g=$(echo $GL|cut -d' ' -f$((gi%NG+1))); LISTS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in $GL; do mkdir -p "$WORK/gpu$g"; printf '[%s]' "${LISTS[$g]%,}" > "$WORK/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$WORK/gpu$g/multi.json" --chains_to_design "C" --bias_AA "C:-100.0" \
      --out_folder "$WORK/gpu$g" --batch_size 6 --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
echo "[long] MPNN done: $(find "$WORK" -path '*/backbones/*.pdb'|wc -l)"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
   --target data/targets/1TNF_AB.pdb --target-chains B --binder-chain C --out "$AF2IN" 2>&1 | grep -v VIRTUAL_ENV | tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus 0-2 --mode fast > runs/long/af2ig_run.log 2>&1
echo "=== LONG-BINDER CAMPAIGN PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++; if($4+0<18)p18++} END{printf "scored=%d PASS=%d pae<10=%d pae<15=%d pae<18=%d\n",n,p+0,pl+0,p15+0,p18+0}' "$AF2IG/out.sc"
(head -1 "$AF2IG/out.sc"; awk 'NR>1 && $4!="nan"' "$AF2IG/out.sc"|sort -k4 -n|head -15) 
echo "LONG_CAMPAIGN_DONE"
