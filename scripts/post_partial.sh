#!/bin/bash
# After partial diffusion (binder=chain A, target=chain B): vanilla-MPNN design chain A (context B, 0-Cys),
# rebuild full-atom target B, af2ig, report pass rate. Usage: post_partial.sh [gpus 0-3] [nseq]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
GPUS=${1:-0-3}; NSEQ=${2:-8}
echo "[pp] waiting for partial diffusion to finish $(date -u +%T)"
for i in $(seq 1 180); do grep -q '\[partial\] done' runs/partial_tnf.log 2>/dev/null && break; sleep 20; done
NB=$(ls runs/partial/out/*/pd_*.pdb 2>/dev/null | wc -l); echo "[pp] $NB partial backbones"
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
STAGE="$REPO/runs/partial/staged"; WORK="$REPO/runs/partial/mpnn"; AF2IN="$REPO/runs/partial/af2in"; AF2IG="$REPO/runs/partial/af2ig"
rm -rf "$STAGE" "$WORK" "$AF2IN" "$AF2IG"; mkdir -p "$STAGE" "$WORK"
n=0; for p in runs/partial/out/*/pd_*.pdb; do [ -e "$p" ]||continue; b=$(echo "$p"|sed 's#runs/partial/out/##;s#/#_#g'); ln -sf "$REPO/$p" "$STAGE/$b"; n=$((n+1)); done
echo "[pp] staged $n"
IFS=- read lo hi <<< "$GPUS"; GL=$(seq $lo $hi); NG=$(echo $GL|wc -w)
mapfile -t ALL < <(ls "$STAGE"/*.pdb)
declare -A LISTS; gi=0
for p in "${ALL[@]}"; do g=$(echo $GL|cut -d' ' -f$((gi%NG+1))); LISTS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in $GL; do mkdir -p "$WORK/gpu$g"; printf '[%s]' "${LISTS[$g]%,}" > "$WORK/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$WORK/gpu$g/multi.json" --chains_to_design "A" --bias_AA "C:-100.0" \
      --out_folder "$WORK/gpu$g" --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
NM=$(find "$WORK" -path '*/backbones/*.pdb'|wc -l); echo "[pp] MPNN done: $NM"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
   --target data/targets/1TNF_AB.pdb --target-chains B --binder-chain A --out "$AF2IN" 2>&1 | grep -v VIRTUAL_ENV | tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus "$GPUS" --mode fast > runs/partial/af2ig_run.log 2>&1
echo "=== PARTIAL-DIFFUSION PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++; if($4+0<20)p20++} END{printf "scored=%d PASS(pae<10,plddt>80,rmsd<2)=%d pae<10any=%d pae<15=%d pae<20=%d\n",n,p+0,pl+0,p15+0,p20+0}' "$AF2IG/out.sc"
echo "--- top 15 ---"; (head -1 "$AF2IG/out.sc"; awk 'NR>1 && $4!="nan"' "$AF2IG/out.sc"|sort -k4 -n|head -15)
echo "POST_PARTIAL_DONE"
