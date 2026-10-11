#!/bin/bash
# v15 stage-2: MPNN (chain A, context B=human TNF, 0-Cys) the 480 partial backbones, build full-atom
# af2in, af2ig fold vs 2-protomer human TNF -> binding filter (pae_interaction). GPUs 1,3,4,5 (free).
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
GPUS="1 3 4 5"; NSEQ="${1:-4}"
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
BASE="$REPO/runs/partial_v15"; STAGE="$BASE/staged"; WORK="$BASE/mpnn"; AF2IN="$BASE/af2in"; AF2IG="$BASE/af2ig"
rm -rf "$STAGE" "$WORK" "$AF2IN" "$AF2IG"; mkdir -p "$STAGE" "$WORK"
n=0; for p in "$BASE"/out/*/pd_*.pdb; do [ -e "$p" ]||continue; b=$(echo "$p"|sed "s#$BASE/out/##;s#/#_#g"); ln -sf "$p" "$STAGE/$b"; n=$((n+1)); done
echo "[v15s2] staged $n backbones; MPNN nseq=$NSEQ on GPUs [$GPUS] $(date -u +%FT%TZ)"
mapfile -t ALL < <(ls "$STAGE"/*.pdb)
declare -A LISTS; gi=0; NG=$(echo $GPUS|wc -w)
for p in "${ALL[@]}"; do g=$(echo $GPUS|cut -d' ' -f$((gi%NG+1))); LISTS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in $GPUS; do mkdir -p "$WORK/gpu$g"; printf '[%s]' "${LISTS[$g]%,}" > "$WORK/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$WORK/gpu$g/multi.json" --chains_to_design "A" --bias_AA "C:-100.0" \
      --out_folder "$WORK/gpu$g" --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
NM=$(find "$WORK" -path '*/backbones/*.pdb'|wc -l); echo "[v15s2] MPNN done: $NM seqs $(date -u +%FT%TZ)"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
   --target data/targets/1TNF_AB.pdb --target-chains B --binder-chain A --out "$AF2IN" 2>&1 | grep -v VIRTUAL_ENV | tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus "1,3,4,5" --mode fast > "$BASE/af2ig_run.log" 2>&1
echo "=== v15 PARTIAL-DIFFUSION BINDING PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++} END{printf "scored=%d PASS(pae<10,plddt>80,rmsd<2)=%d pae<10any=%d pae<15=%d\n",n,p+0,pl+0,p15+0}' "$AF2IG/out.sc" 2>/dev/null
echo "--- top 20 by pae_interaction ---"; (head -1 "$AF2IG/out.sc"; awk 'NR>1 && $4!="nan"' "$AF2IG/out.sc" 2>/dev/null|sort -k4 -n|head -20) 2>/dev/null
echo "V15_STAGE2_DONE $(date -u +%FT%TZ)"
