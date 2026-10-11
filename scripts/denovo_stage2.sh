#!/bin/bash
# De-novo stage-2 funnel: MPNN (chain A, 0-Cys, context B=human TNF) a flat set of de-novo backbones ->
# full-atom af2in -> af2ig binding filter. Generalized (bpair / any flat des_*.pdb set).
# Usage: denovo_stage2.sh "<input glob>" <base_out> "<gpus>" [nseq=8]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
GLOB="${1:-runs/bpair/v15gen/des_*.pdb}"; BASE="${2:-runs/bpair/funnel}"; GPUS="${3:-2 6 7}"; NSEQ="${4:-8}"
case "$BASE" in /*) ;; *) BASE="$REPO/$BASE";; esac   # MPNN cd's into $L, so $WORK must be absolute
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
STAGE="$BASE/staged"; WORK="$BASE/mpnn"; AF2IN="$BASE/af2in"; AF2IG="$BASE/af2ig"
rm -rf "$STAGE" "$WORK" "$AF2IN" "$AF2IG"; mkdir -p "$STAGE" "$WORK"
n=0; for p in $GLOB; do [ -e "$p" ]||continue; b=$(echo "$p" | sed -E 's#^.*/(out|v15gen)/##; s#/#_#g'); ln -sf "$REPO/$p" "$STAGE/$b"; n=$((n+1)); done
echo "[dn-s2] staged $n backbones -> $BASE; MPNN nseq=$NSEQ on [$GPUS] $(date -u +%FT%TZ)"
mapfile -t ALL < <(ls "$STAGE"/*.pdb 2>/dev/null)
[ ${#ALL[@]} -lt 1 ] && { echo "NO BACKBONES matched $GLOB"; exit 1; }
declare -A LISTS; gi=0; NG=$(echo $GPUS|wc -w)
for p in "${ALL[@]}"; do g=$(echo $GPUS|cut -d' ' -f$((gi%NG+1))); LISTS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in $GPUS; do mkdir -p "$WORK/gpu$g"; printf '[%s]' "${LISTS[$g]%,}" > "$WORK/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$WORK/gpu$g/multi.json" --chains_to_design "A" --bias_AA "C:-100.0" \
      --out_folder "$WORK/gpu$g" --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
NM=$(find "$WORK" -path '*/backbones/*.pdb'|wc -l); echo "[dn-s2] MPNN done: $NM seqs $(date -u +%FT%TZ)"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
   --target data/targets/1TNF_AB.pdb --target-chains B --binder-chain A --out "$AF2IN" 2>&1 | grep -v VIRTUAL_ENV | tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus "$(echo $GPUS|tr ' ' ',')" --mode fast > "$BASE/af2ig_run.log" 2>&1
echo "=== DE-NOVO BINDING PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<15)p15++} END{printf "scored=%d PASS(pae<10,plddt>80,rmsd<2)=%d pae<15=%d\n",n,p+0,p15+0}' "$AF2IG/out.sc" 2>/dev/null
echo "--- top 15 by pae_interaction ---"; (head -1 "$AF2IG/out.sc"; awk 'NR>1 && $4!="nan"' "$AF2IG/out.sc" 2>/dev/null|sort -k4 -n|head -15) 2>/dev/null
echo "DENOVO_STAGE2_DONE $(date -u +%FT%TZ)"
