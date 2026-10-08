#!/bin/bash
# Evaluate RFd3 Region-I (2-protomer groove) designs: cif.gz (binder=A, target=B,C) -> MPNN(design A) ->
# full-atom 2-protomer target (Kabsch-superpose clean 1TNF_AB) -> af2ig -> pose-gated pass rate.
# Usage: eval_rfd3_regionI.sh <cif_glob> <tag> <gpus 0-7> <nseq>
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
CIFGLOB=$1; TAG=$2; GPUS=${3:-2-7}; NSEQ=${4:-8}
PDB="$REPO/runs/${TAG}_pdb"; rm -rf "$PDB"; mkdir -p "$PDB"
uv run python scripts/rfd3_to_pdb.py "$CIFGLOB" "$PDB" 2>&1 | grep -v VIRTUAL_ENV | tail -1
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
W="$REPO/runs/${TAG}_mpnn"; AIN="$REPO/runs/${TAG}_af2in"; AIG="$REPO/runs/${TAG}_af2ig"; rm -rf "$W" "$AIN" "$AIG"; mkdir -p "$W"
IFS=- read lo hi <<< "$GPUS"; GL=$(seq $lo $hi); NG=$(echo $GL|wc -w)
mapfile -t ALL < <(ls "$PDB"/*.pdb); declare -A LS; gi=0
for p in "${ALL[@]}"; do g=$(echo $GL|cut -d' ' -f$((gi%NG+1))); LS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in $GL; do mkdir -p "$W/gpu$g"; printf '[%s]' "${LS[$g]%,}" > "$W/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
     --pdb_path_multi "$W/gpu$g/multi.json" --chains_to_design A --bias_AA "C:-100.0" \
     --out_folder "$W/gpu$g" --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$W/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
echo "[$TAG] MPNN: $(find "$W" -path '*/backbones/*.pdb'|wc -l)"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$W/gpu*/backbones/*.pdb" \
    --target data/targets/1TNF_AB.pdb --target-chains "B,C" --binder-chain A --out "$AIN" 2>&1|grep -v VIRTUAL_ENV|tail -2
uv run python scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus "$GPUS" --mode fast > "$REPO/runs/${TAG}_af2ig_run.log" 2>&1
echo "=== ${TAG} Region-I RFd3->MPNN PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10&&$6+0>80&&$2+0<2)p++; if($4+0<14)p14++; if($4+0<18)p18++} END{printf "scored=%d PASS(pae<10)=%d pae<14=%d pae<18=%d\n",n,p+0,p14+0,p18+0}' "$AIG/out.sc"
(head -1 "$AIG/out.sc"; awk 'NR>1 && $4!="nan"' "$AIG/out.sc"|sort -k4 -n|head -12)
echo "EVAL_${TAG}_DONE"
