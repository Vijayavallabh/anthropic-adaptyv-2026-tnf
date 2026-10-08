#!/bin/bash
# Evaluate Genie3 binder backbones (PDB: binder=chain A, target=chain B) -> vanilla MPNN (design A, 0-Cys)
# -> full-atom target (1TNF_B) -> af2ig -> pass rate + pose. Genie3 is the CRO-validated TNF generator.
# Usage: eval_genie3.sh <backbone_glob> <tag> <gpus 0-7> <nseq> [target_pdb]
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
BBGLOB=$1; TAG=$2; GPUS=${3:-0-7}; NSEQ=${4:-8}; TGT=${5:-data/targets/1TNF_B.pdb}
PDB="$REPO/runs/${TAG}_pdb"; rm -rf "$PDB"; mkdir -p "$PDB"
# genie3 backbones are already PDB (binder A + target B); copy with a flat, unique name
n=0; for p in $BBGLOB; do [ -f "$p" ] || continue; b=$(echo "$p" | sed 's#.*genie3_tnf_##; s#/pdbs/#_#; s#/#_#g'); cp "$p" "$PDB/$b"; n=$((n+1)); done
echo "[$TAG] backbones: $n"
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
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$W/gpu*/backbones/*.pdb" --target "$TGT" --target-chains B --binder-chain A --out "$AIN" 2>&1|grep -v VIRTUAL_ENV|tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus "$GPUS" --mode fast > "$REPO/runs/${TAG}_af2ig_run.log" 2>&1
echo "=== ${TAG} Genie3->MPNN PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10&&$6+0>80&&$2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++; if($4+0<18)p18++} END{printf "scored=%d PASS(pae<10,plddt>80,rmsd<2)=%d pae<10=%d pae<15=%d pae<18=%d\n",n,p+0,pl+0,p15+0,p18+0}' "$AIG/out.sc"
(head -1 "$AIG/out.sc"; awk 'NR>1 && $4!="nan"' "$AIG/out.sc"|sort -k4 -n|head -12)
echo "EVAL_${TAG}_DONE"
