#!/bin/bash
# Genie3 Region-I-groove funnel: CA-trace backbones -> full-backbone (2-protomer) -> ProteinMPNN (chain A,
# 0-Cys, oversample) -> full-atom af2in (target B,C) -> af2ig. af2ig is a WEAK ranker for TNF (Overath:
# ipSAE 1.4x better), so this stage is a WIDE filter (record all pae; flag pae<10 strong, pae<16 candidate)
# — the ipSAE+pose stage (run_boltz_tnf + ipsae.py + pose-gate at the groove) does the real ranking next.
# Usage: funnel_genie3_regI.sh "<gpus 'a b c'>" [nseq=6] [tag=genie3_regI]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
IFS=' ' read -ra GA <<< "${1:-4 5 6 7}"; NSEQ="${2:-6}"; TAG="${3:-genie3_regI}"
GPUS_RANGE="$(echo "${GA[@]}" | tr ' ' ',')"
PULCHRA=/home/prachh/v/.claude-private/jobs/5045f83d/tmp/recon_venv/bin/pulchra
TGT="$REPO/data/targets/1TNF_regI_BC.pdb"
ESMC="$REPO/envs/uplifting-esmc/bin/python"
L="$REPO/third_party/LigandMPNN"; MPY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
FBB="$REPO/runs/${TAG}/fbb"; W="$REPO/runs/${TAG}/mpnn"; AIN="$REPO/runs/${TAG}/af2in"; AIG="$REPO/runs/${TAG}/af2ig"
rm -rf "$W" "$AIN" "$AIG"; mkdir -p "$FBB" "$W"   # keep FBB (reuse already-built complexes)
# 1) full-backbone (2-protomer), PARALLELIZED across NW workers (reuses existing outputs via skip-existing)
NW="${FBB_WORKERS:-16}"
echo "[funnel] fullbb x$NW workers $(date -u +%FT%TZ)"
fpids=()
for i in $(seq 0 $((NW-1))); do
  SHARD_I=$i SHARD_N=$NW "$ESMC" scripts/build_genie3_fullbb_bc.py "$PULCHRA" \
    "$REPO/runs/${TAG}/shard_g*/tnf_regI/pdbs/*.pdb" "$FBB" "$TGT" "B,C" > "$FBB/.w$i.log" 2>&1 &
  fpids+=($!)
done
for p in "${fpids[@]}"; do wait "$p" || true; done
echo "[funnel] fullbb done: $(ls "$FBB"/*.pdb 2>/dev/null | wc -l) complexes"
NFBB=$(ls "$FBB"/*.pdb 2>/dev/null | wc -l); echo "[funnel] fullbb complexes: $NFBB"
[ "$NFBB" -ge 1 ] || { echo "NO FBB COMPLEXES — stop"; exit 1; }
# 2) ProteinMPNN design chain A (0-Cys), sharded across GPUs, oversample NSEQ
mapfile -t ALL < <(ls "$FBB"/*.pdb); NG=${#GA[@]}; declare -A LS; gi=0
for p in "${ALL[@]}"; do g=${GA[$((gi%NG))]}; LS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in "${GA[@]}"; do
  [ -n "${LS[$g]:-}" ] || continue
  mkdir -p "$W/gpu$g"; printf '[%s]' "${LS[$g]%,}" > "$W/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$MPY" run.py --model_type protein_mpnn \
      --checkpoint_protein_mpnn "$CKPT" --pdb_path_multi "$W/gpu$g/multi.json" \
      --chains_to_design A --bias_AA "C:-100.0" --out_folder "$W/gpu$g" \
      --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$W/gpu$g.log" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[funnel] MPNN backbones: $(find "$W" -path '*/backbones/*.pdb'|wc -l)"
# 3) full-atom af2in (2-protomer target B,C) -> af2ig
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$W/gpu*/backbones/*.pdb" \
  --target "$TGT" --target-chains "B,C" --binder-chain A --out "$AIN" 2>&1 | grep -v VIRTUAL_ENV | tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus "$GPUS_RANGE" --mode fast \
  > "$REPO/runs/${TAG}_af2ig_run.log" 2>&1
echo "=== ${TAG} af2ig (WIDE filter; af2ig is weak for TNF) ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10&&$6+0>80&&$2+0<2)p10++; if($4+0<16)p16++} END{printf "scored=%d strong(pae<10,plddt>80,rmsd<2)=%d candidate(pae<16)=%d\n",n,p10+0,p16+0}' "$AIG/out.sc"
(head -1 "$AIG/out.sc"; awk 'NR>1 && $4!="nan"' "$AIG/out.sc" | sort -k4 -n | head -20)
echo "FUNNEL_GENIE3_REGI_DONE $(date -u +%FT%TZ)"
