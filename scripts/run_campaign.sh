#!/bin/bash
# Generic de-novo campaign: gen (env GPUS,LENS,N,HOTSPOTS) -> vanilla MPNN (NSEQ) -> single-protomer af2ig -> pass rate.
# Env: TAG (out subdir), GPUS="3 4 5 6 7", LENS="70 80 90 100 110", N=250, NSEQ=16, HOTS="[B86,B87,B88,B90,B91]"
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
TAG=${TAG:-bf}; GPUS=${GPUS:-"3 4 5 6 7"}; LENS=${LENS:-"70 80 90 100 110"}; N=${N:-250}; NSEQ=${NSEQ:-16}; HOTS=${HOTS:-"[B86,B87,B88,B90,B91]"}
GENOUT="$REPO/runs/${TAG}_rfd"; rm -rf "$GENOUT"; mkdir -p "$GENOUT"
echo "[$TAG] gen GPUS=$GPUS LENS=$LENS N=$N $(date -u +%T)"
# inline gen (per-GPU length, cache isolation)
KIT="$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
unset VIRTUAL_ENV; source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export WEIGHTS="$REPO/weights/rfdiffusion1" RFD_ROOT="$REPO/third_party/RFdiffusion" TMPDIR="${TMPDIR:-/tmp}"
read -r -a GA <<< "$GPUS"; read -r -a LA <<< "$LENS"; i=0; pids=()
for g in "${GA[@]}"; do
  L="${LA[$((i % ${#LA[@]}))]}"; pre="$GENOUT/g${g}_L${L}"; mkdir -p "$pre/sched" "$REPO/.jit/triton_g${g}"
  ( cd "$KIT" && CUDA_VISIBLE_DEVICES="$g" MODEL_OPT_JIT_ROOT="$REPO/.jit/L${L}" TRITON_CACHE_DIR="$REPO/.jit/triton_g${g}/cache" TRITON_HOME="$REPO/.jit/triton_g${g}" \
     bash run.sh design --config h100 --mode fast inference.input_pdb="$REPO/data/targets/1TNF_AB.pdb" \
     "contigmap.contigs=[A6-157/0 B6-157/0 ${L}-${L}]" "ppi.hotspot_res=$HOTS" \
     denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.num_designs="$N" \
     inference.schedule_directory_path="$pre/sched" inference.output_prefix="$pre/des" ) > "$GENOUT/gen_g${g}.log" 2>&1 &
  pids+=($!); i=$((i+1)); sleep 3
done
for p in "${pids[@]}"; do wait "$p"||true; done
NB=$(ls "$GENOUT"/g*/des_*.pdb 2>/dev/null|wc -l); echo "[$TAG] gen done: $NB backbones"
# MPNN
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
W="$REPO/runs/${TAG}_mpnn"; AIN="$REPO/runs/${TAG}_af2in"; AIG="$REPO/runs/${TAG}_af2ig"; rm -rf "$W" "$AIN" "$AIG"; mkdir -p "$W"
mapfile -t ALL < <(ls "$GENOUT"/g*/des_*.pdb); declare -A LS; gi=0
for p in "${ALL[@]}"; do g=${GA[$((gi % ${#GA[@]}))]}; LS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in "${GA[@]}"; do mkdir -p "$W/gpu$g"; printf '[%s]' "${LS[$g]%,}" > "$W/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
     --pdb_path_multi "$W/gpu$g/multi.json" --chains_to_design C --bias_AA "C:-100.0" \
     --out_folder "$W/gpu$g" --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$W/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
echo "[$TAG] MPNN: $(find "$W" -path '*/backbones/*.pdb'|wc -l)"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$W/gpu*/backbones/*.pdb" --target data/targets/1TNF_AB.pdb --target-chains B --binder-chain C --out "$AIN" 2>&1|grep -v VIRTUAL_ENV|tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus "$(echo $GPUS|tr ' ' ,)" --mode fast > "$REPO/runs/${TAG}_af2ig.log" 2>&1 || \
  uv run python scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus 3-7 --mode fast >> "$REPO/runs/${TAG}_af2ig.log" 2>&1
echo "=== ${TAG} PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10&&$6+0>80&&$2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++; if($4+0<18)p18++} END{printf "scored=%d PASS=%d pae<10=%d pae<15=%d pae<18=%d\n",n,p+0,pl+0,p15+0,p18+0}' "$AIG/out.sc"
(head -1 "$AIG/out.sc"; awk 'NR>1 && $4!="nan"' "$AIG/out.sc"|sort -k4 -n|head -15)
echo "CAMPAIGN_${TAG}_DONE"
