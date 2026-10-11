#!/bin/bash
# BOLD thrust (SESSION 11 cont.): Genie3 de-novo generation at the CRO-VALIDATED Region-I GROOVE
# (two-protomer B+C hotspots B69/B70/B92 + C109/C110 = V74/L75/I97 + W114/Y115). Genie3 is the only
# generator (with PXDesign) that produced CRO-measured TNF binders in the Anthropic campaign, AND gave
# ALL the human/mouse cross-reactive ones — but I only ever ran it at K90 (wrong site). This runs it at
# Region I. Shards g3batch.py across GPUs, one process/GPU, distinct seeds, per-GPU triton cache.
# Usage: gen_genie3_regI.sh "<gpu space-list>" [total_N=600] [seed_base=7100] [tag=genie3_regI]
# VERIFY GPUs are free first (never clobber a foreign job).
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
G3="$REPO/third_party/uplifting-biomolecular-modeling/genie3"
PY="$G3/venv_g3/bin/python"; BATCH="$G3/opt/genie3_opt/g3batch.py"
DS="$G3/stock_src/genie3/data/design/binder_design/tnf_regIbench"
IFS=' ' read -ra GPUS <<< "${1:-}"; NG=${#GPUS[@]}
[ "$NG" -ge 1 ] || { echo "usage: gen_genie3_regI.sh \"<gpus>\" [N] [seedbase] [tag]"; exit 2; }
TOTAL="${2:-600}"; SEEDBASE="${3:-7100}"; TAG="${4:-genie3_regI}"
PERSHARD=$(( (TOTAL + NG - 1) / NG ))
echo "[gen-regI] $TOTAL designs across ${NG} GPUs (${GPUS[*]}) = ~$PERSHARD/shard, seedbase $SEEDBASE  $(date -u +%FT%TZ)"
pids=()
i=0
for g in "${GPUS[@]}"; do
  SH="$REPO/runs/${TAG}/shard_g${g}"; mkdir -p "$SH" "$REPO/.jit/g3_triton_g${g}"
  REQ="$SH/request.yaml"
  python3 - "$SH" "$DS" "$((SEEDBASE + i))" "$PERSHARD" "$REQ" <<'PY'
import sys
sh,ds,seed,n,req=sys.argv[1:6]
open(req,"w").write(f"""experiment:
  name: tnf_regI_g{seed}
  seed: {seed}
paths:
  rootdir: {sh}
  dataset: {ds}
generation:
  dataset:
    source: target
    selections: tnf_regI
    n_sample: {n}
  sampler:
    sampler:
      direction_scale: 0.0
  base:
    checkpoint: {ds.split('/stock_src/')[0]}/../../../weights/genie3/pretrained/v1/checkpoints/step=600000.ckpt
    config: {ds.split('/stock_src/')[0]}/../../../weights/genie3/pretrained/v1/config.yaml
""")
print("wrote",req)
PY
  ( CUDA_VISIBLE_DEVICES="$g" TRITON_CACHE_DIR="$REPO/.jit/g3_triton_g${g}/cache" TRITON_HOME="$REPO/.jit/g3_triton_g${g}" \
      OMP_NUM_THREADS=16 "$PY" "$BATCH" --config "$REQ" --outdir "$SH" --batch-size 1 \
      --cuda-graphs --hoist --reuse-graphs 16 --lean-pair --wide-capture --alloc expandable \
      --pt-chunk design --tf32 --trimul fpf --compile --timings "$SH/timings.json" ) > "$SH/design.log" 2>&1 &
  pids+=($!); i=$((i+1))
done
echo "[gen-regI] launched ${#pids[@]} shards; waiting..."
for p in "${pids[@]}"; do wait "$p" || true; done
NB=$(find "$REPO/runs/${TAG}"/shard_g*/ -name "*.pdb" 2>/dev/null | wc -l)
echo "[gen-regI] DONE $(date -u +%FT%TZ)  backbones=$NB  (under runs/${TAG}/shard_g*/)"
echo "GEN_GENIE3_REGI_DONE"
