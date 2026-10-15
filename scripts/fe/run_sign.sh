#!/bin/bash
# Plain-TI SIGN wave, GPU-sharded over a job TSV: name<TAB>his<TAB>leg<TAB>seed.
# Completes direct-FE release-direction (ddG sign) coverage across a submission portfolio: fills any
# member lacking bound+free legs and tops thin members to a target seed count. Config matches the cached
# campaign (EQ 0.5 / PRODB 2.5 / PRODF 2.0 / nwin 11 / sample 1.0) so new seeds POOL with the old ones.
# One process per GPU (lanes pull jobs stride-NG). Resumable: skips any leg whose .json exists.
# Engine: alchemy_ddg.py (terminal-His-safe). Analyze with analyze_sign.py. Findings: docs/SWITCH_FE_findings.md.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
export OPENMM_CUDA_COMPILER=/usr/local/cuda-12.3/bin/nvcc OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
PY=envs/md2/bin/python; FE=$R/scripts/fe
ENG=$FE/alchemy_ddg.py
OUT=$R/runs/fe/out; LOG=$OUT/logs_sign; mkdir -p "$OUT" "$LOG"
TSV=${TSV:?need TSV=path (cols: name<TAB>his<TAB>leg<TAB>seed)}
GPUS=(${GPUS:-0 1 2 3 4 5 6 7}); NG=${#GPUS[@]}
NWIN=${NWIN:-11}; EQ=${EQ:-0.5}; PRODB=${PRODB:-2.5}; PRODF=${PRODF:-2.0}; SAMP=${SAMP:-1.0}

mapfile -t JOBS < <(grep -vE '^#|^[[:space:]]*$' "$TSV")
echo "SIGN wave: ${#JOBS[@]} jobs over $NG GPUs (${GPUS[*]})  nwin=$NWIN eq=$EQ prodB=$PRODB prodF=$PRODF samp=$SAMP"

# Good-citizen guard: never clobber a foreign GPU job (SKIP_GUARD=1 overrides deliberately).
busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | wc -l)
if [ "${SKIP_GUARD:-0}" != 1 ] && [ "$busy" -ne 0 ]; then
  echo "ABORT: $busy foreign GPU compute process(es) present — not launching"; nvidia-smi; exit 1
fi

run_lane(){
  local li=$1 i=$1 gpu=${GPUS[$1]}
  while [ $i -lt ${#JOBS[@]} ]; do
    IFS=$'\t' read -r name his leg seed <<< "${JOBS[$i]}"
    local sys=$R/runs/fe/sys/${name}_${leg}
    local out=$OUT/${name}_${leg}_s${seed}
    local prod=$PRODF; [ "$leg" = bound ] && prod=$PRODB
    if [ -f "${out}.json" ]; then echo "skip $name $leg s$seed"; i=$((i+NG)); continue; fi
    if [ ! -f "${sys}_system.xml" ]; then echo "MISS sys $name $leg (build first)"; i=$((i+NG)); continue; fi
    CUDA_VISIBLE_DEVICES=$gpu $PY $ENG --sys "$sys" --switch-his "$his" --binder-chain A \
       --nwin $NWIN --equil-ns $EQ --prod-ns $prod --sample-ps $SAMP --device 0 --seed $seed \
       --out "$out" > "$LOG/${name}_${leg}_s${seed}.log" 2>&1
    echo "done gpu$gpu $name $leg s$seed rc=$?"
    i=$((i+NG))
  done
}
for li in $(seq 0 $((NG-1))); do run_lane $li & done
wait
echo "ALL SIGN JOBS DONE"
