#!/bin/bash
# FE-screen the broader mature pool (runs/fe/screen_pool.tsv: name<TAB>his) for releasers.
# Stage-1: 2 seeds, bound 2.5ns / free 2.0ns. Idempotent. Parallel per-GPU lanes.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
export OPENMM_CUDA_COMPILER=/usr/local/cuda-12.3/bin/nvcc
PY=envs/md2/bin/python; FE=$R/scripts/fe
GPUS=(${GPUS:-1 2 3 4 5 6 7}); NG=${#GPUS[@]}
NWIN=${NWIN:-11}; EQ=${EQ:-0.5}; PRODB=${PRODB:-2.5}; PRODF=${PRODF:-2.0}; SEEDS=${SEEDS:-"1 2"}
TSV=${TSV:-runs/fe/screen_pool.tsv}
find_cif(){ ls runs/boltz*/out/by_seed/"$1"/s0/"$1"_model_0.cif 2>/dev/null | head -1; }

# load table
NAMES=(); HISS=()
while IFS=$'\t' read -r n h; do [ -z "$n" ] && continue; NAMES+=("$n"); HISS+=("$h"); done < "$TSV"
echo "=== SCREEN start $(date -u) ${#NAMES[@]} candidates GPUS=${GPUS[*]} ==="

build_one(){ local n=$1 his=$2 g=$3 cif; cif=$(find_cif "$n"); [ -z "$cif" ] && { echo "[SKIP $n noCIF]"; return; }
  [ -f runs/fe/sys/${n}_bound_system.xml ] || CUDA_VISIBLE_DEVICES=$g $PY "$FE/build_system.py" "$cif" runs/fe/sys/${n}_bound --switch-his "$his" --equil-ps 150 --device 0 >runs/fe/out/${n}_build_bound.log 2>&1
  [ -f runs/fe/sys/${n}_free_system.xml ]  || CUDA_VISIBLE_DEVICES=$g $PY "$FE/build_system.py" "$cif" runs/fe/sys/${n}_free  --switch-his "$his" --free --box-from runs/fe/sys/${n}_bound_box.txt --equil-ps 150 --device 0 >runs/fe/out/${n}_build_free.log 2>&1
  echo "[built $n]"; }
# build phase (NG at a time)
i=0; while [ $i -lt ${#NAMES[@]} ]; do
  build_one "${NAMES[$i]}" "${HISS[$i]}" ${GPUS[$((i%NG))]} &
  i=$((i+1)); (( i % NG == 0 )) && wait
done; wait
echo "=== builds done $(date -u +%T) ==="

# leg task list
TASKS=()
for idx in "${!NAMES[@]}"; do n=${NAMES[$idx]}; h=${HISS[$idx]}
  for sd in $SEEDS; do TASKS+=("$n|$h|bound|$sd|$PRODB" "$n|$h|free|$sd|$PRODF"); done
done
run_leg(){ IFS='|' read -r n his lg sd prod <<<"$1"; local g=$2 out=runs/fe/out/${n}_${lg}_s${sd}
  [ -f ${out}.json ] && { echo "[have $n $lg s$sd]"; return; }
  CUDA_VISIBLE_DEVICES=$g $PY "$FE/alchemy_ddg.py" --sys runs/fe/sys/${n}_${lg} --switch-his "$his" --binder-chain A \
    --nwin $NWIN --equil-ns $EQ --prod-ns $prod --device 0 --seed $sd --out ${out} >runs/fe/out/${n}_${lg}_s${sd}.log 2>&1
  grep -h "\[RESULT\]" runs/fe/out/${n}_${lg}_s${sd}.log 2>/dev/null || echo "[FAIL $n $lg s$sd]"; }
for lane in $(seq 0 $((NG-1))); do ( idx=$lane
  while [ $idx -lt ${#TASKS[@]} ]; do run_leg "${TASKS[$idx]}" ${GPUS[$lane]}; idx=$((idx+NG)); done ) & done
wait
echo "=== SCREEN legs done $(date -u) ==="
$PY "$FE/screen_analyze.py" 2>&1 | grep -vE "^(Warning|JAX|\*|Either|pip|conda| [A-Z])"
