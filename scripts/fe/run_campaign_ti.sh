#!/bin/bash
# pH-switch alchemical ddG (TI) campaign, parallel across GPUs. Idempotent/resumable.
# Phase 1: build bound+box-matched free per design. Phase 2: TI legs in per-GPU lanes. Phase 3: rank.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
export OPENMM_CUDA_COMPILER=/usr/local/cuda-12.3/bin/nvcc
PY=envs/md2/bin/python; FE=$R/scripts/fe
GPUS=(${GPUS:-1 2 3 4 5 6 7}); NG=${#GPUS[@]}
NWIN=${NWIN:-11}; EQ=${EQ:-0.5}; PRODB=${PRODB:-2.5}; PRODF=${PRODF:-2.0}; SEEDS=${SEEDS:-"1 2"}
mkdir -p runs/fe/sys runs/fe/out

declare -A HIS=(
  [des_36133_g4_sa18]=45,53,75 [des_180_6_H9-14-38]=9,14,38
  [des_820_3_r3_sa4_r2_8]=27,33,47 [des_33012_g4_sa8]=11,23,51,53
  [des_495_6_mh_41_49]=41,49 [des_1559_3_mh_48_75]=48,75 )
# interleave tiers so an early partial read still spans the gradient
ORDER="des_36133_g4_sa18 des_495_6_mh_41_49 des_820_3_r3_sa4_r2_8 des_180_6_H9-14-38 des_33012_g4_sa8 des_1559_3_mh_48_75"
find_cif(){ ls runs/boltz*/out/by_seed/"$1"/s0/"$1"_model_0.cif 2>/dev/null | head -1; }

build_design(){ local n=$1 g=$2 his=${HIS[$1]} cif; cif=$(find_cif "$n")
  [ -z "$cif" ] && { echo "[SKIP $n no CIF]"; return; }
  [ -f runs/fe/sys/${n}_bound_system.xml ] || CUDA_VISIBLE_DEVICES=$g $PY "$FE/build_system.py" "$cif" runs/fe/sys/${n}_bound --switch-his "$his" --equil-ps 200 --device 0 >runs/fe/out/${n}_build_bound.log 2>&1
  [ -f runs/fe/sys/${n}_free_system.xml ]  || CUDA_VISIBLE_DEVICES=$g $PY "$FE/build_system.py" "$cif" runs/fe/sys/${n}_free  --switch-his "$his" --free --box-from runs/fe/sys/${n}_bound_box.txt --equil-ps 200 --device 0 >runs/fe/out/${n}_build_free.log 2>&1
  echo "[built $n $(date -u +%T)]" ; }

run_leg(){ local n=$1 lg=$2 sd=$3 g=$4 prod=$5 his=${HIS[$1]} out=runs/fe/out/${n}_${2}_s${3}
  [ -f ${out}.json ] && { echo "[have $n $lg s$sd]"; return; }
  CUDA_VISIBLE_DEVICES=$g $PY "$FE/alchemy_ddg.py" --sys runs/fe/sys/${n}_${lg} --switch-his "$his" \
    --binder-chain A --nwin $NWIN --equil-ns $EQ --prod-ns $prod --device 0 --seed $sd --out ${out} \
    >runs/fe/out/${n}_${lg}_s${sd}.log 2>&1
  grep -h "\[RESULT\]" runs/fe/out/${n}_${lg}_s${sd}.log 2>/dev/null || echo "[FAIL $n $lg s$sd]" ; }

echo "=== CAMPAIGN(TI) start $(date -u) GPUS=${GPUS[*]} nwin=$NWIN seeds='$SEEDS' ==="
i=0; for n in $ORDER; do build_design "$n" ${GPUS[$((i%NG))]} & i=$((i+1)); done; wait
echo "=== builds done $(date -u +%T) ==="

# task list -> per-GPU lanes (round-robin, each lane sequential, lanes concurrent, no barriers)
TASKS=(); for n in $ORDER; do for sd in $SEEDS; do TASKS+=("$n|bound|$sd|$PRODB" "$n|free|$sd|$PRODF"); done; done
for lane in $(seq 0 $((NG-1))); do
  ( g=${GPUS[$lane]}
    idx=$lane
    while [ $idx -lt ${#TASKS[@]} ]; do
      IFS='|' read -r n lg sd prod <<<"${TASKS[$idx]}"
      run_leg "$n" "$lg" "$sd" "$g" "$prod"
      idx=$((idx+NG))
    done ) &
done
wait
echo "=== legs done $(date -u +%T) ==="
$PY "$FE/combine_ddg.py" 2>&1 | grep -vE "^(Warning|JAX|\*|Either|pip|conda| [A-Z])"
echo "=== CAMPAIGN(TI) end $(date -u) ==="
