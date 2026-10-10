#!/bin/bash
# beta-pairing funnel: wait for RFd scale-up, then vanilla ProteinMPNN (clean binders, NO His bias) ->
# af2ig (fast) on gpu 4-7 with the LCB3 de-novo positive control in-batch -> gate.
# Switch-His install happens LATER (His-scan) on the af2ig+cross-react winners.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
log(){ echo "[bpair-funnel $(date -u +%H:%M:%S)] $*"; }
# 1. wait for RFd scale-up to finish (no run_bpair_rfd processes)
log "waiting for RFd scale-up..."
while pgrep -f "run_bpair_rfd.sh rfd " >/dev/null 2>&1; do sleep 60; done
NB=$(ls runs/bpair/rfd/des_*.pdb 2>/dev/null | wc -l); log "RFd done: $NB backbones"
# 2. vanilla ProteinMPNN, design chain A (binder), fix target B, omit C, 4 seqs/backbone, gpu 4-7
L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
W="$R/runs/bpair/mpnn"; AIN="$R/runs/bpair/af2in"; AIG="$R/runs/bpair/af2ig"; rm -rf "$W" "$AIN" "$AIG"; mkdir -p "$W"
mapfile -t ALL < <(ls runs/bpair/rfd/des_*.pdb)
GL="4 5 6 7"; NG=4; declare -A LS; gi=0
for p in "${ALL[@]}"; do g=$(echo $GL|cut -d' ' -f$((gi%NG+1))); LS[$g]+="\"$(readlink -f "$p")\","; gi=$((gi+1)); done
pids=()
for g in $GL; do mkdir -p "$W/gpu$g"; printf '[%s]' "${LS[$g]%,}" > "$W/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
     --pdb_path_multi "$W/gpu$g/multi.json" --chains_to_design A --omit_AA C \
     --out_folder "$W/gpu$g" --batch_size 4 --number_of_batches 1 --temperature 0.15 ) > "$W/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
NS=$(find "$W" -path '*/backbones/*.pdb'|wc -l); log "MPNN backbones: $NS"
# 3. build af2ig inputs (real full-atom target) + drop in the LCB3 positive control
"$R/.venv/bin/python" -I scripts/build_fullatom_af2in_tnf.py --mpnn "$W/gpu*/backbones/*.pdb" \
  --target data/targets/1TNF_B.pdb --target-chains "B" --binder-chain A --out "$AIN" 2>&1 | tail -2
cp runs/_pc_denovo/in/7JZM_lcb3.pdb "$AIN"/ 2>/dev/null && log "added LCB3 control" || log "WARN no control"
NA=$(ls "$AIN"/*.pdb 2>/dev/null|wc -l); log "af2in complexes: $NA"
# 4. af2ig fast on gpu 4-7
"$R/.venv/bin/python" -I scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus "4-7" --mode fast > runs/bpair/af2ig_run.log 2>&1
# 5. gate + control check
log "=== af2ig gate (pae_int<10 & plddt_binder>85 & rmsd<2) ==="
cat "$AIG"/gpu*/out.sc "$AIG"/out.sc 2>/dev/null | awk '$1=="SCORE:"&&$4!="nan"{n++; if($4+0<10&&$6+0>85&&$2+0<2)p++} END{printf "scored=%d PASS=%d\n",n,p+0}'
echo "--- control (expect pae~9.8) ---"
cat "$AIG"/gpu*/out.sc "$AIG"/out.sc 2>/dev/null | awk '$1=="SCORE:"&&$NF~/lcb3/{print "  LCB3 control: pae_int="$4" plddt_binder="$6" rmsd="$2}'
echo "--- top 30 passers (pae plddt rmsd name) ---"
cat "$AIG"/gpu*/out.sc "$AIG"/out.sc 2>/dev/null | awk '$1=="SCORE:"&&$4!="nan"&&$4+0<10&&$6+0>85&&$2+0<2{print $4,$6,$2,$NF}' | sort -n | awk '!s[$NF]++' | head -30 | tee runs/bpair/af2ig_passers.txt
echo "BPAIR_FUNNEL_DONE $(date -u +%FT%TZ)"
