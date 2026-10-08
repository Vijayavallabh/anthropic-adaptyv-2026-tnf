#!/bin/bash
# New-epitope (K90/Y87/V91) RFd de-novo -> vanilla MPNN (chain C) -> SINGLE-protomer af2ig -> pass rate.
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
echo "[e2] waiting for RFd gen (gen_epitope2.log) $(date -u +%T)"
for i in $(seq 1 180); do grep -q 'sweep done' runs/gen_epitope2.log 2>/dev/null && break; sleep 20; done
NB=$(ls runs/rfd_gen/g*/des_*.pdb 2>/dev/null|wc -l); echo "[e2] $NB RFd backbones"
uv run python scripts/rawcontacts2.py "runs/rfd_gen/g*/des_*.pdb" 2>&1 | grep -v VIRTUAL_ENV || true
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
STAGE="$REPO/runs/e2/staged"; WORK="$REPO/runs/e2/mpnn"; AF2IN="$REPO/runs/e2/af2in"; AF2IG="$REPO/runs/e2/af2ig"
rm -rf "$REPO/runs/e2"; mkdir -p "$STAGE" "$WORK"
n=0; for p in runs/rfd_gen/g*/des_*.pdb; do [ -e "$p" ]||continue; b=$(echo "$p"|sed 's#runs/rfd_gen/##;s#/#_#g'); ln -sf "$REPO/$p" "$STAGE/$b"; n=$((n+1)); done
GL="0 1 2"; NG=3; mapfile -t ALL < <(ls "$STAGE"/*.pdb)
declare -A LISTS; gi=0
for p in "${ALL[@]}"; do g=$(echo $GL|cut -d' ' -f$((gi%NG+1))); LISTS[$g]+="\"$p\","; gi=$((gi+1)); done
pids=()
for g in $GL; do mkdir -p "$WORK/gpu$g"; printf '[%s]' "${LISTS[$g]%,}" > "$WORK/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$WORK/gpu$g/multi.json" --chains_to_design "C" --bias_AA "C:-100.0" \
      --out_folder "$WORK/gpu$g" --batch_size 8 --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
echo "[e2] MPNN done: $(find "$WORK" -path '*/backbones/*.pdb'|wc -l)"
uv run python scripts/build_fullatom_af2in_tnf.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
   --target data/targets/1TNF_AB.pdb --target-chains B --binder-chain C --out "$AF2IN" 2>&1 | grep -v VIRTUAL_ENV | tail -1
uv run python scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus 0-2 --mode fast > runs/e2/af2ig_run.log 2>&1
echo "=== NEW-EPITOPE (K90/Y87) PASS RATE ==="
awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<10)pl++; if($4+0<15)p15++; if($4+0<20)p20++} END{printf "scored=%d PASS=%d pae<10any=%d pae<15=%d pae<20=%d\n",n,p+0,pl+0,p15+0,p20+0}' "$AF2IG/out.sc"
echo "--- top 12 ---"; (head -1 "$AF2IG/out.sc"; awk 'NR>1 && $4!="nan"' "$AF2IG/out.sc"|sort -k4 -n|head -12)
echo "PIPE_E2_DONE"
