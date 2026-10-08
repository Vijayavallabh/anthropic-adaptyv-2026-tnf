#!/bin/bash
# WAVE 25: His-ENRICHED soluble MPNN (per-residue His bias at the neutral-acceptor interface) -> full-atom
# target af2ig inputs -> AF2-ig fold filter. De novo (no fixed motif); His biased toward interface positions
# near conserved neutral acceptors (Ser464/Gln432/His433/Phe412) so designs can form a multi-His collective
# pH switch. Cys forbidden globally (0-Cys). Usage: multihis_mpnn_af2ig.sh [staged_dir] [gpus lo-hi] [nseq]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr
cd "$REPO"
STAGED=${1:-runs/wave25_multihis/staged}
GPUS=${2:-0-7}
NSEQ=${3:-4}
BASE=$(dirname "$STAGED")
WORK=$BASE/mpnn; AF2IN=$BASE/af2in; AF2IG=$BASE/af2ig
L="$REPO/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/solublempnn_v_48_020.pt"
MDPY="$REPO/envs/md/bin/python"
IFS=- read -r lo hi <<< "$GPUS"; GPULIST=$(seq "$lo" "$hi"); NG=$(echo $GPULIST | wc -w)

echo "=== [multihis] His-enriched soluble MPNN  $(date -u +%FT%TZ) ==="
mkdir -p "$WORK"
"$MDPY" - "$STAGED/pdbs.json" "$WORK" "$NG" <<'PY'
import json,sys,os
pdbs=json.load(open(sys.argv[1])); work=sys.argv[2]; ng=int(sys.argv[3])
for i in range(ng):
    os.makedirs(f"{work}/gpu{i}",exist_ok=True)
    json.dump(pdbs[i::ng], open(f"{work}/gpu{i}/shard.json","w"))
print(f"sharded {len(pdbs)} pdbs across {ng} gpus")
PY
pids=(); gi=0
for g in $GPULIST; do
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py \
      --model_type soluble_mpnn --checkpoint_soluble_mpnn "$CKPT" \
      --pdb_path_multi "$REPO/$WORK/gpu$gi/shard.json" \
      --bias_AA_per_residue_multi "$REPO/$STAGED/bias_AA_per_residue_multi.json" \
      --chains_to_design A --bias_AA "C:-100.0" \
      --out_folder "$REPO/$WORK/gpu$gi" \
      --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$gi.log" 2>&1 &
  pids+=($!); gi=$((gi+1))
done
for p in "${pids[@]}"; do wait "$p" || true; done
NBB=$(find "$WORK" -path '*/backbones/*.pdb' 2>/dev/null | wc -l)
echo "MPNN done: $NBB designed complexes ($(date -u +%FT%TZ))"
[ "$NBB" -gt 0 ] || { echo "FATAL: no MPNN backbones"; exit 1; }

echo "=== [multihis] full-atom-target af2ig inputs  $(date -u +%FT%TZ) ==="
"$MDPY" scripts/build_fullatom_af2in.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
  --target data/targets/6ARU_domainIII.pdb --target-chain A --out "$AF2IN"
NIN=$(find "$AF2IN" -name '*.pdb' 2>/dev/null | wc -l)
echo "af2in built: $NIN ($(date -u +%FT%TZ))"
[ "$NIN" -gt 0 ] || { echo "FATAL: no af2in inputs"; exit 1; }

echo "=== [multihis] AF2-ig fold filter  $(date -u +%FT%TZ) ==="
"$MDPY" scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus "$GPUS" --mode fast \
  || { echo "FATAL: run_af2ig.py exited non-zero"; exit 1; }
NPRED=$(find "$AF2IG" -name '*_af2pred.pdb' 2>/dev/null | wc -l)
[ "$NPRED" -gt 0 ] || { echo "FATAL: af2ig produced no *_af2pred.pdb"; exit 1; }
echo "MULTIHIS_MPNN_AF2IG_DONE ($NPRED af2pred) $(date -u +%FT%TZ)"
