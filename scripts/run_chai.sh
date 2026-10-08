#!/bin/bash
# Chai-1 single-seq cofold (binder chain A : target DIII chain B, no MSA) for a submission CSV.
# INDEPENDENT 3rd structure predictor (AF2-ig + Boltz-2 + Chai-1 binding consensus on the exact
# submitted sequence). Batches designs across GPUs via per-GPU items.json (model loads once/GPU).
# Usage: run_chai.sh <csv> <tag> [target_seq_file] [gpus]
#   target_seq_file default = data/targets/human_domainIII.seq ; gpus default = 0-2 (user ≤3-GPU cap)
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
CSV=${1:?usage: run_chai.sh <csv> <tag> [target_seq_file] [gpus]}; TAG=${2:-chai}
TGTF=${3:-data/targets/human_domainIII.seq}; GPUS=${4:-0-2}   # default 0-2: user ≤3-GPU cap (2026-10-02); pass arg 4 to override
MD="$REPO/envs/md/bin/python"
[ -s "$CSV" ]  || { echo "run_chai.sh: submission CSV missing/empty: $CSV" >&2; exit 2; }
[ -s "$TGTF" ] || { echo "run_chai.sh: target seq file missing/empty: $TGTF" >&2; exit 2; }
TGT=$(tr -d '\n ' < "$TGTF")
[ -n "$TGT" ]  || { echo "run_chai.sh: empty target sequence from $TGTF" >&2; exit 2; }
IN="$REPO/runs/chai_${TAG}/in"; OUT="$REPO/runs/chai_${TAG}/out"; mkdir -p "$IN" "$OUT"
"$MD" - "$CSV" "$IN" "$TGT" "$GPUS" <<'PY'
import sys, csv, os, json
csvf, ind, tgt, gpus = sys.argv[1:5]
G = list(range(int(gpus.split("-")[0]), int(gpus.split("-")[1])+1)) if "-" in gpus else [int(x) for x in gpus.split(",")]
rows = list(csv.DictReader(open(csvf)))
buckets = {g: [] for g in G}
n = 0
for i, r in enumerate(rows):
    nm = (r.get("name") or r.get("design") or "").strip()
    seq = (r.get("sequence") or "").split(":")[0].strip()
    if not nm:  sys.exit(f"run_chai fasta-gen: row {i} has no name/design value")
    if not seq: sys.exit(f"run_chai fasta-gen: {nm} has empty sequence")
    fa = os.path.join(ind, nm + ".fasta")
    open(fa, "w").write(f">protein|name=A\n{seq}\n>protein|name=B\n{tgt}\n")
    _seeds = [int(x) for x in os.environ.get("CHAI_SEEDS", "0").split(",") if x.strip() != ""]
    buckets[G[i % len(G)]].append({"fasta": fa, "seeds": _seeds or [0]}); n += 1
if not n: sys.exit("run_chai fasta-gen: no rows in CSV")
for g, items in buckets.items():
    json.dump({"items": items}, open(os.path.join(ind, f"items_gpu{g}.json"), "w"))
print("wrote", n, "fastas; per-gpu:", {g: len(v) for g, v in buckets.items()})
PY
rc=$?; [ "$rc" -eq 0 ] || { echo "run_chai.sh: FASTA/items generation failed (rc=$rc)" >&2; exit "$rc"; }
unset VIRTUAL_ENV 2>/dev/null || true
source "$REPO/envs/uplifting-chai1/bin/activate"
set -a; source "$REPO/envs/uplifting.env"; set +a
cd "$REPO/third_party/uplifting-biomolecular-modeling/chai1"
pids=()
shopt -s nullglob
ITEMSLIST=("$IN"/items_gpu*.json)
[ "${#ITEMSLIST[@]}" -gt 0 ] || { echo "run_chai.sh: no items_gpu*.json produced -> nothing to run" >&2; exit 2; }
for ITEMS in "${ITEMSLIST[@]}"; do
  G=${ITEMS##*items_gpu}; G=${G%.json}
  # skip empty buckets (fewer designs than GPUs): launching chai on 0 items can error -> false rc=1
  nit=$("$MD" -c "import json,sys; print(len(json.load(open(sys.argv[1])).get('items',[])))" "$ITEMS" 2>/dev/null || echo 0)
  [ "$nit" -gt 0 ] || { echo "[chai] gpu=$G: 0 items, skipping"; continue; }
  echo "[chai] gpu=$G items=$ITEMS $(date -u +%FT%TZ)"
  CUDA_VISIBLE_DEVICES=$G bash run.sh pred --config h100 --mode fast --input "$ITEMS" \
    --out_dir "$OUT" > "$REPO/runs/chai_${TAG}/chai_gpu${G}.log" 2>&1 &
  pids+=($!)
done
echo "[chai] launched ${#pids[@]} gpu procs tag=$TAG"
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "CHAI_${TAG}_DONE rc=$rc $(date -u +%FT%TZ) | pred dirs: $(find "$OUT" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | wc -l)"
