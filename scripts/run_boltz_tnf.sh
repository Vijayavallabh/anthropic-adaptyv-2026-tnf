#!/bin/bash
# Boltz-2 orthogonal binding confirmation: cofold binders (CSV name,sequence) vs a target seq file -> iptm.
# Independent oracle vs af2-initial-guess. Use human TNF (data/targets/1TNF_B.pdb seq) for binding,
# mouse TNF (data/targets/mouse_TNF_P06804.fasta) for a cross-reactivity iptm delta.
# Usage: run_boltz_tnf.sh <csv name,sequence> <target_seq_file> <tag> <gpu>
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
CSV=$1; TGTF=$2; TAG=${3:-tnf}; GPU=${4:-0}
MD="$REPO/envs/md2/bin/python"
# target seq: if a fasta, strip header lines; else take all non-space chars
TGT=$(grep -v '^>' "$TGTF" 2>/dev/null | tr -d '\n ' || tr -d '\n ' < "$TGTF")
IN="$REPO/runs/boltz_${TAG}/in"; OUT="$REPO/runs/boltz_${TAG}/out"; mkdir -p "$IN" "$OUT"
"$MD" - "$CSV" "$IN" "$TGT" <<'PY'
import sys, csv, os
csvf, ind, tgt = sys.argv[1], sys.argv[2], sys.argv[3]
n=0
for r in csv.DictReader(open(csvf)):
    nm = r.get("name") or r.get("design") or r.get("design_id"); seq = r["sequence"].split(":")[0]
    open(os.path.join(ind, nm + ".yaml"), "w").write(
        "version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: %s\n      msa: empty\n"
        "  - protein:\n      id: B\n      sequence: %s\n      msa: empty\n" % (seq, tgt))
    n += 1
print("wrote", n, "YAMLs ->", ind)
PY
export CUDA_VISIBLE_DEVICES=$GPU
unset VIRTUAL_ENV 2>/dev/null || true
source "$REPO/envs/uplifting-boltz2/bin/activate"
set -a; source "$REPO/envs/uplifting.env" 2>/dev/null || true; set +a
cd "$REPO/third_party/uplifting-biomolecular-modeling/boltz2"
echo "[boltz] gpu=$GPU input=$IN out=$OUT target=$TGTF $(date -u +%FT%TZ)"
bash run.sh pred --config h100 --mode fast --input "$IN" --out_dir "$OUT"
echo "BOLTZ_${TAG}_DONE $(date -u +%FT%TZ) | confidence jsons: $(find "$OUT" -name 'confidence_*.json' 2>/dev/null | wc -l)"
