#!/bin/bash
# Boltz-2 cofold of binders (CSV name,sequence) vs an ARBITRARY target seq file -> iptm.
# Binder-specificity / decoy-TARGET control (generalizes run_boltz.sh, whose target is fixed
# to DIII). Fold the same binders against off-target domains (EGFR DI/II/IV) or decoys to test
# that binding is DIII-specific, not non-specific (e.g. electrostatic stickiness of acidic binders).
# Usage: run_boltz_vs.sh <csv> <target_seq_file> <tag> <gpu>
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
CSV=$1; TGTF=$2; TAG=${3:-spec}; GPU=${4:-1}
MD="$REPO/envs/md/bin/python"
TGT=$(tr -d '\n ' < "$TGTF")
IN="$REPO/runs/boltz_${TAG}/in"; OUT="$REPO/runs/boltz_${TAG}/out"; mkdir -p "$IN" "$OUT"
"$MD" - "$CSV" "$IN" "$TGT" <<'PY'
import sys, csv, os
csvf, ind, tgt = sys.argv[1], sys.argv[2], sys.argv[3]
n=0
for r in csv.DictReader(open(csvf)):
    nm = r.get("name") or r.get("design"); seq = r["sequence"].split(":")[0]
    open(os.path.join(ind, nm + ".yaml"), "w").write(
        "version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: %s\n      msa: empty\n"
        "  - protein:\n      id: B\n      sequence: %s\n      msa: empty\n" % (seq, tgt))
    n += 1
print("wrote", n, "YAMLs ->", ind)
PY
export CUDA_VISIBLE_DEVICES=$GPU
unset VIRTUAL_ENV 2>/dev/null || true
source "$REPO/envs/uplifting-boltz2/bin/activate"
set -a; source "$REPO/envs/uplifting.env"; set +a
cd "$REPO/third_party/uplifting-biomolecular-modeling/boltz2"
echo "[boltz] gpu=$GPU input=$IN out=$OUT target=$TGTF $(date -u +%FT%TZ)"
bash run.sh pred --config h100 --mode fast --input "$IN" --out_dir "$OUT"
echo "BOLTZ_${TAG}_DONE $(date -u +%FT%TZ) | confidence jsons: $(find "$OUT" -name 'confidence_*.json' 2>/dev/null | wc -l)"
