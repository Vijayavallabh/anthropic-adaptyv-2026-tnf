#!/bin/bash
# Boltz-2 single-seq cofold (binder:target DIII, msa:empty) for a submission-style CSV (name,sequence).
# Independent 2nd binding predictor -> iptm. Usage: run_boltz.sh <csv> <tag> <gpu>
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
CSV=$1; TAG=${2:-fin}; GPU=${3:-1}
MD="$REPO/envs/md/bin/python"
TGT=$(tr -d '\n ' < data/targets/human_domainIII.seq)
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
echo "[boltz] gpu=$GPU input=$IN out=$OUT $(date -u +%FT%TZ)"
bash run.sh pred --config h100 --mode fast --input "$IN" --out_dir "$OUT"
echo "BOLTZ_${TAG}_DONE $(date -u +%FT%TZ) | confidence jsons: $(find "$OUT" -name 'confidence_*.json' 2>/dev/null | wc -l)"
