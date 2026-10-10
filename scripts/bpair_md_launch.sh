#!/bin/bash
# Brutal two-state switch MD on the 8 selected bpair His-switch candidates (runs/bpair_mature/md_manifest.json),
# one per GPU, SAME gate params as the validated des_59/des_76 confirmations (replicas 3, prod-ps 3000,
# variants HIP,HIE,GLN). Then audit_brutal all -> SURVIVE/CUT. envs/md2 (OpenCL/CUDA openmm).
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
PY="$R/envs/md2/bin/python"
MANI="${1:-$R/runs/bpair_mature/md_manifest.json}"
OUT="${2:-$R/runs/bpair_mature/md_out}"; mkdir -p "$OUT"
REP="${3:-3}"
PROD="${4:-3000}"
mapfile -t LINES < <("$R/.venv/bin/python" -I -c "import json;[print(x['variant'],x['his_pos'],x['pdb'],x['gpu']) for x in json.load(open('$MANI'))]")
echo "[bpmd $(date -u +%FT%TZ)] launching ${#LINES[@]} brutal MDs"
pids=()
for ln in "${LINES[@]}"; do
  set -- $ln; V=$1; HP=$2; PDB=$3; G=$4
  od="$OUT/$V"; mkdir -p "$od"
  echo "  gpu$G  $V  his$HP"
  nohup "$PY" scripts/md_release_tnf.py --pdb "$PDB" --his-pos "$HP" --out "$od" \
      --binder-chain A --variants HIP,HIE,GLN --device "$G" --replicas "$REP" --prod-ps "$PROD" \
      > "$od/md.log" 2>&1 &
  pids+=($!)
done
echo "[bpmd] ${#pids[@]} MDs running; waiting..."
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[bpmd $(date -u +%FT%TZ)] all MDs done; auditing"
"$R/.venv/bin/python" -I scripts/audit_brutal.py "$OUT"/*/ 2>&1 || "$PY" scripts/audit_brutal.py "$OUT"/*/
echo "BPMD_DONE $(date -u +%FT%TZ)"
