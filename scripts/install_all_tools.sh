#!/bin/bash
# Best-effort installer for the Anthropic uplifting kit envs + weights, plus PyRosetta.
# Each kit -> its own isolated uv venv (scripts/setup_uplifting.sh <kit> --weights).
# Runs sequentially, continues past failures, logs per-kit to logs/, summary in
# logs/install-summary.log. Designed to run in the background (can take hours; several GB/kit).
cd "$(dirname "$0")/.."
mkdir -p logs weights .jit envs
export MODEL_OPT_JIT_ROOT="$PWD/.jit"
SUMMARY=logs/install-summary.log
: > "$SUMMARY"
log(){ echo "[$(date -u +%FT%TZ)] $*" | tee -a "$SUMMARY"; }

# Kits to install: pass a list as args, else the full default set.
# Order: clean pip stacks first, then JAX/finicky, then LM.
KITS="${*:-boltz2 chai1 proteinmpnn caliby af2ig boltzgen rfdiffusion1 esmc}"
log "START uplifting kit installs: $KITS"
for k in $KITS; do
  log "$k: installing (venv + weights)..."
  if timeout 5400 bash scripts/setup_uplifting.sh "$k" --weights > "logs/uplifting-$k.log" 2>&1; then
    log "$k: OK"
  else
    rc=$?
    log "$k: FAILED (rc=$rc) -> logs/uplifting-$k.log"
  fi
done

# PyRosetta (academic / non-commercial, best-effort) for the stage-4 two-pH ddG scoring.
log "pyrosetta: installing (academic/non-commercial, best-effort)..."
if uv venv --clear --seed --managed-python --python 3.11 envs/pyrosetta > logs/pyrosetta.log 2>&1 \
   && envs/pyrosetta/bin/python -m pip install pyrosetta-installer >> logs/pyrosetta.log 2>&1 \
   && envs/pyrosetta/bin/python -c "import pyrosetta_installer as p; p.install_pyrosetta()" >> logs/pyrosetta.log 2>&1; then
  log "pyrosetta: OK"
else
  log "pyrosetta: FAILED -> logs/pyrosetta.log (may need manual license acceptance)"
fi

log "DONE. Summary above; per-kit logs in logs/. Verify each with:"
log "  source envs/uplifting-<kit>/bin/activate && (cd third_party/uplifting-biomolecular-modeling/<kit> && bash run.sh check)"
