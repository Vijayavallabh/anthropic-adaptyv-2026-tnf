#!/bin/bash
# Rebuild envs/pyrosetta with PyRosetta + propka (+ pKAI) for the pH-switch scoring (objective 1).
set -uo pipefail
cd "$(dirname "$0")/.."
echo "[switch-env] $(date -u +%T) creating venv"
uv venv --clear --seed --managed-python --python 3.11 envs/pyrosetta > logs/pyrosetta.log 2>&1
echo "[switch-env] installing pyrosetta-installer"
envs/pyrosetta/bin/python -m pip install pyrosetta-installer >> logs/pyrosetta.log 2>&1
echo "[switch-env] installing pyrosetta (large download)"
envs/pyrosetta/bin/python -c "import pyrosetta_installer as p; p.install_pyrosetta()" >> logs/pyrosetta.log 2>&1 \
  && echo "[switch-env] pyrosetta OK" || echo "[switch-env] pyrosetta FAILED (see logs/pyrosetta.log)"
echo "[switch-env] installing propka + numpy pin"
envs/pyrosetta/bin/python -m pip install 'propka==3.5.1' pandas >> logs/pyrosetta.log 2>&1 \
  && echo "[switch-env] propka OK" || echo "[switch-env] propka FAILED"
echo "[switch-env] verify imports:"
envs/pyrosetta/bin/python -c 'import pyrosetta,propka;print("  pyrosetta",pyrosetta.version() if hasattr(pyrosetta,"version") else "?","propka",propka.__version__)' 2>&1 | tail -2
echo "SWITCH_ENV_DONE $(date -u +%FT%TZ)"
