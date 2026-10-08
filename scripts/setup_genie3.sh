#!/bin/bash
set -uo pipefail
G=/home/prachh/v/anthropic-adaptyv-2026-tnf/third_party/uplifting-biomolecular-modeling/genie3
W=/home/prachh/v/anthropic-adaptyv-2026-tnf/weights/genie3
cd "$G"
echo "[genie3-build] $(date -u +%FT%TZ) creating uv venv py3.10"
uv venv --seed --managed-python --python 3.10 venv_g3 2>&1 | tail -3
. venv_g3/bin/activate
echo "[genie3-build] installing pinned stack (torch 2.7.1+cu126 + stock genie3 VCS); this is the long pole"
pip install -q -r environment/requirements.lock --src "$G/stock_src" 2>&1 | tail -15
echo "[genie3-build] rc_lock=$? ; running run.sh install (kit editable + weights fetch)"
bash run.sh install --weights "$W" 2>&1 | tail -20
echo "[genie3-build] rc_install=$?"
export GENIE3_WEIGHTS="$W/pretrained/v1"
echo "[genie3-build] check:"
bash run.sh check --config h100 --mode fast 2>&1 | tail -20
echo "GENIE3_BUILD_DONE rc=$? $(date -u +%FT%TZ)"
