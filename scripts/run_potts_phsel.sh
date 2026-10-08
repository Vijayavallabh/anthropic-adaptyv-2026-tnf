#!/usr/bin/env bash
# WAVE 16 - run the Proton-PottsMPNN His433 pH-switch corroboration (objective 1) on the 20 complexes.
# CPU-only (zero GPU). Requires the ProtonPottsMPNN venv + a built hbplus binary.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
PPM="$REPO/third_party/ProtonPottsMPNN"
export ADAPTYV_REPO="$REPO"
export HBPLUS_PATH="${HBPLUS_PATH:-$REPO/third_party/HBPLUS/hbplus/hbplus}"
export CUDA_VISIBLE_DEVICES=""          # CPU only; never touches a GPU
[ -x "$PPM/.venv/bin/python" ] || { echo "[FATAL] ProtonPottsMPNN venv missing; run $PPM/install.sh"; exit 1; }
[ -x "$HBPLUS_PATH" ]          || { echo "[FATAL] hbplus not found at $HBPLUS_PATH"; exit 1; }
exec "$PPM/.venv/bin/python" -P "$REPO/scripts/obj1_potts_phsel.py" "$@"
