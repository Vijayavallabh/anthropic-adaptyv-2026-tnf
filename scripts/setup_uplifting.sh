#!/bin/bash
# Set up one Anthropic uplifting-biomolecular-modeling kit in an ISOLATED uv venv (Route C,
# per the kit's STOCK.md):
#   uv venv --clear --seed --managed-python
#   -> install the kit's pinned stack   (environment/requirements*.lock, --no-deps; several GB)
#   -> install the bundled stock wheel  (stock/*.whl, --no-deps) if present
#   -> bash run.sh install [--weights DIR]   (fills the weights cache; prints its cache env var)
#   -> export + record that cache var (envs/uplifting.env) -> bash run.sh check
#
#   scripts/setup_uplifting.sh <kit> [--python X.Y] [--weights]
#
# Not `set -e`: each step logs and continues so one failure doesn't hide the rest.
set -uo pipefail
unset VIRTUAL_ENV 2>/dev/null || true
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"

KIT="${1:-}"
[ -z "$KIT" ] && { echo "usage: scripts/setup_uplifting.sh <kit> [--python X.Y] [--weights]"; exit 2; }
shift
PY=""; DO_WEIGHTS=0
while [ $# -gt 0 ]; do
  case "$1" in
    --python) PY="$2"; shift 2 ;;
    --weights) DO_WEIGHTS=1; shift ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done
# Python pin per kit STOCK.md: most kits 3.11; caliby and esmc need 3.12.
if [ -z "$PY" ]; then case "$KIT" in caliby|esmc) PY=3.12 ;; *) PY=3.11 ;; esac; fi

UP="third_party/uplifting-biomolecular-modeling"
KITDIR="$UP/$KIT"
[ -d "$KITDIR" ] || { echo "kit '$KIT' not found under $UP -- run ./init.sh --tools to clone"; exit 1; }
VENV="$HERE/envs/uplifting-$KIT"
export MODEL_OPT_JIT_ROOT="$HERE/.jit"
ENVFILE="$HERE/envs/uplifting.env"
mkdir -p weights "$HERE/.jit" envs logs
touch "$ENVFILE"

echo "== [$KIT] uv venv (python $PY) =="
uv venv --clear --seed --managed-python --python "$PY" "$VENV" || { echo "[$KIT] venv FAILED"; exit 1; }
# shellcheck disable=SC1091
source "$VENV/bin/activate"

LOCK="$(ls "$KITDIR"/environment/requirements*.lock 2>/dev/null | head -1)"
if [ -n "$LOCK" ]; then
  echo "== [$KIT] installing pinned stack ($(basename "$LOCK"), --no-deps; may be several GB) =="
  uv pip install --no-deps -r "$LOCK" || echo "[$KIT] WARN: stack install reported errors (continuing)"
fi

shopt -s nullglob
whl=("$KITDIR"/stock/*.whl)
if [ ${#whl[@]} -gt 0 ]; then
  echo "== [$KIT] installing bundled stock wheel: $(basename "${whl[0]}") =="
  uv pip install --no-deps --reinstall "${whl[0]}" || echo "[$KIT] WARN: stock wheel install failed"
fi
shopt -u nullglob

INSTALL_LOG="$HERE/logs/.install-$KIT.out"
echo "== [$KIT] run.sh install $( [ $DO_WEIGHTS -eq 1 ] && echo '--weights' ) =="
if [ "$DO_WEIGHTS" -eq 1 ]; then
  ( cd "$KITDIR" && bash run.sh install --weights "$HERE/weights/$KIT" ) 2>&1 | tee "$INSTALL_LOG"
else
  ( cd "$KITDIR" && bash run.sh install ) 2>&1 | tee "$INSTALL_LOG"
fi

# Discover + apply the kit's weights-cache env var (run.sh prints "export VAR=<dir>").
HINT="$(grep -hoE "export [A-Z][A-Z0-9_]*=[^ )]*weights/$KIT[^ )]*" "$INSTALL_LOG" 2>/dev/null | tail -1)"
if [ -n "$HINT" ]; then
  eval "$HINT"
  grep -qF "$HINT" "$ENVFILE" 2>/dev/null || echo "$HINT" >> "$ENVFILE"
  echo "[$KIT] cache var recorded: ${HINT#export }"
fi

# CUDA-13 wheel libs on LD_LIBRARY_PATH for later runs (best-effort; per STOCK.md).
CU13="$("$VENV/bin/python" -c 'import site,os; p=os.path.join(site.getsitepackages()[0],"nvidia","cu13","lib"); print(p if os.path.isdir(p) else "")' 2>/dev/null)"
if [ -n "$CU13" ] && ! grep -q "nvidia/cu13/lib" "$VENV/bin/activate" 2>/dev/null; then
  echo "export LD_LIBRARY_PATH=$CU13:\${LD_LIBRARY_PATH:-}" >> "$VENV/bin/activate"
fi

echo "== [$KIT] run.sh check =="
CHECK_RC=1
for args in "--mode fast" "--mode exact" "--mode fast --allow-partial" "--mode exact --allow-partial" ""; do
  if ( cd "$KITDIR" && bash run.sh check --config h100 $args ); then CHECK_RC=0; break; fi
done
[ $CHECK_RC -ne 0 ] && echo "[$KIT] note: run.sh check nonzero (may be partial activation on this config; see above)"
echo "[$KIT] done (check_rc=$CHECK_RC). Weights: $HERE/weights/$KIT"
exit $CHECK_RC
