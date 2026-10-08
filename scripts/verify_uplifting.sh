#!/bin/bash
# Verify installed uplifting kits WITHOUT reinstalling: discover each kit's weights-cache env
# var from its install log, export it, and run `run.sh check` (fast -> exact -> plain). Records
# all cache vars (+ MODEL_OPT_JIT_ROOT) into envs/uplifting.env for stage runners to source.
#
#   scripts/verify_uplifting.sh [kit ...]      # default: every envs/uplifting-* present
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
export MODEL_OPT_JIT_ROOT="$HERE/.jit"
UP="third_party/uplifting-biomolecular-modeling"
ENVFILE="$HERE/envs/uplifting.env"
mkdir -p envs
: > "$ENVFILE"
echo "export MODEL_OPT_JIT_ROOT=$HERE/.jit" >> "$ENVFILE"

KITS="${*:-}"
[ -z "$KITS" ] && KITS="$(ls -d envs/uplifting-*/ 2>/dev/null | sed 's#envs/uplifting-##; s#/##')"

for k in $KITS; do
  venv="envs/uplifting-$k"
  [ -d "$venv" ] || { echo "$k: NO VENV"; continue; }
  hint="$(grep -hoE "export [A-Z][A-Z0-9_]*=[^ )]*weights/$k[^ )]*" "logs/uplifting-$k.log" "logs/.install-$k.out" 2>/dev/null | tail -1)"
  [ -n "$hint" ] && echo "$hint" >> "$ENVFILE"
  rc=1
  for args in "--mode fast" "--mode exact" ""; do
    if ( source "$venv/bin/activate"; [ -n "$hint" ] && eval "$hint"
         cd "$UP/$k" && bash run.sh check --config h100 $args ) >/dev/null 2>&1; then rc=0; break; fi
  done
  if [ $rc -eq 0 ]; then
    echo "$k: ACTIVE ${hint:+(${hint#export })}"
  else
    echo "$k: check nonzero (see logs/uplifting-$k.log)"
  fi
done
echo "--- cache vars recorded -> $ENVFILE (source it before running kits)"
