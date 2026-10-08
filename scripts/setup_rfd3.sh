#!/bin/bash
# Native Route C setup for the RFdiffusion3 (rfd3) uplifting kit — the one kit the generic
# scripts/setup_uplifting.sh cannot do (it needs Python 3.12.1 EXACTLY, two interpreters,
# the foundry stock tarball, and a locally-compiled NVIDIA apex wheel). Implements the exact
# sequence of rfdiffusion3/STOCK.md "Route C" and environment/Dockerfile, entirely in USER
# SPACE (no sudo, no Docker): the host driver already supports CUDA 13.0 (nvidia-smi), torch's
# cu130 wheels carry the CUDA runtime, and the only thing that needs nvcc 13.0 is building the
# apex wheel once — provisioned here via a throwaway conda CUDA-13 toolkit prefix.
#
#   scripts/setup_rfd3.sh [--arch LIST] [--jobs N] [--skip-weights] [--roles "stock kit"]
#
# Everything lives on /home (3.2 TB free); nothing writes to / (97% full). Not `set -e`: each
# phase gates explicitly so a failure is named, not hidden.
set -uo pipefail
unset VIRTUAL_ENV 2>/dev/null || true
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"

ARCH="9.0+PTX"          # H100 (sm90) only; this box is 8x H100. add "8.0;" for A100.
JOBS="$(nproc --all)"; [ "$JOBS" -gt 16 ] && JOBS=16
DO_WEIGHTS=1
ROLES="stock kit"
while [ $# -gt 0 ]; do
  case "$1" in
    --arch) ARCH="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --skip-weights) DO_WEIGHTS=0; shift ;;
    --roles) ROLES="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

KITDIR="$HERE/third_party/uplifting-biomolecular-modeling/rfdiffusion3"
COMMON="$HERE/third_party/uplifting-biomolecular-modeling/common"
[ -d "$KITDIR" ] || { echo "FATAL: kit not found at $KITDIR (run ./init.sh --tools)"; exit 1; }

# --- /-safe caches + build scratch, all on /home ---
BUILD="$HERE/.rfd3build"
export TMPDIR="$BUILD/tmp"
export UV_CACHE_DIR="${UV_CACHE_DIR:-$HOME/.cache/uv}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$HOME/.cache/pip}"
export MODEL_OPT_JIT_ROOT="$HERE/.jit"
mkdir -p "$TMPDIR" "$HERE/.jit" "$HERE/envs" "$HERE/weights/rfd3" "$HERE/logs" "$KITDIR/stock/wheels"

CONDA="$HOME/miniconda3/bin/conda"
CUDA13="$BUILD/cuda13"          # throwaway CUDA-13 toolkit prefix (BUILD-TIME ONLY: nvcc + dev headers)
PINS="$KITDIR/stock/PINS.json"

say(){ echo ""; echo "=========================================================================="; echo "== $*"; echo "=========================================================================="; }
die(){ echo "FATAL: $*" >&2; exit 1; }

# Read pins (apex commit/wheel, foundry version) straight from PINS.json.
read -r APEX_COMMIT APEX_WHEEL FOUNDRY_VER < <(python3 - "$PINS" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
a=p["pinned_stack"]["apex"]
print(a["commit"], a["wheel"], p["upstream"]["foundry"]["version"])
PY
)
[ -n "${APEX_COMMIT:-}" ] || die "could not read apex commit from $PINS"
echo "pins: apex=$APEX_COMMIT wheel=$APEX_WHEEL foundry=$FOUNDRY_VER  | arch=$ARCH jobs=$JOBS roles='$ROLES'"
df -h /home | tail -1

# =========================================================================
say "PHASE 1 — CUDA 13.0 build toolkit (user-space conda prefix)"
# =========================================================================
if [ -x "$CUDA13/bin/nvcc" ] && "$CUDA13/bin/nvcc" --version 2>/dev/null | grep -q "release 13"; then
  echo "CUDA-13 toolkit already present at $CUDA13 — skipping"
else
  [ -x "$CONDA" ] || die "conda not found at $CONDA"
  echo "creating $CUDA13 (cuda-toolkit 13.0 from the nvidia channel; a few GB, one-time)..."
  "$CONDA" create -y -p "$CUDA13" -c nvidia -c conda-forge \
      "cuda-version=13.0" cuda-nvcc cuda-cudart-dev cuda-nvrtc-dev cuda-nvml-dev \
      cuda-profiler-api cuda-crt cuda-cccl libcublas-dev libcurand-dev libcusparse-dev libcusolver-dev \
      > "$HERE/logs/rfd3_cuda13_create.log" 2>&1 \
    || { echo "--- tail conda log ---"; tail -30 "$HERE/logs/rfd3_cuda13_create.log"; die "conda CUDA-13 create failed"; }
fi
"$CUDA13/bin/nvcc" --version | tail -2 || die "nvcc not runnable in $CUDA13"
# Locate CUDA headers (conda may place them under targets/x86_64-linux/include).
CU_INC=""
for d in "$CUDA13/include" "$CUDA13/targets/x86_64-linux/include"; do
  [ -f "$d/cuda_runtime.h" ] && CU_INC="$d" && break
done
[ -n "$CU_INC" ] || die "cuda_runtime.h not found under $CUDA13 (dev headers missing)"
echo "CUDA headers: $CU_INC"

# =========================================================================
say "PHASE 2 — build the pinned stack into the venv(s)"
# =========================================================================
# lock minus the two non-index lines (installed separately): apex (compiled) + rc-foundry (tarball)
LOCK="$KITDIR/environment/requirements.lock"
STACK="$BUILD/stack.txt"
[ -s "$LOCK" ] || die "requirements.lock missing/empty: $LOCK (kit half-extracted or path wrong)"
grep -v -E '^(#|apex==|rc-foundry==)' "$LOCK" > "$STACK"
[ -s "$STACK" ] || die "filtered pin stack is empty: $STACK (requirements.lock unreadable?)"
# build backends for the foundry sdist (Dockerfile pins these)
BUILDPINS="$BUILD/build-pins.txt"
printf 'hatchling==1.32.0\nhatch-vcs==0.4.0\nsetuptools-scm==10.2.1\nvcs-versioning==2.3.1\n' > "$BUILDPINS"

declare -A VENV
VENV[stock]="$HERE/envs/uplifting-rfd3-stock"
VENV[kit]="$HERE/envs/uplifting-rfd3-kit"

for role in $ROLES; do
  [ -n "${VENV[$role]:-}" ] || die "unknown --roles value '$role' (valid: stock kit)"
done

for role in $ROLES; do
  V="${VENV[$role]}"
  say "PHASE 2.$role — venv $V (Python 3.12.1 exact)"
  uv venv --clear --seed --managed-python --python 3.12.1 "$V" || die "[$role] uv venv failed"
  PY="$V/bin/python"
  "$PY" --version | grep -q "3.12.1" || die "[$role] venv python is not 3.12.1: $("$PY" --version)"
  echo "[$role] installing pinned stack (torch 2.13.0+cu130, triton 3.7.1, ...; several GB, --no-deps)"
  uv pip install --python "$PY" --no-deps -r "$STACK" \
      > "$HERE/logs/rfd3_stack_$role.log" 2>&1 \
    || { tail -25 "$HERE/logs/rfd3_stack_$role.log"; die "[$role] stack install failed"; }
  echo "[$role] installing stock RFdiffusion3 (rc-foundry $FOUNDRY_VER from foundry tarball, --no-deps)"
  SETUPTOOLS_SCM_PRETEND_VERSION="$FOUNDRY_VER" PIP_CONSTRAINT="$BUILDPINS" \
    "$PY" -m pip install --no-cache-dir --no-deps "$KITDIR/stock/foundry-4010e3e.tar.gz" \
      > "$HERE/logs/rfd3_foundry_$role.log" 2>&1 \
    || { tail -25 "$HERE/logs/rfd3_foundry_$role.log"; die "[$role] foundry install failed"; }
  "$PY" -c "import torch,triton,rfd3; print('[$role] torch',torch.__version__,' triton',triton.__version__,' rfd3 import OK')" \
    || die "[$role] post-install import (torch/triton/rfd3) failed"
done

# =========================================================================
say "PHASE 3 — build the NVIDIA apex wheel (once, needs nvcc 13.0)"
# =========================================================================
WHEEL_PATH="$KITDIR/stock/wheels/$APEX_WHEEL"
if [ -f "$WHEEL_PATH" ]; then
  echo "apex wheel (pinned $APEX_WHEEL) already present — skipping build"
else
  OTHER_WHL="$(ls "$KITDIR"/stock/wheels/apex-0.1-*.whl 2>/dev/null | head -1)"
  [ -n "$OTHER_WHL" ] && echo "WARN: apex wheel $OTHER_WHL present but != pinned $APEX_WHEEL; rebuilding to the pin (ABI/CUDA-arch safety)" >&2
  # the stock venv's torch is the ABI target for the wheel
  BV="${VENV[stock]}"; [ -d "$BV" ] || BV="${VENV[kit]}"
  BPY="$BV/bin/python"
  echo "[apex] installing ninja==1.13.0 into $BV"
  uv pip install --python "$BPY" --no-deps ninja==1.13.0 >/dev/null 2>&1 || echo "[apex] WARN: ninja pin install reported errors"
  SRC="$BUILD/apex-src"
  rm -rf "$SRC"; mkdir -p "$SRC"
  TARBALL="$BUILD/apex-$APEX_COMMIT.tar.gz"
  echo "[apex] fetching source @ $APEX_COMMIT"
  if ! curl -fsSL -o "$TARBALL" "https://github.com/NVIDIA/apex/archive/$APEX_COMMIT.tar.gz"; then
    echo "[apex] tarball fetch failed; falling back to git clone"
    git clone https://github.com/NVIDIA/apex.git "$SRC/apex" >/dev/null 2>&1 \
      && git -C "$SRC/apex" checkout "$APEX_COMMIT" >/dev/null 2>&1 \
      || die "[apex] could not obtain source"
    APEX_DIR="$SRC/apex"
  else
    tar -xzf "$TARBALL" -C "$SRC" || die "[apex] extract failed"
    APEX_DIR="$SRC/apex-$APEX_COMMIT"
  fi
  [ -f "$APEX_DIR/setup.py" ] || die "[apex] no setup.py in $APEX_DIR"
  echo "[apex] compiling (CUDA_HOME=$CUDA13, arch=$ARCH, MAX_JOBS=$JOBS, system gcc $(gcc -dumpversion)) — ~10 min"
  ( cd "$APEX_DIR" && \
    PATH="$CUDA13/bin:/usr/bin:$HOME/.local/bin:$PATH" \
    CUDA_HOME="$CUDA13" CPATH="$CU_INC${CPATH:+:$CPATH}" \
    CC=/usr/bin/gcc CXX=/usr/bin/g++ NVCC_PREPEND_FLAGS="-ccbin /usr/bin/g++" \
    APEX_CPP_EXT=1 APEX_CUDA_EXT=1 MAX_JOBS="$JOBS" TORCH_CUDA_ARCH_LIST="$ARCH" \
    "$BPY" -m pip wheel --no-deps --no-build-isolation --no-cache-dir -w "$KITDIR/stock/wheels" . \
  ) > "$HERE/logs/rfd3_apex_build.log" 2>&1 \
    || { echo "--- tail apex build log ---"; tail -40 "$HERE/logs/rfd3_apex_build.log"; die "[apex] wheel build failed (see logs/rfd3_apex_build.log)"; }
  WHEEL_PATH="$(ls "$KITDIR"/stock/wheels/apex-0.1-*.whl 2>/dev/null | head -1)"
  [ -n "$WHEEL_PATH" ] || die "[apex] build finished but no wheel in stock/wheels/"
fi
echo "[apex] wheel: $WHEEL_PATH"
sha256sum "$WHEEL_PATH" | awk '{print "[apex] sha256",$1}'

# =========================================================================
say "PHASE 4 — install apex + the kit (run.sh install) in each interpreter"
# =========================================================================
for role in $ROLES; do
  V="${VENV[$role]}"; PY="$V/bin/python"
  say "PHASE 4.$role"
  echo "[$role] installing apex wheel (--no-deps)"
  "$PY" -m pip install --no-cache-dir --no-deps --force-reinstall "$WHEEL_PATH" \
      > "$HERE/logs/rfd3_apexinst_$role.log" 2>&1 \
    || { tail -20 "$HERE/logs/rfd3_apexinst_$role.log"; die "[$role] apex install failed"; }
  "$PY" -c "import apex, importlib; importlib.import_module('apex.normalization.fused_layer_norm'); import fused_layer_norm_cuda; print('[$role] apex fused_layer_norm_cuda OK')" \
    || die "[$role] apex extension import failed"
  echo "[$role] run.sh install (editable opt_core + rfdiffusion3_opt, then pin check)"
  if [ "$role" = kit ] && [ "$DO_WEIGHTS" -eq 1 ]; then
    ( cd "$KITDIR" && PATH="$V/bin:$PATH" bash run.sh install --weights "$HERE/weights/rfd3" ) 2>&1 | tee "$HERE/logs/rfd3_install_$role.log"
  else
    ( cd "$KITDIR" && PATH="$V/bin:$PATH" bash run.sh install ) 2>&1 | tee "$HERE/logs/rfd3_install_$role.log"
  fi
  [ "${PIPESTATUS[0]}" -eq 0 ] || die "[$role] run.sh install failed (see logs/rfd3_install_$role.log)"
done

# =========================================================================
say "PHASE 5 — apply the kit overlay (KIT interpreter only)"
# =========================================================================
if echo "$ROLES" | grep -qw kit; then
  KV="${VENV[kit]}"
  ( cd "$KITDIR" && "$KV/bin/python" -m rfdiffusion3_opt.install apply && "$KV/bin/python" -m rfdiffusion3_opt.install check-only ) \
      2>&1 | tee "$HERE/logs/rfd3_overlay.log" \
    || echo "WARN: overlay apply/check reported nonzero (see logs/rfd3_overlay.log)"
fi

# record LD_LIBRARY_PATH hint for cu13 libs (per setup_uplifting convention)
for role in $ROLES; do
  V="${VENV[$role]}"
  CU13LIB="$("$V/bin/python" -c 'import site,os;p=os.path.join(site.getsitepackages()[0],"nvidia","cu13","lib");print(p if os.path.isdir(p) else "")' 2>/dev/null)"
  [ -n "$CU13LIB" ] && ! grep -q "nvidia/cu13/lib" "$V/bin/activate" 2>/dev/null \
    && echo "export LD_LIBRARY_PATH=$CU13LIB:\${LD_LIBRARY_PATH:-}" >> "$V/bin/activate"
done

# =========================================================================
say "PHASE 6 — verify (run.sh check)"
# =========================================================================
STOCK_PY="${VENV[stock]}/bin/python"
for role in $ROLES; do
  V="${VENV[$role]}"
  if [ "$role" = kit ]; then MODES="fast exact"; else MODES="off"; fi
  for m in $MODES; do
    if [ "$role" = kit ] && [ ! -x "$STOCK_PY" ]; then
      echo "[$role/$m] SKIP check: stock interpreter not built ($STOCK_PY); re-run with --roles 'stock kit'" >&2
      continue
    fi
    echo "--- [$role] run.sh check --config h100 --mode $m ---"
    ( cd "$KITDIR" && PATH="$V/bin:$PATH" RFDIFFUSION3_STOCK_PYTHON="$STOCK_PY" bash run.sh check --config h100 --mode "$m" ) \
      && echo "[$role/$m] CHECK OK" || echo "[$role/$m] check nonzero (see output above)"
  done
done

say "DONE — rfd3 Route C. venvs: ${VENV[stock]} , ${VENV[kit]} ; weights: $HERE/weights/rfd3"
echo "Design (kit, fast):   cd $KITDIR && PATH=${VENV[kit]}/bin:\$PATH RFDIFFUSION3_STOCK_PYTHON=${VENV[stock]}/bin/python RFD3_CKPT=$HERE/weights/rfd3/rfd3_latest.ckpt bash run.sh design --config h100 --mode fast inputs=<spec.json> out_dir=<dir>"
