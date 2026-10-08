#!/bin/bash
# adaptyv-tnf harness bootstrap + verification (uv-only).
#   ./init.sh          core env + smoke tests (fast, CPU-only) -- the verification gate
#   ./init.sh --tools  additionally scaffold isolated uv venvs for heavy design tools
set -euo pipefail
unset VIRTUAL_ENV 2>/dev/null || true
cd "$(dirname "$0")"

echo "=== adaptyv-tnf :: core environment (uv sync) ==="
uv sync

echo "=== compile check ==="
uv run python -m compileall -q src scripts

echo "=== tests (pytest) ==="
uv run pytest -q

echo "=== submission validator self-check (Track 3) ==="
uv run python scripts/validate_submission.py submissions/example_track3.csv --track 3

if [[ "${1:-}" == "--tools" ]]; then
  echo ""
  echo "=== design tool stacks: Anthropic uplifting-biomolecular-modeling (isolated uv venvs) ==="
  mkdir -p envs third_party weights .jit
  UP="third_party/uplifting-biomolecular-modeling"
  if [ -d "$UP/.git" ]; then
    echo "uplifting-biomolecular-modeling present ($(du -sh "$UP" 2>/dev/null | cut -f1))."
  else
    echo "cloning uplifting-biomolecular-modeling (~1.1 GB, no weights)..."
    git clone --depth 1 https://github.com/anthropics/uplifting-biomolecular-modeling "$UP"
  fi
  echo "export MODEL_OPT_JIT_ROOT=$PWD/.jit   # persistent compile cache"
  echo "Install a kit's env (uv Route C), e.g.:"
  echo "  scripts/setup_uplifting.sh boltz2            # venv only"
  echo "  scripts/setup_uplifting.sh boltz2 --weights  # + weights (several GB)"
  echo "Stage->kit map: configs/uplifting.yaml.  pH scoring (PyRosetta/PypKa): docs/TOOLS.md."
  echo "(Weight installs are NOT auto-run: they are heavy and some are license/gated.)"
fi

echo ""
echo "=== Verification Complete -- environment is clean and restartable ==="
echo "Next steps:"
echo "  1. Read AGENTS.md, progress.md (current status), and docs/PLAN_TNF.md."
echo "  2. The active submission is runs/submission/final_portfolio.{csv,fasta} + METHODS.md;"
echo "     scripts/README.md maps the pipeline. Record evidence before claiming a change done."
echo "  3. Before every upload: uv run python scripts/validate_submission.py <csv> --track 3"
