#!/usr/bin/env python3
"""Dispatch one pipeline stage.

Stage 0 runs in the core uv env. Stages 1-6 need the tool stacks installed by
``./init.sh --tools`` (see docs/TOOLS.md) and are embarrassingly parallel across the
8x H100 (one generator/predictor process per GPU).

Usage:
    uv run python scripts/run_stage.py 0
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

STAGES = {
    0: "stage0_target_prep",
    1: "stage1_generate",
    2: "stage2_sequence",
    3: "stage3_fold_score",
    4: "stage4_ph_crossreact",
    5: "stage5_validate",
    6: "stage6_select",
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Run an adaptyv-tnf pipeline stage.")
    ap.add_argument("stage", type=int, choices=sorted(STAGES), help="pipeline stage 0-6")
    a = ap.parse_args(argv)
    mod = importlib.import_module(f"adaptyv_tnf.{STAGES[a.stage]}")
    mod.main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
