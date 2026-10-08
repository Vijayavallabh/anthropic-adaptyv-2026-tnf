#!/usr/bin/env python3
"""Aggregate generator design metrics into one unified table (stage 2/3 feeder).

Walks a run root for BoltzGen's per-job final_ranked_designs/all_designs_metrics.csv,
tags each row with its source spec/seed/generator, concatenates, and writes one CSV.
Ranking/selection is deliberately NOT done here (see rank/select stage): BoltzGen's own
quality_score favours well-folded natural-like scaffolds (e.g. ubiquitin collapse), which
we must down-rank on novelty. Run with: uv run python scripts/collect_designs.py
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]

UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def ident_to(a: str, b: str) -> float:
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / n if n else 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runroot", default=str(REPO / "runs" / "campaign1"))
    ap.add_argument("--out", default=str(REPO / "results" / "designs_all.csv"))
    ap.add_argument("--generator", default="boltzgen")
    a = ap.parse_args()

    csvs = sorted(glob.glob(
        os.path.join(a.runroot, "**", "final_ranked_designs", "all_designs_metrics.csv"),
        recursive=True))
    if not csvs:
        sys.exit(f"no all_designs_metrics.csv under {a.runroot}")

    frames = []
    for c in csvs:
        try:
            df = pd.read_csv(c)
        except Exception as e:
            print("skip", c, e)
            continue
        rel = Path(c).relative_to(a.runroot).parts
        df["source_spec"] = rel[0] if len(rel) > 0 else "?"
        df["source_seed"] = rel[1] if len(rel) > 1 else "?"
        df["generator"] = a.generator
        df["run_dir"] = str(Path(c).parent.parent)
        frames.append(df)

    alld = pd.concat(frames, ignore_index=True, sort=False)
    seqcol = "designed_chain_sequence" if "designed_chain_sequence" in alld.columns else "designed_sequence"
    alld["binder_seq"] = alld[seqcol].astype(str)
    alld["binder_len"] = alld["binder_seq"].str.len()
    alld["ident_to_ubiquitin"] = alld["binder_seq"].map(lambda s: round(ident_to(s, UBQ), 3))

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    alld.to_csv(a.out, index=False)

    uniq = alld.drop_duplicates(subset=["binder_seq"])
    novel = alld[alld["ident_to_ubiquitin"] < 0.5]
    print(f"{len(csvs)} runs | {len(alld)} designs | {uniq.shape[0]} unique | "
          f"{len(novel)} non-ubiquitin (id<0.5) -> {a.out}")
    if "binder_len" in alld.columns:
        d = alld["binder_len"].describe()
        print(f"length: min={d['min']:.0f} mean={d['mean']:.0f} max={d['max']:.0f}")
    for col in ["design_to_target_iptm", "plip_hbonds_refolded", "plip_saltbridge_refolded",
                "delta_sasa_refolded", "HIS_fraction", "GLU_fraction", "LYS_fraction",
                "designfolding-filter_rmsd", "liability_score"]:
        if col in alld.columns:
            print(f"  {col}: mean={alld[col].mean():.3f} max={alld[col].max():.3f}")
    if "source_spec" in alld.columns:
        print("per-spec counts:", alld["source_spec"].value_counts().to_dict())


if __name__ == "__main__":
    main()
