#!/usr/bin/env python
"""
WAVE 16 merge + correlation: fuse the two Proton-PottsMPNN pH signals into the submission manifest
and test their orthogonality/corroboration against the existing pH-selectivity legs.

Inputs  : results/wave16_potts_phsel.csv  (decoder backend: delta_his433, complex-apo)
          results/wave16_potts_gbind.csv   (Potts-energy backend: ddg_bind_pH two-state)
          results/submission_manifest_v2_wave15.csv  (existing legs: his433_dpKa, pkai_his433_dpKa,
                                                       ens_his433_pka_mean, ddg_pH, ...)
Outputs : results/submission_manifest_v2_wave16.csv  (+ w16_* columns)
          prints Spearman(new Potts signals vs existing legs) + consensus counts.

Run with any python that has pandas+numpy (e.g. the ProtonPottsMPNN venv).
"""
from __future__ import annotations
import os, re, sys
import numpy as np
import pandas as pd

REPO = os.environ.get("ADAPTYV_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = f"{REPO}/results"
PHSEL = f"{R}/wave16_potts_phsel.csv"
GBIND = f"{R}/wave16_potts_gbind.csv"
MANI = f"{R}/submission_manifest_v2_wave15.csv"
OUT = f"{R}/submission_manifest_v2_wave16.csv"


def idx_of(name: str):
    m = re.search(r"(\d{1,2})\s*$", str(name).strip())
    return int(m.group(1)) if m else None


def rankdata(a):
    a = np.asarray(a, float)
    order = a.argsort(kind="mergesort")
    ranks = np.empty(len(a), float)
    ranks[order] = np.arange(1, len(a) + 1)
    # average ties
    _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt)); np.add.at(sums, inv, ranks)
    return (sums / cnt)[inv]


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    if len(x) < 3:
        return float("nan"), int(len(x))
    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan"), int(len(x))
    return float(np.corrcoef(rx, ry)[0, 1]), int(len(x))


def main():
    phsel = pd.read_csv(PHSEL)
    gbind = pd.read_csv(GBIND)
    phsel["idx"] = phsel["name"].map(idx_of)
    gbind["idx"] = gbind["name"].map(idx_of)
    w16 = phsel.merge(gbind[["idx", "gbind_hisP", "gbind_hisS", "ddg_bind_pH", "binder_propensity_hisP"]],
                      on="idx", how="outer")

    mani = None
    if os.path.exists(MANI):
        mani = pd.read_csv(MANI)
        mani["idx"] = mani["name"].map(idx_of)
        if mani["idx"].notna().sum() < len(mani):   # fall back to positional join if names lack an index
            mani = mani.sort_values("name").reset_index(drop=True)
            mani["idx"] = np.arange(1, len(mani) + 1)

    # ---- correlation of the NEW Potts signals vs the EXISTING pH legs ----
    print("==== WAVE16 Potts signals vs existing pH-selectivity legs (Spearman, n) ====")
    existing = {}
    if mani is not None:
        j = w16.merge(mani, on="idx", how="left", suffixes=("", "_m"))
        legs = [c for c in ("his433_dpKa", "pkai_his433_dpKa", "ens_his433_pka_mean", "ddg_pH", "pkai_ddg_pH") if c in j.columns]
        for sig in ("delta_his433", "ddg_bind_pH"):
            for leg in legs:
                rho, n = spearman(j[sig], j[leg])
                print(f"  {sig:14s} vs {leg:22s} rho={rho:+.3f} (n={n})")
            existing[sig] = j[sig]
        # agreement between the two Potts backends (decoder vs potts-energy): expect opposite signs
        rho, n = spearman(j["delta_his433"], j["ddg_bind_pH"])
        print(f"  delta_his433   vs ddg_bind_pH          rho={rho:+.3f} (n={n})  [expect NEGATIVE: both=ON]")
        out = j
    else:
        out = w16
        print("  (manifest not found; skipping cross-leg correlation)")

    # ---- consensus counts ----
    n = len(w16)
    dec_on = w16.delta_his433 > 0
    ene_on = w16.ddg_bind_pH < 0
    w16["w16_potts_consensus"] = np.where(dec_on & ene_on, "ON2",
                                  np.where(dec_on | ene_on, "ON1", "OFF"))
    print("\n==== WAVE16 pH-switch consensus ====")
    print(f"delta_his433 > 0 (decoder: binder shifts His433 -> protonated/ON): {int(dec_on.sum())}/{n}")
    print(f"ddg_bind_pH  < 0 (energy:  binds better protonated = tumour-ON):   {int(ene_on.sum())}/{n}")
    print(f"BOTH Potts backends agree ON (ON2): {int((dec_on & ene_on).sum())}/{n}")
    print(f"at least one backend ON (ON1+ON2): {int((dec_on | ene_on).sum())}/{n}")
    print(f"backends agree on sign (both ON or both OFF): {int(((dec_on & ene_on) | (~dec_on & ~ene_on)).sum())}/{n}")

    # ---- per-design side-by-side: new Potts leg vs existing physics/MD legs ----
    print("\n==== per-design: Potts (new) vs physics legs ====")
    view_cols = ["idx", "delta_his433", "ddg_bind_pH", "w16_potts_consensus"]
    if mani is not None:
        vt = w16.merge(mani, on="idx", how="left")
        for c in ("his433_dpKa", "pkai_his433_dpKa", "ens_his433_pka_mean", "ph_verdict"):
            if c in vt.columns and c not in view_cols:
                view_cols.append(c)
    else:
        vt = w16
    vt = vt.sort_values("delta_his433", ascending=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(vt[[c for c in view_cols if c in vt.columns]].to_string(index=False))

    # ---- write enriched manifest (rename to w16_*) ----
    ren = {"delta_his433": "w16_potts_dHis433_decoder",
           "gap_complex_his433": "w16_potts_gap_complex",
           "gap_apo_his433": "w16_potts_gap_apo",
           "his433_minus_ctrlmean": "w16_potts_his433_vs_ctrl",
           "his433_rank_by_delta": "w16_potts_his433_rank",
           "ddg_bind_pH": "w16_potts_ddg_bind_pH",
           "gbind_hisP": "w16_potts_gbind_hisP",
           "gbind_hisS": "w16_potts_gbind_hisS"}
    keep_cols = ["idx"] + [k for k in ren if k in w16.columns] + (["w16_potts_consensus"] if "w16_potts_consensus" in w16.columns else [])
    keep_idx = w16[keep_cols].rename(columns=ren)
    if mani is not None:
        enriched = mani.merge(keep_idx, on="idx", how="left")
    else:
        enriched = keep_idx
    enriched.to_csv(OUT, index=False)
    print(f"\nwrote {OUT} ({len(enriched)} rows, +{len([v for v in ren.values()])} w16 columns)")


if __name__ == "__main__":
    main()
