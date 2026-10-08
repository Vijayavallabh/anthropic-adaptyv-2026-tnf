#!/usr/bin/env python3
"""Merge pH scores + design metrics, rank by competition priority (pH-sel > cross-react > affinity),
select a diverse top-N, and write a submission CSV (name,sequence,molecule_class), best-first.

Inputs:
  --filtered  results/filtered_designs.csv            (from filter_designs.py: design metadata)
  --variants  stage2 manifest csv                     (from stage2_mutate.py: design,tag,variant_pdb,...)
  --ph65 / --ph74  stage4 csvs at pH 6.5 / 7.4        (from stage4_ph_score.py: pdb,path,ddg,dsasa,his433_state)
Output: a submission CSV; validate afterwards with scripts/validate_submission.py.

pH_margin = ddg(6.5) - ddg(7.4); MORE NEGATIVE = binds better at low pH = tumor-ON (objective 1).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
D3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}


def binder_seq_from_pdb(pdb: str) -> str:
    """Smaller chain's CA sequence (the binder)."""
    chains: dict[str, list] = {}
    try:
        for line in open(pdb):
            if line[:4] != "ATOM":
                continue
            if line[12:16].strip() != "CA":
                continue
            ch = line[21]
            chains.setdefault(ch, []).append(D3.get(line[17:20].strip(), "X"))
    except OSError:
        return ""
    if not chains:
        return ""
    binder = min(chains, key=lambda c: len(chains[c]))
    return "".join(chains[binder])


def base(p):
    return os.path.basename(str(p))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--filtered", default=str(REPO / "results" / "filtered_designs.csv"))
    ap.add_argument("--variants", default=str(REPO / "results" / "variants.csv"))
    ap.add_argument("--ph65", default=str(REPO / "results" / "ph_6.5.csv"))
    ap.add_argument("--ph74", default=str(REPO / "results" / "ph_7.4.csv"))
    ap.add_argument("--out", default=str(REPO / "results" / "submission_track3.csv"))
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--name-prefix", default="egfr_ph")
    a = ap.parse_args()

    var = pd.read_csv(a.variants)
    p65 = pd.read_csv(a.ph65)
    p74 = pd.read_csv(a.ph74)
    for d in (var, p65, p74):
        d["key"] = d[[c for c in d.columns if c in ("variant_pdb", "pdb", "path")][0]].map(base)
    p65 = p65.rename(columns={"ddg": "ddg65", "dsasa": "dsasa65", "his433_state": "his65"})
    p74 = p74.rename(columns={"ddg": "ddg74"})
    m = var.merge(p65[["key", "ddg65", "dsasa65", "his65"]], on="key", how="left")
    m = m.merge(p74[["key", "ddg74"]], on="key", how="left")
    m["ddg65"] = pd.to_numeric(m["ddg65"], errors="coerce")
    m["ddg74"] = pd.to_numeric(m["ddg74"], errors="coerce")
    m["ph_margin"] = m["ddg65"] - m["ddg74"]
    m["binder_seq"] = m[[c for c in m.columns if c in ("variant_pdb", "pdb")][0]].map(binder_seq_from_pdb)
    m["binder_len"] = m["binder_seq"].str.len()

    # merge design-level metrics (cross-react, novelty)
    try:
        filt = pd.read_csv(a.filtered)
        keep = [c for c in ["design", "n_conserved_contacts", "n_divergent_contacts",
                            "glu496_contact", "ident_to_ubiquitin"] if c in filt.columns]
        m = m.merge(filt[keep], on="design", how="left")
    except Exception:
        pass

    # gates: valid sequence, actually binds at 6.5 (ddg65<0), a working switch (margin<0)
    m = m[m["binder_seq"].str.len().between(10, 250)]
    m = m[m["binder_seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", na=False)]
    cand = m[(m["ddg65"] < 0) & (m["ph_margin"] < 0)].copy()
    pool = cand if len(cand) >= a.n else m.copy()   # fall back to full pool if too few pass gates

    # rank: pH-sel (margin asc) > affinity (ddg65 asc) > cross-react (conserved desc)
    cons = pool["n_conserved_contacts"] if "n_conserved_contacts" in pool.columns else 0
    pool = pool.assign(_cons=cons).sort_values(
        by=["ph_margin", "ddg65", "_cons"], ascending=[True, True, False])

    # diversify: one best variant per design, then unique sequences
    pool = pool.drop_duplicates(subset=["design"]).drop_duplicates(subset=["binder_seq"])
    top = pool.head(a.n)

    rows = []
    for i, (_, r) in enumerate(top.iterrows(), 1):
        rows.append({"name": f"{a.name_prefix}_{i:02d}_{r['design'][:30]}",
                     "sequence": r["binder_seq"], "molecule_class": "protein"})
    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)

    print(f"candidates passing gates: {len(cand)} | pool: {len(pool)} | selected: {len(out)} -> {a.out}")
    show = top[["design", "tag", "ph_margin", "ddg65", "ddg74", "binder_len"]].head(a.n)
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
