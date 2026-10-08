#!/usr/bin/env python3
"""Filter + rank generated designs (from analyze_designs.py output) and export top-N complex PDBs.

Gates (cross-reactivity + developability + novelty), then ranks by a composite that favours
epitope-focused, cross-reactive, well-packed, novel binders. Converts the top-N design CIFs to
clean PDBs (gemmi) for the His-switch (stage 2) + pH scoring (stage 4) pipeline.

Run: uv run python scripts/filter_designs.py --analyzed results/designs_analyzed.csv \
       --topn 150 --outdir runs/filtered --manifest results/filtered_designs.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import gemmi
import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyzed", default=str(REPO / "results" / "designs_analyzed.csv"))
    ap.add_argument("--topn", type=int, default=150)
    ap.add_argument("--outdir", default=str(REPO / "runs" / "filtered"))
    ap.add_argument("--manifest", default=str(REPO / "results" / "filtered_designs.csv"))
    ap.add_argument("--min-len", type=int, default=55)
    ap.add_argument("--max-len", type=int, default=130)
    a = ap.parse_args()

    df = pd.read_csv(a.analyzed)
    n0 = len(df)

    # Gates
    df = df[df["binder_len"].between(a.min_len, a.max_len)]
    df = df[df["his433_contact"] == True]                      # must engage the pH sensor
    df = df[df["n_conserved_contacts"] >= 3]                   # cross-reactive epitope core
    df = df[df["n_divergent_contacts"] <= 1]                   # avoid human/mouse-divergent contacts
    df = df[df["ident_to_ubiquitin"] < 0.60]                   # novelty (de novo, not ubiquitin)
    df = df[df["frac_C"] <= 0.08]                              # developability: limit cysteines
    df = df.drop_duplicates(subset=["binder_seq"])
    n_gated = len(df)

    # Composite rank (higher = better). Obj2-heavy (cross-react) + interface quality + novelty.
    df = df.copy()
    df["score"] = (
        2.0 * df["n_conserved_contacts"]
        - 2.0 * df["n_divergent_contacts"]
        + 0.05 * df["n_interface_target"]
        + 3.0 * df["glu496_contact"].astype(float)            # bonus: also near the conserved acidic partner
        + 2.0 * (1.0 - df["ident_to_ubiquitin"])              # novelty bonus
        - 5.0 * df["frac_C"]                                  # penalise cysteines
    )
    df = df.sort_values("score", ascending=False)

    # Diversify: cap per source_spec so we don't over-pick one variant
    top = df.head(a.topn)

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    conv_fail = 0
    for _, r in top.iterrows():
        cif = r["cif"]
        stem = f"{r['source_spec']}_{r['source_seed']}_{r['design_id']}"
        pdb = outdir / f"{stem}.pdb"
        try:
            s = gemmi.read_structure(cif)
            s.setup_entities()
            pdb.write_text(s.make_pdb_string())
        except Exception as e:
            conv_fail += 1
            continue
        rows.append({
            "design": stem, "pdb": str(pdb), "cif": cif,
            "binder_chain": r["binder_chain"], "target_chain": r["target_chain"],
            "binder_len": r["binder_len"], "binder_seq": r["binder_seq"],
            "n_conserved_contacts": r["n_conserved_contacts"],
            "n_divergent_contacts": r["n_divergent_contacts"],
            "his433_contact": r["his433_contact"], "glu496_contact": r["glu496_contact"],
            "ident_to_ubiquitin": r["ident_to_ubiquitin"], "score": round(r["score"], 2),
            "source_spec": r["source_spec"],
        })
    man = pd.DataFrame(rows)
    Path(a.manifest).parent.mkdir(parents=True, exist_ok=True)
    man.to_csv(a.manifest, index=False)

    print(f"analyzed={n0} -> gated={n_gated} -> exported={len(man)} PDBs (conv_fail={conv_fail}) -> {a.outdir}")
    print(f"manifest: {a.manifest}")
    if len(man):
        print("per-spec counts:", man["source_spec"].value_counts().to_dict())
        print(f"binder_len: min={man['binder_len'].min()} med={int(man['binder_len'].median())} max={man['binder_len'].max()}")
        print(f"mean conserved={man['n_conserved_contacts'].mean():.2f} divergent={man['n_divergent_contacts'].mean():.2f}")


if __name__ == "__main__":
    main()
