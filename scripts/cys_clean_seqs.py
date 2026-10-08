#!/usr/bin/env python3
"""Sequence-level cysteine clean (C->S) for refill candidates, so novelty + developability can be
re-scored on the cleaned sequence BEFORE committing GPU to complex refold + MD. C->S is the standard
surface-Cys fix (small, polar, isosteric-ish); the free-Cys survivors get full complex revalidation.

In:  results/refill_candidates.csv (name,sequence,d_occ,occ_HIE,...)
Out: results/refill_candidates_cs.csv (name=<orig>__cs, sequence=C->S, carry d_occ/occ_HIE/lineage, n_cys)
"""
import csv

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
IN = f"{REPO}/results/refill_candidates.csv"
OUT = f"{REPO}/results/refill_candidates_cs.csv"

rows = []
for r in csv.DictReader(open(IN)):
    seq = r["sequence"].upper()
    ncys = seq.count("C")
    rows.append({
        "name": r["name"] + ("__cs" if ncys else ""),
        "orig_name": r["name"],
        "sequence": seq.replace("C", "S"),
        "n_cys_removed": ncys,
        "d_occ": r["d_occ"], "occ_HIE": r["occ_HIE"],
        "clean": r["clean"], "strong": r["strong"], "new_lineage": r["new_lineage"],
        "lineage": r["lineage"],
    })
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["name", "orig_name", "sequence", "n_cys_removed",
                                       "d_occ", "occ_HIE", "clean", "strong", "new_lineage", "lineage"])
    w.writeheader(); w.writerows(rows)
print(f"wrote {OUT}: {len(rows)} cleaned seqs ({sum(1 for r in rows if r['n_cys_removed'])} had Cys)")
