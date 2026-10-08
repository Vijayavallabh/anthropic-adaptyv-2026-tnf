#!/usr/bin/env python3
"""Guaranteed-all-pass fallback: the 13 designs the Adaptyv portal ALREADY scored 3/4 (2026-10-05),
byte-identical to their W19_FINAL sequences. Every submitted design is portal-confirmed novel, so this
file satisfies "all submitted designs >= 3/4 novelty" with ZERO portal iteration.

This is the zero-risk floor. The recommended upside file is the 20-design W21 (adds 7 oracle-safe-harbor
refills that still need a free portal novelty check). Track 3 allows up to 20; a 13-design all-confirmed
submission is valid and strong (8 de novo + 4 strong-switch partial-diffusion + 1 knockout control).
"""
import csv, subprocess

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
W19 = f"{REPO}/results/adaptyv_challenge1_resubmit_W19_FINAL.csv"
OUT = f"{REPO}/results/adaptyv_challenge1_resubmit_SAFE13.csv"

# the 13 portal-confirmed (3/4) designs, in best-of order: strong+clean switches first,
# severe-expression ones in the tail, knockout last.
ORDER = [
    "egfr_phsw_29",  # d_occ 0.927, CLEAN, chai 0.808 — confirmed, clean, binds
    "egfr_phsw_23",  # d_occ 0.957, chai 0.821
    "egfr_phsw_22",  # d_occ 0.962
    "egfr_phsw_27",  # d_occ 0.978 (top switch)
    "egfr_phsw_04",  # locked de novo, d_occ 0.89, chai 0.781
    "egfr_phsw_08",  # locked de novo, d_occ 0.57, CLEAN
    "egfr_phsw_19",  # locked de novo, d_occ 0.48, CLEAN
    "egfr_phsw_12",  # locked de novo, d_occ 0.335, CLEAN
    "egfr_phsw_16",  # locked de novo, d_occ 0.195, CLEAN
    "egfr_phsw_02",  # locked de novo, d_occ 0.155, CLEAN
    "egfr_phsw_06",  # locked de novo, SEVERE expr -> tail
    "egfr_phsw_09",  # locked de novo, SEVERE expr -> tail
    "egfr_phsw_32",  # knockout control (E66Q of 23) — last
]
PORTAL_CONFIRMED = set(ORDER)  # exactly the 13 at 3/4 in the 2026-10-05 novelty check


def main():
    w19 = {r["name"]: r["sequence"].strip().upper() for r in csv.DictReader(open(W19))}
    assert PORTAL_CONFIRMED <= set(w19), "a confirmed design is missing from W19_FINAL"
    assert len(ORDER) == 13 and len(set(ORDER)) == 13
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["name", "sequence", "molecule_class"])
        for n in ORDER:
            s = w19[n]
            assert "C" not in s, f"{n} has a cysteine"
            w.writerow([n, s, "single_chain"])
    print(f"wrote {OUT}: 13 portal-confirmed (3/4) designs, all 0-Cys, guaranteed all-pass")
    print("\n=== validate_submission.py (track 3) ===")
    subprocess.run([f"{REPO}/envs/md/bin/python", f"{REPO}/scripts/validate_submission.py", OUT, "--track", "3"])


if __name__ == "__main__":
    main()
