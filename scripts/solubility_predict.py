#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solubility_predict.py -- Protein-Sol scaled-solubility predictor (offline, sequence-based).

WHAT THIS IS
    A faithful, dependency-free Python reimplementation of the published Protein-Sol
    sequence-based "scaled solubility" algorithm. It adds an ORTHOGONAL,
    experimental-data-derived solubility/expression signal that is independent of
    structure predictors such as ESMFold (it uses only amino-acid composition and
    biophysical property scales fit to an E. coli expression-solubility dataset).

MODEL / PAPER
    Hebditch M, Carballo-Amador MA, Charonis S, Curtis R, Warwicker J.
    "Protein-Sol: a web tool for predicting protein solubility from sequence."
    Bioinformatics 33(19):3098-3100 (2017).  doi:10.1093/bioinformatics/btx345
    Official software: https://protein-sol.manchester.ac.uk  (download tab,
    "protein-sol-sequence-prediction-software", GPLv3, Warwicker & Hebditch 2017).

TRAINING DATA (experimental, orthogonal to ESM)
    eSOL dataset of E. coli proteins solubilised in a *cell-free* expression system,
    Niwa T. et al., PNAS 106:4201-4206 (2009). Protein-Sol uses 3049 proteins (a
    further 653 integral membrane proteins are excluded) and defines a high-solubility
    set (TOP, N=145, mean %solubility 113.2) and a low-solubility set (LOW, N=125,
    mean %solubility 5.2); population mean %solubility = 53.3.

ALGORITHM (transcribed verbatim from the official Perl source:
    server_prediction_seq_export.pl, seq_compositions_perc_pipeline_export.pl,
    seq_reference_data.txt, ss_propensities.txt)
    35 whole-sequence features are computable, but only 10 carry non-zero weight in
    the published fit (the others' weights are 0). For each of the 10 used features f
    with whole-sequence value x, with reference averages xL (LOW set) and xT (TOP set):
        1. clamp x into the interval spanned by [xL, xT] (either orientation);
        2. map linearly to a predicted %solubility:
               pred_f = LOW_PROP + (TOP_PROP - LOW_PROP) * (x_clamped - xL)/(xT - xL)
           with LOW_PROP = 5.208 (LOW mean %sol), TOP_PROP = 113.241 (TOP mean %sol).
    The per-feature predictions are combined by a weighted mean, weights = |z-score
    difference| between the TOP and LOW sets (REF_FEATURES below):
        percent_sol = sum_f( w_f * pred_f ) / sum_f( w_f )
    and reported as a 0-1 "scaled solubility":
        scaled_sol  = (percent_sol - LOW_PROP) / (TOP_PROP - LOW_PROP)

SCALE / INTERPRETATION  (solubility_score = scaled_sol)
    Range 0-1, HIGHER = MORE soluble. Reference population average =
        (53.347 - 5.208) / (113.241 - 5.208) = 0.4456  (~0.45, the published value).
    A score > 0.45 => predicted more soluble than the average eSOL protein; < 0.45 =>
    less soluble. NOTE: this predicts *recombinant expression* solubility, NOT amyloid
    aggregation kinetics -- amyloidogenic peptides (e.g. Abeta, polyQ) can still score
    high, and very short peptides are inflated by the strong length term.

REPORTED PERFORMANCE
    Hebditch et al. (2017) report the predictor separates the high- vs low-solubility
    eSOL populations (scaled-solubility population means ~1.0 vs ~0.0 by construction;
    the fit is a linear model over the 10 features). Protein-Sol is a standard,
    widely-cited sequence solubility baseline. It is a biophysical/compositional model
    (no deep network), which is exactly what makes it orthogonal to ESM-based tools.

FAITHFULNESS / VALIDATION
    This reimplementation reproduces the official Perl tool's `percent-sol`, `scaled-sol`,
    `population-sol` and `pI` outputs to 3 decimal places on all tested controls and on
    the 20 project designs (see scripts companion check; population-sol = 0.446 exactly).

USAGE
    python -I scripts/solubility_predict.py INPUT.csv [-o results/solubility_TAG.csv]
                                            [--tag TAG] [--selftest]
    INPUT.csv must have columns: name,sequence  (extra columns are ignored).
    Run with `-I` when INPUT may be an untrusted directory (this module imports only
    the Python standard library, so `-I` is safe).

Author: generated for the anthropic-adaptyv-2026-tnf project (orthogonal solubility QC).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

MODEL_NAME = "Protein-Sol_scaled_solubility_v2017 (Hebditch et al. 2017, Bioinformatics; eSOL/Niwa E. coli expression data)"

# --- population / set reference %-solubility (seq_reference_data.txt: POP/LOW/TOP AVG) ---
POPULATION_AVG = 53.347   # POP NIWA AVG percentage-solubility
LOW_PROP       = 5.208    # LOW NIWA AVG percentage-solubility (N=125)
TOP_PROP       = 113.241  # TOP NIWA AVG percentage-solubility (N=145)
POPULATION_SCALED = (POPULATION_AVG - LOW_PROP) / (TOP_PROP - LOW_PROP)  # = 0.4456 ~ 0.45

# --- Kyte-Doolittle hydropathy (JMB 1982 157:105), used normalised (x+4.5)/9 in FoldIndex ---
KD = {
    'A':  1.8, 'C':  2.5, 'D': -3.5, 'E': -3.5, 'F':  2.8, 'G': -0.4, 'H': -3.2,
    'I':  4.5, 'K': -3.9, 'L':  3.8, 'M':  1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5,
    'R': -4.5, 'S': -0.8, 'T': -0.7, 'V':  4.2, 'W': -0.9, 'Y': -1.3,
}
# --- net-charge contribution per residue (protein_Q): D,E=-1; K,R=+1; His counted 0 ---
CHG = {'D': -1.0, 'E': -1.0, 'K': 1.0, 'R': 1.0}

AA20 = "ACDEFGHIKLMNPQRSTVWY"

# --- pI calculation constants (protein_pI subroutine) ---
PI_PK  = {'D': 4.0, 'E': 4.4, 'H': 6.3, 'K': 10.4, 'R': 12.0, 'Y': 10.1}
PI_FAC = {'D': -1.0, 'E': -1.0, 'H': 1.0, 'K': 1.0, 'R': 1.0, 'Y': -1.0}

# --- the 10 features that carry non-zero weight in the published fit ---
#   name : (xL = LOW-set average, xT = TOP-set average, weight = |z-score diff|)
#   values transcribed from seq_reference_data.txt (LOW/TOP 'average' rows and the
#   ZDF 'zscore-diff' row; only the 10 flagged use_for_prob=y features appear here).
REF_FEATURES = {
    "K_minus_R_pct": (-2.651, 1.042, 1.096),   # 100*(nK-nR)/L
    "length_naa":    (437.096, 188.055, 1.523), # sequence length
    "His_pct":       (2.923, 2.259, 0.484),
    "Leu_pct":       (10.676, 9.202, 0.557),
    "Val_pct":       (6.001, 7.112, 0.522),
    "AspGlu_pct":    (11.395, 13.657, 0.876),   # 100*(nD+nE)/L
    "aromatic_FWY_pct": (8.918, 6.447, 0.931),  # 100*(nF+nW+nY)/L
    "abs_charge":    (0.016, 0.047, 1.105),     # |mean net charge per residue|
    "fold_index":    (0.150, 0.089, 0.711),     # 2.785*KDnorm_mean - abs_charge - 1.151
    "entropy":       (4.113, 3.995, 1.382),     # Shannon entropy (base 2) of aa composition
}

# The official pipeline writes whole-sequence features to an intermediate CSV at fixed
# printf precision (seq_compositions_perc_pipeline_export.pl) and the prediction step
# (server_prediction_seq_export.pl) reads those ROUNDED values back. We replicate that
# storage precision so the score is byte-identical to the Perl tool:
#   percentages / K-R / D+E / F+W+Y -> 2 dp (%6.2f/%7.2f); length -> integer (%5.0f);
#   abs-charge / fold-index / entropy -> 3 dp (%8.3f).
STORAGE_DECIMALS = {
    "K_minus_R_pct": 2, "length_naa": 0, "His_pct": 2, "Leu_pct": 2, "Val_pct": 2,
    "AspGlu_pct": 2, "aromatic_FWY_pct": 2, "abs_charge": 3, "fold_index": 3, "entropy": 3,
}


def _store(value: float, decimals: int) -> float:
    """Round exactly as C/Perl printf("%.<decimals>f") would (round-half-to-even)."""
    return float(f"{value:.{decimals}f}")


def clean_sequence(seq: str) -> str:
    """Uppercase, drop whitespace. Returns only standard 20-aa characters kept in order.
    (The project sequences are all standard; non-standard letters are dropped, matching
    the composition subroutines which only tally the 20 standard amino acids.)"""
    seq = "".join(seq.split()).upper()
    return "".join(c for c in seq if c in AA20)


def compute_features(seq: str) -> dict:
    """Compute the 10 whole-sequence Protein-Sol features (+ intermediate values)."""
    L = len(seq)
    if L == 0:
        raise ValueError("empty sequence after cleaning")
    counts = {aa: seq.count(aa) for aa in AA20}

    # composition percentages
    def pct(*aas):
        return 100.0 * sum(counts[a] for a in aas) / L

    k_minus_r = 100.0 * (counts['K'] - counts['R']) / L
    his_pct   = pct('H')
    leu_pct   = pct('L')
    val_pct   = pct('V')
    de_pct    = pct('D', 'E')
    fwy_pct   = pct('F', 'W', 'Y')

    # mean net charge per residue -> abs charge
    q_mean = sum(CHG.get(aa, 0.0) * counts[aa] for aa in AA20) / L
    abs_charge = abs(q_mean)

    # normalised Kyte-Doolittle mean, (KD+4.5)/9 per residue
    kdnorm_mean = sum(((KD[aa] + 4.5) / 9.0) * counts[aa] for aa in AA20) / L
    fold_index = 2.785 * kdnorm_mean - abs_charge - 1.151

    # Shannon entropy (base 2) of the composition
    entropy = 0.0
    for aa in AA20:
        if counts[aa]:
            f = counts[aa] / L
            entropy -= f * (math.log(f) / math.log(2.0))

    return {
        "K_minus_R_pct": k_minus_r,
        "length_naa": float(L),
        "His_pct": his_pct,
        "Leu_pct": leu_pct,
        "Val_pct": val_pct,
        "AspGlu_pct": de_pct,
        "aromatic_FWY_pct": fwy_pct,
        "abs_charge": abs_charge,
        "fold_index": fold_index,
        "entropy": entropy,
    }


def compute_pI(seq: str) -> float:
    """Reproduce Protein-Sol's protein_pI subroutine (descending-pH scan, step 0.5,
    linear interpolation of the zero-charge crossing)."""
    aanum = {aa: 0 for aa in PI_PK}
    ionisable = 0
    for aa in seq:
        if aa in aanum:
            aanum[aa] += 1
            ionisable += 1

    def charge_at(pH):
        q = 0.0
        for aa, n in aanum.items():
            if n:
                fac = PI_FAC[aa]
                q += n * (fac * (1.0 / (1.0 + 10 ** (fac * (pH - PI_PK[aa])))))
        return q

    loop = 0
    pH = 14.0
    pH_dec = 0.5
    protein_charge = -100.0
    protein_charge_last = protein_charge
    while protein_charge < 0 and pH > 0.9:
        loop += 1
        protein_charge_last = protein_charge
        protein_charge = charge_at(pH)
        pH -= pH_dec

    if ionisable == 0:
        return 7.0
    if loop == 1:
        return 14.0
    if loop == 27:
        return 1.0
    return (pH + pH_dec) + pH_dec * (protein_charge / (protein_charge - protein_charge_last))


def predict(seq: str) -> dict:
    """Return percent_sol, scaled_sol (solubility_score), pI and the feature dict.

    Features are rounded to the official tool's intermediate-storage precision before
    being fed to the weighted-interpolation prediction (see STORAGE_DECIMALS)."""
    raw = compute_features(seq)
    feats = {k: _store(v, STORAGE_DECIMALS[k]) for k, v in raw.items()}
    pred_sum = 0.0
    wt_sum = 0.0
    for name, (xL, xT, w) in REF_FEATURES.items():
        x = feats[name]
        # clamp to the interval spanned by [xL, xT] (handle either orientation)
        lo, hi = (xL, xT) if xL <= xT else (xT, xL)
        xc = min(max(x, lo), hi)
        pred_f = LOW_PROP + (TOP_PROP - LOW_PROP) * ((xc - xL) / (xT - xL))
        pred_sum += w * pred_f
        wt_sum += w
    percent_sol = pred_sum / wt_sum
    scaled_sol = (percent_sol - LOW_PROP) / (TOP_PROP - LOW_PROP)
    return {
        "percent_sol": percent_sol,
        "scaled_sol": scaled_sol,
        "pI": compute_pI(seq),
        "features": feats,
    }


# feature columns emitted as sub-scores (length is also emitted as sequence_length)
SUBSCORE_COLS = [
    "K_minus_R_pct", "His_pct", "Leu_pct", "Val_pct", "AspGlu_pct",
    "aromatic_FWY_pct", "abs_charge", "fold_index", "entropy",
]


def read_input(path: str):
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        if rdr.fieldnames is None:
            raise SystemExit(f"error: {path} is empty")
        cols = {c.lower().strip(): c for c in rdr.fieldnames}
        if "name" not in cols or "sequence" not in cols:
            raise SystemExit(f"error: {path} must have columns name,sequence (got {rdr.fieldnames})")
        rows = []
        for r in rdr:
            name = (r[cols["name"]] or "").strip()
            raw = (r[cols["sequence"]] or "").strip()
            if not name and not raw:
                continue
            rows.append((name, raw))
    return rows


def run(input_csv: str, out_csv: str):
    rows = read_input(input_csv)
    results = []
    for name, raw in rows:
        seq = clean_sequence(raw)
        p = predict(seq)
        rec = {
            "name": name,
            "sequence_length": len(seq),
            "solubility_score": round(p["scaled_sol"], 3),
            "percent_solubility": round(p["percent_sol"], 3),
            "population_reference": round(POPULATION_SCALED, 3),
            "pI": round(p["pI"], 2),
            "model_name": MODEL_NAME,
        }
        for c in SUBSCORE_COLS:
            rec[c] = round(p["features"][c], 3)
        results.append(rec)

    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    fieldnames = ["name", "sequence_length", "solubility_score", "percent_solubility",
                  "population_reference", "pI", "model_name"] + SUBSCORE_COLS
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for rec in results:
            w.writerow(rec)

    # ranked table (worst -> best solubility), printed to stdout
    ranked = sorted(results, key=lambda r: r["solubility_score"])
    ref = round(POPULATION_SCALED, 3)
    print(f"\nModel: {MODEL_NAME}")
    print(f"Population reference (avg eSOL protein) scaled solubility = {ref}  "
          f"(>{ref} = more soluble than average)\n")
    print(f"{'rank':>4}  {'name':<34} {'len':>4} {'scaled_sol':>10} {'%sol':>7} {'pI':>5}  flag")
    print("-" * 78)
    for i, r in enumerate(ranked, 1):
        flag = "LOW" if r["solubility_score"] < ref else ""
        print(f"{i:>4}  {r['name']:<34} {r['sequence_length']:>4} "
              f"{r['solubility_score']:>10.3f} {r['percent_solubility']:>7.2f} {r['pI']:>5.2f}  {flag}")
    print(f"\nwrote {out_csv}  ({len(results)} sequences)")
    return results


# ---- built-in self-test: reproduces the published population reference and a few
#      hand-checkable invariants (full Perl cross-validation is done separately). ----
def selftest():
    ok = True
    exp = 0.446
    got = round(POPULATION_SCALED, 3)
    print(f"[selftest] population_scaled = {got} (expected {exp})  "
          f"{'PASS' if got == exp else 'FAIL'}")
    ok &= (got == exp)
    # GB1 domain should be clearly soluble (>0.45); known official value 0.899
    gb1 = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
    s = round(predict(gb1)["scaled_sol"], 3)
    print(f"[selftest] GB1 scaled_sol = {s} (official Perl 0.899)  "
          f"{'PASS' if s == 0.899 else 'FAIL'}")
    ok &= (s == 0.899)
    print("[selftest]", "ALL PASS" if ok else "FAILURE")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Protein-Sol scaled-solubility predictor (offline).")
    ap.add_argument("input_csv", nargs="?", help="input CSV with columns name,sequence")
    ap.add_argument("-o", "--out", help="output CSV path (default results/solubility_<tag>.csv)")
    ap.add_argument("--tag", help="tag for default output filename")
    ap.add_argument("--selftest", action="store_true", help="run built-in self-test and exit")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest()
    if not args.input_csv:
        ap.error("input_csv is required (or use --selftest)")

    out = args.out
    if not out:
        tag = args.tag or os.path.splitext(os.path.basename(args.input_csv))[0]
        out = os.path.join("results", f"solubility_{tag}.csv")
    run(args.input_csv, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
