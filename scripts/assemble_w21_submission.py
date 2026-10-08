#!/usr/bin/env python3
"""Assemble the WAVE-21 recovery submission after 7/20 W19 designs failed the portal novelty check.

Final 20 = 12 portal-CONFIRMED novel real designs + 7 oracle-novel, MD-selective, 0-Cys refills
(names 33-39) + 1 knockout control (32). The 7 failed designs (21,24,25,26,28,30,31) are dropped.

Refill sequences are the cys-cleaned (C->S) sequences that were refolded (af2ig) and re-MD'd, so
submit == folded == MD seq. Order is best-of (strongest/cleanest/most-novel + diverse lineages first,
knockout last); an EVT/DPP reorder can refine it once all metrics are joined.
"""
import csv, sys, subprocess

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
W19 = f"{REPO}/results/adaptyv_challenge1_resubmit_W19_FINAL.csv"
CS = f"{REPO}/results/refill_candidates_cs.csv"
OUT = f"{REPO}/results/adaptyv_challenge1_resubmit_W21.csv"

# new refill name -> parent pdnov name (cleaned seq pulled from refill_candidates_cs.csv).
# All 7 are 0-Cys, oracle-novel, MD-selective (d_occ re-measured on the SUBMITTED cleaned seq).
REFILL = {
    "egfr_phsw_33": "pdnov_egfr_phsw_16__des_184_1",  # d_occ 0.932 (cleaned, UP from 0.685), lin16, safe TM 0.782
    "egfr_phsw_34": "pdnov_egfr_phsw_06__des_248_2",  # d_occ 0.927 (cleaned, UP from 0.633), lin06, safe TM 0.773
    "egfr_phsw_35": "pdnov_egfr_phsw_09__des_148_2",  # d_occ 0.883, lin09, 0-Cys native, safe TM 0.778, bind 0.868
    "egfr_phsw_36": "pd2nov_pdnov_egfr_phsw_03__des_33_1__des_232_1",  # DEEP-SAFE lin03 upgrade (near-twin of portal-confirmed egfr_phsw_29, same 03__33_1 family, oracle TM 0.645 ~ 29's 0.646): 0-Cys, d_occ 0.952, occ_HIE 0.043, af2ig pae 6.24, co-fold 0.835, CLEAN expr, mouse~human. Dominates the prior 09__47 (TM 0.791/d_occ 0.773/bind 0.806) on every axis; 09__47 retained as a bench backup.
    "egfr_phsw_37": "pd2nov_pdnov_egfr_phsw_01__des_196_1__des_168_2",  # DEEP-SAFE lin01 upgrade: 0-Cys, d_occ 0.557, occ_HIE 0.167, oracle TM 0.649 (deep harbor), af2ig pae 5.2, co-fold 0.677, CLEAN expr, mouse~human — improves diversity (lin09 3->2, +lin01) + margin + switch over the prior 09__42. (Backups: pdnov_egfr_phsw_09__des_42_2 C107S safe TM 0.773; deeper-safe bench b09=09__53_1__des_168_2 d_occ 0.962/TM 0.676, b06=06__65_2__des_55_1 d_occ 0.675/TM 0.615; lineage-15 alt pdnov_egfr_phsw_15__des_139_2 TM 0.848.)
    "egfr_phsw_38": "pdnov_egfr_phsw_12__des_122_1",  # d_occ 0.55 (cleaned), lin12, safe TM 0.788
    "egfr_phsw_39": "pdnov_egfr_phsw_16__des_102_1",  # d_occ 0.485, lin16, 0-Cys native, safe TM 0.765
}

# best-of order: strongest/cleanest + diverse (lineage) switches first; knockout 32 last. The one
# novelty-RISKY pick (37 = new-lineage 15, oracle TM 0.848) is parked at rank 15 (tail, outside the
# likely-tested band) so a portal-novelty miss costs nothing. An EVT/DPP reorder can refine this.
ORDER = [
    "egfr_phsw_29",  # confirmed, d_occ 0.927, CLEAN, chai 0.808 — clean+confirmed+binds, lin03
    "egfr_phsw_35",  # refill, d_occ 0.883, bind 0.868, lin09
    "egfr_phsw_23",  # confirmed, d_occ 0.957, chai 0.821, lin06
    "egfr_phsw_33",  # refill, d_occ 0.932, lin16
    "egfr_phsw_27",  # confirmed, d_occ 0.978 (top switch), lin08  -> diverse top-5: lin 03/09/06/16/08
    "egfr_phsw_22",  # confirmed, d_occ 0.962, lin01
    "egfr_phsw_34",  # refill, d_occ 0.927, lin06
    "egfr_phsw_36",  # refill, d_occ 0.773, bind 0.806, lin09
    "egfr_phsw_04",  # locked, d_occ 0.89, MODERATE, chai 0.781, lin04
    "egfr_phsw_38",  # refill, d_occ 0.55, lin12, CLEAN
    "egfr_phsw_08",  # locked, d_occ 0.57, CLEAN, lin08
    "egfr_phsw_39",  # refill, d_occ 0.485, bind 0.611, lin16, CLEAN
    "egfr_phsw_19",  # locked, d_occ 0.48, CLEAN, lin19
    "egfr_phsw_12",  # locked, d_occ 0.335, CLEAN, lin12
    "egfr_phsw_37",  # refill, NEW lineage 15, d_occ 0.578, bind 0.856 — novelty-RISKY (portal-verify), parked
    "egfr_phsw_16",  # locked, d_occ 0.195, CLEAN, lin16
    "egfr_phsw_02",  # locked, d_occ 0.155, CLEAN, lin02
    "egfr_phsw_06",  # locked, d_occ 0.595, SEVERE expr -> deep tail, lin06
    "egfr_phsw_09",  # locked, d_occ 0.255, SEVERE expr -> deep tail, lin09
    "egfr_phsw_32",  # knockout control (E66Q of 23) — last (best-of: non-selective control costs nothing)
]


def main():
    w19 = {r["name"]: r["sequence"].strip().upper() for r in csv.DictReader(open(W19))}
    cs = {r["orig_name"]: r["sequence"].strip().upper() for r in csv.DictReader(open(CS))}
    # deep-safe bench designs (0-Cys, no cleaning needed) pulled from the W20 switch pool screen
    bench = {r["name"]: r["sequence"].strip().upper()
             for r in csv.DictReader(open(f"{REPO}/results/w20_bench_in.csv"))}

    seqs = {}
    for name in ORDER:
        if name in REFILL:
            src = REFILL[name]
            seqs[name] = cs.get(src) or bench.get(src)
            if not seqs[name]:
                sys.exit(f"[FATAL] no sequence for refill {name} <- {src}")
        elif name in w19:
            seqs[name] = w19[name]
        else:
            sys.exit(f"[FATAL] no sequence for {name}")

    # hard validation
    assert len(ORDER) == 20, f"expected 20, got {len(ORDER)}"
    assert len(set(ORDER)) == 20, "duplicate names"
    assert len(set(seqs.values())) == 20, "duplicate sequences"
    dropped = {"egfr_phsw_21", "egfr_phsw_24", "egfr_phsw_25", "egfr_phsw_26",
               "egfr_phsw_28", "egfr_phsw_30", "egfr_phsw_31"}
    assert not (set(ORDER) & dropped), f"a dropped novelty-fail is still present: {set(ORDER)&dropped}"
    for n, s in seqs.items():
        if "C" in s:
            sys.exit(f"[FATAL] {n} contains a cysteine (not 0-Cys): {s}")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for name in ORDER:
            w.writerow([name, seqs[name], "single_chain"])
    print(f"wrote {OUT}: 20 designs, all 0-Cys, {len(set(seqs.values()))} distinct")
    print(f"  refills (new): {sorted(REFILL)}")
    print(f"  dropped novelty-fails: {sorted(dropped)}")

    # sanity: official validator
    print("\n=== validate_submission.py (track 3) ===")
    subprocess.run([f"{REPO}/envs/md/bin/python", f"{REPO}/scripts/validate_submission.py", OUT, "--track", "3"])


if __name__ == "__main__":
    main()
