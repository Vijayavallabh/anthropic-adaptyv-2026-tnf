#!/usr/bin/env python3
"""Assemble the refill candidate pool: MD-selective WAVE-19 designs that are NOT already in
the live submission, with sequences, to be scored by the calibrated novelty oracle.

Context: 7 of the 20 submitted designs failed the portal novelty check (2/4). All 7 are
partial-diffusion children that stayed too close to their parent backbone. We need MD-selective
replacements that ALSO clear novelty. The W19 MD screen ran 93 designs (d_occ, occ_HIE, clean,
strong flags); only 11 were promoted. This pulls the rest as refill candidates.

Output: results/refill_candidates.csv (name, sequence, d_occ, occ_HIE, clean, strong,
new_lineage, lineage, pae_interaction, plddt_binder) sorted by d_occ desc.
"""
import csv, os

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
SCREEN = f"{REPO}/results/md_screen_hits_w19.csv"
FASTA = f"{REPO}/runs/novelty_w18.fasta"
SUBMISSION = f"{REPO}/results/adaptyv_challenge1_resubmit_W19_FINAL.csv"
OUT = f"{REPO}/results/refill_candidates.csv"
DOCC_FLOOR = 0.40  # keep anything plausibly selective; novelty filter + final rank come later

# The 11 parents already promoted to submitted designs 21-31 (exclude: do not re-pick).
EXCLUDE = {
    "pdnov_egfr_phsw_11__des_340_1", "pdnov_egfr_phsw_14__des_44_1", "pdnov_egfr_phsw_20__des_66_2",
    "pdnov_egfr_phsw_11__des_97_1", "pdnov_egfr_phsw_03__des_33_1", "pdnov_egfr_phsw_01__des_196_1",
    "pdnov_egfr_phsw_06__des_65_2", "pdnov_egfr_phsw_14__des_171_2", "pdnov_egfr_phsw_08__des_149_1",
    "pdnov_egfr_phsw_06__des_157_2", "pdnov_egfr_phsw_07__des_387_1",
}


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"):
            name = l[1:].strip().split()[0]
        elif name:
            out[name] = out.get(name, "") + l.strip()
    return out


def main():
    seqs = read_fasta(FASTA)
    submitted_seqs = {r["sequence"].strip().upper() for r in csv.DictReader(open(SUBMISSION))}

    rows, missing_seq, dup_seq = [], 0, 0
    for r in csv.DictReader(open(SCREEN)):
        n = r["name"]
        if n in EXCLUDE:
            continue
        try:
            docc = float(r["d_occ"])
        except (ValueError, KeyError):
            continue
        if docc < DOCC_FLOOR:
            continue
        s = seqs.get(n, "").upper()
        if not s:
            missing_seq += 1
            continue
        if s in submitted_seqs:
            dup_seq += 1
            continue
        rows.append({
            "name": n, "sequence": s,
            "d_occ": round(docc, 4),
            "occ_HIE": round(float(r["occ_HIE"]), 4),
            "clean": r["clean"], "strong": r["strong"], "new_lineage": r["new_lineage"],
            "lineage": r["lineage"],
            "pae_interaction": r["pae_interaction"], "plddt_binder": r["plddt_binder"],
        })

    rows.sort(key=lambda x: -x["d_occ"])
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "sequence", "d_occ", "occ_HIE", "clean",
                                           "strong", "new_lineage", "lineage",
                                           "pae_interaction", "plddt_binder"])
        w.writeheader(); w.writerows(rows)

    nstrong = sum(1 for r in rows if r["strong"] == "True")
    nclean = sum(1 for r in rows if r["clean"] == "True")
    nnewlin = len({r["lineage"] for r in rows if r["new_lineage"] == "True"})
    print(f"wrote {OUT}: {len(rows)} candidates (excluded {len(EXCLUDE)} submitted parents; "
          f"{missing_seq} missing seq, {dup_seq} dup-of-submission)")
    print(f"  clean={nclean} strong={nstrong} | lineages={sorted({r['lineage'] for r in rows})}")
    print(f"  top by d_occ:")
    for r in rows[:14]:
        print(f"    {r['name']:34s} d_occ={r['d_occ']:.3f} occHIE={r['occ_HIE']:.3f} "
              f"clean={r['clean']:5s} strong={r['strong']:5s} newlin={r['new_lineage']:5s} {r['lineage']}")


if __name__ == "__main__":
    main()
