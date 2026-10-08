#!/usr/bin/env python3
"""Build results/metrics_w21_final.csv + docs/DESIGN_TABLE_W21.md for the WAVE-21 recovery panel.

20 = 12 portal-CONFIRMED novel real designs (from metrics_w19_final.csv) + 7 oracle-novel, MD-selective,
0-Cys refills (metrics measured this wave) + 1 knockout (32). The 7 W19 novelty-fails are dropped.

Refill metrics are the values verified this wave: d_occ re-measured on the SUBMITTED (cys-cleaned) seq
(salvage MD), novelty = calibrated ESMFold oracle TM (safe harbor <0.80 = 7/7 portal-pass on ground truth),
bind = min(chai_iptm_human, boltz2_iptm), xreact = core-4A mechanism severity (0 = mouse-compatible).
"""
import csv

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
SUB = f"{REPO}/results/adaptyv_challenge1_resubmit_W21.csv"
W19M = f"{REPO}/results/metrics_w19_final.csv"
OUT = f"{REPO}/results/metrics_w21_final.csv"
TABLE = f"{REPO}/docs/DESIGN_TABLE_W21.md"

# refill metrics verified this wave (name -> dict). oracle_tm: ESMFold+TMalign on the submitted seq.
# novelty: "pred-pass" (TM<0.80, 7/7 ground-truth) or "RISKY" (TM in the 0.80-0.86 band -> portal-verify).
REFILL = {
    "egfr_phsw_33": dict(parent="pdnov_egfr_phsw_16__des_184_1 (C53S)", lin="egfr_phsw_16", docc=0.932, occhie=0.02,  tm=0.782, bind=0.534, expr="CLEAN",    nov="pred-pass"),
    "egfr_phsw_34": dict(parent="pdnov_egfr_phsw_06__des_248_2 (C68S)", lin="egfr_phsw_06", docc=0.927, occhie=0.073, tm=0.773, bind=0.452, expr="MODERATE", nov="pred-pass"),
    "egfr_phsw_35": dict(parent="pdnov_egfr_phsw_09__des_148_2",        lin="egfr_phsw_09", docc=0.883, occhie=0.103, tm=0.778, bind=0.868, expr="MODERATE", nov="pred-pass"),
    "egfr_phsw_36": dict(parent="pd2nov_pdnov_egfr_phsw_03__des_33_1__des_232_1", lin="egfr_phsw_03", docc=0.952, occhie=0.043, tm=0.645, bind=0.835, expr="CLEAN", nov="pred-pass (near-twin of portal-confirmed 29)"),
    "egfr_phsw_37": dict(parent="pd2nov_pdnov_egfr_phsw_01__des_196_1__des_168_2", lin="egfr_phsw_01", docc=0.557, occhie=0.167, tm=0.649, bind=0.677, expr="CLEAN", nov="pred-pass"),
    "egfr_phsw_38": dict(parent="pdnov_egfr_phsw_12__des_122_1 (C42S)", lin="egfr_phsw_12", docc=0.55,  occhie=0.38,  tm=0.788, bind=0.436, expr="CLEAN",    nov="pred-pass"),
    "egfr_phsw_39": dict(parent="pdnov_egfr_phsw_16__des_102_1",        lin="egfr_phsw_16", docc=0.485, occhie=0.205, tm=0.765, bind=0.611, expr="CLEAN",    nov="pred-pass"),
}
# all 7 refills are in the oracle safe harbor (TM < 0.80, 7/7 portal-pass on the W19 ground truth).
BACKUP_NOTE = ("All 7 refills are 0-Cys and in the novelty safe harbor (oracle TM 0.645-0.788; six at 0.765-0.791, "
               "egfr_phsw_37 deep-safe at 0.649). Validated deeper-safe BENCH swaps (0-Cys, MD-selective, all from "
               "the W18/W20 switch-pool oracle screen) if the portal fails any refill: b09 = "
               "pd2nov_pdnov_egfr_phsw_09__des_53_1__des_168_2 (d_occ 0.962, TM 0.676, but co-fold binding weak 0.16), "
               "b06 = pd2nov_pdnov_egfr_phsw_06__des_65_2__des_55_1 (d_occ 0.675, TM 0.615, bind 0.444). A lineage-15 "
               "high-upside alternative (oracle TM 0.848, risky band, portal decides): pdnov_egfr_phsw_15__des_139_2.")


def main():
    order = [r["name"] for r in csv.DictReader(open(SUB))]
    sub = {r["name"]: r["sequence"] for r in csv.DictReader(open(SUB))}
    w19 = {r["name"]: r for r in csv.DictReader(open(W19M))}

    cols = ["ship_rank", "name", "tier", "length", "lineage", "d_occ", "occ_HIE",
            "novelty", "oracle_TM", "bind_consensus", "xreact_mech_sev", "expr_risk", "free_cys", "parent", "sequence"]
    rows = []
    for i, name in enumerate(order, 1):
        seq = sub[name]
        if name in REFILL:
            m = REFILL[name]
            rows.append(dict(ship_rank=i, name=name, tier="B_refill_w21", length=len(seq), lineage=m["lin"],
                             d_occ=m["docc"], occ_HIE=m["occhie"], novelty=m["nov"], oracle_TM=m["tm"],
                             bind_consensus=m["bind"], xreact_mech_sev=0, expr_risk=m["expr"], free_cys=0,
                             parent=m["parent"], sequence=seq))
        else:
            r = w19.get(name, {})
            tier = "C_control" if name == "egfr_phsw_32" else "A_portal_confirmed"
            rows.append(dict(ship_rank=i, name=name, tier=tier, length=len(seq), lineage=r.get("lineage", ""),
                             d_occ=r.get("d_occ", ""), occ_HIE=r.get("occ_HIE", ""), novelty="portal 3/4",
                             oracle_TM="", bind_consensus=r.get("chai_human", ""), xreact_mech_sev="",
                             expr_risk=r.get("expr_risk", ""), free_cys=0, parent=r.get("notes", ""), sequence=seq))

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT}: {len(rows)} designs")

    # markdown table
    lines = ["# WAVE 21 design table - novelty-recovery panel (20 designs)", "",
             "12 portal-confirmed novel real designs + 7 oracle-novel 0-Cys MD-selective refills + 1 knockout.",
             "The 7 W19 designs that scored 2/4 on the portal novelty check (21,24,25,26,28,30,31) are removed.", "",
             "| rank | name | tier | len | lineage | d_occ | occ_HIE | novelty | oracle TM | bind | mouse sev | expr | parent |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['ship_rank']} | {r['name']} | {r['tier']} | {r['length']} | {r['lineage']} | "
                     f"{r['d_occ']} | {r['occ_HIE']} | {r['novelty']} | {r['oracle_TM']} | {r['bind_consensus']} | "
                     f"{r['xreact_mech_sev']} | {r['expr_risk']} | {r['parent']} |")
    lines += ["", "## Novelty",
              "Novelty = the calibrated local oracle (`scripts/novelty_oracle.py`): fold the sequence with ESMFold, "
              "Foldseek-TMalign max TM vs PDB+AFDB, exactly as the Adaptyv portal does with ESMFold2. Calibrated on "
              "the 20 W19 portal labels: TM < 0.80 was 7/7 portal-pass (safe harbor); TM >= 0.88 reliably fails. "
              "All seven refills sit in the safe harbor (oracle TM 0.645 to 0.788).",
              "", "## Backup", BACKUP_NOTE,
              "", "## Refill validation (per design, on the SUBMITTED cys-cleaned sequence)",
              "- Switch: 2-state explicit-protonation MD (His433 charged vs neutral), 3 replicas x 2 ns; d_occ = occ_HIP - occ_HIE.",
              "- Binding: AF2-initial-guess interface PAE < 10 + co-fold consensus min(chai_iptm_human, boltz2_iptm).",
              "- Cross-reactivity: core-4A mechanism severity 0 (all tight divergent human/mouse contacts conservative).",
              "- Developability: 0 free Cys (cys-cleaned C->S/C->A where needed, then re-folded + re-MD'd), expression risk CLEAN/MODERATE."]
    open(TABLE, "w").write("\n".join(lines) + "\n")
    print(f"wrote {TABLE}")


if __name__ == "__main__":
    main()
