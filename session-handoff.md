# Session handoff — TNF Challenge 2

Updated 2026-10-10 (SESSION 10). Deadline **2026-10-11 23:59 AoE** (2026-10-12 11:59 UTC).

## Current best submission (validated, code-reviewed)

- **`runs/push/submission_v9.csv`** (md5 `05a344444cf740907e79a177d34c73f5`) — 20 designs, **all 20
  portal-novelty ≥3/4 (3 at 4/4)**, FINAL-TEST pass (20 unique names AND sequences). **Recommended upload.**
- **Fallbacks (both fully portal-verified):** `submission_v8.csv` (md5 `8566b133`) and `submission_v7.csv`
  (md5 `0e84b392`, 20 designs / original floor).
- **Upload is USER-GATED** — do not auto-upload any submission. CSVs live under `runs/` (gitignored); the
  `scripts/integrate_v*.py` scripts are the committed source of truth. Rebuild v9:
  `.venv/bin/python -I scripts/integrate_v9.py`.
- Portal CSV format: columns `name,sequence,molecule_class`; de-novo minis are `single_chain`; length 10–250;
  unique names AND sequences. Account `j_v_v_07`.

## What v9 is

**7 de-novo multi-His pH-switches** (Baker 2025 SOTA recipe: ≥2–3 SEPARATED His each repelling a DISTINCT
conserved target cation; triple-confirmed sign PROPKA `dpKa`<0 + GFN2-xTB `dE`>0) across **5 lineages**, plus
**13 portal-verified diverse v7 floor** = **16 lineages**. Switch core (confirmed-His / min-ipSAE / Δ_hm /
pLDDT / novelty):
- `des_3050_6_c39` 4 / 0.666 / 0.063 / 0.88 / **4/4** — strongest switch (wave-4b His-biased)
- `des_4098_1_b14` 3 / 0.587 / 0.106 / 0.90 / 3/4
- `des_4098_1_c44` 2 / 0.861 / 0.024 / 0.82 / 3/4 — highest binding + balance
- `des_180_6_H9-14-38` 2 / 0.734 / 0.058 / 0.81 / 3/4
- `des_818_4_H16-74` 2 / 0.713 / 0.105 / 0.88 / 3/4
- `des_3050_6_H52-54-79` 2 / 0.631 / 0.095 / 0.79 / 3/4
- `des_875_5_H3-57-80` 1 / 0.600 / 0.016 / 0.85 / **4/4** — most novel + best balance

Provenance + metrics for all 20: `runs/push/master_v9.csv`; portal scores: `runs/push/portal_novelty_v9_full.json`.

## Open items

- **Upload decision (user):** upload `submission_v9.csv` (or v8/v7). Nothing is uploaded.
- Optional: draft the submission methods write-up (offered). The high-value compute is exhausted (strong
  balanced-scaffold pool ≈10 lineages fully mined; 24+ triple-confirmed multi-His switches; cap-2 diversity).

## Key invariants (see AGENTS.md + docs/METHODS.md §8)

- **pH-switch = ≥2–3 SEPARATED His on DISTINCT conserved target cations** (not clustered). Prove the sign per
  His by triple-confirm (PROPKA `dpKa`<0 AND xTB `dE`>0; design OK if ≥2 confirmed), on the exact submitted
  sequence. Functional k_off is wet-lab-only — d_iface≈0 in unbiased MD is EXPECTED, not a failure; the
  thermodynamic proton-linkage is the honest in-silico ceiling. pH 6.0 is harder than Baker's 5.4.
- **Novelty is a hard structural gate (<3/4 = DQ)**; portal-check the final set before locking. MPNN re-seqs
  on an already-verified backbone inherit its novelty (backbone-determined). Local foldseek proxy is
  relative-rank only (no AFDB).
- **Cross-react** by folding BOTH orthologs + balanced min-ipSAE/Δ_hm (not a single folder, not footprint count).
- **Expression** ranks on ESMFold MONOMER pLDDT (primary) + SAP + Glu%; ignore ProtParam instability + ML
  seq-solubility tools (OOD). Consensus-drop only on ≥2 orthogonal axes.
- Commit source/docs/harness only; `runs/ results/ envs/ third_party/ tools/ weights/` are local-only.
