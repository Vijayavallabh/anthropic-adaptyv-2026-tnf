# scripts/ — manifest

This repo was bootstrapped from the challenge-agnostic conditional-binder toolkit (the `anthropic-adaptyv-2026-egfr`
repo). As a result `scripts/` holds three kinds of file:

1. **The active TNF Challenge-2 pipeline** (listed below) — the scripts that produce `runs/submission/final_portfolio.{csv,fasta}`.
2. **Reusable utilities** — generic helpers (PDB cleanup, mutation, renumbering, pKa, novelty, cross-reactivity) used across challenges.
3. **Legacy EGFR Challenge-1 exploration** — mechanism/theozyme/MD scripts targeting the EGFR His433 switch. These are **not** part of the TNF submission path; they are kept for reference and are recognizable by `His433`/EGFR naming. (The 20 EGFR final-selection/ranking/metrics scripts with hardcoded `egfr_phsw_NN` design lists were removed — they live in the `anthropic-adaptyv-2026-egfr` repo.)

All scripts resolve their repo root dynamically (from `__file__` / `$0`), so none depends on an absolute path or a sibling repo.

## Active TNF Challenge-2 pipeline

**Generation (de-novo, all-atom):** `run_gen_tnf.sh`, `partial_tnf.sh`, `partial_cb1div.sh`, `partial_family.sh`,
`run_campaign.sh`, `run_long_campaign.sh`, `pipe_epitope2.sh`, `eval_rfd3.sh`, `eval_regionI.sh` (RFdiffusion3 at the
conserved cb1 epitope → ProteinMPNN).

**Binding validation:** `seq_af2ig_tnf.sh` (AF2 initial-guess), `run_boltz_tnf.sh` (Boltz-2 co-fold),
`build_fullatom_af2in_tnf.py`, `cb1ep2_select_scan.py`.

**pH-switch install + validation (the centerpiece — 4 method-classes):**
`scan_switch.py` (install His at a conserved target Lys) · `md_release_tnf.py` (explicit two-state MD + Gln-isostere control) ·
`ph_release_score.py` + `ph_release_consensus.py` (PROPKA–Wyman + pKAI linkage; reuses `obj1_pkai.py`) ·
`switch_contact_gate.py` (ON-state contact gate) · `qm_switch_fair.py` + `run_qm_switch.py` (GFN2-xTB buried-limit energy) ·
`merge_signals.py`, `finalize_cb1ep_batch.py`, `cb1ep_md_manifest.py` (verdicts → `confirmed.tsv`).

**Cross-reactivity:** `thread_mouse_tnf.py`, `make_mouse_complex.py`, `mouse_mutate.py` (mouse TNF model + per-ortholog Boltz/ipSAE).

**Assembly + novelty:** `assemble_v2.py` (earlier contact-gated, cluster-deduped assembler → `final_portfolio.{csv,fasta}`),
`fold_diversity.py` (Foldseek TM-clustering), `build_core8.py` (cb1 family floor), `tnf_novelty.sh`.

**De-novo MULTI-His pH-switch campaign (v8/v9 — CURRENT submission path, Baker 2025 SOTA; see METHODS §8):**
`wave3_hisscan_build.py` (Track A: point-mutate ≥2–3 His combos, each on a DISTINCT conserved cation) ·
`wave3b_prep.py` + `wave3b_mpnn.sh` + `wave3b_analyze.py` (Track B: LigandMPNN redesign the interface shell
AROUND fixed His so it adapts) · `wave4b.sh` + `wave4b_analyze.py` + `wave4b_triple.py` (His-biased deeper run) ·
`wave3_fold.sh` (generic Boltz human+mouse fold + ipSAE, args: rundir tagpfx) · `wave3_contacts2.py`
(His→distinct-cation 8Å gate) · `wave3_triple.py` / `wave3b_triple.py` (per-His PROPKA `dpKa` + GFN2-xTB `dE`
sign confirm) · `esmfold_monomer.py` (single-seq MONOMER pLDDT expression ranker) · `sap_score.py`
(aggregation) · **`integrate_v9.py`** (assemble switch core + diverse portal-verified floor →
`runs/push/submission_v9.csv`, the current best).

**Submission contract:** `validate_submission.py` + `src/adaptyv_tnf/submission.py` (shared Adaptyv 2026 upload format;
TNF = Track 3, ≤20 designs — caps in `configs/filters.yaml`). The final upload is to the Proteinbase portal (user-gated).
