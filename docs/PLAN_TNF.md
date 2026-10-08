# Challenge 2 (TNF-α) — campaign plan & outcome

Deadline 2026-10-11 23:59 AoE (= 2026-10-12 11:59 UTC). Track 3, ≤20 designs, Claude-selected on
quality + novelty + method-novelty. Objective priority: **pH-selective (bind 7.4 / release 6.0) > mouse
cross-reactive > affinity.** GPU budget: 3 reserved (0,1,2) + any verified-free others.

## Target & epitope

Human TNF-α soluble trimer, PDB **1TNF** (chains A–C, res 6–157), UniProt P01375; mouse ortholog P06804
(soluble domains ~79.5% id). The bindable, cross-react-friendly site is the **cb1** conserved
aromatic/hydrophobic pocket (W114/Y115/Y141 + a conserved interface Lys) at the protomer interface — all
anchors conserved human↔mouse (`configs/target.yaml`). The flat subunit-interface face (Region-I) and the
K90 face **wall** de-novo binding; cb1 is the one epitope that breaks the wall.

## What actually ran (as executed)

1. **Generation — RFdiffusion3 (all-atom).** RFd1 + BoltzGen floored at af2ig PAE ~22 on this flat polar
   target (×2 epitopes, ~6000 designs). RFd3 all-atom at cb1 → ProteinMPNN (0-Cys) broke it: af2ig
   pae_interaction 3.6–4.8 at ~7% hit. (Re-MPNN RFd3 — self-seq folds poorly in AF2.)
2. **Binding validation.** AF2 initial-guess on the full-atom two-protomer target (pae_int<10, pLDDT_binder>80,
   rmsd<2) + Boltz-2 co-fold ipTM 0.89–0.94. Full-atom target is mandatory (backbone-only → false 0%).
3. **pH-switch install.** Interface His next to a **conserved** target Lys → His⁺ repels Lys⁺ at pH 6 →
   release. The cb1 family also carries a **cooperative 2-His** trio (Baker ≥2-His rule → strongest linkage, ddg to −2.34).
4. **Switch validation — 4 method-classes, contact-gated.** Explicit two-state MD (OpenMM, Amber14+GBn2) with a
   Gln-isostere causal control; a **contact gate** (ON-state His–Lys < 8 Å, d_switch > 2·max SEM incl. GLN) that
   rejects non-contact drift artifacts; PROPKA–Wyman linkage; pKAI linkage; and GFN2-xTB buried-limit interaction
   energy (on the six single-His designs). (Rosetta `pH_mode`/`fa_elec` is dead — it does not titrate; see METHODS.)
5. **Cross-reactivity.** Matched human (1TNF) / mouse (P06804-threaded) Boltz-2 co-fold + worst-case-over-species
   **ipSAE > 0.5**; switch anchors are conserved Lys (cross-react-preserving).
6. **Expression / developability.** 0-Cys, instability/GRAVY, SolubleCaliby; no N-glyc/deamidation sequons.
7. **Diversity & novelty.** Foldseek TM-clustering → one design per true shape (backbone-id ≠ shape; the bindable
   surface clusters to ~7 portal-novel folds — the honest ceiling). Final novelty = the Proteinbase portal.
8. **Select + submit.** `scripts/assemble_v2.py` (contact-gated, cross-react-robust, portal-gated, cluster-deduped).

## Outcome

**`runs/submission/final_portfolio.{csv,fasta}` = 9 designs / 6 distinct folds** (byte-verified FASTA==CSV):
the **cb1** family (4 sequences, 1 fold: TNFpH01 single-His dSw 2.68 + a coop-2His trio) + **5 cb1-epitope**
de-novo binders across folds C1/C2/C5/C8/C11. All nine: 3-method switch (MD + PROPKA + pKAI; +xTB on the six
single-His), ipSAE > 0.5, 0-Cys, **portal-novel ≥3/4 (one 4/4)**. `runs/submission/METHODS.md` is the write-up.

**Upload is USER-GATED** (Proteinbase, account j_v_v_07).
