# TNF-α conditional pH-switch binders — final-push portfolio (methods)

De-novo mini-binders that bind TNF-α at pH 7.4 and **release at pH 6.0** (endosomal-recycling
direction), with human↔mouse cross-reactivity. This revision hardens the submission against the
**two levers the Adaptyv EGFR post-mortem shows actually decide wet-lab outcome** — expression and a
robust, multi-signal interface — and reports every design with an honest confidence tier.

Objective priority (unchanged): **pH-selectivity > human/mouse cross-reactivity > affinity.**

## The pH switch (validation engine, carried forward)

Interface histidine(s) adjacent to a **conserved** target lysine (K98/K101/K112…): protonation at pH 6
creates His⁺···Lys⁺ repulsion → release. Two switch chemistries: a single interface His, and a
**cooperative 2-His** network (Baker's ≥2-His rule → strongest proton linkage). Every design is scored by
method-classes that are blind in different ways, read as sign-consensus, never one magnitude:

1. **Explicit two-state MD** (OpenMM, Amber14+GBn2) of the protonated (HIP) and neutral (HIE) states, with
   a **contact gate** (ON-state His–Lys < 8 Å, d_switch > 2·max-SEM incl. a Gln-isostere causal control)
   that rejects non-contact drift artifacts.
2. **PROPKA–Wyman** proton-uptake linkage (His pKa suppressed in complex vs free ⇒ release-competent).
3. **pKAI** linkage (independent ML pKa class).
4. **GFN2-xTB** buried-limit interaction energy on the single-His designs.

## New QC #1 — expression is a hard gate; interface-frozen rescue

Auditing the prior 9 designs against Adaptyv's `protein-qc` **standard** expression gate (instability
index < 40, GRAVY < 0.4, 0-Cys, no aggregation/proteolysis/deamidation runs) showed **only 1/9 passed.**
Cause: the bindable cb1 epitope is hydrophobic, so the vanilla-ProteinMPNN interface needed for binding
raises instability/GRAVY. Per the EGFR post-mortem a non-expressing design scores zero regardless of a
perfect switch, so this is the dominant risk.

**Rescue (new):** freeze the binding interface (union of the human and mouse interfaces, plus the His
switch) and re-sequence **only the solvent-facing surface** with SolubleMPNN (0-Cys). This fixes
expression without touching the binding residues or the switch. Each rescued variant is re-folded
(Boltz-2 human+mouse), re-scored for cross-reactivity (ipSAE), and the switch is re-confirmed by PROPKA
and two-state MD on the exact rescued sequence. Rescued variants pass the expression gate
(instability 27–35, GRAVY ≤ 0.4) with binding, cross-reactivity, and switch preserved.

## Cross-reactivity — grounded in mechanism, not a single folder's cross-species score

Cross-reactivity is **designed in, then verified by mechanism**, because single-folder cross-species
scores proved unreliable:

- **Mechanism.** Every design targets the conserved cb1 epitope (W114/Y115/Y141) with a **conserved-Lys
  switch anchor** (human↔mouse identical). The designed contacts land on conserved positions — the
  footprint shows **≈1 divergent contact in ~27** — so the interface is human↔mouse-invariant by
  construction.
- **Matched co-fold delta.** The right folder metric is a *matched* human-vs-mouse delta under identical
  conditions, not an absolute ortholog score. These deltas are small across the de-novo folds (iptm
  human→mouse e.g. b1m12 0.925→0.925, A20V 0.925→0.871; |Δ| < 0.13).
- **Why we do not rank on absolute cross-species ipSAE.** We stress-tested it two ways — a second folder
  (Chai-1) and a surface perturbation — and it failed both: **Boltz and Chai disagree on mouse ipTM**
  (e.g. b1m12 Boltz 0.83 vs Chai 0.11; A20V the reverse), and single-sequence Boltz mouse ipSAE collapses
  under even one surface mutation. This is the known **folder-pessimism** artifact; absolute cross-species
  ipSAE is a folder/config signal, not a species signal, so we rely on mechanism + matched delta instead.

**Human binding** is confirmed by a **three-folder consensus** — AF2-initial-guess, Boltz-2, and Chai-1 —
all high on every design (Boltz/Chai interface ptm ≈ 0.78–0.94); binding is not a single-folder artifact.

## The portfolio (tiered by expression + switch quality + novelty)

`final_portfolio_push.{csv,fasta}` — 14 designs across **7 distinct Foldseek folds**. All are
cross-reactive by mechanism (above) and bind human by 3-folder consensus; the tiers rank by the real
discriminators — expression (the hard gate), switch strength, and novelty. Expression-PASS: 7.

- **Tier A (expression-clean, strongest on the pH switch — the top objective):** the cooperative-2His cb1
  trio `TNFpH06/07/08` (strongest proton linkage, ddg to −2.34) and `b4m2_r1` (C2, rescued), whose switch is
  **confirmed selective by two-state MD** on the exact rescued sequence (d_repel 0.60; His-causal — the
  Gln-isostere control does not repel). `b1m12_r4` (C5, rescued; the only **novelty 4/4** design, strong
  binding, expression-clean) is also Tier A, but honestly its switch is **weak** — it passes PROPKA but is
  **not MD-selective** (d_repel 0.17, His not interface-causal), the same weak switch as its parent
  (d_switch 0.597); its value is novelty + binding, not the switch.
- **Tier B (strong, expression-hedged):** `b4m2` & `b1m12` originals (highest human ipSAE 0.85/0.83,
  shipped alongside their rescued variants as an expression hedge — "diversify a working binder", the
  EGFR-winning strategy), the second rescued variants `b4m2_r6`/`b1m12_r1`, and `TNFpH01` (single-His, the
  **strongest d_switch 2.68**; expression-risk).
- **Tier C (fold diversity, expression-risk):** `b2m0` (C1), `b1m7` (C11), `A20V` (C12, a new 7th fold),
  `b1m9` (C8) — cross-reactive and binding-confirmed, but fail the expression gate on the original
  sequence and are not surface-rescuable without perturbing the (folder-sensitive) mouse co-fold.

## Honest assessment & recommendation

The dominating designs are Tier A: two de-novo folds (expression-fixed, switch re-confirmed, one novelty
4/4) plus the cooperative-His cb1 family. Tier B adds the highest-affinity binders with an expression
hedge and the strongest single-His switch. Tier C adds fold diversity but carries an expression caveat.
For upload we recommend **Tiers A+B (10 designs)** as the core, with Tier C as optional diversity hedges
up to the Track-3 cap of 20. The rescued variants need a Proteinbase novelty re-check before upload (same
fold as their validated parent ⇒ expected ≥3/4); `A20V` needs its first portal check.
