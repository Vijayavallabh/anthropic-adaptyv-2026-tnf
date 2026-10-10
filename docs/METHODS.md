# Conditional de-novo binder design: methods playbook

A challenge-agnostic record of the pipeline and the validation philosophy for
designing de-novo protein binders with a conditional activity (bind only under a
trigger: a pH window, a cofactor, a conformational state). Swap in the target,
epitope, and trigger per project. The specifics below are stated as transferable
principles, not fixed numbers.

> **This repo's application (TNF-α Challenge 2):** the executed pipeline, thresholds, and the 9-design
> outcome are in [`PLAN_TNF.md`](PLAN_TNF.md); the submission write-up is `runs/submission/METHODS.md`
> and the script map is [`../scripts/README.md`](../scripts/README.md). The TNF pH switch is validated by
> four method-classes (two-state MD + Gln-isostere control, PROPKA–Wyman, pKAI, GFN2-xTB) — **not** Rosetta
> `fa_elec`/`pH_mode`, which does not titrate and is not used.

## 1. Objective ordering

Rank the objectives before designing, and design to the top one. For a
conditional binder the hardest and highest-value property is almost always the
conditionality (selectivity between ON and OFF), not raw affinity. A workable
default order:

1. conditional selectivity (ON vs OFF)
2. any required cross-reactivity (one sequence binds several orthologs)
3. affinity in the ON state

Novelty, when the host scores it, is a gate below a threshold, not a score to
maximise.

## 2. Pipeline

1. **Target prep.** Fetch the structure, define the epitope (the sub-patch to
   block or engage), and if cross-reactivity matters, mark conserved vs
   divergent positions across the orthologs. Build a full-atom target: co-fold
   scoring against a backbone-only target collapses the predicted target plDDT
   and reads as "no binding."
2. **Backbone generation.** De-novo backbones from a diffusion generator.
   Partial diffusion of an already-validated binder passes downstream filters at
   roughly an order of magnitude higher rate than de-novo hotspot sampling, so
   use it to amplify a working lead. To diversify, fold more *original*
   generator sequences rather than re-sequencing one backbone: re-sequencing a
   backbone whose own sequence binds usually folds well but kills binding
   (binding is sequence-specific, not backbone-specific).
3. **Sequence design.** Soluble-MPNN or ProteinMPNN. Aim for zero free cysteines
   and a non-hydrophobic surface from the start; fixing an expression liability
   later costs a full re-validation.
4. **Fold and score binding.** AF2 initial-guess against the full-atom target
   (interface-PAE as the pass gate), plus an orthogonal co-folder (for example
   Boltz-2 and Chai-1). A good fold or ipTM is necessary, not sufficient.
5. **Install and score the switch.** See section 3. This is the decisive step.
6. **Cross-reactivity (if required).** Score by mechanism, not a single
   co-folder. Check the divergent-position contacts are chemically conservative,
   and confirm a *matched* co-fold delta against the ortholog (ortholog score
   minus the repacked-host control), not an absolute ortholog score.
7. **Validate.** Re-run every metric on the exact final sequence, with
   submit == folded == simulated. Strengthening a switch can break binding or
   the geometry, so re-validate the mutated artifact end to end.
8. **Rank and select.** By the objective order. Cap per lineage and per family,
   and guarantee every confirmed backbone family is represented: a diverse
   portfolio beats near-clones of one winning seed. Penalise developability
   liabilities (unpaired Cys, large exposed hydrophobic patches).
9. **Novelty gate.** A fail-loud sequence and structure novelty check with a
   known-positive control, below the host's threshold.

## 3. Scoring a conditional switch (the decisive, dual-confirmed engine)

The central lesson: **a positive static score proves the switch DIRECTION, not
its SELECTIVITY.** A constitutive interaction that sits in the same state at
both trigger values scores positive on a single-structure method yet is active
in both states. So:

- **Do not trust one method or its magnitude.** For a titratable residue, no
  single-structure pKa predictor reliably beats a trivial null. Trust
  cross-method sign-consensus, not any one magnitude.
- **Score the thermodynamic linkage** for a pH switch: integrate the
  proton-uptake difference across the pH window from the predicted pKas of
  complex, target, and binder. The sign gives the direction. Use at least two
  independent pKa methods and clamp obvious buried-charge artifacts.
- **Gate selectivity on explicit-state simulation.** Run the two states (for a
  pH switch, the protonated and neutral forms) under MD and measure the
  interaction-occupancy difference between states. Require the ON-state
  occupancy high AND the OFF-state low (an occupancy delta, not just a positive
  linkage). This is the only test that separates a real switch from a
  constitutive interaction.
- **Confirm** with a two-state binding free-energy delta (MM-GBSA) of the
  correct sign, and an isostere-ablation control: mutate out the switching
  chemistry and the effect should collapse. A knockout design (one mutation that
  removes the switching charge) is a cheap built-in wet-lab control.
- **Validate the whole engine on a known natural switch** (a positive control)
  so the scorer fails loud: it must recover the correct sign on the known case
  and stay silent on a null.

Dead ends to rule out first (verify, do not assume):

- Fixed-charge Rosetta pH mode does not titrate on ref2015 (identical energy
  across pH); the protonated-His type is absent from the residue set.
- Constant-pH MD is fragile to set up and was a NO-GO in practice.
- Burying the titratable group makes its pKa unpredictable to within a couple of
  units, wider than a useful window. Prefer a residue with a naturally in-window
  pKa (histidine for a pH switch) as the sensor, and keep it solvent-facing.
- The off-state partner chemistry governs constitutivity: a buried salt bridge
  reads as constitutive by chemistry even when the linkage looks positive. Score
  the OFF-state interaction energy and reject the ones that stay engaged.

## 4. Rigor rules that repeatedly paid off

- **Fail-loud gates.** Check exit codes, assert a known-positive control, and
  never read an empty tool output as "all pass." A silent tool failure once
  passed a whole gate on nothing.
- **Verify the live artifact by bytes.** Confirm submit == validated == folded
  by diffing bytes, not by trusting a filename; an overloaded name once denoted
  two different files.
- **Audit every gating axis before declaring an optimum.** Developability,
  expression, and liabilities gate every objective even when affinity looks
  maxed. "All levers closed" is usually a missed axis.
- **A small-n win is noise until bootstrapped.** Test significance and decompose
  true vs false positives before adopting a lead, and add a confound-free
  control that isolates the claimed mechanism. Resampling-robustness is
  necessary but not sufficient.
- **Validate a lever in the target regime,** not an adjacent one. A setting that
  helps an unrestrained or generic case can hurt the restrained, production case.
- **Select every hand-set constant jointly** (or pre-register them). A win that
  holds only at one hand-fixed radius or cutoff is a nuisance-parameter
  cherry-pick.

## 5. Reproduce

Core env via uv (`./init.sh`): CPU-only, light (target prep, numbering, liability
scan, submission build and validation, tests). Heavy models install into
isolated per-tool uv venvs (see `docs/TOOLS.md`). The pipeline lives in
`scripts/` and the stage skeleton plus submission contract in `src/`. The
submission gate is `scripts/validate_submission.py`; the novelty gate is
`scripts/novelty_check.sh` (fail-loud, with a positive control).

## 6. De-novo β-pairing switch maturation (2026-10 — the strongest designs)

Breaks the de-novo wall on the flat/polar TNF target by backbone-directing a binder β-strand to pair with
the target's exposed edge strand, then installs and triple-confirms the pH switch. Pipeline:

1. **Generate:** β-pairing RFdiffusion (SS+ADJ fold-conditioning, `--mode off`, binder-first contig,
   edge-strand hotspots) — `scripts/run_bpair_rfd.sh`; interface tensors via
   `scripts/make_interface_tensor_nopyr.py` (biotite, pyrosetta-free). ~2400 backbones.
2. **Sequence + filter:** vanilla ProteinMPNN → af2ig gate (pae<10, plddt_binder>85, rmsd<2) → unique passers.
3. **Cross-react:** Boltz-2 fold each passer vs human AND mouse mature TNF; score min-ipSAE (worst ortholog).
   Top de-novo min-ipSAE 0.54–0.70 (meets/beats best partial-diffusion).
4. **Install switch:** His-scan each binder interface position whose CB is <~8 Å of a CONSERVED target Lys
   in the af2pred bound pose; re-fold human+mouse; keep His variants that preserve cross-react
   (`scripts/bpair_hisscan_build.py` + `scripts/bpair_mature.sh`).
5. **Confirm (decisive gate):** brutal two-state MD (`scripts/bpair_md_prep.py` + `bpair_md_launch.sh`;
   `md_release_tnf.py` runs HIP/HIE/GLN × 3 replicas; `audit_brutal.py` gates ON-state contact<8 Å +
   d_repel>2·SEM + control-margin vs the Gln isostere + direction>0) → triple-confirm (`bpair_triple.py`:
   PROPKA Wyman `dpKa`<0 + GFN2-xTB `dE_switch`>0). Ship only designs where all three method classes
   agree in sign, on a target Lys conserved in both orthologs.

**Lesson:** the candidate with the strongest raw MD signal (des_962_5_H37, d_repel 2.23) FAILED triple-confirm
(PROPKA dpKa +1.05 — wrong sign); MD alone would have shipped a backwards switch. The orthogonal classes are
not optional. Fully-confirmed de-novo survivors `des_287_5_H37` + `des_2106_2_H43` entered submission_v7.

## 7. Novelty gate (structural; hard, <3/4 = DQ)

Adaptyv/Proteinbase novelty is per-design and STRUCTURAL: max TM-score of the single-sequence (ESMFold2) fold,
>70% coverage, vs PDB + AFDB. De-novo sequences are always ≤30% id, so the structure axis decides:
`<0.5 → 4/4`, `0.5–0.8 → 3/4`, `≥0.8 → 2/4 (DQ)`. The user-run portal is the ONLY authority; the local
PDB-only foldseek proxy (`tools/foldseek_dl/foldseek` + `runs/novelty_struct/pdbdb` + `TMalign`) is
RELATIVE-rank only — it failed to reproduce a known 2/4 (AFDB half missing). Always portal-check the final
set and remove every <3/4 before locking. Expression: rank by SAP + charged composition (Glu% is the top
Adaptyv-benchmark predictor) + GRAVY; do not weight the ProtParam instability index (Pro false-positive).

## 8. De-novo multi-His pH-switch (Baker 2025 SOTA) — the v8/v9 switch core

Single-His was the WEAK end of the published method (Ahn, Coventry, … Baker, bioRxiv 2025.09.29.678932):
every wet-lab-validated pH-sensitive binder there carried ≥2 His-cation contacts, fold-change scales with
the NUMBER of independent His→target-cation repulsions, and TNF-α itself hit 79× (at pH 5.4). Recipe =
install **≥2–3 SEPARATED His, each repelling a DISTINCT conserved TARGET cation**, spread across the epitope
— NOT clustered in the binder (clustering is anti-cooperative, the old dead-end). Conserved human↔mouse TNF
cations (mature-seq index): 1,6,27,39,60,77,85,93,98,107,123; **K6 and K93** are the strong close anchors.

Two build routes, both validated here:
- **Track A** (`wave3_hisscan_build.py`): point-mutate His COMBINATIONS (2–3 distinct-cation positions) onto
  already-binding balanced groove scaffolds. Fast; reuses the Boltz x-react pipeline.
- **Track B** (`wave3b_prep.py` → `wave3b_mpnn.sh` → `wave3b_analyze.py`): pre-set His at the chosen
  distinct-cation positions, then LigandMPNN (vanilla `proteinmpnn_v_48_020`) REDESIGNS the interface shell
  AROUND the fixed His so the interface packs to accommodate them (backbone unchanged → novelty = parent's).
  Higher yield AND can FIX a wrong switch sign (des_4098_1 failed as a point-mutant, confirmed after
  shell-adaptation). `wave4b.sh` = His-biased (`--bias_AA H:1.5`) deeper run → designs with 4 confirmed His
  on 4 distinct cations (e.g. `des_3050_6_c39`).

Gates (both tracks): fold vs human AND mouse (Boltz) → **dual-strong** = min-ipSAE≥0.55 AND Δ_hm≤0.13
(balanced = genuine cross-react) AND keeps binding vs the `_wt` scaffold AND ≥2 His within 8Å of DISTINCT
target cations in the bound pose → **triple-confirm EACH His** (PROPKA Wyman `dpKa`<0 AND GFN2-xTB
`dE_switch`>0; design confirmed if ≥2 His pass; `wave3_triple.py` / `wave3b_triple.py`). Expression: ESMFold
single-seq MONOMER pLDDT (`esmfold_monomer.py`, AUROC~0.72) + SAP. Assembly: `integrate_v9.py` (switch core +
diverse portal-verified floor, cap ~2/lineage). All de-novo switch lineages were portal-novelty-confirmed ≥3/4.

Honest ceiling: functional (kinetic) release is seconds-scale and NOT observable in unbiased MD (d_iface≈0
is EXPECTED, not a failure); the thermodynamic proton-linkage (PROPKA + xTB, cross-checked by His-cation
contact count) is the honest in-silico selectivity evidence. pH 6.0 (vs Baker's 5.4) is intrinsically harder
(His pKa ~6.5), so maximize the number/placement of His-cation contacts to push the predicted linkage.
