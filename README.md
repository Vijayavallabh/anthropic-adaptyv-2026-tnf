# anthropic-adaptyv-2026-tnf

De-novo conditional (pH-switch) binders for **TNF-α** — the Anthropic × Adaptyv 2026 Protein Design
Competition, **Challenge 2**. The designs bind TNF-α at **pH 7.4 and release at pH 6.0** (the
endosomal-recycling direction), with **human↔mouse cross-reactivity**, and the switch is proven in
silico before any wet-lab test.

Built on the reusable conditional-binder toolkit (generalized from the Challenge-1 EGFR project); the
TNF-specific target, epitope, and switch live in `configs/` and `docs/`.

## Objective priority

1. **pH-selectivity** — bind ON at pH 7.4, release OFF at pH 6.0 (highest)
2. **cross-reactivity** — one sequence binds human and mouse TNF-α
3. **affinity** in the ON state

Novelty (distance from known proteins) is a hard gate, not a score — final novelty is the Proteinbase portal.

## The submission

**The deliverable** is `submissions/submission_v27.csv` — **18 designs**, validator PASS,
`molecule_class=single_chain`. The method, the design procedure, and the per-design free-energy evidence are in
`docs/METHODOLOGY_v27.md`, `docs/DESIGN_METHOD_v27.md`, and `docs/SWITCH_FE_EVIDENCE_v27.md`.

**What v27 is.** The **#1-pure** portfolio: **0** designs whose switch runs backwards, **5** with a
free-energy-confirmed pH-release switch, every design sequence- and structure-novel. Switch direction is audited
portfolio-wide by direct alchemical free energy — charge the interface His⁰→His⁺ in the bound complex and in the
free binder, then subtract. Close calls are run to 6 → 12 → 24 independent seeds and judged by the z-test and the
distribution-free sign-test, not a fixed-seed "all-but-one" heuristic (a falsifiable control shows that heuristic
breaks at high seed count). The audit caught two designs whose switch is thermodynamically backwards (they grip
*harder* in acid), `des_84_3` and `des_34125`, and both were cut. The free-energy engine is in `scripts/fe/`; the
full audit is `docs/SWITCH_FE_findings.md`.

The set is **18, not the 20** the track allows: two members had to be removed as confirmed anti-switches, and the
built pool held no novelty-safe, non-cousin releaser to backfill the slots without adding portal-novelty risk, so
the honest #1-pure set is 18.

**What is honest, and what is not.** "Dominate every metric" is *not* in-silico-provable, and this write-up does
not claim it. Switch *magnitude* is thermodynamically unrankable (confirmed by both alchemical FE and
enhanced-sampling metadynamics — the bound-leg slow mode is interface reorganization, not the His coordinate);
binding *magnitude* is wet-lab-decided (Anthropic's own data put in-silico affinity ranking near chance, AP≈0.16);
novelty is gated by the Proteinbase portal. What *is* maximized are the provable axes: switch-*direction* purity
(0 confirmed anti-switch), human↔mouse coverage, lineage diversity, and novelty.

Every design carries:

- a **pH-release switch** — ≥2 separated interface histidines, each repelling a distinct conserved target cation
  (K6/K60/K93…) so His⁺ breaks the grip at pH 6 — validated by explicit two-state MD (with a Gln-isostere causal
  control), PROPKA–Wyman and pKAI linkage, and GFN2-xTB cluster energetics;
- **cross-reactivity** by conserved-epitope mechanism (89–95% human↔mouse interface identity, basic cations
  present in both);
- a confident AF2-initial-guess binder, **0 Cys**, expression-screened, and sequence-novel; structural novelty
  (the hard portal ≥3/4 gate) is re-checked on every design before upload.

Upload to the Proteinbase portal is **user-gated**.

## Pipeline (what actually ran)

Most submitted designs come from the **β-pairing RFdiffusion long-shot** at the conserved Region-I edge →
ProteinMPNN (0-Cys) → **AF2 initial-guess binding** (pae_int<10, the load-bearing filter) → partial-diffusion
amplification of the cross-reactive new shapes → install **≥2 SEPARATED His, each repelling a DISTINCT
conserved target cation** (K6/K93/K60…) → **two-state MD** switch validation (+ Gln-isostere control) → PROPKA
+ pKAI linkage → conserved-epitope cross-reactivity → novelty (foldseek TM + mmseqs; portal is authoritative).
**LigandMPNN shell-adaptation** (`scripts/shelladapt_shell.py`): repack the His microenvironment with
the His FIXED so BOTH independent folders (Boltz∩Chai) agree the His engages its conserved cation (fold-robust),
re-confirming binding per design. `scripts/README.md` maps every stage to its scripts.

## Environment & compute

**uv** only for the core env; heavy design models install into isolated per-tool uv venvs under `envs/`
(see `docs/TOOLS.md`). Compute is local multi-GPU (run CUDA directly, one process per GPU). All scripts
resolve their repo root dynamically — no hardcoded paths.

## Quickstart

```bash
./init.sh                 # uv sync + compile + pytest + validator self-check (the gate)
uv run pytest -q          # tests (schema, length, numbering, liabilities)
uv run python scripts/validate_submission.py <designs.csv> --track 3   # submission gate
```

## Layout

- `submissions/submission_v27.csv` — the final 18-design portfolio (the deliverable)
- `docs/METHODOLOGY_v27.md` · `docs/DESIGN_METHOD_v27.md` · `docs/SWITCH_FE_EVIDENCE_v27.md` — the v27 method, design procedure, and free-energy switch evidence
- `configs/` — the target/epitope/switch contract (`target.yaml`), filter thresholds (`filters.yaml`), tool matrix (`uplifting.yaml`)
- `scripts/` — the pipeline implementation (`scripts/README.md` is the manifest: active TNF pipeline vs reusable utils)
- `src/adaptyv_tnf/` — the stage skeleton and the submission contract/validator
- `docs/PLAN_TNF.md` — the campaign plan and outcome · `docs/METHODS.md` — the validation playbook · `docs/TOOLS.md` — the tool stack
