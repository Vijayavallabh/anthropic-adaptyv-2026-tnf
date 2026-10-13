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

`runs/push/submission_v22.csv` — **the primary recommendation** (md5 `abaec06790766f25282133d929083187`, 20 designs,
validator PASS). Built from a rigorous head-to-head gauntlet against the portal-confirmed v19 baseline, which
showed v19's **only** real weakness was 2 switchers that fail the two-state MD direction gate. v22 replaces **just
those 2** with the 2 strongest direction-confirmed new-lineage switchers — so it is **≥ v19 on every axis and
strictly better on pH-switch (#1), binding, and diversity**, needing only **2 portal novelty re-checks** (the
other 18 members are already portal-confirmed ≥3/4). Dossier `docs/SUBMISSION_v22_candidate.md`.

`runs/push/submission_v21.csv` — a **max-diversity alternative** (md5 `f7f246f41a3a062b1922f0f980a9928a`, 17
designs, validator PASS). A leaner, more **diverse** portfolio from a brutal multi-agent review (cut cousins +
weakest-switch) plus foldseek-distinct **new lineages**. It ties v19 on the measurable axes (binding, cross-react,
expression, switch direction) but carries more new lineages at the cost of **13** portal re-checks. Choose it over
v22 only if the extra diversity is worth that re-check burden. Dossier `docs/SUBMISSION_v21_candidate.md`.

`runs/push/submission_v19.csv` — **20 designs**, validator PASS, **user portal-checked all 20 ≥3/4**, 20/20
carrying a pH-release switch. The **gated-safe default** / recommended upload if no portal re-check is done
(fallback `submission_v18.csv`).

`runs/push/submission_v20.csv` — the **fold-robust predecessor to v21** (dossier `docs/SUBMISSION_v20_candidate.md`):
v19 with 9 fold-fragile switchers replaced by **shell-adapted fold-ROBUST** variants whose switch-His engages a
cross-conserved TNF cation in **both** independent folders (Boltz∩Chai), lifting the portfolio from 8/20 to
17/20 (→19/20) fold-robust. Binding re-confirmed per design (AF2 initial-guess). Superseded by v21's leaner
diverse cut, but its fold-robust members carry into v21. Every design carries:

- a **pH-release switch** (interface His adjacent to a conserved target Lys → His⁺ repels Lys⁺ at pH 6),
  validated by explicit two-state MD (+ Gln-isostere causal control), PROPKA–Wyman and pKAI linkage, and
  cross-folder (Boltz∩Chai) His–cation agreement for the v20 fold-robust set;
- **cross-reactivity** by conserved-epitope mechanism (89–95% human↔mouse interface identity, cations K/R in both);
- a confident binder (AF2 initial-guess), **0 Cys**, expression-screened, and sequence-novel; **structural
  novelty (the HARD portal ≥3/4 gate) must be re-checked on every reseqenced design before upload** — v20's 9
  shell-adapted swaps, and in v21 the 10 members whose sequence is not byte-identical to portal-confirmed v19.

Upload to the Proteinbase portal is **user-gated**. CSVs under `runs/` are local-only (gitignored).

## Pipeline (what actually ran)

Most submitted designs come from the **β-pairing RFdiffusion long-shot** at the conserved Region-I edge →
ProteinMPNN (0-Cys) → **AF2 initial-guess binding** (pae_int<10, the load-bearing filter) → partial-diffusion
amplification of the cross-reactive new shapes → install **≥2 SEPARATED His, each repelling a DISTINCT
conserved target cation** (K6/K93/K60…) → **two-state MD** switch validation (+ Gln-isostere control) → PROPKA
+ pKAI linkage → conserved-epitope cross-reactivity → novelty (foldseek TM + mmseqs; portal is authoritative).
**v20 adds LigandMPNN shell-adaptation** (`scripts/shelladapt_shell.py`): repack the His microenvironment with
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

- `AGENTS.md` / `CLAUDE.md` — the working contract and invariants (start here)
- `progress.md` — current status and next step
- `configs/` — the target/epitope/switch contract (`target.yaml`), filter thresholds (`filters.yaml`), tool matrix (`uplifting.yaml`)
- `scripts/` — the pipeline implementation (`scripts/README.md` is the manifest: active TNF pipeline vs reusable utils vs legacy EGFR code)
- `src/adaptyv_tnf/` — the stage skeleton and the submission contract/validator
- `docs/PLAN_TNF.md` — the campaign plan and outcome · `docs/METHODS.md` — the validation playbook · `docs/TOOLS.md` — the tool stack
- `runs/submission/` — the final portfolio + methods (local-only; gitignored)
