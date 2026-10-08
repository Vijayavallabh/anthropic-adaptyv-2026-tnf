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

`runs/submission/final_portfolio.{csv,fasta}` — **9 designs across 6 distinct folds** (byte-verified
FASTA==CSV), with `runs/submission/METHODS.md` as the write-up. Every design carries:

- a **pH-release switch** (interface His adjacent to a conserved target Lys → His⁺ repels Lys⁺ at pH 6),
  confirmed by **three method-classes on all nine** — explicit two-state MD (+ Gln-isostere causal control),
  PROPKA–Wyman linkage, pKAI linkage — plus a **GFN2-xTB** fourth class on the six single-His designs;
- **cross-reactivity** by matched human/mouse Boltz-2 co-fold and worst-case-over-species **ipSAE > 0.5**;
- a confident binder (AF2 initial-guess + Boltz-2), **0 Cys**, expression-screened, portal-novel (≥3/4; one 4/4).

Upload to the Proteinbase portal is **user-gated**.

## Pipeline (what actually ran)

RFdiffusion3 (all-atom) at the conserved **cb1** aromatic/hydrophobic epitope → ProteinMPNN (0-Cys) →
AF2 initial-guess + Boltz-2 binding → His pH-switch install at a conserved Lys → 4-method switch validation →
ipSAE cross-reactivity → Foldseek TM-clustering for true-shape diversity → `assemble_v2.py` (contact-gated,
one design per shape). `scripts/README.md` maps every stage to its scripts.

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
