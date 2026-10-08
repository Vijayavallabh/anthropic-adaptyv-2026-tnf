# conditional-binder-design

A reusable pipeline for designing de-novo protein binders with a *conditional*
activity (bind only under a trigger such as a pH window, a cofactor, or a
specific conformational state) and for proving that conditionality in silico
before any wet-lab test.

Environment: **uv** only. Compute: local multi-GPU, run CUDA directly. Set the
target, epitope, and trigger per project in `configs/target.yaml`.

## What "conditional" means here

The hard, high-value property of a conditional binder is the switch itself:
tight binding in the ON state and weak binding in the OFF state. A good fold or
a high docking score does not prove a switch, because the common structure
predictors are blind to the trigger. Proving selectivity needs explicit
two-state physics (see `docs/METHODS.md`).

## Objective ordering

Rank objectives up front and design to the top one. Default order for a
conditional binder:

1. conditional selectivity, ON vs OFF (highest)
2. any required cross-reactivity (one sequence binds several orthologs)
3. affinity in the ON state

Treat novelty (distance from known proteins) as a gate, not a score, when the
host scores it.

## Quickstart

```bash
./init.sh                 # build the core uv env + run smoke tests (the gate)
uv run pytest             # tests (schema, length, numbering, liabilities)
uv run python scripts/validate_submission.py <designs.csv> --track <N>   # submission gate
bash scripts/novelty_check.sh <designs.csv>                              # novelty gate (fail-loud)
```

Heavy design models install into isolated per-tool uv venvs; see `docs/TOOLS.md`.

## Layout

- `AGENTS.md` / `CLAUDE.md`: working contract and invariants (start here).
- `configs/`: the target/epitope/trigger contract (`target.yaml`, a template),
  filter thresholds (`filters.yaml`), and the tool matrix (`uplifting.yaml`).
- `scripts/`: the pipeline implementation (generation, sequence design,
  folding and binding scoring, switch install plus two-state scoring,
  cross-reactivity, MD, novelty, selection, submission build).
- `src/`: the stage-based skeleton (`stageN_*`) and the submission contract.
- `docs/METHODS.md`: the method and validation playbook.
- `docs/TOOLS.md`: the tool stack and install.

The code under `scripts/` and `src/` is the working implementation carried over
from a prior conditional-binder project. Reuse the patterns and the validation
engine; swap the target-specific parts (the epitope anchors, the switching
residue) for the new project.
