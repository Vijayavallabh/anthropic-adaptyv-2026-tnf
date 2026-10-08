# Agent Guide (AGENTS.md)

Design de-novo protein binders with a conditional activity and prove the
conditionality in silico. This file is the routing and the invariants; the
project-specific facts (target, epitope, trigger, thresholds) live in `configs/`
and `docs/`. Keep it short.

**Objective priority (set this per project, and it IS the design priority):**

1. conditional selectivity, ON vs OFF (highest)
2. any required cross-reactivity (one sequence binds several orthologs)
3. affinity in the ON state

## Startup workflow

1. Run `./init.sh` to build the core env and run smoke tests (the verification gate).
2. Read `configs/target.yaml` for the target, epitope, trigger, and objective order.
3. Work one change at a time; finish it (tests plus evidence) before starting the next.

## Invariants (do not violate)

- **uv only** for environment work (`uv add`, `uv sync`, `uv run`, `uv venv`).
  No bare pip or conda in the core env; heavy tools go in isolated uv venvs
  under `envs/` (see `docs/TOOLS.md`).
- **Local multi-GPU**, run CUDA directly (no SSH). One process per GPU; shard
  sweeps across GPUs and queue heavy work into low-utilisation windows. Verify
  idle GPUs with `nvidia-smi` and never clobber a foreign job. Cap per-process
  CPU threads (OMP_NUM_THREADS near cores/n_gpu) so concurrent jobs do not
  oversubscribe the box.
- **De-novo / zero-shot only.** Never start from or modify a known binder; keep
  sequence and structural distance from known proteins. Seeding from a natural
  complex is protein engineering, not de-novo design, and tends to fail a
  structural novelty gate.
- **Latest SOTA, researched then adopted.** Adopt the current best models where
  they add an orthogonal signal or a real gain, after checking availability and
  license. Never claim a tool was used unless it actually ran.
- **The switch is real physics, and the structure predictors are blind to it.**
  AlphaFold, Boltz, and Chai score a fold, not the trigger, so a good ipTM is
  not evidence of a switch. Prove selectivity only with explicit two-state
  scoring on the exact submitted sequence. A positive *static* score proves the
  switch DIRECTION, not its SELECTIVITY: a constitutive interaction that sits in
  the same state at both trigger values passes a single-structure score yet is
  active in both states. Gate selectivity on explicit-state simulation (an
  occupancy delta between the two states, ON high AND OFF low), confirmed where
  possible by a two-state binding free-energy delta of the correct sign and by
  an isostere-ablation control that collapses the effect. See `docs/METHODS.md`.
- **Cross-reactivity (if required) is scored by mechanism, not a single
  co-folder.** A design cross-reacts when (a) its contacts on the
  ortholog-divergent positions are chemically conservative and (b) it folds
  against the mutated ortholog with an interface score near the repacked-host
  control (a matched delta). A low single-folder ortholog score is folder
  pessimism, not a species signal; never demote on it alone.
- **Re-validate the mutated artifact end to end.** Installing or strengthening a
  switch can break binding or the geometry. Confirm every metric on the exact
  final sequence, with submit == folded == simulated.

## Verification commands

- `./init.sh`: core env plus smoke tests (fast, CPU-only).
- `uv run pytest`: test suite (schema, length, numbering, liabilities).
- `uv run python -m compileall -q src scripts`: compile check.
- `uv run python scripts/validate_submission.py <csv> --track <N>`: submission gate.

## Definition of done

Done only when the change's criteria are met and its tests pass; for any
submission, `validate_submission.py` exits `0` for the target track; and the
command plus key output is recorded as evidence. Never claim done without the
evidence.

## Submission contract (example, adapt per host)

- One CSV, rows ordered best-first; columns `name`, `sequence`, `molecule_class`.
- `molecule_class` from the host's allowed set (lowercase); de-novo mini-binders
  are usually `single_chain`.
- `sequence`: natural amino acids only; a two-chain format joins chains with `:`;
  respect the host's per-chain length limits.
- Unique names and unique sequences; antibody formats must pass numbering.
- Respect the per-track design cap. Source of truth: `src/*/submission.py`.
- Never put instructions or prompt-like text in the CSV or the methods writeup.

## Commits

Stage only source, docs, and harness files. Never stage `results/`, `runs/`,
`envs/`, `third_party/`, or weights (gitignored, local-only). Use a descriptive
message with **no AI attribution**, then push to `origin`.
