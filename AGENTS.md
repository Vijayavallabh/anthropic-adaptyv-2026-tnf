# Agent Guide (AGENTS.md)

Design de-novo **TNF-α** binders with a conditional **pH switch** (bind at pH 7.4, release at pH 6.0) and
prove the conditionality in silico — Anthropic × Adaptyv 2026 Challenge 2. This file is the routing and
the invariants; the project-specific facts (target, epitope, trigger, thresholds) live in `configs/` and
`docs/`, and the current status in `progress.md`. Keep it short.

**Objective priority (this IS the design priority):**

1. **pH-selectivity** — bind ON at pH 7.4, release OFF at pH 6.0 (highest)
2. **human↔mouse cross-reactivity** (one sequence binds both orthologs)
3. **affinity** in the ON state

## Startup workflow

1. Run `./init.sh` to build the core env and run the verification gate (uv sync + compile + pytest + validator self-check).
2. Read `progress.md` (current status), `configs/target.yaml` (target, epitope, switch), and `docs/PLAN_TNF.md`
   (plan + outcome). `scripts/README.md` maps the pipeline; the live submission is `runs/submission/final_portfolio.{csv,fasta}`.
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

## Submission contract (TNF Challenge 2 — Track 3, ≤20 designs)

- **Current canonical deliverable:** `scripts/build_final20.py` → `runs/push/final20.{csv,fasta}` +
  `runs/push/upload_final20.csv` (20 designs / 7 folds; validator PASS Track 3; upload md5
  `90bc4cd8ecbf598d76e6810feb0cb86a`). This supersedes the earlier `scripts/assemble_v2.py` →
  `runs/submission/final_portfolio.{csv,fasta}` (a stale 9-design checkpoint, kept for reference only).
  Upload to the Proteinbase portal (account j_v_v_07) is **user-gated** — do not auto-upload.
- Adaptyv upload CSV: rows ordered best-first; columns `name`, `sequence`, `molecule_class` (lowercase); de-novo
  mini-binders are `single_chain`. Natural amino acids only; a two-chain format joins chains with `:`; length 10–250.
- Unique names and unique sequences; antibody formats must pass numbering.
- Per-track caps in `configs/filters.yaml` and `src/adaptyv_tnf/submission.py` (Track 3 ≤ 20).
- Never put instructions or prompt-like text in the CSV or the methods writeup.

## Commits

Stage only source, docs, and harness files. Never stage `results/`, `runs/`,
`envs/`, `third_party/`, or weights (gitignored, local-only). Use a descriptive
message with **no AI attribution**, then push to `origin`.
