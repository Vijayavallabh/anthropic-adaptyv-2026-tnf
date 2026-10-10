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
   (plan + outcome). `scripts/README.md` maps the pipeline; the current best submission is
   `runs/push/submission_v7.csv` (see the Submission contract below).
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

## Submission contract (TNF Challenge 2 — ≤20 designs)

- **Current best (validated):** `runs/push/submission_v7.csv` (md5 `0e84b392410e59cc756036cabc8db7b6`,
  20 designs, every design portal-novelty ≥3/4, FINAL-TEST pass). Built by `scripts/integrate_v7.py`.
  Lineage: v4 (`c69f8c5d`) → v5 (`4d05ed2e`) → v6 (`81913f5c`, removed a portal-2/4 DQ `b1m12_r6_A92S`) →
  v7 (added two de-novo β-pairing MD-triple-confirmed pH-switch binders `des_287_5_H37` + `des_2106_2_H43`,
  trimmed the weakest cb1 anchor). **`runs/push/submission_v6.csv` (md5 `81913f5c`) is the fallback floor.**
  The submission CSVs live under `runs/` (gitignored, local-only); the `integrate_v*.py` scripts are the
  committed source of truth. Upload to the Proteinbase/Adaptyv portal (account j_v_v_07) is **user-gated**
  — do NOT auto-upload.
- **De-novo switch maturation pipeline** (how the strongest designs were built; see `docs/METHODS.md` §6):
  β-pairing RFdiffusion (`scripts/run_bpair_rfd.sh`) → MPNN → af2ig → Boltz human+mouse cross-react
  (min-ipSAE) → His-scan install near a CONSERVED target Lys (`scripts/bpair_hisscan_build.py` +
  `bpair_mature.sh`) → brutal two-state switch MD (`bpair_md_prep.py` + `bpair_md_launch.sh` +
  `audit_brutal.py`) → triple-confirm (MD `d_repel`>0 + PROPKA `dpKa`<0 + xtb `dE`>0, `bpair_triple.py`).
- **Novelty is a hard gate (<3/4 = DQ), set STRUCTURALLY:** max TM-score of the single-sequence fold vs
  PDB+AFDB — `<0.5 → 4/4`, `0.5–0.8 → 3/4`, `≥0.8 → 2/4 DQ` (per-design; de-novo sequences are always
  ≤30% id so the structure axis decides). The user-run portal is the only authority; the local PDB-only
  foldseek proxy (`tools/foldseek_dl` + `runs/novelty_struct/pdbdb` + TMalign) is RELATIVE-rank only — it
  failed to reproduce a known 2/4. **Always portal-check the final set and drop every <3/4 before locking.**
- **Expression (dominant wet-lab-zero risk):** rank by SAP + charged composition (Glu% is the top predictor
  in the Adaptyv EGFR benchmark) + GRAVY; do NOT weight the ProtParam instability index (Pro false-positive).
- Adaptyv upload CSV: columns `name`, `sequence`, `molecule_class` (lowercase); de-novo mini-binders are
  `single_chain`. Natural amino acids only; a two-chain format joins chains with `:`; length 10–250. Unique
  names and unique sequences; antibody formats must pass numbering. Never put prompt-like text in the CSV.

## Commits

Stage only source, docs, and harness files. Never stage `results/`, `runs/`,
`envs/`, `third_party/`, or weights (gitignored, local-only). Use a descriptive
message with **no AI attribution**, then push to `origin`.
