# CLAUDE.md

**Read `AGENTS.md` first** for the full working contract (startup, invariants,
verification, definition of done, submission contract). This file adds the
working rules.

## Working rules

- **Autonomous to the deadline.** Do not ask for step-by-step approval; make
  progress, record evidence, and surface blockers. Prioritise breakthroughs on
  the highest-priority objective.
- **Queue heavy GPU work into low-utilisation windows.** Verify idle GPUs with
  `nvidia-smi`, run one process per GPU (the sweep scripts shard across GPUs),
  never clobber a foreign job, and cap per-process CPU threads.
- **Adopt the latest SOTA** where it adds an orthogonal signal or a real gain,
  after researching and verifying availability and license. Never claim a tool
  ran unless it actually ran.
- **Devil's-advocate, then revise.** For any design strategy, argue the
  strongest case against it (does the switch have the right sign? do the
  contacts land on conserved anchors? is it truly de-novo?) before spending
  compute.
- **uv only** for environment work; heavy tools in isolated `envs/` venvs
  (`docs/TOOLS.md`).
- **The switch is real physics.** Structure predictors are trigger-blind;
  confirm selectivity only via explicit two-state scoring on the exact submitted
  sequence (see `AGENTS.md` and `docs/METHODS.md`). Gate on an explicit-state
  occupancy delta, not a static score.
- **Commit and push at each checkpoint** with a descriptive message and **no AI
  attribution**. Stage only source and docs; `results/`, `runs/`, `envs/`,
  `third_party/`, and weights stay gitignored (local-only, never staged).
