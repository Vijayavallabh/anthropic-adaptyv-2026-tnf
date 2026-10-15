# FE pipeline — direct alchemical pH-switch free energy (ΔΔG of His⁰→His⁺, bound vs free)

Measures the *functional-release* magnitude of a pH-switch binder as a thermodynamic double difference:

```
ΔΔG_switch = ΔG_prot(bound complex) − ΔG_prot(free binder)      # >0 ⇒ releases at low pH
```

The double difference cancels the proton self-energy, His intramolecular force-field terms, and the PME
net-charge finite-size artifact (legs built in matched boxes), leaving the physical reorganization — so the
*ranking* can be trustworthy even though absolute binding FE is not. Method notes + results:
`docs/SWITCH_FE_findings.md`. Design notes: `PLAN.md`.

## Environment (important)
Run with **`envs/md2`** (OpenMM 8.4) and the system CUDA-12.3 compiler — the `envs/md` (OpenMM 8.6.1) CUDA
platform is PTX-broken on this driver:

```bash
export OPENMM_CUDA_COMPILER=/usr/local/cuda-12.3/bin/nvcc
PY=envs/md2/bin/python          # has openmm 8.4 + openmmtools + pymbar + biopython
```
~770 ns/day on the ~55k-atom solvated systems (HMR + 4 fs).

## Scripts
- `build_system.py CIF OUT --switch-his 45,53,75 [--free --box-from OUT_bound_box.txt]` — PDBFixer →
  disulfides → His variants (switch-His built HIP) → TIP3P/0.15 M → minimize+equilibrate → saves
  `OUT.{pdb,_system.xml,_state.xml,_box.txt}`. Free leg uses `--box-from` to match the bound box.
- `alchemy_ddg.py --sys OUT --switch-his ... --nwin 11 --equil-ns E --prod-ns P --seed S --out J` — TI by
  finite-difference dU/dλ over λ windows (His⁰→His⁺ via NonbondedForce parameter offsets on a global `lam`);
  writes `J.json` with `dG_kcal ± err_kcal`.
- `combine_ddg.py` — ΔΔG = dG_bound − dG_free per design over seeds + the 2σ-gap ranking rule.
- `run_campaign_ti.sh` — build + run both legs × seeds across GPUs (per-GPU lanes). Env: `GPUS`, `NWIN`,
  `EQ`, `PRODB`, `PRODF`, `SEEDS`.
- `run_screen.sh` — same, driven by a `name<TAB>his` TSV (`TSV=…`), for broader-pool screening.
- `screen_analyze.py` — stage-1 FILTER (multiplicity-aware; flags hits, confirms nothing).
- `stage2_analyze.py` — stage-2 confirmation on FRESH seeds only (split-sample), des_820 gold bar
  (all same-sign, |mean| > 5 kcal, |z| > 3; Bonferroni-safe for the screen size).

## Guardrails (learned)
- Split-sample: a 2-seed pass only *filters*; judge each hit on independent fresh seeds. Screening many
  candidates with the noisy (basin-trapped) bound-leg invites multiple-testing flukes.
- The bound leg is basin-trapped — only *reproducible* calls (all seeds same sign, large |z|) are trusted;
  switch magnitude below ~5 kcal is unrankable.
- Terminal His (position 1 or last residue) trips the net-charge assertion in `alchemy_ddg.py` (internal
  template charges don't sum to +1); such designs are not profiled by this pipeline as-is.
- All outputs go under `runs/fe/` (gitignored).
