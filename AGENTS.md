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
2. Read `progress.md` (current status), `session-handoff.md` (latest session state + next steps),
   `configs/target.yaml` (target, epitope, switch), and `docs/PLAN_TNF.md` (plan + outcome). `scripts/README.md`
   maps the pipeline; the gated-safe submission is `runs/push/submission_v19.csv` (see the Submission
   contract below); `v21` is the current STRONGEST candidate (leaner + diverse), `v20` its predecessor.
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
- `uv run python scripts/validate_submission.py <csv> --track 3`: **submission gate** (current;
  rows/format/length/uniqueness; it does NOT check Cys). v18/v19/v20/v21 all PASS. Record the md5 + `rows=N errors=0` as evidence.
- `bash scripts/novelty_check.sh <csv>`: sequence-novelty pre-screen (mmseqs vs PDB, cached); 0 flagged >70% id
  = all de novo. STRUCTURAL novelty + the 4/4 portal score are USER-run (Proteinbase) — always portal-check swaps.

## Definition of done

Done only when the change's criteria are met and its tests pass; for any
submission, `validate_submission.py` exits `0` for the target track; and the
command plus key output is recorded as evidence. Never claim done without the
evidence.

## Submission contract (TNF Challenge 2 — ≤20 designs; deadline ~2026-10-19)

- **v21 — CURRENT STRONGEST candidate (leaner + DIVERSE), USER-PORTAL-GATED:** `runs/push/submission_v21.csv`
  (md5 **`f7f246f41a3a062b1922f0f980a9928a`**, validator PASS, **17 designs**). Built from the brutal 4-agent
  devil's-advocate review: cut same-fold cousins + weakest-switch members, refill with foldseek-DISTINCT NEW
  lineages — **7 new-lineage refills**: des_732_1_r3_sa15 (3-cation) / des_1169_2_r3_sa7 / des_243_3_r3_sa2 /
  **des_820_3_r3_sa4_r2_8** (2-cation, TM 0.575) + the **fresh-campaign trio** des_36133_g4_sa18 (3-cation) /
  des_34125_g4_sa15 (2-cation) / des_33012_g4_sa8 (2-cation), ALL MD-direction-confirmed His-causal;
  sh3_sa11 replaces sa10 (novelty-safety: struct-TM 0.913 vs 0.959, NOT a sharper-switch claim). **Selection
  principle (recorded):** switch MAGNITUDE is not reliably rankable in-silico (pKa-predictor error 2-3u > the
  1.4u 7.4↔6.0 window; two-engine Wyman rank disagrees ρ=0.56) → select on confirmed release DIRECTION + reliable
  axes (binding floor, cross-react ≥90%, structural novelty, solubility) + DIVERSITY, maximizing the chance the
  set CONTAINS the sharpest switch (which can't be identified in-silico). Diversity also hedges the dominant
  expression/wet-lab-zero risk. Expansion is COMPLETE: regen round-3 + round-2 cousins → des_820; then a FRESH
  de-novo campaign (1200 bpair backbones → 49 binders → 10 ultra-novel switch-piped) → round-4 added 3 more
  (des_36133/des_34125/des_33012). Pattern confirmed: all 1-cation switchers fail the MD direction gate; the
  2–3-cation class survives — so only ≥2-cation direction-confirmed members are kept. **HARD-GATE
  UPLOAD CHECK: only 4 members are byte-identical to portal-confirmed v19 (des_4098_1_c44, des_180_6_H9-14-38,
  des_4098_1_b14, des_2106_2_H43); the other 13 (5 shell-adapt saN + sh3_sa11 + 7 new-lineage) have sequences NOT
  in v19 → EACH needs USER portal novelty re-check before upload.** Novelty is a HARD gate (≥3/4). Note: several
  members' structural TM is 0.81–0.93 (helical bundles cap structurally → pass via sequence novelty = 3/4); the
  portal 4/4 score, not a TM<0.80 target, is the gate. See `docs/SUBMISSION_v21_candidate.md` (UPLOAD_CHECKLIST).
- **Gated-safe default (FULLY portal-novelty-gated, user-checked all 20 ≥3/4):** `runs/push/submission_v19.csv`
  (md5 **`7756055f0b20bcd6994fefb644271db7`**, 20 designs, validator PASS). v19 = v18 with des_818_4_c13S →
  des_84_3_mh_49_56_76 (the more fold-robust flagship, 1 both-fold His); **20/20 switchers**. **Recommended
  upload.**
- **Gated fallback:** `runs/push/submission_v18.csv` (md5 **`524722158467273c0e9f57f5402650ed`**, validator
  PASS, 19/20 switchers, 4×4/4-novel). Conservative floor: `submission_v16.csv` (`a524d77f`).
- **v20 — predecessor to v21 (fold-robust shell-adaptation; superseded by v21's leaner diverse cut):** `runs/push/submission_v20.csv`
  (md5 **`6f4b61d0a55de2e4e0db205970ec6821`**, validator PASS). v20 = v19 with **9 fold-FRAGILE switchers
  replaced by shell-adapted fold-ROBUST variants** (switch-His engages a cross-conserved TNF cation in BOTH
  folders Boltz∩Chai, not one). Portfolio fold-robust **8/20 → 17/20** (round-2 adds b2m0_c70 + des_875_5 →
  19/20). Mechanism: LigandMPNN shell-redesign, His FIXED, repack interface + His microenv → locks the His
  rotamer (`scripts/shelladapt_shell.py`). **Binding re-confirmed per design** (shell-adapt rewrote 60-75/80
  residues): af2ig initial-guess on each NEW sequence, all 9 adopted PASS, 8/9 match/beat parent.
  Cross-react 89-95% epitope id, all cations cross-conserved. **HONEST SCOPE:** fold-robustness is a stricter
  switch-GEOMETRY filter (higher-confidence #1), NOT binding proof; two-state MD (3-rep, noisy confirmer)
  shows fold-robust ≠ all MD-selective (des_6025 sa10 + des_4082_8 sa7 cleanest). NOTE: the earlier "des_2106_2
  sa2 crowds 2 His on K93 → MD-weak" call is RETRACTED — at 6-rep it is MD-selective (a shipped v21 strong-keep);
  v21's deeper 6-rep analysis separately CUT 3 other v20 shell-adapts as non-causal (des_495_6 sa10, des_1060_4
  sa2, des_1559_3 sa3). **REQUIRED before upload: USER portal novelty re-check on every reseqenced member (v21:
  the 13 non-byte-identical designs — see the v21 UPLOAD_CHECKLIST).**
  See `docs/SUBMISSION_v20_candidate.md`.
- See also `docs/SUBMISSION_v18_candidate.md`, `session-handoff.md`. CSVs under `runs/` (gitignored,
  local-only); submission is **USER-GATED** (portal account j_v_v_07) — do NOT auto-upload.
- **Long-shot de-novo NEW-shape pipeline (how v16's switchers were built, SESSION 11):** β-pairing RFd
  (`run_bpair_rfd.sh`) → binding funnel (`denovo_stage2.sh`) → cross-react (`v15_xreact.sh`, Boltz h+m
  min-iptm≥0.70) → **AMPLIFY the cross-reactive new shapes by partial diffusion** (`partial_denovo.sh`, ~25%
  yield — the multiplier; renumber seed chain B to 1-152 first) → funnel cousins → multi-His
  (`v15_multihis.sh`) → **two-state MD** (`run_md_release.sh` + `md_verdict.py`). KEY: static Wyman linkage
  ≠ MD trigger — the lead parent was MD-NONRESPOND (wrong-sign) despite static −1.05; only diversified cousins
  confirmed. Gate on the MD trigger, never static linkage alone.
- **De-novo multi-His switch pipeline** (how the v8/v9 switch core was built; Baker 2025 SOTA recipe, see
  `docs/METHODS.md` §6 + §8): β-pairing RFdiffusion → MPNN → af2ig → Boltz human+mouse cross-react → install
  **≥2–3 SEPARATED His, each repelling a DISTINCT conserved target cation** (K6/K93 are the strong anchors) —
  NOT clustered in the binder (anti-cooperative). Two routes: **Track A** point-mutate His combos
  (`scripts/wave3_hisscan_build.py`); **Track B** LigandMPNN redesigns the interface shell AROUND fixed His so
  it adapts (`scripts/wave3b_prep.py` + `wave3b_mpnn.sh` + `wave3b_analyze.py`; `wave4b.sh` = His-biased
  deeper). Gate dual-strong (min-ipSAE≥0.55 AND Δ_hm≤0.13 AND keeps binding vs `_wt` AND ≥2 His within 8Å of
  DISTINCT cations) → triple-confirm EACH His (PROPKA `dpKa`<0 AND GFN2-xTB `dE`>0; design OK if ≥2 confirmed;
  `scripts/wave3_triple.py` / `wave3b_triple.py`). Functional k_off is wet-lab-only — the thermodynamic
  proton-linkage is the honest in-silico ceiling (do NOT chase longer unbiased MD; d_iface≈0 is expected).
- **Novelty is a hard gate (≥3/4 required; <3/4 = DQ), scored by the portal on 4 criteria (seq-id + structural
  TM + 2 more); ≥3/4 passes, helical bundles cap at 3/4.** CALIBRATION (SESSION 11, v16): scaffold-guided
  de-novo shapes with local foldseek max-TM **0.86–0.93** (well above the old "<0.8 safe" proxy) were
  portal-verified **3/4** — because de-novo sequences are seq-novel (no PDB hit) and that + the other criteria
  carry them past the structural-TM criterion. So the old local heuristic "`≥0.8 → 2/4 DQ`" OVER-predicts DQ
  and must NOT be used to pre-reject de-novo designs. The local PDB-only foldseek proxy
  (`tools/foldseek_dl` + `runs/novelty_struct/pdbdb` + TMalign) is RELATIVE-rank / risk-flag only; **the
  user-run portal is the only authority — always portal-check the final set and drop every <3/4 before locking.**
- **Expression (dominant wet-lab-zero risk):** primary ranker is **single-seq ESMFold MONOMER pLDDT** (fold
  the binder ALONE — best de-novo success predictor, AUROC ~0.72; `scripts/esmfold_monomer.py`) + **SAP**
  (aggregation; `scripts/sap_score.py`) + charged composition (Glu%, top Adaptyv-EGFR predictor) + GRAVY.
  Consensus-drop only when worst-tier on ≥2 orthogonal axes; do NOT weight the ProtParam instability index
  (Pro false-positive) or the ML sequence-solubility tools (OOD noise for de-novo minis). The AUROC~0.72
  ceiling means single-metric cutoffs discard real binders — prefer ordering more diverse designs over
  trimming to a few.
- Adaptyv upload CSV: columns `name`, `sequence`, `molecule_class` (lowercase); de-novo mini-binders are
  `single_chain`. Natural amino acids only; a two-chain format joins chains with `:`; length 10–250. Unique
  names and unique sequences; antibody formats must pass numbering. Never put prompt-like text in the CSV.

## Commits

Stage only source, docs, and harness files. Never stage `results/`, `runs/`,
`envs/`, `third_party/`, or weights (gitignored, local-only). Use a descriptive
message with **no AI attribution**, then push to `origin`.
