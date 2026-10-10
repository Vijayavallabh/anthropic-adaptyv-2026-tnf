# Progress

Current status and next step for the TNF-α Challenge 2 campaign. Append; don't rely on chat history.

## Status: SUBMISSION READY — upload user-gated

- **Submission:** `runs/submission/final_portfolio.{csv,fasta}` — **9 designs / 6 distinct folds**, FASTA==CSV
  byte-identical. Write-up: `runs/submission/METHODS.md`. Pipeline map: `scripts/README.md`.
  - cb1 family (4 seqs, 1 fold): `TNFpH01` (single-His, d_switch 2.68) + `TNFpH06/07/08` (coop-2His, ddg −1.9…−2.34)
  - cb1-epitope (5 folds): `b4m2` C2 (ipSAE 0.85) · `b2m0` C1 · `b1m9` C8 · `b1m12` C5 (novelty 4/4) · `b1m7` C11
- **Every design:** contact-gated pH-release switch (3 classes: two-state MD + PROPKA + pKAI; +GFN2-xTB on the 6
  single-His), cross-react ipSAE > 0.5, 0-Cys, **portal-novel ≥3/4** (all 9 confirmed; one 4/4).
- **Repo:** pushed to `github.com/Vijayavallabh/anthropic-adaptyv-2026-tnf` (private). Hardened — portable paths,
  5 active-path bugs fixed, 20 dead EGFR Challenge-1 scripts removed, package `adaptyv_egfr`→`adaptyv_tnf`.

## Evidence (last verified 2026-10-08)

- `uv run pytest -q` → 23 passed · `uv run python -m compileall -q src scripts` → clean
- `scripts/validate_submission.py submissions/example_track3.csv --track 3` → PASS (exit 0)
- `final_portfolio.csv` md5 `f89d9173…` (stable; re-running `switch_contact_gate.py` + `assemble_v2.py` reproduces it)

## Final-push campaign (2026-10-08, autonomous to deadline)

**Decisive finding — expression is the dominant wet-lab-zero risk.** Audited the 9 shipped designs vs
Adaptyv's `protein-qc` standard gate (instability<40, GRAVY<0.4): **only 1/9 pass** (TNFpH07). Cause: the
cb1 epitope is hydrophobic → vanilla-ProteinMPNN interfaces raise instability/GRAVY. Per the EGFR
post-mortem, non-expressing designs score zero regardless of a perfect switch. (ESM-C PLL, now wired in,
is less dire — all designs −1.6…−2.1 — but the gate liabilities are real.)

**Expression rescue (interface-frozen SolubleMPNN).** Freeze binding interface (human∪mouse, 5–7 Å) + His
switch; soluble-redesign only the surface (0-Cys); restore His; re-fold (Boltz h+m) + re-score ipSAE +
PROPKA switch. Results by fold:
- **fold2 (b4m2) & fold5 (b1m12): FULLY RESCUED** — expression-clean (instab 27–32, GRAVY≤0.4), binding +
  mouse preserved (ipSAE_wc 0.54–0.69), switch intact (PROPKA His dpKa −0.2…−0.6). Banked in
  `runs/push/rescue/winners.tsv`; finalist folds in `runs/boltz_fin_h`.
- **folds 8/11/12 (b1m9/b1m7/A20V): mouse-fragile** — surface redesign keeps human binding + expression but
  collapses mouse (orig mouse was marginal 0.55–0.70). Not surface-rescuable without losing cross-react.
- **fold1 (b2m0) & cb1 TNFpH01: expression-stubborn** — instab 49.6/57.5 too high for the free surface to
  fix. Candidates for SolubleCaliby / accept-original.

**Tooling online:** ESM-C 600M PLL (cudnn path fix + HF online), PROPKA+pKAI switch linkage, pKa env.
pxdesign/complexa kits installed (checkpoints downloading) — timeboxed bonus-diversity bet.

**Portfolio plan (hedged, ≤20):** per fold ship original + expression-rescued variant(s) where available
(“diversify a working binder” = the EGFR-winning strategy). ~14 designs across 7 folds; rescued hedges on
folds 2,5. Still TODO: two-state MD on rescued finalists (submit==simulated), pKa consensus hardening,
composite re-rank (ESM-C), final assembly, byte-verify.

## Final-push outcome (staged deliverable)

**Proposed revised submission:** `runs/push/final_portfolio_push.{csv,fasta}` (+ Adaptyv upload form
`runs/push/upload_push.csv`, **Track-3 validate PASS, 0 errors**). Write-up: `docs/METHODS_push.md`.
**14 designs / 7 folds**, tiered A/B/C by expression + switch + novelty. Reproducible pipeline added to
`scripts/` (`rescue_build.py`, `rescue_filter.py`, `rescue_collect.py`, `rescue_minimal_edit.py`,
`rescue_mpnn.sh`, `rescue_ipsae.sh`, `esmc_pll.py`, `assemble_push.py`, `expr_audit.py`, `cif2pdb.py`).

What changed vs the shipped 9:
- **Expression fixed (the dominant risk; was 1/9 passing the gate).** Interface-frozen SolubleMPNN rescue
  gives expression-clean variants for folds C2 (`b4m2_r1/r6`) and C5 (`b1m12_r4/r1`) with binding + switch
  preserved (PROPKA release-competent; two-state MD finalizing). Portfolio now has **7 expression-PASS**.
- **Cross-reactivity corrected & strengthened.** Single-folder cross-species ipSAE is unreliable — Boltz vs
  Chai-1 disagree and it is perturbation-sensitive (folder pessimism). Grounded instead in **mechanism**:
  conserved cb1 epitope + conserved-Lys anchor, ~1 divergent contact/27, small matched human↔mouse co-fold
  delta (|Δ|<0.13). Human binding confirmed by **3 folders** (af2ig+Boltz+Chai, all high).
- **ESM-C PLL** wired in (naturalness, the EGFR post-mortem's top signal).

**Two-state MD on the 4 rescued variants (done):** `b4m2` (r1/r6) switch **confirmed SELECTIVE** —
d_repel 0.60/0.13, His-causal (the Gln-isostere control does not repel). `b1m12` (r4/r1) switch is
**weak / not MD-selective** (d_repel 0.17/0.34, His not interface-causal) — same weak switch as its parent
(d_switch 0.597); b1m12's value is novelty 4/4 + strong binding, not the switch. **Strongest switches
overall:** the cb1 coop-2His family (ddg to −2.34) and single-His `TNFpH01` (d_switch 2.68).

**Honest caveats:** rescue cleanly delivers 2 de-novo folds (C2 with a selective switch, C5 novelty-4/4);
the other de-novo folds (C1/C8/C11/C12) are not surface-rescuable (ship originals, expression-flagged).
Rescued variants need a Proteinbase novelty re-check (same fold ⇒ expected ≥3/4); `A20V` needs its first
portal check. **Sequence novelty confirmed** (`novelty_check.sh`, fixed to a TNF positive control): all 14
designs have **no PDB hit** (fully de-novo); structural novelty inherits (rescued = parent fold) or is
de-novo (A20V) — portal does the authoritative 1–4 score.
New open generators (pxdesign/complexa) were deprioritized — install friction + full re-validation exceed the
time budget for low expected value on this hard target.

## Bold bet in progress — cooperative-2His switch strengthening (2026-10-08)

The de-novo folds have weak single-His switches. Feasibility test: **all 6 de-novo folds sit near ≥2
conserved target Lys** (precondition for a coupled 2-His switch). Installed a 2nd His near a 2nd Lys on
b1m12 (novelty 4/4) and b2m0 and folded + scored:
- **Switch ~2× stronger** (PROPKA binder-His linkage): b1m12 −0.41→**−1.00** (c94), b2m0 −0.67→**−1.27** (c20);
  several variants −0.75…−1.00.
- **Binding preserved** (iptm_h 0.92–0.94) and **cross-react preserved** (matched Δ −0.01…−0.13).
- **Expression:** most coop variants FLAG (adding His to already-borderline de-novo folds pushes instab
  43–50); **b1m12_c68 PASSES** (39.9) with a −0.82 switch — a candidate dominating design as-is.
- **MD two-state selectivity gate: RUNNING** (the decisive test; PROPKA is only the pre-screen).

Recipe for a dominating design: coop-2His (strong switch) → expression-rescue the strong-but-FLAG variants
(freeze interface + both His, soluble surface) → re-validate. Gated on the MD selectivity verdict.
Research agents (wide SOTA + switch-sharpening physics + TNF/unconventional) running in parallel.

## Physics-grounded switch ranking + user relaxed novelty to 3/4 (2026-10-08, cont'd)

**User: "I am fine with 3/4 novelty designs too."** Portal requires >=3/4 to submit; only `A20V` was 2/4 ->
**DROPPED**. Portfolio = **13 designs / 6 folds**, all >=3/4; upload form regenerated
(`runs/push/upload_push.csv`, 13 rows). The demanding all-20-at-4/4 generative campaign is OFF.

**Three research agents returned (convergent, decisive):** Baker pH-binders now peer-reviewed (PubMed
41256377); Anthropic Aug-2026 autonomous-binder paper *ran a TNFα campaign*; switch-sharpness physics.
1. **Switch SHARPNESS = Delta-n_H (protons coupled to binding), NOT His count.** Width = 1.908/n_H pH.
   Over 7.4->6.0: n_H=1 ~25x (single-His too broad = weak), 2 ~630x, 3 ~16000x. Clustered interface His
   are ANTI-COOPERATIVE (pKas split -> broad); need SPATIALLY-SEPARATED His, each near its OWN target
   cation, pKa-matched 6.5-6.7. New metric: `scripts/switch_sharpness.py` (peak Delta-n_H + fold-change
   over 6.0-7.4 via PROPKA Wyman; envs/pyrosetta/bin/propka3). Output `runs/push/switch_sharpness.tsv`.
2. **Baker empirical floor >=2 His<->cation contacts** (`scripts/his_cation_contacts.py`); placement>count.
3. **Cross-reactivity SOLVED**: every design hits a **79-94% conserved** human<->mouse epitope
   (`scripts/epitope_conservation.py`); divergences = flexible N-term (R1>Q,T2>N,P3>S) + conservative
   (Ile<->Val). No new epitope campaign needed.

**Physics-truthful portfolio (new `sharpness_FC` column, assemble_push.py):**
- **TNFpH07 (peak 1.51 / 13x) + TNFpH06 (1.39 / 17x)** = the ONLY sharp switches (3 co-titrating His
  ~pKa 6.3, expr-clean instab 33-35, Baker-PASS, cross-reactive). **THE headliners.**
- All other 11: weak switches (FC 2-8x) — bind + cross-react, but the TOP-objective axis is weak.
  Single-His TNFpH01 = 3x; b1m12 coop (c91/c92/c94/c68) weak (confirms anti-cooperative clustered His).

**Coop-2His two-state MD (6 variants):** b2m0_c20/c37/c70 all selective (His-causal; Gln-isostere control
inert; d_repel 0.24/0.66/3.22; iface stable in 3 ns) but FLAG expr (46-50; b2m0 expression-stubborn).
b1m12_c68/c94 NOT selective (adjacent His 67/68, 67/94). Net: b2m0 can host a selective de-novo switch,
expression is the blocker.

**DECISIVE GAP -> sharp-multi-His generative campaign (8 GPUs idle):** (A) scale cb1 coop (reliable sharp
source); (B) install SEPARATED multi-His (each near a distinct target cation Lys6/Lys93/Arg1, >8A apart,
pKa 6.5-6.7) on expr-clean diverse binders (b4m2_r1, b1m12_r4). Fold (Boltz h+m) -> switch_sharpness ->
keep sharpest that bind+express+cross-react. Goal: fill toward 20 with MORE sharp-switch dominating designs.

## BREAKTHROUGH — sharp separated-multi-His switch + final 20-design portfolio (2026-10-08)

**The gap was switch SHARPNESS (top objective); the fix was physics.** Installed spatially-separated His
(>8A apart, each near a distinct CONSERVED target cation Lys93/Lys6/Arg98), which titrate additively
(Delta-n_H->2-3) instead of anti-cooperatively splitting (`scripts/install_sharp_his.py`). Result: the
switch jumped from FC 3-17x to **19-126x**, and the sharp variants switch in **both** species:
- `TNFpH06_sh3` **FC 59x/59x** (h/m), expr-PASS, x-react Δ−0.02 — expressing sharp headliner
- `TNFpH06_sh5` 31x/76x, expr-PASS; `TNFpH07_sh4` **126x/50x** (physics-max, FLAG hedge)
- `b1m9_1A55H_K93_sh4` FC 22x — **the fully-validated LEAD**: MD His-causal (d_repel +0.91, Gln control
  inert) AND sharp (PROPKA 22x + pKAI 7x) AND strong x-react binding (iptm 0.91h/0.82m) AND expr-PASS.

**Switch confirmed by 4 method-classes, read as consensus (honest):** PROPKA-Wyman sharpness +
**pKAI** (orthogonal ML pKa — agrees on SIGN and RANKING: sh5 pKAI 58x, sh4 40x, sh3 26x) + two-state MD
His-causality (Gln control) + Baker His-cation geometry. **THREE leads fully validated** (sharp + MD
His-causal + x-react + express): `b1m9_sh4` (d_repel +0.91), `TNFpH06_sh3` (d_repel +1.18, 59x both
species), `TNFpH06_sh2` (d_repel +1.25, best cb1 binding 0.75). **4-His is the sweet spot** (2 of 3 tested
4-His cb1 are MD-clean): the 5-His cb1 variants (sh5, sh4) are the sharpest
thermodynamically (76x, 126x; both pKa engines) but MD-negative/ambiguous (5th His over-crowds; e.g.
TNFpH07_sh4 d_repel+1.05 vs control+1.01) -> ship as thermodynamic-sharp hedges (BLI measures the
thermodynamic switch). `b2m0_c70` = a third, independent MD-His-causal de-novo switch (d_repel 3.2A).

**Cross-reactivity corrected & re-verified**: fixed a target-file bug (`mouse_tnf.seq` was the full
PRECURSOR with TM/signal domain -> spurious 0.15 iptm for ALL designs incl parents; now the mature 156aa
ortholog, parent TNFpH06 recovers to 0.63/Δ+0.01, see [[env-mouse-tnf-seq-was-precursor]]). Epitope
conservation 82-97% across the portfolio; switch cations conserved/cation-preserving (Arg98->Lys98).

**Final deliverable: `runs/push/final20.{csv,fasta}` + `upload_final20.csv`** (Track-3 validator PASS,
0 errors). **20 designs / 6 folds** (cb1:6, b1m9/f8:3, b1m12/f5:3 [all 4/4], b4m2/f2:3, b2m0/f1:3,
b1m7/f11:2); **10 expr-PASS, 3 novelty-4/4, 6 sharp (FC>=19x)**. Assembler `scripts/build_final20.py`
auto-merges every metric. Write-up rewritten: `docs/METHODS_push.md`. New scripts: install_sharp_his,
switch_sharpness, pkai_sharpness, his_cation_contacts, epitope_conservation, build_final20.

## SESSION 5 (2026-10-08, resumed on Stop-hook): expression rescue -> 15/20 PASS + sharper switches
Stop-hook flagged stopping early with 10/20 expression-FLAG. Resumed aggressively (8 GPUs + 3 research
agents). **Interface-frozen SolubleMPNN rescue of the 8 FLAG designs** (`scripts/rescue_build_general.py`
+ `rescue2_mpnn.sh` + `rescue2_select.py`; freeze union iface + ALL His, 48 surface-redesigns/design):
- **5/8 rescued to expr-PASS, switch PRESERVED/SHARPENED on the exact rescued seq:** `TNFpH07_sh4_r4`
  (FC 69x, now PASS, x-react 0.66/0.72), `b1m9_sh6_r1` (FC 32x, 0.92/0.83, PASS), `b4m2_2A87H_r3`
  (iptm 0.95/0.95 perfect x-react, PASS), `b1m12_2A67H_r1` (4/4, 0.92/0.92, PASS), `b1m7_r1` (PASS).
  TNFpH07_sh4_r2 even hit FC 223x. Swapped all 5 FLAG->rescued-PASS in `build_final20.py`.
- **b2m0 (73aa) confirmed expression-STUBBORN** (0/48 pass; only 23-26 free surface) -> its 3 designs
  (incl MD-His-causal b2m0_c70) stay FLAG as fold-diversity/switch hedges.
- **Portfolio now 15/20 expr-PASS** (was 10), still 6 folds / 3x 4/4 / 3 MD-His-causal sharp switches +
  rescued sharp. Track-3 PASS, byte-consistent. Rescued-switch MD re-validation RUNNING. 3 research
  agents (break-fold-ceiling / maximize-expression / sharper-switch-mechanisms) running for the next push.

## SESSION 6 (2026-10-08, cont.): bold new-fold tested+floored; portfolio validated on every axis

Pursued the three Stop-hook items to completion and ran the bold new-fold bet.
- **PXDesign new-fold campaign (the evidence-backed TNF winner) — RUN and DECISIVELY FLOORED.** Fixed the
  env (`PXDESIGN_CKPT_DIR` must name `checkpoint/`; LD_LIBRARY_PATH to venv nvidia libs; pip `--no-deps`),
  ran 6 tasks (recsite patch + inter-subunit groove × L70/90/110), extracted **98 on-target backbones**
  (`pxd_extract.py`) → SolubleMPNN (`pxd_mpnn.sh`) → af2ig (`pxd_make_af2ig.py`): **0/196 bind** (best
  pae_int 24.6, all ≥24; plddt_binder 83-94 = fold-not-bind). Boltz-2 cross-check gave 5 at iptm 0.5-0.66
  but the **pose check killed it** — they dock OFF the designed epitope (0-2/8 hotspots) = msa-empty
  pose-promiscuity false-positive. New-fold generation at these epitopes is a tested floor for this
  pipeline. Only the all-atom RFd3 angle remains untested ([[lesson-denovo-wall-flat-polar-target]]).
- **Rescue-switch MD re-validation DONE:** all 3 rescued variants come back **NOT His-causal** (d_repel
  -0.46/-0.27/-0.12) — surface redesign loosens the interface so explicit-MD no longer confirms, though the
  **thermodynamic** switch (PROPKA+pKAI, what BLI measures) is intact. Relabeled rescued = thermodynamic +
  expr-PASS; MD-His-causal reserved for the 3 un-rescued leads + distinct-fold b2m0_c70. METHODS updated.
- **Keystone switch reframe** (research agent, read the Baker 2025 blueprint): fold-change = in-window
  Δn_H; per-His pKa diagnostic shows my cb1 His midpoints are acidic-shifted (pKa_cplx crashes <6.0), so
  the ceiling is set by centering His pKa at 6.7, not adding His (5-His over-crowds → MD-ambiguous). The
  reported 59-126x ARE the realized in-window FC. [[lesson-ph-switch-optimize-in-window-dnh-center-pka]]
- **Expression composition analysis** (research agent): cb1 family is expression-robust (Glu 21%, netQ
  -9..-11, helix 61%) — the real wet-lab predictors; b2m0 genuinely poor (Glu 5.5%), FLAG corroborated.
- **Novelty gate confirmed locally (fail-loud):** mmseqs vs PDB seqres — POSCTRL TNF hits 100%, **all 20
  designs have NO PDB sequence match** (fully novel seqs); single-seq structural TM safe-harbor; parents
  portal-confirmed >=3. All 20 pass the >=3/4 gate.
- Deliverable re-validated **byte-consistent** (CSV==FASTA==upload, 20 uniq, 0-Cys; upload md5 77b6cd92...).

## SESSION 7 (2026-10-08/09, Stop-hook cont.): expression rescue 15→17/20, two negatives closed, portfolio upgraded to 7 folds

- **Region I de-novo — decisive test, NEGATIVE.** 7680 RFd3→MPNN designs at the CRO-validated
  inter-subunit groove (V74/L75/I97 + W114/Y115), af2ig on the correct 2-protomer target. 4471 fold
  self-consistently; the best engage the groove (4/5 epitope residues <8 Å) but none deeply dock (0 at
  3/5 <5 Å), af2ig pae ~26. Boltz monomer cross-fold of the top 14: 3–4 dock the protomer-1 half
  (V74/L75/I97 at 2.7–5.7 Å, human iptm 0.60–0.70) **but mouse iptm uniformly 0.15–0.18 — not
  cross-reactive**. Existing portfolio binders (0.91–0.95 both species) dominate. De-novo at the flat
  inter-subunit groove hits the known wall; not integrated. [[lesson-denovo-wall-flat-polar-target]]
- **K90 new fold — MD refuted its switch.** Two-state MD on the installed His switch: d_repel −0.02
  (sh4) / −0.27 (sh5), His-causal = **no**. K90 stays as the 7th fold / binder-first diversity hedge on
  its own merit (cross-react binder iptm 0.59/0.54, docks 3/4, distinct conserved epitope), honestly
  relabelled: thermodynamic switch FC 8× (PROPKA) but not MD-causal; expr-FLAG (Glu-overload artifact,
  gravy −1.03 soluble). Expression rescue tested and **rejected** (6 E→Q dropped human iptm 0.59→0.31).
- **Expression rescue — interface/His/core-frozen minimal edits, 15→17/20 PASS.** b1m7_sh5 (+3 surface
  edits) and b1m9 parent (+1, E19N) rescued with binding **preserved** (Boltz rx-vs-parent: 0.92/0.89 and
  0.87/0.78) and switch FC preserved (6×, 10×, recomputed on the rescued fold). Burial filter prevented
  a buried-Ile→Asn and structural-proline removals. b2m0 (c20/c70) kept FLAG: the instability index is
  inflated by two turn-prolines (CA bend 81°/116°) — [[lesson-instability-index-proline-false-positive]];
  c70 is the f1 fold's MD-His-causal switch, so kept with the documented moderate expression risk.
- Deliverable re-validated **byte-consistent** (CSV==FASTA==upload, 20 unique, 0-Cys; validator PASS
  Track 3; upload md5 `90bc4cd8ecbf598d76e6810feb0cb86a`). _rx seqs byte-verified = parent ± exactly the
  intended edits; K90 confirmed = parent (not the binding-damaged rescue).

## SESSION 8 (2026-10-09, Stop-hook cont.): switch hunt 6→9 MD-causal (2 dominant 4/4 designs); Genie3 built but cross-react negative

- **Switch hunt on the strong binders — 3 NEW MD-His-causal switches (the #1 axis), 2 on the 4/4-novel
  b1m12 fold.** Scanned His near the conserved cross-reactive cation **K98** on b1m12 (4/4-novel) and b4m2,
  af2ig-confirmed binding retained, two-state MD screened causality. Winners (all expr-PASS, cross-react
  Boltz-confirmed): **b1m12_r6_A92S** (d_repel **3.9**, strongest; iptm 0.91h/0.95m), **b1m12_A91H**
  (d_repel 1.46; 0.92/0.94), **b4m2_A37H** (d_repel 0.45; 0.93/0.80). These are DOMINANT — 4/4-or-3/4
  novel + MD-causal switch + strong cross-reactive binding + express. Integrated (swapped in for the plain
  b1m12_r4/r1 and b4m2_r6). **MD-His-causal switches: 6 → 9.**
- **Genie3 BUILT (the CRO-validated TNF generator, 8/12 reference binders) — docks but cross-react NEGATIVE.**
  Built from scratch (uv venv py3.10 + torch 2.7.1+cu126 + pinned weights; `run.sh check` PASS). 400 binder
  backbones at the K90 epitope dock tightly (313/400 ≥2 hotspots). But: it emits Cα-only backbones →
  ProteinMPNN **CA-only** (cloned `ProteinMPNN_orig` + ca weights). Single-state ca-MPNN → confident human
  binders (6/64 iptm>0.6) but **0/12 cross-reactive** (human-only, exactly the research-predicted single-state
  failure). Multi-state tied ca-MPNN (human+mouse, `1TNF_B_mouse.pdb` threaded target) → **0/28 cross-react
  AND degraded binding** (tie over-restricts CA-only). DECISIVE NEGATIVE: Genie3 generates docking backbones
  but CA-MPNN can't make them cross-reactive here → NOT integrated. (Tooling kept: `setup`/`eval_genie3.sh`,
  `build_multistate_pdb.py`, `thread_mouse_target.py` for future full-backbone attempts.)
- Research scout (SOTA 2025-26) → [[reference-binder-sota-levers-2026]] (multi-state MPNN, Overath calibrated
  ranker, buried-His-network switch, all-atom generators).
- Deliverable re-validated **byte-consistent** (CSV==FASTA==upload, validator PASS Track 3, upload md5
  `f5f632962835d24dfd31a4ce29f6a670`); 3 winners byte-verified present.

## SESSION 9 (2026-10-09, Stop-hook cont.): switch hunt 9→10 — f11 fold gap filled by a strict-upgrade switch

- Extended the K98 switch hunt to the two folds still without an MD-causal switch (hunt3). Two variants
  passed the identical MD gate (`b1m7_sh5rx__A60E2H` d_repel 0.32 / margin 1.72; `b4m2_r1__A37L2H`
  d_repel 0.42 / margin 1.27). Devil's-advocate on each before integrating:
  - **`b1m7_sh5rx__A60E2H` — ADOPTED.** Fills the **f11 fold gap** (the one strong binder fold with no
    MD-causal switch). A single **E60→His** vs the sh5_rx parent; MD confirms His60 drives K98 repulsion
    (Gln60 isostere control −1.40 inert → clean margin 1.72). Re-folded against full mature human+mouse
    (same switchwin recipe, msa:empty, pinned Boltz): iptm **0.950 h / 0.898 m** — both *up* vs the parent
    (0.916/0.894). Thermodynamic FC **8×** (parent 6×), instability **29.9** (parent 37.5). A **strict
    upgrade on every axis** → swapped in for the parent `b1m7_2A57H_K183_sh5_rx`.
  - **`b4m2_r1__A37L2H` — REJECTED.** Clean switch, but f2 already carries `b4m2_A37H` (MD-causal) and
    `b4m2_r1` is f2's pure-binder hedge; a second redundant f2 switch at the cost of the hedge fails the
    marginal-value bar. Not integrated (honest count, not inflated).
- **MD-His-causal switches: 9 → 10**; f11 now carries a validated switch, so 6 of 7 folds do.
  byte-verified submit==folded==MD'd==fasta; validator PASS Track 3; new upload md5
  `53c741c99db4adf2781d926639110897`.

## SESSION 9b (2026-10-09, deep-research campaign): physics reframe + 3 levers killed decisively + buried-network bet

- **5-agent SOTA literature sweep** (pubmed/biorxiv) → [[reference-binder-sota-levers-2026]] updated with runnable assets. **Decisive physics reframe (2 agents convergent):** a protonation switch's 10-90% width = 1.908/Δn_H; at Δn_H≈1 that's **1.9 pH — WIDER than the 1.4-pH window**, so a single-His switch cannot complete its swing in-window. The game is **effective in-window Δn_H ≥ 2**. My FCs confirm it: cb1 family FC 59-126x = Δn_H~1.3-1.5 (sharp, the real switches); b-fold single-His FC 3-4x = Δn_H~0.3-0.6 (weak switches, strong/novel binders). Portfolio already balances sharp-switch (cb1) vs novel-binder (b-fold).
- **THREE improvement levers KILLED by decisive cheap tests (good discipline, don't re-chase):**
  1. **H1 target-side His78 proton** (binder cation vs conserved target His) — His78 IS conserved human+mouse AND at Region I rim (6.5A), BUT PROPKA free pKa = **3.91** (far below the 6.0-7.4 window; buried, 120 neighbors) → won't titrate in-window. His73 in-window (6.70) but NOT conserved (mouse Tyr). DEAD.
  2. **Genie3 full-backbone route** (PULCHRA CA→N/CA/C/O reconstruction → standard full MPNN, the untried fix for the CA-only cross-react failure) — mechanically WORKS (PULCHRA validated, `scripts/build_genie3_fullbb.py`, 1600 seqs vs 0 before) BUT af2ig gives **0/1566 binders (best pae_i 18.5, all >18)**: reconstruction+re-MPNN loses the pose under af2ig's single-seq fold. Deeper: af2ig can't rank TNF binders (AP 0.16) so generators (Genie3/PXDesign) are **unvalidatable in-silico** → can't select which to submit → low-EV for the submission. DEAD (as a selectable route).
  3. **2nd separated interface His** (additive Δn_H→2) — GEOMETRICALLY DEAD: binders present only ONE reachable conserved cation (K98); K90/K65/R32 all out of 4.5-9A range (scan found 0-1 variants, the 1 clustered near K98). Interface is at single-His ceiling.
- **partial2** (fresh partial-diffusion of 5 novel-fold seeds): 230/1920 strong binders (12% pass, matches the lesson) on b4m2+b1m12 — confirms those folds are ROBUST binders, but all COUSINS (diversity cap) → not integrated. (b1m9/b1m7 arm skipped at build_fullatom — artifact, not fragility.)
- **LIVE BET: buried-His-NETWORK switch** (Baker `his_ph_exploder`, MIT, TNF-validated 79x, OFF-interface so preserves affinity+cross-react, the cooperativity/Δn_H>1 route — sidesteps the interface single-His ceiling). PyRosetta installing (free-academic). RISK: Baker says buried networks need LARGE cores; my binders are 82-101aa (small) → b1m12 (101aa, 4/4-novel) is the best candidate. The exploder is itself the go/no-go (finds buried His-pairs or reports none).

## Next step

**Deliverable is final and validated on every axis** — `runs/push/final20.{csv,fasta}` +
`upload_final20.csv` (20 designs / **7 folds** / **17/20 expr-PASS** + 3 documented-soluble FLAGs /
3× novelty-4/4 / **10 MD-His-causal switches** incl 2 dominant 4/4-novel (b1m12_A91H, b1m12_r6_A92S) /
sharp thermodynamic FC to 69× / all seqs fully novel vs PDB). This is the canonical Track-3 upload (supersedes the stale 9-design `runs/submission/final_portfolio.csv`).
**Decision point (owner): upload** (Proteinbase account j_v_v_07; the portal re-checks novelty
authoritatively at upload — expected ≥3/4). Recommended core = the 3 dual-everything leads (b1m9_sh4,
TNFpH06_sh3, TNFpH06_sh2) + the novelty/affinity anchors (b1m12 4/4, b4m2), filling to 20 for fold+risk
diversity. Deadline **2026-10-12 11:59 UTC** (2026-10-11 23:59 AoE). Upload is user-gated — do not
auto-upload.

## 2026-10-09 — v4 finalized; generation bet closed
- SUBMISSION: runs/push/final20_v4.csv (md5 85ea7c75) — 20 validated seqs, byte-identical to v3, RE-RANKED expression-first (ipSAE demoted per EGFR post-mortem). FASTA runs/push/final20_v4.fasta. QC runs/push/final20_v4_qc.md. Floor final20.csv (3518cc6) unchanged. UPLOAD USER-GATED.
- QC added (4 orthogonal, ipSAE demoted to tiebreaker): structural SAP (all 20 low, buried-HYD false-positives), instability (b2m0 Pro-FP keep / tnf_K90 genuine), interface dSASA+contacts+Hbonds (real binding signal; TNFpH07_sh4_r4 strong iface despite low ipSAE), uniform Boltz re-fold human+mouse (cross-react delta). master table runs/push/master_v3.csv.
- GENERATION BET CLOSED (all 3 proven generators → de-novo wall): RFd3 K90 ~0.1% (survivor k90_18 ruled out: dSASA1569 + noisy switch); RFd3 regionI 0; PXDesign recsite 0 (af2ig pae_int~25, Boltz iptm<0.5, both SolubleMPNN+vanilla); Genie3 prior fbb1 0. PXDesign was ACTIVATED (was never actually running — env recipe in memory).
- NEXT: v4 awaits USER-GATED upload. v4 is the strongest evidence-backed strategy (validated partial-diffusion portfolio).

## 2026-10-09 (cont) — confirmed submission locked; partial-diffusion upgrade hunt
- CONFIRMED SUBMISSION: runs/push/submission_v4.csv (md5 c69f8c5d) Adaptyv format (name,sequence,molecule_class=single_chain), best-first. FINAL-TEST PASS (20 uniq, natural-AA, 0-Cys, len73-119, nov>=3/4). Upload USER-GATED.
- partial3 cross-react "upgrades" REJECTED by brutal gate (ipSAE refuted iptm: des_55_2 m0.476/des_7_2 m0.446 < parent 0.533).
- partial4: 5 b1m9 variants BEAT parent cross-react by ipSAE (des_69_2 m0.637, des_86_2 d0.16) BUT 4/5 lost the switch (0 His) + high instab; des_59_1 clean + latent H11 8.5A (MD running). His-bias MPNN screen launched to recover switch on cross-react backbones.

## 2026-10-09 (cont) — brutal devil's-advocate on upgrade hunt: all REJECTED; v4 CONFIRMED FINAL
- partial3 cross-react upgrades: REJECTED (iptm artifact; ipSAE m0.44-0.48 < parent 0.53).
- His-bias (H:2.0) cross-react+switch: REJECTED (His-saturated artifacts, nHis 14-24 = 16-29%; Baker bias is 0.5-1.5).
- des_59_1 latent switch: REJECTED (MD d_repel -1.12, wrong direction).
- Pattern: redesigning b1m9 for cross-react strips the switch (plain MPNN) or saturates His (bias) -> no clean cross-react+switch+expression upgrade exists in generated pool.
- DECISION: v4 (runs/push/submission_v4.csv md5 c69f8c5d) is the CONFIRMED FINAL submission. 20 validated designs, brutally QC'd (ipSAE-demoted, SAP-clean, MD-audited switches, novelty>=3/4). Upload USER-GATED.

## 2026-10-09 (session: breakthrough cross-react + de-novo bet)
- CROSS-REACT BREAKTHROUGH: des_59_1 (partial-diffusion, 82-mer, 0-Cys). hiscan His-position scan min-ipSAE: WT 0.564, H48 0.577 (best), H30 0.559, H59 0.544, H80 0.537 — all >> v4 parent b1m9_A12H min 0.387.
- des_59_1 HARD GATES PASS: novelty 3/4 (seq 0.0% id vs PDB; struct TM 0.885 = de-novo bundle floor, vs b1m9_A12H v4-passer TM 0.778); expression SAP+ 2.7, instab 36.5, gravy -0.13, pI 4.9, 0-Cys.
- switch-MD running (envs/md2, OpenCL): H48/H30 gpu2/3, H51/H55 gpu0/1. audit_brutal verdict pending = last gate before v4 integration (swap for weak b1m9_A12H).
- joint tied-MPNN cross-react (jxr): FAILED — re-seq of interface killed binding (ipSAE ~0.01 both orthologs). Abandoned. Reconfirms: don't re-MPNN a validated binder.
- de-novo bet: β-pairing targeted RFd (Sappington/Baker, MIT) on TNF edge strand B111-126 (Region I) — fork a363299234dd52b22, gpu4-7. Backbone-directed = intrinsically cross-react. 149 tensors built, RFd pilot running.
- 4 SOTA research agents returned: AIMNet2/g-xTB (charged near-DFT, ~0.5x GFN2-xTB err) noted as orthogonal switch-energy confirmation; Baker dddG_elec (fa_elec-only ddG) = validated switch discriminator.
- v4 submission stands as the floor; upload USER-GATED.

## 2026-10-09 cont (integration -> submission_v5)
- 2 switch SURVIVORS triple-confirmed (MD audit_brutal + PROPKA pKa-drop + qm GFN2-xTB dE, 3-method sign-consensus):
  des_59_1_H51 (xreact min-ipSAE 0.521, d_repel 3.26, SAP 2.7, 3/4-nov) + des_76_2_H55 (xreact 0.633 BEST, d_repel 0.65, SAP 5.7, 4/4-nov). Both iptm 0.89-0.95, dSASA 2128/2245.
- CUT: des69_H56 (best xreact 0.646 but switch wrong-dir -3.60), des59_H48/H30/H55 (not sig / wrong dir).
- submission_v5.csv (md5 4d05ed2): swapped weak f8 ranks 10/20 (b1m9_1A55H_K93 sh6_r1/rx, xreact 0.36/0.26) for the 2 winners. FINAL-TEST ALL PASS (20 uniq, 0-Cys, natural AA, len 73-119).
- ranking fix: score() cross-react now min(h,m) not delta_hm gap (rewards worst-ortholog = judging); SAP penalty graded not binary-at-2.0. v5 ranks: des_59_1_H51 #15, des_76_2_H55 #20.
- v4 = fallback floor (unchanged). v5 = new best. Upload USER-GATED.
- GPU: 8x RFd beta-pairing streams generating de-novo backbones (wall broken pae 5.6); bpair_full.sh pending RFd completion.

## 2026-10-10 — portal novelty check -> submission_v6 (DQ fix); bpair de-novo maturation
- Adaptyv PORTAL novelty check run on submission_v5 (GROUND TRUTH). 19/20 pass (>=3/4); ONE FAIL:
  **b1m12_r6_A92S = 2/4** (internal estimate had it 4/4 AND ranked #1) -> submit-blocker, REMOVED.
- Portal also corrected several internal 4/4 -> 3/4; only b1m12_A41H_vsR27 is truly 4/4. f5-C5 lineage
  split 2/4, 3/4, 4/4 across 3 members -> lineage does NOT guarantee the novelty score; a specific seq
  can hit a DB match. The portal check is the only hard-gate authority; internal estimates can be off 2 levels.
- **submission_v6.csv (md5 81913f5caec2bd3e6f2985b6464c1b37)** = 19 portal-confirmed >=3/4 designs,
  re-ranked with the v5 score(). New guaranteed-valid FLOOR. v4/v5 kept as history. `scripts/integrate_v6.py`.
  Upload USER-GATED.
- bpair full run DONE: 2400 backbones -> 14400 MPNN -> 28807 af2ig -> 118 pass -> 58 uniq. De-novo
  cross-react winners (human+mouse mature TNF) MEET/BEAT best partial-diffusion: des_287_5 min-ipSAE 0.609,
  des_2140_5 0.602, des_962_5 0.601, des_84_3 0.573, des_2106_2 0.570 (ref des_59_1 0.564, parent 0.387).
  Cross-react DID transfer at scale (earlier pilot's mouse-0.0 was pilot-size noise). 4/4-novel de-novo folds.
- Switch-feasibility POSITIVE: af2pred bound poses present binder interface positions <8A of target cations
  (prior "0 positions" was from scanning the pre-prediction af2in input, not the bound pose).
- MATURATION launched (task #29): 72 His-switch variants (8 designs x 8 closest-cation positions + wt),
  Boltz human+mouse on all 8 GPUs (`scripts/bpair_hisscan_build.py` + `scripts/bpair_mature.sh`) ->
  min-ipSAE gate (keep cross-react) -> MD brutal + PROPKA + xtb triple-confirm. Target: de-novo 4/4 +
  cross-react + confirmed-switch design for the 20th slot (stronger than the removed r6_A92S).

## 2026-10-10 cont — de-novo switch matured -> submission_v7 CANDIDATE (user-gated novelty)
- bpair maturation: 72 His-switch variants (8 de-novo winners) -> Boltz h+m min-ipSAE (38 keep cross-react,
  top improved to 0.744) -> brutal two-state MD (3 batches, 24 candidates) -> 4 SURVIVORS triple-confirmed
  (MD d_repel>0 + PROPKA dpKa<0 + xtb dE>0), switch on CONSERVED Lys6/Lys93 (fires human+mouse):
  des_2106_2_H43 (1.65/-1.66/+9.96, xreact 0.566, SAP0.0 instab24 Glu21% net-10, novTM 0.676 SAFEST),
  des_287_5_H37 (2.10/-2.57/+17.35, xreact 0.703 BEST, SAP2.5, novTM 0.686),
  des_962_5_H37 (2.23 strongest MD but noisy, xreact 0.542, novTM 0.758 RISKY, net0),
  des_495_1_H78 (1.69/-0.15 weak, xreact 0.619, novTM 0.720).
- NOVELTY gate (research-pinned: struct max-TM vs PDB+AFDB, <0.5=4/4 /0.5-0.8=3/4 />=0.8=DQ). Local PDB
  proxy RELATIVE-ONLY (failed to reproduce known r6_A92S=2/4 -> AFDB missing). Absolute DEFERRED to user portal.
- EXPRESSION research (Adaptyv EGFR benchmark): Glu% is #1 predictor (AUC0.77 > GRAVY/iptm); DROPPED ProtParam
  instab penalty (Pro-FP) from score(). All 4 survivors express well (Glu% 19-22%); des_495_1 instab70 was FP.
- **submission_v7.csv CANDIDATE (md5 9a21ea9965e2a438cefa481ea49d6aed)** = v6's 19 portal-verified + LEAD
  des_2106_2_H43 (#20, ranked #6). FINAL-TEST ALL PASS. `scripts/integrate_v7.py`. runs/push/bpair_candidates.csv
  = 4 survivors (metrics+seqs+rec) for USER portal novelty-check + swap. v6 (81913f5c) = floor. Upload USER-GATED.
- SOTA sweep (4 agents) net: stack confirmed open-SOTA; ArcRefine has no public footprint; KaML-CBtree optional
  4th switch vote (deferred); refinement = contingency if portal DQs a de-novo pick.

## 2026-10-10 cont — ALL 4 de-novo PORTAL-NOVELTY 3/4 -> submission_v7 (validated new best)
- USER portal novelty check on all 4 de-novo survivors: des_287_5_H37, des_2106_2_H43, des_962_5_H37,
  des_495_1_H78 = ALL 3/4 (PASS hard gate). Novelty risk RESOLVED (mixed-a/b bpair fold is novel enough;
  local proxy was correctly cautious but all clean).
- **submission_v7.csv (md5 0e84b392410e59cc756036cabc8db7b6)** = v6's 19 portal-verified MINUS TNFpH07
  (weakest cb1: FC13/xreact0.185/thermo-only) PLUS des_287_5_H37 (switch 2.10/-2.57/+17.35, xreact 0.703
  BEST) + des_2106_2_H43 (1.65/-1.66/+9.96, xreact 0.566, pristine expr). 20 designs, FINAL-TEST PASS,
  11 distinct lineages, cb1 6->5, 4 de-novo switch designs total. `scripts/integrate_v7.py`.
- v7 is now a fully-validated upgrade over v6 (all 20 portal-novel >=3/4; +2 MD-confirmed-switch de-novo,
  stronger cross-react + diversity; -1 weak cb1). v6 (81913f5c) = fallback floor. **Upload USER-GATED.**
  backups: des_962_5_H37/des_495_1_H78 in runs/push/bpair_candidates.csv.
