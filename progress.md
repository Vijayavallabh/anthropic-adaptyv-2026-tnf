# Progress

Current status and next step for the TNF-α Challenge 2 campaign. Append; don't rely on chat history.

## Status: SUBMISSION READY — upload user-gated. CONFIRMED PRIMARY = v24 (all 20 portal-PASSED, 0 re-checks)

- **CONFIRMED PRIMARY:** `runs/push/submission_v24.csv` (md5 **`34f3c096bf63dd9ef70d05308c6c3f1d`**) — 20 designs,
  validator PASS, **user portal-uploaded 2026-10-13 → ALL 20 PASSED** (novelty + internal dedup). 20 UNIQUE,
  8 both-fold multi-His (2 three-cat + 6 two-cat = the #1-axis ceiling), 19 distinct lineages, all
  direction-confirmed, **0 re-checks**. The MAX in-silico-provable set. Dossier `docs/SUBMISSION_v24_candidate.md`.
- **CONFIRMED ALT:** `runs/push/submission_v23.csv` (md5 **`65bf89bb38f5498ad90522542655eae2`**) — all 20
  portal-PASSED; 5 both-fold, 19 lineages, fewer shell-adapt members than v24. Dossier `docs/SUBMISSION_v23_candidate.md`.
- **ZERO-WORK FALLBACK:** `runs/push/submission_v19.csv` (md5 **`7756055f0b20bcd6994fefb644271db7`**) — 20 designs,
  all portal-confirmed ≥3/4 (but carries 4 cousins → ~16 unique under the portal's internal dedup; dominated by v24/v23).
- **SUPERSEDED (records only):** v22→v23, v21→v24, v20 (fold-robust predecessor), v18 (md5 `524722158467273c0e9f57f5402650ed`).
- See the SESSION 13 cont.1–11 entries at the tail for the full chain (portal dedup → v23/v24; 4-His NO-GO;
  xtb switch-magnitude flat; in-silico frontier reached). The 9-design final_portfolio below is historical/SUPERSEDED.
- **Deadline ~2026-10-19.** Full history is append-only below; the authoritative latest state is the SESSION 12
  tail. (The `runs/submission/final_portfolio` 9-design set was an earlier session-8/9 checkpoint — superseded.)
- **Repo:** `github.com/Vijayavallabh/anthropic-adaptyv-2026-tnf` (private). Upload is USER-GATED — do not auto-upload.

## Evidence (current, 2026-10-12)

- `for v in 18 19 20; do uv run python scripts/validate_submission.py runs/push/submission_v$v.csv --track 3; done` → each PASS (20 rows, 0 errors)
- md5: v18 `524722158467273c0e9f57f5402650ed` · v19 `7756055f0b20bcd6994fefb644271db7` · v20 `6f4b61d0a55de2e4e0db205970ec6821`
- `scripts/novelty_check.sh runs/push/submission_v20.csv` → 0/20 flagged >70% id to PDB (all de novo)

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

## Next step (SESSION 9 snapshot — SUPERSEDED; current submission is v19/v20, deadline ~2026-10-19; see the SESSION 12 tail)

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

## 2026-10-10 cont — adversarial stress-test (8-agent workflow) -> best-20-on-all-metrics campaign
- USER mandate: best 20 on ALL metrics at all costs, willing to discard ALL current designs.
- Stress-test workflow (wl9lmukbb, 8 independent adversarial verifiers reading evidence by bytes): under the
  strict all-metrics bar only 2 current designs are STRONG-ON-ALL: des_287_5_H37, des_2106_2_H43. 3 MARGINAL
  (b2m0_c70 mouse-biased Dhm0.318; b1m12_A41H_vsR27 best-binder+only-4/4 but modest switch; des_59_1_H51 mouse
  0.521+expr2.7). 15 DISCARD incl ALL 5 cb1 (switch thermo-only/MD-CONTRADICTED neg d_repel + weak x-react),
  des_76_2_H55 (SAP5.7 expr risk), and (caught my blind spots) b1m7_A66H (switch DECOUPLED from binding, d_iface
  0.003 = constitutive binder), b4m2_A57H (mouse 0.370), b1m9_A12H (mouse 0.387), b2m0_c20 (97% twin of c70).
- **KEY METHOD CATCH: audit_brutal gates d_repel (His-Lys separation) NOT d_iface (interface release).** Both KEEPs
  have d_iface~0 (des_287_5 0.03, des_2106_2 -0.005) -> switch MECHANISM confirmed (interface His-Lys repulsion,
  His-specific via Gln ctrl) but FUNCTIONAL release unconfirmed in 3ns (unbinding is us-scale); single-His fold-change
  is modest (Baker TNF: 79x only via buried-His-networks). des_287_5/des_2106_2 are interface (not internal/decoupled
  like b1m7) so remain valid, with honest modest-magnitude caveat.
- **BIGGEST RISK: one-deep on des_287_5_H37** (only design strong on BOTH switch AND balanced mouse x-react).
  Mouse x-react across the set is folder-optimistic; gate the new wave on balanced two-species binding (small Dhm /
  matched co-fold delta, NOT raw mouse ipSAE) AND triple-confirmed in-contact switch in the SAME molecule -> >=3-4
  dual-strong anchors.
- REFINED wave-2 gate: dual-strong (min-ipSAE>=0.55 AND Dhm<=~0.12 balanced) + triple-confirm + report d_iface;
  longer-MD (10ns) d_iface functional-switch confirmation for finalists. v7 (0e84b392) = fallback floor (contains
  discarded designs -> NOT the best-20; the new campaign targets v8). Upload USER-GATED.

## 2026-10-10 (SESSION 10): research pivot -> multi-His (Baker SOTA) -> wave-3 -> v8
- RESEARCH SWEEP (5-agent) corrected strategy: single-His = WEAK end of Baker 2025 (bioRxiv 2025.09.29.678932);
  validated recipe = >=2-3 SEPARATED His each repelling a DISTINCT conserved TARGET cation (TNF-alpha 79x).
  d_iface~0 is EXPECTED (dissociation seconds-scale, unobservable in MD); honest ceiling = thermodynamic linkage
  (PROPKA+xTB). Expression frontier AUROC~0.72 = monomer ESMFold pLDDT + SAP (ML seq-solubility tools = OOD noise).
  Cross-react: multi-state tied MPNN puts mouse in the OBJECTIVE (Odin-Multi 3-10x); footprint-count REFUTED (rho0.07).
  ArcRefine REAL but fold-preserving (no novelty). No tool breaks the flat-polar+pH6.0 wall -> execution-bound.
  Killed the futile 15ns d_iface MD; see [[lesson-ph-switch-multi-his-separated-cations-baker-sota]].
- WAVE-3 pipeline (conserved cations from human<->mouse align: idx 1,6,27,39,60,77,85,93,98,107,123; K6/K93 strong):
  TRACK A (wave3_hisscan_build.py) multi-His point-combos on top-10 balanced binders -> 139 variants -> fold h+m
  (wave3_fold.sh) -> 31 dual-strong, 11 multi-His switch (>=2 distinct cations <8A) -> triple-confirm 4 OK.
  TRACK B (wave3b_prep/mpnn/analyze) LigandMPNN redesigns interface SHELL around FIXED multi-His (backbone kept ->
  novelty=parent) -> 160 designs -> 13 switch -> triple-confirm 4 OK (interface-adaptation FIXED des_4098_1 wrong sign).
  8 triple-confirmed multi-His switches / 5 lineages; monomer pLDDT 0.79-0.90, SAP<=1.2, novTM proxy 0.59-0.76 (~3/4,
  calibrated vs portal-verified des_287_5/des_2106_2/des_59_1 at 0.56-0.67). des_4098_1 b1/b3 = HIGH novelty risk (~0.75).
- v8 (integrate_v8.py -> runs/push/submission_v8.csv md5 8566b133): 7 new multi-His switches (5 lineages) + 13 diverse
  portal-verified v7 floor (dropped 3 excess cb1 + des_76_2 SAP + b2m0_c20 twin + b1m7_A66H decoupled + b4m2 3rd).
  16 lineages/20. PENDING user portal novelty on 8 new (new_for_portal.csv); backfill = 7 dropped v7. v7=floor. Upload GATED.
- PORTAL NOVELTY (user, 8/8): ALL >=3/4 -> des_875_5_H3-57-80 **4/4**, other 7 = 3/4. v8 LOCKED (portal_novelty_wave3.json).
  des_4098_1 b1/b3 (HIGH-risk proxy ~0.75) came back 3/4 (AFDB did NOT push >=0.8).
- WAVE-4B (His-biased deeper interface-adapt MPNN on the 5 confirmed lineages, bias_AA H:1.5, 48 samples): 245
  designs -> 31 switch candidates (now 4 DISTINCT cations/design) -> 16 TRIPLE-CONFIRMED (double wave-3's 8).
  His-loaded (8-12%) but monomer pLDDT 0.82-0.88 + SAP<=0.1 (His is polar). Within-lineage UPGRADES (same verified
  backbone/novelty): des_3050_6 c39 (4-conf/4-cat/0.666/pLDDT0.88/SAP0.0), des_4098_1 c44 (0.861/dHM0.024).
- v9 (integrate_v9.py -> submission_v9.csv md5 05a34444): v8 with b9->c39, b1->c44. Switch conf-His [4,3,2,2,2,2,1].
  PENDING user portal re-check on c39,c44 (new_for_portal_v9.csv; expect 3/4; fallbacks b9/b1). v7=floor. Upload GATED.
- PORTAL NOVELTY v9 (user, 20/20): ALL >=3/4 -> 3 at 4/4 (des_3050_6_c39, des_875_5_H3-57-80, b1m12_A41H_vsR27),
  rest 3/4. c39 (4-conf-His switch) AND c44 both pass -> v9 LOCKED. portal_novelty_v9_full.json; byte-verified
  submission_v9.csv md5 05a344444cf740907e79a177d34c73f5 (submit==validated). FINAL best-20. v7=floor. Upload USER-GATED.
- CODE-REVIEW (background /code-review, 7 findings) RESOLVED (commit 2c42ceb): fixed wave4b_analyze.base_of to
  handle {scaffold}_c{i} names (the sed-copy had silently disabled the keeps-binding gate for all wave4b
  designs; re-verified c39/c44 keep binding under the active gate -> no v9 impact; switch candidates 31->26);
  added 20-unique-SEQUENCES assert to integrate_v9; wave4b.sh fail-loud empty-jobs + nvidia-smi GPU-idle guard;
  WAVE-3B->WAVE-4B docstrings; removed dead code. submission_v9.csv bytes UNCHANGED (05a34444). HARNESS SYNCED
  (AGENTS.md + session-handoff.md + docs/METHODS.md §8 + scripts/README.md + docs/PLAN_TNF.md -> v9).

## 2026-10-10 (SESSION 11, Stop-hook /goal bold campaign): v9 verified-robust; 4 bold bets A/B/C/D stress-tested
- GOAL = bold/unconventional to deadline; test rigorously, distinguish speculation from results. Synthesis
  (stress-test agent + own analysis): v9 near-optimal, multi-His separated-cation lever mined out; run cheap
  bounded gates for a diversity hedge / legibility figure, then harden + write-up + USER-GATED upload. v9 UNCHANGED.
- EXPRESSION HARDEN (dominant wet-lab-zero risk): expr_audit20.py over ALL 20 -> primary axis ESMFold monomer
  pLDDT **0.79-0.94 (every member strong)**; orthogonal SAP+/GRAVY/acidic%/netQ. 1 consensus-drop (b1m9_A12H:
  SAP+6.7 top-decile + GRAVY+0.30) but **KEPT** -> it is portal-confirmed >=3/4 so a swap trades soft yield-risk
  for HARD novelty-DQ; no lineage-preserving portal-confirmed replacement exists (v8-not-v9 = same-lineage b9/b1);
  pLDDT0.87 => primary predictor says it expresses. No action improves v9 without DQ risk. master_v9.csv completed
  (all 20 expression cols + novelty corrected: 3 at 4/4). submission md5 05a344 UNCHANGED.
- BET A (buried His-network, Baker method-2 cooperative Hill>1): cheapest probe = pre-existing buried His...cation
  pair in a de-novo core -> NONE (c39/des_180_6/des_4098_1_b14 scanned, 0 buried His-cation pairs <8A). Building one
  de-novo mutates core => risks fold + dominant expression axis for an UNVERIFIABLE cooperativity claim -> **DROPPED
  (no low-risk path)**. Honest negative.
- BET C (multi-state tied MPNN cross-react): **refuted by bytes** earlier (runs/jxr 0/256). Dropped.
- BET D (SMD functional-release, md_pull_tnf.py; Jarzynski pull-work in both His states, implicit solvent):
  **NULL — small-n false positive caught by adding reps.** c39 (His49,50,54,79) at 3 reps dW_release +16.1 kcal
  looked clean (GLN control His-specific), BUT at **9 reps dW=+0.1 +/-~11** (HIE 398.5+/-7.8 ~= HIP 398.5+/-7.1,
  consistent with ZERO) — the +16.1 was noise. c44 (His6,74) +4.88@3rep also within noise. CONCLUSION: steered MD
  at feasible sampling does NOT resolve the pH-release signal (dissociation seconds-scale) -> CONFIRMS the
  thermodynamic ceiling (PROPKA+xTB are the switch evidence, NOT SMD). Honest negative; does not weaken the
  switch-core (gated on PROPKA+xTB). runs/pull/{des_3050_6_c39_pull.json, c39_more/des_3050_6_c39_pull.json,
  des_4098_1_c44_pull.json}. [[lesson-ensemble-win-on-small-n-is-noise]] vindicated.
- EVIDENCE DOC: docs/SUBMISSION_v9_EVIDENCE.md (mechanism + 20-design table + triple-confirm incl SMD + expression +
  novelty + cross-react + portfolio + bold-campaign). Declined a quantitative fold-change (interface pKa predictors
  disagree 2-3u > pH window => false precision); conf-His count [4,3,2,2,2,2,1] is the honest selectivity proxy.
- CHAI-1 ORTHOGONAL CO-FOLD (3rd independent predictor) on the 7 switch core vs human+mouse (runs/push/
  chai_switchcore.tsv): 6/7 CONFIRM strong binding BOTH species with balanced matched delta 0.002-0.031 (c39
  0.788/0.772, b14 0.865/0.870, c44 0.863/0.870, H16-74 0.839/0.841, H52-54-79 0.722/0.708, H3-57-80 0.836/0.805)
  -> orthogonal #1+#2 confirmation + model-risk reduction. FLAG des_180_6_H9-14-38: Chai 0.42/0.37 (both species)
  vs Boltz 0.734 = genuine Boltz/Chai disagreement + highest-SAP switch-core member; KEPT (portal-confirmed, 2-His
  switch, distinct lineage, no DQ-safe swap) but logged lowest-binding-confidence. Selectivity-strength proxy ΣdE
  (honest relative, NOT a fold-change): c44 17.8 ~= b14 17.6 > c39 13.7 strongest by per-His QM; c39 top by breadth.
- CHAI v7-floor (runs/push/chai_v7floor.tsv): 9/13 strong both species; 4 flagged (TNFpH06_sh3 0.72/0.31,
  b4m2_A57H 0.83/0.11, b1m9_A12H 0.78/0.15 = mouse-weak+strong-human = likely Chai ortholog pessimism; tnf_K90
  0.42/0.20 = both-weak Boltz/Chai disagreement). ALL-20 stratification: **15/20 strong on both co-folders**;
  5 lower-confidence (2 both-species: des_180_6/tnf_K90; 3 mouse-only). ALL KEPT (Boltz primary + portal-confirmed,
  no DQ-safe swap). Honest risk disclosure, not a submission change.
- BET E (His-headroom: push confirmed-His past c39's 4). Diagnostic: each switch-core epitope contacts 5-7
  conserved cations but His-targets only 0-3 (headroom). Tested placeable NEW His: c39_E77H(->Arg1), c44_R50H
  (->Lys107), c44_E52H(->Lys60), c44_R50H_E52H. Binding gate (ipSAE vs parent): only c39_E77H + c44_R50H_E52H
  survived (singles broke binding). Triple-confirm on BOLTZ structure: c39_E77H DEGRADED (4->3 His, the mutation
  broke His49's confirm); c44_R50H_E52H looked like 4 His (+binding maintained per Chai 0.864, expr pLDDT0.82/SAP0.6).
  BUT re-triple-confirm on the INDEPENDENT CHAI structure: c44_R50H_E52H confirms only **2 His** (orig 6,74; new
  50/52 are fold-dependent artifacts, His52 flips dpKa +0.05). => His-confirmation is FOLD-DEPENDENT; the Boltz
  4-His was overcounted. **BET E NEGATIVE — REJECT c44_R50H_E52H; lever robustly mined-out.** scripts/
  triple_confirm_mut.py; [[lesson-his-confirm-is-fold-dependent]].
- FOLD-ROBUST HIS (applied the Bet-E lesson to v9 itself; runs/push/foldrobust_his.tsv): re-triple-confirmed the 7
  switch core on independent CHAI structures. Confirmed His Boltz/Chai: c44 2/2, b14 3/2, 3050_6_H52 2/2, 180_6 2/3,
  c39 **4/2** (Boltz-inflated: 2 weak same-cation on Chai), 818_4 2/1, 875_5 1/0. **Robust multi-His (>=2 both folds)
  = 5/7**; STRONGEST c44 & b14 (2 strong distinct-cation His, dE+2.5..6.8, both folds). 818_4/875_5 fold-sensitive
  (single/marginal) — kept as portal-novel hedge binders. Tempers headline counts; submission UNCHANGED. Honest #1
  bets = c44, b14. docs/SUBMISSION_v9_EVIDENCE.md §2b updated.
- OUTCOME: all bold bets (A/B/C/D + E His-headroom) explored, NONE beats v9. v9 VERIFIED-OPTIMAL on switch (robust multi-His 5/7 + selection
  optimal, cousins excluded for diversity), expression (all 0.79-0.94, no DQ-safe swap), novelty (20/20 >=3/4),
  cross-react (2 co-folders agree balanced), diversity (16 lineages). v9 is FINAL. submission md5 05a344 UNCHANGED.
  Upload USER-GATED. docs/SUBMISSION_v9_EVIDENCE.md = full record (§2b Chai, §3 SMD-null, §4 expression, §8 bold).
- BET B (Mosaic de-novo, Escalante gradient seq-space; won Adaptyv Nipah/TREM2): ran GPU0, 8 human-TNF designs ->
  6/8 are 36-47% Ala (floor mode). INDEPENDENT held-out Boltz arbiter (runs/boltz_moshuman): 3 high (mosh_0 iptm
  0.917/ipSAE_min **0.778**, mosh_5 0.861/0.67, mosh_3 0.865/0.64). **Mosaic CAN bind** (mosh_0 ipSAE 0.78 >= most
  v9) BUT fails the DECISIVE axes: (1) EXPRESSION mosh_0 GRAVY **+0.63** / mosh_5 **+0.69** (worse than ANY v9
  member; worst v9 was +0.30) + 36-38% Ala => severe aggregation, near-zero wet-lab; the 0.78 DEPENDS on the Ala
  helix that dooms expression. (2) SWITCH mosh_5=0 His, mosh_0=1 incidental His (@70) -> NO multi-His architecture,
  scores 0 on #1 axis. (3) EPITOPE diffuse 20-26 target residues (non-specific wrap), only incidental Region-I
  {114,115}/cation overlap. High iptm = oracle self-consistency trap (Ala-helix on hydrophobic patch). Salvage =
  rebuild v9 wave-3/4b+triple+portal on an inferior non-switchable seed, no time to validate -> **Bet B CLOSED**.
  JOINT bispecific NOT run (same floor mode; won't fix expression/switch). Confirms de-novo generators floor on the
  flat TNF target => v9's multi-His-on-validated-binders is correct.
- v9 remains THE submission (byte-verified 05a344444cf740907e79a177d34c73f5). Bold-bet outputs are HEDGES ONLY;
  fold in only a demonstrable pass that beats v9 AND is portal-novel. Upload USER-GATED.

## 2026-10-10 (SESSION 12, Stop-hook: resumed improving — fold-robust upgrades + v10 candidate)
- Stop-hook pushed back: v9 unchanged, I held instead of improving. CHALLENGED the "frontier" assumption + the
  self-imposed "can't change the 20 (portal-gated)" limit -> the gate is a STEP not a wall: GENERATE stronger
  designs + stage for user portal-check. Target = v9's fold-robust-weakest #1 members (des_875_5 0-robust-His,
  des_818_4_H16-74 1-robust, c39 Boltz-inflated).
- CAMPAIGN F. **F1 (mine non-v9 wave4b under the new FOLD-ROBUST bar)**: Chai-fold + Chai-triple-confirm the 14
  non-v9 Boltz-confirmed designs. All 14 bind on Chai (>=0.72). Fold-robust GEMS found:
  - **des_818_4_c13**: fold-robust 2 His on DISTINCT cations [Lys6,Lys93] (same positions both folds), dHM 0.008
    (best in pool), pLDDT 0.86, SAP 1.1, Chai bind 0.818. CLEAR UPGRADE for des_818_4_H16-74 (1 robust His, dHM
    0.105, SAP 5.4): wins #1+#2+expression, only #3 slightly lower. novelty proxy 0.885 ~= parent 0.888 (3/4).
  - **des_3050_6_c19**: fold-robust 3 His on 2 cations [Lys93,Lys60], pLDDT 0.83; upgrade for c39 (2 weak His/1
    cation). proxy 0.845 < c39 0.859 (more novel). Binding 0.586 (slightly < c39 0.666).
  - des_3050_6_c34 (fold-robust 2, dHM 0.016), des_4098_1_c39 (fold-robust 2) also solid.
- **v10 CANDIDATE** (runs/push/submission_v10.csv md5 6776d3d0): v9 with des_818_4_H16-74->c13 + des_3050_6_c39->c19.
  2 NEW designs need USER portal-check (new_for_portal_v10.csv); likely inherit parent novelty (proxy-confirmed).
  v9 remains byte-locked FALLBACK (05a344). v10 is stronger on fold-robust #1 if both pass portal.
- **F2 (multi-His on DIVERSE backbones = new lineages)**: 11 diverse v7-floor backbones contact BOTH strong anchors
  Lys93+Lys6 (des_2106_2 1.5/2.9A, des_287_5, des_59_1, b2m0, b1m7...). Ran the proven wave3b method (His-lock at
  anchors + LigandMPNN local-shell redesign) -> 176 variants -> Boltz fold h+m RUNNING (then Chai + fold-robust
  triple-confirm). Goal: strong-robust switches in NEW lineages (strong #1 + diversity). scripts/wave5_prep.py.

## 2026-10-10 (SESSION 12 cont.): F2 YIELDED — v10 candidate with 5 fold-robust upgrades (user-portal-gated)
- F2 multi-His on diverse backbones: 176 variants, 159 bind both species (Boltz); top-2/lineage triple-confirmed
  on BOTH folds. Fold-robust (Boltz∩Chai 2-His/2-cation) WINNERS: b1m7_w5_3 (bind 0.879/0.878 Δ0.001, SAP 1.1),
  b4m2_r1_w5_3 (0.796/0.825). b1m9_w5_3 = 1 robust His (mild). Non-robust (fold-dependent His, dropped): des_2106_2,
  b1m12, des_287_5 (Chai-weak bind), b2m0 (0 His).
- v10 CANDIDATE (submission_v10.csv md5 4e0b4ce6, 0-Cys, 16 lineages): 5 same-lineage fold-robust upgrades vs v9 --
  des_818_4_H16-74->c13S (2 His, dHM~0.008, SAP0.9; Cys->Ser fix of c13), des_3050_6_c39->c19 (3 His/2 cat),
  b1m7->b1m7_w5_3 (2 His), b4m2_r1->b4m2_r1_w5_3 (2 His), b1m9_A12H->b1m9_w5_3 (1 His + better expr). ~DOUBLES strong
  fold-robust multi-His switches (~5->~9) while preserving diversity. 5 NEW need USER portal-check (new_for_portal_v10;
  MPNN reseqs on portal-confirmed parent folds -> novelty inherits, proxy-confirmed). docs/SUBMISSION_v10_CANDIDATE.md.
- **CAUGHT a liability**: F2 MPNN (no omit-Cys) put 1 Cys in c13 -> Cys->Ser (c13S) re-folded+re-confirmed 0-Cys,
  fold-robust 2 His intact, binding 0.895/0.807. v10 is 0-Cys.
- v9 (05a344) remains byte-locked FALLBACK. Upload USER-GATED. If any new design <3/4 portal -> revert that swap.

## 2026-10-10 (SESSION 12 cont.): wave5b/5c exhausted the upgrade search — v10 (5 upgrades) is the improved frontier
- wave5b (more variants of non-yielding F2 lineages) + wave5c (des_875_5, des_180_6 F2 treatment, --omit_AA C):
  all BIND well both folds (Chai 0.84-0.88) but FAIL the fold-robust >=2-His gate. b1m12_2A67H_w5_5 = fold-robust 1
  His (2nd fails on Chai); des_875_5_w5c_1 / des_180_6 = fold-robust 0-1. CONCLUSION: these lineages are
  GEOMETRY-LIMITED -- they contact Lys93+Lys6 but cannot present 2 fold-robust His (a backbone property, not fixable
  by more MPNN). Honest negative; upgrade search EXHAUSTED.
- **FINAL improved candidate = v10** (5 fold-robust upgrades: c13S, c19, b1m7_w5_3, b4m2_r1_w5_3, b1m9_w5_3;
  md5 4e0b4ce6, 0-Cys). ~Doubles strong fold-robust multi-His switches (~5->~9) vs v9, 16 lineages preserved.
  5 NEW need USER portal-check (new_for_portal_v10.csv; novelty inherits from portal-confirmed parents, proxy-OK).
  v9 (05a344) = byte-locked FALLBACK. Upload USER-GATED. docs/SUBMISSION_v10_CANDIDATE.md.
- Members NOT upgradeable (stay as v9/v10): des_875_5 (0 robust His, but 4/4-novel binder), des_180_6, des_2106_2,
  des_287_5, des_59_1, b2m0, b1m12_2A67H, TNFpH06, tnf_K90 -- geometry-limited or already-robust. v10 is the ceiling
  of the fold-robust-switch-upgrade campaign.

## 2026-10-10 (SESSION 12 cont. — Stop-hook persistence): BOLD de-novo RFd3 campaign for NOVEL-fold (4/4) switches
- Hook pushed for SUSTAINED bold exploration, not convergence. Challenged the deepest assumption: that the
  submission must be built from EXISTING (3/4-novel, derived) backbones. BOLD play = de-novo generation of BRAND-NEW
  folds designed for the switch -> potential **4/4-novel AND strong-switch** designs (dominate novelty + #1).
- Prior RFd3 "tnf_K90" run gave 2 de-novo binders pae<15 (model_3_4 pae13 plddt90 rmsd0.67; model_15_8 pae14) BUT
  they contact 0 conserved cations (bind a cation-free non-conserved epitope) -> DEAD END for switch + cross-react.
  Lesson: the de-novo binder must target the CONSERVED-CATION epitope.
- Launched scaled RFd3 de-novo on the REGION-I spec (runs/rfd3_v11, GPUs3-6, ~1024 designs, SEED_BASE=5000): hotspots
  V74/L75/I97/W114/Y115 neighbor Lys6/Lys93 -> binders will be switch-capable. All 4 GPUs healthy. Downstream (hours):
  eval_rfd3_regionI.sh (MPNN soluble 0-Cys + af2ig) -> de-novo binders (plddt_binder>80, pae_int<15) -> add multi-His
  switch (F2-style His-lock at cations) -> fold-robust triple-confirm (Boltz+Chai) -> cross-react + expression +
  local novelty (expect 4/4) -> v11 novel-fold strong switches.
- v10 (5 fold-robust upgrades, md5 4e0b4ce6) remains the current best improved candidate; v9 (05a344) fallback. The
  RFd3 campaign is the sustained high-ceiling play for genuinely novel additions. Upload USER-GATED.

## SESSION 11 (2026-10-11) — v12→v16, the two-state-MD trigger matrix, and the de-novo LONG SHOT

- **Two-state MD His-causality gate adopted as the #1 (pH-selectivity) decider.** dSwitch = d(His-Lys|HIP) −
  d(His-Lys|HIE), CONFIRMED iff >2·SEM AND >0 AND dIso(HIP-GLN)>0. Measures the charge TRIGGER (repulsion),
  not wet-lab release (d_iface≈0 in 3ns is expected). Full 20-design v14 matrix: **13 CONFIRMED / 7 NONRESPOND**
  (2 wrong-sign) — `docs/TRIGGER_MATRIX_v14.md`. STATIC Wyman linkage ≠ MD trigger (reconfirmed repeatedly).
- **v12→v13→v14** (verified axis-complete swaps; v14 `7994cf46` = v13 + fix worst cross-reactor; all portal ≥3/4).
- **Proton-PottsMPNN** (SOTA pH inverse-folder) tested → ANTI-correlated with our MD (AUC 0.167) → NOT adopted.
- **LONG SHOT (user directive "go with the long shot"): de-novo NEW shapes for the weak NONRESPOND slots.**
  β-pairing RFd → 480 new shapes → 3 cross-reactive (des_138_7, des_4082_8, des_2044_8) → **partial-diffuse
  AMPLIFY** (the multiplier, ~25% yield) → cross-react → multi-His → two-state MD → **5 MD-CONFIRMED new-shape
  switchers** (lead parent was MD-NONRESPOND; only diversified cousins confirmed). All portal-novelty **3/4**
  (seq-novel + structural-capped; calibration: TM 0.86-0.93 still 3/4 — refutes the old "≥0.8→2/4" heuristic).
- **v16 `runs/push/submission_v16.csv` md5 `a524d77fa7a63ef33b28346ce3a7d87d`, validator PASS = CURRENT BEST.**
  v14 + 4 new-shape switchers (2 new lineages, ≤2 each) for the 4 weakest v14 NONRESPOND → portfolio MD trigger
  13/20→17/20 CONFIRMED, x-react 0.87-0.91, diversity preserved. Fallback v14 (`7994cf46`). Upload USER-GATED.
- **v17 candidate (RESOLVED, USER-PORTAL-GATED):** de-novo lineage **des_6025 yielded 1 MD-confirmed
  new-lineage switcher** `des_6025_4_T20_pd_64_1_mh_34_50_52_74` (des_6126/des_8xxx = dead ends). Re-validated
  on the exact His-variant: d_repel **+0.372 CONFIRMED** (isostere ctrl +0.308), cross-react **re-folded**
  human 0.962 / mouse 0.938 (min 0.938 PASS), solub 0.87, **seq-novel (no PDB hit)**. ONE open risk: foldseek
  structural TM **0.945** (above v16's 0.86-0.93) → likely 3/4 but portal-UNVERIFIED. `submission_v17.csv`
  md5 `9357b6905ed81179da5137a16dd9dd5a` (validator PASS) = v16 with `b4m2_A57H`→des_6025 (18/20 MD-confirmed,
  3 de-novo lineages). **Decision: default upload v16; if a portal novelty check can be spared, check des_6025 —
  ≥3/4→v17, else v16.** Dossier `docs/SUBMISSION_v17_candidate.md`. De-novo long shot now played out.
- Pipeline drivers promoted to `scripts/` (denovo_stage2, partial_denovo, v15_xreact, v15_multihis, v15_stage2);
  see `scripts/README.md` + `session-handoff.md`. GPUs 0,1,6,7 intermittently foreign (VLLM/scGPT) — never clobber.

## SESSION 11 cont. → v18 / v19 (both fully gated)
- **v18 `runs/push/submission_v18.csv` md5 `524722158467273c0e9f57f5402650ed`** (validator PASS) = v16 with 2
  swaps (des_3050_6→des_6126 [4th lineage, MD d_repel 1.36], b4m2→des_6025 mh_52_74 0.93). 19/20 MD-confirmed,
  4 de-novo lineages, all ≥3/4 novel. Genie3-groove hedge = rigorous NO-GO (af2ig floors; "confirmed" switches
  were pose-instability artifacts, binder drifts 4.5-7.1Å). Pose-consensus binding oracle = NO-GO (AUC 0.39-0.55).
- **v19 `runs/push/submission_v19.csv` md5 `7756055f0b20bcd6994fefb644271db7`** (validator PASS) = v18 with
  des_818_4_c13S → des_84_3_mh_49_56_76 (more fold-robust: 1 both-fold His). **User portal-checked: all 20 ≥3/4,
  FULLY GATED.** 20/20 switchers. = CURRENT BEST. Deadline extended to ~2026-10-19.

## SESSION 12 (2026-10-12) — v20 fold-robust shell-adaptation (17→19/20 fold-robust; the live #1-axis win)
- **Hypothesis VALIDATED:** shell-adaptation (LigandMPNN shell-redesign, His FIXED, repack interface + His
  microenv; `scripts/shelladapt_shell.py`) upgrades fold-FRAGILE switchers to fold-ROBUST (switch-His engages a
  cross-conserved cation in BOTH Boltz∩Chai). sa2 campaign: **9/11 fragile members upgraded**; round-2: +2/3 of
  the stubborn (b2m0_c70, des_875_5). Portfolio **8/20 → 17/20 (→19/20) fold-robust.**
- **v20 `runs/push/submission_v20.csv` md5 `6f4b61d0a55de2e4e0db205970ec6821`** (validator PASS) = v19 with 9
  fold-robust swaps. Dossier `docs/SUBMISSION_v20_candidate.md`; picks `runs/foldrobust/v20_final_picks.tsv`.
- **Binding re-confirmed per design** (shell-adapt rewrote 60-75/80 residues): af2ig initial-guess on each NEW
  sequence, all 9 adopted PASS (8/9 match/beat parent). The most-soluble non-binding variants (sa5/sa3, pae 23/22)
  correctly rejected — binding filter is load-bearing. des_6025 sa10 binding confirmed (4.78/92/0.45).
- **Cross-react** 89-95% epitope id, all cations cross-conserved (`crossreact_interface_conservation.md`).
  **Novelty** sequence 0/20 >70% id (all de novo); structural ~ parent (backbone unchanged) → USER portal re-check
  REQUIRED on the reseqenced swaps. **Filter stringency:** both-fold 28% vs single-fold 51% across 176 variants.
- **Two-state MD (physics corroboration, 3-rep):** fold-robust ≠ all MD-selective — des_6025 sa10 + des_4082_8
  clean (2 His, 2 SEPARATED cations); the bc=3 des_2106_2 MD-weak (2 His crowd K93 = anti-cooperative). Vindicates
  separated-His; MD is a noisy confirmer, not a hard cutter — v20 selection stands on fold-robustness + binding.
- **Honest scope:** v20 is a higher-CONFIDENCE switch-geometry gain on #1 (more members reliably switch), NOT
  in-silico binding proof (TNF docking stays wall-limited; #3 affinity wet-lab-decided). Campaign-2 (1200 fresh
  bpair-RFd → funnel) in flight for new-lineage binders (novelty/diversity). v18/v19 gated; v20 USER-portal-gated.


## SESSION 13 (2026-10-12) — v21 leaner DIVERSE portfolio (brutal cut + new-lineage refills)
```json
{
 "title": "SESSION 13 (2026-10-12) \u2014 v21 leaner DIVERSE portfolio (brutal cut + new-lineage refills)",
 "v21_csv": "runs/push/submission_v21.csv",
 "v21_md5": "30915564027fdcea0da997675d29e904",
 "v21_n": 13,
 "v21_status": "validator PASS, LOCAL-ONLY, upload USER-GATED",
 "v21_vs_v20b": "dropped same-fold cousins + weakest-switch members; +3 NEW-lineage regen refills; sh3_sa11 replaces sa10",
 "dossier": "docs/SUBMISSION_v21_candidate.md (JSON-dump of validation_results)",
 "brutal_review": "4 subagents (redundancy / novelty-DQ / group-A / group-B) -> cut the majority; refill with foldseek-DISTINCT new lineages, not cousins",
 "KEY_LESSON": "switch MAGNITUDE not reliably rankable in-silico (pKa err 2-3u > 1.4u window; PROPKA vs pKAI rank DISAGREE rho=0.56). Only release DIRECTION robust (all designs sign-agree). Select on direction + reliable axes (bind/xreact/novelty/sol) + DIVERSITY, NOT magnitude rank.",
 "sh3_sa11": "kept for novelty-safety (struct-TM 0.913 vs sa10 0.959) + binding, NOT a claimed sharper switch",
 "regen_refills_NEW_lineages": {
  "des_732_1_r3_sa15": "3-cation; MD dir 0.568>>Gln-0.322 OK; xr93 nov0.842 sol0.72",
  "des_1169_2_r3_sa7": "MD dir 0.627>>Gln-0.738 OK; xr97 nov0.934 sol0.83",
  "des_243_3_r3_sa2": "best novelty 0.808; MD dir running; xr92 sol0.815 bind-marginal(plddt83)"
 },
 "expansion_in_flight": "regen round-3 (8 new-lineage substrates, GPUs 1-7) + des_243 MD (GPU0) -> fill remaining slots with more foldseek-distinct direction-confirmed switchers (lottery tickets on #1 + expression hedge)",
 "fallback": "v19 (7756055f, portal-confirmed >=3/4 all 20) remains the gated default"
}
```


## SESSION 13 cont. (2026-10-12) — v21 EXPANSION COMPLETE (14) + code-review doc fixes
```json
{
 "v21_final": "14 designs, md5 87066005a74e38ace5984742b6b35a73, validator PASS",
 "expansion_result": "regen round-3 (8 new-lineage) ALL 1-cation (0 reached >=2); 5 bound af2ig but MD direction killed 4 (3 wrong-sign, 1 non-causal). round-2 cousins: des_1120_2 no-bind, des_732_r2_1 cousin-of-shipped, des_820_3_r3_sa4_r2_8 (2-cation, NEW lineage) binds + MD His-causal. NET +1 -> 14.",
 "des_820": "NEW lineage des_820, 2-cation[6,93], af2ig 5.26/87.5/1.84, MD d_repel 0.728>>Gln-0.555 (strongest gap seen), xreact 90.9%, novelty TM 0.575 (best; maybe 4/4), sol 0.659, esm-mono 0.52 (expression flag, plausibly ESMFold pessimism on a novel fold; offset by Boltz+Chai fold + af2ig plddt87.5)",
 "discipline": "rejected every wrong-sign/non-causal/cousin/non-binder -> NO dilution; 14 distinct-fold validated members; pool exhausted (all 23 campaign-2 binders processed).",
 "code_review_fixes": "/code-review (max) found 14 doc inconsistencies; live CSV verified correct (md5 match, 14 rows, 0 Cys, validator PASS). Fixed: byte-identity CORRECTED to 4 byte-identical-v19 / 10 NEED-recheck (earlier '10 byte-identical' was wrong - shell-adapt saN were never portal-checked; DQ-risk averted); WYMAN/REPLACE records (sh3_sa11 shipped for novelty not switch); des_2106_2 'MD-weak' retracted; README v20 portal caveat restored; novelty-TM wording (portal 3/4 is the gate, not TM<0.80); added authoritative UPLOAD_CHECKLIST + FINAL_COMPOSITION to dossier + superseded-record notice for the raw dumps.",
 "upload": "USER-GATED; 10 members need portal novelty re-check before upload; v19 remains gated-safe fallback."
}
```


## SESSION 13 cont.2 (2026-10-13) — BOLD de-novo campaign: +3 ultra-novel switchers -> v21=17
```json
{
 "v21_final": "17 designs, md5 f7f246f41a3a062b1922f0f980a9928a, validator PASS, LOCAL-ONLY",
 "bold_campaign": "fresh bpair-RFd de-novo (1200 backbones, all 8 GPUs) -> 49 af2ig binders (vs campaign-2's 23, >2x) -> 10 ultra-novel (foldseek TM<0.82) switch-piped = round-4.",
 "round4_shelladapt": "2 three-cation + 3 two-cation both-fold winners (vs round-3's all-1-cation) \u2014 fresh backbones had better cation-reach geometry (3 distinct cations closely engaged).",
 "round4_MD_direction": "3 of 3 >=2-cation binders PASSED His-causal gate: des_36133_g4_sa18 (3-cat[6,60,93] gap+2.30), des_34125_g4_sa15 (2-cat[6,93] gap+1.41), des_33012_g4_sa8 (2-cat[6,93] gap+0.23 thin). Both ultra-novel 1-cation binders FAILED (des_30106 wrong-sign, des_36087 non-causal).",
 "added": "+3 new-lineage ultra-novel 2-3 cation members -> v21 14->17. All bind (af2ig pae 4.5-5.7), xreact 91-93%, esm-mono 0.87-0.90, sol 0.66-0.78, nov TM 0.82-0.85.",
 "lesson_reconfirmed": "2-cation is the RELIABLE switch class; 1-cation unreliable (ALL round-3 + round-4 1-cation failed direction). Select only >=2-cation direction-confirmed.",
 "upload": "USER-GATED; 13 of 17 need portal novelty re-check (only 4 byte-identical to portal-confirmed v19); v19 gated-safe fallback."
}
```


## SESSION 13 cont.3 (2026-10-13) — HONEST v21-vs-v19 head-to-head (rigorous)
```json
{
 "finding": "Ran identical in-silico gauntlet on v19 (20) and v21 (17). v21 is NOT 'better in every metric'.",
 "binding_floor": "TIE (both all-pass af2ig pae_int 4-6)",
 "cross_react": "TIE (~90-97% both)",
 "expression": "~TIE (v21 better floor; cut des_84_3 sol 0.452)",
 "switch_direction_#1": "TIE \u2014 v19 raw multi-His ALREADY 13/15 MD His-causal, strong gaps +0.06..+2.3; v21 also ~all His-causal. H43 MD-strong +1.85 despite weak static Wyman; cut des_84_3 MD-strong (+2.14) - cut for expression/novelty not switch.",
 "novelty": "v19 WINS (all 20 portal-confirmed >=3/4 vs v21 13 pending re-check)",
 "fold_robustness_diversity": "v21 WINS (Boltz\u2229Chai both-fold confidence + 15 lineages)",
 "VERDICT": "v19<->v21 is a TRADEOFF not a domination. v21 only marginally better IF portal re-check clears (same measurable axes + diversity); else v19 is the safe strong pick (confirmed novelty + comparable everything). Corrects the earlier 'v21 strongest' overclaim."
}
```


## SESSION 13 cont.4 (2026-10-13) — v22: pragmatic Pareto-improvement over v19 (see cont.5 for the #1-axis correction)
```json
{
 "v22": "runs/push/submission_v22.csv md5 abaec06790766f25282133d929083187, 20 designs, validator PASS, LOCAL-ONLY",
 "thesis": "v22 = v19 minus its 2 MD-FAILING switchers (b1m9_A12H d_repel -1.307 wrong-sign; des_138_7_pd_24_4 -0.097 non-causal) + 2 strongest direction-confirmed NEW-lineage (des_36133_g4_sa18 3-cat gap+2.30; des_732_1_r3_sa15 3-cat). Constructed FROM the rigorous v21-vs-v19 gauntlet, which located v19's only real weakness.",
 "pareto_dominance_vs_v19": ">= v19 on EVERY axis; STRICTLY > on switch-DIRECTION (fixes 2 dead switchers) + diversity; TIE on binding/cross-react/expression; = novelty if the 2 added pass re-check. Only 2 portal re-checks (vs v21's 13); 18 members already portal-confirmed >=3/4. CORRECTION: earlier 'strictly > on binding' was an overclaim (pae 4.5 vs 4.9/6.0 is within noise, does not predict wet-lab binding) -> TIE.",
 "recommendation": "PRIORITY-DEPENDENT, not a flat v22-primary. v22 = pragmatic strict-improvement over the portal-confirmed v19 baseline at minimal (2) re-check risk. v19 = zero-work fallback. v21 = max-#1-axis + max-diversity (13 re-checks)."
}
```

## SESSION 13 cont.5 (2026-10-13) — SWITCH_FC_ORDINAL: v21 leads the #1 axis (correction to recommendation)
```json
{
 "insight": "Switch MAGNITUDE is unrankable in-silico, BUT the both-fold-confirmed COUNT of separated in-window His is a sound ORDINAL for expected fold-change (Wyman linkage is additive; Baker multi-His is the SOTA pH-switch). It is also a geometric-confidence signal (His engages its cation in BOTH independent folders).",
 "both_fold_counts": {"v21": "3 three-cat + 5 two-cat = 8 members with 2-3 confirmed His (+4 carryover)", "v22": "2 three-cat (+18 v19-derived, MD-direction-OK but NOT both-fold geometry-confirmed)", "v19": "0 both-fold-confirmed (20 raw multi-His, MD-direction-OK)"},
 "implication": "On the #1 judged axis (pH-selectivity / expected FC), v21 is STRONGEST via both-fold His-count depth. So the recommendation is PRIORITY-DEPENDENT: v21 for max #1-axis strength (13 re-checks); v22 pragmatic Pareto-over-confirmed-v19 (2 re-checks); v19 zero-work. NOT a flat 'v22 primary'.",
 "also_corrected": "v22-vs-v19 binding claim: 'strictly better' -> TIE (pae_int 4.5 vs 4.9/6.0 is within af2ig noise; the real v22 gain over v19 is switch-DIRECTION [2 non-switchers -> 2 direction-confirmed] + diversity).",
 "recorded_in": "validation_results.json keys SWITCH_FC_ORDINAL + COMPARISON_v21_vs_v19; memory lesson-switch-magnitude-not-rankable-select-diversity.md (refined)."
}
```

## SESSION 13 cont.6 (2026-10-13) — BOLD 4-His campaign (regen5): NO-GO, both-fold ceiling=3 confirmed
```json
{
 "thesis": "Challenge the 'bold space exhausted' claim on the TOP axis: v21's #1-axis ordinal tops out at 3 both-fold cations; a 4-His shell-adapt campaign had never been run. If a variant locked 4 both-fold + binding + MD direction -> beats v21 on pH-selectivity -> v23.",
 "geom_gate_PASSED": "geom4_screen over 34 af2ig-passing bpair binders: 10 reach 4 separated conserved cations (CB<=9A, seq-sep>=8). chainB=mature+80 numbering VERIFIED (all 11 CONS positions are K/R vs human mature TNF).",
 "campaign": "regen5: 8 distinct backbones (6 NEVER tried at 4-His incl des_30043 tight cluster all d<=4.7A), 4 His each, 256 LigandMPNN shell-adapt variants, fold-confirmed Boltz∩Chai.",
 "result": "MAX both_cons=3 (des_32121_g5 sa26 [6,60,93]); des_30043 tight cluster -> only 1 both-fold; NO variant locked 4.",
 "verdict": "NO-GO. both-fold ceiling=3 CONFIRMED by a 2nd independent test (regen4's 4-His installs also maxed 3; des_30106 4-His->1, MD d_repel -0.07 dead). The 4th His does not lock in both folds even from 4-cation geometry. No step-change over v21; no new member adopted (des_32121_g5 sa26 is only a 3-cat = v21's level). Thread CLOSED.",
 "recorded_in": "validation_results.json REGEN5_4HIS_CEILING; memory lesson-his-confirm-is-fold-dependent.md.",
 "recommendation_unchanged": "v21/v22/v19 stand (priority-dependent). This result re-confirms v21 as the #1-axis frontier."
}
```

## SESSION 13 cont.7 (2026-10-13) — USER portal-uploaded v22; portal cousin-dedup -> v23 (dedup-clean, pragmatic best)
```json
{
 "portal_v22_upload": "USER uploaded v22 (20). Portal parse: 16 PASSED (incl des_36133_g4_sa18 + des_732_1_r3_sa15 -> 2 new members now CONFIRMED novelty+dedup); 4 NOT passed = INTERNAL cousin-dedup ('Similar designs in this upload'): des_4098_1_c44, des_4098_1_c44_w6_3, des_2106_2_H43, des_2106_2_mh_37_43_50_75. => v22 scores only 16 UNIQUE.",
 "lesson": "portal dedups near-cousins within an upload (sequence-similarity). Hits v19/v22 (4 cousins) and v21 (2). des_4082_8 x2 PASSED (distinct enough). Vindicates diverse-not-cousins: cap 1 per too-similar cluster.",
 "v23": "runs/push/submission_v23.csv md5 65bf89bb38f5498ad90522542655eae2, 20 designs, validator PASS, LOCAL-ONLY. = v22 - {c44, c44_w6_3, des_2106_2_H43} + {des_34125_g4_sa15, des_820_3_r3_sa4_r2_8, des_33012_g4_sa8}. 20 UNIQUE scored + #1-axis both-fold depth 5 (vs v22's 2). 19 distinct lineages. Re-check: 3 refills (16 already portal-passed).",
 "recommendation": "v23 = pragmatic BEST (dedup-clean, 3 re-checks). v21 = max #1 (8 both-fold, 13 re-checks, needs own dedup). v19 or v22-minus-4 = zero-new-recheck 16-unique. v22 SUPERSEDED by v23.",
 "also": "4-His campaign (cont.6) = NO-GO, both-fold ceiling 3 confirmed."
}
```

## SESSION 13 cont.8 (2026-10-13) — v23 FULLY portal-confirmed; v24 = max-#1 dedup-clean
```json
{
 "v23_confirmed": "USER uploaded v23 -> ALL 20 PASSED ('20 included, 0 needing fixes'). Novelty + internal-dedup both clear. v23 now FULLY portal-confirmed, 20 UNIQUE, 0 re-checks. md5 65bf89bb38f5498ad90522542655eae2.",
 "v24": "runs/push/submission_v24.csv md5 34f3c096bf63dd9ef70d05308c6c3f1d, 20 designs, validator PASS, LOCAL-ONLY. = v23 + 3 same-lineage swaps to v21's SOLID >=2-cat both-fold saN (des_2106_2 sa2, des_4082_8 sa7 d_repel1.69, des_6025_4 sa10_sh3_sa11 d_repel1.2). 8 both-fold (2 three-cat + 6 two-cat) = v21's MAX #1-depth, dedup-clean, only 3 re-checks (vs v21's 13). EXCLUDED des_6126_3 sa1 (1 both-fold cat only) + b2m0 sa17 (weak) after checking v20 bothfold data.",
 "recommendation": "v24 = max #1-axis (3 re-checks). v23 = fully-confirmed-now (0 re-checks). v19 = fallback. v22->v23, v21->v24 SUPERSEDED.",
 "domination_honest": "NOT in-silico-provable: binding magnitude (#3) wet-lab-decided (Anthropic AP~0.16 near-chance); switch MAGNITUDE thermodynamically unrankable. Provable axes MAXED: novelty 20/20 confirmed (v23); #1 both-fold ordinal=8 ceiling (4-His NO-GO cont.6); diversity 19 lineages; direction all-confirmed."
}
```

## SESSION 13 cont.9 (2026-10-13) — v24 FULLY portal-confirmed = the MAX in-silico-provable submission
```json
{
 "v24_confirmed": "USER uploaded v24 -> ALL 20 PASSED ('20 included, 0 needing fixes'). The 3 saN swaps (des_6025_4_sa10_sh3_sa11, des_2106_2_mh_sa2, des_4082_8_sa7) PASS novelty + dedup. v24 FULLY portal-confirmed, 20 UNIQUE, 0 re-checks. md5 34f3c096bf63dd9ef70d05308c6c3f1d.",
 "v24_is_the_frontier": "v24 maxes every PROVABLE axis AND is portal-confirmed: novelty 20/20; #1-axis both-fold ordinal=8 (the CEILING, per 4-His NO-GO cont.6); diversity 19 lineages; switch DIRECTION all-confirmed. = MAX in-silico-provable submission.",
 "recommendation": "v24 = CONFIRMED PRIMARY (max #1, 0 re-checks). v23 = confirmed alt (5 both-fold, fewer shell-adapt). v19 = fallback. v22->v23, v21->v24 SUPERSEDED.",
 "domination_honest": "binding magnitude (#3) + switch MAGNITUDE remain in-silico-unprovable (Anthropic AP~0.16 near-chance; pKa err>window) -> further gains are WET-LAB only. v24 is the ceiling of what in-silico establishes. Upload is USER-GATED."
}
```

## SESSION 13 cont.10 (2026-10-13) — MEASURED switch magnitude (xtb dE): flat across v19/v23/v24
```json
{
 "probe": "Hook pushed: switch-magnitude never MEASURED as an improvement. Answered directly: aggregated per-His GFN2-xTB dE_switch (kcal) per design (the triple-confirm already computes dE; the both-fold COUNT just thresholds it away). Deterministic QM, not MD/SMD/pKa-predictor.",
 "single_fold_total_dE_mean_kcal": {"v24":5.7,"v23":6.1,"v19":6.6},
 "both_fold_total_dE_mean_kcal": {"v24":3.0,"v23":3.0,"v19":2.3},
 "finding": "RAW switch magnitude is FLAT (v19 marginally highest). v24 does NOT beat v19 on measured switch MAGNITUDE; its edge is both-fold CONFIDENCE (+30%) + dedup (20 vs 16 unique) + diversity + direction. Per-member dE noisy/fold-dependent (4 MD-direction-confirmed members read dE=0) -> not a reliable rank -> a dE-informed swap would chase noise. NO v25.",
 "verdict": "Directly measured the axis the hook flagged; the measurement CONFIRMS switch magnitude is not an improvable/rankable axis in this pool (vindicates lesson-switch-magnitude-not-rankable with a direct QM energy). v24 stands as the confirmed primary. Remaining gains are wet-lab-only."
}
```

## SESSION 13 cont.11 (2026-10-13) — in-silico FRONTIER reached; v24 is the max; every scored axis verified
```json
{
 "scored_axes_on_v24": {
  "#1_pH_selectivity": "direction all-confirmed (MD); both-fold His count=8 (CEILING per 4-His NO-GO); MAGNITUDE measured FLAT (xtb dE v24 5.7 ~ v19 6.6) -> not improvable",
  "#2_cross_react": "all members ~90-97% interface conservation (v24 ties v19); x-species affinity magnitude is wet-lab (co-fold ipTM = folder pessimism)",
  "#3_affinity": "af2ig floor passed by all; RANK walled (Anthropic AP~0.16; ground-truth binders unposeable)"
 },
 "gates": {"novelty":"20/20 portal-PASSED (>=3/4 gate; 4/4 not separately scored)","dedup":"20 UNIQUE, 19 lineages","expression":"all solubility-screened"},
 "conclusion": "v24 maximizes every in-silico-improvable component of every scored axis and is fully portal-confirmed. The MAGNITUDE components (binding, switch, x-species) are wet-lab-decided -> 'dominate every metric in-silico' is not achievable for THIS challenge (established by rigorous test: 4-His NO-GO, xtb flat, AP~0.16). Further in-silico iteration = churn. HOLD at v24; remaining gains require wet-lab. Recorded INSILICO_FRONTIER_REACHED."
}
```
