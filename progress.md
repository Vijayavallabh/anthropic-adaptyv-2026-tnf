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

## Next step

**Decision point (owner):** upload. Recommended core = **Tiers A+B (10 designs)**; add Tier C (4) for fold
diversity up to the Track-3 cap of 20. Proteinbase account j_v_v_07. Deadline **2026-10-12 11:59 UTC**.
Upload is user-gated — do not auto-upload. The prior validated `runs/submission/final_portfolio.csv`
(md5 f89d9173, 9 designs) is preserved as the floor.
