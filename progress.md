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

**Honest caveats:** rescue cleanly delivers 2 de-novo folds (C2, C5); the other de-novo folds (C1/C8/C11/C12)
are not surface-rescuable (ship originals, expression-flagged). `b1m12`'s switch is the weakest (parent
d_switch 0.597) and its rescued MD reads borderline — rely on PROPKA there. Rescued variants need a
Proteinbase novelty re-check (same fold ⇒ expected ≥3/4); `A20V` (new 7th fold) needs its first portal check.
New open generators (pxdesign/complexa) were deprioritized — install friction + full re-validation exceed the
time budget for low expected value on this hard target.

## Next step

**Decision point (owner):** upload. Recommended core = **Tiers A+B (10 designs)**; add Tier C (4) for fold
diversity up to the Track-3 cap of 20. Proteinbase account j_v_v_07. Deadline **2026-10-12 11:59 UTC**.
Upload is user-gated — do not auto-upload. The prior validated `runs/submission/final_portfolio.csv`
(md5 f89d9173, 9 designs) is preserved as the floor.
