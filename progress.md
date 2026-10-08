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

## Next step

**Decision point (owner):** upload the 9-design portfolio to the Proteinbase portal (account j_v_v_07), or adjust
first. Deadline **2026-10-12 11:59 UTC**. Upload is user-gated — do not auto-upload.

Optional, if more diversity is wanted: the fold ceiling at the cb1 epitope is ~7 portal-novel shapes; Region-I and
K90 wall de-novo binding, so new folds are unlikely without a new bindable epitope. The strongest de-novo switch
found (`b1m7` A22H, d_switch 2.52) is **not** shipped — it fails ipSAE (0.50) and expression (instability 65).
