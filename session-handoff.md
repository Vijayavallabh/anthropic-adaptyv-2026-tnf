# Session handoff — TNF Challenge 2

Updated 2026-10-10. Deadline **2026-10-11 23:59 AoE** (2026-10-12 11:59 UTC).

## Current best submission (validated)

- **`runs/push/submission_v7.csv`** (md5 `0e84b392410e59cc756036cabc8db7b6`) — 20 designs, every one
  **portal-novelty ≥3/4**, FINAL-TEST pass. This is the recommended upload.
- **Fallback floor:** `runs/push/submission_v6.csv` (md5 `81913f5c`), 19 portal-verified designs.
- **Upload is USER-GATED** — do not auto-upload any submission. (Submission CSVs live under `runs/`,
  which is gitignored; the `scripts/integrate_v*.py` scripts are the committed source of truth. Rebuild:
  `.venv/bin/python -I scripts/integrate_v7.py`.)

## What v7 is

v6's 19 portal-verified designs **minus** `TNFpH07` (weakest cb1 anchor) **plus** two de-novo β-pairing
MD-triple-confirmed pH-switch binders:
- `des_287_5_H37` — switch d_repel 2.10 / dpKa −2.57 / xtb +17.35; cross-react min-ipSAE 0.703 (best); SAP 2.5.
- `des_2106_2_H43` — switch 1.65 / −1.66 / +9.96; cross-react 0.566; pristine expression (SAP 0.0, Glu 21%).
Both switch on target Lys93 (conserved human+mouse). 11 distinct lineages; 4 de-novo switch designs total.

## Backups (portal-3/4, in `runs/push/bpair_candidates.csv`)

- `des_495_1_H78` — triple-confirmed but weak switch (dpKa −0.15); valid backup.
- `des_962_5_H37` — **FAILED** triple-confirm (PROPKA dpKa +1.05, wrong sign) despite strongest raw MD;
  do NOT use (backwards switch).

## Open items

- **Upload decision** (user): upload v7, adjust composition (keep all 6 cb1 / add a 3rd de-novo), or hold.
- **batch-3 MD** was a bonus exhaustive sweep of lower-cross-react tail candidates — fold any survivor into
  the backup pool; v7 is unaffected.

## Key invariants (see AGENTS.md)

- Novelty is a hard structural gate (<3/4 = DQ); portal-check the final set and drop <3/4 **before** locking.
  Local foldseek proxy is relative-rank only (no AFDB).
- Prove the switch with explicit two-state MD + triple-confirm (MD d_repel>0 + PROPKA dpKa<0 + xtb dE>0),
  on the exact submitted sequence, on a conserved target Lys. MD alone is insufficient (one candidate's
  strongest MD signal was a backwards switch caught only by PROPKA).
- Expression ranks on SAP + charged composition (Glu%) + GRAVY; ignore the ProtParam instability index.
- Commit source/docs/harness only; `runs/ results/ envs/ third_party/ tools/ weights/` are local-only.
