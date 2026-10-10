# Submission v10 — candidate (fold-robust switch upgrades), USER-PORTAL-GATED

**v10** `runs/push/submission_v10.csv` md5 `4e0b4ce6b294a924a416581b5bac2645` — 20 designs, 0-Cys, unique,
`single_chain`, len 73–119. **v9** `runs/push/submission_v9.csv` md5 `05a344444cf740907e79a177d34c73f5` remains
the byte-locked, fully-validated FALLBACK.

v10 = v9 with **5 same-lineage fold-robust switch upgrades** (diversity preserved: still 16 lineages). It is
stronger than v9 on the #1 axis (pH-selectivity): it roughly doubles the number of strong fold-robust multi-His
switches (~5 → ~9) by upgrading v9's fold-robust-weakest members and converting diverse v7-floor binders into
switches. **The 5 new designs need a USER portal novelty check before v10 can be uploaded** (each is an MPNN
reseq on a portal-confirmed parent fold, so novelty should inherit — proxy TM confirms the parents' band — but
the portal is the sole authority). If any new design returns <3/4, revert that one swap (fallback to the v9
member) or keep v9 entirely.

## The 5 upgrades (`runs/push/v10_upgrades.tsv`; `new_for_portal_v10.csv`)

| v9 member (out) | v10 design (in) | lineage | fold-robust His (Boltz∩Chai) | gain |
|---|---|---|---|---|
| des_818_4_H16-74 | des_818_4_c13S | des_818_4 | 2 (Lys6,Lys93) | 1→2 robust His; dHM 0.105→~0.008; SAP 5.4→0.9 |
| des_3050_6_c39 | des_3050_6_c19 | des_3050_6 | 3 (Lys60,Lys93) | 2→3 robust His, 1→2 distinct cations |
| b1m7_2A57H_K183_r1 | b1m7_…_w5_3 | b1m7 | 2 (Lys6,Lys93) | floor→2 robust His; bind 0.879/0.878 (Δ0.001); SAP 1.1 |
| b4m2_r1 | b4m2_r1_w5_3 | b4m2 | 2 (Lys6,Lys93) | floor→2 robust His/2 cations |
| b1m9_A12H | b1m9_A12H_w5_3 | b1m9 | 1 (Lys6) | floor(expr-weak)→1 robust His; SAP 6.7→4.7 |

All 5: bind both species on BOTH co-folders (Boltz+Chai), balanced cross-species delta, expression pLDDT ≥0.83,
0-Cys, novelty proxy in the parent's band. fold-robust His = confirmed (PROPKA dpKa<0 AND xTB dE>0) on BOTH the
Boltz AND the Chai structure — the strict cross-fold bar (a single-structure count overcounts; see §2b of the v9
evidence doc).

## How they were found (Campaign F — in response to the "keep improving" directive)

- **F1 (mine under the fold-robust bar):** Chai-folded + Chai-triple-confirmed the 14 non-v9 wave4b designs. All
  bind on Chai; des_818_4_c13 (→c13S after Cys→Ser) and des_3050_6_c19 emerged as fold-robust upgrades.
- **F2 (multi-His on DIVERSE backbones = new lineages):** the proven wave3b method (His-lock at the strong
  anchors Lys93+Lys6 + LigandMPNN local-shell redesign) extended to 11 diverse v7-floor backbones that contact
  both anchors → 176 variants → 159 bind both species (Boltz) → top-2/lineage triple-confirmed on Boltz+Chai →
  b1m7_w5_3, b4m2_r1_w5_3 (fold-robust 2 His), b1m9_w5_3 (1). `scripts/wave5_prep.py`, `scripts/triple_confirm_mut.py`.

## Recommendation

Upload v10 after portal-checking the 5 new designs (`new_for_portal_v10.csv`). If the user prefers zero
novelty-recheck risk, v9 is the fully-validated fallback. The c39→c19 swap is the only one that replaces a
portal-confirmed 4/4 (c39) with an unconfirmed design (c19, proxy-more-novel) — portal-check it first; if it does
not hold 4/4, keep c39. **Upload remains USER-GATED; nothing is auto-uploaded.**
