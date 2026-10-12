# Submission v18 candidate — 4th de-novo lineage + strongest switcher (USER-PORTAL-GATED)

Updated 2026-10-12 (SESSION 11, continued). Deadline **2026-10-11 23:59 AoE = 2026-10-12 11:59 UTC**.

## TL;DR

v18 = **v16 with two swaps**, both MD-NONRESPOND members replaced by **MD-CONFIRMED new-lineage
cross-reactive switchers** harvested from the bpair/partial-diffusion switch sweep:

| slot (v16 NONRESPOND) | swapped in | dSwitch | lineage |
|---|---|---|---|
| `des_3050_6_H52-54-79` (+0.31, below thr) | **`des_6126_3_T20_pd_45_3_mh_9_56_78`** | **+1.36** | **des_6126 (NEW 4th de-novo lineage)** |
| `b4m2_A57H` (+0.22, below thr) | **`des_6025_4_T10_pd_77_4_mh_52_74`** | **+0.93** | des_6025 (3rd de-novo lineage) |

The third v16 NONRESPOND member, **`des_818_4_c13S` (the only 4/4-novel design), is KEPT** — deliberately
not swapped for a 3rd available des_6025 candidate (a cousin), to preserve lineage diversity + the 4/4
novelty showcase (per "diverse not cousins"; the extra cousin would add only a weak +0.68 switcher while
losing a distinct lineage and the lone 4/4).

- Portfolio two-state-MD trigger **17/20 → 19/20 CONFIRMED** (one NONRESPOND kept on purpose: the 4/4 `des_818_4`).
- De-novo new-lineage switchers **2 → 4** (`des_138_7`, `des_4082_8`, `des_6025`, **`des_6126`**).
- des_6126 is the **strongest switcher in the whole portfolio** (d_repel +1.36, all 4 reps positive).
- **v18 CSV:** `runs/push/submission_v18.csv` (md5 **`524722158467273c0e9f57f5402650ed`**), validator **PASS**
  (20 designs, 0 errors, 0 warnings, 0-Cys, unique names+seqs, all `single_chain`).
- **Fallbacks (unchanged):** v17 `9357b6905ed81179da5137a16dd9dd5a`; v16 `a524d77fa7a63ef33b28346ce3a7d87d` (SAFE default).

## Swap-in 1 (flagship): des_6126_3_T20_pd_45_3_mh_9_56_78

`MTEEESRKHRAEVAAILEAAAEKALEEVGDEETAEKLRAFAEAVRAGTSGAGVSIHGGGALVPVFESLGEALAALNKHLA` (80 aa, 0 Cys, His 9/56/78)

| Axis | Metric | Value | Verdict |
|---|---|---|---|
| **#1 pH-switch** | two-state MD d_repel (HIP−HIE), 4 reps × 3 ns | **+1.36** (2·SEM ≈ 0.40) | **CONFIRMED (strongest)** |
|  | per-replica HIP−HIE | +2.04 / +1.02 / +0.66 / +1.72 | **4/4 reps positive** (not an outlier) |
|  | isostere control (HIP−GLN) | +0.39 | positive (His-causal ablation holds) |
|  | `selective_release` | `true` | — |
|  | binder RMSD / interface frac | ~2.0 Å / 0.84 | complex MD-stable across states |
| **#2 cross-react** | human iptm (Boltz-2) | **0.970** | — |
|  | mouse iptm (Boltz-2) | **0.911** | — |
|  | min(h,m) | **0.911** | **PASS** (≥0.70) |
| **expression** | Protein-Sol scaled solubility | 0.82 | PASS (≥0.50) |
| **novelty seq** | mmseqs vs PDB seqres | **NO PDB hit (novel)** | passes seq criterion |
| **novelty struct** | foldseek max-TM vs PDB | **0.948** | 3/4 (helical-bundle cap; see risk) |

## Swap-in 2: des_6025_4_T10_pd_77_4_mh_52_74

`SEETERLIRELNQILLEEMKETVKKLKELGVEVELNVSENGGVSLIASSESHELALEASQLIKETVRRIEQLLHKYVAER` (80 aa, 0 Cys, His 52/74)

| Axis | Metric | Value | Verdict |
|---|---|---|---|
| **#1 pH-switch** | two-state MD d_repel (HIP−HIE), 4 reps × 3 ns | **+0.93** (2·SEM ≈ 0.36) | **CONFIRMED** |
|  | isostere control (HIP−GLN) | +2.18 | positive |
| **#2 cross-react** | human / mouse iptm (Boltz-2) | 0.950 / 0.936 → **min 0.936** | **PASS** |
| **expression** | Protein-Sol scaled solubility | 0.92 | PASS |
| **novelty** | mmseqs seq / foldseek struct TM | NO PDB hit / **0.870** | 3/4 |

Note: this is a **stronger** des_6025 switcher than v17's member (`mh_34_50_52_74`, +0.372) — v18 upgrades the
des_6025 representative while using an independent variant.

## The one open risk: portal structural-novelty (USER-GATED, same as v17/v16)

Both new members are seq-novel (no PDB hit) de-novo helical bundles with foldseek structural max-TM
0.87–0.95 — the same regime as every de-novo member already in v16/v17 (e.g. des_6025 0.945). By the local
calibration (seq-novel + helical ⇒ **3/4**, which passes the ≥3/4 hard gate), both should pass. But the
portal's structural-novelty criterion cannot be confirmed locally with certainty — **USER must portal-check
the two new members (and ideally all 20) before upload.** If either reads <3/4 on the portal, fall back to
v17 (which has only one portal-unverified new member) or v16 (fully verified).

## Recommendation

- **Safest, zero portal checks:** upload **v16** (`a524d77f`).
- **One portal check (des_6025 member):** **v17** (`9357b690`).
- **Three portal checks (the two new v18 members; confirm ≥3/4):** **v18** (`524722158467273c0e9f57f5402650ed`)
  — strongest portfolio on the #1 (pH-selectivity) axis: 19/20 MD-confirmed switchers, a new 4th de-novo
  lineage, the single strongest switcher (+1.36), all cross-reactive (min-iptm 0.91–0.95) and soluble.

**Upload remains USER-GATED. Nothing has been uploaded.**
