# v13 — strict-metric-respecting rebuild: aggregate-only/MD-refuted members replaced by MD-confirmed switches

v13 answers two reviewer steers: (1) re-rank the core by the **strict** per-His fold-robust switch metric,
not the noisier aggregate linkage; (2) a design earns a slot only if **switch AND makeability both hold** —
more His sharpens the switch but hurts expression, so low-His designs that still switch are preferred.

## What changed: exactly 2 swaps (v12 → v13)

| OUT | strict | MD | min x-react | → IN | His | MD dSwitch | min x-react | expr | axis-check |
|---|---|---|---|---|---|---|---|---|---|
| des_962_2_mh_9_32_35 | 0 | **NONRESPOND** | 0.846 | **des_495_6_mh_41_49** | **2** | **+2.30 (CONFIRMED)** | 0.91 | sol 0.95, 2-fold | **DOMINATES** |
| TNFpH06 | 0 | **NONRESPOND** | 0.648 | **des_1559_3_mh_48_75** | **2** | **+1.32 (CONFIRMED)** | 0.78 | sol 0.81, 2-fold | **DOMINATES** |

Both swaps raise the top axis (#1) AND the next (#2) with no axis traded down (per
[lesson] verify-all-axes-before-swap). des_495_6_mh_41_49 has the **strongest MD switch of any design
tested** (dSwitch +2.30). Both swap-ins are **2-His** (lowest makeability risk) yet switch the hardest —
the 4-His des_1060 variants both NONRESPOND, confirming the His-load/expression tradeoff empirically.

`runs/push/submission_v13.csv` — md5 **`4c682b80ee007c3957833daacdb91bb4`** — validator PASS (20 / 0-Cys /
unique / single_chain). TNFpH06_sh3 (not TNFpH06) stays: its 4th His makes it **MD-CONFIRMED (+0.88)**.

## Two-state MD His-causality matrix (the decisive selectivity gate)

HIP vs HIE vs GLN-isostere, 4 replicas × 3 ns. CONFIRMED = dSwitch > 2·SEM AND > 0 AND dIso > 0.

| design | role | dSwitch | verdict |
|---|---|---|---|
| des_495_6_mh_41_49 | **v13 swap-in** | +2.30 | **CONFIRMED** |
| des_2106_2_mh_37_43_50_75 | v12/v13 core | +1.40 | CONFIRMED |
| des_1559_3_mh_48_75 | **v13 swap-in** | +1.32 | **CONFIRMED** |
| des_1060_4_mh_46_74 | alternate (held) | +1.03 | CONFIRMED |
| TNFpH06_sh3 | v12/v13 core | +0.88 | CONFIRMED |
| TNFpH06 | **removed** | +0.29 | NONRESPOND |
| des_962_2_mh | **removed** | +0.40 | NONRESPOND |
| des_818_4_c13S | core (strict-2) | −0.01 | NONRESPOND (static strict-2 holds; MD-flagged) |
| des_1060_2 / des_1060_6 (4-His) | dropped | +0.33 / −0.04 | NONRESPOND |

The 7 v12 strict≥2 members (c19 + 5 others + des_818_4) are being MD-tested for the complete strict×MD
matrix (`runs/rebuild/coremd/`); des_818_4 already shows strict-2 can be MD-flat, so the rebuild's
MD-confirmed additions are genuinely orthogonal evidence, not a restatement of the static metric.

## Why v13 beats v12 on the judged axes

- **#1 pH-selectivity:** switch-confident (strict≥2 OR MD-causal) **9 → 11/20**, and the two new members
  carry the strongest MD switches found (+2.30, +1.32). The two removed members were MD-refuted.
- **#2 cross-reactivity:** removes the portfolio's weakest cross-reactor (TNFpH06, 0.648) and des_962_2_mh
  (0.846); both replacements bind both species (0.78, 0.91). #2 floor rises.
- **expression (dominant risk):** both swap-ins are 2-His, Tier-A (ESMFold 0.87/0.85, Boltz-mono 0.92/0.89,
  2-fold TM-agree, solubility 0.95/0.81, 0 liabilities). No expression regression.
- **novelty (hard gate):** structural proxy calibrated against known-passing v12 analogues (v12-pass maxTM
  0.807–0.955; des_962_2_mh passed at 0.955). Both swap-ins sit **below** that ceiling (0.896, 0.869) → in
  the ≥3/4 regime. **User must confirm on the Proteinbase portal before upload.**
- **diversity:** +2 distinct new lineages (des_495, des_1559); caps respected (des_1060 held at 0 in final).

## Upload ladder (USER-GATED — nothing auto-uploaded)

1. **v13** (`4c682b80…`) — RECOMMENDED. Requires the user to run the portal novelty check on the 2 new
   members (des_495_6_mh_41_49, des_1559_3_mh_48_75); proxy predicts ≥3/4 for both.
2. **v12** (`e34ca96d…`) — fully portal-verified floor; the fallback if either new member fails the portal.
3. v11 (`ef72d8e0…`) / v9 (`05a34444…`) — prior byte-locked fallbacks.

Deadline 2026-10-12 11:59 UTC.
