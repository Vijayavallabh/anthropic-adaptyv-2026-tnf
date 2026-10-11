# v12 re-ranked by the STRICT per-His fold-robust switch metric

The 12-design brutal cull ranked mostly by the **aggregate** Wyman linkage (collective ΔΔG_pH summed
over all His). That measure is noisier than the **strict** per-His fold-robust check, and four of the
eight "core" designs scored zero on the strict check while two strict-positive designs were cut. This
re-rank redoes keep/cut on the strict check so the strong core does not rest on the shakier measure.

## The two metrics

- **STRICT (fold-robust per-His):** a His counts only if PROPKA ΔpKa<0 AND GFN2-xTB ΔE_switch>0 on the
  **submitted-sequence complex folded by BOTH Boltz-2 AND Chai-1** — the fold-robust count is the
  intersection of confirmed His across the two independent folders. This is the rigorous mechanistic
  claim (`runs/push/strict_rerank_v12.json`).
- **MD (two-state His-causality):** HIP vs HIE vs GLN-isostere control, 4 replicas × 3 ns each. A design
  is MD-confirmed if d_switch = His–Lys(HIP) − His–Lys(HIE) > 2·SEM AND > 0 AND the GLN control is
  passed. This equilibrates geometry and bypasses the raw-folder xTB sign-flip; it is an independent
  line of evidence, not a re-run of the strict check.
- **AGGREGATE (what the cull over-weighted):** collective linkage, both folds negative = correct
  direction. Useful as a screen; insufficient alone to claim a confirmed switch.

## Strict fold-robust result across all 20 v12 members

| design | strict (B∩C) | MD | verdict |
|---|---|---|---|
| des_3050_6_c19 | **4** (H5,49,54,79) | — | **core** — the unique strict flagship |
| des_4098_1_b14 | **2** (H6,74) | — | **core** |
| des_4098_1_c44 | **2** (H6,74) | — | **core** |
| des_4098_1_c44_w6_3 | **2** (H6,10) | — | **core** |
| des_180_6_H9-14-38 | **2** (H9,38) | — | **core** (expression caveat: LIAB) |
| des_818_4_c13S | **2** (H38,74) | pending | **core** — strict-2 holds on direct records despite earlier Chai-aggregate +0.09 |
| des_3050_6_H52-54-79 | **2** (H54,79) | — | **core** |
| des_2106_2_mh_37_43_50_75 | 0 | **CONFIRMED** (dSw +1.40) | **core** — MD-causal |
| des_2106_2_H43 | 1 (H43) | — | partial (single fold-robust His) |
| b1m7_2A57H_K183_r1_w5_3 | 1 (H3) | — | partial |
| b4m2_r1_w5_3 | 1 (H37) | — | partial |
| b4m2_A57H | 1 (H37) | — | partial |
| b1m9_A12H_w5_3 | 1 (H55) | — | partial |
| des_875_5_H3-57-80 | 0 | — | binder / diversity |
| des_59_1_H51 | 0 | — | binder / diversity |
| b1m12_A41H_vsR27 | 0 | — | binder / diversity |
| b2m0_c70 | 0 | — | binder / diversity |
| **TNFpH06** | **0** | pending | **aggregate-only** — user's concern |
| **TNFpH06_sh3** | **0** | pending | **aggregate-only** — user's concern |
| **des_962_2_mh_9_32_35** | **0** | **NONRESPOND** (dSw +0.40<2·SEM) | **aggregate-only; MD-refuted** — swap candidate |

## What the re-rank establishes

- **Switch-confident core = 8/20** on the rigorous bar: 7 designs strict≥2 **plus** des_2106_2_mh
  (MD-causal). This is the honest, conservative count — computed from the direct submitted-sequence
  Boltz/Chai complexes, not from more favorable re-folds.
- **5 designs are strict-1** (one fold-robust His) — a real but single switch element; they also carry
  cross-reactivity and expression value.
- **The user's concern is correct for exactly three members.** TNFpH06, TNFpH06_sh3 and des_962_2_mh are
  strict-0; des_962_2_mh is additionally **MD-NONRESPOND**. The v12 #1 gain claimed on des_962_2_mh's
  aggregate linkage does not survive the rigorous check.

## Decision rule (applied from here)

A design is **core (switch-confident)** iff **strict≥2 OR MD-causal** — ideally both. Aggregate linkage
alone is no longer sufficient to claim a confirmed switch. Under this rule:

- **Keep as core:** the 7 strict≥2 + des_2106_2_mh (MD). (keep-decision MD on TNFpH06, TNFpH06_sh3,
  des_818_4 is running; des_818_4 is already core on strict-2, its MD resolves the old aggregate conflict.)
- **Swap candidate:** des_962_2_mh (MD-refuted) — replace with a rebuild design that is MD-causal AND
  makeable AND novel AND a distinct lineage.
- **TNFpH06 / TNFpH06_sh3:** survive as core only if the keep-decision MD confirms them; otherwise they
  drop to the aggregate/diversity tier (still valid portfolio members, not counted as confirmed switches).

## Fix in flight (rebuild)

The rebuild re-mines multi-His switches on the novel cross-reactive backbones and runs the full gauntlet
— both-fold sign + cross-react + **expression hard co-gate** + two-state MD. It is sourcing MD-causal,
well-expressing, distinct-lineage replacements for the aggregate-only members, so the next submission's
#1 gain rests on confirmed physics, not on the aggregate number alone. v12 (`e34ca96d…`) remains the
standing submission and fallback; upload stays user-gated.
