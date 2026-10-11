# v13 — strict-metric-respecting rebuild: aggregate-only/MD-refuted members replaced by MD-confirmed triggers

v13 answers two reviewer steers: (1) re-rank the core by the **strict** per-His fold-robust switch metric,
not the noisier aggregate linkage; (2) a design earns a slot only if **switch AND makeability both hold** —
more His sharpens the switch but hurts expression, so low-His designs that still switch are preferred.

## What the switch validation shows — and what it does NOT (read this first)

All switch numbers here are **computational proxies for a designed mechanism**, not a demonstration of
pH-dependent release. Be precise when ranking and writing methods:

- **MD `dSwitch` (+2.30 etc.) = a confirmed TRIGGER, not confirmed release.** It measures one thing: on
  protonation (pH7.4→6.0), does the His push *away* from its target TNF cation, beyond replicate noise,
  with a GLN-isostere control proving the motion is charge-caused and correctly directed. That is the
  first step of the mechanism. **Full dissociation is a slow (ms–s) event invisible to 3 ns MD**; the
  interface-contact metric stays ~0 over that window. A bigger dSwitch = a stronger trigger, not a
  guarantee the binder detaches in the lab (which is what the contest scores).
- **Aggregate Wyman linkage (−0.90/−0.93 etc.) = the closest *thermodynamic* estimate of the scored
  quantity** (the pH-dependent binding ΔΔG), but it is PROPKA-based, static and approximate — a
  prediction, not proof.
- **Strict per-His = fold-robust confirmation of linkage DIRECTION only.**

The defensible claim for v13's new members is: a **validated, causal, correctly-directed trigger** that
is *also* backed by a negative both-fold thermodynamic linkage — i.e. the designed release mechanism is
present and points the right way on two orthogonal computational measures. Real release is a wet-lab
question. "Switch-confident" below is shorthand for "trigger + linkage confirmed", never "release proven".

## What changed: exactly 2 swaps (v12 → v13)

| OUT | strict | MD | min x-react | → IN | His | MD dSwitch | min x-react | expr | axis-check |
|---|---|---|---|---|---|---|---|---|---|
| des_962_2_mh_9_32_35 | 0 | **NONRESPOND** | 0.846 | **des_495_6_mh_41_49** | **2** | **+2.30 (CONFIRMED)** | 0.91 | sol 0.95, 2-fold | **DOMINATES** |
| TNFpH06 | 0 | **NONRESPOND** | 0.648 | **des_1559_3_mh_48_75** | **2** | **+1.32 (CONFIRMED)** | 0.78 | sol 0.81, 2-fold | **DOMINATES** |

Both swaps raise the top axis (#1) AND the next (#2) with no axis traded down (per
[lesson] verify-all-axes-before-swap). des_495_6_mh_41_49 shows the **strongest MD trigger response of any
design tested** (dSwitch +2.30) and a strong both-fold linkage (−0.90/−0.93). Both swap-ins are **2-His**
(lowest makeability risk) yet trigger the hardest — the 4-His des_1060 variants both NONRESPOND (no trigger
fires), confirming the His-load/expression tradeoff empirically.

`runs/push/submission_v13.csv` — md5 **`4c682b80ee007c3957833daacdb91bb4`** — validator PASS (20 / 0-Cys /
unique / single_chain). TNFpH06_sh3 (not TNFpH06) stays: its 4th His makes it **MD-CONFIRMED (+0.88)**.

## Two-state MD His-causality matrix (the TRIGGER gate — not a release assay)

HIP vs HIE vs GLN-isostere, 4 replicas × 3 ns. CONFIRMED = the protonation TRIGGER fires: dSwitch > 2·SEM
AND > 0 AND dIso > 0. This gate confirms the charge-driven push, not full detachment (see caveat above).

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

- **#1 pH-selectivity:** mechanism-confirmed members (strict≥2 OR MD-trigger-causal) **9 → 11/20**, and the
  two new members carry the strongest MD trigger responses found (+2.30, +1.32) with correctly-directed
  both-fold linkage. The two removed members were MD-refuted (no trigger fired). This is confirmed-trigger
  evidence that the pH-release mechanism is present and correctly directed — not a demonstration of
  wet-lab release, which these short MD runs cannot show.
- **#2 cross-reactivity:** removes the portfolio's weakest cross-reactor (TNFpH06, 0.648) and des_962_2_mh
  (0.846); both replacements bind both species (0.78, 0.91). #2 floor rises.
- **expression (dominant risk):** both swap-ins are 2-His, Tier-A. **Independence structure of the
  makeability call (reviewer, 2026-10-11):** the consensus ANCHORS on ESMFold pLDDT (ESM-2 family), so the
  genuinely-independent hedges are what matter — and the swap-ins clear all of them: **Boltz-2 monomer
  pLDDT 0.92/0.89 (AF3-style diffusion folder, architecturally NON-ESM)**, **cross-architecture TM
  agreement 0.94/0.89 (ESMFold×Boltz — would flag an ESM-only hallucination)**, **Protein-Sol 0.95/0.81
  (experimental eSOL solubility, a non-model modality)**, and **0 developability liabilities (physics/
  rules)**. The makeability does NOT rest on ESM alone. ESM-C pseudo-LL naturalness was also run
  (des_1559_3 most-natural of the 20, des_495_6 in-distribution), but it is **within-ESM-family
  corroboration, NOT an independent vote** — ESMFold and ESM-C share a lineage and can fail together on an
  out-of-distribution design; it is reported as a sanity check, not counted toward independence. Honest
  caveat: the ANCHOR is still ESM-family + the Boltz/Protein-Sol/physics hedges; makeability remains a
  calibrated prediction, used only as a soft tie-breaker, never to override the #1/#2/novelty ordering.
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
