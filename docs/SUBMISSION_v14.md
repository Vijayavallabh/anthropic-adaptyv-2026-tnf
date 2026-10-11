# v14 — v13 + one axis-complete swap, FULLY PORTAL-VERIFIED (new recommended best)

v14 is v13 with a single dominating swap that fixes the portfolio's clearest weak link (the worst
cross-reactor), surfaced by the brutal devil's-advocate pass. Now portal-verified end-to-end.

## The swap (v13 → v14)

| OUT | why out | → IN | result |
|---|---|---|---|
| **TNFpH06_sh3** | MD-confirmed switch (+0.88) but **cross-react 0.62 — worst in the portfolio** (#2 axis) | **des_1060_4_mh_46_74** | **DOMINATES** |

des_1060_4_mh_46_74: MD trigger **+1.03 (CONFIRMED)**, cross-react **0.84**, 3-His, Tier-A expression
(ESMFold 0.87 / Boltz-mono 0.91 / cross-TM 0.94 / Protein-Sol 0.84 / 0 liabilities), **novelty 4/4**.

Every judged axis holds or improves, nothing traded down:
- **#1 pH-selectivity:** switch +0.88 → **+1.03** (both MD-trigger-confirmed).
- **#2 cross-reactivity:** **0.62 → 0.84** (+0.22; removes the portfolio's weakest cross-reactor).
- **novelty:** that slot **3/4 → 4/4** (des_1060_4 is 4/4; the portfolio now has 5 designs at 4/4).
- **expression:** Tier-A both; His load 4 → 3 (lower = better makeability).
- **diversity:** removes the TNFpH06 lineage, adds des_1060 (new to the portfolio; cap respected).

## Portal verification (2026-10-11)

User ran the Proteinbase novelty check on all 20 v14 designs → **every design ≥3/4**
(`runs/push/portal_novelty_v14.tsv`); the new member des_1060_4_mh_46_74 passed at **4/4**. 5 designs at
4/4: des_818_4, des_875_5, des_1060_4_mh_46_74, b1m12, des_2106_2_mh.

`runs/push/submission_v14.csv` — md5 **`7994cf46f0a55921a986cf73cee5e13d`** — validator PASS (20 / 0-Cys /
unique / single_chain).

## Upload ladder (USER-GATED — nothing auto-uploaded)

1. **v14** (`7994cf46…`) — RECOMMENDED, fully portal-verified, strictly dominates v13 (#2 + novelty up,
   #1 edges up, expression held). Ready for upload.
2. **v13** (`4c682b80…`) — fully portal-verified, retained as the conservative fallback.
3. v12 (`e34ca96d…`) / v11 (`ef72d8e0…`) / v9 (`05a34444…`) — prior verified fallbacks.

## Context / honest caveats carried forward
- Switch numbers are confirmed **triggers** (charge-driven His–cation repulsion, causal + correctly
  directed), not demonstrated wet-lab release — see `docs/SUBMISSION_v13.md`.
- Makeability anchors on the ESM family + the independent Boltz-2 / Protein-Sol / physics legs; a soft
  tie-breaker, not a guarantee.
- SOTA sweep (Zatom-2/Nucleus, ArcRefine, Proton-PottsMPNN) found no runnable ceiling-breaker that
  survived calibration on our data; stack is at the field frontier.
- The full 20-design two-state-MD switch-trigger matrix is still completing (rigor; may upgrade/confirm
  individual members but does not change the v14 composition).

Deadline 2026-10-12 11:59 UTC.
