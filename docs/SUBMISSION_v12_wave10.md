# v12 — wave-10 multi-His switch scan on NOVEL cross-reactive folds → a submission that beats v11

## The campaign (genuinely new, not a re-tread)
Every prior switch campaign (waves 6–9) worked the SAME backbones (c19 / des_4098 family). Wave-10 attacks
an **untapped** region: the 14 **novel, cross-reactive de-novo bpair folds** (session-9 β-pairing campaign)
that were only ever **single-His** scanned or never scanned. A **multi-His separated-cation** switch (the
Baker SOTA form that gives the c19 flagship its 4 His) had never been tried on them.

`scripts/bpair_multihis_build.py` → for each fold, read the af2pred bound pose, find binder positions within
8 Å of a target cation, enumerate His **combinations** mutually separated ≥8 Å on distinct cations (singles,
pairs, triples, quads) → **147 variants (99 multi-His)**. Then: Boltz human linkage screen
(`wave10_screen.py`) → top-8 → **Chai fold + Boltz mouse + PROPKA+xTB both-fold triple-confirm**
(`wave10_confirm.sh`). All 8 GPUs, full triple-confirm rigor.

## What wave-10 established (two honest findings)
1. **Strong multi-His AGGREGATE switches DO install on novel cross-reactive folds.** Several novel folds gave
   3–4 His variants with both-fold Wyman linkage −0.8 to −1.0 (approaching c19's −1.274) AND excellent
   cross-reactivity (iptm 0.76–0.92 both species) — comparable to v11's shipped aggregate switches
   (TNFpH06_sh3 −1.15) but on NEW lineages.
2. **The STRICT per-His fold-robust multi-His switch does NOT reproduce** (c19 stays unique). On the independent
   Chai fold the xTB dE_switch sign flips for individual His (e.g. the 4-His des_2375_4 variant: 4/4 confirmed
   on Boltz, 0/4 on Chai), so the both-fold intersection is 0 — the same fold-dependence that makes c19's
   4 strict-fold-robust His irreducible ([[lesson-his-confirm-is-fold-dependent]]).
   - **The both-fold check caught an anti-switch**: `des_2375_4_mh_38_49_51` is Boltz −1.125 but **Chai +0.406**
     (wrong direction). Boltz-only would have shipped a backwards switch — exactly the des_287_5 failure mode.
   - **Multi-signal expression caught a hidden risk**: the strongest Boltz switch (`des_2375_4_mh_38_45_49_76`,
     4 His, linkage −1.007) has ESMFold monomer pLDDT **0.646** (below bar; the His load destabilised the fold)
     despite Boltz 0.89 — rejected on expression grounds.

## v12 = v11 with the two portfolio LAGGARDS replaced by two clean wave-10 switches
Validated on the identical pipeline (both folders, both-fold linkage, cross-species iptm, expression consensus):

| | OUT (v11 laggard) | IN (wave-10) | why it wins |
|---|---|---|---|
| swap 1 | **tnf_K90…sh4** — #2 laggard: iptm 0.594/0.543, switch −0.521 | **des_2106_2_mh_37_43_50_75** — iptm 0.920/0.886, switch Boltz −1.008 / Chai −0.904 | #1 ↑ (strongest both-fold switch found), #2 ↑↑; novelty-SAFE (des_2106_2 fold already portal 3/4 in v11) |
| swap 2 | **b1m12_2A67H_K112_r1** — #1 laggard: switch −0.208, solubility 0.585 (weakest), 7-member b-lineage | **des_962_2_mh_9_32_35** — switch Boltz −0.795 / Chai −0.565, iptm 0.945/0.846, **NEW de-novo lineage** | #1 ↑↑, expression ↑ (sol 0.83, 2-fold agree), diversity ↑ (new lineage; b-overload 7→6) |

**v12 beats v11 on the judged axes:**
- **#1 pH-selectivity:** two weak switches (−0.521, −0.208) → two strong **both-fold** aggregate switches
  (−1.008/−0.904, −0.795/−0.565), each confirmed correct-direction on BOTH Boltz and Chai. (All four are 0 on
  the *strict* per-His fold-robust metric — like v11's TNFpH06_sh3 — so c19 remains the unique strict flagship;
  the gain is on the aggregate ΔΔG_pH magnitude, the physical pH-selectivity.)
- **#2 cross-reactivity:** removes the portfolio's weakest cross-reactor (tnf_K90 0.54) and the b1m12 slot;
  both replacements bind both species strongly (min iptm 0.846, 0.886). Net #2 depth up.
- **expression (dominant wet-lab-zero risk):** six-signal consensus **improves** — Tier-A 12→**14**, Tier-C
  still **0**, expected expressing **17.1→17.2/20**, P(≥14)=0.985. Both new members are Tier A (2-fold agree,
  soluble, 0 liabilities).
- **novelty:** +1 novel de-novo lineage (des_962_2); des_2106_2_mh inherits portal 3/4.

`runs/push/submission_v12.csv` **md5 `e34ca96d0674cf85872f5b1976f396b5`** — validation gate PASS (20 / 0-Cys /
unique / single_chain).

## Upload ladder (USER-GATED — portal re-checks novelty at upload)
1. **v12** (`e34ca96d…`) — RECOMMENDED. Strongest on #1+#2+expression. **One contingency:** `des_962_2` is a
   de-novo lineage whose novelty is portal-UNCHECKED (expected ≥3/4, but de-novo ≠ guaranteed —
   [[lesson-portal-novelty-structural-calibrate-oracle]]). USER portal-check it before upload.
2. **v12_safe** (`0fbfe60f8c7ff0729bd44d8e669efde2`) — if `des_962_2` fails ≥3/4: the novelty-SAFE 1-swap
   (only tnf_K90→des_2106_2_mh, which inherits portal 3/4; keeps b1m12). Still beats v11 on #1+#2, zero novelty risk.
3. **v11** (`ef72d8e09f11d75d7aef959b4b18a2eb`) — portal-verified floor. v9 (`05a344…`) byte-locked fallback.

Nothing auto-uploaded. Deadline 2026-10-12 11:59 UTC.
