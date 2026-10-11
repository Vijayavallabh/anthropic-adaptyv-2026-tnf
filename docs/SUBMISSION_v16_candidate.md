# v16 candidate — de-novo NEW-shape switchers (the long shot), NOVELTY-GATED

`runs/push/submission_v16_candidate.csv` — md5 **`a524d77fa7a63ef33b28346ce3a7d87d`** — validator **PASS**
(20 / 0-Cys / unique / single_chain). **Upload + novelty are USER-gated; v14 (`7994cf46`) stays the locked floor.**

## What the long shot produced
User directive: "go with the long shot" — genuinely NEW shapes for the 7 weak (MD-NONRESPOND) v14 slots,
not partial-diffusion cousins of the 4 already-confirmed lineages (which would collapse diversity).

De-novo (β-pairing RFdiffusion) → 480 new shapes → 20 binders (15 shapes) → **3 cross-reactive new shapes**
(des_138_7, des_4082_8, des_2044_8) → partial-diffusion amplification of those 3 → 576 cousins → 548 binders
→ 92 cross-reactive → 75 strong multi-His switches → **5 MD-trigger-CONFIRMED new-shape switchers** on 2
new lineages:

| swap-in (new shape) | lineage | dSwitch (MD) | x-react | ESMFold pLDDT | Protein-Sol | seq-novelty |
|---|---|---|---|---|---|---|
| des_138_7_T10_pd_24_4_mh_49_58_77 | des_138_7 | **+1.11 CONF** | 0.87 | 0.87 | 0.71 | no PDB hit |
| des_138_7_T10_pd_39_4_mh_49_77 | des_138_7 | **+1.09 CONF** | 0.87 | 0.90 | 0.76 | no PDB hit |
| des_4082_8_T10_pd_10_4_mh_51_79 | des_4082_8 | **+0.95 CONF** | 0.90 | 0.81 | 0.64 | no PDB hit |
| des_4082_8_T10_pd_9_1_mh_43_55_57_79 | des_4082_8 | **+0.46 CONF** | 0.93 | 0.86 | 0.52 | no PDB hit |
| des_138_7_T20_pd_9_2_mh_36_53_77 (backup) | des_138_7 | +0.63 CONF | 0.71 | 0.87 | 0.87 | no PDB hit |

The lead parent des_138_7_mh_47_53_80 was MD-**NONRESPOND** (wrong-sign −0.78) despite strong static
linkage −1.05 — the static-linkage≠MD-trigger discipline holding again; only backbone-diversified cousins
passed MD.

## The v16 swap (contingent on novelty)
Swap **OUT** the 4 weakest v14 members by the two-state MD trigger:
`b1m7_2A57H_K183_r1_w5_3` (−1.08, wrong-sign), `des_59_1_H51` (−0.62, wrong-sign),
`b1m12_A41H_vsR27` (+0.07, 2·SEM 1.70 — unreliable), `des_3050_6_c19` (+0.25, flat; des_3050 keeps H52-54-79).
Swap **IN** the 4 strongest new-shape MD-confirmed switchers (≤2 per new lineage).
Net: **+4 MD-confirmed on the #1 axis, +2 genuinely new lineages, 0 Cys, validator PASS.**

## NOVELTY gate: PASSED (portal-verified 2026-10-11)
User ran the Proteinbase/Adaptyv novelty check on all 5 new-shape switchers → **every one 3/4 (≥3/4 = PASS)**.
As predicted from the local proxy: sequence-novel (no PDB-seqres hit) + structural-capped (max TM 0.86–0.93,
scaffold-guided β-pairing) → 3/4 — the same profile as the v14 bpair members (des_1060_4 at 4/4). All 4 v16
swap-ins clear the hard novelty gate.

## Verdict: v16 is a verified strict improvement over v14
- **#1 pH-selectivity:** 4 MD-NONRESPOND members (2 wrong-sign) → 4 MD-trigger-CONFIRMED new-shape
  switchers. Portfolio two-state-MD trigger goes **13/20 → 17/20 CONFIRMED** (3 NONRESPOND remain, kept for
  their other merits/diversity).
- **#2 cross-reactivity:** swap-ins min-iptm **0.87–0.91** (strong; no higher axis traded down — swap-outs
  were NONRESPOND on the top axis).
- **novelty:** all swap-ins 3/4 (gate passed). **diversity:** +2 genuinely new lineages (des_138_7,
  des_4082_8); the over-represented des_3050 trimmed to 1.
- **expression:** ESMFold pLDDT 0.81–0.90, Protein-Sol pass, sequence-novel.
- validator PASS, 0-Cys, 20 unique single_chain.

**v16 (`a524d77f`) is the recommended best. Upload remains USER-gated; v14 (`7994cf46`) stays the fallback.**
A 3rd new lineage from the still-running bpair generation could feed a v17 that also displaces the 3
remaining NONRESPOND members.
