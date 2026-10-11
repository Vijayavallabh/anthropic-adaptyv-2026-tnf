# Complete 20-design two-state-MD switch-TRIGGER matrix (v14)

Every member of `submission_v14.csv` (md5 `7994cf46…`) scored on the explicit two-state MD
His-causality gate (HIP/HIE/GLN isostere, 4 reps × 3 ns). This is the full matrix — no member
is un-MD'd anymore.

**What the number is (honest framing).** `dSwitch = d(His–Lys | HIP) − d(His–Lys | HIE)` measures the
charge-driven **trigger**: does protonating the His at pH 6 push it away from a conserved target cation
(Coulomb repulsion)? `dIso = HIP − GLN` is the causal control (charge vs neutral isostere). A positive,
>2·SEM, dIso-consistent ladder CONFIRMS the trigger. It does **not** demonstrate wet-lab release —
full dissociation is ms-scale and invisible in 3 ns (`d_local` ≈ 0 throughout). "CONFIRMED" = the
first mechanistic step is real and correctly directed, nothing more. See `docs/SUBMISSION_v13.md`.

## Matrix (sorted by dSwitch)

| design | HIP | HIE | GLN | dSwitch | dIso | 2·SEM | verdict |
|---|---|---|---|---|---|---|---|
| des_495_6_mh_41_49 | 12.55 | 10.25 | 11.96 | **+2.30** | 0.59 | 0.73 | CONFIRMED |
| b2m0_c70 | 7.99 | 5.77 | 5.05 | **+2.22** | 2.94 | 0.81 | CONFIRMED |
| des_2106_2_H43 | 5.42 | 3.48 | 4.02 | **+1.94** | 1.40 | 0.23 | CONFIRMED |
| des_875_5_H3-57-80 | 9.85 | 8.41 | 9.70 | **+1.44** | 0.15 | 0.89 | CONFIRMED |
| des_2106_2_mh_37_43_50_75 | 10.39 | 8.99 | 8.53 | **+1.40** | 1.86 | 0.34 | CONFIRMED |
| des_1559_3_mh_48_75 | 5.99 | 4.67 | 3.74 | **+1.32** | 2.25 | 0.39 | CONFIRMED |
| des_4098_1_b14 | 6.65 | 5.57 | 4.22 | **+1.08** | 2.43 | 0.31 | CONFIRMED |
| des_1060_4_mh_46_74 | 11.03 | 10.00 | 8.57 | **+1.03** | 2.45 | 0.30 | CONFIRMED |
| des_180_6_H9-14-38 | 5.90 | 5.01 | 4.41 | **+0.89** | 1.49 | 0.19 | CONFIRMED |
| b1m9_A12H_w5_3 | 10.19 | 9.40 | 8.50 | **+0.79** | 1.70 | 0.38 | CONFIRMED |
| des_4098_1_c44_w6_3 | 7.72 | 7.15 | 6.27 | **+0.57** | 1.45 | 0.37 | CONFIRMED |
| des_4098_1_c44 | 5.88 | 5.36 | 3.38 | **+0.51** | 2.50 | 0.28 | CONFIRMED |
| b4m2_r1_w5_3 | 6.70 | 6.30 | 5.34 | **+0.40** | 1.37 | 0.26 | CONFIRMED |
| des_3050_6_H52-54-79 | 7.70 | 7.38 | 6.22 | +0.31 | 1.47 | 0.41 | NONRESPOND |
| des_3050_6_c19 | 8.05 | 7.80 | 7.39 | +0.25 | 0.66 | 0.41 | NONRESPOND |
| b4m2_A57H | 6.83 | 6.61 | 5.02 | +0.22 | 1.81 | 0.39 | NONRESPOND |
| b1m12_A41H_vsR27 | 7.91 | 7.84 | 4.96 | +0.07 | 2.96 | 1.70 | NONRESPOND |
| des_818_4_c13S | 10.70 | 10.71 | 10.46 | −0.01 | 0.24 | 0.41 | NONRESPOND |
| des_59_1_H51 | 7.14 | 7.77 | 7.86 | **−0.62** | −0.72 | 0.26 | NONRESPOND (wrong sign) |
| b1m7_2A57H_K183_r1_w5_3 | 6.80 | 7.89 | 6.00 | **−1.08** | 0.80 | 0.65 | NONRESPOND (wrong sign) |

**13 CONFIRMED · 7 NONRESPOND.**

## Decision implications
- The **13 CONFIRMED** designs carry the #1 (pH-selectivity) axis with a causal, correctly-directed
  trigger. They are the keep-core.
- The **7 NONRESPOND** are the portfolio's weak links on the dominant axis. Two (`des_59_1_H51`,
  `b1m7_2A57H…`) have a **wrong-sign** ladder — the fast MD shows no repulsion trigger. These are the
  priority swap targets for the v15 regeneration campaign (partial-diffused + β-paired new switchers),
  but only a v15 candidate that is itself MD-trigger-CONFIRMED **and** holds cross-react ≥0.70 **and**
  novelty ≥3/4 **and** expression earns a swap — never a switch score alone.
- NONRESPOND is a deprioritiser, not an auto-drop: a NONRESPOND design with a strong static both-fold
  switch (PROPKA+xTB) + 4/4 novelty (e.g. `des_818_4_c13S`) still has independent merit; MD just can't
  confirm its trigger in 3 ns.
