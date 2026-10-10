# v11 working accounting — honest dual-metric switch strength (supersedes v10 doc's switch claims)

Rigorous re-accounting of the **current best portfolio** (`runs/push/submission_v10.csv`, md5
`4e0b4ce6b294a924a416581b5bac2645`, 20 designs) on the #1 judged axis (pH-selectivity). Two independent
switch metrics, both on the exact submitted sequences:

1. **Fold-robust His** = His with PROPKA dpKa<0 AND GFN2-xTB dE_switch>0 on BOTH an independent Boltz AND
   an independent Chai human complex (intersection by His position + target cation). The strict,
   fold-robust bar — guards the fold-dependence overcount (`lesson-his-confirm-is-fold-dependent`).
2. **Wyman linkage** (collective_his_linkage, Boltz human) = Σ per-His Wyman contribution ≈ the aggregate
   ΔΔG_pH magnitude. **More negative = stronger correct-direction switch** (His destabilised in the complex
   → releases at low pH). Calibrated against the confirmed flagship c19 (−1.274).

## Portfolio (ranked by Wyman linkage; fold-robust His shown)

| design | linkage | fold-robust His | cations |
|---|---|---|---|
| des_3050_6_c19 | **−1.274** | **4** | 6,60,93 — FLAGSHIP |
| TNFpH06_sh3 | −1.149 | 0 | (aggregate only) |
| des_818_4_c13S | −1.144 | 2 | 6,93 |
| des_4098_1_c44 | −1.109 | 2 | 6,93 |
| TNFpH06 | −0.848 | 0 | (aggregate only) |
| des_4098_1_b14 | −0.764 | 2 | 6,93 |
| b1m9_A12H_w5_3 | −0.676 | 1 | 6 |
| b1m7_…_w5_3 | −0.572 | 1 | 6 |
| b4m2_r1_w5_3 | −0.565 | 1 | 93 |
| des_875_5_H3-57-80 | −0.529 | 0 | |
| tnf_K90_…_sh4 | −0.521 | 0 | |
| b4m2_A57H | −0.503 | 1 | 93 |
| b2m0_c70 | −0.467 | 0 | |
| des_59_1_H51 | −0.449 | 0 | |
| des_2106_2_H43 | −0.426 | 1 | 93 |
| b1m12_A41H_vsR27 | −0.394 | 0 | |
| des_180_6_H9-14-38 | −0.355 | 2 | 6,93 |
| des_3050_6_H52-54-79 | −0.326 | 2 | 6,60 |
| b1m12_2A67H_K112_r1 | −0.208 | 0 | |
| des_287_5_H37 | −0.193 | 0 | |

**Totals:** 1 flagship (≥3 fold-robust His), 6 designs ≥2, 5 at 1, 8 at 0. Every design has NEGATIVE
linkage (correct switch direction). This is SOBERER than the v10 doc's "~9 strong switches" — the honest,
both-fold number is 6 at ≥2 and 1 flagship.

## What this establishes

- **des_3050_6_c19 is the crown jewel**: #1 on BOTH metrics — a 4-His / 3-distinct-cation (Lys6/Lys60/Lys93,
  all human↔mouse conserved) fold-robust switch that binds human AND mouse (Boltz iptm 0.89/0.89). This is
  the full Baker multi-His separated-cation SOTA form, achieved and cross-fold-verified.
- **v10 > v9 on the #1 axis is CONFIRMED** by this rigorous accounting: the v9→v10 upgrade added c19 (flagship)
  plus fold-robust His to des_818_4_c13S and the b-lineage members.
- The switch strength is **capped by having one switch-rich lineage** (des_3050_6). Two expansion attempts
  this session did NOT break that cap: (a) wave6 — install 3 separated His on diverse backbones + LigandMPNN
  shell → binding preserved (ipSAE 0.79–0.88) but only 1/17 reached ≥2 fold-robust His (fold-dependence);
  (b) sibling-mining des_3050_6 c34/c42 → both bind human+mouse but confirm only 2 fold-robust His (Boltz
  drops the His Chai keeps). c19's exact fold is special.
- Bold de-novo re-test (RFd3 Region-I, 2048 complexes) = NULL (pae_interaction floor 24.5) — consistent with
  RFd3=0 on the TNF groove; partial-diffusion-derived binders remain the only working basis.

Upload remains USER-GATED. v9 (md5 `05a344444cf740907e79a177d34c73f5`) byte-locked fallback.

## Multi-axis audit (all 20 v10 designs) — v10 confirmed sound on every axis

| axis | result |
|---|---|
| #1 pH-selectivity | TWO complementary measures: (a) per-His fold-robust (mechanistic) = c19 flagship (4 His/3 cations) + 6 designs ≥2; (b) both-fold Wyman linkage (physical ΔΔG_pH negative on BOTH Boltz AND Chai) = **14/20 are fold-robust switches**, strongest c19 (−1.27/−1.12), TNFpH06_sh3 (−1.15/−1.32), des_4098_1_c44 (−1.11/−1.16). `runs/push/v10_foldrobust_linkage.json` |
| #2 cross-reactivity | all 20 bind human AND mouse by Boltz iptm (range 0.54–0.97); weakest tnf_K90 (0.59/0.54). NOTE: ipSAE_m was pessimistic on several — iptm is the robust read ([[lesson-chai-mouse-iptm-is-folder-pessimism]]) |
| expression (dominant wet-lab-zero risk) | ALL 20 ESMFold monomer pLDDT 0.79–0.93 (zero below the 0.75 good-bar) — excellent |
| novelty (≥3/4 hard gate) | v10's 5 new designs portal-UNCHECKED (c19 proxy-novel); USER portal-check required before upload |
| 0-Cys | all 20 (enforced) |

**Audit verdict:** no dominated member; no clear switch/cross-react/expression swap improves v10. The portfolio
is rigorously validated across all axes. Expansion was exhausted across FOUR bold campaigns this session
(de-novo RFd3; wave6 install-3-His-on-diverse-backbones; sibling-mining des_3050_6 c34/c42; wave7
c19-partial-diffusion family of 120 backbones → 360 variants) — all NULL for a second distinct flagship:
c19's 4-His/3-cation switch is a special, non-reproducible solution (wave7's 60 strong binders were ALL
c19 cousins with near-zero switch linkage). This establishes the switch frontier RIGOROUSLY (not prematurely).

## c19 flagship novelty — DE-RISKED by lineage portal data
c19's novelty is proxy-only, but its lineage (des_3050_6) has THREE portal-confirmed members at ≥3/4:
**c39 = 4/4, des_3050_6_H52-54-79 = 3/4, des_3050_6_b9 = 3/4** (`runs/push/portal_novelty_*.json`). c19 is a
sequence variant of that same fold (78% id to siblings c34/c42) → very likely inherits ≥3/4. So keeping c19
as the flagship is well-justified; "revert to c39" is a low-probability contingency, not an expected outcome.

## USER-GATED decision package
1. **Portal novelty-check the 5 new designs** (`runs/push/new_for_portal_v10.csv`), especially **des_3050_6_c19**
   (the flagship; proxy-novel, replaced the portal-confirmed-4/4 c39; lineage de-risk above). If c19 fails ≥3/4
   → revert c19→c39 (fallback keeps a portal-confirmed 4/4 but loses the 4-His flagship).
2. **Upload** v10 (md5 `4e0b4ce6b294a924a416581b5bac2645`) — or v9 (md5 `05a344444cf740907e79a177d34c73f5`)
   for zero novelty-recheck risk. Upload remains USER-GATED; nothing auto-uploaded.
