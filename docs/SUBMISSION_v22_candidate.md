# v22 candidate — pragmatic Pareto-improvement over v19 (LOCAL-ONLY, upload USER-GATED)

Built from the rigorous v21-vs-v19 gauntlet: v19's only real weakness was 2 MD-failing switchers.
v22 swaps exactly those 2 for 2 strong direction-confirmed new-lineage switchers -> >= v19 on every axis,
strictly > on switch-DIRECTION + diversity, TIE on binding/cross-react/expression, 2 portal re-checks.
NOTE: v22 does NOT dominate v21 on the #1 axis — v21's both-fold His-count depth (8 members) is higher expected FC.

## STATUS
```json
{
 "n": 20,
 "md5": "abaec06790766f25282133d929083187",
 "validator": "PASS rows=20 errors=0",
 "cys": 0,
 "upload": "USER-GATED"
}
```

## THESIS
```json
"v22 = v19 with ONLY its 2 MD-FAILING switchers replaced by 2 strongest direction-confirmed new-lineage switchers. Constructed from the rigorous v21-vs-v19 gauntlet, which showed v19's weakness is exactly 2 non-switching members (not binding/cross-react/expression, which tie). PARETO-DOMINATES v19 with only 2 portal re-checks (vs v21's 13); keeps 18 v19 members unchanged (already portal-confirmed >=3/4)."
```

## SWAP
```json
{
 "cut_v19_MD_fails": {
  "b1m9_A12H_w5_3": "MD d_repel -1.307 (WRONG-sign switch), bind pae 4.87, sol 0.78",
  "des_138_7_T10_pd_24_4_mh_49_58_77": "MD d_repel -0.097 (non-causal switch), bind pae 6.01, sol 0.71 (also has cousin des_138_7_pd_39 still in set)"
 },
 "added_direction_confirmed_new_lineage": {
  "des_36133_g4_sa18": "3-cation [6,60,93], MD d_repel 1.107 >> Gln -1.190 (gap +2.30, strongest), bind pae 4.48, xreact 92.7%, sol 0.77, esm 0.90, nov_TM 0.842 (fresh de-novo campaign)",
  "des_732_1_r3_sa15": "3-cation, MD d_repel 0.568 >> Gln -0.322, bind pae 4.53, xreact 93%, sol 0.72, nov_TM 0.842 (regen)"
 }
}
```

## PARETO_DOMINANCE_vs_v19
```json
{
 "switch_DIRECTION_#1": "BETTER (robust, DIRECTION not magnitude) \u2014 removed 2 members that FAIL the direction gate (d_repel -1.307 wrong-sign, -0.097 non-causal) and added 2 His-causal (+1.107, +0.568). This is a direction improvement (switching vs not-switching), which is robust; it is NOT a magnitude/sharpness claim (magnitude is unrankable).",
 "switch_FC_ordinal_#1": "NOTE: on expected fold-change via both-fold His-count, v21 > v22 (v21 has 8 both-fold 2-3cat members vs v22's 2). v22's 18 v19-derived members are direction-OK (MD) but NOT both-fold geometry-confirmed. So v22 fixes v19's 2 DEAD switchers but does not match v21's high-FC-confidence depth. See SWITCH_FC_ORDINAL.",
 "binding": "TIE \u2014 all members pass the af2ig floor; pae_int differences (4.5 vs 4.9/6.0) are within noise and do NOT predict wet-lab binding. (Earlier 'strictly better on binding' was an overclaim; corrected.)",
 "cross_react": "TIE \u2014 added 92.7/93% vs cut ~92-93%.",
 "expression": "TIE \u2014 added sol 0.77/0.72 vs cut 0.78/0.71.",
 "diversity": "BETTER \u2014 +2 new lineages, removed a des_138_7 cousin.",
 "novelty": "= IF the 2 added pass portal re-check (nov TM 0.842 helical -> expect 3/4). 18 members already portal-confirmed.",
 "VERDICT": "v22 Pareto-improves v19: >= on EVERY axis, STRICTLY > on switch-DIRECTION (fixes 2 dead switchers) + diversity; TIE on binding/cross-react/expression; = novelty pending 2 re-checks. It does NOT dominate v21 on the #1 axis (v21 has more both-fold multi-His = higher expected FC). v22 = pragmatic strict-improvement over the confirmed baseline at minimal (2) re-check risk."
}
```

## vs_v21
```json
"v22 keeps 18 confirmed-novelty members (vs v21's 4) and needs 2 re-checks (vs 13); it forgoes v21's extra diversity (7 new lineages) for far lower re-check risk + a strict improvement over the portal-confirmed baseline. v21 remains available for max diversity if all 13 re-checks are done."
```

## RECOMMENDATION
```json
"PRIORITY-DEPENDENT, not a flat 'v22 primary'. v22 = the PRAGMATIC pick: a strict improvement over the portal-confirmed v19 baseline (fixes 2 dead switchers) at minimal (2) re-check risk. For MAX #1-axis strength pick v21 instead (8 both-fold 2-3cat His vs v22's 2; cost 13 re-checks). v19 = zero-work fallback."
```

## composition
```json
[
 "des_4098_1_b14",
 "des_4098_1_c44",
 "des_180_6_H9-14-38",
 "des_84_3_mh_49_56_76",
 "des_6126_3_T20_pd_45_3_mh_9_56_78",
 "des_875_5_H3-57-80",
 "des_4098_1_c44_w6_3",
 "des_2106_2_H43",
 "des_1060_4_mh_46_74",
 "b2m0_c70",
 "b4m2_r1_w5_3",
 "des_6025_4_T10_pd_77_4_mh_52_74",
 "des_2106_2_mh_37_43_50_75",
 "des_495_6_mh_41_49",
 "des_1559_3_mh_48_75",
 "des_138_7_T10_pd_39_4_mh_49_77",
 "des_4082_8_T10_pd_10_4_mh_51_79",
 "des_4082_8_T10_pd_9_1_mh_43_55_57_79",
 "des_36133_g4_sa18",
 "des_732_1_r3_sa15"
]
```
