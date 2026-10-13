# v24 candidate — MAX #1-axis, dedup-clean (LOCAL-ONLY, upload USER-GATED)

The max pH-selectivity (#1-axis) option: 8 both-fold-confirmed multi-His members (= v21's depth)
in a 20-UNIQUE dedup-clean set, needing only 3 portal re-checks (vs v21's 13). Built on the
now-fully-confirmed v23 by swapping 3 members for v21's solid >=2-cat both-fold saN variants.
Choose v24 for max #1-axis depth (3 re-checks); v23 if you want it fully confirmed now (0 re-checks).

## V24_MAX1_DEDUP_CLEAN
```json
{
 "csv": "runs/push/submission_v24.csv md5 34f3c096bf63dd9ef70d05308c6c3f1d, 20 designs, validator PASS, LOCAL-ONLY",
 "construction": "v23 with 3 same-lineage swaps to v21's SOLID >=2-cat both-fold saN variants: des_2106_2_mh->sa2 (K93,K60), des_4082_8_pd_10_4->sa7 (K93,K6, d_repel 1.69), des_6025_4_pd_77_4->sa10_sh3_sa11 (K60,K6, d_repel 1.2). EXCLUDED des_6126_3 sa1 (only 1 both-fold cat) + b2m0 sa17 (weak 2-cat d_repel 0.208) -> kept v23 versions.",
 "both_fold_depth": "8 (2 three-cat: des_36133,des_732 + 6 two-cat: des_34125,des_820_3,des_33012,des_2106_2_sa2,des_4082_8_sa7,des_6025_4) = MATCHES v21's 8, the MAX #1-axis depth, in a dedup-clean 20-UNIQUE set",
 "recheck_burden": "3 (the saN swaps, never portal-checked; 17 others already portal-passed in v22/v23 uploads). vs v21's 13 re-checks + v21 still needs its own dedup.",
 "supersedes": "v21 (same 8 both-fold, but dedup-clean + 3 vs 13 re-checks)",
 "composition": [
  "des_4098_1_b14",
  "des_180_6_H9-14-38",
  "des_84_3_mh_49_56_76",
  "des_6126_3_T20_pd_45_3_mh_9_56_78",
  "des_875_5_H3-57-80",
  "des_1060_4_mh_46_74",
  "b2m0_c70",
  "b4m2_r1_w5_3",
  "des_6025_4_T10_pd_77_4_mh_52_sa10_sh3_sa11",
  "des_2106_2_mh_37_43_50_75_sa2",
  "des_495_6_mh_41_49",
  "des_1559_3_mh_48_75",
  "des_138_7_T10_pd_39_4_mh_49_77",
  "des_4082_8_T10_pd_10_4_mh_51_79_sa7",
  "des_4082_8_T10_pd_9_1_mh_43_55_57_79",
  "des_36133_g4_sa18",
  "des_732_1_r3_sa15",
  "des_34125_g4_sa15",
  "des_820_3_r3_sa4_r2_8",
  "des_33012_g4_sa8"
 ]
}
```

## V23_PORTAL_CONFIRMED
```json
{
 "event": "USER uploaded v23 (20) to Proteinbase 2026-10-13",
 "result": "ALL 20 PASSED \u2014 '20 included, 0 needing fixes'. Novelty gate + internal dedup both clear for the FULL set.",
 "resolves": "the 3 refills (des_34125_g4_sa15, des_820_3_r3_sa4_r2_8, des_33012_g4_sa8) PASS novelty; des_2106_2_mh passes dedup alone (H43 removed). ZERO re-checks remaining.",
 "status": "v23 is now FULLY portal-confirmed, 20 UNIQUE, dedup-clean. Strictly > v22 (16 unique): +4 unique scored + #1-axis both-fold depth 5 vs 2. v23 = the confirmed PRIMARY upload-ready submission.",
 "md5": "65bf89bb38f5498ad90522542655eae2"
}
```

## SWITCH_FC_ORDINAL
```json
{
 "insight": "both-fold-confirmed COUNT of separated in-window His is a sound ORDINAL proxy for expected fold-change on the #1 axis (Wyman linkage is additive over titratable groups: more separated His each repelling a distinct conserved cation -> larger |\u0394n_H| -> larger FC; Baker multi-His SOTA). NOT a precise rank (per-His pKa-in-window varies), but a robust ordinal + a geometric-confidence signal (both-fold = Boltz\u2229Chai agree the His engages its cation).",
 "per_submission_both_fold_counts": {
  "v21": "3-cat=3, 2-cat=5, 1-cat=5 (8 members with 2-3 separated His) + 4 carryover",
  "v22": "3-cat=2 (des_36133, des_732) + 18 raw v19 multi-His (direction-OK via MD, NOT both-fold geometry-confirmed)",
  "v19": "0 both-fold + 20 raw multi-His (direction-OK, not both-fold)"
 },
 "IMPLICATION": "On the #1 axis (pH-selectivity), v21 is the STRONGEST: 8 members with 2-3 both-fold-confirmed separated His = highest expected-FC distribution + highest geometric confidence. v22/v19 have direction-confirmation (MD) but little/no both-fold geometry confirmation. This REVISES the earlier v22-primary call."
}
```
