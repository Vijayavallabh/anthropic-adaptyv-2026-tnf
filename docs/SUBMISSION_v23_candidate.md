# v23 candidate — DEDUPED upgrade of v22 (LOCAL-ONLY, upload USER-GATED)

Built in direct response to the USER's v22 portal upload: the portal's internal cousin-dedup
flagged 4 v22 designs as 'Similar designs in this upload' (des_4098_1 x3, des_2106_2 x2), so v22
scored only 16 UNIQUE. v23 resolves each cousin cluster to one and refills with 3 distinct-lineage
direction-confirmed switchers -> 20 UNIQUE scored + higher #1-axis both-fold depth (5 vs 2).
The 2 new members des_36133_g4_sa18 + des_732_1_r3_sa15 PASSED in the v22 upload (now confirmed).

## PORTAL_V22_UPLOAD_RESULT
```json
{
 "what": "USER uploaded v22 (20) to Proteinbase; portal parse returned pass/not-passed per design",
 "passed_16": "all v22 members EXCEPT the 4 cousins; INCLUDES des_36133_g4_sa18 + des_732_1_r3_sa15 (the 2 new members) -> both now portal-CONFIRMED (novelty + dedup)",
 "not_passed_4_internal_dedup": [
  "des_4098_1_c44",
  "des_4098_1_c44_w6_3",
  "des_2106_2_H43",
  "des_2106_2_mh_37_43_50_75"
 ],
 "cause": "portal 'Similar designs in this upload' = INTERNAL cousin dedup (NOT the PDB-novelty gate). des_4098_1 had 3 cousins (b14 kept), des_2106_2 had 2 (both flagged). des_4082_8 x2 PASSED (distinct enough). => v22 scores only 16 UNIQUE.",
 "lesson": "portal penalizes cousins within an upload -> cap 1 per too-similar cluster; vindicates diverse-not-cousins even harder"
}
```

## V23_DEDUP
```json
{
 "csv": "runs/push/submission_v23.csv md5 65bf89bb38f5498ad90522542655eae2, 20 designs, validator PASS, LOCAL-ONLY",
 "construction": "v22 - {des_4098_1_c44, des_4098_1_c44_w6_3, des_2106_2_H43} + {des_34125_g4_sa15, des_820_3_r3_sa4_r2_8, des_33012_g4_sa8} (3 distinct-lineage direction-confirmed 2-cat refills)",
 "why_better_than_v22": "20 UNIQUE scored (vs v22 effective 16 after dedup) AND #1-axis both-fold depth 5 (2 three-cat + 3 two-cat) vs v22's 2; 19 distinct lineages (only des_4082_8 x2, portal-passed)",
 "recheck_burden": "3 (the refills des_34125/des_820_3/des_33012, never portal-checked) + re-confirm des_2106_2_mh alone passes dedup (low risk, novelty already v19-confirmed). 16 members already portal-passed in the v22 upload.",
 "zero_recheck_fallback": "v22 minus the 4 flagged cousins = 16 UNIQUE, ALL already portal-passed (if the user prefers 16 confirmed over 20 with 3 re-checks)"
}
```

## STATUS
```json
{
 "n": 20,
 "md5": "65bf89bb38f5498ad90522542655eae2",
 "validator": "PASS rows=20 errors=0",
 "upload": "USER-GATED"
}
```

## composition
```json
[
 "des_4098_1_b14",
 "des_180_6_H9-14-38",
 "des_84_3_mh_49_56_76",
 "des_6126_3_T20_pd_45_3_mh_9_56_78",
 "des_875_5_H3-57-80",
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
 "des_732_1_r3_sa15",
 "des_34125_g4_sa15",
 "des_820_3_r3_sa4_r2_8",
 "des_33012_g4_sa8"
]
```
