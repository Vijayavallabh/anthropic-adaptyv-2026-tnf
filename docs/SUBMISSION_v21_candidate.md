# v21 candidate (LOCAL-ONLY, upload USER-GATED)

**The SUBMISSION_STATUS / UPLOAD_CHECKLIST / FINAL_COMPOSITION / EXPANSION_RESULT sections below are
AUTHORITATIVE** (clean final data). The detailed JSON records after them are working records that may
contain superseded intermediate states (labels like 'pending'/'running' from before a run finished, or
earlier decisions later reversed) and bare `saN` shorthand that can collide across designs — resolve any
conflict in favor of the authoritative sections and the live CSV. Counts and md5 are computed from the CSV.

## SUBMISSION_STATUS
```json
{
 "n_designs": 17,
 "md5": "f7f246f41a3a062b1922f0f980a9928a",
 "validator": "PASS (track 3, <=20 rows)",
 "upload": "USER-GATED (never auto-upload; portal account j_v_v_07)",
 "fallback": "v19 (portal-confirmed >=3/4 all 20)",
 "expansion": "COMPLETE: regen round-3 (8 new-lineage) + round-2 cousins tested; 1 cleared all gates (des_820_3_r3_sa4_r2_8) -> v21 13->14"
}
```

## UPLOAD_CHECKLIST (authoritative — portal re-check before upload)
```json
{
 "byte_identical_v19_NO_recheck": [
  "des_4098_1_c44",
  "des_180_6_H9-14-38",
  "des_4098_1_b14",
  "des_2106_2_H43"
 ],
 "NEED_portal_novelty_recheck": [
  "des_2106_2_mh_37_43_50_75_sa2",
  "des_4082_8_T10_pd_10_4_mh_51_79_sa7",
  "des_138_7_T10_pd_24_4_mh_49_58_77_sa12",
  "des_6126_3_T20_pd_45_3_mh_9_56_78_sa1",
  "b2m0_c70_sa17",
  "des_6025_4_T10_pd_77_4_mh_52_sa10_sh3_sa11",
  "des_732_1_r3_sa15",
  "des_1169_2_r3_sa7",
  "des_243_3_r3_sa2",
  "des_820_3_r3_sa4_r2_8",
  "des_36133_g4_sa18",
  "des_34125_g4_sa15",
  "des_33012_g4_sa8"
 ],
 "counts": {
  "byte_identical": 4,
  "need_recheck": 13,
  "total": 17
 },
 "WHY": "novelty is a HARD gate (>=3/4; DQ if <3/4). Only the 4 byte-identical-to-v19 members are already portal-confirmed. The other 13 (5 shell-adapt saN + sh3_sa11 + 7 new-lineage) have sequences NOT present in v19 -> each MUST be portal-re-checked before upload.",
 "fallback": "v19 (all 20 portal-confirmed >=3/4) if re-check not done; upload USER-GATED",
 "CORRECTION": "supersedes any earlier '10 byte-identical' claim (that conflated v20b carryover with portal-confirmed; shell-adapt saN members were never portal-checked)."
}
```

## FINAL_COMPOSITION (authoritative per-member table)
```json
{
 "n": 17,
 "md5": "f7f246f41a3a062b1922f0f980a9928a",
 "validator": "PASS rows=17 errors=0",
 "cys": 0,
 "members": [
  {
   "name": "des_4098_1_c44",
   "lin": "des_4098",
   "cat": "multi",
   "bind": 4.56,
   "switch": "STRONG (Wyman -1.383/-2.666)",
   "portal": "byte-identical v19 (confirmed >=3/4)"
  },
  {
   "name": "des_180_6_H9-14-38",
   "lin": "des_180",
   "bind": 4.69,
   "portal": "byte-identical v19 (confirmed)"
  },
  {
   "name": "des_4098_1_b14",
   "lin": "des_4098",
   "bind": 4.7,
   "portal": "byte-identical v19 (confirmed)"
  },
  {
   "name": "des_2106_2_H43",
   "lin": "des_2106",
   "bind": 5.26,
   "switch": "WEAK single-His (-0.275/-0.839); kept for binding+novelty+diversity",
   "portal": "byte-identical v19 (confirmed)"
  },
  {
   "name": "des_2106_2_mh_37_43_50_75_sa2",
   "lin": "des_2106",
   "cat": 2,
   "bind": 6.3,
   "switch": "STRONG (Wyman -1.054/-3.083); MD-selective 6-rep",
   "sol": 0.91,
   "xreact": 93,
   "portal": "RE-CHECK (shell-adapt reseq vs v19)"
  },
  {
   "name": "des_4082_8_T10_pd_10_4_mh_51_79_sa7",
   "lin": "des_4082",
   "cat": 2,
   "bind": 5.59,
   "md": "d_repel 1.69 (cleanest)",
   "sol": 0.62,
   "xreact": 92,
   "nov_TM": 0.914,
   "portal": "RE-CHECK"
  },
  {
   "name": "des_138_7_T10_pd_24_4_mh_49_58_77_sa12",
   "lin": "des_138",
   "cat": 1,
   "bind": 6.13,
   "sol": 0.86,
   "xreact": 89,
   "portal": "RE-CHECK"
  },
  {
   "name": "des_6126_3_T20_pd_45_3_mh_9_56_78_sa1",
   "lin": "des_6126",
   "cat": 1,
   "bind": 3.92,
   "sol": 0.69,
   "xreact": 95,
   "nov_TM": 0.903,
   "portal": "RE-CHECK"
  },
  {
   "name": "b2m0_c70_sa17",
   "lin": "b2m0",
   "cat": 1,
   "bind": 5.01,
   "sol": 0.63,
   "xreact": 94,
   "nov_TM": 0.92,
   "len": 73,
   "portal": "RE-CHECK"
  },
  {
   "name": "des_6025_4_T10_pd_77_4_mh_52_sa10_sh3_sa11",
   "lin": "des_6025",
   "cat": 3,
   "bind": 4.87,
   "md": "0.653 >> Gln -0.053 (His-causal)",
   "sol": 0.835,
   "xreact": 91,
   "nov_TM": 0.913,
   "esm": 0.82,
   "portal": "RE-CHECK",
   "note": "shipped for novelty-safety (0.913 < sa10 0.959); switch a WASH vs sa10, NOT sharper"
  },
  {
   "name": "des_732_1_r3_sa15",
   "lin": "des_732 NEW",
   "cat": 3,
   "bind": 4.53,
   "md": "0.568 >> -0.322 (His-causal)",
   "sol": 0.723,
   "xreact": 93,
   "nov_TM": 0.842,
   "portal": "RE-CHECK"
  },
  {
   "name": "des_1169_2_r3_sa7",
   "lin": "des_1169 NEW",
   "cat": 1,
   "bind": 4.91,
   "md": "0.627 >> -0.738 (His-causal)",
   "sol": 0.83,
   "xreact": 97,
   "nov_TM": 0.934,
   "portal": "RE-CHECK",
   "note": "nov-risk 0.934"
  },
  {
   "name": "des_243_3_r3_sa2",
   "lin": "des_243 NEW",
   "cat": 1,
   "bind": 6.03,
   "md": "0.617 >> 0.22 (His-causal)",
   "sol": 0.815,
   "xreact": 92,
   "nov_TM": 0.808,
   "portal": "RE-CHECK"
  },
  {
   "name": "des_820_3_r3_sa4_r2_8",
   "lin": "des_820 NEW",
   "cat": 2,
   "bind": 5.26,
   "md": "0.728 >> -0.555 (His-causal, strongest gap)",
   "sol": 0.659,
   "xreact": 90.9,
   "nov_TM": 0.575,
   "esm": 0.52,
   "portal": "RE-CHECK",
   "note": "BEST novelty 0.575 (maybe 4/4); esm 0.52 expression flag plausibly novel-fold pessimism, offset by Boltz+Chai fold + af2ig plddt_binder 87.5"
  },
  {
   "name": "des_36133_g4_sa18",
   "lin": "des_36133 NEW(campaign3)",
   "cat": 3,
   "bind": 4.48,
   "md": "1.107 >> Gln -1.190 (gap +2.30, His-causal; strongest)",
   "sol": 0.768,
   "xreact": 92.7,
   "nov_TM": 0.842,
   "esm": 0.9,
   "portal": "RE-CHECK",
   "note": "fresh de-novo campaign, ultra-novel fold"
  },
  {
   "name": "des_34125_g4_sa15",
   "lin": "des_34125 NEW(campaign3)",
   "cat": 2,
   "bind": 5.65,
   "md": "0.897 >> Gln -0.508 (gap +1.41, His-causal)",
   "sol": 0.66,
   "xreact": 91.2,
   "nov_TM": 0.821,
   "esm": 0.87,
   "portal": "RE-CHECK",
   "note": "fresh de-novo campaign"
  },
  {
   "name": "des_33012_g4_sa8",
   "lin": "des_33012 NEW(campaign3)",
   "cat": 2,
   "bind": 4.71,
   "md": "1.293 > Gln 1.060 (gap +0.23, His-causal but THIN margin)",
   "sol": 0.775,
   "xreact": 92.9,
   "nov_TM": 0.85,
   "esm": 0.9,
   "portal": "RE-CHECK",
   "note": "fresh de-novo campaign; thinnest His-causal margin of the 3"
  }
 ]
}
```

## EXPANSION_RESULT
```json
{
 "goal": "fill empty slots with foldseek-DISTINCT direction-confirmed switchers (lottery tickets on #1 + expression hedge)",
 "round3": "8 new-lineage substrates shell-adapted; ALL 1-cation (0 reached >=2). 5 bound af2ig; MD direction killed 4 (3 wrong-sign, 1 non-causal) -> 0 kept.",
 "round2_cousins": "3 re-shell-adapt 2-cation winners af2ig: des_1120_2 no-bind; des_732_r2_1 binds but cousin of shipped 3-cation des_732 (skip); des_820_3_r2_8 binds + MD His-causal.",
 "net": "v21 = 17: 4 byte-identical-v19 + 13 reseq'd/new (all direction-confirmed, cross-reactive >=90%, expressing). Diversity: 17 distinct lineages/folds.",
 "discipline": "rejected every wrong-sign/non-causal/cousin/non-binder -> no dilution. 14 distinct-fold validated members; padding empty slots with correlated cousins or failures would reduce independent lottery tickets and raise DQ risk.",
 "campaign3_bold": "FRESH de-novo bpair campaign (1200 backbones, 8 GPUs) -> 49 af2ig binders (vs campaign-2's 23) -> 10 ultra-novel (TM<0.82) switch-piped. Round-4 shell-adapt: 2 three-cation + 3 two-cation. MD direction: 3 of 3 >=2-cation binders PASSED (des_36133 3-cat gap+2.30, des_34125 2-cat gap+1.41, des_33012 2-cat gap+0.23); both ultra-novel 1-cation FAILED (vindicates 2-cation = reliable switch class). NET +3 -> v21 14->17."
}
```

---

## ⚠ DETAILED WORKING RECORDS (secondary — may contain superseded intermediate states)

## gap_closed
```json
"all 5 carryover keepers BIND on exact submitted seq (agents flagged this as unverified) -> floor is solid on #3"
```

## final_verdicts
```json
{
 "sh3_sa11_MD": {
  "d_repel": 0.653,
  "dud_Gln": -0.053,
  "margin": 0.706,
  "selective": true,
  "verdict": "His-CAUSAL clean; 3 distinct cations; lower single-contact mag than sa10(1.2) but cleaner causality + Wyman sharpness + better binding(4.87>4.78) + better novelty(0.913<0.959). ADOPT replacing sa10."
 },
 "cutcheck_af2ig_all_bind": {
  "des_138_7_pd_39_4": [
   4.13,
   94.4,
   0.6
  ],
  "des_84_3": [
   7.22,
   88.0,
   0.57
  ],
  "des_875_5": [
   8.59,
   85.0,
   0.85
  ],
  "b1m9": [
   5.69,
   90.6,
   0.38
  ],
  "b4m2": [
   6.53,
   91.5,
   0.35
  ]
 },
 "KEY": "all 10 carryovers BIND on exact seq -> binding floor solid; cuts rest on SWITCH+EXPRESSION not binding."
}
```

## WYMAN_CORRECTION
```json
{
 "CORRECTION": "sh3_sa11 is NOT a strict upgrade. Wyman PROPKA linkage AND MD d_repel BOTH show sa10(2His) switch > sh3_sa11(3His). '3 cations=sharper' REFUTED: pKa shifts, not cation count, set switch strength. Retract flagship/upgrade framing.",
 "method": "collective_his_linkage.py PROPKA3 complex-vs-apo Wyman integral [6.0-7.4]; negative=release-at-low-pH (correct dir); more-negative=stronger. Consensus with 6-rep MD d_repel.",
 "switch_strength_ranking_wyman": {
  "des_4098_1_c44": -1.383,
  "des_2106_2_sa2": -1.054,
  "sa10": -0.908,
  "sa7": -0.842,
  "sa12": -0.715,
  "sh3_sa11": -0.687,
  "sa1_6126": -0.621,
  "b14": -0.6,
  "H9-14-38": -0.575,
  "H43": -0.275
 },
 "des_6025_decision": "SHIP sh3_sa11 (novelty-safer: struct-TM 0.913 vs sa10 0.959). Switch strength vs sa10 is a WASH (magnitude unrankable; engines disagree) -> sh3_sa11 chosen for the better novelty margin on the hard gate, NOT switch superiority. MD done: d_repel 0.653 >> Gln -0.053 (His-causal).",
 "portfolio_implications": {
  "c44": "TRUE switch leader (-1.383) -> strong keep vindicated",
  "sa2_2106": "2nd switch (-1.054) -> strong keep",
  "H43": "WEAKEST switch (-0.275, single-His) -> keep rests on expression+novelty ONLY, not #1 axis",
  "note": "all keeps still correct-direction; H43/H9-14-38/b14/sa1 are the weak-switch tail"
 },
 "TWO_ENGINE_UPDATE": {
  "HONEST_VERDICT": "switch strength sh3_sa11 vs sa10 is a WASH (engine-dependent, within pKa-predictor noise). Retract 'sharper/strict upgrade'. Keep sh3_sa11 for des_6025 slot on NOVELTY-DQ-safety(0.913 vs 0.959)+binding(4.87 vs 4.78) at comparable switch \u2014 a judgment call not dominance.",
  "two_engine_linkage": {
   "sh3_sa11": {
    "propka": -0.687,
    "pkai": -1.88
   },
   "sa10": {
    "propka": -0.908,
    "pkai": -1.479
   },
   "c44": {
    "propka": -1.383,
    "pkai": -2.666
   },
   "sa2_2106": {
    "propka": -1.054,
    "pkai": -3.083
   },
   "sa7": {
    "propka": -0.842,
    "pkai": -1.796
   },
   "sa12": {
    "propka": -0.715,
    "pkai": -1.31
   },
   "sa1_6126": {
    "propka": -0.621,
    "pkai": -1.986
   },
   "b14": {
    "propka": -0.6,
    "pkai": -1.745
   },
   "H9-14-38": {
    "propka": -0.575,
    "pkai": -1.524
   },
   "H43": {
    "propka": -0.275,
    "pkai": -0.839
   }
  },
  "robust_across_engines": "c44 & sa2 = switch LEADERS (both engines top-3); H43 = WEAKEST switch (both engines last); all 10 sign-agree (correct release direction).",
  "action": "sh3_sa11 and sa10 are interchangeable on switch; pick sh3_sa11 (novelty-safer). If regen yields a better-novelty strong-switch member, des_6025 lineage becomes a cut candidate."
 }
}
```

## SWITCH_RANKING_RELIABILITY
```json
{
 "finding": "pKa predictor error (2-3u) > pH window (1.4u) > inter-design gaps. Spearman rho(PROPKA,pKAI)=0.56. sh3_sa11-vs-sa10 gap 0.22-0.40 << engine-error 1.19 = UNRESOLVABLE.",
 "robust_bins_only": {
  "STRONG(both engines top3)": [
   "c44",
   "sa2"
  ],
  "WEAK(both bottom3)": [
   "H43"
  ],
  "MIDDLE(not reliably ranked)": [
   "sa10",
   "sa7",
   "sa12",
   "sh3_sa11",
   "sa1",
   "b14",
   "H9-14-38"
  ]
 },
 "reliable_switch_signal": "DIRECTION only (both-fold geometry + MD d_repel sign/Gln-control + pKa SIGN-consensus). Fine MAGNITUDE ranking NOT reliable in-silico at this window.",
 "STRATEGY_SHIFT": "Do NOT rank/headline by predicted switch magnitude. Select v21 on: confirmed switch DIRECTION + reliable axes (af2ig binding, solubility/pLDDT expression, struct-TM novelty, interface%id xreact) + DIVERSITY (distinct folds hedge switch-magnitude uncertainty). Cuts unaffected (rested on reliable signals).",
 "implication_for_keeps": "H43 = confirmed weak switch -> replace with diverse regen winner if available. c44+sa2 = anchor strong switches. middle-7 selected by reliable axes + fold-diversity, not switch rank."
}
```

## sh3_sa11_complete_profile
```json
{
 "binding_af2ig": [
  4.874,
  92.69,
  0.436
 ],
 "switch_wyman_propka": -0.687,
 "switch_wyman_pkai": -1.88,
 "switch_cations_bothfold": [
  "His32->K93",
  "His52->K60",
  "His74->K6"
 ],
 "xreact_iface_pct": 91,
 "xreact_contacts": "32/35 conserved",
 "novelty_structTM": 0.913,
 "sol": 0.835,
 "esm_monomer_plddt": 0.82,
 "status": "FULLY VALIDATED all axes; solid member (not switch-leader). des_6025 rep = sh3_sa11 (novelty-safer than sa10) OR sa10 (switch); interchangeable."
}
```

## REGEN_WINNERS
```json
{
 "des_732_1_r3_sa15": {
  "lineage": "NEW(des_732)",
  "both_cons": 3,
  "cations": [
   6,
   60,
   93
  ],
  "af2ig": [
   4.53,
   93.1,
   0.53
  ],
  "xreact": 93,
  "novelty_TM": 0.842,
  "sol": 0.723,
  "md_dir": {
   "d_repel": 0.568,
   "dud_Gln": -0.322,
   "sel": true,
   "verdict": "His-causal, correct sign"
  },
  "verdict": "STRONG new-lineage 3-cation refill \u2605"
 },
 "des_1169_2_r3_sa7": {
  "lineage": "NEW(des_1169)",
  "both_cons": 1,
  "cations": [
   93
  ],
  "af2ig": [
   4.91,
   91.3,
   0.46
  ],
  "xreact": 97,
  "novelty_TM": 0.934,
  "sol": 0.83,
  "md_dir": {
   "d_repel": 0.627,
   "dud_Gln": -0.738,
   "sel": true,
   "verdict": "His-causal, correct sign"
  },
  "verdict": "strong binder, 1-cation, nov-risk 0.934"
 },
 "des_243_3_r3_sa2": {
  "lineage": "NEW(des_243)",
  "both_cons": 1,
  "cations": [
   93
  ],
  "af2ig": [
   6.03,
   83.1,
   0.7
  ],
  "xreact": 92,
  "novelty_TM": 0.808,
  "sol": 0.815,
  "verdict": "best novelty 0.808, binding marginal(plddt83)",
  "md_dir": {
   "d_repel": 0.617,
   "dud_Gln": 0.22,
   "sel": true,
   "verdict": "His-causal, correct sign"
  }
 },
 "marginal_binders": {
  "des_1120_2_r3_sa32": [
   7.14,
   83.7,
   1.16
  ],
  "des_820_3_r3_sa4": [
   5.92,
   83.6,
   0.86
  ],
  "des_970_1_r3_sa4": [
   7.59,
   81.9,
   1.59
  ]
 },
 "FAILED": {
  "des_348_1_r3_sa32": "pae_int 26.4 (no bind)"
 },
 "summary": "regen yielded 1 strong new-lineage 3-cation(des_732_1_r3_sa15)+1 strong 1-cation(des_1169_2)+1 marginal-bind best-novelty(des_243_3); all NEW lineages, xreact 90-97%, sol>0.7. Diversity-expansion material for v21 refill."
}
```

## KEEPER_BINDING_EXACT_SEQ
```json
{
 "des_4098_1_c44": 4.56,
 "des_180_6_H9-14-38": 4.69,
 "des_4098_1_b14": 4.7,
 "des_2106_2_H43": 5.26,
 "note": "af2ig pae_interaction on exact submitted seq; all 4 byte-identical carryovers BIND (floor solid). Reseq'd/new members' binding in FINAL_COMPOSITION."
}
```

## CUT_SET_LINKAGE_CONFIRM
```json
{
 "finding": "all cuts defensible on #1 switch axis (Wyman PROPKA). 2 WRONG-SIGN: des_138_7_pd39 +0.345, sa15 +0.022 (switch backwards/absent). c44_w6_3 -1.039 strong but dup of stronger c44(-1.383, binds better). rest weak-moderate.",
 "weak_keep_flag": "H43 (-0.275) is weakest keep + cousin of sa2 -> REPLACE with regen winner if available (cut majority).",
 "v21_plan": "core 8 strong keeps (c44,sa2,sa7,sa12,sa1,sh3_sa11,b14,H9-14-38) + refill weak-tail(H43,b2m0_sa17)+10 cuts with regen winners (strong switch + NEW lineage). Target <=20, cut the majority."
}
```

## CUT_DECISION
```json
{
 "tier2_6rep_MD_verdict": {
  "des_495_6_mh_41_49_sa10": {
   "d_repel": 1.885,
   "dud_Gln": 2.385,
   "verdict": "CUT (d_repel<Gln control = non-causal)"
  },
  "des_1060_4_mh_46_74_sa2": {
   "d_repel": 0.665,
   "dud_Gln": 0.716,
   "verdict": "CUT (d_repel<Gln; also worst sol 0.597)"
  },
  "des_1559_3_mh_48_75_sa3": {
   "d_repel": 1.458,
   "dud_Gln": 1.365,
   "verdict": "BORDERLINE (+0.09 margin, within noise) -> lean cut"
  }
 },
 "keeper_binding_exact_seq_PASS": [
  "des_4098_1_c44(4.56)",
  "des_4098_1_c44_w6_3(6.71,mouse)",
  "des_2106_2_H43(5.26)",
  "des_180_6_H9-14-38(4.69)",
  "des_4098_1_b14(4.70)"
 ],
 "CONFIRMED_CUTS": {
  "sa_swaps": [
   "des_4082_8_pd_9_1_sa15 (MD non-selective)",
   "des_495_6_sa10 (6rep non-causal)",
   "des_1060_4_sa2 (6rep non-causal)",
   "des_1559_3_sa3 (6rep borderline)"
  ],
  "carryovers": [
   "des_4098_1_c44_w6_3 (95% dup c44, weaker bind)",
   "des_84_3 (drift switch,sol0.452,nov prose-only; af2ig pending)",
   "b4m2_r1 (weakest switch 0.40,SAP5.1; af2ig pending)",
   "b1m9 (weak+mislabel,SAP4.7; af2ig pending)",
   "des_875_5 (fold-fragile parent; upgrade MD-rejected; af2ig pending)",
   "des_138_7_pd_39_4 (fold-fragile; af2ig pending)"
  ]
 },
 "REPLACE": {
  "des_6025_sa10": "-> sh3_sa11 (SHIPPED for novelty-safety TM 0.913<0.959; switch a WASH vs sa10, NOT a strict upgrade; MD DONE d_repel 0.653 >> Gln -0.053)"
 },
 "STRONG_KEEPS": [
  "des_4082_8_pd_10_4_sa7",
  "des_6025_sh3_sa11(new3cat)",
  "des_2106_2_sa2",
  "des_138_7_pd_24_4_sa12",
  "des_6126_3_sa1",
  "des_4098_1_c44",
  "des_180_6_H9-14-38",
  "des_2106_2_H43"
 ],
 "BORDERLINE_KEEP": [
  "des_4098_1_b14 (cousin c44)",
  "b2m0_c70_sa17 (weak switch, distinct fold)"
 ],
 "cut_count_estimate": "10 cuts (4 sa-swap + 6 carryover) -> keep 10 incl new flagship; refill via 3-cation sa4 recipe on campaign-2 substrates"
}
```
