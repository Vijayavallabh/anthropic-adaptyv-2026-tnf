# v25 candidate — direct-FE #1-axis upgrade (b14→c44), portal-confirmed (LOCAL-ONLY, upload USER-GATED)

**`runs/push/submission_v25.csv`** — md5 `41d5e12b4d89a353c44c7c3eb1289064`, 20 designs, validator PASS.
**Portal-parsed 2026-10-15: ALL 20 PASSED, 0 needing fixes, no cousin/dedup flags.**

v25 = **v24 with exactly one swap**: `des_4098_1_b14 → des_4098_1_c44` (same des_4098 lineage). Every other
member is identical to v24. This is the single change the direct free-energy moonshot justified.

## Why the swap — direct alchemical binding free energy

The #1 axis (pH-selectivity) is the top-weighted objective, and its *magnitude* had never been attacked by a
direct binding free energy (only trigger proxies: both-fold His count, xTB, 3 ns distances). The moonshot
computed **ΔΔG_switch = ΔG_prot(bound complex) − ΔG_prot(free binder)** for the His⁰→His⁺ protonation
(TI, finite-diff dU/dλ, 11 λ-windows, box-matched legs so the proton self-energy / FF / PME-Δq² artifacts
cancel; OpenMM, ff14SB/TIP3P, HMR+4 fs). ΔΔG > 0 ⇒ protonation disfavors binding ⇒ releases at low pH.

- **`des_4098_1_b14` (was in v24) is a reproducible ANTI-switch**: all 6 seeds negative, ΔΔG = −3.9 kcal/mol,
  z = −6.3. Its direction is rock-solid — it leans toward binding *tighter* at low pH, the wrong way for a
  release switch.
- **`des_4098_1_c44` is a confirmed strong RELEASER**: fresh-seed (split-sample) ΔΔG = +8.0 kcal/mol, z = +6.5,
  all seeds positive — cleared the gold bar used to confirm des_820. Mechanistically, **c44 = b14 + His63 + His66**
  (K63H, L66H): two extra interface histidines that drive the release.

The screen that found it was multiplicity-guarded (27 mature-pool candidates; a cheap 2-seed pass only *filtered*,
and each hit was judged on 4 *fresh independent* seeds — 4 of 10 stage-1 hits regressed to noise, as expected).
See `docs/SWITCH_FE_findings.md` for the full method, the basin-trapping finding, and the detection floor.

## Gate — every axis cleared

| Axis | b14 (removed) | c44 (added) | Verdict |
|---|---|---|---|
| **#1 pH-switch (top)** | anti-switch −3.9 (z −6.3) | **releaser +8.0 (z +6.5)** | **major upgrade** |
| **#2 cross-reactivity** | 95% footprint conserved (N-term artifact only) | 94% conserved (same T1→N, P2→S) | tie (clean) |
| **Novelty (DQ hard gate)** | pass | pass — portal ≥3/4, byte-identical in v19; v25 parse Passed | pass |
| **Binding floor** | pass | co-fold ipTM 0.974, bind 4.56 | pass |
| **Dedup** | 0.82 TM to des_1060 (portal-distinct) | 0.85 TM (0.99 ≡ b14, inherits niche) | portal 20-unique ✓ |

Result: the set holds **4 FE-confirmed releasers** (des_820, b4m2, des_6025_sh3_sa11, c44) and **zero confirmed
anti-switch**.

## Composition (20)

```
des_4098_1_c44                    <- SWAPPED IN (was des_4098_1_b14)
des_180_6_H9-14-38
des_84_3_mh_49_56_76
des_6126_3_T20_pd_45_3_mh_9_56_78
des_875_5_H3-57-80
des_1060_4_mh_46_74
b2m0_c70
b4m2_r1_w5_3
des_6025_4_T10_pd_77_4_mh_52_sa10_sh3_sa11
des_2106_2_mh_37_43_50_75_sa2
des_495_6_mh_41_49
des_1559_3_mh_48_75
des_138_7_T10_pd_39_4_mh_49_77
des_4082_8_T10_pd_10_4_mh_51_79_sa7
des_4082_8_T10_pd_9_1_mh_43_55_57_79
des_36133_g4_sa18
des_732_1_r3_sa15
des_34125_g4_sa15
des_820_3_r3_sa4_r2_8
des_33012_g4_sa8
```

## Pick & caveats

**v25 is the recommended primary** — a strict improvement on the top axis over v24 (removes the one reproducible
anti-switch, adds a confirmed releaser) with no cost on any other axis, now portal-confirmed. **v24 is the
conservative alternative** (the exact earlier-blessed set) for anyone who prefers not to rely on the in-silico
FE swap. v23 / v19 are fallbacks. Honest caveats: in-silico-motivated swap; b14's anti-switch direction is
rock-solid but moderate in magnitude (−3.9, under the +5 gold cutoff for a *gold-confirmed* anti-switch); switch
*magnitude* in general remains unrankable (only the reproducible calls — b14, c44 — are trusted, and both are).
Final submit is USER-GATED.
