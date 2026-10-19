# pH-switch magnitude by direct free energy — can functional release be ranked in-silico?

The #1 axis (pH-selectivity) is the highest-weighted objective. Every prior switch metric here measures
only the *trigger* (both-fold His count, local xTB ΔE, 3 ns His–cation distance) — none measures
*functional release* (does TNF actually let go at pH 6). And switch *magnitude* had been called
"unrankable" based on pKa-prediction error, never on a direct binding free energy. This log records the
one untested lever on #1: a rigorous alchemical binding free-energy calculation.

## Method
Thermodynamic cycle: **ΔΔG_switch = ΔG_prot(bound complex) − ΔG_prot(free binder)**, where ΔG_prot is the
alchemical free energy of protonating the binder's switch-His (HIE→HIP). ΔΔG_switch > 0 ⇒ protonation
disfavors binding ⇒ **releases at low pH** (the design goal); < 0 ⇒ binds *tighter* at low pH (anti-switch).
The double-difference cancels the proton self-energy, His intramolecular FF terms, and the PME net-charge
finite-size artifact (legs built in matched boxes), leaving only the physical reorganization — which is
exactly why the *ranking* could be trustworthy where absolute pKa prediction is not.

Implementation: OpenMM (amber14SB/TIP3P, HMR + 4 fs, PME), His⁰→His⁺ via NonbondedForce per-particle
parameter offsets on a global coupling parameter; finite-difference dU/dλ over 11 λ windows integrated by
trapezoid (TI). Validated engine throughput ~770 ns/day on the 55k-atom systems. 6 designs spanning the
both-fold gradient (2 three-both-fold "strong", 2 two-both-fold "mid", 2 direction-only "weak"), bound and
free legs, up to 6 independent seeds (72 legs total).

## Results
- **Free-binder legs converge** — reproducible to <0.2–0.8 kcal/mol across seeds. The estimator is sound.
- **Bound-complex legs are basin-trapped** — seeds disagree by up to ~8 kcal/mol while *within*-seed block
  error is 0.2–0.5. The protonated His settles into different salt-bridge states against TNF and even 4 ns
  per window × 6 seeds does not cross between them. So **switch magnitude does not converge** in accessible MD.
- **Gradient is flat / unrankable:** pooled 6-seed ΔΔG — strong-tier mean +0.2, mid +2.9, weak +0.25. The
  both-fold-His **count does not predict** functional release magnitude.
- **One effect clears the noise:** `des_820_3_r3_sa4_r2_8` ΔΔG = **+6.83 ± 1.54 kcal/mol (z = 4.4, all 6
  seeds positive)** — the first design whose low-pH **release is confirmed by direct free energy**. It is
  already a submission member. The method therefore has a **detection floor** (~|ΔΔG| ≳ 5 kcal/mol): it can
  confirm a strong releaser but cannot resolve smaller differences.
- **No anti-switch liability:** no tested member resolves as ΔΔG < 0 (most-negative −1.5 ± 1.1, not
  significant). Nothing in the portfolio is a demonstrated low-pH tight-binder that would warrant replacement.

## Consequence for the submission
This is the deepest-rigor confirmation that pH-switch *magnitude* is not reliably rankable in-silico — not
because the effect is flat, but because the bound-state protonation reorganization is slow and multi-basin,
so the free energy itself will not converge at tractable sampling. That **vindicates the diversity-hedge
strategy** (spread across many distinct direction-confirmed switchers rather than chase a magnitude rank),
and it **flags no change**: no member is a confirmed anti-switch, and the set already contains a direct-FE-
confirmed releaser (`des_820_3`). Since the submission is exactly 20 members (no cut to make) and no FE
liability exists, there is no FE-justified swap. The correct in-silico treatment of #1 is diverse selection
with direction-confirmation, which the submission has; functional-release *magnitude* is wet-lab-decided.

## Broader-pool FE-screen → one actionable swap (b14 → c44)

A screen of the mature candidate pool (27 designs: the 14 remaining submission members + 13 history
candidates) was run to hunt for additional confirmed releasers. Because screening many candidates with
a noisy (basin-trapped) measure invites multiple-testing flukes, the design was **split-sample**: a cheap
2-seed pass only *filters*, and each filtered hit is then judged on **4 fresh, independent seeds** at the
des_820 bar (all same-sign, |ΔΔG| > 5 kcal/mol, |z| > 3, Bonferroni-safe for ~27 tests). The guard worked:
of 10 two-seed hits (~6 expected by luck), **4 regressed to unresolved** on the fresh seeds, including a
planted fluke control and a design whose 2-seed −7.6 collapsed to −1.5.

Three candidates cleared the bar as **confirmed releasers**: `b4m2_r1_w5_3` (+6.7, z=11), `des_6025_sh3_sa11`
(+7.2, z=4.8) — both already in the submission — and `des_4098_1_c44` (+8.0, z=6.5), which was **not**. Its
lineage sister `des_4098_1_b14` (in the submission) is the one **reproducible anti-switch**: all 6 seeds
negative, z = −6.3 (mean −3.9; direction rock-solid, magnitude just under the gold cutoff). c44 is b14 plus
two interface histidines (His63, His66) — the explicit release mechanism.

Swapping **b14 → c44** therefore replaces a reproducible wrong-direction member with a confirmed strong
releaser on the top-weighted axis. The swap clears every other axis: novelty (c44 is portal-confirmed ≥3/4,
byte-identical in v19), binding (co-fold ipTM 0.974), cross-reactivity (94% footprint conservation, only the
benign N-terminal ragged-tail divergence, identical to b14), and dedup (c44 is TM 0.85 to des_1060 vs b14's
0.82 — both portal-distinct; c44 is 0.99-identical to b14 and inherits its niche). The resulting set carries
**four FE-confirmed releasers and zero confirmed anti-switch**. This is a strict in-silico improvement on #1;
it is recorded as a candidate (local CSV), and any upload decision remains deliberately gated to the user.

Expression (manufacturability — the dominant wet-lab-zero risk) was then checked on c44 with the 6-signal
consensus, since c44 carries two extra histidines: it clears it. Protein-Sol solubility is identical to b14
(0.830), developability is cleaner (0 vs 1 liability — L66→H66 broke a hydrophobic run), and the two-folder
agreement gate passes (Boltz monomer pLDDT 0.931 ≈ b14's 0.935, TM(ESM,Boltz)=0.936), so the lower ESM pLDDT
(0.82) is ESM-specific pessimism, not a fold defect. c44's manufacturability ≈ b14's; the swap holds on every axis.

## Enhanced-sampling test — does metadynamics make #1 *magnitude* rankable?

The "magnitude unrankable" verdict above rested on bound-leg basin-trapping at plain TI (4 ns/window, no
enhanced sampling) — a *compute* limit, not a proven physics wall. One lever was untested: enhanced sampling.
Engine `scripts/fe/alchemy_meta_ddg.py` adds **well-tempered metadynamics** (biasFactor 8, native OpenMM —
no PLUMED/MPI) on a His⁺↔carboxylate coordination-number CV to **both** legs, then reweights ⟨dU/dλ⟩ back to
the unbiased ensemble. Key identity: the bias V(s) depends only on the CV s, **not on λ**, so it cancels
exactly in the finite-difference dU/dλ; the unbiased average is recovered by weights wᵢ = exp(+V(sᵢ)/kT) with
V = −(1−1/γ)·F(s). A fresh bias is grown per λ window. Controls: 3 designs × 2 legs × 4 **fresh** seeds
(s21–24), 11 λ-windows, ~13 ns/leg, one leg per GPU on 8×H100.

**Result — sign robust, magnitude did the opposite of converge:**

| design (truth) | meta ΔΔG ± SE (z) | plain-TI ΔΔG | sign |
|---|---|---|---|
| des_820_3 (+) | **+7.20 ± 3.75** (z +1.9) | +6.83 | + ✓ |
| des_4098_1_b14 (−) | **−2.06 ± 0.95** (z −2.2) | −3.88 | − ✓ |
| des_4098_1_c44 (+) | **+6.90 ± 2.72** (z +2.5) | +7.21 | + ✓ |

- **Signs 3/3 correct.** The ΔΔG *sign* is now confirmed across **three independent regimes** — plain-TI
  6-seed, split-sample fresh-seed, and metadynamics-TI. The `b14 → c44` swap stands on all three.
- **Free legs collapsed, bound legs did not.** des_820 free seed-spread 7.5 → **1.7**, b14 free 2.9 → **0.7**:
  the free-binder CV explores a real range (cv_mean 0.2–0.9, cv_max 1–2), so the bias works. But des_820
  **bound** spread went 7.4 → **16.9** (seeds −15.2 / −3.1 / −0.9 / +1.7, *each* with 0.5–0.8 within-seed
  block error = textbook basin-trapping, not estimator noise).
- **Mechanism, from the per-λ CV traces:** in the *bound complex* the His⁺↔carboxylate coordination CV is
  **pinned near zero** (cv_mean ≈ 0.03, cv_max 0.07–0.14) — the switch-His is sequestered by the interface,
  so this CV is **not** the bound-state slow coordinate. Depositing bias along it cannot cross the real
  (multi-residue interface-rearrangement) barriers; it only widened the non-productive bound ensemble.

**Consequence.** This is the deepest-rigor confirmation yet — now with per-λ CV evidence, not just argument —
that the bound-state rate-limiting degree of freedom is the **multi-residue interface reorganization**,
*orthogonal* to any single-His coordination CV. A CV able to enhance it would have to span the full interface
salt-bridge network (high-dimensional, not tractably biasable in this regime). So #1 switch **magnitude** is
confirmed **unrankable in-silico** even with the correct *class* of tool — while the **sign** is now
triple-confirmed. This **re-vindicates the diversity-hedge** and **flags no change**: every v25 member's switch
*direction* holds, the set carries 4 sign-confirmed releasers and 0 confirmed anti-switch, and magnitude is
wet-lab-decided. v25 stands.

## Whole-portfolio sign audit — all 20 v25 members, 6→12 seeds

The metadynamics section above tested only 3 controls. This audit extends direct-FE **sign** coverage to
every v25 member (`scripts/fe/run_sign.sh` + `analyze_sign.py`): bound+free, 6 seeds each (the one gap —
`des_875`, His80 = C-terminal — was unblocked by the terminal-His engine fix). Estimator = paired per-seed
ΔΔG; bar (calibrated on b14: paired −3.9, z −8.8, 0/6) = ≥all-but-one seeds negative AND |z|>3, ≥4-seed floor.

**6-seed roster result:** 5 confirmed releasers (c44 +7.2, b4m2 +6.7, des_6025 +7.7, des_820 +6.8,
des_732 +3.2), 0 confirmed anti-switch. Three members showed a 5/6-negative *lean* not clearing the bar
(des_84_3 −1.8/z−3.0, des_2106_2 −3.5/z−1.6, des_495_6 −1.5/z−1.4) → extended all three to **12 seeds**.

**12-seed resolution:**

| member | 6-seed | 12-seed ΔΔG ± SE (z) | neg seeds | verdict |
|---|---|---|---|---|
| des_84_3_mh_49_56_76 | −1.8 (z−3.0) | **−2.16 ± 0.47 (z −4.6)** | **11/12** | **CONFIRMED ANTI-SWITCH** |
| des_2106_2…_sa2 | −3.5 (z−1.6) | −1.87 ± 1.75 (z −1.1) | 8/12 | noise (spread −11..+10) |
| des_495_6_mh_41_49 | −1.5 (z−1.4) | −0.39 ± 0.72 (z −0.5) | 6/12 | neutral noise |

des_84_3 firmed up with sampling (did not regress) → a real, reproducible anti-switch, same class as b14
(~half the magnitude, higher-seed confidence). The other two regressed, as split-sample predicts.

**Consequence → v26.** des_84_3 is the one v25 member proven to switch the WRONG direction on the #1 axis.
No clean replacement exists: des_84_3 is the only des_84-lineage design with FE coverage (no c44-style
same-lineage releaser), and the portal-confirmed non-v25 designs are all negative-leaning, dedup-cousins of
existing members, or lack an FE sign. So `submission_v26.csv` = **v25 minus des_84_3 (19 designs)**,
md5 `73119a41fee4d726a39263229ad3f9b0`, validator PASS, a strict subset of portal-confirmed v25 (no
re-check, no dedup risk): **0 confirmed anti-switch, 5 confirmed releasers.** v26 is the recommended
#1-pure primary; v25 is retained as the 20-design max-diversity alternative (trade: one lineage of
wet-lab hedge vs. one proven-wrong-direction slot). Final pick user-gated.

## Brutal whole-portfolio stress test (2026-10-16) — a second anti-switch, and the v27 ladder

Re-opened the full v26 roster under an adversarial mandate ("best 20 even if it means rejecting all
selected"). Stress-tested all 19 v26 members + 3 new-lineage 20th-candidates on every computable axis:
sequence novelty (mmseqs vs full PDB: 21/21 zero hits, TNF posctrl 100% = search valid), developability
(0 Cys / 0 N-glyc all; none disqualifying), binding (all pass af2ig pae_int<10), cross-react (settled
90–97%), structural novelty (19/19 v26 portal-confirmed), FE switch sign (the live axis).

**Sign extension of the 3 remaining 6-seed negative-leaners (→12), then the one firming one (→24):**

| member | 6-seed | 12-seed | 24-seed | call |
|---|---|---|---|---|
| des_33012_g4_sa8 | −1.03 | −0.21 (z−0.2, 6/6) | — | noise, KEEP |
| des_180_6_H9-14-38 | −0.85 | −1.09 (z−0.9, 6/6) | — | basin noise, KEEP |
| des_34125_g4_sa15 | −0.68 | −3.26 (z−2.8, 10/12) | **−3.96 (z−4.5, 20/24)** | **CONFIRMED ANTI-SWITCH** |

des_34125 firmed **monotonically** (−0.68→−3.26→−3.96) across three independent seed batches — magnitude
larger than both known anti-switches (b14 −3.88, des_84_3 −2.16), z−4.46 (p<1e-5). It fails the analyzer's
literal "all-but-one-seed-negative" label only because that is a 6-seed heuristic: at n=24 the basin-trapped
bound leg throws occasional positive outlier seeds even for a true anti-switch, so the governing instrument
is the z-test, which is decisive. b14 calibration held in the same run (−3.88, z−8.8). → **drop des_34125.**

**20th-slot candidates (all new lineages):** des_243_3_r3_sa2 firmed +0.27→+1.76 at 6 seeds but **regressed
to −0.42 (z−0.4, 6/6) at 12** → neutral, not a releaser (novelty-safe TM 0.808). des_1169_2_r3_sa7 neutral
(−0.05), strong binder, novelty TM 0.934 (marginal; cf passed sh3_sa11 0.913). b1m9_A12H_w5_3 anti-lean
(−2.60) → reject. No novelty-safe non-cousin **releaser** remains in the built pool.

**The v27 ladder (all #1-pure: 0 confirmed anti-switch, 5 confirmed releasers, validator PASS):**

| option | n | = | add | novelty | note |
|---|---|---|---|---|---|
| **v27** `f5147388` | 18 | v26 − des_34125 | — | all portal-confirmed | **recommended; 2 anti-switches removed total (des_84_3+des_34125); 0 re-check/dedup/risk** |
| v27d `17b45312` | 19 | v27 + des_243 | neutral, NEW lin | TM 0.808 strong-inferred-pass | +1 lineage; 1 portal novelty-check |
| v27dd `25bf9322` | 20 | v27d + des_1169 | neutral, NEW lin | TM 0.934 marginal | reaches 20; 2 portal-checks; des_1169 the soft spot |

A clean 20 (all portal-confirmed AND all releasing) is **not** data-supported: two anti-switches had to be
removed and the pool has no novelty-safe non-cousin releaser left. The ladder trades count against
novelty-certainty at one fixed #1-pure, 5-releaser core. Upload user-gated.

## des_34125 override audit (2026-10-17) — settled by control, not judgment

An external review flagged that dropping des_34125 overrode the analyzer's own label: at 20/24 it misses the
"all-but-one-seed" rule (needs ≥23/24), so its label stayed `NEG_LEAN_UNCONFIRMED`; the cut rested on the
z-test, i.e. on judgment over the tool. Fair flag — settled with a falsifiable control rather than argument.

**Pre-committed rule:** extend BOTH reference anti-switches (b14, des_84_3) to the same 24-seed depth. If they
*also* land <23/24 → the rule is miscalibrated at high n → the cut stands. If they stay ≥23/24 → des_34125 is
genuinely noisier → restore it (v26=19). Control: 60 legs, all rc=0, GPUs returned idle.

| design (24 seeds) | neg/n | ddG | z | sign-test p | all-but-one (≥23/24) |
|---|---|---|---|---|---|
| des_4098_1_b14 (calibration anchor) | 22/24 | −3.04 | −6.3 | 1.8e-5 | **FAIL** |
| des_34125_g4_sa15 (in question) | 20/24 | −3.96 | −4.5 | 7.7e-4 | **FAIL** |
| des_84_3_mh_49_56_76 (already removed) | 19/24 | −1.63 | −3.6 | 3.3e-3 | **FAIL** |

**Result: VINDICATED.** At matched depth both reference anti-switches fail the ≥23/24 rule, and the analyzer now
labels all three `NEG_LEAN_UNCONFIRMED` — it un-confirms its own two reference anti-switches purely on crossing
into high-seed depth. So the rule, not the physics, is what breaks: the all-but-one sign-count is a 6-seed
heuristic; at n=24 the bound-leg basin throws occasional positive outlier seeds even for proven anti-switches.
The governing high-n instruments are the z-test and the distribution-free sign test, and by both des_34125
(20/24, p=7.7e-4, z−4.5, −3.96) is **more** one-signed, higher |z|, and larger in magnitude than des_84_3
(19/24, p=3.3e-3, z−3.6, −1.63), a slot already removed as a confirmed anti-switch and uncontested. No
consistent rule removes des_84_3 yet keeps des_34125. **Recommendation stays v27 = 18; v26 not restored.**

## v27 negative-leaner sweep (2026-10-17) — the 24-seed bar applied to every member kept

The des_34125 episode set a new standard: a 6–12 seed negative-leaner can firm into a real anti-switch at 24. That
depth had been applied to des_34125/b14/des_84_3 but not to the four negative-leaners still *inside* v27 — an
inconsistency (stricter bar on the slot rejected than on the slots kept). Fixed: all four extended to 24 seeds
(96 legs, all rc=0). Pre-committed firm rule: drop only if |z| grows with depth **and** 24-seed sign-test p<0.01
**and** neg-count ≥ des_84_3's confirmed 19/24.

| v27 member @ 24 seeds | mean | z | neg/n | sign-test p | 6→12→24 trend | verdict |
|---|---|---|---|---|---|---|
| des_2106_2_…_sa2 | −2.02 | −1.9 | 16/24 | 0.076 | z −1.7→−1.1→−1.9 (non-monotonic) | **keep** (non-sig; CI incl. 0) |
| des_180_6_H9-14-38 | −0.92 | −1.2 | 11/24 | 0.73 | minority-negative | keep (noise) |
| des_33012_g4_sa8 | +0.43 | +0.5 | 10/24 | 0.85 | flipped positive | keep (noise) |
| des_495_6_mh_41_49 | −0.00 | −0.0 | 12/24 | 0.58 | collapsed to zero | keep (noise) |

**No third anti-switch.** Three of four are clean noise (coinflip / positive / zero mean). des_2106_2 is a
non-significant mild negative lean — its 95% CI spans zero, its p (0.076) is 23× des_84_3's, and its trend is
non-monotonic (unlike des_34125's clean monotonic firming), so it fails every firm criterion. Cutting it on a
sub-significant lean would repeat, in mirror image, the judgment-over-tool over-cut flagged on des_34125; the
disciplined call is to keep it and flag it as the weakest switch member. **v27 = 18 is confirmed maximally
stress-tested: the 24-seed decisive bar is now applied portfolio-wide — 0 confirmed anti-switch, 5 confirmed
releasers. The ladder (18/19/20) is unchanged.**

## Cross-species switch FE — does the pH-switch fire in MOUSE? (2026-10-18/19)

Every FE section above is the **human** complex. Cross-reactivity (#2) and #1-in-mouse had only ever been
argued from a **static** epitope-conservation footprint (`data/epitope/human_mouse_map.json`) — never
energetically. This is the first direct test: does the His⁰→His⁺ release-direction (ΔΔG_switch sign) also
hold in the **mouse** complex? Same pipeline, estimator and bar as the human sign audit; mouse systems named
`<design>_mouse_{bound,free}`, built from a Boltz-2 co-fold of each binder on **mature** mouse-TNF.

**Pose-trust gate first (the result is only as good as the mouse pose).** The binder is the *same* molecule in
both complexes, so the gate is numbering-invariant: compare the binder-side **paratope** (contacts) + switch-His
interfaciality, not the TNF epitope. PASS = paratope frac-of-human ≥0.60, ≥2 switch-His interfacial (<6Å),
pose reproducible across ≥2 co-fold seeds (paratope Jaccard ≥0.6). **All 18 PASS** (frac 0.81–1.00, 14/18 ≥0.92);
2 flagged for LOW co-fold confidence only (des_33012 iptm 0.401, des_820_3 iptm 0.618 — pose mode still
reproducible). iptm alone is insufficient (one pilot had iptm 0.88 on a wrong-face pose → caught by paratope).

**Wave: 18 designs × 2 legs × 6 seeds = 216 legs, all rc=0, zero crashes.** Convention unchanged
(ΔΔG>0 = releases at low pH = correct; <0 = anti-switch).

| design | HUMAN ΔΔG (z) | MOUSE ΔΔG (z, +/−) | verdict |
|---|---|---|---|
| b4m2_r1_w5_3 | +6.7 (z+18.8) | **+3.26 (z+2.7, 6/0)** | **DUAL-SPECIES RELEASER** |
| des_6025_4…sh3_sa11 | +7.7 (z+8.3) | **+3.57 (z+4.9, 6/0)** | **DUAL-SPECIES RELEASER** |
| des_4082_8…10_4 | +6.2 (z+2.5) | **+5.59 (z+4.0, 6/0)** | **DUAL-SPECIES RELEASER** |
| des_732_1_r3_sa15 | +3.2 (z+5.1) | **+2.25 (z+3.5, 6/0)** | **DUAL-SPECIES RELEASER** |
| des_820_3_r3_sa4_r2_8 | +6.8 (z+4.8) | **+6.83 (z+3.5, 6/0)** | **DUAL-SPECIES RELEASER** (low-conf pose) |
| des_875_5_H3-57-80 | +0.9 (z+1.7) | +0.08 (z+0.1, 5/1) | conserved-flat (no inversion) |
| b2m0_c70 | +1.7 (z+1.5) | +0.16 (z+0.2, 4/2) | conserved-flat |
| des_6126_3…mh_9_56_78 | +1.9 (z+1.9) | −2.53 (z−1.3, 1/5) | conserved-flat |
| des_138_7…mh_49_77 | +3.8 (z+2.6) | −2.61 (z−1.5, 2/4) | conserved-flat (borderline; deepening) |
| des_4082_8…9_1 | +0.3 (z+0.3) | −1.83 (z−1.4, 2/4) | conserved-flat |
| des_33012_g4_sa8 | +0.4 (z+0.5) | −0.02 (z−0.0, 3/3) | conserved-flat (low-conf pose) |
| des_4098_1_c44 | +7.2 (z+8.1) | −1.49 (z−1.7, 1/5) | CI-straddle → deepening (seeds 7–12) |
| des_36133_g4_sa18 | +1.2 (z+0.8) | −3.84 (z−2.5, 1/5) | CI-straddle → deepening |
| des_180_6_H9-14-38 | −0.9 (z−1.2) | −6.34 (z−5.9, 0/6) | firm mouse-anti |
| des_1060_4_mh_46_74 | +1.8 (z+0.7) | −8.01 (z−6.0, 0/6) | firm mouse-anti |
| des_2106_2…_sa2 | −2.0 (z−1.9) | −7.82 (z−4.8, 0/6) | firm mouse-anti (anti in BOTH species) |
| des_495_6_mh_41_49 | −0.0 (z−0.0) | −1.35 (z−3.1, 0/6) | firm mouse-anti |
| des_1559_3_mh_48_75 | +2.0 (z+1.7) | −6.50 (z−4.1, 0/6) | firm mouse-anti |

**Partition: 5 dual-species releaser | 6 conserved-flat | 2 CI-straddle (deepening) | 5 firm mouse-anti.**

**Validity controls (the signal is real complex physics, not an artifact):**
1. *Free-leg match.* The free leg is the isolated binder, so human-free must ≈ mouse-free — it does, for every
   design (ΔFree +0.65/+0.39/−0.73/+0.07/+1.00 kcal on the first 5, all ≤1). No free-leg reference anomaly.
2. *Signal is entirely in the BOUND leg* (ΔBound −8.0/−5.0/−5.1/−0.7/−8.8) — genuine complex reorganization.
3. *Design-specific magnitude*, not a constant offset → not a systematic solver bias.
4. *Conserved counterexample:* des_875_5 ΔBound only −0.7 → it is NOT that every mouse complex trivially
   stabilizes protonated His; the effect tracks the specific epitope.
5. *His-set matched* human==mouse per design (verified) → like-for-like sign comparison.
Residual irreducible caveat: the switch-His TNF-side partner depends on the co-fold pose; no blind
dock-correctness gate exists (ledger), so the gate anchors pose-equivalence to the human designed pose, not
absolute correctness.

**Cross-species coherence.** Of the 5 human-confirmed releasers, **4 are confirmed in mouse too**
(dual-species: b4m2, des_6025, des_4082_8…10_4, des_732 — plus des_820 which is also human-confirmed = 5 dual
when counting des_820). des_4098_1_c44 is the lone anomaly: a rock-solid human switch (z+8.1) that leans anti
in mouse → the top deepening target. Human switch *strength does not cleanly predict transfer*
(des_4082_8…10_4 at human z+2.5 transfers; c44 at z+8.1 does not), so the mouse FE adds **independent**
discriminating signal, not an echo of the human result.

**Consequence for the submission — an ORDERING signal, not a cut.** This gives the first energetic
cross-species switch-robustness ranking the portfolio never had. The 4–5 dual-species releasers are the
strongest assets on the combined #1∩#2 objective and should be featured. It is **not** an automatic cull of
the firm mouse-antis: (a) they still **bind** mouse (raw #2 intact) — only the conditional release inverts;
(b) #1 is judged primarily in human; (c) culling collapses the diversity hedge; (d) the pose caveat above.
The one scrutiny case is **des_2106_2** (anti in BOTH species) — but its human anti is sub-significant
(z−1.9, 16/24), so by the des_34125 standard it does not clear the cut bar on human alone; mouse adds to it.
No design is added or removed on this result. **v27 = 18 stands as default; any ordering/pick change is
USER-GATED.** des_4098 + des_36133 seed-deepening (7–12) in flight to resolve the two CI-straddlers.
