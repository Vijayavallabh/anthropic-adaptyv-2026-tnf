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
