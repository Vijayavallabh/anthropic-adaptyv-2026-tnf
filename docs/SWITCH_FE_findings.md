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
