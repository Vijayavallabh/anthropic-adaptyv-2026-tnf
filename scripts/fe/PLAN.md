# Moonshot: alchemical pH-switch ΔΔG — can functional RELEASE be ranked in-silico?

## Thesis
The #1 axis (pH-selectivity) rests only on *trigger* evidence (both-fold His count, local xTB dE,
3ns d_switch). None measures *functional release* (does TNF actually let go at pH6). Switch
*magnitude* was declared "unrankable" based on pKa-prediction error (2–3u > 1.4u window) — never on
a direct binding free-energy calc. The one observable that can rank where absolute pKa cannot:

  ΔΔG_switch = ΔG_prot(bound) − ΔG_prot(free)   [= 2.303·RT·(pKa_free − pKa_bound)]

A thermodynamic double-difference: proton self-energy, His intramolecular FF, and PME Δq² finite-size
artifact all CANCEL (box matched between legs). Only real physics (protonation-triggered
reorganization) survives → absolute FE unreliable, but the RANKING is trustworthy (HREX + seeds).
ΔΔG>0 ⇒ protonation disfavors binding ⇒ releases at low pH. Larger = sharper switch.

Validated as the only lever: SOTA scan (2026-10-13) found NO off-the-shelf method ranks pH-switch
magnitude (ArcRefine=no pH axis; no ddG predictor models His0→His+; Baker Sep-2025 ranks with
heuristics, drives selection to wet lab). Physics is the frontier.

## Method (method-research-validated)
Alchemical HIE→HIP on ALL binder His (chain A); TNF His fixed neutral (common-mode background).
NonbondedForce per-particle offsets on global 'lam': ring charges HIE→HIP, HD1 appears (q+eps 0→HIP).
Engine: envs/md2 OpenMM 8.4 CUDA (OPENMM_CUDA_COMPILER=cuda-12.3 nvcc; md/8.6.1 PTX-broken) + HMR/4fs
+ openmmtools HREX over 12 λ + MBAR (pymbar4). Box-matched bound/free (artifact cancel). ~770 ns/day.

## Test set (switch-strength gradient)
STRONG (3 both-fold): des_36133_g4_sa18 (H45,53,75); des_180_6_H9-14-38 (H9,14,38)
MID (2 both-fold):    des_820_3_r3_sa4_r2_8 (H27,33,47); des_33012_g4_sa8 (H11,23,51,53)
WEAK (dir-only):      des_495_6_mh_41_49 (H41,49); des_1559_3_mh_48_75 (H48,75)
2 seeds each; 12 win; equil 1ns + prod 5ns(bound)/4ns(free).

## Decision criteria
- DISCRIMINATES if STRONG > WEAK with (ΔΔG_i−ΔΔG_j) > 2·√(SEM_i²+SEM_j²). Then: re-rank all 20 v24
  by ΔΔG; a principled swap (central-His 2bf beating peripheral-His 3bf) = genuine #1-axis gain.
- FLAT (all within 2σ) ⇒ highest-rigor confirmation switch magnitude is unrankable ⇒ v24 diversity
  hedge vindicated. Either outcome decisive; v24 untouched unless a swap is 2σ-justified AND verified.

## Guards
Box-matched legs; assert Σdq=+N; 2σ gap rule (not raw order); cross-seed SEM; HREX mixing check;
sham control (protonate a non-interface His ⇒ expect ΔΔG≈0). Lessons: not SMD/Jarzynski (equilibrium
MBAR); re-validate every metric on the mutated artifact; small-n needs seeds+2σ.
