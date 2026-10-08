# Challenge 2 (TNF-alpha) campaign plan

Deadline 2026-10-11 23:59 AoE (= 2026-10-12 11:59 UTC). Track 3, <=20 designs, minibinder-first.
Objective priority: pH-selective (bind 7.4 / release 6.0) > mouse cross-reactive > affinity.
GPU budget: 3 reserved (0,1,2) + any free others (verify nvidia-smi, never clobber a foreign job).

## Stage 0 - setup
- [x] workspace from the generalized toolkit; core uv env
- [x] fetch 1TNF (chains A-C, res 6-157) + mouse P06804
- [x] human<->mouse map: 79.5% id, anchors 65/67/135/146 conserved, 32 divergent (configs/target.yaml)
- [ ] toolchain build (background: proteinmpnn, rfdiffusion1, af2ig, boltz2, chai1, boltzgen, caliby, esmc, pyrosetta)
- [ ] ESMFold novelty-oracle env (envs/esmfold2-pinned) + foldseek/mmseqs
- [ ] build full-atom TWO-protomer target (chains A+B) for af2ig; define hotspots (conserved anchors at the A-B cleft)
- [ ] build mouse two-protomer model (thread P06804 on 1TNF) for the cross-react co-fold

## Stage 1 - generate (RFdiffusion1)
De-novo hotspot binders at the A-B receptor-binding cleft + partial diffusion to amplify. Minibinder 40-100 aa lead; a few microbinders (<40) and larger (>100) for category coverage. One process per free GPU.

## Stage 2 - sequence (soluble-MPNN)
0-Cys. Produce plain + His-biased variants (for both switch mechanisms).

## Stage 3 - fold/score binding
af2ig on the full-atom two-protomer target (interface-PAE gate) + Boltz-2/Chai co-fold. Full-atom target is mandatory (backbone-only target = false 0% yield).

## Stage 4 - install + score the pH switch (release at low pH)
Mechanism A: binder His next to CONSERVED target K65 -> protonation@6 -> repulsion -> release.
Mechanism B: buried His charged H-bond network in the binder core -> acid destabilizes fold (target-independent, safest for cross-reactivity).
Score: Rosetta fa_elec ddG_elec(pH5)-ddG_elec(pH7) >= 0 (repulsion rises at low pH) + PROPKA linkage + EXPLICIT-STATE MD occupancy gate (His-neutral interface intact @7.4 / His-protonated broken @6.0), tuned so the transition sits at pH 6.0. Re-validate every metric on the exact engineered sequence (submit==folded==simulated).

## Stage 5 - cross-reactivity (mouse)
Mechanism severity (contacts must avoid the 32 divergent positions; hard-avoid D140K charge reversal) + matched mouse-vs-human co-fold delta. A switch via mechanism A must anchor on a CONSERVED basic (K65), not a divergent one (R31/R131/R138).

## Stage 6 - novelty + developability
ESMFold oracle TM safe-harbor + sequence novelty; explicitly check vs the Baker 2025 pH-sensitive designs (bioRxiv 2025.09.29.678932). 0-Cys, SAP patch.

## Stage 7 - select + submit
Rank by priority order; diverse lineages; <=20; methods writeup + metrics; validate_submission.py PASS; CSV best-first. Upload USER-GATED.
