# Abstract / Theoretical-Physics Concepts for a pH-Switchable, Cross-Reactive TNF-alpha Binder

Scope: a dense catalog of physics frameworks rarely/never transferred into *de novo* binder design, each with (a) core idea, (b) a concrete operational hook in our pipeline (RFdiffusion / soluble-ProteinMPNN / AF2-ig+Boltz-2+Chai-1 / PyRosetta+PROPKA+pKAI / OpenMM MD + GFN2-xTB / foldseek+ESMFold, 8xH100), (c) risk/feasibility in ~4 days, (d) truly-untried vs emerging. Priority weighting throughout: switch KINETICS (fast koff at pH 6.0) > sharp OFF near pH 6.0 > human/mouse cross-reactivity > affinity at pH 7.4.

**Baseline to beat (the "state of the art" our novelties must exceed).** The most direct precedent is Ahn et al. 2025, *Computational design of pH-sensitive binders* (bioRxiv 2025.09.29.678932; PMC12621740). They design one-sided pH-dependent binders with **histidine-bias ProteinMPNN** (bias 0.5/1.0/1.5) + AF2 + Rosetta FastDesign, and score pH with a single **equilibrium electrostatics delta**: Rosetta `fa_elec` ddG at pH 5 vs pH 7 (`dddG_elec`), plus a hand "pH score" weighting interface H-bonds by burial/protonation. They use **NO PROPKA, NO explicit pKa calc, NO MD, NO kinetic/koff scoring, NO non-equilibrium thermodynamics, NO cooperativity/multi-site modeling.** Results include **TNFalpha: 79-fold weaker binding at pH 5.4 vs 7.4**; TNFR2: 122-fold (6/42 designs >2-fold pH-sensitive); EphA2: up to 1000-fold; IL-6: 6-fold (8/36). Transitions were "binary rather than gradual," and they report no cooperativity engineering. — [Ahn et al. 2025, PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740/). Implication: every framework below is "truly untried" *in the de novo binder loop* even where it is standard physics elsewhere; and the two highest-value, least-crowded gaps are (1) scoring the **kinetic** release (koff) directly and (2) engineering **cooperative sharpness** of the OFF transition — neither exists in the baseline.

## Non-equilibrium / stochastic thermodynamics (Jarzynski, Crooks) for pH-dependent binding free energy and unbinding work

### Takeaway
Jarzynski/Crooks let you recover equilibrium binding free energies and full unbinding-work distributions from short, massively-parallel steered-MD (SMD) pulls — a far more rigorous, trajectory-based alternative to single-structure PROPKA for ranking the 7.4-vs-6.0 difference, and it directly yields the mechanical "rupture" signal we care about. Feasible on 8xH100 in 4 days for ranking.

### Cited Findings
- Jarzynski's equality gives the theoretical basis for computing a potential of mean force (an equilibrium PMF) from non-equilibrium SMD pulls; SMD mimics AFM by pulling atoms at constant velocity, and is a standard enhanced-sampling route to ligand-receptor PMFs. — [kbbox SMD toolbox, HITS](https://kbbox.h-its.org/toolbox/methods/molecular-simulation/steered-molecular-dynamics/)
- Jarzynski+SMD "nonequilibrium binding free energy" ranked affinities across beta-lactamase, thrombin, factor Xa, HIV-1 protease, Mcl-1, CDK2 with experimental correlation **R = 0.80-0.89**, and is recommended for *ranking* binding affinities as accurate and robust. — [Non-equilibrium binding free energy study, PMC7590978](https://pmc.ncbi.nlm.nih.gov/articles/PMC7590978)
- Caveat from the same literature: Jarzynski can suffer relaxation/slow-convergence problems for complex protein-ligand systems and may need long pulls for reliable PMFs (exponential average is dominated by rare low-work tails). — [PMC7590978](https://pmc.ncbi.nlm.nih.gov/articles/PMC7590978)
- The rigorous pH route (next section's engine) integrates proton-linkage; Jarzynski supplies the *reference* absolute dG that linkage theory shifts with pH.

### Inferences
- **Operational hook:** Run SMD unbinding pulls of the designed binder:TNF complex in OpenMM at **two fixed protonation states** representing pH 7.4 (HIE/HID neutral His) and pH 6.0 (HIP protonated His at the switch positions, chosen from PROPKA/pKAI). Collect N pulls/state (cheap, embarrassingly parallel across 8 GPUs), compute the Jarzynski free-energy estimate and, better, use **Crooks fluctuation theorem** with bidirectional (pull+push) work to get a Bennett-acceptance-ratio (BAR) dG with error bars. The *selectivity score* is dG_unbind(7.4) - dG_unbind(6.0); the *rupture-work* and *max-force* distributions give a mechanical proxy for how abruptly the interface fails when His is charged.
- Crooks/BAR is more convergence-robust than one-directional Jarzynski; prefer it. This is a direct, trajectory-based upgrade over the baseline's single-structure Rosetta `fa_elec` delta.
- Because His protonation is set explicitly per state, this cleanly separates the pH effect from everything else — the exact variable we must engineer.

### Gaps
- No source found applying Jarzynski/Crooks *specifically to a pH-triggered protein-protein release* (the combination is novel); risk that constant-protonation SMD misses titration *during* unbinding — mitigate by also running the constant-pH route below.
- Convergence in 4 days for a protein-protein (vs protein-small-molecule) interface is uncertain; use as a *ranking* filter on a shortlist, not an absolute-dG oracle.

## Kramers' rate theory / TST / Markov state models / enhanced-sampling kinetics for koff (the RELEASE objective)

### Takeaway
The "release" objective is fundamentally kinetic (fast koff at pH 6.0), and there is a mature, directly-usable toolkit to *score koff from MD*: infrequent (frequency-adaptive) metadynamics is the fastest and most 4-day-feasible; weighted ensemble (WESTPA) and Markov state models (MSM/transition-path theory) are the rigorous options; Kramers/Grote-Hynes gives an analytic koff from a PMF + friction. None of these is used in the baseline — this is the single biggest untapped lever for our top objective.

### Cited Findings
- koff (dissociation/residence time) is "a parameter of fundamental importance"; MSMs give the metastable states and inter-state transition rates from MD and are a standard route to koff; transition-path theory on the MSM yields the ensemble of (un)binding pathways and rate-limiting steps. — [Frontiers Mol Biosci 2022 review](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.899805/full); [MSM of protein-protein encounters, arXiv 2105.02767](https://ar5iv.labs.arxiv.org/html/2105.02767); [Noe et al., TPT in protein folding, FU-Berlin PNAS](https://publications.imp.fu-berlin.de/826/1/09_PNAS_NoeSchuetteVandenEijndenReichWeikl_TPTpinWW.pdf)
- **Infrequent metadynamics (Tiwary-Parrinello):** accelerates unbinding and recovers *unbiased* koff from the mean residence time; benchmarked vs brute-force MD it reproduces the same dG profile, pathway, and residence time at **6-50x lower cost**. — [Infrequent metadynamics kinetics, PMC5451192](https://pmc.ncbi.nlm.nih.gov/articles/PMC5451192); [Frequency-adaptive metadynamics, arXiv 1802.04182](https://ar5iv.labs.arxiv.org/html/1802.04182)
- **Weighted ensemble (WESTPA):** a statistical-mechanics path-sampling strategy that returns rate constants and pathways using a fraction of brute-force MD time, by cloning/merging trajectories; explicitly applied to protein-ligand (un)binding koff and large conformational transitions. Open-source, GPU-friendly, documented tutorials. — [Weighted ensemble koff, PMC8825686](https://pmc.ncbi.nlm.nih.gov/articles/PMC8825686); [WESTPA 2.0 tutorials, LiveCoMS](https://livecomsjournal.org/index.php/livecoms/article/view/v5i1e1655)
- **Kramers rate theory:** rate of diffusive barrier crossing; in the strong-friction (diffusive) limit relevant to proteins in water, koff ~ (1/gamma)*exp(-dG‡/kT), so lowering the barrier *or* raising friction at pH 6.0 both speed release; protein barrier crossing is confirmed to be in the diffusive limit. — [Kramers/friction in protein folding, RSC c8sc01319e](https://pubs.rsc.org/en/content/articlehtml/2018/sc/c8sc01319e); [Memory-dependent friction in folding, arXiv 2208.13842](https://arxiv.org/pdf/2208.13842)

### Inferences
- **Operational hook (cheapest, do first):** **Infrequent metadynamics in OpenMM/PLUMED** on the complex, biasing a binder-TNF separation CV (e.g., interface contact number or COM distance), run once per protonation state. Report koff(6.0)/koff(7.4). Select designs maximizing this ratio. 6-50x speedup makes a ~10-design shortlist tractable in 4 days on 8xH100.
- **Rigorous confirm (if time):** WESTPA with the same CV for the top 2-3 designs — gives koff with CIs and the unbinding pathway (tells you *which* His-mediated contact ruptures first, feeding back into design).
- **Analytic fast-screen:** Combine a cheap umbrella-sampling/SMD PMF with a Kramers prefactor to get a relative koff *without* a full MSM — a quick triage before committing GPU to metadynamics/WE.
- **TST framing for the sharp OFF:** the pH 6.0 state should present a *low, broad* barrier (fast escape); the pH 7.4 state a *high* barrier (kinetic trapping = slow koff = effective binding). Design target = maximize dG‡(7.4) - dG‡(6.0).

### Gaps
- koff convergence for a real protein-protein interface in 4 days is the chief risk; mitigate by using these as *relative* (7.4 vs 6.0, design-vs-design) scores, not absolute s^-1.
- No published case of koff-based scoring inside a *de novo design* loop found — genuinely untried; CV choice (what coordinate the release proceeds along) is the main methodological unknown.

## Statistical mechanics of cooperativity (Wyman linkage, Ising/Potts, MWC, Hill) for a SHARP OFF at pH 6.0

### Takeaway
A sharp, switch-like OFF near pH 6.0 requires *cooperative* titration — a single His gives a shallow Henderson-Hasselbalch curve (Hill n=1, only ~1 pH unit from 10%->90%). Wyman-Tanford linkage says multiple coupled titratable sites are *necessary* for highly pH-dependent binding, and coupling them (Ising/Potts electrostatic interaction, or an MWC two-state conformational switch) can push the effective Hill coefficient n>1 to compress the transition into the ~0.4 pH-unit window we want. This is the most theoretically principled route to "sharp OFF" and is entirely absent from the baseline.

### Cited Findings
- Wyman's binding-polynomial / linked-functions formalism unifies all equilibrium coupling (cooperativity, allostery, pH linkage); it is the rigorous language for proton-linked binding. — [Wyman biographical memoir, NAS](https://www.nasonline.org/wp-content/uploads/2024/06/wyman-jeffries-1.pdf)
- **The linkage requirement:** proton-linkage thermodynamics dictate that *multiple* ionizable groups undergoing a pKa change on binding are *necessary* for highly pH-dependent binding; incorporating multiple histidines makes a protein highly sensitive to small pH decreases. — [Wyman linkage / histidine pH-dependence, web synthesis of PMC sources]; [Analysis of pH-dependent folding/stability via Tanford-Wyman linkage, Stony Brook](https://researchconnect.stonybrook.edu/en/publications/analysis-of-the-ph-dependent-folding-and-stability-of-histidine-p/)
- Histidine titration can be *anomalously sharp*, and multiple interacting histidines show cooperativity (sometimes negative) explained by direct electrostatic His-His coupling; interaction (Ising-like) models fit the data better than independent-site models. — [Sharp/cooperative His titration, PMC426825](https://pmc.ncbi.nlm.nih.gov/articles/PMC426825)
- **MWC in de novo design is now real:** Baker lab designed allosterically switchable assemblies "inspired by the classic Monod-Wyman-Changeux model," with ligand-binding cooperativity "comparable to classic natural systems such as haemoglobin," via rigid-body coupling of switchable hinge modules. — [De novo allosterically switchable assemblies, bioRxiv 2023.11.01.565167 / Nature 2024](https://www.biorxiv.org/content/10.1101/2023.11.01.565167.full.pdf); [IPD PDF s41586-024-07813-2](https://www.ipd.uw.edu/publication-pdfs/323/93bad0a8e360a61c7ad5034bc5291fb8/s41586-024-07813-2.pdf)
- Engineered ionizable residues buried in a protein core convert it into a pH-dependent allosteric (two-state) switch — the ensemble/MWC model of allostery as a design strategy. — [Engineered pH allostery, PMC4704092](https://pmc.ncbi.nlm.nih.gov/articles/PMC4704092)

### Inferences
- **Operational hook (direct, computable):** From constant-pH MD or a pKa panel (PROPKA+pKAI+xtb), compute the net proton uptake n(pH) on binding and fit the **effective Hill coefficient** of the binding-vs-pH curve; *select designs with n_eff > 1*. By linkage, steepness d(dG_bind)/d(pH) = -2.303 RT * (delta n_protons) — so maximize the number of switch histidines whose protonation state *flips on binding* (not just present, but coupled to the interface).
- **Ising/Potts hook for coupling:** place **2-4 histidines close enough to electrostatically couple** (a His cluster). Model the cluster as a few-site Ising system (states = protonation microstates, couplings J_ij from Coulomb/PB between sites) and solve the small partition function in Python (2^4=16 states) to predict the *collective* titration midpoint and sharpness. Design toward *positive* cooperativity (aligned couplings) so the cluster charges up abruptly just below 6.5. This is a ~50-line computation feasible immediately.
- **MWC hook for an on/off conformational switch:** couple His protonation to a two-state (bound-competent vs bound-incompetent) conformational equilibrium of the binder (e.g., a hinge/loop). The MWC allosteric constant L amplifies a modest pKa shift into a steep, switch-like OFF — the sharpest achievable mechanism. Buildable with RFdiffusion hinge motifs (precedent exists), but higher design risk.
- Clear design rule emerging: *one His = gradual; a coupled His cluster or a His-gated two-state switch = sharp.* The baseline used up to 8-11 His but did NOT engineer their coupling/cooperativity — so there is headroom to convert "binary-looking" into genuinely sharp-and-tuned.

### Gaps
- Quantitative J_ij couplings from PB are approximate; cooperativity sign (His-His repulsion often gives *negative* cooperativity, broadening the curve) must be checked — placement geometry decides sign and is the design risk.
- No found example of an explicit Ising/Potts titration model used as a *selection filter* in de novo binder design — untried.

## Linear-response theory / fluctuation-dissipation for fast mutation sensitivity (switch tuning + cross-reactivity)

### Takeaway
Linear-response (LRA) and the fluctuation-dissipation theorem (FDT) let you predict how a point mutation shifts a pKa or a binding free energy *from a single unperturbed simulation's fluctuations* — a cheap sensitivity scan to (a) tune the switch midpoint and sharpness and (b) find mutations that preserve binding across human and mouse TNF simultaneously. Standard in electrostatics/pKa work; untried as a cross-reactivity co-optimization tool.

### Cited Findings
- **LRA for pKa:** the PDLD/S-LRA (linear-response approximation) method computes pKa by averaging electrostatic energies over trajectories of both charge states — a workhorse semi-macroscopic electrostatics method. — [LRA/PDLD-S pKa, PMC1305621](https://pmc.ncbi.nlm.nih.gov/articles/PMC1305621/)
- A residue's protonation propensity is set by the electrostatic environment of all partial charges (incl. waters); continuum-electrostatics/PB methods and now ML predictors (trained on continuum data, ~1000x faster) estimate pKa shifts — i.e., the LRA observable is cheaply differentiable w.r.t. perturbations. — [ML-accelerated pKa predictors](https://par.nsf.gov/servlets/purl/10382876)
- **FDT = fluctuations encode response:** near equilibrium, the linear response of an observable to a perturbation equals an equilibrium correlation/covariance — "fluctuations and response encode the same information." For receptor binding, sensitivity is bounded by a Hill function whose coefficient can be enhanced by driving. — [FDT modern perspective, arXiv 2312.10134](https://ar5iv.labs.arxiv.org/html/2312.10134); [FDT bound on receptor sensitivity, PRL 115.118102](https://link.aps.org/accepted/10.1103/PhysRevLett.115.118102)

### Inferences
- **Operational hook (pKa tuning):** FDT/LRA estimate of a mutation's effect on a His pKa: from one constant-pH or fixed-charge MD of the wild design, the shift in a site's electrostatic potential upon a nearby mutation is approximated by the *equilibrium covariance* of that potential with the perturbing charge — no new simulation per mutant. Rank candidate mutations that move the switch His pKa into the [6.0, 6.5] sweet spot.
- **Operational hook (cross-reactivity):** Linear-response/thermodynamic-perturbation estimate of d(dG_bind) for the *same* binder mutation against human vs mouse TNF from the covariance of interaction energy with the mutation site, computed on human and mouse complexes in parallel. Select mutations with small |d(dG)| on *both* species (dual-robust) — a principled alternative to brute-force double mutant scans. The 79% human/mouse identity means the conserved-epitope interactions dominate; LRA identifies which contacts are species-invariant.
- This is the cheapest framework here: the covariances come almost for free from trajectories we already run for the kinetics/thermo scores.

### Gaps
- LRA assumes small, linear perturbations; a His->something charge change is large and may break linearity — use LRA for *triage/ranking*, confirm top hits with explicit mutant MD/xtb.
- No found precedent for FDT-covariance mutation scanning as a cross-reactivity co-design filter — untried transfer.

## Renormalization-group / coarse-graining / effective-field-theory for multiscale binder energetics

### Takeaway
RG/coarse-graining is the rigorous bridge between all-atom detail and the few collective coordinates that actually govern the switch; its practical, 4-day-usable face is relative-entropy / force-matching coarse-graining and Mori-Zwanzig memory kernels (which also feed Kramers friction). As a *design principle* it says: identify the few relevant order parameters (the switch reaction coordinate) and integrate out the rest. Most speculative of the "classical" entries; low direct ROI vs the kinetics/cooperativity levers.

### Cited Findings
- Repeated coarse-graining *is* the core idea of RG; the Mori-Zwanzig formalism projects many-body dynamics onto a few coarse coordinates (and yields the memory/friction kernel used in generalized-Langevin Kramers rates). — [Multiscale CG of protein energy landscape, PMC2891700](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2891700/)
- **Relative-entropy coarse-graining:** a statistical-mechanical variational principle — minimize the information lost on coarse-graining — gives a universal route to CG interaction potentials; force-matching atomistic sims reproduces folding/aggregation energetics. — [PMC2891700](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2891700/); [Relative-entropy coarse-graining, AIChE](https://www.aiche.org/conferences/aiche-annual-meeting/2011/proceeding/paper/37c-coarse-graining-proteins-using-relative-entropy-theory)

### Inferences
- **Operational hook:** use a memory-kernel/Mori-Zwanzig estimate of *friction along the unbinding CV* to supply the Kramers prefactor (ties directly into the koff section) — the one concrete, high-value use of this machinery for us. The friction may itself be pH-dependent (a protonated, more rigid interface vs a loosened one), giving a second koff lever beyond the barrier height.
- As *selection*: project designs onto a low-dimensional order-parameter space (e.g., via TICA/VAMP from the MD) and keep designs whose pH 6.0 and 7.4 ensembles are *well-separated* in the few relevant coordinates — a clean, RG-flavored separability score.
- Honest assessment: full RG on a binder energy landscape is not a 4-day deliverable; the *coarse-coordinate / friction* fragments are, and they are better pursued inside the MSM/Kramers workstream than as a standalone.

### Gaps
- No precedent for RG/EFT *per se* in binder design (as opposed to CG force fields); the "effective field theory of an interface" is an analogy, not an established method — highest-novelty, lowest-maturity entry.

## Topological-physics analogies (topological protection, Berry phase) for robust cross-reactivity

### Takeaway
Topological protection — a property robust to local perturbation unless a qualitative (gap-closing) change occurs — is a compelling *conceptual* frame for "protected" binding that survives the 21% human->mouse sequence change (cross-reactivity as robustness to disorder). There is real theory connecting topological protection to non-equilibrium biochemical networks, but no operational protein-*design* method; treat as inspiration/framing, not a 4-day scoring hook. Most exotic entry.

### Cited Findings
- Topologically protected modes exist in *non-equilibrium stochastic* (biochemical) systems: non-equilibrium fluxes can support protected modes of topological origin, robust and insensitive to local perturbation, analogous to topological boundary modes in electronic/mechanical systems; invoked for robust sensory adaptation and enzymatic specificity despite thermal noise. — [Topologically protected modes in non-equilibrium stochastic systems (Murugan & Vaikuntanathan), PMC5234070](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5234070/); [Topology in biology, MPI-DS](https://www.ds.mpg.de/3763345/210723_topology)
- An adiabatic geometric **Berry phase** arises from nonlinear normal modes; its quantized integer invariant protects nonlinear topological boundary modes from disorder — a mechanical/vibrational analog potentially mappable to protein normal modes. — [Topological Berry phase in nonlinear normal modes, Springer Nature Communities](https://communities.springernature.com/posts/topologically-invariant-berry-phase-defined-in-nonlinear-normal-modes-for-general-and-strongly-nonlinear-dynamics)

### Inferences
- **Weak operational hook (robustness-as-selection):** interpret cross-reactivity as *spectral robustness* — compute the binder's interface elastic-network / normal modes (free in PyRosetta/Python) and keep designs whose binding-relevant soft modes are *insensitive* to the human->mouse residue substitutions (small spectral shift, no "gap closing"/mode reordering). This operationalizes "topological protection" as a mode-stability score without needing literal invariants.
- **Framing value:** design for binding mediated by a *distributed, conserved-backbone* contact pattern (robust, "protected") rather than a few side-chain-specific contacts (fragile) — consistent with the 79%-identity conserved epitope strategy and with our prior lesson that distributed binding behaves differently from hotspot binding.
- Berry-phase-per-se has no credible binder-design hook in 4 days; the elastic-network robustness proxy is the salvageable, computable piece.

### Gaps
- No source applies topological invariants to protein-protein *binding affinity/specificity* design — purely analogical; high risk of over-claiming. Report as inspiration, with the normal-mode-robustness proxy as the only concrete derivative.

## Preorganized electric fields and vibrational Stark effect for interface electrostatics and protonation tuning

### Takeaway
The vibrational Stark effect (VSE) framework — quantify the local electric field a bond/charge experiences, and design the field to be "preorganized" — is a mature enzyme-catalysis paradigm (Boxer/KSI) that transfers cleanly to tuning a switch His pKa and interface electrostatics. It gives a single scalar (projected field in MV/cm) that is cheap to compute and directly controls both the His pKa and the strength/direction of the pH-dependent interaction. Strong, under-used, 4-day-feasible transfer.

### Cited Findings
- VSE spectroscopy measured the electric field inside ketosteroid isomerase's active site and found it "nearly perfectly oriented to stabilize the transition state"; the field is **preorganized** — the active site is set up with respect to the functionally relevant field changes. — [Preorganized electric field in KSI (Boxer lab), PMC7474534](https://pmc.ncbi.nlm.nih.gov/articles/PMC7474534); [Boxer lab VSE PDF](https://boxerlab.stanford.edu/_files/ugd/4006c9_bdbbac6fba0b497392ce44c88736a589.pdf)
- VSE provides a direct experimental measure of the electric field a substrate feels in an active site and is compelling evidence for a dominant electrostatic contribution to catalysis — establishing "electric field" as a designable, quantifiable observable. — [PMC7474534](https://pmc.ncbi.nlm.nih.gov/articles/PMC7474534)
- (Foundational: Fried, Bagchi & Boxer, *Science* 2014, 346:1510, quantified the KSI field at ~-144 MV/cm and linked it to rate — the canonical field-catalysis result.)

### Inferences
- **Operational hook (pKa tuning via field):** compute the **electrostatic field projected onto the switch His** imidazole from the rest of the designed binder+target (trivial from partial charges / APBS / xtb). A field that *disfavors* the His+ cation at pH 7.4 keeps His neutral (binding ON); the environment should *permit* protonation as pH drops to 6.0. Rank/tune designs by this projected field — a physically-grounded, single-number pKa dial that complements PROPKA/pKAI.
- **Operational hook (interface strength):** VSE logic says the pH-dependent interaction's strength is set by the field across the key interface contact. Select for a *large field change* between the His-neutral and His-protonated states at the interface — i.e., maximize the electrostatic *contrast* that drives the OFF switch. GFN2-xTB with --alpb water (already in our stack) can give polarized fields more faithfully than fixed-charge FF.
- Preorganization principle for cross-reactivity: a *preorganized* binder (interface field already set in the apo state, minimal reorganization on binding) pays less reorganization penalty and is more tolerant of the human/mouse substitutions — select for small apo->bound field reorganization.

### Gaps
- VSE/field-catalysis is established for *enzyme TS stabilization*, not for *pH-switch binder* tuning — the transfer is novel though low-risk. No calibration of "field -> His pKa shift" magnitude for our specific geometry; would need a quick internal check against PROPKA/pKAI on a few designs.

## Polymer-physics / phase-behavior theory for an acid-destabilized-fold ("conformational collapse") switch

### Takeaway
An alternative to a local interface-breaking His is a *global* mechanism: engineer the binder (or a subdomain) to sit near a coil-globule / fold-unfold phase boundary so that His protonation at pH 6.0 tips it across into a molten-globule/unfolded state that cannot bind — an intrinsically *sharp* (near-first-order) transition with fast kinetics. Polymer/phase theory provides the order parameter and the sharpness; buildable in principle but the highest design-risk of the switch mechanisms. Untried for binders.

### Cited Findings
- Coarse-grained simulations capture pH-dependent folding/unfolding of coiled coils; a pH drop protonates acidic side chains and raises overall hydrophobicity, driving conformational change. — [pH-dependent coiled coils (CG MD), arXiv 1304.7241](https://arxiv.org/pdf/1304.7241); [pH-driven coil-helix at interface, arXiv 1405.2886](https://arxiv.org/pdf/1405.2886)
- The coil-globule transition is an entropy(coil)->energy(globule) switch and can be continuous *or discontinuous* (first-order-like, i.e. sharp) depending on conditions; at low pH many proteins occupy a compact **molten-globule** state that is globular but lacks defined structure. — [Globule-coil mean-field theory, arXiv cond-mat/0612256](https://arxiv.org/pdf/cond-mat/0612256); [Coil-globule lattice framework, APS 2023](https://meetings-archive.aps.org/mar/2023/s04/12)
- De novo pH switches at the fold level are real: a designed peptide switches between a parallel pentamer (pH 8) and a trimer of antiparallel dimers (pH 6) driven purely by Glu protonation — demonstrating designable, protonation-triggered quaternary/fold reconfiguration. — [De novo pH-responsive assemblies (web synthesis; see also OSTI pH-responsive filaments)](https://www.osti.gov/pages/biblio/2332902-de-novo-design-ph-responsive-self-assembling-helical-protein-filaments)

### Inferences
- **Operational hook:** place titratable groups (His, or buried Glu/Asp) in the binder *core* so their protonation at pH 6.0 introduces buried charge / electrostatic frustration that destabilizes the fold. Score the fold stability at both protonation states with xtb/MD and ESMFold/AF2 pLDDT + an explicit **dG_fold(6.0) - dG_fold(7.4)**; select designs poised so that only the pH-6.0 state unfolds/molten-globules. The coil-globule order parameter (radius of gyration, contact order) from short MD gives the sharpness readout.
- **Kinetic bonus:** an unfolding-coupled release can be *faster* than a local contact breaking (global collapse removes the whole interface at once) — directly serving the koff objective; and a near-first-order collapse is intrinsically sharp, serving the sharp-OFF objective.
- Our memory flags buried-acid pKa as unpredictable (PROPKA vs pKAI disagree 2-3 units) — so prefer a **buried-His** core-destabilizer (natural in-window pKa) over buried Glu/Asp for this mechanism.

### Gaps
- Designing a binder that is folded+binding at 7.4 but *reliably* unfolds at 6.0 without being marginally stable (aggregation/expression risk) is hard in 4 days; the fold-stability predictors (ESMFold/AF2) are not calibrated for marginal/molten-globule states. High design risk; treat as a bold secondary track, not the primary switch.
- No precedent found for a *binder* using an acid-unfolding core switch (vs antibodies using local His) — genuinely untried, correspondingly risky.

## Cross-cutting synthesis (for the report-writer)

### Takeaway
Ranked by (impact on top objectives x 4-day feasibility): **(1) koff scoring via infrequent metadynamics / Kramers** (directly scores the kinetic release, the #1 objective, and is absent from all prior work); **(2) Wyman-linkage + Ising/MWC cooperativity** (the principled route to a sharp OFF; absent from baseline); **(3) constant-pH MD with Wyman-Tanford linkage** as the rigorous pH-dG engine (replaces single-structure PROPKA); **(4) preorganized-field / VSE** His-pKa tuning (cheap, physically grounded); **(5) LRA/FDT** sensitivity scans for cross-reactivity (nearly free); then the bolder/riskier **(6) acid-unfolding polymer-physics switch**, **(7) RG/friction** fragments, and **(8) topological-robustness** framing.

### Cited Findings
- The rigorous pH-dG engine: constant-pH MD + **Wyman-Tanford linkage**, d(dG)/d(pH) = -2.303 RT (n_monomer - n_dimer), integrated via cubic Hermite splines, with **per-residue decomposition** dG_i = integral of -2.303 RT (n_i,M - n_i,D) d(pH) identifying which titratable sites drive the pH dependence; d(n)/d(pH) = -2.303 var(n) and per-site = -2.303 cov(n_i, n). Tools: GROMACS+MEAD+PETIT (we can substitute OpenMM constant-pH / our pKa stack). — [da Rocha et al. 2022, JCTC 18:1982, PMC9775224](https://pmc.ncbi.nlm.nih.gov/articles/PMC9775224)
- The variance/covariance identities above *are* the fluctuation-dissipation statement for proton binding — the cooperativity (var(n) = sharpness) and per-residue-linkage (cov) and LRA sensitivity frameworks are mathematically the same object computed from one constant-pH trajectory. This is the key unifying insight: **one constant-pH MD ensemble per design, per protonation window, yields simultaneously the pH-dG (linkage), the sharpness/Hill (variance), the per-residue drivers (covariance), and the mutation sensitivities (LRA).**
- Precedent anchor for the whole effort: pH-dependent antibody engineering relies on His (pKa ~6.0, protonated below 6.5) to give *faster dissociation at acidic pH* for endosomal release/recycling — our exact direction (bind 7.4, release 6.0) and confirmation that His is the right sensor. — [pH-dependent antibody histidine engineering, PMC4622719](https://pmc.ncbi.nlm.nih.gov/articles/PMC4622719); [PMC4623423](https://pmc.ncbi.nlm.nih.gov/articles/PMC4623423)

### Inferences
- **Minimum-regret 4-day plan:** run one short **constant-pH (or two-fixed-state) OpenMM ensemble** per shortlisted design; from it extract (A) Wyman-linkage dG_bind(pH) and per-residue drivers, (B) variance-based Hill/sharpness, (C) LRA cross-reactivity sensitivities; then run **infrequent metadynamics** on the top few for koff(6.0)/koff(7.4); add the **projected-field (VSE)** scalar and the **Ising-cluster** cooperativity prediction as cheap Python post-processing. This single pipeline operationalizes six of the eight frameworks from essentially one simulation per design.
- Against the Ahn et al. baseline (single-structure Rosetta `fa_elec` delta, no kinetics, no cooperativity), this is a strict superset that adds the two things the baseline provably lacks and that our objectives most reward: the **kinetic** release score and the **cooperative-sharpness** score.

### Gaps
- Every framework's 4-day convergence on a real protein-protein interface is the shared risk; all should be used as *relative/ranking* filters over a shortlist, validated by our existing orthogonal oracles, not as absolute predictors.
- I did not find any single published pipeline combining non-equilibrium kinetics + cooperativity + linkage for *de novo pH-switch binder* design — the integration itself is the novel contribution.
