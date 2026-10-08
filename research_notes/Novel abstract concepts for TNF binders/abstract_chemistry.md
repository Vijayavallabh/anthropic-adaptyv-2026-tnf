# Novel Abstract / Theoretical-Chemistry Concepts for a pH-Switchable, Cross-Reactive De Novo TNF-alpha Binder

Scope note for the report writer: this catalog targets a de novo PROTEIN binder (20 natural AAs only) against human TNF-alpha trimer (1TNF), binding the receptor-binding groove at the two-protomer interface (anchors K65/Q67/E135/E146). Priorities: (1) sharp COOPERATIVE OFF-switch near pH 6.0 (bind pH 7.4, release pH 6.0); (2) mouse cross-reactivity (79% id); (3) affinity at pH 7.4. Toolchain: RFdiffusion, soluble-ProteinMPNN, AF2-initial-guess + Boltz-2 + Chai-1, PyRosetta/PROPKA/pKAI, APBS, OpenMM MD + GFN2-xTB, foldseek/ESMFold. Baker precedent to go BEYOND: interface His-next-to-Arg/Lys repulsion, and buried-His charged-H-bond-network fold destabilization (bioRxiv 2025.09.29.678932). Each concept is tagged with an operational hook, a 4-day feasibility/risk read, and a truly-untried-vs-emerging verdict.

## Q1: Charge-regulation theory (Kirkwood-Tanford, Linderstrom-Lang, Lund-Jonsson) — can we design an interface whose net electrostatic PPI turns repulsive at pH 6.0, and score it?

### Takeaway
Charge regulation gives a single, computable design observable — the protein "charge capacitance" C = <Z^2> - <Z>^2 (equivalently the titration-curve slope dZ/dpH) — that peaks exactly where titratable groups ionize. Engineering the interface so its capacitance peaks near pH ~6.5 while the APBS mean-field interaction energy crosses from attractive (pH 7.4) to repulsive (pH 6.0) is a near-untried, directly scoreable switch criterion that goes beyond the Baker single-His-repulsion motif.

### Cited Findings
- Charge capacitance is defined as the variance of the net charge, C = <Z^2> - <Z>^2, and is "an intrinsic property on par with the net charge and the dipole moment"; it is obtainable experimentally from the derivative of the titration curve (dZ/dpH) or from simulation — [Charge regulation in biomolecular solution, Quarterly Reviews of Biophysics](https://www.cambridge.org/core/journals/quarterly-reviews-of-biophysics/article/abs/charge-regulation-in-biomolecular-solution/E4145C91436EA0EA2874A5A20233330D).
- Charge regulation is the perturbation of titratable-site proton equilibria when a protein is approached by "other proteins, surfaces and membranes," generating a charge response fully described by the capacitance; the capacitance "directly enters the free energy expression for intermolecular interactions" — [QRB review](https://www.cambridge.org/core/journals/quarterly-reviews-of-biophysics/article/abs/charge-regulation-in-biomolecular-solution/E4145C91436EA0EA2874A5A20233330D).
- The fluctuation-induced (Kirkwood-Shumaker) interaction scales as R^-2 and coincides with the classical Kirkwood-Shumaker result when the capacitance definition is substituted; it is a thermodynamic perturbation around the non-interacting state accounting for proton charge fluctuations — [Kirkwood-Shumaker interactions and general thermal fluctuation forces, arXiv:1401.5252](https://arxiv.org/pdf/1401.5252).
- Induced ionization on approach can add "several kT" of attraction (charge-induced-charge), i.e., charge regulation is a sizeable, pH- and salt-tunable force — [arXiv:1401.5252](https://ar5iv.arxiv.org/html/1401.5252); Lund-Jonsson corpus — [On the charge regulation of proteins (Lund group)](https://portal.research.lu.se/en/publications/on-the-charge-regulation-of-proteins/), [enhanced protein adsorption due to charge regulation](https://portal.research.lu.se/en/publications/enhanced-protein-adsorption-due-to-charge-regulation/).
- Historical frameworks the review builds on: Kirkwood & Shumaker (1952) proton-fluctuation forces, and Linderstrom-Lang (1924) "smeared charge" model treating ionizable groups as a continuous surface charge rather than discrete points — [QRB review](https://www.cambridge.org/core/journals/quarterly-reviews-of-biophysics/article/abs/charge-regulation-in-biomolecular-solution/E4145C91436EA0EA2874A5A20233330D).

### Inferences
- OPERATIONAL HOOK (directly implementable): add a "capacitance-peak + sign-flip" score. (i) Run PROPKA/pKAI on the bound complex to get per-residue titration curves, sum to Z_binder(pH) and Z_interface(pH), and compute C(pH) = -(1/ln10)*dZ/dpH numerically; reward designs whose interface capacitance is MAXIMIZED at pH ~6.3-6.6 (that is exactly what packing several His at the interface does, because capacitance peaks when pH ~ pKa). (ii) Run APBS (nonlinear Poisson-Boltzmann) twice, at protonation states for pH 7.4 vs pH 6.0, and compute the binding electrostatic free energy dG_elec(pH) = G_complex - G_binder - G_target for each; the design criterion is dG_elec(7.4) < 0 (attractive) AND dG_elec(6.0) > 0 (net repulsive). This is the mean-field complement to PROPKA Wyman linkage they already use, but it scores the SIGN of the force, not just the pKa shift — the missing selectivity axis.
- Because the Kirkwood-Shumaker capacitance force is ALWAYS attractive, the OFF-switch must come from the MONOPOLE-MONOPOLE + dipole (mean-field) term flipping sign, not from the fluctuation term. Design goal = make the binder+interface net charge of the contact patch go from near-complementary (attractive) at 7.4 to like-signed (repulsive) at 6.0. Concretely: cluster 2-3 interface His so protonation adds +2/+3 to a patch that faces a conserved basic TNF anchor (K65) — collective repulsion, not a single His-Arg clash.
- A smeared-charge (Linderstrom-Lang/Tanford-Kirkwood) surrogate is cheap enough to screen thousands of RFdiffusion backbones BEFORE folding: approximate each patch as a disc of smeared charge Z(pH)/area and score the cross-patch Coulomb sign flip. Use it as a pre-filter feeding AF2-ig.
- Salt dependence is a free robustness knob: charge-regulation forces weaken with ionic strength, so test the sign flip at 150 mM (physiological) in APBS to avoid a switch that only works at low salt.

### Gaps
- No open source gave the exact closed-form two-protein charge-regulation energy with the correct electrostatic prefactor for a finite-size, non-spherical interface; the R^-2 far-field scaling and product-of-capacitances dependence are established but a docked-interface geometry needs the full PB calc (APBS), not the analytic monopole formula.
- Whether PROPKA/pKAI titration curves are accurate enough to compute dZ/dpH (capacitance) for a de novo interface is unverified; capacitance is a derivative and will amplify pKa error (see Q8 for xtb cross-check).

## Q2: Coupled-pKa networks / proton tautomerism / proton wires for a COOPERATIVE multi-His switch (Hill n>1) sharp at pH 6.0

### Takeaway
Cooperativity (Hill n>1, a sharp OFF near pH 6.0) is a SOLVED, protein-encodable problem: Baker/Boyken designed buried His hydrogen-bond-NETWORKS (3 His each) whose protonation is cooperative, and showed transition pH and steepness are tunable by the number of His networks vs hydrophobic layers. Porting this from a fold-switch to an INTERFACE-disruption switch (networks that bridge the binder-TNF interface rather than the binder core) is the key, lightly-explored transfer.

### Cited Findings
- A general strategy preorganizes His residues in buried hydrogen-bond networks to build homotrimers/heterodimers stable above pH 6.5 that undergo cooperative, large-scale conformational change when pH is lowered and electrostatic/steric repulsion builds as network His protonate — [Science 364:6441, aav7897 (2019)](https://www.ipd.uw.edu/publication-pdfs/4/517f778cf57bbe00f3bcd20be228243a/science.aav7897.pdf); [PMC7072037](https://pmc.ncbi.nlm.nih.gov/articles/PMC7072037).
- Transition-pH set point is the balance of two terms: intrinsic neutral-pH stability (increases with n = hydrophobic layers) vs the magnitude of pH-dependent destabilization (increases with m = His-network layers); replacing ONE His network with hydrophobic packing shifts the transition pH from above 5 down to ~3.5 — [PMC7072037](https://pmc.ncbi.nlm.nih.gov/articles/PMC7072037).
- Transition STEEPNESS (cooperativity) is set by the total number of protonatable His (increases with m); the one-His-network variant pRO-2.5 was "the least cooperative design tested"; each network held three His that protonate cooperatively, and protonation of one His disrupts the network, making the remaining His more accessible (positive cooperativity mechanism) — [PMC7072037](https://pmc.ncbi.nlm.nih.gov/articles/PMC7072037).
- Network LOCATION matters: His networks near chain termini gave higher transition pH and reduced cooperativity due to local flexibility allowing earlier protonation — [PMC7072037](https://pmc.ncbi.nlm.nih.gov/articles/PMC7072037).
- Precedent for anomalously SHARP single-His titration: 19F-NMR on fluorinated hemoglobin showed His beta-146 titrates abnormally sharply, especially with DPG, linked to cooperative ligand binding — [Huestis & Raftery, PNAS 1972 (Caltech)](https://resolver.caltech.edu/CaltechAUTHORS:HUEpnas72).

### Inferences
- OPERATIONAL HOOK: use RFdiffusion + soluble-ProteinMPNN to place a 3-His buried hydrogen-bond network that STRADDLES the binder-TNF interface (one or two His contributed by the binder, H-bonded to a polar TNF anchor such as Q67 or the backbone near E135/E146). At pH 7.4 the neutral network is a buried polar H-bond zipper holding the complex; at pH 6.0 collective protonation injects +2/+3 into a buried, low-dielectric pocket -> huge desolvation + mutual His-His repulsion -> interface breaks. Because the His sit in a low-dielectric buried environment, their effective pKa and the cooperativity both rise (burial sharpens the transition) — this is the single most promising "beyond-Baker" motif because it moves the proven cooperative network from fold-destabilization to interface-destabilization.
- Scoring the Hill coefficient: run OpenMM constant-pH MD (or a titration sweep with PROPKA/pKAI at a ladder of pH) on the complex, fit occupancy of the interface vs pH to a Hill equation, and SELECT on fitted n>1 and midpoint in [6.0, 6.5]. This gives a direct, rankable cooperativity metric rather than a single pKa.
- Proton-wire / tautomerism refinement: His has two neutral tautomers (HID/HIE) and the charged HIP; a buried His whose Ndelta/Nepsilon H-bond donor/acceptor roles are preorganized will titrate more sharply. Enforce a specific tautomer in the H-bond network (MPNN + explicit HID/HIE assignment) so the network is a cooperative unit, not independent sites. Independent His give n~1 (gradual); H-bond-COUPLED His give n>1 (sharp) — this is the mechanistic lever for the "sharp OFF at 6.0" prize.
- Two coupled His sharing one buried carboxylate (Asp/Glu) create an anti-cooperative-then-cooperative "proton pair": the first protonation is easy (neutralizes the acid), the second is hard until pH drops further, which can be tuned to place a steep step right at 6.0. This is classic coupled-pKa behavior and is designable with one buried acid + two His.
- FEASIBILITY (4 days): HIGH. The motif is a direct reuse of an existing, validated design principle; the pipeline (RFd + MPNN + His-network H-bond constraints + PROPKA/pKAI titration sweep + Hill fit) is all already in their stack. Main risk = getting a buried 3-His network to also be an INTERFACE contact without killing fold stability (the Baker networks were in the core).

### Gaps
- No source quantifies the maximum Hill n achieved (the paper reports relative cooperativity, not absolute n values for each design); the exact n needed for a "clean" OFF at 6.0 vs ON at 7.4 (likely n>=2 to get ~90% release over a 1.4-unit window) is an inference, not a cited number.
- Interface (vs core) His networks for binder RELEASE are not reported in the retrieved literature — this is the untried part and carries buildability risk.

## Q3: Conceptual DFT (Fukui functions, electrophilicity, hardness/softness) to design/score the His-partner contacts

### Takeaway
Conceptual-DFT descriptors (condensed Fukui functions f+/f-, local softness, electrophilicity omega = mu^2/2eta) have been used for ligand docking, active-site detection, and protein folding, but essentially NEVER as a side-chain-pair compatibility score in de novo binder design. They offer a cheap, principled way to rank which natural partner (carboxylate vs aromatic vs amide) best accepts a proton/charge from a switching His — an emerging-to-untried transfer.

### Cited Findings
- The Fukui function is a local reactivity descriptor (change in electron density on adding/removing an electron) related to local softness; its usefulness has been explored for three biological problems: ligand docking, active-site detection, and protein folding — [Conceptual DFT Fukui-function study of biological problems, PMC2848499](https://pmc.ncbi.nlm.nih.gov/articles/PMC2848499).
- Global electrophilicity omega = mu^2 / (2*eta) (mu = electronic chemical potential, eta = chemical hardness) is a reliable descriptor of biological activity; local philicity omega_k = omega * f_k partitions it to atoms — [conceptual DFT review, PMC9505550](https://pmc.ncbi.nlm.nih.gov/articles/PMC9505550); [Nature Index: chemical reactivity descriptors in DFT](https://www.nature.com/nature-index/topics/l4/chemical-reactivity-descriptors-in-density-functional-theory).
- Local electrophilicity governs charge transfer between a toxin and receptor; condensed-to-atom Fukui functions and local philicity are standard for ranking site reactivity — [PMC2848499](https://pmc.ncbi.nlm.nih.gov/articles/PMC2848499).

### Inferences
- OPERATIONAL HOOK (cheap, GFN2-xTB-compatible): for each candidate His + partner contact, extract the capped side-chain pair, run GFN2-xTB (which gives atomic charges and can be post-processed for condensed Fukui functions via the finite-difference N, N+1, N-1 densities), and compute: (i) hardness eta = (E_LUMO - E_HOMO)/2 of the partner, (ii) the nucleophilic Fukui f- on the proton-accepting atom. The "soft-soft / hard-hard" HSAB principle predicts a protonated-His (hard cation) prefers a HARD base (carboxylate O) for a strong, pH-switchable electrostatic bridge, while a neutral-His cation-pi prefers a SOFT aromatic (Trp > Phe). Rank partner identity by HSAB match to the DESIRED switch direction. This directly informs the choice their memory already wrestles with (salt bridge vs cation-pi; Trp>Phe>Tyr).
- Electrophilicity/chemical-potential equalization gives a scalar "charge-transfer propensity" per contact that complements MD d_occ and xtb dE_switch — use it as a tie-breaker, not a primary gate.
- FEASIBILITY (4 days): MEDIUM. Condensed Fukui via three xTB single-points per side-chain pair is a ~50-line Python wrapper; the risk is that CDFT descriptors are calibrated for covalent reactivity/charge-transfer, and their relevance to a purely electrostatic/H-bond proton switch is an analogy, not established — treat as a secondary ranking signal.
- NOVELTY: using Fukui/HSAB to pick the His PARTNER residue in a de novo switch appears genuinely untried; conceptual DFT in proteins has been docking/folding analysis, not generative side-chain selection.

### Gaps
- No source demonstrates CDFT descriptors predicting pH-switch strength or direction specifically; the HSAB-to-switch-direction mapping is a reasoned transfer, not a validated result.

## Q4: Preorganized electric fields (Warshel) and vibrational Stark effect (Boxer) to tune a protonation-dependent interface

### Takeaway
The electrostatic-preorganization framework (enzyme active sites exert large, oriented intrinsic electric fields that stabilize charge-separated states) is directly transferable: design the binder to project a strong oriented field onto the His-partner bond so that protonation is energetically sharply pH-coupled, and SCORE the field with APBS/MD. This is standard in enzyme design but essentially untried as a pH-SWITCH tuning criterion in binder design.

### Cited Findings
- Vibrational Stark effect (VSE) spectroscopy measured the intrinsic electric field inside ketosteroid isomerase; the field is nearly perfectly oriented to stabilize the transition state and is "preorganized," adjusting the substrate so minimal structural change occurs on activation — [Fried, Bagchi, Boxer; PMC7474534](https://pmc.ncbi.nlm.nih.gov/articles/PMC7474534); [Boxer lab](https://boxerlab.stanford.edu/_files/ugd/4006c9_bdbbac6fba0b497392ce44c88736a589.pdf).
- A direct proportionality exists between activation free energy and the electric field PROJECTED onto the bond undergoing charge reorganization — a quantitative, computable catalysis/stabilization metric — [PMC7474534](https://pmc.ncbi.nlm.nih.gov/articles/PMC7474534).
- Active-site fields deviate markedly from solvent fields (preorganization is the difference between enzyme electrostatics and bulk solvent) — [PMC7474534](https://pmc.ncbi.nlm.nih.gov/articles/PMC7474534).

### Inferences
- OPERATIONAL HOOK: treat His protonation at the interface as a "charge reorganization" along the His-N...H+ axis. Compute the electrostatic field (APBS gives the potential; take its gradient, or use OpenMM to sum Coulomb field from all binder atoms) projected onto the His imidazole N that gains the proton. A design that preorganizes a LARGE field OPPOSING protonation at pH 7.4 (field from nearby binder carboxylates/dipoles stabilizing neutral His) and whose field environment REVERSES or desolvates on partner rearrangement amplifies the pH sensitivity. Select on |field projection| and its sign at the switching atom. This is the microscopic "why" behind a sharp pKa and is cheaper than full QM.
- Helix-macrodipole trick (fully protein-encodable): point an alpha-helix N-terminus (positive end) or C-terminus (negative end) at the switching His to bias its pKa by preorganized field — a classic, natural-AA-only way to tune pKa by ~1 unit and sharpen the transition.
- FEASIBILITY (4 days): MEDIUM-HIGH. Field projection from APBS/OpenMM is straightforward Python; mapping a target field magnitude to a pKa shift quantitatively is the hard part (needs calibration against PROPKA/pKAI/xtb on a few designs). Use as a mechanistic selection/explanation layer.
- NOVELTY: field-projection scoring is standard in enzyme design (Warshel/Boxer) but applying it to tune a binder pH-switch is untried.

### Gaps
- No retrieved source maps electric-field magnitude to a His pKa shift in a protein interface; the field->pKa->switch-sharpness chain is a reasoned inference requiring in-house calibration.

## Q5: Supramolecular / host-guest cooperativity and allosteric-chemistry principles

### Takeaway
Positive cooperativity (the engine of a sharp Hill n>1 OFF-switch) is the protein analog of host-guest/allosteric cooperativity; the actionable transfer is to make protonation events COUPLED (one event changes the affinity of the next) rather than independent. The Baker His-network result IS this principle instantiated in protein; explicit supramolecular/MWC-KNF framing mainly clarifies how to design the coupling.

### Cited Findings
- In the Baker networks, protonation of one His "will likely disrupt the surrounding network, making the remaining His residues more accessible" — a textbook positive-cooperativity (induced) mechanism driving a sharp collective transition — [PMC7072037](https://pmc.ncbi.nlm.nih.gov/articles/PMC7072037).
- Cooperativity (steepness) scales with the number of coupled protonatable His; isolated His give a gradual transition — [PMC7072037](https://pmc.ncbi.nlm.nih.gov/articles/PMC7072037).

### Inferences
- OPERATIONAL HOOK: design the interface so binding is CONDITIONAL on a coupled set (e.g., two His + one shared carboxylate + one cation-pi Trp all in one contiguous patch). Model it as a 2-4 site coupled-equilibrium (classic Wyman linkage / MWC two-state), fit the pH-release curve, and select on the apparent Hill n. Coupling comes from shared electrostatics (two charges in one low-dielectric pocket repel -> their pKas split -> cooperative collapse) and shared H-bond partners.
- Allosteric transfer: instead of the switch living IN the contact, put it ADJACENT (a buried His network one layer behind the interface whose protonation shifts the backbone and peels the paratope) — an allosteric OFF-switch. Their memory already flags the monolithic-fold buildability problem (PivotLid NO-GO), so prefer coupling WITHIN the contact patch (electrostatic frustration) over a mechanical hinge.
- FEASIBILITY (4 days): HIGH for in-contact electrostatic coupling (reuses Q2 pipeline); LOW for mechanical allostery (buildability, per their PivotLid lesson).

### Gaps
- No protein-specific host-guest cooperativity citation beyond the Baker result was retrieved; the MWC/KNF framing is standard textbook theory applied by analogy (not separately cited here).

## Q6: Specific-ion / Hofmeister effects, water structuring, desolvation thermodynamics as switch levers

### Takeaway
Desolvation is a real, exploitable OFF-switch lever: burying a charge (protonated His) in a low-dielectric interface costs large desolvation energy, so a buried His that protonates at pH 6.0 pays a desolvation penalty that breaks the interface. Hofmeister/specific-ion effects are a robustness/environment caveat (the assay salt matters) more than a designable natural-AA lever.

### Cited Findings
- Hofmeister effects: effective interactions between charged/neutral objects in water depend on salt TYPE, not just concentration, affecting solubility, crystallization, and protein-protein interactions — [Ion-specific effects on phase transitions in protein solutions, arXiv:0901.1599](https://ar5iv.arxiv.org/html/0901.1599).
- Direct ion-backbone (amide) interactions are crucial to salt-induced protein destabilization; water structure near the interface and ion hydration/dispersion must be modeled to capture specific-ion behavior — [arXiv:0901.1599](https://ar5iv.arxiv.org/html/0901.1599); [specific ion effects review, PMC3694758](https://pmc.ncbi.nlm.nih.gov/articles/PMC3694758).

### Inferences
- OPERATIONAL HOOK (desolvation): compute the Born/PB desolvation penalty of the switching His in the BOUND (buried, low local dielectric) vs FREE (solvent-exposed) states via APBS; a large bound-state desolvation penalty for the PROTONATED His is exactly what makes protonation break the interface. Select for interface His that are buried in the complex (high desolvation on protonation) but solvent-accessible in the free binder (so the free binder doesn't pay the penalty and stays folded). This reframes the Baker "buried charged network" as a desolvation-driven OFF-switch and gives a clean APBS scoring criterion.
- The desolvation lever is MECHANISTICALLY the same one their own notes already circle (QM dE_switch, APBS continuum NET, the "buried-His desolvation vs salt-bridge strength" frontier) — Hofmeister theory just names why continuum-NET is the right probe and why ECC/fixed-charge MD mis-weights it.
- Hofmeister as robustness check: validate the switch at physiological 150 mM NaCl and note the wet-lab buffer; a switch that depends on specific-ion water structuring is fragile. Not a primary design lever for natural-AA protein design.
- FEASIBILITY: HIGH for the desolvation APBS criterion (already in-stack); specific-ion modeling is NOT worth building in 4 days.

### Gaps
- No source quantifies a target desolvation penalty (kcal/mol) that guarantees interface rupture; must be calibrated in-house against the APBS dG_elec sign-flip (Q1).

## Q7: Marcus theory / proton-coupled electron-transfer (PCET) analogies for proton-coupled binding

### Takeaway
PCET/Marcus theory is a poor literal fit (no electron transfer in a pH-switch binder) but offers one genuinely useful conceptual import: treat binding + protonation as COUPLED along two coordinates (a binding/geometric coordinate and a proton coordinate), and design so the two are tightly coupled (proton transfer is contingent on the bound geometry). This "reorganization energy" framing is a lens, not a scoring function here.

### Cited Findings
- PCET theory extends Marcus electron-transfer theory: the proton is treated quantum-mechanically and reactions proceed via nonadiabatic transitions between mixed electron-proton VIBRONIC states, described by reorganization of the solvent/environment — [Hammes-Schiffer PCET theory, PMC2783471](https://pmc.ncbi.nlm.nih.gov/articles/PMC2783471); [PMC3005854](https://pmc.ncbi.nlm.nih.gov/articles/PMC3005854).
- Marcus theory represents reactant/product states as harmonic free-energy curves; driving force connects to rate via the Marcus relation — [PMC2783471](https://pmc.ncbi.nlm.nih.gov/articles/PMC2783471).

### Inferences
- OPERATIONAL HOOK (weak): model the switch as two coupled free-energy parabolas (bound-neutral-His vs released-protonated-His) along a reaction coordinate; the "reorganization energy" = structural rearrangement cost on protonation. A design with LOW reorganization energy (protonation directly destabilizes the bound state without needing a big conformational move) gives a sharper, more reliable switch. This can be probed by a short OpenMM steered-MD / umbrella sampling along the His-protonation + interface-separation coordinate, but that is expensive.
- FEASIBILITY (4 days): LOW as a quantitative method; MEDIUM as a design heuristic ("make protonation directly repulsive in the bound geometry, minimize required rearrangement"). Not recommended to build PCET machinery.
- NOVELTY: literal PCET is untried AND largely inapplicable; the two-coordinate coupling heuristic is a useful framing but not a differentiated lever.

### Gaps
- No evidence of PCET/Marcus formalism applied to pH-switch protein binding; this is an analogy with limited actionable payoff — flag as low priority for the report.

## Q8: QM/MM or constant-pH / xtb-based pKa and free-energy methods beyond PROPKA

### Takeaway
Two upgrades over single-structure PROPKA are feasible in 4 days and directly improve switch scoring: (1) constant-pH MD (CpHMD) to get titration curves + cooperativity under conformational sampling, and (2) GFN2-xTB microstate/tautomer pKa (ALPB water) as an independent, physics-based pKa leg. Both are stronger than PROPKA for buried, coupled His networks where PROPKA is known to be the harshest/least reliable proxy.

### Cited Findings
- Constant-pH MD requires no a priori assumption of which residues titrate or their states; discrete CpHMD couples MD with Monte-Carlo protonation moves, continuous CpHMD assigns each site a lambda in [0,1] evolving dynamically; a pH-based replica-exchange protocol accelerates protonation sampling — [A Guide to Continuous Constant-pH MD in Amber and CHARMM, LiveCoMS v4i1e1563](https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1563); [Dobrev 2017 JCTC](https://www.em.mpg.de/634541/Dobrev_2017_JCTC.pdf).
- CpHMD predicts pKa values AND proton-coupled conformational dynamics for enzymes/transporters — i.e., it natively captures the cooperative coupling a switch needs — [LiveCoMS guide](https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1563).
- GFN2-xTB first-principles pKa: enumerate protonation microstates/tautomers, optimize in gas phase + ALPB implicit water, derive macroscopic pKa from charge-state partition functions; calibrated SAMPL6 MAE 1.74 / RMSE 2.18 pKa units — [SAMPL6 first-principles pKa](https://www.k-dense.ai/use-cases/sampl6-first-principles-pka-prediction).
- Implicit-solvent GFN2-xTB pKa can reach r^2 ~ 0.95 on a congeneric set (C-PCM); accuracy ceiling is GFN2's limited accuracy for relative stability of protonation-site isomers — [k-dense SAMPL6 report](https://www.k-dense.ai/use-cases/sampl6-first-principles-pka-prediction); [RSC d0cp05281g](https://pubs.rsc.org/en/content/articlehtml/2021/cp/d0cp05281g); newer g-xTB improves pKa — [Rowan g-xTB pKa](https://rowansci.substack.com/p/g-xtb-pka-and-website-redesign).

### Inferences
- OPERATIONAL HOOK #1 (CpHMD, HIGH value): run OpenMM-based constant-pH MD (or Amber/CHARMM CpHMD) on the top complexes across a pH ladder (5.5-7.5). Outputs give (i) titration curves -> capacitance (Q1), (ii) Hill n of interface disruption (Q2), (iii) whether the switch survives conformational averaging (not just a single relaxed structure). This is the single best method upgrade: it is the ONLY approach that natively yields the cooperative, multi-His, conformationally-averaged OFF-transition that is the actual objective. OpenMM constant-pH exists; budget 1-2 days to stand up.
- OPERATIONAL HOOK #2 (xtb microstate pKa, HIGH value, cheap): for each switching-His cluster, extract the pocket, enumerate HID/HIE/HIP microstates (and the partner's protonation), optimize with GFN2-xTB --alpb water, and compute macroscopic pKa from the microstate partition function. This is their existing qm_switch_fair.py philosophy extended from dE_switch to an actual pKa, giving a 4th/5th independent pKa leg (vs PROPKA/pKAI/JustHISpKa) that is burial-aware and coupling-aware. Calibrate the systematic offset on a few residues with known pKa.
- Combine: xtb/CpHMD give the midpoint and sharpness; APBS (Q1) gives the force sign; the switch passes only if ALL agree on midpoint in [6.0,6.6] and n>1. This multi-method consensus is exactly the "triple-confirm" discipline in their memory, upgraded for cooperativity.
- FEASIBILITY: xtb microstate pKa = HIGH (extends existing script). CpHMD = MEDIUM (setup cost, but the payoff is the cooperativity metric nothing else gives). Both beat single-structure PROPKA for buried/coupled His.

### Gaps
- GFN2-xTB's ~1.7-2.2 pKa-unit MAE is marginal for resolving a midpoint to the ~0.3-unit precision a sharp 6.0 switch needs; it is a consensus leg, not a sole oracle (consistent with their buried-acid pKa NO-GO lesson). Microstate pKa for strongly coupled multi-His networks (where macroscopic pKa is not a simple sum) may need a full linkage/titration treatment, not independent-site xtb.

## Q9: Bioisostere / chemical-robustness logic for cross-reactivity (human->mouse, 79% id)

### Takeaway
Cross-reactivity by chemistry = anchor the paratope on RESIDUES AND INTERACTION TYPES that survive the human->mouse substitutions. The actionable logic: (a) map which of K65/Q67/E135/E146 and neighbors are conserved vs substituted human->mouse, (b) design contacts to conserved positions, and (c) where a position differs, make the contact tolerant via conservative-substitution/bioisostere reasoning (contact the backbone, or use an interaction robust to the specific side-chain swap).

### Cited Findings
- Antibodies can be cross-reactive across species (e.g., tissue-factor antibodies binding both human and mouse TF); amino-acid substitutions at cross-reactive epitope residues can reduce/ablate cross-reactivity — i.e., cross-reactivity hinges on contacting positions that are conserved OR substitution-tolerant — [US Patent 11447566, anti-TF antibodies](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11447566).
- Conservative-substitution groups (side chains of similar chemistry, generally functionally interchangeable): Val-Leu-Ile, Phe-Tyr, Lys-Arg, Ala-Val, Glu-Asp, Asn-Gln — [US Patent (conservative substitution definition)](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6783936).

### Inferences
- OPERATIONAL HOOK (do this FIRST, it's cheap and gates everything): align human (1TNF) vs mouse TNF-alpha sequences, map the 21% differences onto the receptor-binding groove surface, and CLASSIFY each anchor contact: (1) identical -> safe to make a specific side-chain contact; (2) conservative swap (e.g., Glu<->Asp, Lys<->Arg) -> design an interaction robust to it (a salt bridge to a carboxylate works for both Glu and Asp; an H-bond/charge contact to a cation works for both Lys and Arg); (3) non-conservative swap -> AVOID the side chain, contact the BACKBONE (main-chain H-bonds are species-invariant) or contact a nearby conserved residue instead. Bake this into the RFdiffusion hotspot selection and an MPNN interface-residue filter.
- "Bioisostere" in a natural-AA protein = choose the interaction GEOMETRY that both human and mouse side chains satisfy: carboxylate bidentate contacts satisfy Asp or Glu; a pi/cation contact tolerant to Lys/Arg; a hydrophobic packing tolerant to Val/Leu/Ile. Score each design by "fraction of interface energy contributed by species-invariant contacts" (backbone + conserved + conservative-tolerant), and select for high invariant fraction.
- Cross-check with the co-fold delta they already use (per memory): fold the SAME binder sequence onto human and mouse TNF (Boltz-2/Chai-1/AF2-ig) and require low ipTM/pae delta — but drive the DESIGN with the conservation map, not just post-hoc folding.
- IMPORTANT tension for the report: the pH-switch (Q1-Q6) wants to EXPLOIT a specific charged anchor (e.g., K65), while cross-reactivity wants robustness to substitution. CHECK that the switch-partner anchor is conserved human->mouse; if K65/the switch anchor differs, the switch and cross-reactivity objectives collide and the switch must be rebuilt on a conserved anchor. This is a concrete go/no-go the report must surface.
- FEASIBILITY: HIGH. Sequence alignment + contact classification + invariant-fraction scoring is a day of Python, reusing PyRosetta per-residue interface energies.

### Gaps
- The human-vs-mouse conservation status of the SPECIFIC anchors K65/Q67/E135/E146 and the switch partner was not verified in this research pass (no sequence alignment run); this must be confirmed from UniProt/1TNF vs mouse TNF before committing switch chemistry — flagged as the top prerequisite.
- Peer-reviewed de novo "designed cross-reactivity by conserved-contact selection" precedent was not found (patent/antibody literature only); the invariant-fraction score is a reasoned transfer, not a cited validated method.

## Cross-cutting synthesis for the report writer (priorities)
- HIGHEST payoff, lowest risk, beyond-Baker: (1) interface-straddling COOPERATIVE 3-His hydrogen-bond network as a desolvation+repulsion OFF-switch (Q2+Q6); (2) APBS dG_elec SIGN-FLIP + capacitance-peak score for pH-6.0 repulsion (Q1); (3) CpHMD + xtb-microstate pKa to measure midpoint AND Hill n (Q8); (4) conservation-map-driven invariant-contact design for cross-reactivity, run FIRST as a gate (Q9).
- MEDIUM, mechanistic/selection-layer: preorganized electric-field projection + helix-macrodipole pKa tuning (Q4); Fukui/HSAB partner selection (Q3); coupled-equilibrium/Wyman cooperativity modeling (Q5).
- LOW / deprioritize: literal PCET-Marcus (Q7); explicit Hofmeister specific-ion modeling (Q6, keep only as a 150 mM robustness check).
- The genuinely UNTRIED transfers (the "bold" asks): charge-capacitance-peak-at-pH-6.5 as a design objective; APBS force-SIGN-flip as the switch gate; Fukui/HSAB side-chain-partner selection; electric-field-projection as pH-switch tuning; interface (not core) cooperative His networks for binder release. The cooperative-His-network and desolvation ideas are proven in protein (emerging), the rest are abstract-chemistry imports that have not been used in de novo binder design.
