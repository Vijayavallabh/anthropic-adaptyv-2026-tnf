# Abstract / pure-mathematics concepts for a pH-switchable, cross-reactive TNF-alpha binder pipeline

Scope: a catalog of mathematical frameworks under-used (or never used) in de novo binder design, each with (a) core idea, (b) a concrete pipeline hook on *our* toolchain outputs (sequences, PDBs, PROPKA/pKAI pKas, OpenMM/xTB trajectories, AF2-ig/Boltz-2/Chai-1 scores, foldseek, the ~375-design screen, the <=20 ranked submission), (c) 4-day Python feasibility, (d) untried-vs-emerging verdict. Priority order throughout matches the objective: pH-switch sharpness (sharp OFF at pH 6.0) > mouse cross-reactivity > affinity.

Convention below: **Cited Findings** = sourced facts; **Inferences** = my pipeline hook + feasibility + novelty verdict (my synthesis, not sourced claims).

---

## Q1. Topological data analysis / persistent homology (PH)

### Takeaway
PH turns a structure or an MD ensemble into a multiscale topological fingerprint (persistence diagram/landscape/image) that is coordinate- and alignment-free and stable to noise. It is *emerging* in protein-ligand ML but *essentially untried* for the two specific jobs we need: (i) a two-state (pH 7.4 vs 6.0) ensemble-separation score, and (ii) an interface shape-complementarity/robustness score for a de novo binder.

### Cited Findings
- PH builds a filtration (e.g. Vietoris-Rips or alpha complex over atom coordinates at growing radius) and records when topological features (components, loops, voids = Betti_0/1/2) are born and die; features that persist across a wide scale range are real signal, not sampling noise — [PH-based ML survey, arXiv:1811.00252](https://arxiv.org/pdf/1811.00252); [Weighted PH for biomolecular data, arXiv:1903.02890](https://arxiv.org/pdf/1903.02890).
- Element-specific / "topological fingerprint" PH coupled to ML gave top performance on protein-ligand binding-affinity and virtual-screening blind challenges (Cang & Wei) — [Topological fingerprints reveal protein-ligand binding mechanism, arXiv:1703.10982](https://arxiv.org/pdf/1703.10982); [Representability of algebraic topology for biomolecules in ML scoring/VS, arXiv:1708.08135](https://arxiv.org/pdf/1708.08135).
- "The majority of active-site and allosteric-pathway residues lie near the most persistent loop in the filtered Vietoris-Rips complex" — PH captures binding-relevant dynamics — [PH + dynamical distances for protein binding, arXiv:1412.1394](https://arxiv.org/pdf/1412.1394).
- Persistence *landscapes* detected the open vs closed conformational change of maltose-binding protein directly from structure/MD — i.e. PH discriminates two functional states of the *same* protein — [arXiv:1412.1394](https://arxiv.org/pdf/1412.1394).
- PH features from MD reveal stable, misfolded, and transition states and compress long-MD data to low-dimensional kinetic embeddings — [Learning Topological Representations for MD, arXiv:2506.14737](https://arxiv.org/pdf/2606.14737); [Protein-folding analysis using PH features, arXiv:1905.05942 / PMC7300307](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7300307/).
- PH connects to the geometry/topology of conformational spaces and energy landscapes — [Topology & geometry of molecular conformational spaces and energy landscapes, arXiv:1907.07770](https://arxiv.org/pdf/1907.07770).

### Inferences
- **Pipeline hook #1 (TOP PRIORITY — pH switch):** From OpenMM trajectories of the bound complex at pH 7.4 (His neutral, HIE/HID) and pH 6.0 (His protonated, HIP), compute the persistence diagram of the *interface atoms only* (both protomers' epitope shell + binder paratope) per frame. Define a **topological switch score = Wasserstein/bottleneck distance between the pH-7.4 and pH-6.0 persistence-diagram distributions** (Betti_0 = contact-cluster fragmentation; Betti_1/2 = pocket loops/voids opening as the interface disassembles). This is a *protonation-state ensemble summary* that is independent of RMSD alignment and complementary to the existing MD d_occ gate — it should catch "the interface topologically falls apart at 6.0" even when a single salt bridge metric is ambiguous.
- **Pipeline hook #2 (shape complementarity / robustness):** Alpha-complex PH of the paired epitope groove (between adjacent protomers, where TNFR1/TNFR2 and adalimumab bind) gives a pocket-void fingerprint; score binder paratope-vs-pocket shape complementarity as similarity of their PH fingerprints (fills voids = low surviving Betti_2). Binding-mode robustness = bottleneck stability of the diagram across MD frames (a persistent-across-frames fingerprint = a stable binding mode).
- **Pipeline hook #3 (selection):** vectorize persistence diagrams as persistence images and use them as cheap features to rank/cluster the 375-design screen.
- **Feasibility (4 days): LOW RISK.** GUDHI and Ripser are pip-installable; alpha/Rips persistence on a few-hundred-atom interface is milliseconds/frame; `persim`/`gudhi` give persistence images, landscapes, and bottleneck/Wasserstein diagram distances out of the box. The only real work is choosing the atom subset and the filtration.
- **Novelty verdict: EMERGING for affinity ML, UNTRIED for pH-two-state switch scoring and for de novo binder interface robustness.** The maltose-binding-protein two-state result (1412.1394) is the proof-of-concept that it *can* separate conformational states; applying it to a *designed protonation switch* is the novel transfer.

### Gaps
- No source directly demonstrates PH applied to a *protonation-state* (pH-titrated) ensemble; the two-state evidence is for ligand-driven/solvent-driven conformational change. The transfer is well-motivated but unproven.

---

## Q2. Optimal transport / Wasserstein geometry (+ Gromov-Wasserstein, + DPP for diversity)

### Takeaway
OT gives a principled distance between *distributions* (two conformational ensembles, two oracles' output distributions, a design portfolio). Three distinct, implementable hooks: Wasserstein distance between pH-7.4/6.0 ensembles (switch magnitude), **Gromov-Wasserstein** for alignment-free cross-species and two-state structure comparison, and **Wasserstein-barycenter ensembling** of folding oracles. OT-for-ensemble-comparison is *emerging* (IDP field); GW for cross-species epitope matching and W-barycenter oracle fusion are *untried* in binder design. For the <=20 portfolio, the right tool is a **determinantal point process (DPP)** (quality x diversity), not OT.

### Cited Findings
- Wasserstein distance is the natural metric on probability distributions from OT; a conformational ensemble can be treated as a probability distribution and two ensembles compared by Wasserstein distance — WASCO compares intrinsically-disordered-protein ensembles as distributions on R^3 (global) and the 2D flat torus (local, backbone dihedrals) — [WASCO, HAL hal-04056662](https://hal.laas.fr/hal-04056662).
- A root-mean Wasserstein distance (RMWD) has been proposed that generalizes all-atom RMSD to an ensemble-vs-ensemble distributional measure of positional accuracy (reported in the conformational-ensemble-generation literature) — [TEMPO ensemble generation, arXiv:2511.05510](https://arxiv.org/pdf/2511.05510).
- **Gromov-Wasserstein (GW)** compares two structures by their *internal* residue-residue distance matrices, minimizing distortion of intra-molecular distances; it is invariant to rigid transforms (no alignment needed), approximable in near-linear time, and GWProt uses it to "identify conformational switches in individual proteins," "detect functional domains shared among evolutionarily distant proteins," and reveal topological rearrangements — [GWProt, bioRxiv 2025.11.03.686256](https://www.biorxiv.org/content/10.1101/2025.11.03.686256.full.pdf); [GWProt docs](https://gwprot.readthedocs.io/).
- **Wasserstein-barycenter model ensembling** finds the OT consensus of several models' output distributions; it balances confidence and semantic agreement better than arithmetic/geometric mean of outputs — [Wasserstein Barycenter Model Ensembling, ICLR 2019, arXiv:1902.04999](https://arxiv.org/pdf/1902.04999).
- **DPPs** select subsets that are simultaneously high-quality and diverse via a repulsive point process; greedy MAP approximations with guarantees beat ad-hoc quality+diversity heuristics for batch selection/active learning — [Batch Active Learning Using DPPs, arXiv:1906.07975](https://arxiv.org/pdf/1906.07975).

### Inferences
- **Hook #1 (pH switch, alignment-free):** GW distance between the residue-distance matrix of the interface at pH 7.4 vs 6.0 = a single scalar "how much did the interface geometry rearrange" that needs no superposition and is robust to rigid-body drift of one protomer. GWProt literally advertises "identify conformational switches" — this is the closest published match to our exact need.
- **Hook #2 (cross-reactivity, TOP-2 priority):** GW between the human-TNF-complex interface and the mouse-TNF-complex interface (same binder, 79% identical target) measures *epitope-geometry conservation under the binding mode*; low GW distortion at the conserved anchors (K65/Q67/E135/E146) = the binder sees a geometrically-consistent epitope across species. Also usable to pre-screen which conserved-anchor subset is most geometrically invariant human->mouse before committing designs.
- **Hook #3 (oracle fusion):** AF2-ig, Boltz-2, Chai-1 disagree; instead of averaging ipTM/pae, treat each oracle's per-residue contact/PAE map or pose distribution as a measure and take the **Wasserstein barycenter** as the consensus pose/score. More robust to one outlier oracle than the mean (ties to the memory lesson that a single co-folder's cross-species ipTM is folder pessimism).
- **Hook #4 (the <=20 submission):** DPP over the PASS designs with quality = objective score (pH-selectivity > cross-react > affinity composite) and a similarity kernel from sequence/fold/epitope features; greedy MAP gives a diverse-yet-strong 20. Directly operationalizes the memory directive "pick for diversity, not cousins."
- **Feasibility (4 days): LOW RISK.** `POT` (Python Optimal Transport) provides `emd2`/Sinkhorn (Wasserstein), `gromov_wasserstein2`, and free-support `barycenter`. GW on ~50-residue interface distance matrices is seconds. DPP: `dppy`, or a ~15-line greedy Cholesky MAP. All pip-installable, all CPU.
- **Novelty verdict:** Wasserstein ensemble comparison = EMERGING (IDP field, WASCO/RMWD). GW for *cross-species epitope consistency* and for *pH two-state switch magnitude* = UNTRIED in binder design (GWProt is 2025, protein-structure-general, not binder-design). W-barycenter oracle fusion = UNTRIED in protein design. DPP for binder-portfolio selection = UNTRIED here (standard elsewhere).

### Gaps
- Exact primary citation for RMWD (root-mean Wasserstein) not pinned — it appears in the ensemble-generation literature (AlphaFlow/str2str lineage) but I could not confirm the originating paper within the search budget.

---

## Q3. Dynamical systems / bifurcation & catastrophe theory (the SHARP pH-6.0 OFF transition)

### Takeaway
A sharp cooperative OFF transition is mathematically a steep sigmoid; the right abstract lens is a **fold (saddle-node) bifurcation / cusp catastrophe as "organizing center"** plus thermodynamic linkage (binding polynomial) to compute the Hill slope. The abstract framing (design *toward* a bifurcation to maximize slope) is untried in binder design, but the *operational* hook (maximize the Hill coefficient of the pH-titration of the interface, pin the midpoint at 6.0) is cheap and directly buildable from PROPKA/pKAI output.

### Cited Findings
- The **cusp** is the codimension-2 "organising center" governing switch/bistable regimes; its normal form is dx/dt = a + b·x - x^3, and it is the intersection of a (super)critical pitchfork and a fold (saddle-node) bifurcation — [Monostability and bistability of biological switches, arXiv:2104.04227](https://arxiv.org/pdf/2104.04227).
- The **Hill coefficient** quantifies cooperative multi-ligand binding; steeper (higher-n) responses require multiple cooperating sites; n>=2 is the minimum that can induce bistability in a feedback circuit — [arXiv:2104.04227](https://arxiv.org/pdf/2104.04227).
- Multiple buried ionizable residues give steeper pH sensitivity over a narrower pH range than a single residue; de novo histidine hydrogen-bond-network proteins were designed to be stable above pH 6.5 and undergo **cooperative** conformational change when His protonate at lower pH, with transition pH and cooperativity *tunable* via the number of His networks and surrounding hydrophobicity — [Baker lab, De novo design of tunable pH-driven conformational changes (Janelia listing)](https://www.janelia.org/publication/de-novo-design-of-tunable-ph-driven-conformational-changes).
- Inserting charge into a hydrophobic interior is an established strategy for engineering pH-sensitive switches (anomalous buried pKa) — [US patent 10138471, pH-sensitive switches via buried charge](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10138471).

### Inferences
- **The right model for OUR switch:** our pH switch is a *monotonic* sharp OFF transition near pH 6.0, not a hysteretic toggle, so the relevant object is the **fold/threshold** and the **Hill slope of the titration curve**, not a full bistable cusp with hysteresis (no positive feedback loop in a static binder). Concretely: from PROPKA/pKAI pKas of the complex, target-alone, and binder-alone, build the Wyman-linkage proton-binding curve Delta_n_H(pH) and the implied binding free energy DeltaG_bind(pH); fit a Hill function; **two selection criteria: (i) midpoint pH_1/2 as close to 6.0 as possible, (ii) maximize Hill slope n (the "sharpness").**
- **Bifurcation as a design objective (the bold transfer):** the Hill slope rises as the system approaches a fold/criticality; stacking *multiple coupled* His at the interface (cooperatively coupled protonation, like the Baker His-network cooperativity) pushes the effective n up. So the abstract recommendation is: **select/engineer for >=2 coupled His whose protonation jointly collapses the interface, and rank designs by proximity-to-fold (steepness), not just by sign of DeltaG_pH.** This upgrades the existing memory lesson "far-better needs Hill n>1" into an explicit, scored objective.
- **Feasibility (4 days): LOW RISK for the Hill-slope objective** (it is post-processing of pKa outputs already in the pipeline + a `scipy.optimize.curve_fit`). HIGHER RISK / optional for the full dynamical treatment (building an ODE/transition model and locating the bifurcation numerically) — not needed for the deadline.
- **Novelty verdict: UNTRIED as an explicit design objective in de novo binder design.** Cooperativity/Hill is standard biochemistry and the Baker His-network cooperativity is published for *standalone switch proteins*, but framing *binder-interface* pH-off-sharpness as proximity-to-a-fold-bifurcation and selecting on it is novel.

### Gaps
- No source directly designs a *binder-interface* pH switch by bifurcation analysis; evidence is for standalone switch proteins (Baker) and for generic bistable circuits. The transfer is sound but unvalidated for binders.
- Could not pin the exact journal/year of the Baker tunable-pH paper beyond the Janelia listing within budget.

---

## Q4. Information theory (two-state discrimination + few-shot selection)

### Takeaway
The bind@7.4-vs-not@6.0 problem is literally a **1-bit channel-capacity / discrimination** problem, and few-shot selection is an **information-gain / rate-distortion** problem. Information theory is mature in cell-signaling fidelity (channel capacity of signaling networks) and *emerging* in protein design (information-bottleneck generators), but using channel capacity as a *pH-switch discriminability score* and info-gain as the *<=20 acquisition rule* is untried here.

### Cited Findings
- Ligand-receptor binding can be modeled as a finite-state Markov channel (the "BIND channel") with computable capacity for binary output/binary state — [Thomas & Eckford, IEEE Trans. Inf. Theory 2016](https://case.edu/math/thomas/ThomasEckford2016IEEE-Trans-IT.pdf).
- Channel capacity of a signaling network measures its fidelity in distinguishing inputs; capacity = 1 bit means the system can at best distinguish two input levels (e.g. present/absent) — many mammalian pathways sit near ~1 bit — [eLife reviewed-preprint 87747, sensing heterogeneously distributed in cell populations](https://elifesciences.org/reviewed-preprints/87747v1); capacity-achieving input distribution characterization — [Thomas & Eckford 2016](https://case.edu/math/thomas/ThomasEckford2016IEEE-Trans-IT.pdf).
- The **information bottleneck** extracts a minimal representation that is maximally predictive of a target by trading off I(X;Z) down and I(Z;Y) up; PAC-Bayesian IB quantifies per-sample information density to cope with scarce protein-ligand complex data — [IBEX, arXiv:2508.10775](https://arxiv.org/pdf/2508.10775).
- Concept-bottleneck protein language models (CB-pLM) enforce an interpretable information bottleneck and have been scaled to 3B parameters for generative protein design — [CB-pLM, arXiv:2411.06090](https://arxiv.org/pdf/2411.06090).

### Inferences
- **Hook #1 (pH switch as a channel):** treat the pH state (7.4 vs 6.0) as a binary channel input and a binding readout (predicted DeltaG_bind, ipTM, or a contact statistic) as the output; estimate the **mutual information / 1-bit capacity** between pH-state and predicted-binding across the MD/oracle ensemble. A perfect switch approaches 1 bit; a constitutive binder has ~0 bits. This is a *single principled scalar* that fuses "binds at 7.4" and "doesn't at 6.0" into the exact quantity we care about, and it penalizes designs whose two-state distributions overlap (noise-aware, unlike a mean-difference).
- **Hook #2 (cross-reactivity):** extend to a joint channel over {species} x {pH} and require high conditional MI for binding on both human and mouse at 7.4 — an information-theoretic encoding of "same sequence, both species."
- **Hook #3 (few-shot selection, the ~375 screen -> <=20):** pick designs by **expected information gain / mutual information** about the top-objective (max-entropy / BALD-style acquisition), which naturally balances informativeness and reduces redundancy; complements the DPP diversity selection (Q2).
- **Feasibility (4 days): LOW-MEDIUM.** MI/capacity between two 1-D score distributions is a Blahut-Arimoto or a k-NN MI estimator (`sklearn`/`npeet`), a few lines. The IB/CB-pLM *generators* are NOT 4-day builds (skip). The scoring/selection hooks are the buildable part.
- **Novelty verdict:** channel capacity for cell signaling = MATURE; IB *generators* for protein/molecule design = EMERGING; **channel capacity as a pH-switch discriminability score and info-gain as the binder-selection acquisition = UNTRIED** in this setting.

### Gaps
- No source applies channel capacity to a *designed* two-state binder; the capacity literature is descriptive (natural signaling), not a design objective. Classic Cheong et al. Science 2011 ("information transduction capacity of noisy biochemical signaling networks", the ~1-bit result) was not fetched directly within budget.

---

## Q5. Morse theory / discrete Morse theory (switch-relevant saddle points)

### Takeaway
A transition state between the bound and the pH-disassembled state is, by definition, an **index-1 saddle** (Hessian with exactly one negative eigenvalue) — pure Morse theory. Sublevel-set persistence / merge trees summarize the energy-vs-pH landscape cheaply; high-index saddle dynamics can *find* the saddles. Morse theory is standard in energy-landscape physics but untried as a *binder pH-switch transition-barrier* score.

### Cited Findings
- On a potential energy surface the key points are minima and the transition states joining them; a transition state is a stationary point whose Hessian has exactly one negative eigenvalue, whose eigenvector is the reaction coordinate; the **index** = number of negative Hessian eigenvalues — [energy-landscape / saddle search results summarized in arXiv:cond-mat/9904304 (energy landscape of a model protein)](https://arxiv.org/pdf/cond-mat/9904304).
- Minima linked by higher-index saddles are also linked by lower-energy true (index-1) transition states — [arXiv:2002.10690](https://arxiv.org/pdf/2002.10690).
- **High-index saddle dynamics (HiSD)** systematically computes saddle points of any index and maps the whole "solution landscape" (the connectivity of minima and saddles) — [Searching the solution landscape by generalized HiSD, arXiv:2002.10690](https://arxiv.org/pdf/2002.10690); [Improved HiSD, arXiv:2502.03694](https://arxiv.org/pdf/2502.03694).
- The topology/geometry of molecular energy landscapes (Morse-theoretic structure: minima, saddles, basins) can be analyzed directly — [arXiv:1907.07770](https://arxiv.org/pdf/1907.07770).

### Inferences
- **Hook #1 (pH switch barrier, TOP priority):** the quantity that makes an OFF transition *sharp and complete* is a low barrier from bound-protonated(HIP) to unbound and a high barrier the other way near pH 6.0. Run a string/NEB or steered path in OpenMM (or along a GW/PH reaction coordinate from Q1/Q2) between the bound@7.4 (HIE/HID) and bound@6.0 (HIP) states; the **index-1 saddle energy = the switch transition barrier.** Rank designs by a *low* barrier to disassembly at 6.0 and a *high* barrier at 7.4 — a Morse-theoretic complement to the thermodynamic DeltaG_pH (kinetics, not just equilibrium).
- **Hook #2 (cheap landscape summary):** build the scalar field "interface interaction energy as a function of pH (and of a 1-D disassembly coordinate)" from a pH-scan of xTB/MD energies; compute its **merge tree / sublevel-set persistence** (1-D discrete Morse) — the number and depth of basins and the saddle between them *is* the switch cooperativity signature (one deep basin that collapses = clean switch; two shallow basins = leaky). This is minutes of compute.
- **Feasibility (4 days): MEDIUM.** Index-1 Hessian eigenvalue check and 1-D merge-tree/sublevel persistence are easy (`numpy`, `gudhi`). A full NEB/HiSD transition-state search in OpenMM is doable but is the highest-effort item here; the 1-D energy-landscape merge-tree is the lazy, high-value version. HiSD reference code exists but integrating it to all-atom in 4 days is risky.
- **Novelty verdict:** Morse/saddle analysis of energy landscapes = MATURE in biophysics; **applying it as a scored pH-switch *transition-barrier* selection criterion for a de novo binder = UNTRIED.**

### Gaps
- No source computes a transition-state barrier specifically for a *protonation-driven binder-interface disassembly*; evidence is for folding/landscape physics generally.

---

## Q6. Robust / distributionally-robust optimization (DRO) for cross-reactivity

### Takeaway
Cross-reactivity (same sequence must bind human AND mouse TNF, 79% identical) is exactly a **robustness-to-distribution-shift** problem: optimize worst-case binding over an ambiguity set containing both species (and perturbations). Wasserstein-DRO is a mature OR/ML framework with out-of-sample guarantees but — as best I can find — has **never been applied to protein/binder design**, making this one of the boldest genuinely-untried transfers.

### Cited Findings
- DRO minimizes the worst-case risk over an **ambiguity set** of distributions; the **Wasserstein ambiguity set** is a ball of radius epsilon (in Wasserstein distance) around the empirical distribution, contains both continuous and discrete distributions, and yields finite-sample out-of-sample performance guarantees; epsilon tunes conservativeness (larger epsilon = more robust to sampling error) — [Data-driven DR-MPC using the Wasserstein metric, arXiv:2105.08414](https://arxiv.org/pdf/2105.08414); [Data-driven chance constraints under Wasserstein ambiguity, arXiv:1805.06729](https://arxiv.org/pdf/1805.06729).
- Protein surfaces (where epitopes/paratopes live) are more mutationally robust than cores; such robustness principles are proposed as useful for designing sequences robust to point mutations — [BMC Biology s12915-020-00870-9](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-020-00870-9).

### Inferences
- **Hook (cross-reactivity, TOP-2 priority):** define the ambiguity set as the finite set {human-TNF complex, mouse-TNF complex, + a few in-silico perturbations: AF2/Boltz samples, point-mutant targets at non-conserved positions, side-chain repacks}. **Select/rank binders by their worst-case (min over the set) predicted affinity/ipTM at pH 7.4**, i.e. a max-min (robust) objective instead of the human-only score. This is a *tractable sampled-ambiguity-set DRO* — no DRO solver needed: just take the min over the sampled targets as the score, which is exactly the empirical worst-case risk. A Wasserstein radius epsilon can be emulated by how aggressively you perturb the non-conserved positions.
- **Second-order use:** the DRO view says a binder anchored on the *conserved* residues (K65/Q67/E135/E146) is the one whose score is flat across the ambiguity ball — so DRO formally justifies the existing "design against conserved anchors" heuristic and gives a number (worst-case score) to rank on.
- **Feasibility (4 days): LOW RISK** as a *selection criterion* (min over an ensemble of targets = a few lines around existing oracle calls). Implementing a genuine Wasserstein-DRO *training/generation* loop is NOT a 4-day build — but it is also not needed; the robust-selection version captures the value.
- **Novelty verdict: GENUINELY UNTRIED in protein design** (no source found applying DRO to binder design). High upside, low implementation cost — strong candidate for a headline "novel method" in the dossier.

### Gaps
- No protein-design application of DRO found at all (could be a true gap in the literature, or simply not indexed) — I searched "distributionally robust optimization protein design" and found only generic mutational-robustness work. Report as "not found," which is itself the signal.
- Canonical Wasserstein-DRO theory paper (Mohajerin Esfahani & Kuhn, Math. Programming 2018) not fetched directly; cited via downstream data-driven-DRO papers.

---

## Q7. Category theory / sheaves / operads (multi-chain consistency + oracle fusion)

### Takeaway
**Sheaf theory** is the mathematically correct framework for *fusing disagreeing local observations into a global-consistency score* — exactly the multi-oracle (AF2/Boltz/Chai) and multi-chain (two-protomer/trimer) consistency problem. The **consistency radius / sheaf cohomology** give a computable "do my oracles/chains agree" diagnostic. Operads/ologs for hierarchical protein assembly exist but are more conceptual. Sheaf fusion is mature in sensor networks, essentially untried in protein design.

### Cited Findings
- "Sheaves are the canonical data structure for sensor integration"; a sheaf over a simplicial/cell complex assigns data to local pieces with restriction maps, and the **consistency radius** quantifies how far the local data are from gluing into one global section; data fusion = finding the globally-consistent section by constrained optimization — [Robinson, arXiv:1603.01446](https://arxiv.org/pdf/1603.01446).
- **Sheaf cohomology** detects problems in the fusion process and quantifies the impact of each data relationship on overall consistency; the workflow is formalize-system -> quantify local consistency -> derive global-consistent data -> measure relationships via cohomology — [sheaf sensor-fusion summary, arXiv:1603.01446](https://arxiv.org/pdf/1603.01446); applied to heterogeneous geolocation uncertainty — [PMC7349656](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7349656).
- Category theory (ologs) and **operads** model hierarchical protein materials: operads "control how a higher-level system is formed as an arrangement of interacting components," and protein-material ologs share category-theoretic structure with other complex systems — [Giesa, Spivak & Buehler, Category-Theoretic Analysis of Hierarchical Protein Materials, PLoS ONE 2011, PMC3169555](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3169555/).

### Inferences
- **Hook #1 (oracle fusion / trust filter):** build a tiny sheaf whose cells are {AF2-ig, Boltz-2, Chai-1} predictions for a design and whose restriction maps encode "these two oracles should agree on the pose/interface-contacts." The **consistency radius is a per-design trust score**: low radius = oracles agree (safe to submit), high radius = disagreement (demote, as in the memory lesson about folder pessimism). This is a principled upgrade over ad-hoc "require all three to pass."
- **Hook #2 (multi-chain / multi-protomer consistency):** TNF is a trimer and the epitope spans two protomers; make a sheaf over the chains/interfaces and require the binder's predicted contacts to glue consistently across the two-protomer groove (and across the 3 symmetric interfaces) — the consistency radius flags designs that only "work" on one protomer copy or break the trimer symmetry.
- **Hook #3 (conceptual, lower priority):** operad/olog formalization of "paratope module + switch module + cross-reactive-anchor module" composing into a binder — useful for *organizing* the design logic and the dossier narrative, not a 4-day scoring tool.
- **Feasibility (4 days): MEDIUM.** `pysheaf` exists but is finicky; the consistency-radius computation for a small, finite sheaf is implementable by hand (assign global values minimizing max local disagreement under restriction maps = a small optimization, ~30-50 lines). Operads/ologs are not a coding deliverable in 4 days (conceptual only).
- **Novelty verdict: sheaf fusion = GENUINELY UNTRIED in protein design** (mature only in sensor fusion); operads/ologs for protein *materials* = EMERGING/conceptual, untried for binder *design/selection*.

### Gaps
- No source applies sheaves to folding-oracle fusion or multi-chain structural consistency; the transfer is by analogy to sensor fusion. Strong conceptual fit, zero prior art found.

---

## Q8. Algebraic / geometric invariants for structural novelty beyond TM-score

### Takeaway
TM-score is a Cα-superposition similarity and misses topological/contact-arrangement differences. Three alignment-free invariants give *orthogonal* novelty signals: **persistent-homology (persistent Betti) fingerprints**, **circuit topology** (the arrangement of intra-chain contacts: parallel/series/cross, with polynomial invariants), and **Gromov-Wasserstein distortion** (Q2). PH/GW novelty is emerging-general; circuit topology as a de novo *binder-novelty* score is untried.

### Cited Findings
- Persistent-homology **persistent Betti functions (PBF)** reduce a complex biomolecular structure to simple 1-D functions for quantitative, resolution-robust structure comparison (a "persistent similarity") — [A quantitative structure comparison with persistent similarity, arXiv:1707.03572](https://arxiv.org/pdf/1707.03572); PH extracts phylogenetic/structure information robustly from noisy structures — [arXiv:1811.00252](https://arxiv.org/pdf/1811.00252).
- **Circuit topology** classifies the arrangement of intra-chain contacts within a folded linear chain (parallel/series/cross relations), enabling classification/comparison of chains; **polynomial invariants of molecular circuit topology** have been constructed (efficient, handling soft/hard contacts), and knot-theory tools (quandles) give further invariants — [Polynomial invariant of molecular circuit topology, arXiv:2109.02391](https://ar5iv.labs.arxiv.org/html/2109.02391); [Circuit topology for bottom-up engineering of molecular knots, arXiv:2106.03925](https://ar5iv.org/html/2106.03925).
- **Gromov-Wasserstein** reveals "topological rearrangements in homologous folds" and "recurrent short structural motifs," alignment-free — [GWProt, bioRxiv 2025.11.03.686256](https://www.biorxiv.org/content/10.1101/2025.11.03.686256.full.pdf).
- TM-score quantifies topological similarity from aligned Cα distances with length-independent normalization; >=0.5 = same fold — [TM-score man page / Zhang](https://manpages.org/tm-score).

### Inferences
- **Hook (novelty beyond TM-score, feeds the Adaptyv novelty requirement):** compute a **multi-invariant novelty vector** for each binder vs the PDB/AFDB nearest neighbors: (i) PBF/persistence-image distance, (ii) circuit-topology fingerprint distance (contact-arrangement histogram), (iii) GW distortion. A design can be TM-score-similar yet circuit-topology- or PH-distinct (novel contact wiring) — these catch novelty TM-score blurs, and can be combined with the existing ESMFold-oracle novelty calibration as orthogonal evidence. Circuit topology additionally characterizes *which contacts form the switch* (the contacts that must break at pH 6.0 are specific circuit relations).
- **Feasibility (4 days): LOW-MEDIUM.** PBF/persistence images = `gudhi`/`persim` (same stack as Q1). Circuit topology = contact map + a combinatorial relation classifier over contact pairs (parallel/series/cross is a simple ordering test), ~40-60 lines; polynomial invariant is more work but the raw fingerprint histogram is enough for a novelty distance. GW = `POT` (same as Q2).
- **Novelty verdict:** PH/GW structural novelty = EMERGING (general structure comparison, not binder-design novelty scoring). **Circuit-topology as a de novo binder novelty/switch descriptor = UNTRIED.**

### Gaps
- No source uses circuit topology or PBF specifically as a *de novo binder novelty filter*; evidence is for general structure comparison/classification. Whether these invariants correlate with the Adaptyv portal's ESMFold+TMalign novelty metric is unknown (would need calibration against portal labels, per the existing memory lesson).

---

## Cross-cutting notes for the report-writer

- **Highest value-to-risk (ship these first):** (1) GW distance for pH-two-state switch magnitude + cross-species epitope consistency [Q2]; (2) PH Wasserstein distance between pH-7.4/6.0 interface ensembles [Q1]; (3) robust max-min selection over a {human,mouse,+perturbation} target ambiguity set = sampled Wasserstein-DRO [Q6]; (4) Hill-slope/proximity-to-fold objective on the PROPKA/pKAI titration curve [Q3]; (5) sheaf consistency-radius oracle-trust filter [Q7]; (6) DPP for the diverse <=20 [Q2]. All are LOW/MEDIUM risk, all pip-installable (`gudhi`, `persim`, `POT`, `dppy`, `scipy`), all operate on outputs the pipeline already produces.
- **Boldest genuinely-untried transfers (dossier headline candidates):** DRO for cross-reactivity [Q6] and sheaf fusion for oracle/multi-chain consistency [Q7] — no prior protein-design art found for either.
- **Common implementation stack:** GUDHI/Ripser/persim (TDA), POT (Wasserstein/GW/barycenter), dppy or hand-rolled greedy MAP (DPP), scipy (Hill fit, MI, merge trees), existing OpenMM/xTB/PROPKA/pKAI/AF2-ig/Boltz/Chai outputs as inputs. No new heavy dependencies, all CPU-side post-processing of existing GPU outputs.
- **Priority alignment:** every TOP-priority (pH-switch) hook is in Q1/Q2/Q3/Q5; cross-reactivity hooks in Q2/Q6; selection/novelty in Q2/Q4/Q7/Q8. Affinity is deliberately served last (barely a dedicated hook), matching the stated objective ordering.
