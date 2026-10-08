# Conditional / pH-switchable / environmentally-responsive / allosteric protein binder design + validation methods (2024–2026, beyond Baker 2025)

Scope anchor: our #1 objective is a de-novo binder that BINDS human TNF-α at pH 7.4 and RELEASES at pH 6.0 (endosomal-recycling direction). Current mechanism = interface His adjacent to a conserved target Lys (His+ repels Lys+ at low pH → release), validated by two-state MD + PROPKA Wyman linkage, with cooperative 2-His variants. This catalog lists methods we may have MISSED, flagging for each: mechanism, direction-vs-SELECTIVITY, code/license/uv+https installability, relevance to a pH-6.0 TNF release switch + structural-carryover refinement, and DOI/date.

---

## Q1. 2025–2026 methods for pH-dependent / conditional / logic-gated / state-switching binder design beyond Baker 2025

### Takeaway
Beyond Baker's pH-binder paper there is a small but distinct 2024–2026 toolkit: (a) protease-activated "peptide mask" conditional miniprotein binders (fully de-novo, large on/off, orthogonal to pH); (b) de-novo two-state hinge proteins and allosteric switchable assemblies (ligand/peptide-triggered conformational switches); (c) LOCKR / Co-LOCKR / LucCage logic-and-biosensor switches; (d) emerging multistate generative models (ProDiT, two-state DL, RFdiffusion2). Critically, Baker 2025 itself scores only *direction* (electrostatic ddG, ad-hoc pH-score), not two-state SELECTIVITY, and installs His with a PyRosetta script — so our two-state-MD + Wyman-linkage pipeline is already more selectivity-aware than the reference method.

### Cited Findings
- **Baker 2025 (the reference) scores DIRECTION only, not two-state selectivity.** It uses RFdiffusion + ProteinMPNN (histidine-biased) + AF2 filtering (pae_interaction/pLDDT) + Rosetta FastRelax; pH-selectivity is scored by the **Rosetta "pH module"** (protonation-state sampling at pH 0 vs pH 7), a custom ad-hoc **"pH score"** (H-bond patterns weighted by solvent accessibility), and **"dddG_elec"** (Rosetta `fa_elec` scoreterm only). No PROPKA, no constant-pH MD, no Wyman/multistate linkage. Histidines are installed via **a PyRosetta script** (license blocker for us). Two mechanisms: (1) interface His as H-bond acceptor next to a cationic target Arg/Lys → repulsion on protonation; (2) buried His H-bond networks destabilizing the monomer fold at low pH. Targets include **TNFR2 and TNFα** plus EphA2, IL-6, PCSK9, Neo2. — [bioRxiv 10.1101/2025.09.29.678932, 2025-09-29](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740/)
- **Protease-activated de-novo peptide masks for conditionally-active miniprotein binders** (Escobar-Rosales, Baker, Oller-Salvia). Extends the binder C-terminus with a protease-cleavable linker + a de-novo masking helix that sterically blocks the interface; tumor-protease cleavage restores binding. ~half of 20 designs gave >100-fold affinity reduction (best >1000-fold EGFR); binding restored in 19/20 on cleavage; micromolar binder–mask affinity suffices for robust OFF + fast ON. Targets EGFR domI/III, FGFR2, IL7Rα. This is a SELECTIVITY (two-state OFF/ON) conditional mechanism, orthogonal to pH. — [bioRxiv 10.1101/2025.02.25.639922, 2025-02-25](https://pmc.ncbi.nlm.nih.gov/articles/PMC12703739/); published JACS 2025, 147(49) — [americanpeptidesociety.org summary](https://americanpeptidesociety.org/research/conditional-miniprotein-activation/)
- **De-novo stimulus-responsive two-state HINGE proteins** (Praetorius, Leung et al.). Proteins populate designed state 1 (apo) and designed state 2 (ligand-bound); explicit sculpting of an energy landscape with two minima. Conformational + binding equilibria are thermodynamically coupled; validated by X-ray, cryo-EM, DEER. Note: Science **2023** (not 2024). Follow-up crosslinks hinges to azobenzene for light-driven cycling. — [Science 2023, PMC10697137](https://pmc.ncbi.nlm.nih.gov/articles/PMC10697137/)
- **De-novo allosterically switchable protein assemblies** — peptide-switchable hinge modules rigid-body-coupled to interfaces to direct alternative oligomeric states (MWC-cooperativity-inspired); cyclic rings incorporate/eject subunits on peptide binding; up to 40-fold activation. — [bioRxiv 10.1101/2023.11.01.565167, 2023-11](https://www.biorxiv.org/content/10.1101/2023.11.01.565167.full.pdf)
- **LOCKR** (de-novo bioactive protein switches): 5-helix "cage" whose interface binds either an intramolecular "latch" or an intermolecular "key" peptide; functional motifs on the latch act only when the key displaces it. — [Langan et al., Nature 2019, PubMed 31341284](https://pubmed.ncbi.nlm.nih.gov/31341284/)
- **Co-LOCKR** (colocalization-dependent switches): AND/OR/NOT Boolean logic, activating by conformational change only when co-localized on a target cell surface; used to gate T-cell targeting on antigen combinations. — [Lajoie et al., Science 2020, aba6527](https://www.ipd.uw.edu/?p=5553)
- **LucCage/LucKey biosensors**: cage+latch / key + analyte-binder; target binding thermodynamically flips OFF→ON, reconstituting split luciferase (e.g., lucCagePTH, 10 nM LOD). A de-novo allosteric readout for binding. — [Quijano-Rubio et al., Nature 2021, PMC8085813](https://pmc.ncbi.nlm.nih.gov/articles/8085813/)
- **Schrödinger pH-dependent antibody workflow** (industry): physics-based FEP+ **combined with constant-pH MD**, His substitutions to tune pH sensitivity; experimentally confirmed strong binding at neutral / weak at acidic pH. Proprietary (Schrödinger license). — [AET 2025, Sunidhi Lenka](https://www.schrodinger.com/life-science/resources/event/antibody-engineering-therapeutics-2025/)
- **mAbs 2025**: rational design of pH-dependent antigen binding from structural insights into broadly-neutralizing anti-α-neurotoxin antibodies (histidine-centric). — [prelights/mAbs 2025](https://prelights.biologists.com/highlights/computational-design-of-ph-sensitive-binders/)

### Inferences
- We are already "beyond Baker" on the axis that matters for a *release* switch: Baker reports fold-change (direction) and uses `fa_elec`-only dddG; our two-state MD occupancy + Wyman linkage targets SELECTIVITY. Baker's cheap `dddG_elec` (electrostatic-only ddG between protonation states) is, however, a fast open-reimplementable *filter* we could add (OpenMM/pdb2pqr electrostatics) without PyRosetta.
- Protease-mask conditionality is the strongest fallback if the His-pH ΔΔG proves too weak at the TNF interface: it is fully de-novo, gives >100–1000× OFF/ON, and is pH-independent. But it changes the product profile (protease-gated, not endosomal-recycling) — note for the writer that it does not satisfy the pH-6.0 recycling spec, it is an alternative conditional modality.
- Allosteric (hinge/LOCKR) switches give true two-state SELECTIVITY but our own prior lesson (allosteric hinged-lid binders are ~1/512 buildable; diffusion+MPNN+folder give monolithic folds) makes them a low-yield route for a TNF mini-binder; keep as a validation-readout idea (LucCage) rather than the primary binder mechanism.

### Gaps
- No published **generative** method as of 2026-10 explicitly designs a pH switch at a binder interface (i.e., a diffusion/flow model conditioned on "bind at pH 7.4, release at pH 6.0"). Baker 2025 remains the only dedicated pH-binder design paper; everything else is either conformational-switch design or conditional via masking/colocalization.

---

## Q2. His/titratable-residue switch design; constant-pH / two-state scoring; pKa predictors beyond PROPKA

### Takeaway
There is a rich ML pKa-predictor landscape beyond PROPKA (pKAI/pKAI+, DeepKa, PypKa, pKALM, pKa-ANI, pKAANI) — and the key, directly-actionable finding is that **PROPKA3 is the WORST method for histidine pKa** in head-to-head benchmarks (His R²≈0.05), exactly the residue our switch depends on. Use a sign-consensus across method-classes (empirical ML + physics-based + PLM + QM-environment). For two-state scoring, open constant-pH-MD exists but is effectively NO-GO in a uv/https env (unreleased or Java), vindicating our explicit two-state MD; Rosetta multistate design is capable but PyRosetta-license-blocked.

### Cited Findings
- **Benchmark (JCIM 2023, 10.1021/acs.jcim.3c00165):** compared PROPKA3, DeepKa, pKAI, pKAI+ (empirical) and DelPhiPKa, MCCE2, H++ (physics). **For the 41 histidines: DeepKa RMSE 0.89 (best empirical, R²0.43); H++ 0.94 (R²0.26); pKAI+ 0.94 (R²0.27); PROPKA3 1.20 (R²0.05, "no clear advantage").** A simple **DeepKa+PROPKA3 arithmetic consensus** gave the best overall (RMSE 0.76, R²0.45). All seven were installed and run locally. — [PMC10466379](https://pmc.ncbi.nlm.nih.gov/articles/PMC10466379/)
- **pKAI / pKAI+** — first DL predictor of pKa shifts trained on continuum-electrostatics data; >1000× faster than physics methods; pKAI+ tuned to predict experimental pKa from a single conformation. Pip-installable (`pip install pkai`), open. — [pypi.org/project/pkai](https://pypi.org/project/pkai); [Reis et al., JCTC 2022, PMC9369009](https://pmc.ncbi.nlm.nih.gov/articles/PMC9369009)
- **PypKa** — open-source Python module for Poisson–Boltzmann pKa (proton tautomerism), CPU-parallel, runs on PDB or MD frames, also predicts pI and assigns protonation states; simple API/CLI (`pip install pypka`). — [pypka.readthedocs.io](https://pypka.readthedocs.io/)
- **DeepKa** — DL predictor trained/validated on **continuous constant-pH MD of 279 soluble proteins**; grid-charge representation of electrostatics; accuracy close to CpHMD; best single empirical His performer above. ACS Omega 2021; a DeepKa-Web server exists (standalone-code availability not confirmed in sources). — [ChemRxiv 61570a1b](https://chemrxiv.org/engage/chemrxiv/article-details/61570a1bcc53d94e30b901bd)
- **pKALM** — protein-language-model (sequence-only) pKa predictor; 8 ionizable groups (Asp/Glu/His/Lys/Cys/Tyr + termini); RMSE 0.8321; 4,961 pKa/sec; GitHub + license linked from the server. JCTC 2025. — [bioRxiv 10.1101/2024.09.16.613101](https://www.biorxiv.org/content/10.1101/2024.09.16.613101.full.pdf); [onodalab server](https://onodalab.ees.hokudai.ac.jp/pkalm)
- **pKa-ANI** (Gokcan & Isayev) — deep representation learning using ANI-2x atomic-environment vectors + learned QM representation; coordinates-only input; predicts all five titratable types. An environment/burial-aware (QM-flavored) class distinct from continuum-electrostatics ML. — [PMC8864681](https://pmc.ncbi.nlm.nih.gov/articles/PMC8864681/)
- **pKAANI** — additional pip-installable ML pKa tool. — [pypi.org/project/pkaani](https://pypi.org/project/pkaani/)
- **PROPKA3** — empirical, open (jensengroup), the method we currently use; benchmark-worst for His. — [github.com/jensengroup/propka](https://github.com/jensengroup/propka)
- **His-network / titratable-switch DESIGN (buried-His, two-state).** De-novo "tunable pH-driven conformational changes": buried His in preorganized H-bond networks, stable above pH 6.5, cooperative large-scale conformational change as network His protonate; transition pH + cooperativity tuned by number of His networks and surrounding hydrophobic strength. pH-responsive filaments from subunits with 6 or 9 buried His: sharp reversible assembled↔disassembled transition over **0.3 pH units**, disassembly <1 s. This is the design+validation precedent for Baker's "monomer-destabilization" (buried-His) mechanism and for cooperative multi-His switches. — [Janelia/IPD; Nat. Nanotechnol. 2024, 10.1038/s41565-024-01641-1](https://link.springer.com/article/10.1038/s41565-024-01641-1); PDB 6MSQ/6MSR
- **Two-state / multistate SCORING tools.** Rosetta `mpi_msd` (explicit negative states) and `ProteinInterfaceMS` (target state = bound complex; competitor states = unbound/unfolded partners) implement multistate design with explicit negative design — but require Rosetta/PyRosetta (license blocker). — [Rosetta ProteinInterfaceMSMover docs](https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/ProteinInterfaceMSMover); generic MSD program [Leaver-Fay, PLoS ONE 2011, PMC3130737](https://pmc.ncbi.nlm.nih.gov/articles/PMC3130737/)
- **Constant-pH MD (open, for two-state validation).** `protons` (Chodera lab, OpenMM) implements CpHMD for protonation/tautomer sampling but has **no official release** (dev-stage). A `molecular_simulations` module offers OpenMM constant-pH with pH replica exchange. Force Field X (FFX) implements polarizable (AMOEBA) CpHMD but is **Java**. Amber (pmemd CpHMD) and CHARMM (PHMD) are the mature implementations but are licensed C/Fortran stacks. — [protons.readthedocs.io](https://protons.readthedocs.io/); [molecular-simulations cph_simulation](https://molecular-simulations.readthedocs.io/en/stable/api/molecular_simulations.simulate.cph_simulation.html)

### Inferences
- **Immediate upgrade to our pKa stack:** our Wyman linkage currently rests on PROPKA, the *weakest* His predictor (R²≈0.05). Add DeepKa (best empirical His), pKAI+/H++ (next best), and an orthogonal class — pKALM (PLM, sequence-only) or pKa-ANI (QM-environment) — then read **sign-consensus** on ΔpKa(His) between bound and unbound/low-pH states. This operationalizes our existing lesson "pH-switch pKa landscape: method-classes for sign-consensus" with concrete, uv-installable tools.
- The benchmark's winning move (DeepKa+PROPKA3 consensus) suggests a cheap ensemble ΔΔG_pH: compute the Wyman integral with each predictor and keep sign-consensus; disagreement on His > ~0.5 pKa = flag (matches our buried-carboxylate undesignability lesson).
- Open CpHMD being NO-GO in uv/https (protons unreleased; FFX Java; Amber/CHARMM licensed) confirms our explicit fixed-protonation two-state MD + xTB is the right open path; there is no drop-in open CpHMD to switch to.
- Rosetta multistate (ProteinInterfaceMS) is conceptually exactly our "bound vs unbound-protonated" scoring but is license-blocked; our MD occupancy delta is the open equivalent and (per our lessons) catches selectivity that static linkage misses.

### Gaps
- No source gives a His-specific pKa benchmark that includes pKALM or pKa-ANI head-to-head with DeepKa on *interface* histidines (benchmarks use mostly buried/soluble-protein His) — interface/surface His accuracy for our exact geometry is unquantified.
- DeepKa standalone (non-web-server) code location/license not confirmed in the sources (web server exists; a GitHub repo was not surfaced).

---

## Q3. Allosteric / conformational-switch binder design (hinge, lid, order-disorder)

### Takeaway
The open de-novo allosteric toolkit is Baker/IPD-centric: two-state hinges, peptide-switchable allosteric assemblies, LOCKR/Co-LOCKR latch-cage-key, and LucCage biosensors — all give genuine two-state SELECTIVITY but are hard to realize as a compact TNF mini-binder (monolithic-fold failure mode we already hit). For VALIDATING a conformational switch cheaply and open-source, AF-Cluster (MSA-clustering + AF2) predicts alternative/fold-switch states and is sensitive to single point mutations.

### Cited Findings
- **Two-state hinges** (apo vs ligand-bound, thermodynamically coupled) — see Q1; the canonical de-novo "lid/hinge" conformational switch, atom-accurate in both states. — [Science 2023, PMC10697137](https://pmc.ncbi.nlm.nih.gov/articles/PMC10697137/)
- **Allosterically switchable assemblies** — peptide-triggered oligomeric-state switching, 40× activation, MWC design logic. — [bioRxiv 2023.11.01.565167](https://www.biorxiv.org/content/10.1101/2023.11.01.565167.full.pdf)
- **LOCKR / Co-LOCKR / LucCage** — latch–cage–key allostery; Boolean logic; split-reporter biosensing (Q1 citations).
- **Two-state protein design via deep learning** (UW thesis line): de-novo backbone design of *two* states + sequence design using per-position amino-acid probabilities + simulated annealing; two designed sequences expressed and showed **temperature-dependent** switching between states — a generative route to order/order two-state proteins. — [UW ResearchWorks](https://digital.lib.washington.edu/researchworks/items/bc4faf1c-f575-4180-93cb-bdc5233d941a/full)
- **AF-Cluster (validation).** Clustering an MSA by sequence similarity lets AF2 sample alternative states of metamorphic proteins (KaiB, RfaH, Mad2) with high confidence; **sensitive to point mutations** — the authors designed+verified a 3-mutation set that flips KaiB's state. Open; runs on AF2/ColabFold (uv-installable via localcolabfold). A cheap first-pass check that a designed sequence can populate two states. — [Wayment-Steele et al., Nature 2023; bioRxiv 2022.10.17.512570, PMC10808063](https://pmc.ncbi.nlm.nih.gov/articles/PMC10808063/)

### Inferences
- For TNF, an allosteric/order-disorder binder is a poor primary bet (our own buildability lesson: monolithic folds, ~1/512), but AF-Cluster is a near-free addition to our validation stack: run it on the His-switch binder sequence to check whether the low-pH (protonated) state has any propensity to adopt a distinct backbone state — useful negative control that our switch is electrostatic (same fold, weaker interface) rather than a fold-switch artifact.
- LucCage is the natural way to turn our TNF binder into a pH biosensor for wet-lab readout of the switch, if an orthogonal confirmation of two-state behavior is wanted.

### Gaps
- None of the allosteric design methods ship a pH-triggered variant; they are peptide/ligand/temperature/light-triggered. Converting them to a proton trigger (buried-His network inside the hinge) is implied by the Q2 buried-His work but not demonstrated as a *binder*.

---

## Q4. Methods that combine a switch with affinity maturation while preserving the pose (switch + structural carryover — "ArcRefine-like")

### Takeaway
Two pose-preserving maturation engines are current and relevant — PPIFlow in-silico maturation (flow-matching) and RFdiffusion partial diffusion — but BOTH optimize affinity only; neither imposes or preserves a pH switch. The practical recipe is maturation-with-constraints: freeze the switch residues (His + partner) and the interface pose, mature the periphery, then re-score the switch two-state. No published tool does switch+maturation jointly.

### Cited Findings
- **PPIFlow** — flow-matching generative model + a dedicated **in-silico affinity-maturation stage: interface rotamer enrichment + partial flow refinement to optimize interfacial packing**; effective for de-novo VHHs at pM/nM affinity. Licensed CC-BY-NC 4.0; code release not stated in the abstract. Explicitly **affinity-only — no conditional / pH / two-state**. — [bioRxiv 10.64898/2026.01.19.700484, 2026-01-22, Mingchen Chen / Changping Lab](https://doi.org/10.64898/2026.01.19.700484)
- **Partial diffusion (RFdiffusion)** — treat a suboptimal/validated binder as the starting structure, add controlled noise and denoise ("back-and-forth") to explore nearby structure space **while preserving the binding mode**; used to generate picomolar binders by refining existing designs. The canonical pose-carryover operator; open-source (RFdiffusion, RosettaCommons). — [Watson et al., Nature 2023 (RFdiffusion)](https://ideas.repec.org/a/nat/nature/v626y2024i7998d10.1038_s41586-023-06953-1.html); [partial-diffusion usage, tamarind docs](https://docs.tamarind.bio/tools/rfdiffusion); [protflow binder tutorial](https://protflow.readthedocs.io/en/latest/_sources/tutorials/binder_design.rst.txt)
- **ProteinMPNN fixed-position redesign** — standard open route to re-sequence everything except pinned switch/interface residues while holding the backbone; the open, PyRosetta-free analogue of Rosetta FastRelax carryover. (Implicit across the binder-design toolchains above; our own lessons already use fixed-position MPNN.)
- **ODesign** — "world model for biomolecular interaction design" (generative framework spanning interaction design), 2025; potential maturation/refinement engine, conditional-switch support unverified. — [arXiv 2510.22304](https://arxiv.org/pdf/2510.22304)

### Inferences
- "ArcRefine-like" = exactly PPIFlow-style maturation (freeze critical interface residues, noise the neighbors, partial-refine) OR RFdiffusion partial diffusion at low noise, with the SWITCH residues added to the frozen set. The missing capability is a maturation objective that *also* conserves ΔΔG_pH; today that must be bolted on as a post-filter (re-run two-state MD + multi-pKa consensus on every matured design), not optimized inside the generator.
- Because PPIFlow and partial diffusion both re-pack/perturb the interface, they RISK erasing the His···Lys geometry — our prior lesson "AF2 binding is sequence-specific; re-seq of a binding backbone can kill binding" applies doubly here: mature by perturbing periphery only, keep His/partner + their rotamers fixed, and re-validate the switch, not just affinity.
- A defensible open carryover stack with no license blockers: RFdiffusion partial diffusion (pose) → ProteinMPNN fixed-position (sequence, His pinned) → AF2/af2ig (fold+binding) → two-state MD occupancy + multi-pKa Wyman consensus (switch). This mirrors Baker's generate→MPNN→AF2 loop but replaces their PyRosetta His-install + `dddG_elec` with our open switch scorer.

### Gaps
- No method optimizes affinity and a pH/conditional switch jointly while preserving pose; "switch + carryover" is an unfilled niche as of 2026-10. PPIFlow code availability/license-for-use is unconfirmed (CC-BY-NC = non-commercial).

---

## Q5. Which methods are open/runnable with no PyRosetta blockers (uv/pip + https)

### Takeaway
The switch-scoring and pose-carryover pieces we need are available open-source and pip-installable; the PyRosetta-dependent pieces (Baker's His-install script, Rosetta multistate `ProteinInterfaceMS`, Rosetta pH-module dddG) are the only hard blockers, and each has an open substitute.

### Cited Findings
- **Open / pip (uv + https) — pKa:** pKAI/pKAI+ (`pip install pkai`), PypKa (`pip install pypka`), pKAANI (`pip install pkaani`), PROPKA3 (github jensengroup, pip), pKALM (GitHub + license, PLM), pKa-ANI (repo, ANI-2x). — citations in Q2.
- **Open — design/carryover:** RFdiffusion + partial diffusion, ProteinMPNN, AF2/ColabFold, AF-Cluster — all open-source, no PyRosetta. — citations in Q1/Q3/Q4.
- **Blocked / licensed:** Baker 2025 His-install **PyRosetta** script; Rosetta `mpi_msd` / `ProteinInterfaceMS` multistate (Rosetta license); Rosetta pH-module `dddG_elec`/pH-score; Schrödinger FEP++CpHMD workflow (commercial); Amber/CHARMM CpHMD (licensed); FFX CpHMD (open but **Java**, not uv/pip); `protons` OpenMM CpHMD (open but **no release**). — citations in Q1/Q2.

### Inferences
- Every license-blocked capability maps to an open substitute we can run: PyRosetta His-install → ProteinMPNN fixed-position + manual His placement; Rosetta multistate → our explicit two-state MD occupancy; Rosetta `dddG_elec` → pdb2pqr/APBS or OpenMM electrostatic ΔΔG between protonation states (cheap open reimplementation of Baker's filter); PROPKA-only Wyman → multi-class pKa consensus (DeepKa/pKAI+/pKALM/pKa-ANI).
- The one genuinely unavailable open capability is mature constant-pH MD in uv/https — reaffirming explicit fixed-protonation two-state MD as the open workhorse.

### Gaps
- Could not confirm from sources whether PPIFlow or DeepKa expose a usable open code release (PPIFlow CC-BY-NC, no repo cited; DeepKa web server only). Treat both as "verify install first-hand" (our standing lesson) before relying on them.
