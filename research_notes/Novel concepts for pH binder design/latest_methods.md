# Latest (2024–2026) computational methods for de-novo protein binder design, with a focus on hard/flat polar cytokine targets (TNF-α) and pH-conditional binders

Scope note: compiled 2026-10-08. Every method is dated. Project target = human TNF-α (flat, polar, β-jellyroll homotrimer, PDB 1TNF); goals = pH-conditional binding (bind 7.4 / release 6.0), mouse cross-reactivity, affinity. Existing stack = RFdiffusion1/3, ProteinMPNN, AF2-initial-guess, Boltz-2, BoltzGen, partial diffusion.

---

## Key Question 1 — State-of-the-art generators (2025–2026): capabilities, cytokine/flat-PPI success, experimental hit rates, availability

### Takeaway
The frontier has split into (a) AF2-hallucination pipelines (BindCraft) and (b) all-atom diffusion/flow co-design generators (RFdiffusion2/3, BoltzGen, Chai-2, La-Proteina). All-atom generators that model side chains/ligands directly are the newest and most relevant for flat polar surfaces; BoltzGen and RFdiffusion3 are open (MIT / Foundry) and already in or adjacent to the project stack, while AlphaProteo and Chai-2 — the two with the most striking hit rates — are proprietary. Critically, TNF-α is a documented failure case for one SOTA generator (AlphaProteo) yet a documented success for others (classical Rosetta minibinders; BoltzGen benchmark), so generator choice and epitope choice both matter.

### Cited Findings

**BindCraft (EPFL/MIT, Pacesa, Ovchinnikov, Correia; Nature 2025; preprint/record)**
- Open-source automated pipeline; backpropagates through AlphaFold2-multimer (MSA-free/single-sequence) to hallucinate binder sequences, guided by multiple loss functions, then ProteinMPNN + AF2 filtering — [PMC11642741](https://pmc.ncbi.nlm.nih.gov/articles/PMC11642741); [GitHub martinpacesa/BindCraft](https://github.com/martinpacesa/BindCraft)
- Reported experimental success rates of **10–100%**, nanomolar affinity, "one-shot": hits in nM range from ≤10 designs; validated against cell-surface receptors, allergens, CRISPR-Cas9 — [PMC11642741](https://pmc.ncbi.nlm.nih.gov/articles/PMC11642741); [Rosetta Commons summary](https://rosettacommons.org/2025/09/25/one-shot-design-of-functional-protein-binders-with-bindcraft/)

**AlphaProteo (Google DeepMind; Zambaldi et al., Sept 2024)**
- 9–88% experimental success across 7 targets; 3–300× better affinity than prior methods; one round of medium-throughput screening; targets = BHRF1 (88%), SARS-CoV-2 RBD, IL-7Rα, PD-L1, TrkA, IL-17A, VEGF-A — [arXiv 2409.08022](https://arxiv.org/abs/2409.08022); [DeepMind blog](https://deepmind.google/discover/blog/alphaproteo-generates-novel-proteins-for-biology-and-health-research/)
- **Explicitly FAILED on TNF-α**: "despite successes across seven other targets, the authors were unable to obtain binding hits for TNFα," attributed to a "flat, highly polar binding site" at the interface between two subunits of the homotrimer — [arXiv 2409.08022 §2.1](https://arxiv.org/abs/2409.08022); [emergentmind summary](https://www.emergentmind.com/open-problems/de-novo-binders-to-tnfa)
- Availability: **proprietary** — DeepMind has not released code or weights (no public repo/API found in search).

**RFdiffusion2 (Baker lab/IPD; preprint Apr 2025, Nature Methods Dec 2025)**
- Flow-matching model (trained from random init) for *unindexed atomic motif scaffolding* — scaffolds minimally-defined catalytic sites ("theozymes") without preset residue indices/rotamers; models joint distribution of rotamers, sequence indices, scaffolds from atom-level site descriptions — [biorxiv 2025.04.09.648075](https://www.biorxiv.org/content/10.1101/2025.04.09.648075); [IPD](https://www.ipd.uw.edu/2025/04/introducing-rfdiffusion2); [Nature Methods s41592-025-02975-x](https://www.ipd.uw.edu/wp-content/uploads/2025/12/s41592-025-02975-x.pdf)
- Primary demo is enzyme active sites, not binders, but the flow-matching + atomic-motif machinery underlies RFdiffusion3.

**RFdiffusion3 (Baker lab/IPD; preprint Sept 2025, released Dec 3 2025)**
- **Atom-level diffusion of backbone AND side-chain atoms** (atoms, not residues, are the diffused unit); explicitly models side chains, ligands, nucleic acids, metal ions during generation; ~**10× faster than RFdiffusion2** — [biorxiv 2025.09.18.676967](https://www.biorxiv.org/content/10.1101/2025.09.18.676967); [IPD "RFdiffusion3 now available"](https://www.ipd.uw.edu/2025/12/rfdiffusion3-now-available/)
- Availability: training + inference code and weights via **Rosetta Commons Foundry** (Dec 2025) — [Rosetta Commons](https://rosettacommons.org/2025/12/22/rfdiffusion3-is-now-available-in-foundry/). (Project memory already records RFd3 all-atom breaking the TNF-α wall → af2ig PAE 13.0 real binder.)

**BoltzGen (MIT, Stärk et al.; preprint Nov 20 2025) — already in project stack**
- All-atom diffusion model unifying structure prediction + design; built on Boltz-2 architecture with a **geometry-based residue encoding** (virtual atoms on backbone positions encode amino-acid identity, keeping everything in continuous space); SOTA folding too — [PMC12697729](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729); [biorxiv 2025.11.20.689494](https://www.biorxiv.org/content/10.1101/2025.11.20.689494)
- Modalities: nanobodies, globular proteins (80–180 aa), linear peptides (5–50), disulfide-bonded cyclic peptides, helicons (stapled), AMPs — [PMC12697729](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729)
- Experimentally validated across **26 targets / 8–9 campaigns**. For **5 benchmark targets that include TNF-α** (PD-L1, TNFα, PDGFR, IL-7Rα, InsulinR): nanobodies 80% target-level success (4/5), proteins 80% (4/5). For 9 novel targets (no >30% PDB identity): 66% success, best affinity ~6–9 nM — [PMC12697729](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729)
- Design-specification language: covalent bonds (disulfide/ligation/staple), structure conditioning (fixed scaffold, motif, pairwise distance constraints), binding sites (positive AND "not-binding"/negative residue annotations), secondary-structure (helix/sheet/coil) — [PMC12697729](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729)
- Availability: **MIT License**, weights + training/inference code + all designs — [GitHub HannesStark/boltzgen](https://github.com/HannesStark/boltzgen); [MIT jclinic](https://jclinic.mit.edu/boltzgen/)

**Chai-2 (Chai Discovery; technical report July 2025; follow-up Nov 2025)**
- Fully de-novo **antibody** design: **~16% hit rate** (close to 20%), >100× over prior computational methods (which were <0.1%); ≤20 designs per target across 52 targets; ≥1 hit for **50% of targets** in a single round — [biorxiv 2025.07.05.663018](https://www.biorxiv.org/content/10.1101/2025.07.05.663018); [Chai intro](https://www.chaidiscovery.com/news/introducing-chai-2)
- **Miniprotein** design: **68% wet-lab success**, routinely picomolar; 5/5 miniprotein targets validated — [biorxiv 2025.07.05.663018](https://www.biorxiv.org/content/10.1101/2025.07.05.663018)
- Hard target classes: ≥1 binder per target for **six GPCRs** (10–73 designs each); specific pMHC binders to KRAS G12V with no cross-binding to variants; follow-up "Drug-like antibody design against challenging targets" — [biorxiv 2025.11.29.691346](https://www.biorxiv.org/content/10.1101/2025.11.29.691346)
- Availability: **proprietary/closed** — partner-only (e.g., Eli Lilly, Jan 2026), early-access signup; no public API/weights — [search summary of silicon.co.uk / businesswire](https://www.chaidiscovery.com/news/introducing-chai-2)

**Boltz-1 / Boltz-2 (MIT) — structure+affinity predictors, used as generators' oracle** (see Q3).

**La-Proteina (NVIDIA; preprint Jul 2025, ICLR 2026)**
- All-atom co-design via *partially latent* flow matching (explicit coarse backbone + per-residue latent variables for sequence + atoms); SOTA on all-atom co-designability/diversity/validity; scales to **800 residues** where baselines collapse — [arXiv 2507.09466](https://arxiv.org/abs/2507.09466); [NVIDIA](https://research.nvidia.com/index.php/publication/2026-01_la-proteina-atomistic-protein-generation-partially-latent-flow-matching)
- Baselines it compares against (the all-atom co-design field): P(all-atom)/Pallatom, PLAID, ProteinGenerator, APM, Protpardelle — [arXiv 2507.09466](https://arxiv.org/abs/2507.09466). Primarily a monomer generator; binder conditioning not the headline.

**ProteinGenerator (Baker lab; Nat Biotech 2024)**
- RoseTTAFold *sequence-space* diffusion: jointly denoises sequence+structure, can bias amino-acid composition and guide with experimental sequence-activity data; demonstrated **multistate design** (one sequence folding to different supersecondary structures) and peptide caging (melittin) — [biorxiv 2023.05.08.539766](https://www.biorxiv.org/content/10.1101/2023.05.08.539766); open source (RosettaCommons/protein_generator). Multistate capability is conceptually relevant to conditional/switchable binders.

**EvoBind2 (Bryant/SciLifeLab; preprint 2024)**
- In-silico directed evolution with AF2/AF2-multimer for **linear + cyclic peptide** binders from target *sequence only* (no binding-site or scaffold needed; blind) — [PMC12280060](https://pmc.ncbi.nlm.nih.gov/articles/PMC12280060); [biorxiv 2024.06.20.599739](https://www.biorxiv.org/content/10.1101/2024.06.20.599739)
- Success rates **46% (linear), 75% (cyclic)**, affinities 0.26 nM – 7.9 µM — [PMC12280060](https://pmc.ncbi.nlm.nih.gov/articles/PMC12280060). Available on neurosnap/tamarind; open code.

**Germinal (Stanford/Arc Institute; preprint Sept 2025)**
- De-novo *epitope-targeted* antibodies/nanobodies: structure predictor + antibody-specific PLM co-optimize CDRs onto a user framework; **4–22% experimental success** across 4 targets testing only 43–101 designs; nM affinity, mammalian expression — [biorxiv 2025.09.19.677421](https://www.biorxiv.org/content/10.1101/2025.09.19.677421)

**TNFR-superfamily target-conditioned diffusion (preprint Sept 2024)**
- Target-conditioned diffusion generated potent TNFR-superfamily agonists AND antagonists — directly in the TNF/TNFR structural family — [biorxiv 2024.09.13.612773](https://www.biorxiv.org/content/10.1101/2024.09.13.612773)

### Inferences
- For a flat polar cytokine, the **all-atom co-design generators (RFdiffusion3, BoltzGen)** are the right class because they place explicit side chains/rotamers at the interface rather than relying on a coarse backbone + downstream MPNN — this matches the project's own observation that RFd3 broke the wall where RFd1/BoltzGen-on-wrong-epitope floored.
- The TNF-α contradiction (AlphaProteo fail vs BoltzGen/Rosetta success) is best explained by **epitope choice**, not raw generator power (see Q4) — the "unbindable" label attaches to the polar subunit-interface patch, not to all of TNF-α.
- BoltzGen's **negative/"not-binding" residue annotations** and secondary-structure conditioning are unusual levers worth exploiting to steer away from the polar patch and toward an edge β-strand.

### Gaps
- No dedicated "thermodynamic/physics-guided diffusion" binder method surfaced with experimental validation; the physics signal in the frontier is mostly post-hoc (Boltz-2 affinity, FEP — Q3) rather than guiding generation. Flow-matching for binders is represented by RFd2/3 and La-Proteina (monomer-focused).
- Could not confirm Chai-1 (the open predecessor) license specifics from search; believed open (chaidiscovery/chai-lab, code + weights) but unverified here.
- BoltzGen's TNF-α result is reported at the 5-benchmark-target level (80% = 4/5); the source summary did not isolate whether TNF-α specifically was among the 4 hits — verify in the paper's per-target table.

---

## Key Question 2 — Best practices for HARD / flat polar targets

### Takeaway
Four levers dominate for hard targets: (1) pick a *compact, contactable, ideally hydrophobic* hotspot/epitope rather than a dispersed polar patch; (2) reduce the target to a *discontinuous epitope-only template* (big in-silico win on hard targets); (3) *partial diffusion* of validated binders rather than pure de-novo; (4) bias sequence design to *satisfy interface polar atoms* and engage *β-strand/edge* conformations. Oracle filtering should be *consensus* (AF2-IG ipAE/ipTM + newer predictors) because single-metric pae_interaction is necessary but not sufficient.

### Cited Findings

**Hotspot / epitope selection**
- For TNF-α specifically, the successful Rosetta campaign chose **five hydrophobic hotspot residues in "Region I"** (the receptor-binding cleft): Val74, Leu75, Ile97 on one monomer + Trp114, Tyr115 on the adjacent monomer — selected for "hydrophobic feature and comparatively flat surface morphology," conserved across human/canine — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665). This is the opposite choice from AlphaProteo's failed polar subunit-interface patch.
- General principle (project memory, consistent with literature): hotspots must be a **compact contactable cluster (<~12 Å)**, verified by raw contacts, not scattered/elongated.

**Discontinuous epitope-only templates ("EpitopeCraft"; NeurIPS 2025)**
- Retain only the *discontinuous surface residues surrounding the binding site* (selected by distance to hot spots that contact known binders); AF2's `residue_index` feature re-positions the fragments. Rationale: protein-folding NNs are **"local-first"** — sensitive to local interactions, insensitive to global foldability — [arXiv 2509.25479](https://arxiv.org/abs/2509.25479)
- **Up to 80% in-silico success improvement for difficult targets** (TBLR1); ClpP (390 aa, dimeric) and ALS3 (326 aa) moved from "undesignable" to designable; **up to 44.2× faster** per successful design; built on/modified from BindCraft — [arXiv 2509.25479](https://arxiv.org/abs/2509.25479)
- Companion optimization steps: **Monte-Carlo evolution** (simulated annealing, accepts unfavorable mutations; targets low-pLDDT residues) gave **+130% to +364%** success for difficult targets; **biased MPNN redesign** with a *hydrophilic-occupancy penalty* (penalizes unpaired polar residues on the PPI interface) and a *surface-charge penalty* (keeps pI 7.5–9; raised eligible designs 10%→20%) — [arXiv 2509.25479](https://arxiv.org/abs/2509.25479)
- Related: "Programming co-folding to design binders for intrinsically disordered epitopes" (ICLR 2025) and "Harnessing protein-folding algorithms to drug intrinsically disordered epitopes" ([biorxiv 2025.11.11.687846](https://www.biorxiv.org/content/10.1101/2025.11.11.687846)).

**Extended / larger interfaces for flat polar surfaces**
- RFdiffusion was adapted to "generate backbones with a more extended contact surface that can engage dispersed surface hydrophobic residues" for flat polar targets (the TNFR-class problem) — [emergentmind summary](https://www.emergentmind.com/open-problems/de-novo-binders-to-tnfa)
- Minibinders to the TNFR1 ectodomain (another extended flat polar surface) reached micromolar to sub-micromolar — [search summary, PMC de novo therapeutics thesis](https://digital.lib.washington.edu/researchworks/items/12ac4931-49de-4d29-8268-b47112065358)

**β-strand / edge-strand binders (directly relevant to a β-jellyroll)**
- RFdiffusion can diffuse binders into **β-strand conformations** to engage targets: for the IDR IL-2RG, binders were designed to β-strand conformations of the target with **Kd 10–100 nM**; broader IDP set (amylin, C-peptide, VP48, BRCA1) Kd 3–100 nM — [search summary of Nature s41586-025-09248-9](https://www.ipd.uw.edu/wp-content/uploads/2025/08/s41586-025-09248-9-compressed.pdf). This is the de-novo capability closest to β-sheet augmentation against an edge strand of the TNF β-jellyroll.

**Partial diffusion / noising validated binders**
- Baker pH-binder pipeline uses partial diffusion as a refinement step after RFdiffusion (see Q4) — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740).
- Project memory (own data): partial diffusion of validated binders ~11% pass vs ~0.3% de-novo; re-fold ORIGINAL generator sequences rather than MPNN-reseq to amplify (binding is sequence-specific).

**Interface metrics beyond pae_interaction / consensus filtering**
- Bennett et al. (Nat Commun 2023) established **AF2-initial-guess pae_interaction** as the key discriminator of binders vs non-binders; outputs pLDDT, PAE, pae_interaction, RMSD — [Nat Commun s41467-023-38328-5](https://ideas.repec.org/a/nat/natcom/v14y2023i1d10.1038_s41467-023-38328-5.html); [GitHub nrbennet/dl_binder_design](https://github.com/nrbennet/dl_binder_design)
- Meta-analysis of **3,766 experimentally characterised binders** (preprint Aug 2025): AF2 confidence metrics (pLDDT, **ipAE, ipTM**) are predictive of in-vitro binding and **outperform physics-based Rosetta energies**; newer predictors (AF3, Boltz-1/2, Chai-1) have emerged but their in-vitro-binding predictivity is **not yet systematically validated**; generalizability of alternative metrics across diverse datasets is "unclear" — [biorxiv 2025.08.14.670059](https://www.biorxiv.org/content/10.1101/2025.08.14.670059)
- Practical metric set in use for modern pipelines: ipTM, complex_pLDDT, **complex_pde** (mean inter-chain PAE), Rosetta **shape complementarity** (computed on the Boltz-2 target-MSA complex), interface H-bonds/salt-bridges, ΔSASA, **ddG**, **SAP** score — [search summary; BoltzGen metrics guide](https://neurosnap.ai/blog/post/interpreting-boltzgen-metrics-and-filtering-in-protein-design/691cccde8b9522d6ffeff2d9); open scoring repo [DigBioLab/de_novo_binder_scoring](https://github.com/DigBioLab/de_novo_binder_scoring)

### Inferences
- The strongest actionable recipe for TNF-α flat-polar: **(a)** target the Region-I hydrophobic cluster (Trp114/Tyr115/Val74/Leu75/Ile97) or an edge β-strand, **not** the polar subunit interface; **(b)** reduce to a discontinuous epitope-only template; **(c)** use an all-atom generator (RFd3/BoltzGen) or partial-diffuse validated binders; **(d)** add a hydrophilic-occupancy penalty in MPNN so interface polar atoms get satisfied; **(e)** filter by consensus, not pae_interaction alone.
- "Flat + polar" is partly mitigable by **larger/extended interfaces** and **β-strand engagement**, both now demonstrated de-novo.

### Gaps
- The 3,766-binder meta-analysis PDF exceeded fetch size; I captured its qualitative conclusions but not the concrete AUC/enrichment numbers or recommended numeric thresholds per metric — worth a direct read for threshold calibration.
- No head-to-head benchmark found that ranks consensus oracles (AF2-IG vs Boltz-2 vs Chai-1 vs AF3) *specifically on de-novo mini-protein binders to cytokines*.

---

## Key Question 3 — Binding/affinity prediction & ranking for de-novo binders: reliability

### Takeaway
For small-molecule affinity, Boltz-2's affinity head is a genuine advance (approaches FEP at ~1000× lower cost). But for **protein–protein** affinity — the project's actual regime — the Boltz-2 affinity head does NOT transfer well and is beaten by sequence-based models; de-novo protein-binder ranking should still rest on **structural confidence (ipAE/ipTM/pae_interaction)** plus physics rescoring, not on a learned protein-protein affinity number.

### Cited Findings
- **Boltz-2 affinity (small molecule)**: first AI model to approach **FEP** accuracy on protein–ligand affinity; FEP+ 4-target benchmark **Pearson R ≈ 0.66** (beats all ML baselines, nears OpenFE physics); **>1000× faster** than FEP; wins CASP16 affinity track; full complex + energetics in <20 s — [biorxiv 2025.06.14.659707](https://www.biorxiv.org/content/10.1101/2025.06.14.659707); [PMC12262699](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699); [boltz.bio/boltz2](https://boltz.bio/boltz2)
- Availability: Boltz-1 and Boltz-2 are **open, MIT License** — [boltz.bio](https://boltz.bio/boltz2) (GitHub jwohlwend/boltz).
- **Boltz-2 for protein–protein affinity is weak** (Dec 2025): fine-tuned "Boltz-2-PPI" on TCR3d and PPB-affinity *underperforms sequence-based alternatives* in both small and large data regimes despite high structural accuracy; "current structure-based representations are not primed for performant [protein-protein] affinity prediction"; hybridizing Boltz-2-PPI embeddings with sequence embeddings gives complementary gains — [arXiv 2512.06592 (King et al., Synteny)](https://arxiv.org/abs/2512.06592)
- **AF2-IG pae_interaction** remains the validated primary ranking signal for de-novo binders (Q2); meta-analysis confirms AF2 ipAE/ipTM > Rosetta energies, newer predictors not yet systematically validated — [Nat Commun 2023](https://ideas.repec.org/a/nat/natcom/v14y2023i1d10.1038_s41467-023-38328-5.html); [biorxiv 2025.08.14.670059](https://www.biorxiv.org/content/10.1101/2025.08.14.670059)
- **AlphaFold3** as an oracle: structure code is open (non-commercial) and weights are available to academics on request (non-commercial only) — [DeepMind/AF3 open-source for academia](https://www.maginative.com/article/google-deepmind-makes-alphafold-3-open-source-for-academia/)
- **FEP / MD for binders**: the canine/human TNF-α minibinder work used MM/GBSA and MM/PBSA (GROMACS/AMBER22) binding-energy calculations as an orthogonal rescore of designed interfaces — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665).

### Inferences
- Rank de-novo TNF-α binders by **AF2-IG pae_interaction + ipTM/complex_pde** as the primary gate, with **Rosetta ddG/SC/SAP** and **MM-GBSA** as orthogonal rescoring; treat any learned *protein-protein affinity* score (incl. Boltz-2's) as low-trust until validated in-house.
- For the *pH-switch magnitude* specifically, a learned affinity model is unlikely to be reliable; explicit-state physics (the project's MD occupancy / QM switch readouts, PROPKA linkage) remains the honest metric.

### Gaps
- No reliable, experimentally-validated AI predictor of *absolute* protein–protein Kd for de-novo binders was found — the field consensus is that structural confidence + physics rescoring is still the state of practice.
- Reliability of AF3/Chai-1 ipTM vs AF2-IG pae_interaction *for de-novo binder discrimination* is unquantified in the sources found.

---

## Key Question 4 — 2025–2026 papers specifically on TNF-α / cytokine binders & antibody-alternative scaffolds

### Takeaway
TNF-α is simultaneously a famous failure (AlphaProteo, polar interface) and a 2025 success (sub-nM de-novo minibinders to a hydrophobic cleft; BoltzGen benchmark; Baker pH-switch TNF-α). The decisive variable is which surface you target. The Baker pH-sensitive-binders paper is the single most relevant prior art and already includes a TNF-α and a TNFR2 pH-switch.

### Cited Findings

**De-novo minibinders to human + canine TNF-α (Communications Biology 2025)**
- Classical Rosetta pipeline (RIF generation + PatchDock docking + Rosetta interface design + AF3 complex prediction + GROMACS/AMBER MD + MM-GBSA/PBSA); scaffolds = pre-validated inert de-novo **triple-helix bundles, 56–58 aa** (no RFdiffusion/ProteinMPNN) — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665)
- Epitope = **"Region I"** hydrophobic cluster at the receptor-binding cleft between monomers (Val74/Leu75/Ile97 + Trp114/Tyr115), conserved human/canine — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665)
- Scale: 169,975 docked configs from 1,000 scaffolds → top 5 → 4/5 expressed (~80%); **TMB-1 Kd 0.55 nM (canine), ~1.5 nM (human)**; TMB-2–4 1.0–2.5 nM; robust thermostability, stable in simulated GI fluids; cross-species L929 IC50 ~2–6 nM both species — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665)
- Structures: **PDB 9JEC**, EMDB EMD-61412 (cTNFα–TNF30K, 3.1 Å). DOI 10.1038/s42003-025-09030-7.
- Note: this is the explicit counter-example to AlphaProteo's TNF-α failure; same molecule, different (hydrophobic) epitope, different (Rosetta) method.

**Computational design of pH-sensitive binders (Baker lab; preprint Sept 29 2025 — THE project blueprint)**
- Two principles: **(1) interface electrostatic repulsion** — place His adjacent to positively charged residues (His/Lys/Arg) at the interface so protonation at low pH repels and weakens binding; **(2) buried His-containing charged H-bond networks** in the binder core so protonation destabilizes the monomer — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740); [biorxiv 2025.09.29.678932](https://www.biorxiv.org/content/10.1101/2025.09.29.678932)
- Pipeline: RFdiffusion (5,000–10,000 backbones) → **ProteinMPNN with His bias (0.5/1.0/1.5)** → AF2 filter (**pae_interaction <8, pLDDT >88**) → partial diffusion + FastRelax → select on ddG, SAP, custom pH-scoring; a **PyRosetta script installs His networks at designated burial depths** + local redesign — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740)
- Targets and pH 7.4→5.4 fold-changes: **EphA2 up to 1,000× (interface); TNFR2 122× (interface); TNFα 79× (network); IL-6 6× (network); PCSK9 3.5× (network); Neo2 >2× (network, 7.4→6.0)** — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740)
- Success rates: TNFR2 6/43, IL-6 8/36, PCSK9 7/40, Neo2 2/48 — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740)
- Availability: **no code/data release stated**; provisional patent filed — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740)

**Other cytokine de-novo binders**
- RFdiffusion on **IL-7Rα**: success ~2 orders of magnitude above the original campaign (1 order from RFdiffusion, 1 from AF2 filtering) — [emergentmind / Watson-Bennett summary](https://www.emergentmind.com/open-problems/de-novo-binders-to-tnfa)
- IL-2: selective de-novo βγ-binding variants via RFdiffusion + ProteinMPNN + Boltz-1 (open "no-code" pipeline) — [Superbio IL-2 writeup](https://superbioai.notion.site/Designing-a-Selective-IL-2-Variant-with-Superbio-s-Protein-Design-Toolkit-PPI-Validation-with-B-1f26c2eab984806f8b67c3c7927402bb)
- AlphaProteo cytokine successes: IL-7Rα, IL-17A, VEGF-A, PD-L1 (binders obtained) — [arXiv 2409.08022](https://arxiv.org/abs/2409.08022)

**Antibody-alternative scaffolds for cytokines**
- Minibinders (helical bundles) are the dominant antibody-alternative and work on TNF-α (above). Germinal (nanobody/antibody, 4–22%) and Chai-2 (antibody 16%, miniprotein 68%) are the generative antibody routes; EvoBind2 gives cyclic/linear peptides. DiffPepBuilder is a target-specific de-novo peptide binder generator — [arXiv 2405.00128](https://arxiv.org/pdf/2405.00128).

### Inferences
- The project's 7.4→6.0 release window is **harder than the Baker paper's main 7.4→5.4 demonstrations**: His pKa ≈ 6, so by pH 6.0 only ~half the His is protonated, giving a weaker switching signal. The only Baker example at 7.4→6.0 (Neo2) achieved just **>2×**. Expect to need the stronger **interface-repulsion** principle (which gave the largest switches, EphA2 1000×/TNFR2 122×) rather than core-network destabilization to hit a usable switch at 6.0.
- TNFR2 (122×, interface) is the closest validated analog in the family and used the interface principle — a template for a TNF-α-side switch on the Region-I hydrophobic epitope.
- Given that classical Rosetta already reached sub-nM on Region I, an all-atom generator (RFd3/BoltzGen) seeded on the same cluster plus an interface His···Lys/Arg repulsion pair is a plausible route to affinity + switch simultaneously.

### Gaps
- No 2025–2026 paper found that combines **TNF-α + de-novo + pH-switch + mouse cross-reactivity** in one (the project's exact spec); the Baker paper (TNF-α pH-switch) and the Comm Biol paper (human/canine cross-reactive, no pH) are the two closest, non-overlapping precedents. Mouse (vs canine) TNF-α cross-reactivity specifically was not found.

---

## Key Question 5 — Availability / license summary (open vs proprietary)

### Takeaway
The best open stack for this project is already largely in hand: RFdiffusion3 + BoltzGen + Boltz-2 + ProteinMPNN + AF2-IG, all open. The two highest-hit-rate generators (AlphaProteo, Chai-2) are proprietary and effectively unavailable; AF3 is usable as an oracle but non-commercial-only.

### Cited Findings
- **Open / runnable**: RFdiffusion1 (open), **RFdiffusion2 & 3** (code+weights via Rosetta Commons Foundry, Dec 2025) — [IPD](https://www.ipd.uw.edu/2025/12/rfdiffusion3-now-available/); **BindCraft** (open) — [GitHub](https://github.com/martinpacesa/BindCraft); **ProteinMPNN** (MIT, open); **AF2-initial-guess / dl_binder_design** (open) — [GitHub nrbennet/dl_binder_design](https://github.com/nrbennet/dl_binder_design); **Boltz-1/2** (MIT) — [boltz.bio](https://boltz.bio/boltz2); **BoltzGen** (MIT) — [GitHub HannesStark/boltzgen](https://github.com/HannesStark/boltzgen); **ProteinGenerator** (open); **EvoBind2** (open; neurosnap/tamarind); **DiffPepBuilder** (open) — [arXiv 2405.00128](https://arxiv.org/pdf/2405.00128).
- **Semi-open**: **AlphaFold3** — code open non-commercial, weights to academics on request, non-commercial only — [maginative](https://www.maginative.com/article/google-deepmind-makes-alphafold-3-open-source-for-academia/); **La-Proteina** (NVIDIA; research release) — [arXiv 2507.09466](https://arxiv.org/abs/2507.09466).
- **Proprietary / closed**: **AlphaProteo** (DeepMind, no code/weights) — [arXiv 2409.08022](https://arxiv.org/abs/2409.08022); **Chai-2** (partner-only, no public API/weights) — [Chai](https://www.chaidiscovery.com/news/introducing-chai-2).
- **Not released**: the Baker **pH-sensitive binders** code (provisional patent) — [PMC12621740](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740); "EpitopeCraft" discontinuous-epitope pipeline (no repo stated) — [arXiv 2509.25479](https://arxiv.org/abs/2509.25479).

### Inferences
- Everything needed to reproduce the Baker pH pipeline is itself open (RFdiffusion, ProteinMPNN His-bias, AF2 filter, PyRosetta, FastRelax) even though their exact code/weights are not released — the project can re-implement the His-bias + His-network-installation steps.
- Chai-1 (open) can substitute for Chai-2's folding/ranking even though Chai-2's generator is closed.

### Gaps
- Germinal and La-Proteina code-release status not explicitly confirmed in search (Germinal likely open — Arc/Stanford; verify repo).
- Chai-1 precise license terms unverified here.

---

## Cross-cutting actionable synthesis (for the report-writer)
1. **Re-target the epitope.** TNF-α's "unbindable" label is specifically the polar subunit-interface patch (AlphaProteo fail). The receptor-binding **Region I hydrophobic cluster** (Trp114/Tyr115/Val74/Leu75/Ile97) yielded sub-nM de-novo minibinders (PDB 9JEC), and is conserved — the strongest lead for both affinity and cross-species.
2. **Use an all-atom generator** (RFd3 already works here; BoltzGen open, benchmark includes TNF-α) and/or **partial-diffuse validated binders**; both beat coarse-backbone de-novo on this class.
3. **Shrink the target to a discontinuous epitope-only template** and add a **hydrophilic-occupancy penalty** in MPNN to satisfy interface polar atoms (direct remedy for "flat polar").
4. **pH switch at 6.0 is hard** (His pKa≈6): prefer the Baker **interface electrostatic-repulsion** principle (His adjacent to Lys/Arg; gave 122–1000× in-family) over core-network destabilization (gave only >2× at 6.0 for Neo2).
5. **Rank by structural confidence + physics, not learned PP-affinity**: AF2-IG pae_interaction/ipTM as primary; ddG/SC/SAP/MM-GBSA orthogonal; Boltz-2 affinity head is unreliable for protein–protein.
