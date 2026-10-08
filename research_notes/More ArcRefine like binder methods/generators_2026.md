# Newest de-novo protein binder GENERATORS (2025-2026): dense catalog + pose-carryover / flat-polar suitability

Scope: NEW or possibly-MISSED de-novo binder generators relevant to a hard, flat/polar, antibody-only target (human TNF-alpha), small 60-120 aa 0-Cys mini-binders, with emphasis on (a) hit rate on flat polar epitopes and (b) ability to carry over / refine an existing validated pose (partial diffusion, template/motif conditioning, affinity maturation). Already in the team's stack (not re-cataloged in depth): RFdiffusion1, RFdiffusion3/all-atom, ProteinMPNN, AF2 initial-guess, Boltz-2, BoltzGen.

Note on "ArcRefine": no public tool by that exact name surfaced in 2025-2026 literature or GitHub this session. The closest documented analogs for "carry over an existing pose and affinity-mature it" are RFdiffusion/ProteinGenerator **partial diffusion**, **AffinityFlow** (alternating structure↔inverse-folding maturation), **ODesign** (partial diffusion of suboptimal binders), **SALAD** (editable per-step denoising), and the Oct-2025 analysis "Mapping the evolution of computationally designed protein binders." See Gaps under Q3.

---

## Q1. What binder generators appeared or updated in 2025-2026 (one-line method / hit rate / pose carryover / availability / paper)

### Takeaway
The 2025-2026 wave splits into four families: (1) structure-prediction-**inversion** binder designers (BoltzDesign1, PXDesign/Protenix, ODesign) that generalize AF3-class models to binder design and now report 20-73% nanomolar hit rates on standard targets; (2) **antibody/CDR** generators (Chai-2, Germinal, RFantibody, IgDesign) that are the natural fit for flat/polar antibody-only epitopes but output Cys-containing Ab/VHH formats, not 0-Cys mini-binders; (3) **all-atom backbone** generators (RFdiffusion3, La-Proteina, Pallatom, SALAD, Genie3, Proteus, Chroma) for the mini-binder format; and (4) explicit **pose-carryover / affinity-maturation** methods (AffinityFlow, partial diffusion, ODesign). The two strongest "missed" open items for a mini-binder campaign are **PXDesign** (open, beats AlphaProteo) and **BoltzDesign1** (open, MIT). The two strongest overall but closed are **Chai-2** and **AlphaProteo**.

### Cited Findings

**Family 1 — Structure-prediction inversion (AF3/Boltz-class → binder design)**

- **BoltzDesign1** — inverts the all-atom Boltz-2 structure-prediction model for generalized biomolecular binder design (protein, and non-protein targets). MIT license; code at github.com/yehlincho/BoltzDesign1 with a Colab notebook; authors Cho, Pacesa, Zhang, Correia, Ovchinnikov; bioRxiv 2025. Open weights via Boltz. — [HF/BoltzDesign1](https://huggingface.co/JoeyCh/BoltzDesign1); [Colab](https://colab.research.google.com/github/yehlincho/BoltzDesign1/blob/main/Boltzdesign1.ipynb)
- **PXDesign (Protenix-design, ByteDance)** — open-source binder-design suite built on the Protenix (AF3-reproduction) stack; offers BOTH a diffusion generator (PXDesign-d) and a hallucination approach (PXDesign-h), then filters by structure prediction + confidence. Reports **20-73% nanomolar hit rates across five of six diverse targets, surpassing AlphaProteo**. Code at github.com/bytedance/PXDesign, benchmark at github.com/bytedance/PXDesignBench, webserver protenix-server.com, sequences at protenix.github.io/pxdesign. — [ByteDance Seed: pxdesign](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders); [bioRxiv node](https://www.biorxiv.org/node/5109812)
- **Protenix-v2** — ByteDance's updated structure-prediction + biomolecular-design model (superset context for PXDesign). Team memory notes Protenix-v2 is proprietary/403 for weights; PXDesign is the open path. — [ByteDance Seed: protenix-v2](https://seed.bytedance.com/zh/public_papers/protenix-v2-broadening-the-reach-of-structure-prediction-and-biomolecular-design)
- **ODesign** — all-atom generative "world model" for all-to-all biomolecular interaction design; lets you **specify epitopes on arbitrary targets** and generate diverse binder classes with fine-grained control; multimodal (protein/nucleic-acid/small-molecule); claims 2-4 orders-of-magnitude higher design throughput and outperforms modality-specific baselines across 11 tasks. Published Oct 25 2025 (arXiv 2510.22304); no-code webserver at odesign.lglab.ac.cn. Explicitly frames **partial diffusion of "suboptimal binders"** (perturb-and-correct while preserving binding mode) as a use case. — [arXiv 2510.22304](https://arxiv.org/abs/2510.22304v2)

**Family 2 — Antibody / nanobody / CDR generators (natural fit for flat polar epitopes; Cys-containing formats)**

- **Chai-2 (Chai Discovery)** — zero-shot/foundation model that designs fully de-novo CDRs (and miniproteins) from epitope + target structure, no antibody templates. **16-20% experimental hit rate** validated across 52 previously unaddressed antigens at only 20 candidates/target; ≥1 hit for 50% of targets in a single round; **binders for 5/5 miniprotein targets**; ~100x over prior methods (field baseline <0.1%). Announced June 30 2025. Proprietary — web/partnership platform, **no code or weights released**. — [Chai: introducing Chai-2](https://www.chaidiscovery.com/news/introducing-chai-2); [BusinessWire](https://www.businesswire.com/news/home/20250630307418/en/Chai-Discovery-Unveils-Chai-2-Breakthrough-Achieving-Fully-De-Novo-Antibody-Design-With-AI)
- **Germinal** — open generative pipeline for epitope-targeted de-novo nanobodies and scFvs: ColabDesign hallucination → selective AbMPNN redesign → cofolding filter, onto a **user-specified antibody framework**. **4-22% experimental success across four targets at only 43-101 designs each**, nanomolar affinities, good mammalian expression. Open code at github.com/SantiagoMille/germinal (an "OpenGerminal" variant exists). bioRxiv 2025.09.19.677421 (Sep 19 2025). — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.09.19.677421.full.pdf); [GitHub](https://github.com/SantiagoMille/germinal)
- **RFantibody (IPD/Baker)** — structure-based de-novo antibody/VHH/scFv pipeline: antibody-finetuned RFdiffusion (backbone) → ProteinMPNN (sequence) → antibody-finetuned RoseTTAFold2 (in-silico filter), starting from target structure + desired epitope. Free for academic/personal/**commercial** use; code at github.com/RosettaCommons/RFantibody. Published in Nature, Nov 2025. — [IPD: RFantibody in Nature](https://www.ipd.uw.edu/2025/11/rfantibody-in-nature); [GitHub](https://github.com/RosettaCommons/RFantibody)
- **IgDesign (Absci)** — antibody CDR (incl. HCDR3) design by inverse folding (IgMPNN) conditioned on antigen + framework; experimentally validated. GitHub absci/igdesign (training-knowledge; see Gaps for DOI). — see Gaps

**Family 3 — All-atom / backbone generators (the 60-120 aa mini-binder format)**

- **RFdiffusion3 (all-atom)** — already in the team's stack and the method that broke the TNF-alpha wall (all-atom generator → MPNN → af2ig pae 13.0 at ~0.4% hit, per team memory). Noted here as the in-house baseline to beat. — (team memory; not a web source)
- **La-Proteina (NVIDIA)** — partially-latent flow matching: models backbone Cα explicitly, captures sequence + atom detail in per-residue latent variables; generates all-atom structure + sequence unconditionally OR does **atomistic motif scaffolding** (builds a protein around a predefined functional site). SOTA on all-atom co-designability/diversity/validity. arXiv 2507.09466 (Jul 2025); code on GitHub (NVIDIA-Digital-Bio/la-proteina) + NGC, released Sep 10 2025; ICLR 2026. Likely NVIDIA non-commercial research license (verify). — [NVIDIA Research](https://research.nvidia.com/labs/genair/publication/laproteina2025/); [arXiv 2507.09466](https://arxiv.org/abs/2507.09466v1)
- **Pallatom** — all-atom generative model learning the joint distribution of structure + sequence via a dual-track (token + atom) framework with an "atom14" representation; strong designability/diversity/novelty. ICML 2025 (Qu et al.). Follow-up **Pallatom-Ligand** (all-atom diffusion for ligand-binding proteins) at ICLR 2026. Public code repo not confirmed this session (see Gaps). — [ICML 2025 poster](https://icml.cc/virtual/2025/poster/43525); [Pallatom-Ligand ICLR 2026](https://iclr.cc/virtual/2026/poster/10006900)
- **SALAD** — family of "sparse all-atom denoising" diffusion models with sub-quadratic complexity; smaller/faster than SOTA, generates up to 1000 aa; supports designable backbones with specified shapes, **scaffolds for functional motifs**, repeat proteins, multi-state proteins, and **rapid prototyping by editing model outputs at each denoising step** (a per-step carryover handle). bioRxiv 2025.01.31.635780 (Jan 2025); code + params on Zenodo (record 14711580, salad-0.1.0). — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.01.31.635780.full.pdf); [Zenodo 14711580](https://zenodo.org/records/14711580)
- **Genie3** — fast all-atom SE(3)-equivariant diffusion; SOTA on unconditional generation, **motif scaffolding, AND binder design**, retaining equivariant efficiency; Columbia (AlQuraishi) + UC Berkeley + Puxano BV + IT4Innovations. Successor to Genie2 (arXiv 2405.15489, 2024; github.com/aqlaboratory/genie — motif scaffolding at structural-universe scale). Genie3 standalone DOI not confirmed this session (see Gaps). — [Genie2 arXiv 2405.15489](https://arxiv.org/html/2405.15489v1); [EuroHPC: Genie3](https://www.eurohpc-ju.europa.eu/genie-30-optimized-denoising-diffusion-probabilistic-model-protein-design_en)
- **Proteus** — deep diffusion backbone generator using graph-based triangle methods + multi-track interaction, with NO structure-prediction pretraining dependency; high designability/efficiency for de-novo backbones. bioRxiv 2024.02.10.579791; ICML 2024. — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.02.10.579791.full.pdf)
- **Chroma (Generate Biomedicines)** — programmable generative model: composable building blocks → diverse all-atom structures via diffusion + equivariant GNNs + conditional random fields; conditioners for shape/symmetry/substructure. Code at github.com/generatebio/chroma (weight download gated by a free registration token). Nature 2023. — [GitHub: chroma](https://github.com/generatebio/chroma)

**Family 4 — Sequence-space generators & explicit affinity-maturation / carryover methods**

- **AffinityFlow** — guided flow matching for **antibody affinity maturation from sequence only**: alternates (1) sequence-fixed, affinity-guided structure generation with (2) inverse folding to propose mutations, plus a co-teaching module using noisy biophysical energies; SOTA on affinity-maturation benchmarks. arXiv 2502.10365, ICML 2025 (Chen et al.). Code to be open-sourced "after acceptance" (verify availability). The clearest published "take an existing binder and mature it" analog. — [arXiv 2502.10365](https://arxiv.org/abs/2502.10365); [ICML 2025](https://proceedings.mlr.press/v267/chen25y.html)
- **ProteinGenerator (Baker lab, Lisanza et al.)** — RoseTTAFold-based diffusion operating in **sequence space**, jointly generating sequence + structure; supports partial diffusion, motif/sequence inpainting, and amino-acid **composition biasing** (useful for 0-Cys / polar-surface constraints). bioRxiv 2023.09.11.556673; code RosettaCommons/protein_generator (training-knowledge for repo). — [bioRxiv 2023.09.11.556673](https://biorxiv.org/cgi/content/short/2023.09.11.556673v1)
- **EvoDiff (Microsoft)** — evolutionary-scale discrete diffusion in **sequence space**; controllable generation, scaffolds functional motifs, and can design regions inaccessible to structure-based models (e.g., disordered). Open: github.com/microsoft/evodiff + PyPI `evodiff` (MIT). Weak fit for structured flat-epitope binders but useful for sequence-level diversification/inpainting. eLife reviewed preprint 112029 (orig. bioRxiv 2023). — [Microsoft Research](https://www.microsoft.com/en-us/research/publication/protein-generation-with-evolutionary-diffusion-sequence-is-all-you-need/); [PyPI](https://pypi.org/project/evodiff/)
- **BindCraft** — one-shot de-novo binder design by AlphaFold2 backpropagation/hallucination; experimental success 10-100%, nanomolar, even without known binding site; does de-novo AND fixed-backbone/target modes. Open source; bioRxiv 2024.09.30.615802; published Nature 2025. — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.09.30.615802.full.pdf); [RosettaCommons blog](https://rosettacommons.org/2025/09/25/one-shot-design-of-functional-protein-binders-with-bindcraft/)
- **Protein Hunter** — exploits structure hallucination WITHIN a diffusion structure-prediction model to design binders from out-of-distribution inputs; supports cyclic-peptide binders and **partial redesign**. bioRxiv 2025.10.10.681530 (Oct 10 2025), new/possibly-missed. — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.10.10.681530.full.pdf)
- **EvoBind2 / AfCycDesign** — AF2/AF-Multimer two-stage sequence optimization with a cyclic-offset trick for (cyclic) peptide binders; ~3x higher success vs prior; peptide-scale, not mini-protein. — [AfCycDesign README (RosettaCommons)](https://huggingface.co/datasets/RosettaCommons/AfCycDesign/blob/main/README.md)
- **AlphaProteo (Google DeepMind)** — generative binder designer; **9-88% success on 7/8 targets**, IL-17A ~700x over prior, BHRF1 88%, best binders at picomolar KD on 4 targets; binding affinity 3-300x over prior methods. Announced Sep 5 2024. **Closed** — academic/industry access only via a limited Trusted Tester Program; commercial work at Isomorphic Labs; no open weights. — [DeepMind/Gigazine coverage](https://gigazine.net/gsc_news/en/20240906-google-deepmind-alphaproteo/)
- **"Mapping the evolution of computationally designed protein binders"** — Chang Liu, bioRxiv 2025.10.04.680454 (Oct 4 2025); an analysis/review of designed-binder development (likely a cross-method benchmark and the best single overview of the "evolution/maturation" theme). Full text not extractable this session (binary PDF). — [DOI 10.1101/2025.10.04.680454](https://doi.org/10.1101/2025.10.04.680454)

### Inferences
- The field's headline hit rates (Chai-2 16-20%, AlphaProteo 9-88%, PXDesign 20-73%) are reported on each group's own benchmark targets, which are generally NOT flat/polar antibody-only targets like TNF-alpha; they should be read as ceilings on "nice" targets, not as TNF-alpha expectations.
- PXDesign and BoltzDesign1 are the most actionable NEW/MISSED **open** mini-binder generators for this campaign: both are structure-prediction-inversion designers (same conceptual family as the team's BoltzGen/Boltz-2), open-licensed, and PXDesign explicitly claims to beat AlphaProteo.
- The all-atom generators (La-Proteina, Pallatom, SALAD, Genie3) are the natural reinforcement for the RFdiffusion3 all-atom approach that already broke the TNF wall; they add motif-scaffolding/editable-denoising handles the team can use for the hotspot-cluster strategy noted in memory.

### Gaps
- Pallatom public code repo/license not confirmed this session (ICML/ICLR pages only).
- Genie3 standalone preprint DOI/date and license not confirmed (only Genie2 arXiv 2405.15489 + aqlaboratory/genie verified).
- IgDesign exact DOI/repo not re-verified this session (training knowledge: Absci, Shanehsazzadeh et al., bioRxiv 2023, github absci/igdesign).
- ODesign code/weights/license (vs webserver-only) not confirmed; appears webserver-first.
- Full benchmark numbers from the Chang Liu "Mapping the evolution…" paper could not be extracted (binary PDF); worth a manual read for cross-method hit-rate comparisons.

---

## Q2. Which are strongest on FLAT / polar / antibody-only epitopes like TNF-alpha?

### Takeaway
For a flat/polar antibody-only epitope, the antibody/CDR generators are mechanistically best-suited (they design concave CDR-loop paratopes that engage flat surfaces the way real antibodies do) — Chai-2 (closed, best hit rate), RFantibody (open), Germinal (open), IgDesign — but they output Cys-containing VHH/scFv formats, not 0-Cys 60-120 aa mini-binders. For the required mini-binder format, the evidence points to all-atom generators (the in-house RFdiffusion3 route, reinforced by La-Proteina/Pallatom/SALAD) plus the inversion designers PXDesign and BoltzDesign1; AlphaProteo's IL-17A result is the main public proof that a generator can crack a hard inflammatory-cytokine target.

### Cited Findings
- Chai-2 designs fully de-novo CDRs from epitope + structure with no templates and hits 16-20% across 52 antigens — i.e., it is explicitly built for the epitope-targeted antibody regime that flat polar targets demand. — [Chai](https://www.chaidiscovery.com/news/introducing-chai-2)
- Germinal is epitope-targeted onto a user-specified framework and succeeded (4-22%) on four diverse protein targets with nanomolar nanobodies — antibody-format, CDR-loop engagement. — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.19.677421.full.pdf)
- RFantibody designs VHHs/scFvs from a target structure + a **desired epitope**, i.e., you can point it at the TNF epitope directly. — [IPD](https://www.ipd.uw.edu/2025/11/rfantibody-in-nature)
- AlphaProteo achieved a ~700x success-rate improvement on **IL-17A** (an inflammatory cytokine, a genuinely hard target class), with picomolar binders on four targets — the strongest public evidence a mini-binder generator can crack a hard cytokine. — [Gigazine/DeepMind](https://gigazine.net/gsc_news/en/20240906-google-deepmind-alphaproteo/)
- PXDesign reports 20-73% nanomolar hit rates across 5/6 diverse targets and claims to surpass AlphaProteo, using both diffusion and hallucination generators with structure-confidence filtering. — [ByteDance Seed](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders)
- ODesign allows specifying epitopes on **arbitrary** targets and generating diverse binder classes — aimed exactly at controllable epitope targeting. — [arXiv 2510.22304](https://arxiv.org/abs/2510.22304v2)

### Inferences
- The team's own result (pure de-novo RFd1/BoltzGen floored at af2ig pae ~22 on TNF; only RFdiffusion3 all-atom broke it) strongly implies that for flat polar epitopes, **all-atom generation matters more than the specific codebase** — so the highest-value additions are other all-atom generators (La-Proteina, Pallatom, SALAD, Genie3) and all-atom inversion designers (PXDesign, BoltzDesign1), not more backbone-only RFd1-class tools.
- Antibody/CDR methods (Chai-2/RFantibody/Germinal/IgDesign) likely give the highest raw binding hit rate on TNF-alpha, but they conflict with the 0-Cys 60-120 aa mini-binder constraint (VHH ≈ 120 aa with 1-2 conserved disulfides; scFv larger). Use them as a mechanistic hedge or to mine binding modes/epitope poses, not as drop-in mini-binder generators.
- IL-17A (AlphaProteo) and TNF-alpha are both trimeric inflammatory cytokines with largely flat/polar surfaces; AlphaProteo's IL-17A success is the closest public precedent and argues the target class is crackable by a strong generator + heavy filtering.

### Gaps
- No source reports a binder-generator hit rate specifically on TNF-alpha (the team's own runs are the only TNF-specific data).
- Whether Chai-2's miniprotein mode (5/5 miniprotein targets) can be pointed at a flat polar cytokine epitope at mini-binder scale is not documented; miniprotein targets in their report are unspecified here.

---

## Q3. Which support carrying over / refining an existing pose (vs fresh generation only)?

### Takeaway
Pose carryover is well-supported across the open stack, via three mechanisms: (1) **partial diffusion** — noise an existing binder a little, re-denoise, re-sequence (RFdiffusion1/3 [in-house], ProteinGenerator, ODesign, SALAD editable denoising, La-Proteina/Genie/Pallatom motif scaffolding); (2) **template / framework conditioning** (RFantibody framework+epitope, IgDesign framework+antigen, BoltzDesign1/BoltzGen template, Chai-2 epitope conditioning); and (3) **explicit affinity maturation of a known binder** (AffinityFlow, EvoBind2 sequence optimization, and the theme of the Chang-Liu "Mapping the evolution…" paper). AffinityFlow is the closest published "ArcRefine-style" method but is antibody/sequence-only and code may be unreleased.

### Cited Findings
- ODesign explicitly frames partial diffusion of "suboptimal binders": treat an existing binder as the initial molecule, introduce controlled perturbations/corrections within the diffusion/inverse-folding framework to explore optimization space **while preserving the binding mode**. — [arXiv 2510.22304](https://arxiv.org/abs/2510.22304v2)
- AffinityFlow matures an existing antibody by alternating sequence-fixed affinity-guided structure generation with inverse-folding mutation proposals (co-teaching on biophysical energies) — i.e., it carries over a pose and optimizes affinity. — [arXiv 2502.10365](https://arxiv.org/abs/2502.10365)
- SALAD enables rapid prototyping by **editing model outputs at each denoising step**, and scaffolds functional motifs / multi-state proteins — a direct handle for seeding from and refining an existing structure. — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.31.635780.full.pdf)
- La-Proteina performs **atomistic motif scaffolding** (build around a predefined functional site), the all-atom analog of grafting a known interface. — [NVIDIA Research](https://research.nvidia.com/labs/genair/publication/laproteina2025/)
- Genie3/Genie2 support **motif scaffolding** (steer generation to contain a desired functional substructure). — [Genie2 arXiv](https://arxiv.org/html/2405.15489v1)
- RFantibody conditions on a specified framework + epitope (carry over the framework/epitope), and Germinal designs CDRs onto a **user-specified structural framework** — both carry over scaffold while regenerating the interface. — [IPD](https://www.ipd.uw.edu/2025/11/rfantibody-in-nature); [Germinal bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.19.677421.full.pdf)
- BindCraft can run de-novo OR fixed-backbone/target design, giving control to preserve an existing backbone. — [RosettaCommons](https://rosettacommons.org/2025/09/25/one-shot-design-of-functional-protein-binders-with-bindcraft/)
- Protein Hunter supports **partial redesign** (plus cyclic peptides) via diffusion hallucination. — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.10.10.681530.full.pdf)
- ProteinGenerator does sequence-space diffusion with motif/sequence inpainting and composition biasing (partial/conditioned regeneration from an existing design). — [bioRxiv 2023.09.11.556673](https://biorxiv.org/cgi/content/short/2023.09.11.556673v1)

### Inferences
- The team already owns the single most practical carryover lever: **RFdiffusion partial diffusion**, which team memory records as ~11% pass vs ~0.3% de-novo when partial-diffusing VALIDATED binders. The new tools above are complements, not replacements: ODesign/SALAD/ProteinGenerator add perturb-and-preserve modes in different model families; AffinityFlow adds an explicit affinity objective.
- For "ArcRefine-like" refinement of a VALIDATED TNF mini-binder pose, the lowest-friction route is RFdiffusion3 partial diffusion + ProteinMPNN re-seq + af2ig (already in stack), optionally cross-checked with BoltzDesign1/PXDesign partial/template modes. AffinityFlow's objective is attractive but it targets antibody sequence space, not 0-Cys mini-binders, and code availability is uncertain.

### Gaps
- "ArcRefine" is not a locatable public tool/paper (no GitHub/preprint match this session); it may be an internal/codename or a very obscure method. Functional analogs listed above.
- AffinityFlow code release status ("after acceptance") unconfirmed — may not be runnable yet.
- Whether BoltzDesign1 / PXDesign expose a first-class "seed from existing binder" (partial) mode vs only template conditioning was not confirmed from primary docs this session.

---

## Q4. Open + runnable without special licenses (uv + https-only Linux/H100)?

### Takeaway
Open and installable: BoltzDesign1 (MIT), PXDesign (open), Germinal, RFantibody (free incl. commercial; Docker/PyRosetta friction), BindCraft, ProteinGenerator, EvoDiff (MIT/PyPI), SALAD, Genie2/3, Proteus, Chroma (free token to fetch weights). License-flagged / friction: La-Proteina (NVIDIA, likely non-commercial research license; NGC/HF). Closed / not locally installable: Chai-2, AlphaProteo, Protenix-v2 weights, ODesign (webserver-first), AffinityFlow (code possibly unreleased).

### Cited Findings
- BoltzDesign1: MIT license, GitHub + Colab, built on open Boltz weights. — [HF/BoltzDesign1](https://huggingface.co/JoeyCh/BoltzDesign1)
- PXDesign: full pipeline open at github.com/bytedance/PXDesign + benchmark repo + webserver + public sequences. — [ByteDance Seed](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders)
- Germinal: open code + full computational/experimental protocols at github.com/SantiagoMille/germinal (plus an OpenGerminal variant). — [GitHub](https://github.com/SantiagoMille/germinal)
- RFantibody: shared freely for academic, personal, AND commercial use; github.com/RosettaCommons/RFantibody. (Pipeline uses RF2/PyRosetta components — Docker-based install, a wrinkle for a pure uv environment.) — [IPD](https://www.ipd.uw.edu/2025/11/rfantibody-in-nature)
- EvoDiff: github.com/microsoft/evodiff + PyPI `evodiff`. — [PyPI](https://pypi.org/project/evodiff/)
- Chroma: github.com/generatebio/chroma (weights behind a free registration token). — [GitHub](https://github.com/generatebio/chroma)
- La-Proteina: NVIDIA code on GitHub + NGC (research model; NVIDIA model licenses are typically non-commercial — verify before use). — [NVIDIA Research](https://research.nvidia.com/labs/genair/publication/laproteina2025/)
- Chai-2 (web/partnership only), AlphaProteo (Trusted Tester Program only; Isomorphic for commercial): no open weights/code. — [Chai](https://www.chaidiscovery.com/news/introducing-chai-2); [Gigazine/DeepMind](https://gigazine.net/gsc_news/en/20240906-google-deepmind-alphaproteo/)
- ODesign: available as a no-code webserver (odesign.lglab.ac.cn); downloadable code/weights not confirmed. — [arXiv 2510.22304](https://arxiv.org/abs/2510.22304v2)

### Inferences
- For a uv + https-only H100 setup, the frictionless adds are **BoltzDesign1, PXDesign, BindCraft, ProteinGenerator, EvoDiff, SALAD, Genie2/3, Proteus** (pip/conda + weights over HTTPS from HF/Zenodo/GitHub). All should run on a single H100.
- **RFantibody** is open but ships as a Docker image with PyRosetta/RF2 dependencies — budget container/pyrosetta setup time; not a clean `uv pip install`.
- **La-Proteina** is open but gated behind NVIDIA licensing and NGC; confirm the license permits competition use before investing.
- **Chai-2 and AlphaProteo** — the two highest-reported hit rates — are NOT locally runnable; access is via platform/partnership or a trusted-tester program. Treat them as unavailable for an autonomous in-house pipeline.

### Gaps
- Exact licenses unconfirmed this session for: La-Proteina (NVIDIA terms), Pallatom (repo/license), Genie3 (repo/license), ProteinGenerator and IgDesign repos (training-knowledge).
- ODesign and AffinityFlow: downloadable-code availability (vs webserver / unreleased) not confirmed.
- Did not test actual uv-installability of any package this session (no installs run); claims are about stated availability, not verified builds on this H100.
