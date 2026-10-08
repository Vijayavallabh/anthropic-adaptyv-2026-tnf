# In-silico Affinity Maturation and Pose-Preserving Binder Optimization / Ranking Methods (2024–2026)

Scope: methods to mature the affinity of an existing de-novo protein mini-binder **without losing the docked pose** (structural-carryover, "ArcRefine" use-case). Emphasis on (a) ML ddG / affinity predictors for ranking interface mutations on a fixed backbone, (b) in-silico directed-evolution / active-learning loops, (c) what is open-source + PyRosetta-free + uv/pip/https-installable, and (d) published correlation to measured Kd. Context: user already runs AF2-initial-guess (pae_interaction 4–8), Boltz-2 ipTM, ProteinMPNN; wants MISSED/NEWER tools; needs the His pH-switch preserved; cross-reactivity (human + mouse) ranking a plus.

---

## KQ1: 2024–2026 ML affinity / ddG predictors usable for ranking interface mutations on a fixed backbone

### Takeaway
There is now a mature set of fast ML binding-ddG predictors that operate on a **fixed complex backbone** (exactly the carryover use-case): the strongest PyRosetta-free options are **StaB-ddG** (ProteinMPNN-based, matches FoldX, CC-BY, pip-installable) and **GearBind** (geometric GNN, wet-lab validated up to 17-fold affinity gains). Crucially, the **Boltz-2 affinity head is small-molecule-focused and is weak on protein-protein affinity** (Pearson 0.15 on TCR3d), so it is not a reliable ranker for a protein mini-binder.

### Cited Findings
- **StaB-ddG** predicts mutational ΔΔG_bind from interface structure+sequence on a **fixed backbone**; it is pure deep learning built on ProteinMPNN (structure-conditioned log-likelihoods as a folding-energy proxy), **requires neither PyRosetta nor FoldX**, enforces thermodynamic antisymmetry and path-independence, and is the first DL predictor to match force-field FoldX on binding ΔΔG — [arXiv:2507.05502 (deposited 7 Jul 2025)](https://arxiv.org/html/2507.05502v1).
- StaB-ddG per-interface Spearman on SKEMPI v2 test split = **0.45** (not statistically different from FoldX or Flex ddG); **ensemble StaB-ddG+FoldX = 0.53 Spearman (SOTA)**; speed ~**0.2 s/mutation on GPU, ~1000× faster than FoldX** (~260 CPU-h). Code: **github.com/LDeng0205/StaB-ddG, CC BY 4.0** — [arXiv:2507.05502](https://arxiv.org/html/2507.05502v1).
- **GearBind**: pretrainable geometric graph NN (multi-relational, multi-level message passing on GearNet, contrastive pretraining on CATH, fine-tuned on SKEMPI v2 with antisymmetric predictor) for ΔΔG_bind; surpasses prior SOTA on SKEMPI v2 (Pearson, MAE, RMSE). Wet-lab validated: increased **CR3022 affinity vs SARS-CoV-2 Omicron spike up to 17-fold** and a human single-domain antibody (UdAb) vs 5T4 **up to 5.6-fold** — [Nat. Commun. 15, 2024, s41467-024-51563-8](https://doaj.org/article/326a856cdb144984ae87d5a5ddcdd218); [GearBind biorxiv PDF](https://www.biorxiv.org/content/10.1101/2023.08.10.552845.full.pdf). Code: **github.com/DeepGraphLearning/GearBind** (TorchDrug/PyTorch, Mila) — [GitHub](https://github.com/DeepGraphLearning/GearBind).
- **Boltzmann-Aligned inverse-folding ddG (BA-DDG)**: transfers an inverse-folding model to ΔΔG via Boltzmann distribution + Bayes, **explicitly modeling the unbound state** in the thermodynamic cycle. SKEMPI v2 Spearman **0.3201 (unsupervised) / 0.5134 (supervised)**, beating prior SOTA (0.2632 / 0.4324). Demonstrated on binding-energy prediction, docking, and **antibody optimization** — [ICLR 2025, arXiv:2410.09543](https://arxiv.org/html/2410.09543v1).
- **Boltz-2** (MIT license, open, co-folds structure + predicts affinity) approaches FEP accuracy for **small-molecule** Kd/IC50 (~0.6 correlation on FEP+ benchmark; beat all CASP16 affinity-track entries) at ~18 s/GPU, ≥1000× faster than FEP — [Boltz-2 biorxiv 2025.06.14.659707](https://www.biorxiv.org/content/10.1101/2025.06.14.659707.full.pdf); [PMC12262699](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699/).
- **Boltz-2 is weak for protein-protein affinity**: Boltz-2-PPI on TCR3d gave **Pearson 0.153 / Spearman 0.091**, *worse* than an ESM2-650M sequence baseline (0.239) — [PMC12262699](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699/).
- **ThermoMPNN**: transfer-learning ddG predictor using frozen ProteinMPNN embeddings + the ~300-protein "Megascale" stability dataset; predicts **monomer stability ΔΔG** (not binding) for all point mutations in seconds–minutes; ThermoMPNN-D handles double mutants. PyTorch, no PyRosetta. Code: github.com/Kuhlman-Lab/ThermoMPNN — [biorxiv 2023.07.27.550881](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10402116/); [UNC Kuhlman lab](https://klab.web.unc.edu/predicting-mutation-effects-on-protein-stability/).
- **ProteinMPNN-ddG** (Peptone): aligns ProteinMPNN outputs for stability prediction using a second term derived from the model itself, **no retraining/no extra data**, improving ProteinMPNN's success at picking stabilizing variants by **11%**; has a fast JAX implementation. Monomer stability focus. Code: **github.com/PeptoneLtd/proteinmpnn_ddg** — [biorxiv 2024.06.15.599145](https://www.biorxiv.org/content/10.1101/2024.06.15.599145v2); [NeurIPS 2024](https://nips.cc/virtual/2024/102469).
- **mCSM-PPI2**: graph-based ML (inter-residue network signatures + evolutionary + energetic terms) for single-point PPI ΔΔG; ranked 1st of 27 on CAPRI blind tests; **webserver only** at biosig.unimelb.edu.au/mcsm_ppi2 — [NAR 2019, PMC6602427](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602427/).
- **DDMut-PPI**: graph convolutional deep-learning model for single + multiple-point PPI ΔΔG_bind (Siamese + graph signatures); **webserver** (biosig) — [UoM listing](https://findanexpert.unimelb.edu.au/scholarlywork/1891744); mmCSM-PPI handles multiple mutations — [PMC8262703](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8262703/).
- **Twin Peaks**: dual-head, **structure-free** (sequence-only) prediction of PPI binding affinity and mutation effects — [arXiv:2509.22950 (Sep 2025)](https://arxiv.org/pdf/2509.22950).

### Inferences
- For a **protein** TNF-α mini-binder, the dedicated binding-ddG models (StaB-ddG, GearBind, BA-DDG) are the right rankers; the Boltz-2 affinity head should **not** be used as the affinity signal here (it is tuned to small-molecule Kd and under-performs ESM2 on PPI). Keep using Boltz-2 only for ipTM/pose confidence.
- StaB-ddG is the single best fit for ArcRefine: fixed-backbone, PyRosetta-free, CC-BY, GPU-fast, and its FoldX-ensemble (0.53 Spearman) is the current accuracy ceiling for cheap binding-ddG ranking. GearBind adds independent, wet-lab-validated signal — an **ensemble/consensus of StaB-ddG + GearBind + FoldX** is the robust ranking strategy.
- ThermoMPNN / ProteinMPNN-ddG are **stability** predictors, not binding — use them as a *developability/foldability guardrail* (reject interface mutations that destabilize the monomer) rather than as the affinity signal.

### Gaps
- No confirmed public code repo found for BA-DDG (arXiv:2410.09543) or for Twin Peaks (arXiv:2509.22950) — availability/license unverified.
- Did not find head-to-head numbers for StaB-ddG vs GearBind on an identical split.

---

## KQ2: In-silico directed evolution / active-learning / Bayesian-optimization loops for binder affinity

### Takeaway
Two complementary loop families are available and PyRosetta-free: (a) **gradient/MCMC sampling over sequence** driven by any differentiable fitness (EvoProtGrad), and (b) **active-learning / learning-to-rank** loops that call a ddG or co-fold "oracle" iteratively (ALLM-Ab, antibody AL). AlphaFold-as-oracle directed evolution (EvoBind2) reports validated nanomolar–picomolar peptide binders, establishing the co-folder-oracle loop as real.

### Cited Findings
- **EvoProtGrad**: plug-and-play directed evolution via **gradient-based discrete MCMC**; wraps *any* neural network as a black-box fitness predictor (product-of-experts in discrete sequence space) with **no fine-tuning**, using gradients to propose beneficial mutations; efficiently reaches high-likelihood/high-activity variants multiple mutations from WT. Code: github.com/NREL/EvoProtGrad — [Mach. Learn.: Sci. Technol. 2023 / arXiv:2212.09925](https://arxiv.org/abs/2212.09925v1).
- **EvoBind / EvoBind2**: AlphaFold-based **in-silico directed evolution** of peptide binders from target sequence only, adapting the receptor interface during optimization and scoring with plDDT/AF confidence. EvoBind2 reports success rates **46% (linear), 75% (cyclic)** with affinities **0.26 nM – 7.9 µM**; reusing AF as a scoring function reportedly raised experimental success to ~90% vs 0–5% for physics-based scoring. Code: github.com/patrickbryant1/EvoBind — [EvoBind biorxiv 2022.07.23.501214](https://www.biorxiv.org/content/10.1101/2022.07.23.501214.full.pdf); [EvoBind2 (Neurosnap)](https://neurosnap.ai/service/EvoBind2).
- **ALLM-Ab**: active-learning antibody optimization using **LoRA fine-tuned protein language models + learning-to-rank**; validated offline on deep-mutational-scanning data and in **online active-learning trials minimizing Flex ddG across 15 antigens**, finding high-affinity variants faster than baselines **while preserving developability** — [biorxiv 2025.08.05.668775 (Aug 2025)](https://www.biorxiv.org/content/10.1101/2025.08.05.668775.full.pdf); [PMC12606632](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12606632/).
- **Active learning for antibody affinity**: iterative pool/oracle loop selecting batches for labeling to data-efficiently improve affinity prediction — [arXiv:2406.07263 (Jun 2024)](https://www.arxiv.org/pdf/2406.07263).
- **Bayesian optimization + physics oracle** (small-molecule analogue of the loop): BO surrogate selects a subset each round for MMGBSA/MD, "focusing on promising regions" to get MD-level affinity at scale — [biorxiv 2025.06.22.660936](https://www.biorxiv.org/content/10.1101/2025.06.22.660936.full.pdf).
- **EvoPro**: genetic-algorithm protein optimization using AF2 scoring (Kuhlman lab) — code github.com/Kuhlman-Lab/evopro — [GitHub](https://github.com/Kuhlman-Lab/evopro).

### Inferences
- The highest-leverage, lowest-effort loop for this project: **EvoProtGrad with a binding-ddG model (StaB-ddG/GearBind) as the fitness function**, constrained to interface positions and excluding the His-switch residues — this is pose-preserving by construction (sequence-only moves on the fixed backbone) and PyRosetta-free.
- ALLM-Ab's design (PLM + learning-to-rank + ddG oracle + developability constraint) is a directly transferable template for a cross-reactive, developability-aware maturation loop, though it is antibody-framed and not a turnkey mini-binder tool.
- EvoBind2 shows an AF-oracle loop reaching picomolar peptides, but it co-adapts the receptor interface (not strictly fixed-pose) and is peptide-oriented.

### Gaps
- No benchmark found that uses a **co-folder (AF2-ig/Boltz) ipTM as the oracle inside a Bayesian-optimization loop specifically for de-novo mini-protein affinity** (the closest are EvoBind's AF-plDDT loop for peptides and the MMGBSA-BO for small molecules).

---

## KQ3: Methods specifically for pose-preserving (fixed backbone) affinity maturation

### Takeaway
Pose-preserving maturation = rank/optimize interface **sidechains and sequence on a frozen backbone**, then verify the pose is retained with a co-folder. The field splits into (a) fixed-backbone ddG ranking (KQ1 tools) and (b) "light-touch" generative refinement that keeps the interface and only relaxes neighboring regions (interface-rotamer-enrichment + partial diffusion/flow). Pure fixed-backbone ddG ranking is the safest for keeping both the dock and the His switch.

### Cited Findings
- Partial diffusion combined with ProteinMPNN yields **significantly more designable** proteins than ProteinMPNN sequence redesign alone on the native backbone; partial-diffusion specificity reprogramming produced a **record 10 pM monomeric binder** — [Model-Centric Review of DL for Protein Design, arXiv:2502.19173](https://arxiv.org/pdf/2502.19173).
- A described **in-silico maturation strategy**: interface rotamer enrichment, selection/merging of key interface residues, **fixation of critical residues with noise added only to neighboring backbone**, then partial flow-based refinement of unconstrained regions to optimize interfacial packing — [search synthesis; see PPIFlow below](https://arxiv.org/pdf/2502.19173).
- **PPIFlow + in-silico affinity maturation** (newest, Jan 2026): pairformer flow-matching generator that models backbone rigid-body transforms as continuous flows; maturation "combines **interface rotamer enrichment with partial flow refinement** to optimize energetic packing," prioritizing candidates with **AF3Score** (a score-only AlphaFold3 adaptation). Reports consistent picomolar/nanomolar binders and **sub-nanomolar VHHs**. Note: it **refines the backbone** (not strictly fixed-pose). Changping Lab (Mingchen Chen), posted 22 Jan 2026 — [biorxiv 10.64898/2026.01.19.700484](https://www.biorxiv.org/content/10.64898/2026.01.19.700484v1).
- Fixed-backbone design models (ProteinMPNN family) are noted to **not natively address affinity optimization, off-target reduction, or accurate interface side-chain conformations**, motivating dedicated interface-maturation steps — [arXiv:2502.19173](https://arxiv.org/pdf/2502.19173).
- **FoldX / mutateX** give fixed-backbone, physics ΔΔG for saturation mutagenesis of a complex; mutateX is a Python wrapper automating per-residue saturation (incl. multimers, ensembles, PTMs). FoldX is a free-for-academia binary (not PyRosetta). Code: github.com/ELELAB/mutatex — [mutateX, Brief. Bioinform. 2022 bbac074](https://academic.oup.com/bib/article/23/3/bbac074/6552273); [FoldX suite](https://foldxsuite.crg.eu/node/3337).

### Inferences
- The cleanest ArcRefine recipe that preserves both pose and the His switch: **(1)** enumerate interface point/double mutations on the frozen backbone (exclude His + its salt-bridge partner), **(2)** rank by a binding-ddG ensemble (StaB-ddG + GearBind + FoldX), **(3)** gate on monomer stability (ThermoMPNN) for developability, **(4)** verify the top candidates still dock with AF2-ig/Boltz-2 ipTM (pose carryover check). No single turnkey tool does all four, but each step is PyRosetta-free.
- FoldX has an explicit **pH parameter**, so it can score ΔΔG at pH 7.4 vs pH 6.0 — uniquely useful for not breaking the His pH-switch while maturing affinity (physics, no PyRosetta).
- Generative refinement (PPIFlow / partial diffusion) risks moving the backbone and therefore the pose/His geometry; prefer it only if fixed-backbone ranking plateaus, and re-validate the switch geometry after.

### Gaps
- PPIFlow code availability/license and PyRosetta dependence were **not stated** in the abstract; no measured Kd table was visible in the excerpt.

---

## KQ4: Which are open and runnable WITHOUT PyRosetta, within uv/pip + https-only Linux/H100

### Takeaway
Most modern ML tools are pure PyTorch/JAX and pip/uv-installable; the hard **PyRosetta/Rosetta blockers are Flex ddG and Cartesian ddG**. FoldX is PyRosetta-free but a proprietary (free-academic) binary; mCSM-PPI2 and DDMut-PPI are webserver-only.

### Cited Findings
- **PyRosetta-free, pip/uv-installable (PyTorch/JAX), open license**: StaB-ddG (CC-BY) — [arXiv:2507.05502](https://arxiv.org/html/2507.05502v1); GearBind (PyTorch/TorchDrug) — [GitHub](https://github.com/DeepGraphLearning/GearBind); ProteinMPNN-ddG (JAX) — [GitHub](https://www.biorxiv.org/content/10.1101/2024.06.15.599145v2); ThermoMPNN (PyTorch) — [biorxiv](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10402116/); Boltz-2 (MIT) — [biorxiv](https://www.biorxiv.org/content/10.1101/2025.06.14.659707.full.pdf); EvoProtGrad — [arXiv:2212.09925](https://arxiv.org/abs/2212.09925v1); EvoBind/EvoBind2 — [biorxiv](https://www.biorxiv.org/content/10.1101/2022.07.23.501214.full.pdf); BindCraft (open-source, AF2-based) — [BindCraft biorxiv 2024.09.30.615802](https://www.biorxiv.org/content/10.1101/2024.09.30.615802.full.pdf).
- **ipSAE** (improved interface scoring from AF2/AF3/Boltz PAE — "what's wrong with ipTM and how to fix it"): **pip-installable** (`pypi.org/project/ipsae`) — [PMC11844409](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11844409/); [PyPI](https://pypi.org/project/ipsae/).
- **PyRosetta/Rosetta REQUIRED (blockers for this project)**: Flex ddG (Rosetta backrub ensemble, multiple CPU-hours/mutation) — [Flex ddG, J. Phys. Chem. B 2018, PMC5980710](https://pmc.ncbi.nlm.nih.gov/articles/PMC5980710); Cartesian ddG / cartddg2020 (Rosetta Cartesian minimization) — [Rosetta docs](https://docs.rosettacommons.org/docs/latest/cartesian-ddG).
- **PyRosetta-free but proprietary binary (free academic license)**: FoldX, wrapped by mutateX — [FoldX](https://foldxsuite.crg.eu/node/3337); [mutateX](https://academic.oup.com/bib/article/23/3/bbac074/6552273).
- **Webserver-only (no documented local install)**: mCSM-PPI2 — [biosig.unimelb.edu.au/mcsm_ppi2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602427/); DDMut-PPI (biosig) — [UoM](https://findanexpert.unimelb.edu.au/scholarlywork/1891744).

### Inferences
- A fully PyRosetta-free maturation stack is feasible on the 8×H100 box: StaB-ddG + GearBind + ProteinMPNN-ddG/ThermoMPNN + EvoProtGrad + AF2-ig/Boltz-2 + ipSAE, all pip/uv.
- The common Flex ddG and Cartesian-ddG baselines are **out** (need PyRosetta). The usable physics option is **FoldX/mutateX** (binary download over https; verify academic-license terms) — valuable as the physics leg of a consensus and for pH-dependent ΔΔG.
- Webserver tools (mCSM-PPI2, DDMut-PPI) are awkward for batch/autonomous loops and send sequences off-box; treat as optional cross-checks, not pipeline components.

### Gaps
- Did not verify whether StaB-ddG/GearBind are published to PyPI (install likely via `pip install git+https://github.com/...`); both repos are public so git-pip install over https should work.

---

## KQ5: Published correlation to measured Kd / expression; cross-reactivity-aware ranking

### Takeaway
Binding-ddG ML models have **moderate** correlation to experimental ΔΔG (SKEMPI Spearman ~0.45–0.53), roughly matching force fields; GearBind and ALLM-Ab add real wet-lab affinity gains. Co-fold confidence (ipTM/iPAE) **discriminates binders from non-binders but is a weak quantitative Kd ranker** — which is the gap a dedicated ddG model fills. No tool natively does human+mouse cross-reactivity ranking; the practical route is dual-ortholog consensus.

### Cited Findings
- SKEMPI v2 ΔΔG correlations: StaB-ddG **0.45** (ensemble with FoldX **0.53**) — [arXiv:2507.05502](https://arxiv.org/html/2507.05502v1); BA-DDG **0.5134 supervised / 0.3201 unsupervised** — [arXiv:2410.09543](https://arxiv.org/html/2410.09543v1); Flex ddG Pearson **≈0.65** (but multiple CPU-hours/mutation) — [arXiv:2507.05502 (reporting Flex ddG)](https://arxiv.org/pdf/2507.05502).
- **Data-leakage caveat**: "when data leakage is corrected, the performance of these deep learning predictors lags behind Flex ddG" — i.e., reported ML gains over force fields can be inflated by train/test overlap — [search synthesis of Flex ddG comparisons](https://www.biorxiv.org/content/10.1101/221689.full.pdf).
- **Wet-lab affinity gains**: GearBind — CR3022 **up to 17-fold**, UdAb **up to 5.6-fold** — [Nat. Commun. 2024](https://doaj.org/article/326a856cdb144984ae87d5a5ddcdd218); EvoBind2 — designed peptides **0.26 nM–7.9 µM** — [Neurosnap](https://neurosnap.ai/service/EvoBind2); ALLM-Ab — faster discovery of high-affinity variants validated against DMS data — [biorxiv 2025.08.05.668775](https://www.biorxiv.org/content/10.1101/2025.08.05.668775.full.pdf).
- Boltz-2 affinity correlation to experiment: **~0.6 on FEP+ (small molecule)**; **0.15 Pearson on protein-protein (TCR3d)** — [PMC12262699](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699/).
- Co-fold confidence vs Kd: binders have **higher average ipTM** than weak/non-binders (binary discrimination), but across **55 PDBbind entries with known Kd the iPAE–Kd correlation is "relatively weak"**; AF3 interface ipTM is often ~0.3 and "not yet fully reliable" for distinguishing binders — [arXiv:2506.01995 (GFRAL antibody discovery)](https://arxiv.org/pdf/2506.01995); [ipSAE/PMC11844409](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11844409/).
- **ipSAE** reweights ipTM to only well-predicted interface residue pairs, improving binder/non-binder separation over raw ipTM/iPAE — [PMC11844409](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11844409/).

### Inferences
- Expected real-world behavior: ML binding-ddG rankers will correctly order *large* affinity differences but are noisy at the single-mutation level (Spearman ~0.5) — so **mature in small accumulating steps and rank by a consensus** (StaB-ddG + GearBind + FoldX), not any single model, and confirm wins experimentally in batches (echoes the project's "small-n win is noise until bootstrapped" lesson).
- Because ipTM/pae_interaction discriminate binders but rank Kd weakly, the value-add here is a **dedicated binding-ddG predictor layered on top of the existing AF2-ig/Boltz-2 pose gate** — use co-fold for pose/binary pass, ddG ensemble for affinity ranking, ipSAE to sharpen the pose gate.
- **Cross-reactivity (human + mouse)**: run the binding-ddG predictor (and co-fold) on **both the human-TNF and mouse-TNF complexes** and require a mutation to be non-deleterious (ideally improving) on **both** — i.e., rank by min(ΔΔG_human, ΔΔG_mouse) or a consensus, and reject mutations that help one species but hurt the other. This is a configuration of the above tools, not a dedicated method.
- FoldX pH-dependent ΔΔG lets you check that a maturation mutation does not flatten the pH 7.4 vs 6.0 ΔΔG gap that encodes the His switch — maturation and switch-preservation can be co-scored.

### Gaps
- No method found that **jointly** optimizes affinity across two orthologs (human+mouse) out of the box — it must be assembled from per-complex ddG runs.
- No published direct correlation of StaB-ddG or GearBind predictions to **expression/developability**; developability must be covered separately (ThermoMPNN stability, aggregation/charge heuristics).
- The data-leakage critique is reported second-hand via the StaB-ddG/Flex ddG discussion; did not independently confirm the specific leakage-corrected benchmark numbers.

---

## Quick reference: method catalog (carryover = fixed-backbone sequence/sidechain optimization)

| Method | Fixed-backbone (carryover) or redesign | Affinity signal | Open / PyRosetta-free / uv-pip | Primary source + date | Exp. affinity evidence |
|---|---|---|---|---|---|
| **StaB-ddG** | Carryover (rank interface mutations) | ML binding ΔΔG (ProteinMPNN folding-energy cycle) | Yes / yes / yes (CC-BY, GitHub) | [arXiv:2507.05502, Jul 2025](https://arxiv.org/html/2507.05502v1) | SKEMPI 0.45 (0.53 w/ FoldX); matches FoldX |
| **GearBind** | Carryover | ML binding ΔΔG (geometric GNN) | Yes / yes / yes (GitHub) | [Nat.Commun. 2024](https://doaj.org/article/326a856cdb144984ae87d5a5ddcdd218) | 17-fold & 5.6-fold wet-lab gains |
| **BA-DDG** | Carryover | ML ΔΔG (Boltzmann-aligned inverse folding, unbound-state term) | Open paper / yes / code unverified | [ICLR 2025, arXiv:2410.09543](https://arxiv.org/html/2410.09543v1) | SKEMPI 0.51 supervised |
| **ProteinMPNN-ddG** | Carryover | ML stability ΔΔG (no retraining; JAX) | Yes / yes / yes (GitHub) | [biorxiv 2024.06.15.599145](https://www.biorxiv.org/content/10.1101/2024.06.15.599145v2) | +11% stabilizing-variant picks |
| **ThermoMPNN** | Carryover | ML stability ΔΔG (Megascale) | Yes / yes / yes (GitHub) | [biorxiv 2023.07.27.550881](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10402116/) | SOTA on stability benchmarks |
| **FoldX / mutateX** | Carryover | Physics ΔΔG (incl. pH) | Free-academic binary / yes / wrapper is pip | [mutateX 2022](https://academic.oup.com/bib/article/23/3/bbac074/6552273) | long-standing physics baseline |
| **Flex ddG** | Carryover | Physics ΔΔG (Rosetta backrub ensemble) | **No — PyRosetta/Rosetta** | [J.Phys.Chem.B 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5980710) | Pearson ~0.65 |
| **Cartesian ddG** | Carryover (monomer) | Physics ΔΔG | **No — Rosetta** | [Rosetta docs](https://docs.rosettacommons.org/docs/latest/cartesian-ddG) | stability baseline |
| **mCSM-PPI2 / DDMut-PPI** | Carryover | ML/graph binding ΔΔG | Webserver only | [NAR 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602427/) | CAPRI blind 1st |
| **Boltz-2 affinity head** | Scoring (co-fold) | ML Kd/IC50 (small molecule); weak on PPI | Yes (MIT) / yes / yes | [biorxiv 2025.06.14.659707](https://www.biorxiv.org/content/10.1101/2025.06.14.659707.full.pdf) | FEP+ 0.6 (SM); PPI 0.15 |
| **ipSAE** | Scoring (pose gate) | Reweighted ipTM/PAE | Yes / yes / yes (PyPI) | [PMC11844409](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11844409/) | better binder separation |
| **EvoProtGrad** | Carryover loop (seq MCMC) | Any plug-in fitness (wrap a ddG model) | Yes / yes / yes (GitHub) | [arXiv:2212.09925, 2023](https://arxiv.org/abs/2212.09925v1) | — (framework) |
| **EvoBind2** | Partial (co-adapts receptor) | AF plDDT/confidence directed evolution | Yes / yes / yes | [biorxiv 2022.07.23.501214](https://www.biorxiv.org/content/10.1101/2022.07.23.501214.full.pdf) | peptides 0.26 nM–7.9 µM |
| **ALLM-Ab** | Carryover loop (AL + PLM) | Learning-to-rank + Flex ddG oracle | Yes / mostly (oracle uses Flex ddG) | [biorxiv 2025.08.05.668775](https://www.biorxiv.org/content/10.1101/2025.08.05.668775.full.pdf) | DMS-validated, developability-aware |
| **PPIFlow (+ maturation)** | Redesign (backbone flow) + interface maturation | Flow-matching + AF3Score | Code unverified | [biorxiv 2026.01.19.700484, Jan 2026](https://www.biorxiv.org/content/10.64898/2026.01.19.700484v1) | sub-nM VHHs (no Kd table seen) |
| **BindCraft** | Redesign (de-novo) | AF2 backprop (ipTM/confidence) | Yes / yes / yes | [biorxiv 2024.09.30.615802](https://www.biorxiv.org/content/10.1101/2024.09.30.615802.full.pdf) | 10–100% success, nM binders |
