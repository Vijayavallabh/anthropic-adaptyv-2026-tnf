# De-novo Binder Generators & Pipelines (2024–2026) with Code — ArcRefine-like releases we may have missed

Scope: NEW de-novo mini-protein/binder generators and pipelines released 2024–2026 WITH CODE, biased to indie/lab GitHub/HuggingFace/Zenodo drops, evaluated for a flat polar cytokine (human TNF-α β-jellyroll) with pH-release + mouse cross-reactivity, runnable on 8×H100 in ~2 days.

Already-known (per assignment, not re-cataloged except where a new fact appeared): RFdiffusion1/3, BoltzGen, ProteinMPNN, AF2-initial-guess, Boltz-2, partial diffusion, BindCraft, AlphaProteo, Chai-2, EvoBind2, Germinal, ProteinGenerator, La-Proteina, PXDesign, ArcRefine (fold-conditioned AF2-Multimer hallucination, bioRxiv 2025.07.02.662497).

---

## Key Question 1 — NEW generators/pipelines beyond our list (indie/lab GitHub, HF, Zenodo 2025–2026)

### Takeaway
At least four genuinely new generators we did not have are code-released and relevant: **β-pairing targeted RFdiffusion** (Baker, Nat Commun 2025 — runs on our existing RFdiffusion), **Proteina-Complexa** (NVIDIA/ManifoldBio, ICLR 2026, open weights on HF/GitHub/NGC), **SeedProteo** (ByteDance Seed, Dec 2025 — webserver only), and **BoltzDesign1** (MIT-licensed hallucination on Boltz). Plus several pipelines/agents (ProteinDJ, Minibinder Lab) and an Anthropic agentic-campaign benchmark.

### Cited Findings
- **β-pairing targeted RFdiffusion** — conditions RFdiffusion via SS/ADJ (secondary-structure + block-adjacency) tensors so a binder β-strand (length L) forms an extended β-sheet pairing with a chosen target edge β-strand (identity T); a "simple algorithm" maps user (L,T) onto the conditioning tensors and the rest of the binder is left unassigned. Preprint CC-BY 4.0; peer-reviewed as Nature Communications 2025 (s41467-025-67866-3, IPD UW publication #365) — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf); [IPD Nat Commun ref](https://www.ipd.uw.edu/publication-pdfs/365/4c42c439d6300d56d83cbc25f1b2b04d/s41467-025-67866-3_reference.pdf)
- β-pairing uses the SS/ADJ machinery already in RosettaCommons/RFdiffusion (co-author David E. Kim) — runs on a standard RFdiffusion install — [RFdiffusion GitHub](https://github.com/RosettaCommons/RFdiffusion); [davidekim GitHub](https://github.com/davidekim)
- **Proteina-Complexa** (NVIDIA + ManifoldBio) — fully-atomistic flow-matching generator (encoder+decoder+denoiser, pair-biased transformer; extends La-Proteina) combining large-scale generative pretraining on synthetic data (new "Teddymer" dataset) with **test-time compute** (adaptive inference-time search); ICLR 2026 Oral = "Scaling Atomistic Protein Binder Design with Generative Pretraining and Test-Time Compute" — [alphaXiv](https://alphaxiv.org/abs/2603.27950); [NVIDIA GenAIR page](https://research.nvidia.com/labs/genair/proteina-complexa)
- Proteina-Complexa code/weights: GitHub [NVIDIA-Digital-Bio/Proteina-Complexa](https://github.com/NVIDIA-Digital-Bio/Proteina-Complexa) and mirror [NVIDIA-BioNeMo/Proteina-Complexa](https://github.com/NVIDIA-BioNeMo/Proteina-Complexa); HF [nvidia/NV-Proteina-Complexa-Protein-Target-160M-v1](https://huggingface.co/Iambackup/NV-Proteina-Complexa-Protein-Target-160M-v1) (+ Ligand-Target variant); NGC catalog. License = **NVIDIA Open Model License (planned Apache-2.0)**; runs on **NVIDIA A100, ~9–10 s/design** — [HF model card](https://huggingface.co/Iambackup/NV-Proteina-Complexa-Protein-Target-160M-v1)
- Proteina-Complexa evaluated on a **19-target protein benchmark** with "markedly higher in-silico success rates" than RFdiffusion and test-time search "greatly outperforms" BindCraft hallucination; novelty 0.71–0.80; enzyme AME 41/41; ManifoldBio×NVIDIA campaign tested ~1M designs vs 127 targets (~100M PPIs) — [alphaXiv](https://alphaxiv.org/abs/2603.27950); [ICML Minibinder-Lab listing](https://icml.cc/virtual/2026/67892)
- **SeedProteo** (ByteDance Seed) — diffusion-based de-novo **all-atom** design that repurposes a folding architecture into a generative framework via self-conditioning; inputs = target mmCIF + hotspot residues (optional SS conditioning); arXiv Dec 30 2025 — [arXiv 2512.24192](https://arxiv.org/pdf/2512.24192)
- SeedProteo wet-lab hit rates **70–80%** on PD-L1 (picomolar) and SARS-CoV-2 RBD (sub-nM); SeedProteo-D best novelty in 8/10 targets; >60% unconditional success at 1000 residues. **Code/weights: webserver only (https://seedfold.io/proteinDesign); no open GitHub/HF weights found** (searches only surfaced the unrelated "Seed-Coder" LLM) — [arXiv 2512.24192](https://arxiv.org/html/2512.24192v2)
- **BoltzDesign1** (yehlincho) — inverts the Boltz all-atom structure predictor by gradient ascent on logits (pre/soft/temp/hard stages) = hallucination analogous to BindCraft but on Boltz, then ProteinMPNN redesign; **MIT license**; inputs = target PDB/id + target type (protein/SM/DNA/RNA) + chains + GPU id; **currently single-chain only (multi-chain WIP), NOT experimentally validated** — [github.com/yehlincho/BoltzDesign1](https://github.com/yehlincho/BoltzDesign1)
- **ProteinDJ** (PapenfussLab) — freely-available GitHub pipeline (Nextflow-style) for de-novo binder design orchestrating RFdiffusion-family generation + ProteinMPNN + folding validation for HPC scale-out — [github.com/PapenfussLab/proteindj](https://github.com/PapenfussLab/proteindj)
- **Minibinder Lab** — agentic framework of 12 specialized agents coordinating generation/analysis/critique/refinement (LLM reasoning + structural/sequence tools); ICML 2026 "The Reliability Gap of Agents for Designing High-Quality Protein Binders" (a cautionary reliability study) — [ICML 2026](https://icml.cc/virtual/2026/67892)
- **"Autonomous Binder Design with Claude"** (Anthropic, ~Aug 2026) — agentic de-novo binder campaigns run 24–48 h across 16 targets; 1,440 designs + binding data + prompts released as a reproducible protocol/benchmark — [Anthropic PDF](https://www-cdn.anthropic.com/30bf50e22a01388bb29bf077ee3f244531594b7a.pdf)
- **PXDesign** (known) new fact: open at [github.com/bytedance/PXDesign](https://github.com/bytedance/PXDesign) with HF mirror [OneScience-Group/PXDesign](https://huggingface.co/OneScience-Group/PXDesign); 20–73% nM hit across 5/6 targets, diffusion (PXDesign-d) + hallucination (PXDesign-h) — [ByteDance Seed](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders)
- **BoltzGen** (known) new fact: open MIT (weights+data+train+inference), 66% nanomolar hit on novel targets — [github.com/HannesStark/boltzgen](https://github.com/HannesStark/boltzgen); [MIT page](https://jclinic.mit.edu/boltzgen/)

### Inferences
- Proteina-Complexa is the best NEW locally-runnable generator for us: open weights, H100-compatible, ~10 s/design → a multi-thousand-design TNF-α campaign fits inside 2 days on 8×H100, and it reportedly beats RFd/BindCraft in-silico on a broad protein benchmark.
- SeedProteo is the strongest reported newcomer numerically but **webserver-only** makes it a poor fit for our local 8×H100, no-outbound-except-https, high-throughput, deadline-bound setting.
- BoltzDesign1 and ProteinDJ are low-cost add-ons (reuse Boltz/RFdiffusion we already run) but add little over our existing stack; BoltzDesign1's single-chain limitation and lack of wet validation are concerns.

### Gaps
- Exact per-design GPU/runtime for SeedProteo and β-pairing RFdiffusion not published.
- Whether SeedProteo will release weights (ByteDance released PXDesign openly, so a later drop is plausible) — unconfirmed.

---

## Key Question 2 — Specialized for HARD/flat cytokine targets, epitope-focused, β-strand/edge-strand, or minibinder scaffolds

### Takeaway
The standout is **β-pairing targeted RFdiffusion**, purpose-built for polar targets with exposed **edge β-strands** (incl. immunoglobulin β-sandwich folds — structurally the same family as TNF-α's β-jellyroll) where standard RFdiffusion fails; it lifted in-silico hit rate ~9× and gave 137 pM–mid-nM binders. A training-free epitope-focused **complementarity dock-and-design** method and a **peptide** diffusion binder (DiffPepBuilder) round out the specialized options.

### Cited Findings
- β-pairing RFdiffusion explicitly targets "therapeutically relevant polar targets" with "beta-sheets with unpaired and exposed beta-strands"; design set spans single exposed edge strands (α-CTX, FCRL5, NRP1) and **immunoglobulin-fold β-sandwich domains (FCRL5, KIT, PDGFRα)** — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf)
- β-pairing hit-rate lift: **>9.2% of β-strand-interface-conditioned designs met quality metrics vs 0.98% for hotspot-conditioned RFdiffusion** (~9×); 5-fold improved in-silico success on exposed-edge targets; β-pairing yielded in-silico successes for **all 7** targets whereas hotspot conditioning failed entirely for ALK-3 and FCRL5 — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf)
- β-pairing measured affinities (SPR): ALK-2bp 164 nM, **ALK-3bp 528 pM**, α-CTXbp 1.9 nM, KITbp 64.7 nM, FCRL5bp 82 nM, NRP1bp 84 nM, **PDGFRαbp-7LBF 137 pM**, PDGFRαbp-AF2model 4 nM (overall 137 pM–mid nM) — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf)
- β-pairing screening scale (yeast display → SPR): KIT best from 1230 designs; FCRL5 best from 4846; PDGFRα from 3780 and 391-member libraries — modest libraries, wet-screened — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf)
- β-pairing specificity: all-by-all SPR showed **no off-target binding even among homologous family members** (ALK-2/ALK-3, KIT/PDGFRα); KITbp:KIT co-crystal (**PDB 9H71**, 2.8 Å) is 1.9 Å all-atom RMSD to the design model; binders express in E. coli, monomeric, thermostable (TM 63 to >95 °C) despite exposed edge strands — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf)
- **Complementarity-based dock-and-design** (Advanced Science 2025) — training-free: retrieves complementary scaffolds from a structure DB to a query epitope via a fast surface fingerprint, then carves a binding site by mutation; **24 candidates → nanomolar binders across 5 folds** vs 3 epitopes on VEGF and **IL-7Rα** (a cytokine receptor); anti-VEGF tumor-inhibiting in vivo; PDBs 8BL5/8BL9 — [PMC12412469](https://pmc.ncbi.nlm.nih.gov/articles/PMC12412469/)
- **DiffPepBuilder** (2024) — diffusion model for target-specific de-novo **peptide** binder design (SS-bonded scaffolds) — [arXiv 2405.00128](https://arxiv.org/pdf/2405.00128)
- **mBER** (bioRxiv 2025.09.26) — controllable de-novo **antibody** design with million-scale experimental screening (antibody modality, platform-tied) — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.09.26.678877.full.pdf)
- Cytokine-binder precedent: Rosetta designed compact (~60 aa) binders for γc/βc/IFNAR1/TSLPR at low-nM; RFdiffusion used where higher success needed (IFNAR2, IL10Rβ); "Novokines" fuse designed binder pairs into cytokine mimics (bioRxiv 2025.10.12) — [PMC12632978](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12632978/)

### Inferences
- β-pairing targeted RFdiffusion is the most mechanistically-matched tool for TNF-α: a β-jellyroll presents exposed edge β-strands, and β-augmentation (strand-pairing) is a **different binding mode that sidesteps the flat-polar-face wall** that has floored our af2ig PAE ~22. Because it rides on the RFdiffusion we already run (SS/ADJ conditioning), it is implementable within the deadline by scripting the (L,T)→tensor mapping.
- The training-free dock-and-design approach is an orthogonal epitope-focused hedge (works on IL-7Rα, a cytokine receptor) but code availability is unconfirmed, so it is lower-priority for a 2-day build.

### Gaps
- No public code repo pinpointed for the complementarity dock-and-design method (Advanced Science); need to check its supplement.
- β-pairing paper's exact conditioning scripts (the make-SS/ADJ helper for strand pairing) not located as a standalone repo in search; likely in the Nat Commun supplement / RFdiffusion helper_scripts / a Zenodo deposit — verify before relying on turnkey scripts.

---

## Key Question 3 — TNF-α-specific or TNFR-superfamily-specific design releases

### Takeaway
A 2025 paper designed **universal minibinders against human + canine TNF-α at sub-nanomolar affinity** and — critically — did it by targeting a **conserved hydrophobic "Region I" cleft** (not the flat polar face), yielding natural human/canine cross-reactivity AND functional activity in **mouse** cells. Separately, the Baker pH-sensitive binder paper (our known blueprint) includes **TNFR2** (TNFR-superfamily). Together these reframe TNF-α as bindable at a specific conserved epitope and give us exact hotspots.

### Cited Findings
- **TNF-α universal minibinder** (PMC12644665, 2025): designed with **Rosetta RIF + PatchDock + Rosetta design (NOT deep learning)**; scaffold = **triple-helix bundle, 56–58 aa** (two binding helices + one stabilizing helix) — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665/)
- Its target epitope = conserved hydrophobic **Region I** at the ligand–receptor interface: **Val74, Leu75, Ile97 (one monomer) + Trp114, Tyr115 (adjacent monomer)**; the design deliberately **avoids the non-conserved D-E loop**, which is what delivers cross-species conservation — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665/)
- Affinity/potency: best design **TMB-1 = 0.55 nM to canine TNFα**, nanomolar to human; cell IC50 **1.87–5.42 nM (canine), 4.92–6.0 nM (human)**; **4/5 designs expressed**; cross-reactive human+canine and **inhibited TNFα-driven osteoclastogenesis in mouse BMDMs** (functional mouse cross-reactivity); cryo-EM **PDB 9JEC / EMD-61412**; no code repo (Rosetta protocol + supplementary data only) — [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665/)
- **pH-sensitive binders** (Baker, bioRxiv 2025.09.29.678932 — our known blueprint): two methods — (i) His adjacent to positively-charged residues at the interface (low-pH electrostatic repulsion) and (ii) buried His-containing charged H-bond networks (acid destabilization); designed pH-dissociating binders against **TNFR2 (TNFR-superfamily)**, IL-6, EphA2, PCSK9, and the IL-2 mimic Neo2, built with RFdiffusion + ProteinMPNN — [PMC12621740](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12621740/); [preLights](https://prelights.biologists.com/highlights/computational-design-of-ph-sensitive-binders/)
- Newest all-atom generators still fail on TNFα in-silico: SeedProteo benchmark Table 3 binder-design successes = **TNFα: SeedProteo-R 1 / BoltzGen 0 / RFdiffusion3 0** (vs PD-L1 380/17/1 and SC2RBD 92/0/0) — TNFα is the hardest of the 10 benchmark targets and only SeedProteo scored a single in-silico success — [arXiv 2512.24192](https://arxiv.org/html/2512.24192v2)

### Inferences
- The TNFα-minibinder result is the single most actionable intel in this report even though it isn't a new generator: it proves TNF-α is de-novo bindable at **sub-nM with a plain 3-helix bundle** when you target the conserved hydrophobic **Region I** inter-monomer cleft rather than the flat polar β-sheet face, and that this epitope choice **automatically buys human/canine/mouse cross-reactivity** (directly serving our mouse-cross-react requirement). We should re-point every generator (RFd1/RFd3/BoltzGen/Proteina-Complexa/β-pairing) at Region I hotspots (V74/L75/I97 + W114/Y115 across the dimer interface).
- SeedProteo Table 3 independently corroborates our "TNF-α is a wall" lesson across the 2026 SOTA, and that an all-atom generator with hotspot+mmCIF conditioning is the regime where the one success appeared — consistent with our RFd3 breakthrough. This argues for (a) all-atom generators and (b) a precise, bindable epitope (Region I), not the flat face.
- pH switch + TNFR-superfamily is already de-risked by the Baker paper (TNFR2 in-scope), so the switch chemistry transfers; our novelty is applying it on a TNF-α (ligand) Region I binder.

### Gaps
- No open code for the TNFα-minibinder (Rosetta RIF protocol only) — but the hotspot definition and the 9JEC structure are directly reusable as a target spec for our DL generators.
- Whether anyone has publicly released a de-novo binder specifically against the TNF-α **trimer flat face** (as opposed to Region I or the receptor side) — no such code release found; the literature consensus is that the flat polar face is the hard/failing site.

---

## Bottom line — TOP 3 for our situation (flat polar TNF-α, pH switch, mouse cross-react, ~3 days, 8×H100)

1. **β-pairing targeted RFdiffusion** (Baker, Nat Commun 2025; CC-BY) — purpose-built for polar targets with exposed edge β-strands incl. Ig β-sandwich folds (≈ TNF-α β-jellyroll); ~9× in-silico hit-rate lift, 137 pM–mid-nM binders, proven specificity + co-crystal. **Runs on the RFdiffusion we already have** via SS/ADJ strand-pairing conditioning. A fundamentally different binding mode that bypasses the flat-polar-face wall. HIGHEST priority. [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.10.11.617496.full.pdf)
2. **TNFα Region I epitope intel** (PMC12644665) — not a generator, but the validated, cross-species-conserved **hydrophobic hotspots (V74/L75/I97 + W114/Y115)** + a 0.55 nM 3-helix scaffold + mouse functional activity. Redirect all our generators here; it directly serves mouse cross-reactivity. [PMC12644665](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644665/)
3. **Proteina-Complexa** (NVIDIA/ManifoldBio, ICLR 2026; open weights, NVIDIA-Open→Apache) — newest open all-atom generator with test-time compute, beats RFd/BindCraft in-silico, **~10 s/design on A100/H100** → large Region-I TNF-α campaign fits in <2 days locally; complements our RFd3. [GitHub](https://github.com/NVIDIA-Digital-Bio/Proteina-Complexa) · [HF](https://huggingface.co/Iambackup/NV-Proteina-Complexa-Protein-Target-160M-v1)

Runner-ups: **SeedProteo** (only method with an in-silico TNFα success, but webserver-only → poor local fit); **BoltzDesign1** (MIT, reuses our Boltz, cheap hallucination add-on, single-chain/unvalidated caveats).
