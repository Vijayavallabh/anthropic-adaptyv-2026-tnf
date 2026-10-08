# Newer design ecosystem (2025–2026) & cross-species / multi-state binder design we may have missed

Scope: (A) new design frameworks / model hubs / pipelines that bundle generation+refinement+scoring, and (B) methods to design/score a SINGLE binder cross-reactive across human+mouse TNF-α (orthologs ~79–80% identical). Complements the sibling note `../ArcRefine like binder methods/ecosystem.md` (FoldArc/ArcRefine, Mosaic, BindCraft, Germinal, Adaptyv skills, Proteinbase) — those are NOT re-catalogued here; this note is the MISSED/NEWER layer. Research date 2026-10-08. Many entries are 2026 preprints (unreviewed); dates/IDs from live search.

---

## KQ1. New design platforms / frameworks / model hubs (2025–2026) — is there a "FoldArc", what do they bundle?

### Takeaway
"FoldArc" is already fully covered in the sibling note (it is Ken Osumi's commercial service; its one method is **ArcRefine**, a Boltz-2/Mosaic structural-carryover refiner with a ready TNF-α example). The biggest *missed* ecosystem item is **Anthropic's own `uplifting-biomolecular-modeling`** — 36 Apache-2.0 inference-optimization kits that effectively enumerate the sanctioned 2025–26 SOTA toolchain (AF3, Boltz-2, Chai-1, RFdiffusion3, BoltzGen, ProteinMPNN, Caliby, pxdesign, Proteina-Complexa, Genie3, ESM-C, E1…) with uv-installable, GPU-optimized wrappers. Other new bundles: NVIDIA **BioNeMo** (binder-design Blueprint), the **Boltz** design suite (BoltzGen, BoltzDesign1, Boltz-2 affinity), **Chai-2** (proprietary), and gradient "invert-a-predictor" tools **TorchCraft** (AF3) and **BoltzDesign** (Boltz-1).

### Cited Findings

**Anthropic `uplifting-biomolecular-modeling` (the key missed hub):**
- Repo `anthropics/uplifting-biomolecular-modeling`, **Apache-2.0** (kit code; upstreams keep own licenses), reference release (~Sept 2026) accompanying the blog "How Claude is uplifting biomolecular modeling"; **explicitly NOT maintained** — [GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- **36 kits**, each wraps one pinned upstream and adds drop-in speed/memory modes: `off` (stock), `exact` (bit-identical, ~2–10×), `fast` (small numeric diffs within seed variation, ~3–15×, default), `big` (lowest peak GPU mem; `--n_gpu P` splits one prediction across GPUs) — [GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- Binder-design kits: **af2ig**, **colabdesign** (BindCraft), **boltzgen** (BoltzGen 0.3.2), **mosaic** (Boltz-2 via joltz), **pxdesign** (Protenix 0.5.0+pxd), **complexa** (Proteina-Complexa 1.1.0, all-atom binder gen), **ef2inv** (ESM-cookbook gradient binder design on ESMFold2), **genie3** (backbone diffusion + binder design) — [GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- Cofolding kits: **af3_jax**, **af3_torch** (xfold), **boltz2** (2.2.1), **chai1** (0.6.1), **colabfold**, **protenix_v1/v2** (2.0.0), **openfold3** / **openfold3_ob0** (OpenBind-0 weights), **opendde**, **atlasfold**, **rosettafold3**, **esmfold2** — [GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- Sequence-design / PLM / scoring kits: **proteinmpnn**, **caliby** (Caliby/SolubleCaliby), **esm_if1**, **esmc** (ESM C 300m/600m/6b), **e1** (Profluent E1 mutant scoring), **progen2**; backbone: **rfdiffusion1**, **rfdiffusion3** (all-atom) — [GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- Install is per-kit (Docker / Apptainer / **uv venv**); `bash run.sh install` then `check` then `<cmd> --mode fast`; weights SHA-256 verified — matches our uv/https-only constraint — [GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)

**Boltz design suite (open, MIT):**
- **Boltz-2** jointly predicts complex structure + binding **affinity**, "first DL model to approach FEP accuracy while ~1000× faster"; **all code+weights MIT (academic+commercial)** — [boltz.bio](https://boltz.bio); DOI 10.1101/2024.11.19.624167 (Boltz-1) / Boltz-2 tech report 2025
- **BoltzGen** — "Toward Universal Binder Design", MIT J-Clinic; general binder generator (we already run it; kit `boltzgen` 0.3.2) — [jclinic.mit.edu/boltzgen](https://jclinic.mit.edu/boltzgen/)
- **BoltzDesign1** — inverts Boltz-1 (open AF3 reproduction) to design PPIs/complexes without finetuning; Colab + conda setup; **not yet wet-lab validated** — [github.com/yehlincho/BoltzDesign1](https://github.com/yehlincho/BoltzDesign1)

**NVIDIA BioNeMo:**
- Framework + **NIM microservices** + **Blueprints**; the **"generative protein binder design" Blueprint** bundles RFdiffusion + ProteinMPNN + AlphaFold2 (and DiffDock 2.0 for docking) as a GPU-accelerated reference workflow — [NVIDIA dev blog](https://developer.nvidia.com/blog/?p=94851); [NVIDIA blog](https://blogs.nvidia.com/blog/bionemo-large-language-models-drug-discovery/)
- BioNeMo was **open-sourced** by NVIDIA to broaden access — [Future Medicine](https://www.futuremedicine.com/articles/nvidia-open-sourced-bionemo-to-advance-computational-biology-for-global-biopharma-and-science)

**Chai-2 (proprietary, strong):**
- "Zero-shot antibody design in a 24-well plate" — **~16% de novo antibody hit rate** (double-digit, 16–20%), >100× over prior compute methods; 52 novel antigens, **≤20 candidates/target**, scFv + **VHH**, bound ~50% of targets; 2-week design→wet cycle — DOI **10.1101/2025.07.05.663018** (2025-07-05) — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.07.05.663018.full.pdf); [Chai](https://www.chaidiscovery.com/news/introducing-chai-2)
- Availability: **not open** — API/partnership only (no weights/code) — [Chai](https://www.chaidiscovery.com/news/introducing-chai-2)

**TorchCraft (new unified gradient designer on AF3):**
- "TorchCraft: Unified binder design by inverting an all-atom structure predictor" — optimizes sequence logits through **frozen AlphaFold3** weights (confidence+contact+geometric+seq-prior objectives); produced **minibinders and VHHs with experimentally measured binding**, plus cyclic-peptide/ligand benchmarks — arXiv **2609.19770** (2026-09-17) — [arXiv](https://arxiv.org/abs/2609.19770)

### Inferences
- `uplifting-biomolecular-modeling` is the single highest-leverage ecosystem find: it is both a *map* of the current SOTA (if Anthropic optimized a kit for it, it is in-scope for this competition) and a *speed multiplier* for our existing H100 pipeline (fast/big modes), installable the uv way we already use. Priority adds from its list we are likely NOT running: **Caliby/SolubleCaliby** (expression), **ef2inv** and **pxdesign/complexa/genie3** (orthogonal binder generators), **E1/ESM-C** (sequence-plausibility scorers Adaptyv's post-mortem ranked above iPAE).
- For the pH-switch multi-objective, **Mosaic** (already identified) remains the right differentiable substrate; TorchCraft/BoltzDesign are the AF3/Boltz analogues if we want a second, orthogonal gradient generator.

### Gaps
- Could not confirm Proteina-Complexa's standalone paper/DOI or BioNeMo's exact OSS license text (Apache vs other) from primary pages; both are evidenced via the uplifting-kit list and NVIDIA blogs only.
- Did not verify whether BioNeMo's binder Blueprint is runnable without NVIDIA's hosted NIM stack (may not be pure uv/https).

---

## KQ2. Multi-target / multi-state design — design ONE sequence satisfying human AND mouse TNF

### Takeaway
There is now a dedicated literature for exactly our problem. **Chamaileon** unifies multi-state (MS, conformations) and **multi-target (MT, distinct proteins)** design and optimizes ONE sequence across contexts — the clearest on-target method. **Odin-Multi** optimizes one binder against multiple complexes while *discouraging* off-targets. **DynamicMPNN** is a trained multi-state inverse-folding model (conformations, not yet orthologs). The oldest/simplest option, **ProteinMPNN multi-state (MSD)** = averaging per-position logits across states, works but has low experimental success. The pragmatic route available today is an **averaged/worst-case AF2/Boltz loss over the two ortholog complexes** inside Mosaic/BindCraft/ColabDesign.

### Cited Findings
- **Chamaileon** — "Cross-Context Binder Design with Contextualized Modeling and Mixed Sampling": formulates design as **cross-context binding-landscape modeling**, unifying **multi-state AND multi-target** so one sequence binds multiple distinct targets/states; mechanisms = **I3CD** (In-Context Complex Co-Design; decoupled noise schedules for binder seq vs structure to keep sequence consistent across contexts) + **MoPS** (Mixture-of-Paths Sampling at inference, optimizes a single sequence across contexts despite scarce multi-conformation data); evaluated on a new **CROSS** benchmark; code "on GitHub" (URL not retrievable); arXiv **2607.23518**, DOI 10.48550/arXiv.2607.23518 (2026-07-26, rev 2026-09-13) — [arXiv abs](https://arxiv.org/abs/2607.23518); [PDF](https://arxiv.org/pdf/2607.23518)
- MS vs MT definitions (useful framing): "multi-state design maintains interface compatibility across the target's conformational landscape; multi-target design is multi-specificity, where a single sequence is optimized to engage **similar epitopes on distinct proteins** with high affinity" — [arXiv 2607.23518 via search](https://arxiv.org/pdf/2607.23518); [Model-Centric Review of DL for Protein Design, arXiv 2502.19173](https://arxiv.org/pdf/2502.19173)
- **Odin-Multi** — "optimises one binder sequence against **multiple complexes simultaneously**, encouraging desired interactions while **discouraging specified off-target** interactions" (exactly a human+mouse-on / off-target-off objective) — surfaced via search + the DL-design review; primary paper/URL not pinned — [Model-Centric Review, arXiv 2502.19173](https://arxiv.org/pdf/2502.19173)
- **DynamicMPNN** — "Multi-state Protein Design with DynamicMPNN": **first explicit multi-state inverse-folding model**, learns joint p(Y|X₁…Xₘ) over conformations (independent encoders → Deep-Set pooling → autoregressive decode; +DSS cross-attention variant); trained on **46,033 conformational pairs** (75% of CATH); beats ProteinMPNN-MSD by **+25% decoy-normalized RMSD** (4.30 vs 5.84 Å) and **+12% sequence recovery** (42.7 vs 38.0%); authors Abrudan, Joshi, Greenig, Khmelinskaia, Meiler, Vendruscolo, Knowles; CC BY 4.0; arXiv **2507.21938v2** (2025-10-15); ICML 2026 — [arXiv HTML](https://arxiv.org/html/2507.21938v2)
- **ProteinMPNN-MSD** (baseline multi-state): inference-time **averaging of per-position logits** from independent single-state ProteinMPNN runs; historically **low experimental success** (e.g., only 46 2-state hinge designs expressed), reflecting weak multi-state datasets/benchmarks — [arXiv 2507.21938v2](https://arxiv.org/html/2507.21938v2)
- **MotifCraft** — AlphaFold-Multimer **hallucination** workflow for motif scaffolding with binder-design applications (ICML 2026) — [ICML 2026](https://icml.cc/virtual/2026/70782)
- **HalluDesign** — iterative structure-hallucination + sequence design for optimization and de novo design; bioRxiv **2025.11.08.686881** (2025-11-08) — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.11.08.686881.full.pdf)
- **BindCraft** composite loss is already a single differentiable AF2-Multimer objective (interface pLDDT/PAE + contacts + helicity/compactness), so adding a second target term is natural — [BindCraft, DOI 10.1101/2024.09.30.615802](https://www.biorxiv.org/content/10.1101/2024.09.30.615802)

### Inferences
- Ranked options for human↔mouse TNF, by readiness: (1) **cheapest today** — run our existing Mosaic/BindCraft gradient design with a **two-complex loss** = `L(human) + L(mouse)` or `max(L(human), L(mouse))` (worst-case) over the conserved-epitope footprint; this is the "averaged-gradient / min-over-targets hallucination" the user asked about, and Mosaic composes arbitrary model terms so it is directly expressible. (2) **Purpose-built** — try **Chamaileon** (MT mode) and **Odin-Multi** (on/off-target) once their repos are located; they are literally designed for "one sequence, two targets." (3) **DynamicMPNN** is for conformations, not orthologs, so treat ortholog use as experimental — but it is the right tool if we also want the *same* binder to tolerate TNF's trimer-interface flexibility.
- Worst-case (`max`/`min`-over-orthologs) objectives are safer than averaged ones for cross-reactivity: averaging lets a strong human term mask a failing mouse term (the classic multi-state failure mode ProteinMPNN-MSD shows).

### Gaps
- **Could NOT verify** the specific "multitarget AF2 optimization → cross-reactive **KIT (mouse+human)** and **HA stalk (H1/H3/H5)**" claim: it appeared only in AI-generated WebSearch summaries and did **not** appear in any primary source I fetched (AlphaProteo arXiv 2409.08022, BindCraft 2024.09.30.615802, and N. Bennett's UW thesis "Deep Learning Tools for Protein Binder Design" 2023 all confirmed to NOT contain it). Likely a benchmark case of Chamaileon/Odin-Multi mis-attributed by the summarizer — report-writer should verify against Chamaileon's CROSS benchmark before citing KIT/HA as precedent.
- Chamaileon and Odin-Multi **GitHub URLs, licenses, and experimental validation** not confirmed (Chamaileon says "code on GitHub"; Odin-Multi primary source not located). DynamicMPNN paper gives **no code/weights URL** (gap for installability).

---

## KQ3. Scoring cross-species cross-reactivity more reliably than a single co-folder ipTM

### Takeaway
No single magic metric exists, but three concrete upgrades over "one Boltz-2 co-fold delta": (1) replace ipTM with **ipSAE_min** (interface-focused, documented to separate true/false complexes better); (2) take the **worst-case (min) interface score across orthologs**, not an averaged delta, and fold it into the design objective (Chamaileon/Odin-Multi) rather than only scoring; (3) **ensemble multiple orthogonal co-folders** to cancel single-folder pessimism. The deeper validity point from the antibody literature: **cross-reactivity is governed by epitope conservation**, which confirms our conserved-epitope strategy and says the most reliable "score" is a structural check that the designed contacts land on residues identical between human and mouse TNF.

### Cited Findings
- **ipSAE** (interprotein Score from Aligned Errors) is interface-focused and "separates true from false predicted complexes more reliably than ipTM"; AF3 **ipSAE_min gave ~1.4× average-precision gain** over the ipAE score RFdiffusion pipelines commonly filter on; thresholds ipSAE_min >0.61 / >0.70 — [Adaptyv ipsae SKILL.md](https://cdn.jsdelivr.net/gh/adaptyvbio/protein-design-skills@main/skills/ipsae/SKILL.md) (also in sibling `ecosystem.md`)
- Adaptyv's own post-mortem: **iPAE and ipTM are weak predictors (ROC ≈ 0.64)**; composite scoring (fold confidence + interface + sequence plausibility) is required, "not one in silico metric" — DOI 10.1101/2025.04.17.648362 — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- Cross-reactivity follows conserved epitopes (antibody precedent): anti-IL-6 binders that bind an epitope **conserved between human and cynomolgus but divergent in mouse/rat/dog** cross-react with cyno but not the others — [US Patent 8,198,414](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8198414); highly cross-reactive CXC-chemokine binders recognize **residues conserved across chemokines** — [US Patent 11,485,781](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11485781)
- Multi-target design is itself the strongest "scorer": a model that must satisfy both ortholog contexts in one objective (Chamaileon MoPS; Odin-Multi) internalizes cross-reactivity rather than scoring it post hoc — [arXiv 2607.23518](https://arxiv.org/abs/2607.23518)

### Inferences
- Concrete recipe to replace our single co-fold delta: for each design, co-fold vs **human TNF and mature mouse TNF** with ≥2 independent folders (e.g., Boltz-2 + Chai-1 + AF3/Protenix-v2 from the uplifting kits), compute **ipSAE_min per ortholog**, and rank by `min(ipSAE_min_human, ipSAE_min_mouse)` plus a hard check that the interface footprint contacts human/mouse-identical positions. This directly addresses the "folder-pessimism" lesson in memory ([[lesson-chai-mouse-iptm-is-folder-pessimism]]): a *matched* host-repack delta per folder, then take the consensus/worst-case, instead of trusting one folder's absolute cross-species ipTM.
- Because the deciding factor is epitope conservation, the cheapest high-value "score" is sequence-level: overlay the designed contact map on the human/mouse TNF alignment and reject any design whose energetic hotspots sit on non-conserved positions — this is deterministic and folder-independent.

### Gaps
- No published, experimentally-calibrated metric that *predicts* human→mouse cross-reactivity probability for de novo minibinders (as opposed to antibodies) was found; the min-ipSAE/ensemble recipe above is a synthesis, not a validated benchmark.

---

## KQ4. Developability / expression predictors (2024–2026) beyond GRAVY / instability / ESM

### Takeaway
The headline new tool is **Caliby / SolubleCaliby** — an ensemble-conditioned inverse-folding model that **beats SolubleMPNN**: designs SolubleMPNN called "undesignable" are designable under SolubleCaliby, and it is already packaged in Anthropic's uplifting kit (so uv-runnable). For charge-driven developability (viscosity/aggregation, and nonspecific binding that matters for the two-pH assay) there is **ZetaDial** (dial net charge at inference). Multi-objective inverse folding that *co-optimizes* solubility + evolutionary plausibility exists (**property-driven preference-aligned IF**). PLM likelihood scorers (**ESM-C, E1, ProGen2**) remain the best single sequence-plausibility signal (Adaptyv ranked ESM log-likelihood ROC 0.72, above iPAE/ipTM).

### Cited Findings
- **Caliby** — "Ensemble-conditioned protein sequence design": designs sequences for **non-idealized backbones** out of reach of current models by prioritizing broadly generalizable structural signals over phylogenetic bias; **superior to SolubleMPNN (SMPNN)** on binder quality (high predicted binder pLDDT) — DOI **10.1101/2025.09.30.679633** (2025-09-30); Tianyu Lu et al. (Stanford) — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.09.30.679633.full.pdf); [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12621727/)
- **SolubleCaliby** (soluble-only variant): **Protpardelle-1c binder designs deemed undesignable by SolubleMPNN are designable under SolubleCaliby**, exposing SolubleMPNN filter limitations — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.09.30.679633.full.pdf); hosted on [Tamarind](https://www.tamarind.bio/tools/caliby); **L-Caliby** (ligand-conditioned) at ICLR 2026 — [ICLR 2026](https://iclr.cc/virtual/2026/10014509); packaged as kit `caliby` — [uplifting GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- **ZetaDial** — "dialing net charge of protein binders at inference time for therapeutic developability": a **post-sampling, per-protein secant controller around fixed-backbone ProteinMPNN** (feedback loop: sample→measure charge→correct), lowers charge-target MAE vs static biases; targets viscosity/aggregation; BindCraft sweeps explored; arXiv **2609.05451** (~2026-08/09); Syed & Dinneen — [arXiv](https://arxiv.org/abs/2609.05451)
- **Property-driven Protein Inverse Folding with Multi-Objective Preference Alignment** — fine-tunes IF to **jointly improve solubility + evolutionary plausibility** while matching designability; arXiv **2603.06748** (2026) — [arXiv PDF](https://arxiv.org/pdf/2603.06748)
- **Composition regularizers raise expression**: minibinder campaigns adding amino-acid-composition regularizers (match residue distribution to a reference) showed **higher expression than unregularized rounds** — consistent with Adaptyv's finding that Glu/Lys composition predicts expression (Glu-comp ROC 0.77) — [search synthesis]; [Adaptyv EGFR post-mortem, 10.1101/2025.04.17.648362](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- PLM scorers in the uplifting kit for sequence plausibility/expression proxy: **esmc** (ESM C 300m/600m/6b), **e1** (Profluent E1 150/300/600m mutant scoring), **progen2** — [uplifting GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)

### Inferences
- **Swap/augment SolubleMPNN with SolubleCaliby** as the final expression check: it recovers good binders SolubleMPNN wrongly rejects, which matters because expression is the Adaptyv hard gate and our minibinder pool is small. It is uv-runnable via the `caliby` kit, satisfying our constraints.
- **ZetaDial** is a useful secondary knob: net charge affects both nonspecific binding in BLI (Adaptyv sets buffer pH ±1 from pI) and the His-salt-bridge pH-switch electrostatics; dialing net charge to a developable band without disturbing the engineered His is a reversible, cheap pass.

### Gaps
- Neither ZetaDial nor the property-driven IF paper reported wet-lab expression/solubility validation in the abstracts fetched, and **no code/license/URL** was confirmed for either (installability unknown).
- No dedicated 2024–26 predictor specifically validated for **E. coli vs mammalian de novo minibinder expression** was found beyond composition rules + SolubleMPNN/SolubleCaliby; antibody-specific tools (TAP, AbNatiV) exist but are out-of-domain for a non-antibody miniprotein (not pursued here).

---

## KQ5. Which tools are open + runnable in uv/pip under https-only?

### Takeaway
Most of the high-value items are open and uv/https-friendly; the notable closed ones are Chai-2 and AlphaProteo. Anthropic's uplifting kits document a uv-venv install route explicitly, which is the lowest-friction way to get the whole 2025–26 toolchain running the way we already work.

### Cited Findings
- **uv/pip + https, confirmed open**: Boltz-2/BoltzGen/BoltzDesign1 (**MIT**, pip/conda) — [boltz.bio](https://boltz.bio), [github.com/yehlincho/BoltzDesign1](https://github.com/yehlincho/BoltzDesign1); Mosaic (JAX, pip) — [github.com/escalante-bio/mosaic](https://github.com/escalante-bio/mosaic); ProteinMPNN/SolubleMPNN, ESM-IF1, ESM-C, ProGen2 — all via the **uplifting kits (Apache-2.0, uv venv route)** — [uplifting GitHub](https://github.com/anthropics/uplifting-biomolecular-modeling)
- **Open, heavier / container-leaning**: NVIDIA **BioNeMo** (open-sourced, but binder Blueprint leans on NIM/Docker) — [Future Medicine](https://www.futuremedicine.com/articles/nvidia-open-sourced-bionemo-to-advance-computational-biology-for-global-biopharma-and-science); **Caliby/SolubleCaliby** (kit `caliby`; also Tamarind-hosted) — [Tamarind](https://www.tamarind.bio/tools/caliby)
- **Open but code/weights URL unconfirmed** (installability TBD): **Chamaileon** (arXiv 2607.23518, "code on GitHub"), **DynamicMPNN** (arXiv 2507.21938, no URL in paper), **Odin-Multi** (no primary located), **ZetaDial** (arXiv 2609.05451), **TorchCraft** (arXiv 2609.19770), property-driven IF (arXiv 2603.06748)
- **NOT open**: **Chai-2** (API/partnership only) — [Chai](https://www.chaidiscovery.com/news/introducing-chai-2); **AlphaProteo** (DeepMind, "no ML methods released due to biosecurity/commercial considerations") — arXiv 2409.08022 — [ar5iv](https://ar5iv.labs.arxiv.org/html/2409.08022)

### Inferences
- The fastest path to adopting these with our uv/https setup is to pull the relevant **uplifting kits** (caliby, mosaic, pxdesign, complexa, genie3, ef2inv, esmc, e1, boltz2, protenix_v2, rfdiffusion3) rather than installing each upstream from scratch — they ship pinned, digest-verified, uv-venv-installable stacks and GPU speedups for our 8×H100.

### Gaps
- Exact licenses for Chamaileon, DynamicMPNN, Odin-Multi, ZetaDial, TorchCraft not confirmed (arXiv default is non-exclusive distribution of the *paper*, which says nothing about code license).

---

## One-paragraph bottom line for the writer
The biggest missed ecosystem item is Anthropic's own **`uplifting-biomolecular-modeling`** (Apache-2.0, 36 uv-installable, H100-optimized kits) — it both enumerates the sanctioned 2025–26 SOTA (AF3, Boltz-2, Chai-1, RFdiffusion3, BoltzGen, pxdesign, Proteina-Complexa, Genie3, Caliby, ESM-C, E1, …) and accelerates our pipeline; from it we should actually adopt **SolubleCaliby** (beats SolubleMPNN on expression — the Adaptyv hard gate) and the extra gradient generators. For the human↔mouse cross-reactivity objective, the purpose-built methods are **Chamaileon** (unifies multi-target + multi-state, one sequence across contexts) and **Odin-Multi** (one sequence vs multiple complexes with explicit off-target repulsion); the cheapest route today is a **worst-case (min/max) two-complex AF2/Boltz loss in Mosaic/BindCraft** over our conserved epitope. To *score* cross-reactivity better than one co-folder's ipTM, use **ipSAE_min per ortholog across an ensemble of folders and rank by the worst ortholog**, plus a deterministic check that contacts land on human/mouse-identical residues (the antibody literature confirms cross-reactivity = epitope conservation). Caveats: the widely-quoted "multitarget-AF2 KIT/HA-stalk" precedent could NOT be verified in any primary source (likely an AI-summary artifact — verify against Chamaileon's CROSS benchmark); and several 2026 preprints (Chamaileon, Odin-Multi, DynamicMPNN, ZetaDial, TorchCraft) lack confirmed code URLs/licenses, so installability is unproven.
