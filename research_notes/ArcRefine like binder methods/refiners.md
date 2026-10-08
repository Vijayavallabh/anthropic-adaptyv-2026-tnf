# ArcRefine family: gradient / differentiable-folding sequence-optimization & refinement for de-novo binders (2024–2026)

> Scope: methods that improve/optimize a binder **sequence** by running a structure predictor (Boltz, AF2/3, OpenFold, ESMFold, Protenix, RF) *in the loop* — by gradients (hallucination), by directed evolution, or by MPNN↔fold relax-refine — with emphasis on code-released, Boltz-2/Mosaic/ProteinMPNN-compatible tools usable in ~2 days to improve existing human-TNFα de-novo pH-switch binders.

> **Provenance note on the anchor (ArcRefine):** verified real but preliminary. GitHub repo `ken-osumi/ArcRefine` exists (GitHub API `id=1402815951`, desc "improve de novo protein binder designs through structural carryover. Built on Mosaic and Boltz-2"). Zenodo preprint 10.5281/zenodo.23115832 exists, dated **2026-10-05 (3 days before this research)**, **single author Ken Osumi, FoldArc Co., Ltd.** (a commercial binder-design vendor, foldarc.com), **not peer-reviewed**. All wet-lab/ipSAE claims are self-reported in that preprint. Treat as adoptable code + an unvalidated claim, not established result.

---

## KQ1: Methods that optimize a binder SEQUENCE by gradients through a differentiable structure predictor

### Takeaway
There is a well-populated 2024–2026 "hallucination / differentiable-folding" family. The ones on the user's exact stack are **Mosaic** (differentiable Boltz-1/2/AF2/Protenix/ESMFold, MIT — already installed) and **ArcRefine** (Mosaic + "structural carryover", the anchor). Adjacent code-released siblings: **BindCraft** (AF2 backprop, MIT), **BoltzDesign1** (invert Boltz, no clear license), **PXDesign-h** (Protenix hallucination, Apache-2.0), **ColabDesign/AfDesign** (AF2 gradient framework underneath BindCraft & Germinal), **EvoBind2** (AF2 directed-evolution, non-gradient), **Germinal** (ColabDesign hallucination for Ab/nanobody, Apache-2.0).

### Cited Findings
- **ArcRefine** (anchor): "Refine an existing protein binder design while preserving its fold and binding mode… uses structural carryover: Boltz-2's single and pair representations pass between sequence updates. The resulting sequence is then predicted independently with Boltz-2." Built on Mosaic, "with structural carryover added to its optimization loop." — [ArcRefine README (raw)](https://raw.githubusercontent.com/ken-osumi/ArcRefine/main/README.md)
- ArcRefine mechanism detail: starting-complex prediction initializes the carried state; sequence probabilities start random and are optimized with Mosaic's structural losses + soluble ProteinMPNN; default **50/50/15 updates** (phases), tunable via `--steps`, disable with `--carryover off`; final independent prediction uses 3 recycles / 25 diffusion steps. — [ArcRefine README](https://raw.githubusercontent.com/ken-osumi/ArcRefine/main/README.md)
- ArcRefine preprint self-reported results: modified ~2/3 of residues while retaining binding interface in **45/60 PD-L1 and 23/60 TNFα** designs; median ipSAE **0.49→0.69** across six predictors; "six BLI-confirmed binders against TNFα" at ~$700 campaign cost. — [Zenodo 10.5281/zenodo.23115832](https://zenodo.org/doi/10.5281/zenodo.23115832) (user's brief adds "6/10 BLI, 13 nM" — the 6-of-10 denominator and 13 nM are **not** in the abstract I read; unverified, see Gaps)
- ArcRefine preprint frames the family: "Hallucination methods such as **BindCraft and Mosaic** update sequences using gradients from structure predictors and achieve high success rates, but they build binders from random sequences and cannot refine a given design while keeping its binding mode." — [Zenodo abstract](https://zenodo.org/doi/10.5281/zenodo.23115832)
- **Mosaic** (`escalante-bio/mosaic`): "gradient-based optimization over a continuous, relaxed sequence space"; supports **Boltz-1, Boltz-2, AlphaFold2, OpenFold3, Protenix (mini/tiny/base/v1/v2), ESMFold**; integrates ProteinMPNN for sequence recovery; configurable recycling steps; MIT license. Reported 8/10 (PD-L1) and 7/10 (IL7Ra) experimental validation. — [Mosaic GitHub](https://github.com/escalante-bio/mosaic); license/metadata via GitHub API (MIT, pushed 2026-10-06, 375★)
- **BindCraft** (`martinpacesa/BindCraft`): backpropagation through AF2-Multimer + (Sol)MPNN + PyRosetta; iterative soft/temporary/hard/greedy stages minimizing a multi-term loss (binder pLDDT, interface pAE, contacts, radius of gyration, helicity). MIT, 1239★. — [BindCraft GitHub](https://github.com/martinpacesa/BindCraft); [BindCraft wiki](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft)
- **BoltzDesign1** (`yehlincho/BoltzDesign1`): "Inverting all-atom structure prediction model for generalized biomolecular binder design" (Cho, Pacesa, Zhang, Correia, Ovchinnikov, MIT/EPFL). Inverts Boltz via logits→softmax→temperature-anneal→hard phases; supports **partial redesign with interface residues (<4 Å) fixed**; uses LigandMPNN/ProteinMPNN; AF3 for validation. — [BoltzDesign1 GitHub](https://github.com/yehlincho/BoltzDesign1)
- **PXDesign** (`bytedance/PXDesign`): both diffusion (PXDesign-d) and **hallucination (PXDesign-h)** through Protenix; pipeline = diffusion/hallucination → ProteinMPNN → AF2-IG → Protenix eval; reports 17–82% nanomolar hit rate on 6/7 targets, "surpassing AlphaProteo." Apache-2.0, 266★. — [PXDesign GitHub](https://github.com/bytedance/PXDesign); [ByteDance Seed page](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders)
- **ColabDesign / AfDesign** (`sokrypton/ColabDesign`): AF2-based framework for fixed-backbone design, hallucination, **binder design**, and **partial scaffolding/loop redesign** (can fix positions / use templates); the gradient engine underneath BindCraft and Germinal. 946★. License = **NOASSERTION** (no standard LICENSE file detected by GitHub API). — [ColabDesign GitHub](https://github.com/sokrypton/ColabDesign)
- **EvoBind2** (`patrickbryant1/EvoBind`): AF2/AF-Multimer **in-silico directed evolution** (iterative mutate→fold→score, MCMC-style, *not* gradient) of linear & cyclic peptide binders from target **sequence only**; adversarial false-positive checks; 46% (linear) / 75% (cyclic) success, affinities 0.26 nM–7.9 µM. No license file (GitHub API: NONE), 280★. — [EvoBind GitHub](https://github.com/patrickbryant1/EvoBind); [EvoBind2 preprint](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12280060/)
- **Germinal** (`SantiagoMille/germinal`, Stanford): 3-step de-novo **antibody/nanobody** pipeline = ColabDesign hallucination → AbMPNN selective redesign → cofolding; 4–22% experimental success testing 43–101 designs/antigen, nanomolar affinities. Apache-2.0, 297★. — [Germinal GitHub](https://github.com/SantiagoMille/germinal); [Germinal preprint](https://www.biorxiv.org/content/10.1101/2025.09.19.677421.full.pdf)
- Related 2025–2026 hallucination work (code status less certain): **MotifCraft** (AF2-Multimer hallucination for motif scaffolding, ICML 2026) — [ICML listing](https://icml.cc/virtual/2026/70782); **Proteina-Complexa** (unifies generative + hallucination, fully atomistic binder generation) and **"Scaling Atomistic Protein Binder Design with Generative Pretraining and Test-Time Compute"** (ICLR 2026) — [ICLR poster](https://www.iclr.cc/virtual/2026/poster/10007211); **"Fold-Conditioned De Novo Binder Design via AF2-Multimer Hallucination"** — [bioRxiv 2025.07.02.662497](https://www.biorxiv.org/content/10.1101/2025.07.02.662497.full.pdf)

### Inferences
- ArcRefine is not a new model — it is a **thin fork of Mosaic** ("the `structural-carryover` command names are retained for compatibility"; "includes the required Mosaic modules"). The entire novelty is: (a) seed the optimizer from an existing design + target complex, and (b) **carry Boltz-2's single/pair embeddings across sequence-update steps** instead of recomputing from scratch — a warm-start that biases every update toward the already-found binding mode.
- Because the user already has Mosaic + Boltz-2 + ProteinMPNN, they can either run ArcRefine directly **or** reproduce its core trick inside their own Mosaic install (seed logits from the current TNFα design; reuse/freeze the Boltz-2 trunk representation between steps; keep interface residues biased) — the latter avoids ArcRefine's noncommercial license.
- "Gradient-through-folder sequence optimization" and "MPNN↔fold relax-refine" are two different mechanisms that are often chained: BindCraft/Mosaic/BoltzDesign1/PXDesign-h all end with an MPNN re-sequence + an independent AF2-IG/Protenix/Boltz check, matching the user's existing ProteinMPNN + AF2-initial-guess loop.

### Gaps
- The "6/10" BLI denominator and "13 nM" affinity in the user's brief are not in the Zenodo abstract I could read; they may be in the manuscript body/supplement (not fetched). Flagged as unverified.
- Could not fetch per-step GPU/time numbers for most tools (see KQ3).

---

## KQ2: Embedding-carryover / recycling / binding-mode-preserving refinement tricks

### Takeaway
"Binding-mode preservation" is achieved in three code-released ways: (1) **feed the prior/target structure as a template or AF2 "initial guess"** (dl_binder_design/AF2-IG — the user already has this; BindCraft `predict_initial_guess`/`predict_bigbang`; Boltz/BoltzDesign `template_pdb`); (2) **warm-start the folder's internal state across optimization steps** (ArcRefine's single/pair-embedding carryover; AF2/Boltz recycling); (3) **fix/constrain interface residues during redesign** (BoltzDesign1 <4 Å fix; ColabDesign position masks). ArcRefine's specific contribution is #2 applied to Boltz-2 *between sequence updates*, which is not otherwise packaged as a product.

### Cited Findings
- **ArcRefine structural carryover** = reuse Boltz-2 single & pair representations between sequence updates, initialized from the starting-complex prediction, so the optimizer keeps the existing binding mode; the final sequence is then re-predicted *without* carried state/templates as an honesty check. — [ArcRefine README](https://raw.githubusercontent.com/ken-osumi/ArcRefine/main/README.md)
- **AF2 "initial guess"** (Bennett et al. 2023, Nat Commun; `nrbennet/dl_binder_design`, MIT, 419★): AF2 modified to accept the designed complex coordinates as an initial structural guess; raises design success ~10×. This is the binding-mode-preserving *evaluation* trick the user already runs (AF2-IG). Pipeline = ProteinMPNN → optional PyRosetta FastRelax → AF2-initial-guess. — [dl_binder_design GitHub](https://github.com/nrbennet/dl_binder_design)
- **BindCraft** exposes `predict_initial_guess` (binder atom positions as a bias) and `predict_bigbang` (atom-position bias initialization), recommended for post-MPNN / large-complex cases — i.e., template the binder pose to preserve the mode. — [BindCraft wiki](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft)
- **Mosaic** supports configurable recycling steps and ProteinMPNN sequence recovery inside the loop — the primitives ArcRefine builds carryover on top of. — [Mosaic GitHub](https://github.com/escalante-bio/mosaic)
- **BoltzDesign1 / ArcRefine** both accept `template_pdb` + `template_chain` for target chains (relative paths in the JSON), letting you pin the target (and in ArcRefine the whole starting complex) — [ArcRefine README](https://raw.githubusercontent.com/ken-osumi/ArcRefine/main/README.md); [BoltzDesign1 GitHub](https://github.com/yehlincho/BoltzDesign1)

### Inferences
- The user's existing AF2-initial-guess is the *prediction-side* version of the same idea ArcRefine applies on the *optimization-side*. The cheap, license-clean experiment: inside their current Mosaic run, (a) seed from the existing TNFα binder sequence, (b) raise `recycling_steps`, (c) keep interface residues fixed/biased — approximating carryover without ArcRefine's code.
- For the **pH-switch**: none of these tools model protonation/pH directly, so binding-mode-preserving refinement is complementary to (not a substitute for) the user's existing PROPKA/xTB/MD pH-selectivity gating. Refine the pH-7.4 bound state's confidence/interface with a carryover method, then re-run the pH-selectivity gate on the refined sequence.

### Gaps
- No evidence any released tool carries embeddings across steps *other than* ArcRefine; standard recycling (within a single prediction) is the closest widely-available analog. The claim "only ArcRefine does cross-step carryover" is plausible but I found no independent confirmation.

---

## KQ3: Code availability, license, inputs, GPU/time per design

### Takeaway
All core candidates are code-released. Licenses: **MIT** (Mosaic, BindCraft, dl_binder_design), **Apache-2.0** (PXDesign, Germinal), **PolyForm Noncommercial 1.0.0** (ArcRefine — free for universities/public research, commercial needs permission), and **no standard license detected** (ColabDesign=NOASSERTION, EvoBind=NONE, BoltzDesign1=NONE → treat as all-rights-reserved until the LICENSE file is checked). Per-design GPU/time is rarely published; expect seconds–minutes/step on an H100 for Boltz-2-sized gradient steps, minutes–tens-of-minutes per full design.

### Cited Findings (code · license · inputs · hardware)
- **ArcRefine** — [github.com/ken-osumi/ArcRefine](https://github.com/ken-osumi/ArcRefine). License **PolyForm Noncommercial 1.0.0** (commercial → ken@foldarc.com; universities/public research exempt from separate patent license). Inputs: `design.json` with `parent_sequence` (starting binder) + one entry per target chain (sequence, optional `template_pdb`/`template_chain`, `use_msa`). Install: Linux + NVIDIA CUDA GPU, `uv venv --python 3.12`, `uv pip sync requirements-linux-cuda12.lock`. Boltz-2 weights auto-download; soluble ProteinMPNN weights bundled; bundles its own Mosaic (use a separate env). Optional Modal backend = **1× H200**, 40-min launcher timeout, `--smoke` = 1 update/phase smoke test. Ships `examples/tnfa/config.json`. — [ArcRefine README](https://raw.githubusercontent.com/ken-osumi/ArcRefine/main/README.md)
- **Mosaic** — [github.com/escalante-bio/mosaic](https://github.com/escalante-bio/mosaic). **MIT**. Inputs: Python API; target/binder as sequences, composable loss pytrees; Boltz-1/2, AF2, OpenFold3, Protenix, ESMFold backends; ProteinMPNN integration. Hardware: JAX/CUDA GPU; "first JAX-compiled call is slow," no per-design benchmark published. Actively developed (pushed 2026-10-06). — GitHub API + [repo](https://github.com/escalante-bio/mosaic)
- **BindCraft** — [github.com/martinpacesa/BindCraft](https://github.com/martinpacesa/BindCraft). **MIT**. Inputs: target PDB (trimmed), hotspot residues or null (AF2 picks), binder length range. Hardware: CUDA GPU, **≥32 GB VRAM recommended**; "a few hundred to thousands of trajectories" per target; per-accepted-design time not stated. Colab exists but limited for large complexes. De-novo from scratch (not a refiner), but `predict_initial_guess` biases toward a pose. — [BindCraft wiki](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft)
- **BoltzDesign1** — [github.com/yehlincho/BoltzDesign1](https://github.com/yehlincho/BoltzDesign1). License **NOT detected by GitHub API (NONE)** — a secondary summary claimed "MIT" but that is **unverified; check the LICENSE file before relying on it**. Inputs: FASTA or PDB (or PDB ID), target = protein/small-molecule/NA/metal/PTM; partial redesign fixes interface residues <4 Å. Deps: Boltz weights (auto), LigandMPNN/ProteinMPNN, AF3 (separate) for validation, optional PyRosetta, conda/py3.10. GPU via `--gpu_id`; no time benchmark. "Under active development, NOT experimentally validated." — [repo](https://github.com/yehlincho/BoltzDesign1)
- **PXDesign** — [github.com/bytedance/PXDesign](https://github.com/bytedance/PXDesign). **Apache-2.0**. Inputs: target structure; pipeline diffusion/hallucination → ProteinMPNN → AF2-IG → Protenix. Web server at protenix-server.com; benchmark suite `bytedance/PXDesignBench`. Hardware not specified in results. — [ByteDance Seed](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders)
- **ColabDesign/AfDesign** — [github.com/sokrypton/ColabDesign](https://github.com/sokrypton/ColabDesign). **NOASSERTION** (no standard LICENSE). Inputs: target PDB + binder length (binder mode), position masks/templates for partial design; AF2 weights. Colab-first; GPU required. — GitHub API + repo
- **EvoBind2** — [github.com/patrickbryant1/EvoBind](https://github.com/patrickbryant1/EvoBind). License **NONE**. Inputs: target **sequence only** (optionally receptor interface residues); AF2/AF-Multimer weights. Directed evolution = many fold calls/design; peptide-scale. Also hosted on Tamarind/Neurosnap. — GitHub API + [EvoBind2 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12280060/)
- **Germinal** — [github.com/SantiagoMille/germinal](https://github.com/SantiagoMille/germinal). **Apache-2.0**. Inputs: target structure + user-specified epitope residues; outputs nanobody/scFv. Deps: ColabDesign (AF2), AbMPNN, cofolder. — GitHub API + [preprint](https://www.biorxiv.org/content/10.1101/2025.09.19.677421.full.pdf)
- **dl_binder_design (AF2-initial-guess)** — [github.com/nrbennet/dl_binder_design](https://github.com/nrbennet/dl_binder_design). **MIT**. Inputs: binder–target complex backbone → ProteinMPNN → optional PyRosetta FastRelax → AF2-IG (outputs pLDDT, PAE, interface PAE, RMSD). This is the user's existing AF2-IG relax-refine loop. — [repo](https://github.com/nrbennet/dl_binder_design)

### Inferences
- **2-day runnability verdict (on 8× H100, Boltz-2 v2.2.1 + Mosaic + ProteinMPNN + AF2-IG already present):**
  - ArcRefine — **YES, fast path.** `git clone` + `uv pip sync` + adapt `examples/tnfa/config.json` with the user's current binder + 1TNF target. Bundles Mosaic+ProteinMPNN; Boltz-2 weights auto-download. Main risk = env/lockfile friction and the noncommercial license.
  - Mosaic DIY-carryover — **YES.** Already installed; write a small driver that seeds logits from the current design and raises recycling / fixes interface. MIT, license-clean.
  - BoltzDesign1 partial-redesign — **YES, moderate.** Separate conda env + Boltz weights; interface-fixed redesign is the closest built-in "refine existing binder" mode. License unclear.
  - BindCraft — **YES but it's de-novo**, not a refiner; useful only if generating *fresh* TNFα binders, not improving existing ones. Needs AF2 params + ≥32 GB VRAM (have it).
  - PXDesign-h — **MAYBE.** Needs Protenix env; Apache-2.0 is clean; more setup than the above.
  - EvoBind2 — **NO for this goal** (peptide-only, de-novo, no existing-binder refine).
  - Germinal — **NO** (antibody/nanobody only; user designs mini-protein binders).

### Gaps
- Hard per-design wall-clock/VRAM numbers are unpublished for Mosaic, BoltzDesign1, PXDesign, ArcRefine; only BindCraft (≥32 GB, hundreds–thousands of trajectories) and ArcRefine's Modal config (1× H200, 40-min cap) are concrete.
- BoltzDesign1 / ColabDesign / EvoBind actual license terms need a direct LICENSE-file read before any non-research use.

---

## KQ4: Direct compatibility with Boltz-2 / Mosaic / ProteinMPNN (what the user already has)

### Takeaway
Three tools are natively on the user's stack: **ArcRefine** (= Mosaic + Boltz-2 + ProteinMPNN, by construction), **Mosaic** (Boltz-2 is a first-class backend; ProteinMPNN integrated), and **BoltzDesign1** (Boltz + ProteinMPNN). BindCraft/EvoBind/Germinal are AF2-based (compatible with the user's AF2-IG weights but not Boltz). PXDesign is Protenix-based (Mosaic can also drive Protenix, so there's an indirect path).

### Cited Findings
- ArcRefine: "Built on Mosaic… structural carryover added to its optimization loop," uses **Boltz-2** + **soluble ProteinMPNN** (bundled). — [ArcRefine README](https://raw.githubusercontent.com/ken-osumi/ArcRefine/main/README.md)
- Mosaic: explicit **Boltz-1 and Boltz-2** backends + ProteinMPNN sequence recovery; MIT. — [Mosaic GitHub](https://github.com/escalante-bio/mosaic)
- BoltzDesign1: inverts **Boltz**, uses **ProteinMPNN/LigandMPNN**. — [BoltzDesign1 GitHub](https://github.com/yehlincho/BoltzDesign1)
- PXDesign uses **ProteinMPNN + AF2-IG** and **Protenix** (which Mosaic also supports). — [ByteDance Seed](https://seed.bytedance.com/en/public_papers/pxdesign-fast-modular-and-accurate-de-novo-design-of-protein-binders)

### Inferences
- Lowest-friction, highest-ArcRefine-fidelity action: run **ArcRefine** on the current best TNFα binder(s) using `examples/tnfa/config.json` as a template, with the pH-7.4 bound complex as `template_pdb`; then re-run the user's pH-selectivity gate (PROPKA/xTB/MD) on the refined sequences. Carryover should tighten interface confidence (ipSAE/PAE) without moving the pose.
- License-clean equivalent if noncommercial is a problem: reproduce carryover in the user's MIT Mosaic by warm-starting logits from the existing design + high recycling + interface bias.

### Gaps
- Whether Boltz-2 **v2.2.1** specifically matches ArcRefine's pinned Boltz-2 (via its `requirements-linux-cuda12.lock`) is unverified; a version mismatch could force a separate env (ArcRefine already recommends an isolated env).

---

## TOP 3 most ArcRefine-like + directly adoptable (flagged)

1. **ArcRefine itself** — [github.com/ken-osumi/ArcRefine](https://github.com/ken-osumi/ArcRefine). The exact "refine existing binder, preserve binding mode via Boltz-2 single/pair carryover" method; ships a TNFα example; runs on the user's Mosaic+Boltz-2+ProteinMPNN stack. Caveat: PolyForm **Noncommercial** license (fine for a research competition; not for commercial use) and a 3-day-old single-author unvalidated preprint.
2. **Mosaic with DIY carryover / warm-start** — [github.com/escalante-bio/mosaic](https://github.com/escalante-bio/mosaic). Already installed, **MIT**, Boltz-2 native. Seed from the existing TNFα design, raise recycling, bias/fix interface residues to approximate carryover. Most flexible and license-clean; best if ArcRefine's env or license is a blocker.
3. **BoltzDesign1 partial-redesign** — [github.com/yehlincho/BoltzDesign1](https://github.com/yehlincho/BoltzDesign1). Boltz-based, ProteinMPNN-based, explicit **interface-fixed redesign** mode = refine-existing-binder while holding the contact residues. Runner-up only because the license is undetected and it defaults to Boltz-1. 

Strong runners-up: **BindCraft `predict_initial_guess`** (MIT, pose-biased, but de-novo not refine) and **PXDesign-h** (Apache-2.0, Protenix hallucination + MPNN + AF2-IG).
