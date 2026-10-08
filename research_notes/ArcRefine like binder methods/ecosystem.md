# Non-standard protein-design publishing & Adaptyv competition ecosystem (2024–2026)

Scope: surface actionable, code-available binder-design methods and competition intelligence that standard arXiv/bioRxiv/PubMed searches miss — the class of release that caused us to miss ArcRefine. Emphasis on items usable for a de-novo TNF-α pH-switch mini-binder and on scoring well in the Adaptyv assay. Research date 2026-10-08.

---

## Q1. What is foldarc.com? Other methods it hosts, who is behind it, what else they released

### Takeaway
FoldArc is a one-person commercial binder-design service run by **Ken Osumi** (FoldArc Co., Ltd.) whose only original published method is **ArcRefine** ("structural carryover") — a Boltz-2 + Mosaic gradient refinement that improves de-novo binders while preserving binding mode, with a ready-to-run **TNF-α example** and **6/10 Adaptyv-BLI-confirmed TNFα binders for ~$677** of compute. Code is on GitHub, data/paper on Zenodo; nothing of this is indexed by arXiv/PubMed, which is exactly why we missed it.

### Cited Findings
- FoldArc is a "custom protein design service that creates binders targeting specific protein interaction sites"; delivers ~100-aa (10–20 kDa) sequences + per-design computational evaluation; standard plan = 10 candidate sequences — [FoldArc home](https://foldarc.com)
- The site's one research artifact is the preprint **"ArcRefine improves de novo binder designs while preserving binding modes,"** author **Ken Osumi**, 2026, unreviewed — [FoldArc home](https://foldarc.com); research page [foldarc.com/research/structural-carryover/](https://foldarc.com/research/structural-carryover/)
- **Method**: gradient-based refinement that uses "structural carryover" — Boltz-2's single and pair representations are passed between sequence updates (built atop the **Mosaic** framework), with ProteinMPNN for soluble sequence generation; validated by an *independent* Boltz-2 prediction without carried state ("fresh_prediction.cif") — [ArcRefine GitHub](https://github.com/ken-osumi/ArcRefine); [Zenodo record](https://zenodo.org/records/23115832)
- **TNFα result**: 6 of 10 ArcRefine-optimized designs bound vs 1 of 10 unoptimized; interface confidence median rose 0.49 → 0.69 across multiple predictors; parent designs came from RFdiffusion/BoltzGen; ~US$677 compute (excl. wet-lab) — [foldarc.com/research/structural-carryover/](https://foldarc.com/research/structural-carryover/); [Zenodo](https://zenodo.org/records/23115832)
- **PD-L1 result**: 4 of 10 optimized bound vs 1 of 10 unoptimized (same protocol) — [foldarc.com/research/structural-carryover/](https://foldarc.com/research/structural-carryover/)
- Honest limitation stated by the author: results "do not isolate the effect of structural carryover alone or establish that binding occurs at the predicted site" — [foldarc.com/research/structural-carryover/](https://foldarc.com/research/structural-carryover/)
- **Experimental validation was done by Adaptyv Bio via BLI** on PD-L1 and TNFα — [FoldArc home](https://foldarc.com)
- **Repo specifics**: Linux + NVIDIA CUDA GPU, `uv`, Python 3.12; Boltz-2 weights auto-download, ProteinMPNN bundled; CLI `structural-carryover design.json --output run` then `structural-carryover-predict`; default "50/50/15 updates" (`--carryover off` to ablate); ships `examples/tnfa/config.json` and `examples/pdl1/config.json`; `scripts/evaluate.py` reports contact retention, structural displacement, and **ipSAE**; Modal cloud path (one H200, 40-min timeout) via `modal_run.py` — [ArcRefine GitHub](https://github.com/ken-osumi/ArcRefine)
- **License**: PolyForm Noncommercial 1.0.0 (commercial use needs permission from ken@foldarc.com; academic/public-research use allowed without separate patent license) — [ArcRefine GitHub](https://github.com/ken-osumi/ArcRefine)
- **Zenodo**: version DOI 10.5281/zenodo.23156179 (V2.3, 2026-10-05); concept DOI 10.5281/zenodo.23115832; files = ArcRefine_manuscript_V19.pdf (2.0 MB) + structural-carryover-data.zip (394.6 KB); CC BY 4.0 — [Zenodo record](https://zenodo.org/records/23115832)
- **Ken Osumi's other GitHub repos are all forks, no second original method**: ArcRefine (9★, original) plus forks of RFdiffusion, ColabDesign (sokrypton), Protenix (ByteDance), a `mosaic-fork`, and Prefect (workflow orchestration) — [github.com/ken-osumi](https://github.com/ken-osumi)

### Inferences
- ArcRefine is directly actionable for us *now*: take our existing de-novo TNFα backbones/sequences (e.g., the RFd3-derived binders in memory), drop them into `examples/tnfa/config.json` as the `parent_sequence`, and run structural-carryover to (a) lift hit rate (demonstrated 6× on TNFα) and (b) generate sequence-diverse variants that keep the same interface — because ~two-thirds of residues can change while the binding arrangement is preserved, it doubles as a **novelty generator** that holds the interface (relevant to the structural/sequence novelty gate).
- His fork list reveals the real FoldArc stack and lineage: RFdiffusion/BoltzGen (generate) → ColabDesign/Mosaic (gradient refine) → Boltz-2/Protenix (fold/validate), orchestrated with Prefect. ArcRefine = Mosaic + Boltz-2 carryover wrapper.
- Treat ArcRefine as a hit-rate/novelty booster, not proof of site-specific binding (author's own caveat); keep our independent af2ig/Boltz validation and (for the pH switch) our explicit-state MD occupancy gate downstream of it.

### Gaps
- No standalone pH-switch functionality in ArcRefine; it preserves a single binding mode, so dual-pH selectivity still needs our own machinery layered on top.
- Zenodo data zip contents not inspected byte-level here (sequences of the 6 TNFα binders may be inside structural-carryover-data.zip — worth downloading to seed partial-diffusion/ArcRefine runs from a *validated* TNFα binder).

---

## Q2. Adaptyv competition ecosystem: winners, baselines, assay, and what scored well on TNF-α / cytokines

### Takeaway
Adaptyv's assay is **cell-free expression → label-free kinetics (BLI for EGFR, SPR for Nipah), Kd from a concentration series, ≥2 replicates** — so **expression is a hard gate** and the single biggest, most-ignored lever (winners optimized solubility; most entrants only optimized binding). The published EGFR post-mortem shows **iPAE/ipTM are weak predictors (ROC≈0.64)** while **ESM log-likelihood, interface ΔG, contact count and H-bond fraction** separate real binders; **SolubleMPNN is essential**; and **diversifying/optimizing an existing binder beats pure de-novo** on success rate. The current Anthropic×Adaptyv 2026 program (5 challenges, Proteinbase, results Dec 15 2026) selects designs via an undisclosed **multi-metric Claude workflow, explicitly "not one in silico metric."**

### Cited Findings

**The published post-mortem (primary source):**
- Paper: *Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition*, bioRxiv 2025.04.17.648362 — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362); mirror [Twist](https://twistbioscience.com/resources/publication/crowdsourced-protein-design-lessons-adaptyv-egfr-binder-competition)
- **Assay = cell-free protein expression followed by Bio-Layer Interferometry (BLI), ≥2 replicates per design, EGFR at multiple concentrations to measure Kd** — [bioRxiv full text](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- Expression/binding yields: R1 73% expressed (146/201), 5 novel binders, hit rate 2.5%; R2 95% expressed (378/400), 53 confirmed binders; 60 binders total — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- **Metric → binding-success ROC AUC: iPAE 0.64, ipTM 0.64, pLDDT 0.55–0.66 (all weak); ESM2 log-likelihood 0.72 (best single; ESM3/ESM-C stronger); interface area pDockQ/mpDockQ ~0.68** — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- **Expression predictors: Glu composition ROC 0.77, Lys composition 0.73, α-helix fraction 0.67** — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- Real binders had favorable binding ΔG (p<0.001), more interface contacts (p<0.001), higher H-bond % (p=0.004); iPAE correlated with affinity only at Pearson −0.25; **the top-2 R1 binders ranked 54th and 89th on the iPAE leaderboard** — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full); [Adaptyv substack "Protein Optimization 102"](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons)
- **SolubleMPNN is decisive for expression: 93.1% (27/29) vs ~0% with standard ProteinMPNN; "use SolubleMPNN as a final check before submitting"** — [Adaptyv substack](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons)

**Who won / what they used:**
- R1 1st = **BindCraft** (AF2 hallucination of the interface; loss = interface confidence + residue count + radius of gyration + helicity), Kd 491 nM; open-source Colab — [Adaptyv substack](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons)
- R1 2nd = Khondamir Rustamov: **AF2Seq hallucination + SolubleMPNN** (ranked 54th by iPAE, but SolubleMPNN gave expression) — [Adaptyv substack](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons)
- R1 3rd = Tripp & Kaltenbrunner: **RFdiffusion → LigandMPNN → ESMFold/ColabFold, 3 refinement cycles, ProtFlow** orchestration — [Adaptyv substack](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons)
- R2 winner = **Cradle**: lead-optimized Cetuximab-scFv (mutated framework, kept CDRs) → 1.21 nM, 8.2× over Cetuximab (9.94 nM); only 2 designs beat Cetuximab — [Cradle blog](https://www.cradle.bio/blog/adaptyv-protein-design-competition); [Cradle update](https://cradle.bio/blog/adaptyv2)
- Submission-volume popularity: de-novo+filter (RFdiffusion→ProteinMPNN→AF2 iPAE filter), diversified-binder+filter, AF2 hallucination; **RFdiffusion reported ~19% hit rate** — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)

**Nipah competition (most recent completed; closest "hard target" precedent):**
- >10,000 designs, >650 participants; **RFdiffusion had the highest hit rate**; validated by **SPR**; top-3 teams chosen on experimental binding (classification + Kd) — [WebSearch of Adaptyv/Proteinbase Nipah results](https://manifold.markets/Proteinbase)
- One Nipah challenge: 1,196 tested → 1,028 expressed → 111 bound (~9.3%) — consistent with the 86% expression / 9.3% binding Nipah figures cited in the competition guide — [intuitionlabs guide](https://intuitionlabs.ai/articles/anthropic-adaptyv-protein-design-competition-guide)

**Current Anthropic × Adaptyv 2026 competition:**
- 5 challenges, Sep 28–Oct 31 2026; validation by Nov 30; **results (incl. failures) published openly on Proteinbase Dec 15 2026**; >$1M wet-lab validation + $1M Claude credits; >5,000 designs tested; 3 tracks (Track 1 expert ≤~20 teams, Track 2 academics/industry, Track 3 open-to-all) — [intuitionlabs guide](https://intuitionlabs.ai/articles/anthropic-adaptyv-protein-design-competition-guide); [Dealroom](https://dealroom.co/news/155261-anthropic-and-adaptyv-bio-launch-protein-design-competition-with-5-000-e/); [Scalevise](https://scalevise.com/resources/claude-protein-design-competition-adaptyv-bio/)
- Selection uses a "jointly agreed Claude-based workflow"; **"selection will not rely on one in silico metric"**; exact workflow disclosed only after the competition — [intuitionlabs guide](https://intuitionlabs.ai/articles/anthropic-adaptyv-protein-design-competition-guide)
- Adaptyv's own binding-characterization skill covers **both SPR and BLI**, and explicitly handles **buffer pH** (for charge-mismatch/nonspecific binding, set buffer pH ±1 from the protein pI; add 1 mg/mL BSA) — confirms the assay can run **binding at two pH values**, which is how a pH-switch would be read — [Adaptyv binding-characterization skill](https://agentskills.so/skills/adaptyvbio-protein-design-skills-binding-characterization)

### Inferences
- **Scoring-well playbook for the Adaptyv assay (synthesized):** (1) pass the expression gate first — run SolubleMPNN and Adaptyv's expression filters (instability index <40, GRAVY <0.4, avoid Glu/Lys-poor, favor α-helix); a beautiful non-expresser scores zero. (2) Rank by a *composite*, not iPAE/ipTM alone — add ESM-C/ESM3 log-likelihood (length-normalized), ipSAE_min, interface ΔG/contacts/H-bonds. (3) Hedge de-novo risk by also submitting *optimized/diversified* variants of any binder that already works (ArcRefine is purpose-built for this). (4) For the pH switch, plan BLI at pH 7.4 and pH 6.0 as the actual readout; the His-salt-bridge direction + our MD occupancy/xtb selectivity gates are what the assay will ultimately reward.
- Because final selection is a Claude multi-metric workflow (not a single threshold), a design that is *defensible across many orthogonal metrics* (fold confidence, interface metrics, expression, plausibility, and a clear mechanistic pH story) is favored over one that maxes a single number.
- Proteinbase is a goldmine the report-writer/team should mine directly: every past Adaptyv design + its experimental outcome (expressed? bound? Kd?) is open — i.e., a labeled dataset of what actually binds cytokine-class targets in *this exact assay*.

### Gaps
- The **TNF-α Challenge 2 exact scoring rubric and pH/assay spec are not in public prose** (challenge pages live inside the Proteinbase app; the selection workflow is undisclosed by design). The project's own brief already encodes the requirements (bind pH7.4 / release pH6.0, mouse cross-react, affinity) — treat those as authoritative over anything inferable from web text. Verify the exact pH values, temperature, and buffer against the live Proteinbase challenge page.
- No TNF-α-specific Adaptyv post-mortem exists yet (this challenge is still open). The EGFR paper + Nipah results are the best proxies.

---

## Q3. Non-arXiv release venues (Zenodo / OSF / protocols.io / HuggingFace / Papers-with-Code), 2025–2026

### Takeaway
The actionable binder-design work of 2025–26 increasingly ships *code-first* on GitHub + Zenodo + personal domains, with the paper as a secondary (or absent) arXiv/bioRxiv artifact — ArcRefine (foldarc.com + Zenodo, no arXiv) is the archetype. The highest-leverage single resource is **Adaptyv's own open skills repo**, which is a curated, Modal-runnable catalog of the current actionable toolset with the organizer's QC thresholds baked in.

### Cited Findings
- **`adaptyvbio/protein-design-skills`** — 24 Claude Code skills, MIT-licensed, Modal-GPU-wrapped (via `biomodals`). Categories: **Design (13)**: boltzgen, bindcraft, rfdiffusion, mosaic, germinal, protenix, proteinmpnn, ligandmpnn, solublempnn, chai, boltz, alphafold, esm; **Eval (2)**: protein-qc, ipsae; **Experimental (2)**: cell-free-expression, binding-characterization; **Orchestration (3)**: binder-design, protein-design-workflow, campaign-manager; **Utils (4)**: pdb, uniprot, foldseek, setup — [GitHub](https://github.com/adaptyvbio/protein-design-skills); [skills.md](https://cdn.jsdelivr.net/gh/adaptyvbio/protein-design-skills@main/docs/skills.md)
- **Adaptyv `protein-qc` thresholds (verbatim, standard / stringent):** ipTM >0.50 / >0.60; PAE_interaction <12 Å / <10 Å; **ipSAE_min >0.61 / >0.70**; pLDDT >0.85 / >0.90; shape complementarity >0.50 / >0.62; interface_dG < −10 / < −15; buried-unsat H-bonds ≤4 / ≤2; scRMSD <2.0 / <1.5 Å. Expression: instability index <40 / <30; GRAVY <0.4 / <0.2; ESM2 PLL >0.0 / >0.2; folding ΔG < −2 / < −4 kcal/mol. Liabilities: even Cys count; flag NG/NS/NT (deamidation), ≥3 consecutive K/R (proteolysis), hydrophobic runs ≥6 (aggregation). **"Individual metrics alone are too weak — use composite scoring"** (30% pLDDT + 20% ipTM + 20% PAE + 15% SC + 15% seq-plausibility). Campaign health: >15% pass = excellent, 10–15% good, 5–10% marginal, <5% poor — [protein-qc SKILL.md](https://cdn.jsdelivr.net/gh/adaptyvbio/protein-design-skills@main/skills/protein-qc/SKILL.md)
- **Adaptyv `protein-design-workflow` pipeline (verbatim):** target prep → backbone (RFdiffusion or BindCraft, 100–500 designs) → ProteinMPNN (8 seq/backbone) → structure validation (AlphaFold2 or Chai) → filter to 50–200 by composite score — [protein-design-workflow SKILL.md](https://cdn.jsdelivr.net/gh/adaptyvbio/protein-design-skills@main/skills/protein-design-workflow/SKILL.md)
- **Adaptyv `ipsae` skill** — ipSAE = interprotein Score from Aligned Errors; interface-focused, "separates true from false predicted complexes more reliably than ipTM"; AF3 ipSAE_min gave **1.4× average-precision gain over the ipAE score RFdiffusion pipelines commonly filter on**; inputs = PAE (AF2/3 JSON or Boltz1 NPZ) + structure + PAE/dist cutoffs 10–15 Å; thresholds ipSAE_min >0.61/>0.70, LIS >0.35/>0.45, pDockQ >0.5/>0.6; works for AF2/AF3/Boltz1 — [ipsae SKILL.md](https://cdn.jsdelivr.net/gh/adaptyvbio/protein-design-skills@main/skills/ipsae/SKILL.md)
- **Mosaic** (`github.com/escalante-bio/mosaic`) — JAX, multi-objective design by continuous relaxation / gradient descent over soft sequence (lineage: ColabDesign → RSO → BindCraft); composes multiple expert models (Boltz-1/2, ESM-2, ProteinMPNN, AF2, Protenix) into one differentiable loss; experimentally validated binders vs PD-L1, IL7Ra, Nipah. This is the engine ArcRefine extends — [Mosaic GitHub](https://github.com/escalante-bio/mosaic)
- **ArcRefine** itself is the Zenodo exemplar of the non-arXiv pattern (GitHub + Zenodo DOI + personal domain, no arXiv/PubMed index) — [Zenodo](https://zenodo.org/records/23115832)

### Inferences
- The Adaptyv skills repo is effectively the organizer telling us the sanctioned toolchain and thresholds; mirroring its QC (especially ipSAE_min >0.61 and the expression filters) aligns our in-silico pre-screen with how designs get judged.
- Mosaic is the right substrate for a *custom multi-objective* TNFα pH-switch loss (bind@pH7.4 term + repulsion@pH6.0 term + SolubleMPNN/expression term in one differentiable objective) — more flexible than ArcRefine's fixed single-mode carryover for encoding the switch directly.

### Gaps
- Did not separately sweep HuggingFace model hub, OSF, or protocols.io; the GitHub+Zenodo+Adaptyv-skills axis already covers the current actionable toolset, but a dedicated HuggingFace sweep (e.g., ESM-C, Boltz-2, Protenix-v2 weights) could add model-availability detail.
- `campaign-manager` and `binder-design` orchestration skills not read in full (likely encode the end-to-end automation Adaptyv expects).

---

## Q4. Community hubs & lab-group drops (X/Twitter, Discord, blogs, GitHub) that announce code-released tools

### Takeaway
The signal lives on GitHub (often via personal/forked accounts), Substack post-mortems, and non-standard aggregators (skills registries, prediction markets) rather than arXiv. Recurring, actionable lab/tooling drops relevant to us: **BindCraft (Pacesa/Nickel), Mosaic (escalante-bio), Germinal/OpenGerminal (de-novo antibody/nanobody), and Adaptyv's skills + Substack.**

### Cited Findings
- **Adaptyv Substack** is the primary post-mortem channel (design strategies, assay lessons, winner interviews) — [adaptyvbio.substack.com](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons); [archive](https://adaptyvbio.substack.com/archive)
- **Cradle** publishes method write-ups of competition wins (Cetuximab scFv lead-optimization) on its blog — [cradle.bio/blog](https://www.cradle.bio/blog/adaptyv-protein-design-competition)
- **Germinal** (`github.com/SantiagoMille/germinal`) — de-novo antibody/nanobody/scFv vs user-specified epitopes; 3 steps: ColabDesign hallucination → **AbMPNN** redesign → cofold; **4–22% experimental success testing only 43–101 designs/antigen**; Sept-2025 preprint; **OpenGerminal** open implementation; runnable on Tamarind — [Germinal GitHub](https://github.com/SantiagoMille/germinal); [Twist pub](https://www.twistbioscience.com/resources/publication/efficient-generation-epitope-targeted-de-novo-antibodies-germinal); [Tamarind OpenGerminal](https://www.tamarind.bio/tools/opengerminal)
- **Tamarind Bio** hosts browser-runnable versions of these tools (Mosaic hallucination, Germinal, OpenGerminal) — useful if GPU setup is a bottleneck — [tamarind.bio](https://www.tamarind.bio/tools/mosaic-hallucination)
- **Claude-skill registries** (skills.cat, claudskills.com, smithery.ai, skills.sh, agentskills.so) re-index the Adaptyv skills and others — a non-standard discovery surface for method drops — [skills.cat/u/adaptyvbio](https://skills.cat/u/adaptyvbio)
- **Proteinbase** (proteinbase.org; [LinkedIn](https://linkedin.com/company/proteinbase)) — open DB of all Adaptyv competition designs + experimental outcomes; also runs **Manifold prediction markets** on competition outcomes (e.g., Nipah "which method wins") as a crowd-signal of what scores well — [Manifold/Proteinbase](https://manifold.markets/Proteinbase)
- **BindCraft** repeatedly surfaces as the de-novo workhorse for small-budget campaigns (AF2 interface hallucination), open Colab — [Adaptyv substack](https://adaptyvbio.substack.com/p/protein-optimization-102-lessons)

### Inferences
- The detective heuristic going forward: watch GitHub (incl. personal accounts and forks), Zenodo DOIs, Substack post-mortems, Tamarind's tool list, Claude-skill registries, and Proteinbase — not arXiv listings — for code-released binder methods. ArcRefine would have been caught by monitoring Zenodo "binder" deposits and foldarc.com, or by reading Adaptyv's validated-design list (it was validated by Adaptyv's own BLI).
- Proteinbase's open per-design outcomes + the Manifold markets are an underused, competition-specific intelligence source for "what scored well" that no literature search returns.

### Gaps
- Did not directly scrape X/Twitter or Discord (not reliably web-fetchable here); specific handle-level method drops (e.g., @martin_pacesa, ColabDesign/Sergey Ovchinnikov, Boltz/MIT, IPD) are inferred as the originators of the tools above but not quoted from their own posts.
- No evidence found of a public TNF-α-specific community design thread for the current challenge yet (challenge is open; discussion likely in the competition's private channels).

---

## One-paragraph bottom line for the writer
ArcRefine is the exact method class we missed: code-first (GitHub `ken-osumi/ArcRefine`) + Zenodo (10.5281/zenodo.23115832) + personal domain (foldarc.com), no arXiv/PubMed footprint, and it is *directly on-target* — a Boltz-2/Mosaic structural-carryover refiner with a ready TNFα example that turned 1/10 into 6/10 Adaptyv-BLI-confirmed TNFα binders for ~$677, while changing ~2/3 of residues (so it doubles as an interface-preserving novelty generator). Its engine is Mosaic (escalante-bio), which also supports encoding a custom pH-switch multi-objective loss. For scoring well in the Adaptyv assay specifically: expression is the hard gate (SolubleMPNN + Adaptyv's expression filters), iPAE/ipTM are weak (ROC≈0.64) so rank by a composite incl. ESM-C/ESM3 log-likelihood + ipSAE_min (>0.61) + interface ΔG/contacts/H-bonds, mirror Adaptyv's own `protein-qc` thresholds, submit optimized/diversified variants alongside de-novo, and plan the pH switch as a two-pH BLI readout. The richest ongoing intelligence source is Proteinbase's open, per-design experimental outcomes plus Adaptyv's Substack post-mortems and their own `adaptyvbio/protein-design-skills` repo.
