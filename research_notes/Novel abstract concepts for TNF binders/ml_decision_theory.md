# Novel ML / Decision-Theory Concepts for Diverse-Novel Generation and Few-Shot Selection of pH-Switchable, Cross-Reactive TNF-alpha Binders

Scope note for the report-writer: this catalog covers abstract ML-theory, representation/geometry, and decision-theory frameworks, each tied to a concrete hook in our pipeline (RFdiffusion backbones -> soluble-ProteinMPNN sequences -> AF2-initial-guess + Boltz-2 + Chai-1 interface scores -> PROPKA/pKAI/MD pH-switch scores -> novelty filter; ~375 pooled wet-lab screen, submit <=20 best-first; strict priority pH-selectivity > mouse cross-reactivity > affinity; structural novelty is a hard gate). Each concept is labeled "feasibility in ~4 days" and "truly-untried vs emerging in protein design." All claims carry inline source URLs; unsourced items are moved to Gaps.

## GFlowNets for diverse, high-reward binder generation

### Takeaway
GFlowNets sample candidates with probability proportional to a (multi-objective) reward instead of collapsing onto a single optimum, which is exactly the right inductive bias for producing a DIVERSE portfolio of novel binders; the method is established for biological sequences (ICML 2022) and is emerging-but-barely-used for antibody/binder design (2025-2026), so applying it to de novo TNF binders is on the frontier rather than fully untried.

### Cited Findings
- GFlowNets were introduced for de novo biological sequence design with an active-learning loop that uses epistemic-uncertainty estimation and a GFlowNet generator to produce "a diverse batch of useful and informative candidates after each round," generating candidates with probability proportional to reward R(x); the method produced "more diverse and novel batches with high scoring candidates compared to existing approaches" (Jain et al., ICML 2022) — [arXiv:2203.04115](https://arxiv.org/abs/2203.04115).
- A GFlowNet learns a sampler whose terminating-state distribution is proportional to a given reward, so high-reward modes are sampled frequently but many distinct modes are covered rather than one argmax (same source) — [arXiv:2203.04115](https://arxiv.org/pdf/2203.04115).
- GFlowNet has been applied "for the first time to direct binding energy optimization in de novo diffusion-based CDR design" via AbFlowNet (Diffusion+GFlowNet fusion), improving top-1 binding energy by 38.1% and amino-acid recovery by +3.06% over diffusion-only baselines (2025) — [arXiv:2505.12358](https://arxiv.org/abs/2505.12358).
- PG-AbD unites GFlowNets with a pretrained protein language model (as a Products-of-Experts reward with a Potts model) to generate "highly potent, diverse and novel antibody candidates" (AAAI 2025) — [AAAI 2025](https://ojs.aaai.org/index.php/AAAI/article/view/34370).
- Bonobo trains a generator "on the structure prediction loss using a GFlowNet formulation in which the loss is transformed to act as a reward" for library-scale de novo antibody design (ICLR 2026) — [ICLR 2026](https://iclr.cc/virtual/2026/10022777).
- Local Search GFlowNets improve mode-finding by alternating GFlowNet sampling with local backtrack-and-reconstruct search, relevant to refining backbones/sequences (2023) — [arXiv:2310.02710](https://arxiv.org/pdf/2310.02710).

### Inferences
- Operational hook (cheapest, 4-day-feasible): do NOT train a GFlowNet over backbones from scratch. Instead run a GFlowNet over the SEQUENCE step only — fix the RFdiffusion backbone set, and let a small GFlowNet (amino-acid-autoregressive, soluble-MPNN-logit-initialized) sample sequences with reward R = scalarized {pH-switch score, cross-react delta, interface score}, replacing or augmenting soluble-ProteinMPNN sampling. This is the Jain-2022 setting (sequence GFlowNet + proxy reward) and is the only part trainable in days on 8xH100 given our existing per-design oracle scores as the proxy's training labels.
- Alternative lazy hook: skip training entirely and use GFlowNet-style reward-proportional RESAMPLING of the ~375-candidate pool we already have (importance-weight each design by exp(reward/T) and sample a diverse subset) — captures the "sample in proportion to reward, not argmax" benefit without a network.
- Risk: GFlowNet training is finicky (reward temperature, trajectory-balance loss stability); the proxy reward must be cheap (our noisy oracle scores are the natural proxy). Training a backbone-level GFlowNet in 4 days is high-risk; sequence-level or pool-resampling is low-risk.
- Untried-vs-emerging: GFlowNet-for-sequences is established; GFlowNet-for-binder-affinity is emerging (3 papers in 2025-2026, all antibody CDR, none on de novo mini-binders or pH-conditional rewards). A pH-switch-conditioned GFlowNet reward is genuinely untried.

### Gaps
- No source found applying GFlowNets to RFdiffusion-style de novo backbone generation (vs antibody CDR grafting), so backbone-level feasibility is unconfirmed.

## Multi-objective optimization matched to a STRICT priority order (lexicographic / epsilon-constraint / Chebyshev)

### Takeaway
Our strict priority (pH-selectivity > cross-reactivity > affinity) maps cleanly onto lexicographic optimization or the epsilon-constraint method; for ranking/scalarizing, the AUGMENTED Chebyshev (Tchebycheff) scalarization is preferred over linear/weighted-sum because it can reach Pareto-optimal designs in NON-convex regions that any linear scalarization provably misses — a subtle but decisive point for a 3-objective front.

### Cited Findings
- Linear (weighted-sum) scalarization only recovers Pareto-optimal solutions on CONVEX fronts; "some Pareto front solutions never maximize any linear scalarization," which is why preference-conditional multi-objective methods and Chebyshev scalarizations are needed on non-convex fronts (Multi-Objective GFlowNets discussion) — [arXiv:2210.12765](https://ar5iv.labs.arxiv.org/html/2210.12765).
- The Augmented Weighted Tchebycheff method "ensures that all efficient solutions are found, even those that are in the nonconvex region of the Pareto front," by varying the weights, addressing the linear-combination limitation (multi-objective scalarization literature) — [arXiv:2305.11774](https://arxiv.org/pdf/2305.11774).
- Multi-objective GFlowNets (MOGFN) generate "diverse Pareto optimal candidates"; MOGFN-PC is a preference-conditional GFlowNet that scalarizes into single-objective sub-problems, and goal-conditioned GFlowNets (GC-GFNs, Roy et al. 2023) "uniformly explore solutions along the entire Pareto front under challenging objective landscapes" where preference-conditioning fails on non-convex fronts — [arXiv:2210.12765](https://ar5iv.labs.arxiv.org/html/2210.12765), [arXiv:2306.04620](https://arxiv.org/html/2306.04620v2).
- Sample-efficient multi-objective molecular optimization with GFlowNets (hypernetwork over preference weights) targets diverse Pareto-optimal molecules under a limited oracle budget (2023) — [arXiv:2302.04040](https://arxiv.org/pdf/2302.04040).

### Inferences
- Operational hook (4-day-feasible, no training): implement strict priority as an epsilon-constraint ranking over the existing ~375 pool — (1) filter to designs passing a pH-selectivity threshold (primary), (2) within survivors, filter by a cross-reactivity threshold (secondary), (3) rank survivors by affinity (tertiary). This is lexicographic selection and is literally a 20-line pandas sort; it respects the hard priority exactly and never trades primary for tertiary.
- Softer hook for ranking the final <=20: augmented Chebyshev score with weights heavily skewed to the primary objective, plus a small augmentation term to break ties toward designs good on all three — gives a single defensible best-first ordering while still reaching non-convex-front designs that a weighted-sum would silently drop.
- If a GFlowNet is used for generation, condition its reward with preference weights (MOGFN-PC) biased to pH-selectivity; but given 4 days, the pool-level lexicographic/epsilon-constraint rule is the lazy, robust default.
- Untried-vs-emerging: lexicographic/epsilon-constraint selection is standard decision theory, essentially untried as an explicit, pre-registered selection rule in de novo binder campaigns (where ad hoc weighted sums dominate); augmented-Chebyshev ranking of designs is untried in this domain.

### Gaps
- None material; the methods are classical and the scalarization caveat is well-documented.

## Active learning / Bayesian optimal experimental design / batch bandits for choosing <=20 under the ~375 budget

### Takeaway
Choosing <=20 designs to maximize expected winners (or information) under a fixed pooled screen is a batch active-learning / Bayesian optimal experimental design (BOED) problem; batch Bayesian optimization, pool-based active learning, Thompson sampling, and "probability of including the optimum" acquisitions are all directly transferable, established in molecule/protein engineering, and runnable in days on our oracle scores.

### Cited Findings
- Batch Bayesian Optimization is an established method for protein sequence design: build a surrogate of the protein landscape and do model-based population selection to reduce wet-lab burden (2023) — [arXiv:2303.10429](https://arxiv.org/pdf/2303.10429).
- Pool-based active learning accelerates high-throughput virtual screening by iteratively selecting the most informative/promising candidates from a fixed library (MolPAL, Graff et al.) — [arXiv:2012.07127](https://arxiv.org/pdf/2012.07127).
- Batched Thompson Sampling obtains a batch by sampling multiple times from the posterior per round and achieves asymptotic regret equal to the sequential setting; it "naturally generates diverse candidates by independently sampling from the posterior" (Karbasi et al. 2021) — [arXiv:2110.00202](https://arxiv.org/pdf/2110.00202).
- "Probability of including the optimum" (qPO) is a batch-BO acquisition that directly maximizes the chance the selected batch contains the best design — well matched to "maximize expected winners in <=20" (2024) — [arXiv:2410.06333](https://arxiv.org/pdf/2410.06333).
- In BOED, many selection objectives (e.g., mutual information / expected information gain) are SUBMODULAR, so greedy batch construction is near-optimal with a (1-1/e) ~ 63% guarantee of the global optimum (submodular OED) — [arXiv:2602.09285](https://arxiv.org/html/2602.09285v1).
- Active-learning-assisted directed evolution (ALDE) methods maximize end-of-campaign success under "low-N" conditions (only a few dozen measurable mutants) using naturalness warm-starting and diversity-aware batch selection — [FolDE review](https://www.themoonlight.io/tw/review/low-n-protein-activity-optimization-with-folde) (secondary/aggregator source — treat as pointer, verify primary).

### Inferences
- Operational hook (4-day-feasible): fit a cheap probabilistic surrogate (Gaussian process or a small deep ensemble) that maps per-design oracle features (AF2ig ipae/pae, Boltz-2 iptm, Chai-1 iptm, PROPKA/pKAI ddG_pH, MD d_occ, novelty distance) to P(wet-lab winner). Then select the <=20 by either (a) qPO / Thompson sampling to maximize expected winners, or (b) greedy submodular maximization of expected information gain if we also want to learn for a next round. All operate on scores we already have; no new structure prediction.
- Because the screen is ONE-SHOT (no feedback round), the right objective is "maximize expected number of winners in the batch" (qPO-style, batch Thompson), NOT sequential information gain — note this distinction explicitly for the report.
- Risk: surrogate calibration on tiny labeled data; mitigate by using leave-one-out conformal (next section) rather than trusting GP variance. Low implementation risk; main risk is overfitting the surrogate to circular oracle features.
- Untried-vs-emerging: batch BO/active learning is established in protein ENGINEERING (mutant libraries) but rarely applied as the explicit SELECTION rule for a de novo binder pool with a novelty constraint; Thompson/qPO batch selection over co-folding oracle scores is essentially untried.

### Gaps
- FolDE/ALDE details came from an aggregator; the primary ALDE paper URL was not retrieved.

## Conformal prediction / calibrated uncertainty for trustworthy ranking and a joint "bind@7.4 AND not@6.0" guarantee

### Takeaway
Conformal prediction gives distribution-free, finite-sample confidence sets that remain valid even when the test set is CHOSEN based on the model (exactly the design-selection setting, "feedback covariate shift"); conformal risk control / Learn-Then-Test and conformal selection extend this to control a COMBINED risk such as a joint "binds at 7.4 AND not at 6.0" criterion with an explicit error budget — this is the single best-matched framework to our few-shot, trust-calibrated selection problem and is almost entirely untried in de novo binder design beyond one biomolecular-design paper.

### Cited Findings
- "Conformal Prediction Under Feedback Covariate Shift for Biomolecular Design" (Fannjiang et al., PNAS 2022, 119(43) e2204569119) constructs confidence sets with finite-sample guarantees that hold "for any prediction algorithm, even when a trained model chooses the test-time input distribution," precisely the dependence that arises when designed candidates are selected using the model trained on prior data — [arXiv:2202.03613](https://arxiv.org/abs/2202.03613), [PNAS/PMC9618043](https://pmc.ncbi.nlm.nih.gov/articles/PMC9618043).
- That method is recommended "for design algorithm selection," letting practitioners pick hyperparameters that trade off high predicted fitness vs low predictive uncertainty, and was demonstrated on real protein fitness data (AAV5 7-mer packaging) — [arXiv:2202.03613](https://arxiv.org/abs/2202.03613).
- Learn-Then-Test (LTT) calibrates any model so predictions satisfy explicit finite-sample statistical guarantees by reframing risk control as MULTIPLE HYPOTHESIS TESTING, without refitting (Angelopoulos et al.) — [arXiv:2110.01052](https://arxiv.org/pdf/2110.01052v4).
- Conformal Risk Control handles real-valued losses by selecting the least-conservative threshold whose upper confidence bound satisfies a target risk, built on the LTT paradigm — [MAPIE risk-control docs](https://mapie.readthedocs.io/en/latest/content/risk-control/).
- Conformal selective prediction with general risk control adapts conformal p-values to SELECTIVE settings (select a subset while controlling a chosen risk across selections), the formal template for "select <=20 while controlling false-positive rate on a joint criterion" (2026) — [arXiv:2603.24704](https://arxiv.org/pdf/2603.24704).

### Inferences
- Operational hook (4-day-feasible, high value): use split conformal on any held-out calibration set (or leave-one-out on our few labeled designs / known controls) to turn each oracle score into a calibrated predictive interval for (a) binding affinity at pH 7.4 and (b) binding at pH 6.0. Then select only designs whose conformal set satisfies "bind@7.4 = yes" AND "bind@6.0 = no." Use LTT to tune the two per-objective thresholds jointly so the batch-level probability of violating the joint criterion is <= a chosen alpha (Bonferroni-split the budget across the two conditions, or a single combined loss).
- The feedback-covariate-shift correction (Fannjiang) matters because our <=20 are selected BY the oracle scores — naive conformal would be anti-conservative; their likelihood-ratio reweighting restores validity. This is the one paper that already solved our exact statistical pathology; cite it as the methodological backbone.
- Ranking hook: order the final <=20 best-first by lower-confidence-bound on the primary (pH-selectivity) rather than point estimate, so the ranking is risk-aware (favor designs that are confidently selective over designs that are high-mean but uncertain).
- Risk: conformal needs a calibration set; with truly zero labeled TNF binders, use (i) exchangeable proxy labels from oracle cross-validation, or (ii) known positive/negative controls (a real binder and a decoy) to anchor nonconformity scores. Coverage is only as meaningful as the calibration exchangeability assumption — flag this.
- Untried-vs-emerging: conformal prediction is essentially untried in de novo binder SELECTION; only the Fannjiang biomolecular-design paper exists in-domain, and a JOINT pH-conditional conformal guarantee has no precedent.

### Gaps
- No in-domain example of conformal risk control applied to a two-condition (pH-ON/pH-OFF) selection; the combination is a novel transfer we are proposing, not a documented result.

## Submodular optimization / DPPs / Thompson sampling for diverse batch selection under a novelty constraint

### Takeaway
A diverse <=20 batch that respects the hard novelty gate is naturally expressed as submodular maximization under a matroid/knapsack constraint (near-optimal greedy, 1-1/e) or as sampling from a Determinantal Point Process whose kernel encodes design-design similarity; DPP-Thompson-Sampling unifies "high reward" with "repulsive diversity," and all three are drop-in on our existing embeddings/scores.

### Cited Findings
- DPPs model repulsive correlations: more-similar items have lower co-occurrence probability, and a more diverse subset has a larger determinant and higher selection probability — a principled diversity sampler given a similarity kernel — [arXiv:2506.21158](https://arxiv.org/html/2506.21158v1).
- Diverse mini-batch selection via DPPs "substantially improved diversity while still obtaining solutions of high quality" in de novo drug design RL (2025) — [arXiv:2506.21158](https://arxiv.org/abs/2506.21158).
- DPP-Batch Bayesian Optimization (DPP-BBO) and DPP-Thompson-Sampling (DPP-TS) induce batch diversity using DPP repulsion and achieve tighter Bayesian simple-regret bounds than classical batched Thompson sampling (Nava et al., AISTATS 2022) — [arXiv:2110.11665](https://arxiv.org/abs/2110.11665).
- DPPs for mini-batch diversification provide a natural framework to sample a diverse set from specified similarity information (Zhang et al. 2017) — [arXiv:1705.00607](https://arxiv.org/pdf/1705.00607).
- Submodular objectives (e.g., mutual information, facility-location coverage) admit greedy algorithms with a (1-1/e) ~ 63% optimality guarantee, and support cardinality/matroid constraints — directly a "pick 20 most-covering designs" formulation — [arXiv:2602.09285](https://arxiv.org/html/2602.09285v1).
- Batched Thompson sampling independently samples the posterior to form a diverse batch with regret matching the sequential case — [arXiv:2110.00202](https://arxiv.org/pdf/2110.00202).

### Inferences
- Operational hook (4-day-feasible, very high value): build a similarity kernel K over the ~375 pool from structural/sequence embeddings (e.g., TM-score matrix, ESM/ProtT5 cosine, or RFdiffusion latent distances). Then either (a) k-DPP sample 20 designs from a QUALITY-weighted DPP kernel L = diag(q) K diag(q) where q = conformal reward (this is Thompson/DPP hybrid, repels near-clones), or (b) greedy-maximize a submodular facility-location + reward objective under the cardinality-20 constraint. Both give the diverse, high-reward, non-cousin portfolio the campaign memory explicitly wants (the "diverse not cousins" lesson).
- The novelty gate is enforced as a HARD matroid/filter constraint (drop any design failing the structural-distance threshold BEFORE DPP/greedy), while DPP/submodularity handles INTRA-pool diversity among survivors — separating "novel vs known proteins" (gate) from "diverse vs each other" (DPP) is the clean formulation.
- k-DPP with a quality-weighted kernel is the single most elegant fit: it jointly trades reward and diversity in one determinant, is ~30 lines (numpy eigendecomposition + sampling), and needs only a precomputed similarity matrix we can get from existing structures.
- Untried-vs-emerging: DPP/submodular batch selection is used in drug/molecule design and BO; it is essentially untried as the explicit <=20-selection mechanism in de novo binder campaigns, and quality-weighted k-DPP over co-folding-oracle reward is a novel transfer here.

### Gaps
- None material; methods and kernels are standard and cheap.

## Optimal transport / Wasserstein barycenters to ensemble disagreeing folding oracles (AF2 / Boltz-2 / Chai-1)

### Takeaway
When AF2-ig, Boltz-2, and Chai-1 disagree on an interface/binding score, a Wasserstein barycenter of their (possibly weighted) predictive distributions gives a geometry-respecting consensus that preserves modes better than naive averaging; this is proven in multi-model forecast ensembling but is, as far as sources show, untried for folding-oracle fusion.

### Cited Findings
- Multi-model ensembling can be reformulated as the Wasserstein BARYCENTER of discrete predictive distributions; the Wasserstein distance is the optimal-transport cost and "better preserves certain geometric properties of the input ensemble distributions," and the barycenter can give different weights to members (weighted ensemble) (2023/2026, S2S forecasting) — [arXiv:2310.17933](https://arxiv.org/abs/2310.17933).
- Empirically the Wasserstein barycenter and the L2 barycenter share the same ensemble mean but differ in uncertainty representation; the L2 barycenter overspreads (underconfident) while the Wasserstein barycenter "performs better more frequently" on CRPS — [arXiv:2310.17933](https://arxiv.org/abs/2310.17933).

### Inferences
- Operational hook (4-day-feasible): treat each oracle's per-design score (or its bootstrap/seed distribution of scores) as a 1-D empirical distribution; the 1-D Wasserstein barycenter is just the weighted average of QUANTILE functions (closed form, trivially computed with numpy.percentile) — no OT solver needed in 1-D. Weight oracles by their calibrated reliability (e.g., inverse conformal interval width, or correlation with known controls). Output one fused, calibrated binding score per design that respects disagreement geometry (a bimodal AF2-vs-Boltz disagreement is preserved as spread, not washed to the mean).
- This directly addresses our "disagreeing oracles" problem and the campaign lesson that a single co-folder's cross-species delta is folder-pessimism: a Wasserstein-barycenter fuse with reliability weights is a principled upgrade over picking one folder or averaging.
- Risk: very low for 1-D (closed form). Higher-D OT barycenters (over joint score vectors) need POT/ott-jax but are still days-feasible on 8xH100; not necessary for the core win.
- Untried-vs-emerging: truly untried for folding/co-folding oracle fusion in protein design — all located precedent is in weather/forecast ensembling. This is one of the boldest, lowest-risk transfers in the catalog.

### Gaps
- No protein-design source for Wasserstein-barycenter oracle fusion; the transfer rests on forecasting precedent plus the 1-D closed form.

## Causal inference / do-calculus to separate "switch is CAUSED by the designed His contact" from confounds

### Takeaway
The isostere-ablation test we already run is informally an INTERVENTION (do(His -> isostere)); formalizing it with Pearl's causal hierarchy and counterfactual/interventional feature-importance turns "the contact correlates with the switch" into "the contact causes the switch," controlling for confounds like net charge or burial — a rigorous upgrade to our mechanism-validation gate, and untried in protein design as an explicit causal formalism.

### Cited Findings
- Pearl's three-level causal hierarchy places interventions (do-operator) above association, and counterfactuals above interventions; "counterfactuals sit at the top... and subsume interventions, allowing one to assign fully causal explanations" (deep twin networks) — [ar5iv:2109.01904](https://ar5iv.labs.arxiv.org/html/2109.01904).
- Causal mechanisms "cannot be uniquely determined from observations and interventions alone" — two models matching on interventional distributions can disagree on counterfactuals — so counterfactual claims require an explicit structural model (same source) — [ar5iv:2109.01904](https://ar5iv.labs.arxiv.org/html/2109.01904).
- Counterfactual Importance Distribution (CID) ranks feature importance by generating positive/negative counterfactuals and measuring distributional dissimilarity, a causal (not merely associational) feature-attribution method (MLR 2026) — [PMLR v307 Conti et al.](https://proceedings.mlr.press/v307/conti26a.html).
- Causality-oriented interpretability "can identify causes and effects of model decisions" rather than association, unlike traditional interpretable ML (same search cluster) — [PMLR v307](https://proceedings.mlr.press/v307/conti26a.html).

### Inferences
- Operational hook (4-day-feasible): formalize isostere ablation as a do-intervention. For each candidate switch, compute the switch metric (ddG_pH or MD d_occ) under do(His433-contact = native His) vs do(contact = isostere, e.g., His->Gln/Asn that removes titratability but preserves sterics) vs do(remove partner Asp/Glu). A significant drop ONLY under the titratability-removing intervention, holding sterics fixed, is the interventional signature that the switch is caused by the protonation contact, not by a confound (burial, net charge, a constitutive salt bridge). This matches our memory's warning that constitutive salt bridges pass isostere tests yet lack selectivity — a matched-sterics counterfactual (His->Gln) is the discriminating intervention.
- Add a simple structural causal DAG: {design features} -> {His contact geometry} -> {protonation-dependent ddG} -> {pH selectivity}, with net-charge/burial as back-door confounders; use it to justify which ablations block which back-door paths. This is a one-page formalization, not new compute.
- Risk: very low compute; the "causal" content is in experimental design of the ablations, not a heavy estimator. Main risk is over-claiming identifiability from a handful of in-silico interventions.
- Untried-vs-emerging: explicit do-calculus / counterfactual framing of a design mechanism test is untried in de novo binder design; our pipeline already does the intervention, so formalizing it is low-cost and high-rigor.

### Gaps
- No protein-design paper uses do-calculus for switch-mechanism attribution; the transfer is ours, grounded in general causal-inference theory.

## Information-geometric / topological measures of portfolio coverage and robustness

### Takeaway
Persistent homology (TDA) and information-geometric diversity measures quantify how well a <=20 portfolio COVERS the reachable design space and how robust that coverage is to perturbation — turning "diverse portfolio" from a vibe into a number; persistent homology is established for protein CONFORMATIONAL analysis, but using it to measure DESIGN-LIBRARY coverage is a novel transfer.

### Cited Findings
- Persistent homology "embeds geometric information into topological invariants" and has characterized protein structure, flexibility, and folding with "excellent consistence between persistent homology prediction and molecular dynamics simulation" (Xia & Wei 2014) — [arXiv:1412.2779](https://ar5iv.labs.arxiv.org/html/1412.2779).
- A hierarchical-clustering + persistent-homology method detects regions of interest in protein CONFORMATIONAL space from coarse-grained search, "without requiring explicit enumeration of all possible states" — useful for summarizing ensemble/library coverage — (TDA protein analysis) — [arXiv:1412.2779](https://ar5iv.labs.arxiv.org/html/1412.2779).
- DPPs provide an information-theoretic diversity score (log-determinant of the similarity kernel) that doubles as a coverage/robustness metric for a selected set — [arXiv:1705.00607](https://arxiv.org/pdf/1705.00607).

### Inferences
- Operational hook (4-day-feasible): two cheap coverage metrics on the final <=20. (1) Information-geometric: log-det of the DPP kernel L over the selected set = a single "portfolio diversity" scalar; maximize it subject to the reward/novelty constraints (identical machinery to the DPP selection section). (2) Topological: build a Vietoris-Rips/alpha complex over design embeddings (ripser/giotto-tda, minutes to run) and report the persistence diagram / number of persistent H0 clusters as a measure of how many DISTINCT fold-lineages the portfolio spans — guards against submitting 20 variations of one winning seed (the "diverse not cousins" lesson).
- Robustness angle: recompute the persistence diagram after dropping each design (leave-one-out); a portfolio whose topology is stable to removals is robust coverage, one that collapses depends on a single design.
- Risk: low compute (ripser is fast on 375 points); interpretation risk — persistence features on learned embeddings can be metric-dependent, so report alongside the simpler log-det number.
- Untried-vs-emerging: persistent homology is established for protein conformational/folding analysis; applying it (or log-det coverage) to score a DE NOVO DESIGN PORTFOLIO's diversity/robustness is untried.

### Gaps
- No source directly computes TDA coverage of a de novo design LIBRARY (vs conformational ensembles of one protein); the transfer is inferential.

## Equivariance beyond SE(3): gauge / geometric-algebra / Clifford equivariant representations for interfaces

### Takeaway
Clifford/geometric-algebra (Clifford-group) equivariant networks and geometric vector perceptrons generalize SE(3)-equivariance and give compact, frame-independent interface representations; gains reported so far are modest (~2% accuracy with fewer parameters) and these are emerging-not-mainstream in protein ML, so the realistic 4-day win is using an existing geometric-algebra/GVP embedding as a better interface FEATURIZER for our surrogate/DPP kernel, not training a new equivariant model.

### Cited Findings
- Clifford Group Equivariant Neural Network (CGENN) layers operate on multivectors of a Clifford algebra in an E(n)-equivariant way (generalizing beyond SE(3)); for protein structure prediction they increased accuracy "by up to 2.1% with fewer trainable parameters" and gave interpretable geometric outputs (Pepe et al., 2024) — [PMLR v233 Pepe et al.](https://proceedings.mlr.press/v233/pepe24a.html).
- Geometric Vector Perceptrons (GVP-GNN) augment GNNs with layers that do "expressive geometric reasoning over Euclidean vector features," a widely used SE(3)-aware protein-structure representation (Jing et al., ICLR 2021) — [arXiv:2009.01411](https://arxiv.org/abs/2009.01411).
- GLGENN provides parameter-light Clifford/geometric-algebra equivariant networks (2025), continuing the geometric-algebra line — [arXiv:2506.09625](https://arxiv.org/html/2506.09625v1).

### Inferences
- Operational hook (4-day-feasible, pragmatic): do NOT train a Clifford/gauge-equivariant model from scratch. Instead extract a GVP or CGENN embedding of each design's INTERFACE region (the receptor-binding-groove contact residues at the two-protomer interface) using a pretrained geometric encoder, and feed those frame-invariant features into (a) the conformal surrogate (better winner-prediction) and (b) the DPP/TDA similarity kernel (better diversity measurement). Equivariant features are more sample-efficient than raw coordinates, which matters in our few-label regime.
- The reported gains are small (~2%), so equivariance-beyond-SE(3) is a representation refinement, not a step-change — position it as a nice-to-have featurizer, not a core bet, for a 4-day window.
- Risk: moderate — pretrained CGENN weights for interfaces may not be readily available; GVP is more mature and lower-risk. Training any equivariant model in 4 days is not advised.
- Untried-vs-emerging: Clifford/geometric-algebra equivariance is emerging in protein structure prediction (a handful of 2024-2025 papers); gauge-equivariant (beyond global E(n), truly local gauge) interface models are essentially untried in binder design.

### Gaps
- No source quantifies a binding/interface-specific gain from geometric-algebra equivariance (the ~2% is general structure prediction), so the expected benefit for interface scoring is unquantified.

## Cross-cutting synthesis (priority recommendations for the report)

### Takeaway
Ranked by (impact on our two decisive bottlenecks) x (4-day feasibility): the highest-value, lowest-risk transfers are conformal selection under feedback covariate shift, quality-weighted k-DPP / submodular batch selection, 1-D Wasserstein-barycenter oracle fusion, and lexicographic/epsilon-constraint priority enforcement — all operate on our EXISTING per-design outputs with no model training. GFlowNet generation and equivariant representations are higher-risk, higher-ceiling bets better scoped as stretch goals.

### Cited Findings
- Conformal under feedback covariate shift is the one framework already validated on real protein design-selection data (Fannjiang, PNAS 2022) — [arXiv:2202.03613](https://arxiv.org/abs/2202.03613).
- DPP/submodular diverse batch selection is validated in de novo drug design and batched BO with tighter regret than vanilla Thompson sampling — [arXiv:2506.21158](https://arxiv.org/abs/2506.21158), [arXiv:2110.11665](https://arxiv.org/abs/2110.11665).
- Wasserstein-barycenter ensembling preserves distributional geometry better than L2 averaging — [arXiv:2310.17933](https://arxiv.org/abs/2310.17933).
- GFlowNet-for-binders is emerging (3 papers, 2025-2026) and antibody-only, implying de novo mini-binder + pH-conditional reward is frontier, not settled — [arXiv:2505.12358](https://arxiv.org/abs/2505.12358), [AAAI 2025](https://ojs.aaai.org/index.php/AAAI/article/view/34370), [ICLR 2026](https://iclr.cc/virtual/2026/10022777).

### Inferences
- A single coherent 4-day pipeline: (1) hard-gate the ~375 pool on structural novelty; (2) fuse AF2ig/Boltz-2/Chai-1 per-design scores via reliability-weighted 1-D Wasserstein barycenter; (3) calibrate each design's pH-7.4-bind and pH-6.0-no-bind predictions with split/feedback-shift conformal + LTT to control the JOINT error; (4) lexicographically filter by pH-selectivity -> cross-reactivity, then (5) k-DPP-select <=20 from survivors with quality = conformal LCB on selectivity and kernel = interface-embedding similarity; (6) report portfolio coverage via log-det + persistent-homology cluster count; (7) formalize the switch-causality claim with a matched-sterics (His->Gln) do-intervention. Every step runs on existing outputs in Python; none requires retraining a generator.
- GFlowNet (sequence-level or pool-resampling) and equivariant interface featurization are optional accelerators layered on top if time permits.

### Gaps
- Feasibility estimates are reasoned, not benchmarked on our specific 8xH100 setup; actual wall-clock depends on how many embeddings/structures must be recomputed for the similarity kernel.
