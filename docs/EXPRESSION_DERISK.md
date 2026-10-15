# Expression de-risking dossier — treating "will it express" as a quantified prediction, not a solved fact

> **Note (SESSION 15):** this is a v11-era snapshot. The *methodology* (orthogonal multi-signal expression
> consensus) is current and still applied, but the "current best = v11" framing below is superseded — the
> **confirmed primary is `submission_v25.csv`** (v24 with the direct-FE swap b14→c44; portal-parsed 20/20),
> with `submission_v24.csv` the conservative alternative (all 20 portal-PASSED) and `submission_v23.csv`/`v19`
> fallbacks (see AGENTS.md / README.md / session-handoff.md / docs/SWITCH_FE_findings.md).

**Why this doc.** Expression is the dominant wet-lab-zero risk: a design that cannot be manufactured scores
zero regardless of its pH-switch. Our prior confidence rested mainly on **one** predictor — ESMFold monomer
pLDDT, whose published de-novo-success AUROC is ~0.72 (Garcia/Dixit/Rocklin 2026, PMC12817478). A single
~0.72-AUROC predictor saying "all 20 pass" is a hopeful prediction, not proof. This dossier replaces that
binary with a **multi-signal, orthogonal consensus** and an **honest expected-yield with a risk band**.

## The six signals and why they are orthogonal

| # | signal | modality | independent of ESMFold? |
|---|---|---|---|
| 1 | ESMFold monomer pLDDT | ESM2 language-model folder | — (the baseline) |
| 2 | Boltz-2 monomer pLDDT + pTM | AF3-style **diffusion** folder | YES — different architecture & training |
| 3 | cross-folder **TM(ESMFold,Boltz)** | structural agreement of 1↔2 | YES — the Baker/Rocklin gate: two independent folders must agree on the SAME fold |
| 4 | structure-based SAP | exposed-hydrophobic aggregation surface (Shrake-Rupley SASA on the fold) | YES — a biophysical surface property, not a confidence score |
| 5 | Protein-Sol scaled solubility | biophysical/compositional model trained on **experimental** eSOL E. coli cell-free solubility (Hebditch 2017; Niwa 2009) | YES — **no neural net, no ESM**; a measured-solubility model |
| 6 | developability liabilities | net charge / GRAVY / charge-runs / N-glyc / low-complexity | YES — sequence rules |

Signals 1–3 ask *"does it fold?"*; 4–6 ask *"does the folded/charged sequence stay in solution and express?"*
Errors in these are largely independent, so agreement across them is far stronger than any one alone.

Tooling (all read-only, committed): `scripts/expr_consensus.py` (folders + TM + aggregate), `scripts/expr_yield_dist.py`
(Poisson-binomial risk), `scripts/sap_scan.py`/inline SAP, `scripts/solubility_predict.py` (Protein-Sol,
stdlib, `--selftest`; validated **byte-identical** to the official Perl tool, max|Δ|=0.0000 on 9 controls + 20
designs; control ordering PASS: min positive 0.560 > max membrane-negative 0.505), `scripts/developability_scan.py`.

## Headline results (submission_v11 — the retained current best, md5 `ef72d8e09f11d75d7aef959b4b18a2eb`)

- **Two architecturally independent folders agree** on the fold for **18/20** designs at TM ≥ 0.85 (all 20 ≥ 0.745).
  The fold is no longer one predictor's claim; ESM2 and AF3-diffusion concur.
- **All 20 pass the experimental-data solubility check**: Protein-Sol scaled solubility 0.585–1.000, every design
  **above** the 0.446 eSOL population reference (i.e. predicted more soluble than an average E. coli expression
  protein). None is a hard solubility risk.
- **Expected number expressing ≈ 17.1 / 20** (literature-anchored per-design probabilities), band **[15.5, 18.1]**.
  Poisson-binomial over the 20 independent Bernoulli(p_i): **P(≥14)=0.983, P(≥16)=0.848, mode=17**.
- **Tiers: A=12 (clean + 2-fold agree), B=8 (one soft axis), C=0 red flags.**
- **Switch-survival hedge:** P(≥1 of the 7 fold-robust-His switch designs expresses) = **1.0000**. Even under the
  ~0.72 per-predictor uncertainty, the pH-switch mechanism (the #1 judged axis) is essentially guaranteed to be
  represented in the expressed set. The single flagship c19 is a 0.76 bet on its own, but the lineage hedge
  (des_3050_6 + des_4098_1 + des_818_4 + TNFpH0x families) makes switch-survival a near-certainty.

## What the multi-signal consensus CAUGHT (the point of the exercise)

Relying on any single predictor would have misled us. Two concrete catches:

1. **False alarm corrected — tnf_K90.** `developability_scan` flags it **SEVERE** on net charge (−15.1). On that
   flag alone it looked like the weakest member and a swap target. But the orthogonal **experimental** solubility
   model scores it **0.883 (rank 17/20, among the most soluble)** — acidic proteins are typically highly soluble —
   and its **exposed-hydrophobic surface is the cleanest of all 20** (SAP patch 174 Å² vs portfolio median 384).
   The charge "SEVERE" is a magnitude heuristic that the measured-solubility model explicitly clears. tnf_K90 is
   **kept**; swapping it on the single flag would have been the exact predictor-over-trust error. Its only genuine
   weakness is mid-pack cross-reactivity (iptm 0.59/0.54), not expression.
2. **Hidden risk exposed — A20V candidate.** A rejected swap candidate (A20V_1A20H_K197) had the *highest* fold
   confidence of anything tested (ESMFold 0.95 / Boltz 0.97), which pLDDT-alone would read as safest. The
   structure-based SAP revealed the **worst exposed-hydrophobic patch of all** (534 Å², 25% exposed) and the lowest
   solubility (0.63). pLDDT alone would have picked a surface-aggregation-prone design. It was excluded.

**Convergent weak group:** three orthogonal axes (developability hydrophobic flags, structure-SAP, Protein-Sol)
agree the relatively weaker-expressing members are the **b\* binder-scaffold lineage** (b1m12/b4m2/b1m9/b1m7/b2m0),
driven by exposed aromatic/aliphatic stretches — not the charge-flagged tnf_K90. All b\* members still clear the
solubility population reference (0.585–0.78); they are the lowest *within* a portfolio that is soluble throughout.

## Swap candidates evaluated — none strictly improves v11, so v11 is RETAINED

The consensus flagged **b1m12_2A67H_K112_r1** as the apparent weakest member (**weakest solubility 0.585**,
**weakest switch** both-fold linkage −0.208, redundant 7-member b-lineage). I fully validated a replacement
(**TNFpH07**: solubility 0.840, aggregate both-fold switch −0.69/−0.80) on the identical pipeline — folds, Boltz
+ Chai complexes, PROPKA+xTB per-His confirmation, and cross-species iptm. The disciplined verdict: **do not swap.**

| axis | b1m12_2A67H_K112_r1 (incumbent) | TNFpH07 (candidate) | honest read |
|---|---|---|---|
| switch, **strict per-His fold-robust** (PROPKA∧xTB on BOTH folds) | **0** | **0** (xTB confirms His49 on Boltz but His5 on Chai — different His, fold-dependent) | **TIE** — neither is a mechanistic multi-His switch |
| switch, aggregate both-fold Wyman linkage | −0.208 | −0.69 / −0.80 | candidate stronger *on the aggregate metric only* |
| **cross-reactivity** iptm (h / m) | **0.92 / 0.916 (strong)** | 0.642 / 0.627 (mid) | **incumbent much stronger — this is the catch** |
| solubility (Protein-Sol) | 0.585 | 0.840 | candidate better (but both above population ref) |
| diversity | b1m12 (unique variant) | 3rd near-identical TNFpH0x cousin | candidate worse |

**Why no swap:** the portfolio's #1 (pH-selectivity) *rank* is set by its strongest switches — flagship c19
(4 fold-robust His, linkage −1.274) and the fold-robust members — which this marginal slot does not affect. On the
**strict** fold-robust metric the incumbent and candidate **tie at 0**, so the swap buys only a stronger *aggregate*
linkage on a non-flagship member. Against that marginal #1 gain, the swap would **sacrifice a genuinely strong
cross-reactive binder (iptm 0.92/0.916 → 0.63)** — a real loss on the #2 judged axis — and add a third
near-duplicate TNFpH0x fold. That is a lateral-to-negative trade, not a strict improvement, so forcing it would be
churn. **No v11 member is strictly dominated by any validated candidate** (des_495_1_H78: top solubility + strong
x-react but negligible switch −0.068; A20V: strong x-react but the *worst* exposed-hydrophobic patch and lowest
solubility; des_962_5_H37: anti-switch +0.627). v11 stands.

The candidate evidence is retained at `runs/push/rejected/submission_v12_tradeoff_NOT_ADOPTED.csv`
(md5 `e32ec77a68b200cdb18a88cc0c19ccbb`) for provenance — explicitly **not** a submission.

**Bottom line:** the critique asked that expression be treated as a prediction, not solved. This six-signal
consensus does exactly that — and, applied honestly, it (a) *confirms* v11 is well-de-risked (C=0, all 20 soluble,
E≈17/20, switch-survival ≈1.0), (b) *corrects* a false alarm that would have wrongly ejected tnf_K90, and (c) shows
no swap strictly improves the portfolio. v11 (`ef72d8e09f11d75d7aef959b4b18a2eb`) remains current best; v9
(`05a344444cf740907e79a177d34c73f5`) byte-locked fallback. Upload remains USER-GATED.

## Honest caveats

- These are all **in-silico predictions**. None is a wet-lab assay; the strongest single de-novo-success AUROC here
  is ~0.72. The value is in **orthogonal agreement** (independent errors cancel) and in stating risk as a
  **distribution**, not a pass/fail.
- The per-design p_i are **literature-anchored**, not fitted to the Adaptyv assay; the band and the Poisson-binomial
  reflect that. The robust, prior-free headline is the **tier count (A=13/B=7/C=0)** and **all-20-above-population
  solubility**.
- Protein-Sol predicts **recombinant-expression solubility**, not amyloid kinetics; it is used for the former only.
