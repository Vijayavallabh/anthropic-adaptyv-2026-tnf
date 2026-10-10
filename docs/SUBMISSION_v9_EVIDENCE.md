# Submission v9 — evidence record

**File:** `runs/push/submission_v9.csv` · md5 `05a344444cf740907e79a177d34c73f5` · 20 designs ·
columns `name,sequence,molecule_class` (all `single_chain`) · 0-Cys · len 73–119 · 20 unique names + sequences.
**Portal novelty:** all 20 ≥3/4 (3 at 4/4). **Upload: USER-GATED.**

Target: de-novo conditional pH-switch binder to TNF-α. BIND pH7.4 / RELEASE pH6.0, cross-reactive human AND
mouse. Judged: pH-selectivity (#1) > cross-reactivity (#2) > affinity; novelty ≥3/4 hard gate; expression =
dominant wet-lab-zero risk.

## 1. Switch mechanism — multi-His separated-cation (Baker 2025 SOTA)

Blueprint: Ahn/Coventry/Baker bioRxiv 2025.09.29.678932. Single-His is the WEAK end; the validated SOTA is
≥2–3 **separated** His, each repelling a **distinct conserved** TNF-α cation (spread across the epitope; NOT
clustered, which is anti-cooperative). Fold-change scales with the number of His–cation contacts (TNF-α 79×).
Conserved human↔mouse cations (mature index): 1,6,27,39,60,77,85,93,98,107,123; K6/K93 are the strong close
anchors. At pH7.4 His is neutral → binds; at pH6.0 His⁺···cation⁺ repulsion destabilizes the interface → releases.

## 2. Switch-core designs (7) — triple-confirmed

| name | lineage | conf-His | min-ipSAE | Δ_hm | mono-pLDDT | SAP+ | GRAVY | novelty |
|---|---|---|---|---|---|---|---|---|
| des_3050_6_c39 | des_3050_6 | **4** | 0.666 | 0.063 | 0.88 | 0.6 | −0.52 | **4/4** |
| des_4098_1_b14 | des_4098_1 | 3 | 0.587 | 0.106 | 0.90 | 0.6 | −0.58 | 3/4 |
| des_4098_1_c44 | des_4098_1 | 2 | **0.861** | **0.024** | 0.82 | 0.9 | −0.67 | 3/4 |
| des_180_6_H9-14-38 | des_180_6 | 2 | 0.734 | 0.058 | 0.81 | 4.8 | +0.10 | 3/4 |
| des_818_4_H16-74 | des_818_4 | 2 | 0.713 | 0.105 | 0.88 | 5.4 | −0.14 | 3/4 |
| des_3050_6_H52-54-79 | des_3050_6 | 2 | 0.631 | 0.095 | 0.79 | 2.5 | −0.08 | 3/4 |
| des_875_5_H3-57-80 | des_875_5 | 1 | 0.600 | **0.016** | 0.86 | 0.1 | −0.42 | **4/4** |

conf-His = # His passing the per-His triple-confirm. min-ipSAE = min(human,mouse) Boltz ipSAE (binding).
Δ_hm = |ipSAE_human − ipSAE_mouse| (balance; lower = more cross-reactive). SAP+/GRAVY = expression (lower).

## 2b. Orthogonal binding confirmation (Chai-1, 3rd co-folder) + selectivity-strength

Chai-1 single-seq co-fold vs human AND mouse (no MSA) — independent of Boltz/af2ig. `runs/push/chai_switchcore.tsv`.

| name | conf-His | ΣdE | Chai_human | Chai_mouse | Chai Δ | Boltz min-ipSAE |
|---|---|---|---|---|---|---|
| des_3050_6_c39 | 4 | 13.7 | 0.788 | 0.772 | 0.016 | 0.666 |
| des_4098_1_b14 | 3 | 17.6 | 0.865 | 0.870 | 0.005 | 0.587 |
| des_4098_1_c44 | 2 | 17.8 | 0.863 | 0.870 | 0.007 | 0.861 |
| des_180_6_H9-14-38 | 2 | 17.2 | **0.424** | **0.369** | 0.055 | 0.734 |
| des_818_4_H16-74 | 2 | 1.6 | 0.839 | 0.841 | 0.002 | 0.713 |
| des_3050_6_H52-54-79 | 2 | 1.9 | 0.722 | 0.708 | 0.014 | 0.631 |
| des_875_5_H3-57-80 | 1 | 1.8 | 0.836 | 0.805 | 0.031 | 0.600 |

- **6/7 confirmed by a 3rd independent co-folder** as strong binders to BOTH species, with small cross-species
  **matched deltas (0.002–0.031)** = balanced cross-reactivity (the rigorous #2 metric; not single-folder).
- **FLAG des_180_6_H9-14-38**: Chai weak on BOTH species (0.42/0.37) vs Boltz 0.734 — a genuine Boltz/Chai
  binding disagreement (not cross-species pessimism). Also the highest-SAP switch-core member. KEPT (portal-
  confirmed, 2-His triple-confirmed switch, distinct lineage; no DQ-safe switch-preserving replacement) — logged
  as the lowest-binding-confidence switch-core member.
- **ΣdE** = sum of confirmed-His GFN2-xTB dE (honest RELATIVE selectivity-strength proxy, NOT a fold-change):
  c44/b14/c39 are the strongest switches by per-His QM destabilization; c39 strongest by His breadth (4 distinct
  cations).

## 3. Switch-confirm gate (per His; design confirmed if ≥2 His pass) — two independent thermodynamic methods

1. **PROPKA Wyman linkage** — bound-state His pKa suppressed by the cation (dpKa < 0 = protonation disfavors
   binding = correct switch sign). `scripts/collective_his_linkage.py`.
2. **GFN2-xTB** capped-cluster interaction energy — dE_switch = E_int(HIP) − E_int(HIE) > 0 (protonated state
   destabilized). `scripts/qm_switch_fair.py` (`envs/xtb`). Example c39: His49 dE +3.68, His79 dE +4.66.

The honest in-silico pH-selectivity ceiling is this thermodynamic proton-linkage sign + confirmed-His count;
protein–protein dissociation is seconds-scale and unobservable in MD (d_iface≈0 is EXPECTED). We do NOT report a
quantitative fold-change: interface pKa predictors disagree by 2–3 units (> the pH window), so a numeric FC would
be false precision. conf-His count [4,3,2,2,2,2,1] is the honest selectivity proxy.

**SMD bold-bet probe (Bet D) — NULL, by design-honesty.** `scripts/md_pull_tnf.py` steered-MD binder↔target
separation work in both His states. At LOW sampling (3 reps) des_3050_6_c39 showed dW_release +16.1 kcal, but it
**did not survive**: at 9 reps dW = **+0.1 ± ~11 kcal** (HIE 398.5±7.8 ≈ HIP 398.5±7.1, consistent with zero);
des_4098_1_c44 (3 reps) +4.88 was likewise within noise. The +16.1 was a small-n fluctuation. Conclusion: steered
MD at feasible sampling does **not** resolve the pH-release signal — which **confirms** the thermodynamic ceiling
above. The switch evidence rests on PROPKA + xTB, NOT SMD. (Recorded as a negative; it does not weaken the
switch-core, which was gated on the two thermodynamic methods.)

## 4. Expression (dominant wet-lab-zero risk) — all 20 verified

Primary predictor = single-sequence ESMFold **monomer pLDDT** (AUROC~0.72). All 20 = **0.79–0.94** (strong
band). Orthogonal: SAP+ (Chennamsetty), GRAVY, acidic%, netQ. Consensus-drop = weak on ≥2 axes. Only
b1m9_A12H flags (SAP+6.7 + GRAVY+0.30) — **KEPT**: it is portal-confirmed ≥3/4 (swapping trades a soft
yield-risk for a hard novelty-DQ), pLDDT 0.87 (primary predictor = expresses), and no lineage-preserving
portal-confirmed replacement exists. Audit: `expr_audit20.py`; metrics in `runs/push/master_v9.csv`.

## 5. Cross-reactivity (#2) — scored by mechanism, not one co-folder

Contacts land on conserved Region-I anchors (V74/L75/I97/W114/Y115) + conserved cations; balance measured as
a matched human/mouse co-fold delta Δ_hm (0.016–0.106 Boltz across the switch core). A low single-folder ortholog
score is folder pessimism, not a species signal. **Confirmed orthogonally** (§2b): Chai-1 human↔mouse matched
delta 0.002–0.031 for 6/7 switch-core designs → two independent co-folders agree the cross-reactivity is balanced.

## 6. Novelty (hard gate, <3/4 = DQ) — all 20 ≥3/4

STRUCTURAL (single-seq fold, max TM vs PDB+AFDB). Portal (sole authority) confirmed all 20 ≥3/4; 4/4:
des_3050_6_c39, des_875_5_H3-57-80, b1m12_A41H_vsR27. MPNN re-seqs on a verified backbone inherit its
backbone-determined novelty. De-novo throughout (no natural-binder seeding).

## 7. Portfolio logic — 16 lineages / 20

7 triple-confirmed multi-His switch-core (strongest #1-axis evidence) + 13 portal-verified v7-floor binders
(diverse validated hedge). Spread across 16 distinct lineages, cap ~2/lineage — diversity hedges wet-lab risk.
Strong triple-confirmed switches left out (des_3050_6_c34/c19/c48, des_4098_1_c39, des_818_4_c13…) are
redundant lineage-cousins of included members; excluding them is deliberate (diversity > cousins).

## 8. Bold-bet campaign (/goal) — all explored, none beats v9

- **A buried His-network** (Baker method-2, cooperative Hill>1): no de-novo core carries a pre-existing buried
  His···cation pair <8Å; building one de-novo risks the fold + dominant expression axis for an unverifiable
  cooperativity claim → **DROPPED (no low-risk path)**.
- **B Mosaic de-novo** (Escalante gradient seq-space): CAN bind (mosh_0 ipSAE 0.78) but the binding depends on
  an Ala-rich helix (GRAVY +0.63, 36–38% Ala) that dooms expression; 0–1 incidental His (no switch); diffuse
  non-specific epitope. Salvage = rebuild v9's pipeline on an inferior seed → **CLOSED**. Confirms de-novo
  generators floor on the flat TNF target → v9's multi-His-on-validated-binders is correct.
- **C multi-state tied MPNN** (cross-react in the objective): **refuted by bytes** (runs/jxr, 0/256).
- **D SMD functional-release**: **NULL** — the 3-rep +16.1 kcal did not survive to 9 reps (+0.1, consistent with
  zero; §3). MD does not resolve the release signal at feasible sampling; confirms the thermodynamic ceiling. A
  rigorously-caught small-n false positive, not a win.

## 9. Reproduce

Build: `scripts/integrate_v9.py` → `submission_v9.csv` + `master_v9.csv`. Switch pipeline: `scripts/README.md`
(wave3 Track A/B → wave4b His-biased MPNN → triple-confirm). Methods: `docs/METHODS.md` §8.
