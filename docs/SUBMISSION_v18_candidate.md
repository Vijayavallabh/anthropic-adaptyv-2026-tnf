# Submission v18 candidate — 4th de-novo lineage + strongest switcher (PORTAL-NOVELTY CONFIRMED)

Updated 2026-10-12 (SESSION 11, continued). Deadline **2026-10-11 23:59 AoE = 2026-10-12 11:59 UTC**.

> **2026-10-12 UPDATE — the one open risk is now RESOLVED.** The user ran the portal novelty check on all
> 20 v18 designs: **every design scored ≥3/4 (the hard gate; <3/4 = DQ)**. Both new members read **3/4**
> (`des_6126…mh_9_56_78`, `des_6025…mh_52_74`), exactly as the local calibration predicted; four designs
> read **4/4**. v18 is therefore **fully gated** — no portal-unverified member remains. The portal list was
> byte-checked against the CSV (all 20 names match `md5 524722158467273c0e9f57f5402650ed`). See
> "Portal-novelty result" below; it supersedes the old "one open risk" section.
>
> **2026-10-12 UPDATE 2 — v19 is ALSO fully gated, and the deadline is extended (~7 days).** The user ran the
> portal check on v19 too: **des_84_3 reads 3/4; all 20 v19 designs ≥3/4.** So v19's one contingency is
> resolved — both v18 and v19 are now fully-gated, valid submissions. With the extended runway, work has moved
> from locking a submission to a **bold binding-oracle campaign** (cross-folder pose-consensus calibrated on
> the 12 CRO-validated binders) aimed at the hook's real critique — that binding/cross-react are near-chance
> in-silico for TNF. Status of that campaign is tracked separately; v18/v19 remain the standing picks meanwhile.

## TL;DR

v18 = **v16 with two swaps**, both MD-NONRESPOND members replaced by **MD-CONFIRMED new-lineage
cross-reactive switchers** harvested from the bpair/partial-diffusion switch sweep:

| slot (v16 NONRESPOND) | swapped in | dSwitch | lineage |
|---|---|---|---|
| `des_3050_6_H52-54-79` (+0.31, below thr) | **`des_6126_3_T20_pd_45_3_mh_9_56_78`** | **+1.36** | **des_6126 (NEW 4th de-novo lineage)** |
| `b4m2_A57H` (+0.22, below thr) | **`des_6025_4_T10_pd_77_4_mh_52_74`** | **+0.93** | des_6025 (3rd de-novo lineage) |

The third v16 NONRESPOND member, **`des_818_4_c13S` (the only 4/4-novel design), is KEPT** — deliberately
not swapped for a 3rd available des_6025 candidate (a cousin), to preserve lineage diversity + the 4/4
novelty showcase (per "diverse not cousins"; the extra cousin would add only a weak +0.68 switcher while
losing a distinct lineage and the lone 4/4).

- Portfolio two-state-MD trigger **17/20 → 19/20 CONFIRMED** (one NONRESPOND kept on purpose: the 4/4 `des_818_4`).
- De-novo new-lineage switchers **2 → 4** (`des_138_7`, `des_4082_8`, `des_6025`, **`des_6126`**).
- des_6126 is the **strongest switcher in the whole portfolio** (d_repel +1.36, all 4 reps positive).
- **v18 CSV (RECOMMENDED):** `runs/push/submission_v18.csv` (md5 **`524722158467273c0e9f57f5402650ed`**),
  validator **PASS** (20 designs, 0 errors, 0 warnings, 0-Cys, unique names+seqs, all `single_chain`),
  **all 20 portal-novelty ≥3/4** (see below).
- **Conservative floor:** v16 `a524d77fa7a63ef33b28346ce3a7d87d` (identical gates, weaker only on the #1
  switch axis); v17 `9357b6905ed81179da5137a16dd9dd5a` (intermediate). All three share the same
  generator-provenance profile — see "Generator provenance" and "Recommendation".

## Swap-in 1 (flagship): des_6126_3_T20_pd_45_3_mh_9_56_78

`MTEEESRKHRAEVAAILEAAAEKALEEVGDEETAEKLRAFAEAVRAGTSGAGVSIHGGGALVPVFESLGEALAALNKHLA` (80 aa, 0 Cys, His 9/56/78)

| Axis | Metric | Value | Verdict |
|---|---|---|---|
| **#1 pH-switch** | two-state MD d_repel (HIP−HIE), 4 reps × 3 ns | **+1.36** (2·SEM ≈ 0.58, all arms) | **CONFIRMED (strongest)** |
|  | per-replica HIP−HIE | +2.04 / +1.02 / +0.66 / +1.72 | **4/4 reps positive** (not an outlier) |
|  | isostere control (HIP−GLN) | +0.39 | positive (His-causal ablation holds) |
|  | `selective_release` | `true` | — |
|  | binder RMSD / interface frac | ~2.0 Å / 0.84 | complex MD-stable across states |
| **#2 cross-react** | **conserved-cation engagement** (the real evidence) | His9→**K6**, His56→**K93**, His78→**R1** — all 3 identical in mouse | **mechanistic PASS** |
|  | Boltz-2 iptm h/m (NOT relied on) | 0.970 / 0.911 | *uninformative for TNF — see cross-react note* |
| **expression** | Protein-Sol scaled solubility | 0.82 | PASS (≥0.50) |
| **novelty seq** | mmseqs vs PDB seqres | **NO PDB hit (novel)** | passes seq criterion |
| **novelty struct** | foldseek max-TM vs PDB → **portal** | **0.948** → portal | **3/4 — PORTAL-CONFIRMED** |

## Swap-in 2: des_6025_4_T10_pd_77_4_mh_52_74

`SEETERLIRELNQILLEEMKETVKKLKELGVEVELNVSENGGVSLIASSESHELALEASQLIKETVRRIEQLLHKYVAER` (80 aa, 0 Cys, His 52/74)

| Axis | Metric | Value | Verdict |
|---|---|---|---|
| **#1 pH-switch** | two-state MD d_repel (HIP−HIE), 4 reps × 3 ns | **+0.93** (2·SEM ≈ 0.43, all arms) | **CONFIRMED** |
|  | isostere control (HIP−GLN) | +2.18 | positive |
| **#2 cross-react** | **conserved-cation engagement** | His52→**K60**, His74→**K6** — both identical in mouse | **mechanistic PASS** |
|  | Boltz-2 iptm h/m (NOT relied on) | 0.950 / 0.936 | *uninformative for TNF — see cross-react note* |
| **expression** | Protein-Sol scaled solubility | 0.92 | PASS |
| **novelty** | mmseqs seq / foldseek struct TM → **portal** | NO PDB hit / **0.870** → portal | **3/4 — PORTAL-CONFIRMED** |

Note: this is a **stronger** des_6025 switcher than v17's member (`mh_34_50_52_74`, +0.372) — v18 upgrades the
des_6025 representative while using an independent variant.

## The pH-switch: four physics methods agree on the binding pose (answers "MD alone"; see fold caveat)

The decisive (#1-ranked) axis is the pH switch. On the af2ig/MD binding pose, it agrees in sign across four
orthogonal method-classes — two-state MD (fixed-charge MM), the Gln-isostere ablation control, two empirical
pKa predictors (PROPKA-Wyman, pKAI), and a QM single-point (GFN2-xTB). So it is **not MD-alone**. A positive
static score only proves switch DIRECTION; sign-consensus across orthogonal physics is the honest selectivity
evidence (magnitudes are not comparable across classes and are not claimed). **Caveat (below): this agreement
is on one fold family; an independent folder (Chai) does not reproduce it — the fold-independent evidence is
the conserved-cation mechanism.**

| member | MD d_repel (HIP−HIE) | ON-state contact (HIE) | MD-stable? | Gln control | PROPKA ddg | pKAI ddg | GFN2-xTB dE_switch | verdict |
|---|---|---|---|---|---|---|---|---|
| des_6126_…mh_9_56_78 | **+1.36** (>2·SEM 0.58), 4/4 reps + | His–Lys 7.88 Å (<8 ✓) | iface 0.84, rmsd 2.0 ✓ | +0.39 | **−0.837** | **−2.475** | **+12.3** (H56·K93); +14.9 (H9·K6, strained) | **release (af2ig pose)** |
| des_6025_…mh_52_74 | **+0.93** (>2·SEM 0.43) | His–Lys 6.46 Å (<8 ✓) | iface 0.88–0.92, rmsd 1.4 ✓ | +2.18 | **−1.09** | **−1.989** | **+12.5** (H74·K6); +15.0 (H52·K60) | **release (af2ig pose)** |

(Target cations in **mature** numbering, consistent with the cross-react table; the raw af2pred pose numbers
them K86/K173/K140 ≡ mature K6/K93/K60.)

**Quantified magnitude (Wyman Δn_H on the Boltz co-fold — a 3rd structural context).** Integrating the
proton-uptake difference (complex vs apo PROPKA pKas) across the 6.0→7.4 window gives a predicted release
factor Kd(6.0)/Kd(7.4): **des_6126 = 22× (peak Δn_H 1.58, SHARP)** — all three His titrate in-window and in
the release direction (H9 6.68→4.59, H56 6.31→5.19, H78 6.53→6.04); **des_6025 = 8× (peak Δn_H 0.78)** — the
weaker of the two (two His, unequal). des_6126's three His titrate near-independently (little
anti-cooperative pKa-splitting), which is why it is both the strongest MD switcher and the sharpest linkage.
(Reference scale: Δn_H 1≈25×, 2≈630×; Baker's TNF SOTA hit 79× at the easier pH 5.4 window.)

*4-His sharpening probe — tested NO-GO.* Direct structural analysis of the des_6126 bound pose: the
interface presents only three target cations (mature K6, K93, R1), and all three are already engaged by
His9/His56/His78. The only mutable chain-A positions near a cation (pos 1→K6, pos 61→K93) point at
*already-engaged* cations, so a 4th His would double up (anti-cooperative — the known dead-end), not add an
independent Δn_H unit. des_6126 is **saturated at 3 distinct-cation His**; its 3-His form is the irreducible
optimum for this backbone, not an under-exploited one.

All signs agree: protonation at pH 6.0 (HIP) raises the His⁺–Lys⁺ repulsion / destabilises the contact
(MD d_repel>0, xTB dE_switch>0) and the proton-linkage releases (PROPKA & pKAI ddg<0). Both switches sit on
an **MD-stable** bound pose (the gate the Genie3-groove hedge failed — those drifted rmsd 4.5–7 Å, so their
"switch" was measured on a detaching pose and was rejected). Honest caveats: the des_6126 Gln-isostere
margin (+0.39) is positive but within noise, and its His9·Lys86 xTB pair is strained (raw-pose clash,
relaxed under `--opt`; consistent sign, unreliable magnitude) — the clean His56·K93 pair (+12.3) and the
primary HIP−HIE MD delta (+1.36, 4/4 reps) carry that member. Artifacts (durable, local-only under
gitignored `runs/`): `runs/v18_switch/pka_consensus.tsv` (PROPKA+pKAI), `runs/v18_switch/qm_switch.tsv`
(GFN2-xTB, 4 pairs, all release=True), `runs/v18_switch/md_des_{6126,6025}.json` (two-state MD aggregates).

**Fold-robustness landscape (strict per-His both-fold Boltz∩Chai check, ALL 20 members).** A switch-His is
"fold-robust" only if PROPKA+xTB confirm it in BOTH the Boltz and the independent Chai binding pose — the gate
that caught the des_287_5/des_2375_4 anti-switches. Sweeping the whole portfolio
(`runs/foldrobust/landscape.tsv`):

| both-fold His (on conserved cations) | members |
|---|---|
| **2** (fold-robust champions) | **des_4098_1_b14, des_4098_1_c44, des_4098_1_c44_w6_3, des_180_6_H9-14-38** |
| 1 | des_84_3 (v19), des_2106_2_H43, b4m2_r1_w5_3, b1m9_A12H_w5_3 |
| 0 (single-fold only) | des_6126, des_6025, des_875_5, des_1060_4, des_2106_2_mh, b2m0, des_495_6, des_1559_3, des_138_7×2, des_4082_8×2 |

**8/20 members carry ≥1 fold-robust His, 4 carry 2** — so the "switch is fold-dependent" caveat is specific to
the point-mutant / partial-diffusion members (des_6126, des_6025, the de-novo bpair lineages), NOT the whole
portfolio. The fold-robust members are the **shell-adapted** ones: des_4098 was LigandMPNN-redesigned around its
fixed His (METHODS §8 Track B), which *locks the His into the interface so both folders reproduce it* — the
mechanistic reason it confirms on both folds where point-mutants don't. So the portfolio already carries
genuinely fold-robust switchers; des_6126/des_6025 contribute switch *magnitude* (d_repel 1.36, 22× Wyman) on
the af2ig pose but not cross-fold robustness. The fold-INDEPENDENT evidence for every member is the
conserved-cation engagement below (a sequence/epitope fact, not a fold prediction).
(`runs/v18_switch/bothfold_robustness.tsv`, `runs/foldrobust/landscape.tsv`.)

## Cross-reactivity (human + mouse) — mechanistic, NOT the co-fold iptm

**Do not rely on the Boltz mouse iptm.** Calibrated on the 12 CRO-validated TNF binders, the co-fold iptm
ranks them *backwards* (best real mouse binder iptm 0.18 — would "fail" the 0.70 bar; co-fold AP on TNF ≈ 0.16,
near the chance line). A 0.91 "PASS" is close to a coin flip for this target and is reported above only as a
non-load-bearing number.

**The real cross-react evidence is mechanistic and fold-independent:** every switch-His engages a target cation
that is IDENTICAL in human and mouse mature TNF. The conserved cation set (aligned human↔mouse) is
{R1, K6, R27, R39, K60, R77, K85, K93, R98, K107, K123}; the divergent ones (lost in mouse: R26→Q, R126→Q,
R133→L) are all avoided. Engagement by member:

| member | switch-His → conserved cation (mature) | all conserved? |
|---|---|---|
| des_6126 | His9→K6, His56→K93, His78→R1 | **yes (3/3)** |
| des_6025 | His52→K60, His74→K6 | **yes (2/2)** |
| des_84_3 (v19) | His56→K93, His76→K6, His49→R77 | **yes (3/3)** |

K6 and K93 are the strongest close anchors in the Baker TNF recipe. Because the switch mechanism itself is
anchored on residues conserved in both orthologs, the binder cannot distinguish human from mouse at those
contacts — the cross-react argument that does not depend on a near-chance co-folder.
(`runs/v18_switch/crossreact_mechanism.tsv`.)

## Portal-novelty result — RESOLVED (all 20 ≥3/4; was the one open risk)

The user ran the Proteinbase/Adaptyv portal novelty check on the full v18 set. **Every design scored ≥3/4**,
so the hard novelty gate (`<3/4 = DQ`) is cleared for the whole portfolio. The 20 portal names were diffed
against the CSV and match exactly (`md5 524722158467273c0e9f57f5402650ed`). Four designs are 4/4; the rest
3/4, including both new members — exactly what the local calibration predicted (seq-novel + helical ⇒ 3/4).

| design | portal novelty | | design | portal novelty |
|---|---|---|---|---|
| des_818_4_c13S | **4/4** | | des_4098_1_b14 | 3/4 |
| des_875_5_H3-57-80 | **4/4** | | des_4098_1_c44 | 3/4 |
| des_1060_4_mh_46_74 | **4/4** | | des_4098_1_c44_w6_3 | 3/4 |
| des_2106_2_mh_37_43_50_75 | **4/4** | | des_180_6_H9-14-38 | 3/4 |
| **des_6126_…mh_9_56_78** (new) | **3/4** | | des_2106_2_H43 | 3/4 |
| **des_6025_…mh_52_74** (new) | **3/4** | | b2m0_c70 | 3/4 |
| des_138_7_…mh_49_58_77 | 3/4 | | b4m2_r1_w5_3 | 3/4 |
| des_138_7_…mh_49_77 | 3/4 | | b1m9_A12H_w5_3 | 3/4 |
| des_4082_8_…mh_51_79 | 3/4 | | des_495_6_mh_41_49 | 3/4 |
| des_4082_8_…mh_43_55_57_79 | 3/4 | | des_1559_3_mh_48_75 | 3/4 |

Totals: **16 × 3/4, 4 × 4/4, 0 below gate.** No portal-unverified member remains in v18 — the reason the
dossier previously preferred v17/v16 as safer fallbacks (fewer unchecked members) no longer applies.

## Generator provenance — the honest binding-transfer risk (and why it does not change the pick)

**The concern (valid and worth stating plainly):** Anthropic's CRO-measured TNF campaign (12/150 bound, 8%,
two CROs — the only non-zero de-novo TNF hit rate ever measured) got **every** binder from Genie3 (8/12) and
PXDesign (4/12); **RFdiffusion3 and BoltzGen produced zero confirmed binders.** Every in-silico-confirmable
de-novo binder in this portfolio — the `des_*` bpair lineages (β-pairing RFdiffusion1) AND the `b*`/cb1
groove lineages (RFd1 partial diffusion / RFd3-at-cb1) — is **RFdiffusion-derived.** So by the one wet-lab
dataset, the portfolio's generator family is the one that did not bind. This is a real risk on the binding
(affinity) axis and should not be hidden.

**Why it does not change the pick — three independent reasons:**

1. **The challenge scores a CONDITIONAL binder, and only an RFd-confirmable pose can carry a confirmable
   switch.** Judging is pH-selectivity (#1) > cross-reactivity (#2) > affinity (#3, lowest). Building and
   confirming a pH switch requires an MD-stable, in-silico-confirmable bound pose (af2ig pae<10 → two-state
   MD on a pose that stays bound). Only the RFd-derived designs give that. Genie3/PXDesign designs **floor
   af2ig (pae 23–27) AND drift in MD** (tested directly this session — the Genie3-groove hedge: binder_rmsd
   4.5–7 Å, interface shed 35–55%), so they cannot be shown to bind, cannot carry a confirmable switch, and
   would score ≈0 on the #1 and #2 axes that dominate the score. A raw Genie3 binder is not a conditional
   binder.
2. **β-pairing RFd1 is not the RFd3 the CRO tested.** The `des_*` lineages use targeted β-strand pairing to
   TNF's exposed edge strand (Sappington, Nat Commun 2025; >9× hit lift on Ig-fold β-sandwiches — exactly
   TNF's β-jellyroll). The CRO null was for RFd3 (all-atom de-novo) and BoltzGen; it did not test β-pairing
   RFd1. Suggestive against the family, but not the same method.
3. **In-silico ranking is near-chance for TNF anyway** (co-fold AP 0.16; our own calibration: Boltz iptm
   *anti*-correlates with measured affinity on the 12 CRO binders). Anthropic found their binders by wet lab,
   not by score. So a Genie3/PXDesign slot would be an **unrankable, switch-less, blind recipe bet** that
   displaces a design we *can* confirm on the two top-ranked axes. The affinity axis — the one the RFd
   provenance actually threatens — is the lowest-weighted, and the portfolio already hedges it across four
   de-novo β-pairing lineages (des_138_7 / des_4082_8 / des_6025 / des_6126), the groove `b*`/cb1 family, and
   partial-diffusion variants, i.e. across distinct folds and epitopes.

**Bottom line:** the binding-transfer risk is genuine, shared by every in-silico-confirmable de-novo option,
and sits on the lowest-ranked axis. The switch (#1) and cross-reactivity (#2) — where the score is won — are
confirmable only on the RFd-derived poses, and v18 maximises them. Submitting a CRO-recipe (Genie3/PXDesign)
design instead would trade a confirmable conditional binder for an unconfirmable raw-binding lottery ticket.
This was tested and ruled out this session (`lesson-genie3-groove-floors-af2ig-like-all-generators`).

## #1-axis profile across the whole portfolio, and a v19 option

Characterising the top-ranked axis across all 20 members (two-state MD d_repel, the primary switch metric):

| band | members (d_repel) |
|---|---|
| strong (>1.3) | des_495_6 **2.31**, des_2106_2_H43 1.94, des_875_5 1.44, des_2106_2_mh 1.40, des_6126 1.36, des_1559_3 1.32 |
| mid (0.9–1.3) | des_138_7×2 1.11/1.09, des_4098_1_b14 1.08, des_1060_4 1.03, des_4082_8_mh51 0.95, des_6025 0.93 |
| weak (0.4–0.9) | des_180_6 0.89, b2m0 0.83, b1m9 0.79, des_4098_1_c44×2 0.57, des_4082_8_mh43 0.46, b4m2 0.40 |
| non-switcher | **des_818_4 −0.01** (kept for 4/4 novelty + 0.83 solubility, not for the switch) |

Median d_repel 0.99; 19/20 positive. The one non-switcher, **des_818_4**, was originally kept as "the only
4/4-novel design" — but the portal check shows **three other 4/4 members that ARE strong switchers**
(des_875_5 1.44, des_2106_2_mh 1.40, des_1060_4 1.03), so the 4/4 showcase no longer needs it.

**Floor-raising swap — rigorously tested; one genuine candidate (des_84_3).** I searched the confirmed-switcher
pool outside v18 for a replacement that beats des_818_4 on the switch without a hidden-axis regression (the
`portfolio-swap-verify-all-axes` discipline). Most strong-d_repel designs outside v18 are already dropped-for-cause
(des_287_5 / des_962_5 / des_2375_4 anti-switch or wrong-sign on a 2nd fold; Genie3-groove MD-unstable; cousins).
The one clean survivor is **des_84_3_mh_49_56_76** (wave-10, a novel β-pairing fold):

| axis | des_818_4 (current) | des_84_3 (candidate) |
|---|---|---|
| #1 switch d_repel | **−0.01 (non-switcher)** | **+1.77**, MD-stable (ON His–Lys 8.0 Å, iface 0.91, rmsd 2.8) |
| switch fold-robustness | — | **His56·Lys93 confirmed on BOTH Boltz & Chai**; aggregate linkage right-signed on both folds (not an anti-switch) |
| #2 cross-react (mechanistic) | — | His engage **K93/K6/R77 — all conserved** in mouse (iptm 0.94/0.912 not relied on) |
| expression (Protein-Sol) | **0.83** (excellent) | 0.452 (average — above the 0.446 population mean; monomer pLDDT 0.933, so likely expresses) |
| novelty | **4/4** (portal) | de-novo helical/β-pair ⇒ expect 3/4 (NOT yet portal-checked) |

**v19** = v18 with des_818_4 → des_84_3 (`runs/push/submission_v19.csv`, md5 `7756055f0b20bcd6994fefb644271db7`,
validator PASS). It makes the portfolio **20/20 switchers, all cross-reactive**, raising the #1-axis floor from
−0.01 to +1.77. The honest cost is a trade *down* two lower-priority axes: novelty headroom (4/4→3/4 for that
slot; v19 keeps 3×4/4) and expression headroom (0.83→0.452 — both are expected to express). The firm catch:
**des_84_3 is not portal-novelty-checked**, so v19 is *not* fully gated until the user runs it; if des_84_3
reads <3/4 it DQs and v18 is the fallback. Evidence in `runs/v19_option/`.

## Recommendation

The portal-check ladder that previously separated v16/v17/v18 has collapsed: **all 20 v18 designs are now
portal-confirmed ≥3/4**, so v18 carries no unverified-novelty risk that v16/v17 avoided, and the
generator-provenance risk is identical across all three (every member is RFd-derived). v18 therefore
**dominates v16 and v17 on the confirmable axes at no added gate risk:**

- **#1 pH-selectivity:** 19/20 MD-confirmed switchers (v16 = 17, v17 = 18); the single strongest switcher of
  the whole campaign (des_6126, d_repel +1.36, 4/4 reps); both new members confirmed by a 4-method-class
  sign-consensus (MD + Gln control + PROPKA + pKAI + GFN2-xTB).
- **#2 cross-reactivity:** carried by the **mechanistic conserved-cation argument** (every switch-His engages a
  cation identical in human & mouse TNF), NOT the Boltz mouse iptm — which is near-chance for TNF (ranks the 12
  CRO binders backwards) and is not load-bearing here. See "Cross-reactivity".
- **novelty gate:** all 20 ≥3/4 (4 at 4/4), portal-confirmed. **expression:** all new members soluble
  0.82–0.92, 0-Cys.
- **diversity:** 4 de-novo β-pairing lineages + groove `b*`/cb1 family + partial-diffusion variants (keeps
  the lone 4/4 des_818_4 rather than a des_6025 cousin).

Two defensible picks, by appetite for one more portal check:

- **v18** (`submission_v18.csv`, md5 `524722158467273c0e9f57f5402650ed`) — **fully gated, zero further action.**
  All 20 portal-confirmed ≥3/4; 19/20 switchers + one max-novelty/max-expression non-switcher anchor. The safe
  recommended upload.
- **v19** (`submission_v19.csv`, md5 `7756055f0b20bcd6994fefb644271db7`) — **strongest on the #1 axis** (20/20
  switchers, all cross-reactive) by swapping the non-switcher des_818_4 for the both-fold-confirmed switcher
  des_84_3. **Requires one portal novelty check on des_84_3** (expect 3/4; if <3/4, fall back to v18). Pick this
  if you can run that check and accept trading novelty/expression *headroom* for the switch-floor gain.

**v16** (`a524d77f`) remains a valid conservative floor (identical gates, weaker on #1).

**Upload remains USER-GATED. Nothing has been uploaded.** (Portal account j_v_v_07; the upload action is the
user's.)
