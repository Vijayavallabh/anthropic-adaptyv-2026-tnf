# Submission v18 candidate — 4th de-novo lineage + strongest switcher (PORTAL-NOVELTY CONFIRMED)

Updated 2026-10-12 (SESSION 11, continued). Deadline **2026-10-11 23:59 AoE = 2026-10-12 11:59 UTC**.

> **2026-10-12 UPDATE — the one open risk is now RESOLVED.** The user ran the portal novelty check on all
> 20 v18 designs: **every design scored ≥3/4 (the hard gate; <3/4 = DQ)**. Both new members read **3/4**
> (`des_6126…mh_9_56_78`, `des_6025…mh_52_74`), exactly as the local calibration predicted; four designs
> read **4/4**. v18 is therefore **fully gated** — no portal-unverified member remains. The portal list was
> byte-checked against the CSV (all 20 names match `md5 524722158467273c0e9f57f5402650ed`). See
> "Portal-novelty result" below; it supersedes the old "one open risk" section.

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
| **#2 cross-react** | human iptm (Boltz-2) | **0.970** | — |
|  | mouse iptm (Boltz-2) | **0.911** | — |
|  | min(h,m) | **0.911** | **PASS** (≥0.70) |
| **expression** | Protein-Sol scaled solubility | 0.82 | PASS (≥0.50) |
| **novelty seq** | mmseqs vs PDB seqres | **NO PDB hit (novel)** | passes seq criterion |
| **novelty struct** | foldseek max-TM vs PDB → **portal** | **0.948** → portal | **3/4 — PORTAL-CONFIRMED** |

## Swap-in 2: des_6025_4_T10_pd_77_4_mh_52_74

`SEETERLIRELNQILLEEMKETVKKLKELGVEVELNVSENGGVSLIASSESHELALEASQLIKETVRRIEQLLHKYVAER` (80 aa, 0 Cys, His 52/74)

| Axis | Metric | Value | Verdict |
|---|---|---|---|
| **#1 pH-switch** | two-state MD d_repel (HIP−HIE), 4 reps × 3 ns | **+0.93** (2·SEM ≈ 0.43, all arms) | **CONFIRMED** |
|  | isostere control (HIP−GLN) | +2.18 | positive |
| **#2 cross-react** | human / mouse iptm (Boltz-2) | 0.950 / 0.936 → **min 0.936** | **PASS** |
| **expression** | Protein-Sol scaled solubility | 0.92 | PASS |
| **novelty** | mmseqs seq / foldseek struct TM → **portal** | NO PDB hit / **0.870** → portal | **3/4 — PORTAL-CONFIRMED** |

Note: this is a **stronger** des_6025 switcher than v17's member (`mh_34_50_52_74`, +0.372) — v18 upgrades the
des_6025 representative while using an independent variant.

## The pH-switch is confirmed by FOUR method-classes, not MD alone

The decisive (#1-ranked) axis is the pH switch. For both new members it agrees in sign across four
independent method-classes — two-state MD (fixed-charge MM), the Gln-isostere ablation control, two
empirical pKa predictors (PROPKA-Wyman, pKAI), and a QM single-point (GFN2-xTB). A positive static score
only proves switch DIRECTION; sign-consensus across orthogonal physics is the honest selectivity evidence
(magnitudes are not comparable across classes and are not claimed).

| member | MD d_repel (HIP−HIE) | ON-state contact (HIE) | MD-stable? | Gln control | PROPKA ddg | pKAI ddg | GFN2-xTB dE_switch | verdict |
|---|---|---|---|---|---|---|---|---|
| des_6126_…mh_9_56_78 | **+1.36** (>2·SEM 0.58), 4/4 reps + | His–Lys 7.88 Å (<8 ✓) | iface 0.84, rmsd 2.0 ✓ | +0.39 | **−0.837** | **−2.475** | **+12.3** (H56·K173); +14.9 (H9·K86, strained) | **release, 4-class sign-consensus** |
| des_6025_…mh_52_74 | **+0.93** (>2·SEM 0.43) | His–Lys 6.46 Å (<8 ✓) | iface 0.88–0.92, rmsd 1.4 ✓ | +2.18 | **−1.09** | **−1.989** | **+12.5** (H74·K86); +15.0 (H52·K140) | **release, 4-class sign-consensus** |

**Quantified magnitude (Wyman Δn_H on the Boltz co-fold — a 3rd structural context).** Integrating the
proton-uptake difference (complex vs apo PROPKA pKas) across the 6.0→7.4 window gives a predicted release
factor Kd(6.0)/Kd(7.4): **des_6126 = 22× (peak Δn_H 1.58, SHARP)** — all three His titrate in-window and in
the release direction (H9 6.68→4.59, H56 6.31→5.19, H78 6.53→6.04); **des_6025 = 8× (peak Δn_H 0.78)** — the
weaker of the two (two His, unequal). des_6126's three His titrate near-independently (little
anti-cooperative pKa-splitting), which is why it is both the strongest MD switcher and the sharpest linkage.
(Reference scale: Δn_H 1≈25×, 2≈630×; Baker's TNF SOTA hit 79× at the easier pH 5.4 window.)

All signs agree: protonation at pH 6.0 (HIP) raises the His⁺–Lys⁺ repulsion / destabilises the contact
(MD d_repel>0, xTB dE_switch>0) and the proton-linkage releases (PROPKA & pKAI ddg<0). Both switches sit on
an **MD-stable** bound pose (the gate the Genie3-groove hedge failed — those drifted rmsd 4.5–7 Å, so their
"switch" was measured on a detaching pose and was rejected). Honest caveats: the des_6126 Gln-isostere
margin (+0.39) is positive but within noise, and its His9·Lys86 xTB pair is strained (raw-pose clash,
relaxed under `--opt`; consistent sign, unreliable magnitude) — the clean His56·Lys173 pair (+12.3) and the
primary HIP−HIE MD delta (+1.36, 4/4 reps) carry that member. Artifacts (durable, local-only under
gitignored `runs/`): `runs/v18_switch/pka_consensus.tsv` (PROPKA+pKAI), `runs/v18_switch/qm_switch.tsv`
(GFN2-xTB, 4 pairs, all release=True), `runs/v18_switch/md_des_{6126,6025}.json` (two-state MD aggregates).

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

## Recommendation

The portal-check ladder that previously separated v16/v17/v18 has collapsed: **all 20 v18 designs are now
portal-confirmed ≥3/4**, so v18 carries no unverified-novelty risk that v16/v17 avoided, and the
generator-provenance risk is identical across all three (every member is RFd-derived). v18 therefore
**dominates v16 and v17 on the confirmable axes at no added gate risk:**

- **#1 pH-selectivity:** 19/20 MD-confirmed switchers (v16 = 17, v17 = 18); the single strongest switcher of
  the whole campaign (des_6126, d_repel +1.36, 4/4 reps); both new members confirmed by a 4-method-class
  sign-consensus (MD + Gln control + PROPKA + pKAI + GFN2-xTB).
- **#2 cross-reactivity:** all members min-iptm 0.91–0.95 (note: co-fold iptm is a weak TNF signal — the
  mechanistic conserved-cation argument carries the real cross-react claim).
- **novelty gate:** all 20 ≥3/4 (4 at 4/4), portal-confirmed. **expression:** all new members soluble
  0.82–0.92, 0-Cys.
- **diversity:** 4 de-novo β-pairing lineages + groove `b*`/cb1 family + partial-diffusion variants (keeps
  the lone 4/4 des_818_4 rather than a des_6025 cousin).

**Recommended upload: v18** (`runs/push/submission_v18.csv`, md5 `524722158467273c0e9f57f5402650ed`).
**v16** (`a524d77f`) remains a valid conservative floor (identical gates, weaker only on the #1 axis) if a
last-minute reason to minimise changes arises; there is no longer a novelty reason to prefer it.

**Upload remains USER-GATED. Nothing has been uploaded.** (Portal account j_v_v_07; the upload action is the
user's.)
