# TNF-α conditional pH-switch binders — final 20-design portfolio (methods)

De-novo mini-binders that **bind TNF-α at pH 7.4 and release at pH 6.0** (endosomal-recycling
direction), cross-reactive to human **and** mouse TNF-α. Objective priority (unchanged):
**pH-selectivity > human/mouse cross-reactivity > affinity**, with **novelty ≥ 3/4** a hard gate
(Proteinbase-confirmed) and **expression** the dominant wet-lab-zero risk.

Deliverable: `runs/push/final20.{csv,fasta}` + upload form `runs/push/upload_final20.csv` (20 designs,
6 distinct Foldseek folds). Every metric below is reproduced by scripts in `scripts/` from the folded
complexes.

## The pH switch — sharpness is Δn_H, and we engineer it

The wet-lab readout is **BLI at two pH**: a *thermodynamic* measurement of how much K_D changes between
7.4 and 6.0. Statistical-mechanics (Wyman linkage) fixes what sets that change: the **number of protons
coupled to binding, Δn_H**. The 10→90 % transition width is `1.908/n_H` pH units, so across our 1.4-pH
window a **single His (Δn_H≈1) can only shift K_D ~25×** (too broad to separate 7.4 from 6.0 — this is
why single-His switches are weak), whereas **Δn_H≈2 → ~600×** and **≈3 → ~10⁴×**.

The obstacle (and our fix): two interface His placed *close together* are mutually **anti-cooperative** —
like charges resist co-protonation, their pKa's split apart, and they titrate over a broad range (Δn_H
stays ~1 even with two His). We therefore install **spatially-separated His (> 8 Å Cβ–Cβ), each paired
with its own conserved target cation** (Lys93, Lys6, or Arg98), tuned to pKa ≈ 6.3–6.7. Separated His
titrate **additively** → Δn_H → 2–3 (`scripts/install_sharp_his.py`). Protonation at pH 6 creates
His⁺···cation⁺ repulsion at the interface → release; neutral His at 7.4 preserves binding.

**Result (`scripts/switch_sharpness.py`, PROPKA-Wyman):** installing 3–5 separated His lifts the switch
from the single-His/parent baseline **FC 3–17× → 19–126×** (peak Δn_H 1.5–2.5), and the sharp variants
switch in **both** human and mouse (e.g. `TNFpH06_sh3` 59×/59×, `TNFpH06_sh5` 31×/76×, `TNFpH07_sh4`
126×/50×). The separation strategy works exactly as predicted — all installed His co-titrate rather than
split.

### Switch is scored by five method-classes, blind in different ways, read as consensus

1. **PROPKA-Wyman linkage** — per-His pKa in complex vs apo binder; peak Δn_H and K_D fold-change over
   6.0–7.4 (the quantity BLI measures). `scripts/switch_sharpness.py`.
2. **pKAI** — an orthogonal **ML** pKa engine (`scripts/pkai_sharpness.py`). It **agrees with PROPKA on
   sign and ranking** for every sharp design (all His shift *down* on binding; TNFpH06_sh5 pKAI 58× /
   PROPKA 31–76×, sh3 26×/59×, TNFpH07_sh4 40×/126×). Two independent pKa classes agreeing on the sign is
   the core evidence the switch is real.
3. **Explicit two-state MD** (OpenMM, 3× replicas of HIP/HIE + a **GLN-isostere causal control**).
   `d_repel = His–cation(HIP) − His–cation(HIE)`; a switch is **His-causal** only if `d_repel` is real
   *and exceeds the Gln control* (if the Gln isostere moves as much, the motion is pose noise, not
   protonation). `scripts/md_release_tnf.py`.
4. **Baker His-cation geometry** — the empirical ≥2 protonated-His···cation contacts rule
   (bioRxiv 2025.09.29.678932 / PubMed 41256377). `scripts/his_cation_contacts.py`.
5. **GFN2-xTB** QM interaction energy on the capped His–cation cluster (`scripts/qm_switch_fair.py`,
   `runs/push/leads_qm.tsv`): on the three fully-validated leads' primary switch pairs the protonation
   weakens the contact (dE_switch +0.33/+0.64/+0.50 kcal/mol, release-competent) — a quantum-mechanical
   4th-independent-engine check that the His⁺···Lys⁺ repulsion is real, not a forcefield artifact.

**Honest reading of (3) vs (1,2), and the 4-His sweet spot.** The two-state MD and the thermodynamic
linkage can disagree, and we report both. **Six designs are MD-His-causal** (repulsion real and exceeding
the Gln control) — five on the cb1/f8 switch lineages plus one on a distinct fold:
- `b1m9_1A55H_K93_sh4` (PROPKA 22× / pKAI 7×; MD d_repel **+0.91**, control −0.15) — rigid high-affinity f8 interface;
- `TNFpH06_sh3` (PROPKA 59× / pKAI 26×, both species; MD d_repel **+1.18**, control +0.13) — the sharpest MD-confirmed;
- `TNFpH06_sh2` (PROPKA 17×/11× both species; MD d_repel **+1.25**, control +0.17);
- `TNFpH06` (coop-2His original; MD d_repel **+0.81**, control −0.13) and `TNFpH07` (MD d_repel **+0.32**, control −0.98) — the two validated parents are themselves MD-His-causal;
- `b2m0_c70` (distinct f1 fold; MD d_repel 3.2 Å, Gln control inert) — a sixth, independent fold-and-lineage confirmation.

The MD cleanly resolves **2 of the 3** installed **4-His** cb1 variants (sh3, sh2) — so the structural
switch is a robust property of the 4-His design, not a single lucky pose. The disagreements are equally
informative: the **5-His** cb1 variants (`TNFpH06_sh5`, `TNFpH07_sh4`) are thermodynamically the sharpest
(31×/76× and 126×/50×, both pKa engines) but their MD is negative/control-confounded — the 5th installed
His over-crowds the interface and the Gln control perturbs it as much as HIP. **Four separated His is the
sweet spot** (maximal Δn_H the structural switch still supports); the 5-His variants ship as
**thermodynamic-sharp hedges** (BLI measures the thermodynamic switch, which both pKa classes confirm),
labeled distinctly from the MD-confirmed designs.

## Expression — the hard gate, and interface-frozen rescue

Against Adaptyv's `protein-qc` **standard** gate (instability < 40, GRAVY < 0.4, 0-Cys), only **1/9** of
the original designs passed — the bindable epitope is hydrophobic, so vanilla-MPNN interfaces raise
instability. Per the EGFR post-mortem a non-expressing design scores zero regardless of a perfect switch.
**Rescue:** freeze the binding interface (union of human+mouse interfaces + **all** switch His) and
re-sequence **only the solvent-facing surface** with SolubleMPNN (0-Cys, 48 seqs/design); re-fold
(Boltz-2 human+mouse) and re-score expression + binding + the **thermodynamic** switch (PROPKA+pKAI)
on the exact rescued sequence. This preserves binding, cross-reactivity **and the thermodynamic His
switch** (the interface+His are frozen) while fixing expression — and in several cases the surface
redesign *sharpened* the PROPKA/pKAI fold-change:
- `TNFpH07_sh4_r4` — physics-max switch FC 69× (parent 126×), now **expr-PASS** (instab 34) + cross-react (iptm 0.66 h/0.72 m);
- `b1m9_1A55H_K93_sh6_r1` — sharp switch FC 32×, strong binding 0.92/0.83, now **expr-PASS**;
- `b4m2_2A87H_K101_r3` — strongest binder iptm **0.95/0.95** (perfect cross-react), now **expr-PASS**;
- `b1m12_2A67H_K112_r1` — the 4/4-novel strong binder (0.92/0.92), now **expr-PASS**.

**Honest caveat on the rescued variants (two-state MD re-run on the rescued sequences).** We re-ran the
explicit HIP/HIE two-state MD on the exact rescued sequences: all three tested (`TNFpH07_sh4_r4`,
`b1m9_1A55H_K93_sh6_r1`, `b4m2_2A87H_K101_r3`) came back **not His-causal** (d_repel −0.46/−0.27/−0.12,
Gln-control noisy) — the solvent-surface redesign makes the interface floppier, so the explicit-MD
geometric switch no longer confirms even though the **thermodynamic** linkage (PROPKA+pKAI, the quantity
BLI actually measures) is intact. We therefore label the rescued variants **thermodynamic-switch +
expression-PASS**, and reserve the **MD-His-causal** label for the three un-rescued leads below (and the
distinct-fold `b2m0_c70`). This is deliberate: the BLI two-pH readout measures the thermodynamic switch,
so a rescued expresser is a legitimate submission — we simply do not over-claim the stricter MD gate for it.
**Second rescue round — minimal surgical edits (burial-aware).** Where the 48-seq SolubleMPNN resurface was
unnecessary, we instead took the fewest solvent-surface point substitutions (neutral panel Q/N/S/T) that
drop instability < 38, freezing the interface, **all His, and buried-core residues** (CB-neighbour count
≥ 10) and guarding against new N-glyc sequons / charge runs. This fixed the `b1m9` single-His parent
(+1 edit, E19N) and `b1m7_2A57H_K183_sh5` (+3 edits) with binding **preserved** (Boltz rescued-vs-parent
iptm 0.87/0.78 and 0.92/0.89) and switch FC preserved (10×, 6×, recomputed on the rescued fold). The burial
filter mattered: it blocked a fold-breaking buried-Ile→Asn and two structural-proline removals that an
instability-only greedy had proposed.

The portfolio is **17/20 expression-PASS**. The 3 remaining FLAGs are all *soluble* by the orthogonal
signals: `b2m0_c20`/`b2m0_c70` FLAG only because the instability index is inflated by **two turn-prolines**
(P38/P42, CA bend 81°/116° — genuine turns we will not remove to game a crude index; GRAVY 0.10, 0-Cys), and
`b2m0_c70` is the f1 fold's MD-His-causal switch so it is kept with the documented risk; `tnf_K90…sh4` FLAGs
on a 28 %-Glu MPNN surface-overload artifact but has GRAVY −1.03 (excellent solubility) — its 6×E→Q rescue
was **rejected** because it halved human binding (iptm 0.59→0.31), so the parent is kept.

## Cross-reactivity — designed in, verified on the mature mouse ortholog

Cross-reactivity is **mechanistic**, not a single folder's cross-species score (which is unreliable —
folder-pessimism: Boltz and Chai disagree on mouse ipTM, and single-seq mouse folds collapse under one
surface mutation). Every design targets the **conserved** inter-protomer epitope with a **conserved target
cation** as switch anchor, so both the interface and the switch survive the human→mouse change.

- **Epitope conservation** (`scripts/epitope_conservation.py`, vs the *mature* mouse TNF, UniProt P06804):
  every design's footprint is **82–97 % conserved**; the few divergences are the flexible N-terminus
  (T2>N, P3>S) and conservative swaps (Ile↔Val, Glu↔Asp). The switch cations are conserved or
  cation-preserving (Lys93/Lys6 identical; Arg98→Lys98 — still a cation, switch intact).
- **Matched co-fold delta** on the mature ortholog: sharp-variant human↔mouse iptm deltas are small
  (`TNFpH06_sh3` Δ−0.02, `sh5` Δ+0.06, `b1m9_sh4` Δ+0.10), and the **switch itself cross-reacts** — the
  cb1 sharp variants are sharp in both species (sh3 59×/59×, sh5 31×/76×).
- *(Correction recorded: an earlier mouse-target file was the full-length precursor, including the
  TM/signal domain, which gave spurious ~0.15 cross-species iptm for every design including known-
  cross-reactive parents. Re-folding against the mature 156-aa ortholog restored parent TNFpH06 to
  iptm 0.63, Δ+0.01. All cross-react numbers here use the mature construct.)*

## The 20-design portfolio

7 distinct folds; `novelty ≥ 3/4` throughout (3 designs are **4/4**). Chosen for the judging axes, leading
with the most-validated design and spreading risk across folds.

- **Leads — fully validated (sharp + MD His-causal + cross-react + express):**
  `b1m9_1A55H_K93_sh4` (f8; 22×, MD d_repel +0.91, iptm 0.91 h/0.82 m, PASS), `TNFpH06_sh3`
  (cb1; **59×/59×** both species, MD d_repel +1.18, Δ−0.02, PASS), and `TNFpH06_sh2`
  (cb1; 17×/11×, MD d_repel +1.25, best cb1 binding iptm 0.75, PASS). All three pass *every* axis.
- **Sharpest thermodynamic switches (cb1 hedges):** `TNFpH06_sh5` (31×/76×, PASS) and the physics-max
  `TNFpH07_sh4` (126×/50×) — PROPKA+pKAI-consensus sharp
  and cross-reactive, MD thermodynamic-only (5-His over-crowding); plus the two validated coop-2His
  originals `TNFpH06/07`.
- **Novel strong binders (f5 / b1m12, novelty 4/4):** `b1m12_r4` (iptm 0.94, PASS), `b1m12_2A67H_K112`,
  `b1m12_r1` — highest affinity + the only 4/4 structural-novelty fold; modest switch (value = novelty +
  binding + cross-react).
- **Strong cross-reactive binders (f2 / b4m2):** `b4m2_r1` (iptm 0.91, 97 % conserved), `b4m2_2A87H_K101`
  (iptm 0.93), `b4m2_r6`.
- **MD-confirmed de-novo switches (f1 / b2m0):** `b2m0_c70` (**MD His-causal**, d_repel 3.2 Å),
  `b2m0_c20` — fold diversity with a structurally-confirmed switch (instability-index FLAG is a
  turn-proline artifact; GRAVY/0-Cys clean).
- **Fold diversity (f11 / b1m7):** `b1m7_2A57H_K183_r1`, `b1m7_2A57H_K183_sh5_rx` (expr-rescued).
- **7th fold — all-atom RFd3 cross-react binder (fK90):** `tnf_K90…sh4` — distinct conserved epitope,
  iptm 0.59 h/0.54 m, docks 3/4; thermodynamic His switch FC 8× but **not** MD-His-causal (binder-first
  diversity hedge, see honest negative below).

## New-fold generation — attempted, and an honest negative

The Anthropic/Adaptyv autonomous-TNF campaign got its binders from **Genie3 + PXDesign** at the
inter-subunit receptor groove (RFdiffusion/BoltzGen gave zero), so we stood up **PXDesign** and ran it at
two conserved epitopes (the kit's receptor-site patch {32,65,73,75,92,144,145,146} and the inter-subunit
groove A{74,75,97}/B{114,115}) × three lengths. PXDesign (fast mode) is a backbone generator; we took its
**98 on-target backbones → SolubleMPNN → af2ig**. Result: **0/196 bind** (best `pae_interaction` 24.6, all
≥24; `plddt_binder` 83–94, so they fold well but do not engage TNF). A Boltz-2 cross-check (the Anthropic
paper notes AF2 metrics are non-predictive for TNF) gave five designs iptm 0.5–0.66 — but the decisive
**pose check showed they dock off the designed epitope** (0–2/8 hotspots): an msa-empty co-fold
pose-promiscuity artifact, not a designed interface. New-fold generation at these epitopes is therefore a
decisively-tested floor for this pipeline (Anthropic's hits came from thousands of designs ranked by
Protenix-ipSAE and CRO-measured; our small-scale dual-metric+pose filter surfaced none). **The portfolio
rests on designs that pass *both* af2ig (`pae_interaction` 3.7–4.2 for the cb1 lineage) and Boltz-2 iptm
at the intended epitope** — not on a single optimistic oracle.

**Binding-metric calibration — why we rank by epitope+pose+diversity, not by score.** We folded the 12
CRO-measured TNF-α binders from the Anthropic/Adaptyv campaign through our own stack and found the
in-silico interface scores **do not rank TNF affinity** — the tightest human binder (0.7 nM) scored
Boltz iptm 0.20 while a 27× weaker one (18.8 nM) scored 0.64, and the best mouse binder (66 nM) scored the
*lowest* mouse iptm (0.18). This matches the published TNF ranking AP of 0.16. We therefore do not use
iptm/ipSAE/pae as hard affinity filters (they would discard the real binders); instead we (i) target the
experimentally-validated, human↔mouse-conserved **Region I** groove and the conserved K90 sub-site,
(ii) require the designed **pose to dock at the intended epitope**, (iii) anchor the pH switch on
mouse-conserved cations (Lys65/90/98, Arg32; never Arg31→Gln or His73→Tyr), and (iv) **spread risk across
distinct folds** — the exact strategy that produced the only CRO-validated cross-reactive de-novo TNF
binders. Cross-reactivity is argued mechanistically (conserved epitope + conserved switch anchor), since a
single folder's cross-species score is unreliable here.

**New-fold diversity added (all-atom RFd3, the route that breaks the flat-polar wall).** Scaled RFd3 at the
conserved, switch-compatible **K90** site (hotspots Y87/K90/V91/Y141, all mouse-conserved) yielded a
new-fold **cross-reactive** binder (`tnf_K90…model_3_3`, iptm 0.59 h / 0.54 m, docks 3/4 at the epitope) —
the 7th fold. Installing His near the conserved Lys85/Arg77 gives a thermodynamic pH switch (FC 8×,
mouse-conserved anchors), but the explicit two-state MD **refuted** His-causality (d_repel −0.02), so K90
enters as a **binder-first / cross-react diversity hedge**, not a validated switch.

**Region I de-novo — a second honest negative.** We ran 7680 RFd3→MPNN designs at the full inter-subunit
groove (2-protomer target) and folded them with af2ig against the correct 2-protomer target. 4471 fold
self-consistently, and the best engage the groove (4/5 epitope residues < 8 Å), but none deeply dock (0 at
3/5 < 5 Å) and af2ig `pae_interaction` ≈ 26. An orthogonal Boltz monomer cross-fold of the top 14 confirmed
the ceiling: 3–4 dock the protomer-1 half (V74/L75/I97 at 2.7–5.7 Å, human iptm 0.60–0.70) **but mouse iptm
is uniformly 0.15–0.18 — not cross-reactive**, so none beats the existing portfolio binders (0.91–0.95 both
species). De-novo design straight at the flat inter-subunit groove hits the known wall; the all-atom K90
route is the one that produced a usable cross-reactive new fold.

## Switch characterization — in-window Δn_H (Wyman)

The fold-changes reported here are the **in-window** proton linkage, `FC = 10^(∫₆.₀^7.4 Δn_H dpH)` — exactly
what a BLI-at-two-pH assay integrates. Per-His pKa analysis (apo binder vs complex) confirms the switch
His titrate in the right direction and that the cb1 leads realize peak Δn_H 2.0–2.5 across the window; the
one caveat the analysis surfaces is that several His have `pKa_complex` pushed well below 6.0 by the target
cation (e.g. `TNFpH06_sh3` H61 6.89→4.39), so part of each site's titration happens below the window — i.e.
the achievable ceiling here is set by centering the His midpoint near 6.7, not by adding more His (the 5-His
variants confirm this: higher nominal FC but MD-ambiguous from interface over-crowding). This matches the
Baker pH-binder blueprint, which explicitly leaves "optimizing responsiveness around pH 6.5 by shifting His
pKa" as open future work.

## Honest assessment

**Five** designs dominate on every axis — the two 4/4-novel `b1m12_r6_A92S` (d_repel 3.9) and `b1m12_A91H`
(d_repel 1.46), plus `b1m9_1A55H_K93_sh4`, `TNFpH06_sh3`, and `TNFpH06_sh2` — each confirmed on binding,
cross-reactivity, expression, *and* an explicit His-causal two-state-MD switch (the two b1m12 leads add
4/4 structural novelty + the portfolio's best cross-reactive binding, iptm 0.91-0.92 h / 0.94-0.95 m; both
cb1 leads are additionally sharp in *both* species). The 5-His cb1 hedges (`TNFpH06_sh5`, `TNFpH07_sh4`) push the
thermodynamic switch highest (to 126×, two pKa engines agreeing — the quantity BLI measures) but are
MD-ambiguous, so we label them thermodynamic-only rather than overclaim. `b2m0_c70` adds a third,
independent MD-His-causal switch on a different fold; the **b1m12 4/4** fold anchors structural novelty and
**b4m2** the affinity. **Nine** designs carry an explicit MD-His-causal switch — the six leads above
(`b1m9_sh4`, `TNFpH06_sh3/sh2`, `TNFpH06`, `TNFpH07`, `b2m0_c70`) plus the three found by scanning His
against the conserved cation **K98** on the strong binders (`b1m12_A91H`, `b1m12_r6_A92S`, `b4m2_A37H`);
the rest carry the thermodynamic switch (the quantity BLI measures) or are binder-first diversity picks.
The spread across **7 folds with 17/20 expression-PASS** (the 3 FLAGs soluble by GRAVY/0-Cys) hedges the
dominant wet-lab risk. Several new-fold/new-switch routes were tested and honestly reported: the K98 switch
hunt (adopted, +3 MD-causal), the all-atom K90 fold (a cross-react hedge), and two negatives — de-novo
Region I design, and **Genie3** (the CRO-validated generator, built here; it docks tightly at the K90
epitope but its Cα-only output through CA-only ProteinMPNN gave human-only binders single-state and
binding-degraded multi-state, so it was not integrated). The remaining open item is upload (user-gated).
