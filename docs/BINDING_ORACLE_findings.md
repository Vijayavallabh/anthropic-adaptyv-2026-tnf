# Binding / cross-react in-silico oracle — attempts and findings (TNF)

The challenge's #2 (cross-react) and #3 (affinity) axes need a way to tell, in-silico, which designs will
actually bind (and bind mouse). This log records what was tried and why TNF resists it. **Bottom line: every
co-fold-based binding oracle is near-chance for TNF; the only robust in-silico evidence is physics (MD pH-switch
on a pose) + mechanism (switch-His engaging human↔mouse-conserved cations).**

## Established before this session
- **Single co-folder scores (Boltz/Chai iptm, ipSAE) are near-chance for TNF.** Co-fold AP ≈ 0.16 (near the
  12/150 chance line); on the 12 CRO-validated binders, Boltz iptm *anti*-correlates with measured affinity
  (best mouse binder iptm 0.18). Boltz-msa-empty recovers only 3/12 (misses all 8 Genie3-L83, incl the 3
  mouse-cross-reactive). (`reference-anthropic-tnf-dataset`, `lesson-chai-mouse-iptm-is-folder-pessimism`.)

## This session — cross-folder POSE-CONSENSUS (a new oracle, not a score ensemble)
**Hypothesis:** even if each folder's *score* is near-chance, if two independent seq-only co-folders (Boltz-2,
Chai-1) *dock a binder to the same epitope in the same orientation*, that pose AGREEMENT could flag a real
interface. Metric (`scripts/pose_consensus.py`): superpose the TNF target across the two complexes, then
binder-Cα-RMSD (low = same placement) and interface-residue Jaccard (high = same epitope).

**Calibration** (`runs/binding_oracle/pose_consensus_calibration.tsv`): 12 CRO-validated binders + 12
composition-matched scrambles, folded with Boltz+Chai vs human and mouse.

| test | metric | AUC | verdict |
|---|---|---|---|
| binding (real vs scramble, human) | iface Jaccard (real>scram) | **0.389** | near-chance (scrambles agree *more*) |
| binding (real vs scramble, human) | binder-RMSD (real<scram) | **0.549** | near-chance |
| cross-react (mouse-XR vs human-spec, mouse) | iface Jaccard | **0.667** (n=3) | weak, underpowered, not significant |

Real binders show Boltz↔Chai binder-RMSD **~35 Å** (docked to *different faces* of TNF) — indistinguishable
from scrambles. **NO-GO: pose-consensus cannot confirm TNF binding.** The two folders don't even agree on which
face a *known real binder* uses. (The mouse-XR signal 0.667 is suggestive — the 3 real mouse-cross-reactive
binders pose-converge slightly better on mouse — but n=3 is far too small to rely on.)

## Consequence for the submission
Three independent lines now agree that TNF binding/cross-react is **in-silico-walled** (single-folder scores,
strict both-fold switch, and pose-consensus all near-chance). Therefore:
- The **cross-react (#2) claim rests on the mechanistic conserved-cation argument**, not any co-fold score
  (see `SUBMISSION_v18_candidate.md` §Cross-reactivity). This is fold-independent and the right basis.
- The **switch (#1) is the strongest axis** because it is physics (two-state MD + pKa + QM on a pose), and its
  fold-INDEPENDENT core is the conserved-cation chemistry.
- Binding/affinity (#3, lowest-weighted) is genuinely wet-lab-decided for TNF — as Anthropic found (they
  selected by an ipSAE ensemble that still only reached 8% CRO hit rate). No in-silico method available here
  changes that; claiming in-silico binding *dominance* would be unsupported.

## A binding oracle is MOOT for the final submission (2026-10-13)

An untested angle was considered: calibrate a binding oracle on **real PDB TNF-α complexes with published Kd**
(already crystallographically posed — a continuous Kd gradient, no fold-generation noise), rather than on the
de-novo CRO folds above (binary labels, ~35 Å Boltz↔Chai pose disagreement). Two reasons make it moot before
any compute is spent, independent of whether it could calibrate:

1. **No cut to make.** Track 3 caps at 20 designs; the submission (v24) is *exactly 20*. A ranker only adds
   value when candidates must be dropped to fit a cap — here there is nothing to drop.
2. **#3 is the lowest-priority axis.** Any swap justified by a better #3 score would eject a member selected
   for its #1 (switch direction/depth), #2 (cross-react) or novelty strength — trading a higher-ranked axis
   for a lower one, which is churn, not improvement (`lesson-portfolio-swap-verify-all-axes-same-pipeline`).

Also, a PDB-Kd set is dominated by antibody/receptor interfaces (large CDR/loop contacts); transferring such
an oracle to small Region-I minibinders is an out-of-regime extrapolation (`lesson-validate-in-target-regime`,
`lesson-powered-mechanism-can-still-be-nogo-in-target-regime`). So the oracle is both out-of-regime AND moot.

**Consequence:** the correct in-silico treatment of #3 is a **floor** (every member passes af2ig: pae_int
4.1–6.7, plddt_binder 85–94, rmsd < 1 Å) — which v24 has — not a rank. This is a provably-correct design
choice, not a concession. Binding *magnitude* dominance is wet-lab-decided and cannot be manufactured in-silico.

## Cross-reactivity (#2) footprint audit — all 20 members clean (2026-10-13)

The provable sub-metric for #2 is the fraction of the binding footprint on human↔mouse-**identical**
residues (Anthropic's x-react determinant: cross-reactivity tracks # substituted footprint residues). All 20
v24 members were re-folded (binder + human TNF, Boltz-2) and every contacted TNF residue (<5 Å) classified
against the mouse ortholog. Result: **mean footprint conservation 93.2 % (min 90 %); zero genuine
cross-reactivity liabilities** — i.e. no member has a non-conservative, deep (≥6 binder atoms <5 Å),
non-N-terminal divergent core contact. The entire ~7 % divergent residual is the ragged flexible N-terminus
(`T1>N`, `P2>S` — human `RTPSDK…` vs the mouse N-terminal extension; an alignment/flex artifact), conservative
substitutions (`I→V`, `E→D`), or peripheral single-atom brushes. Separately, 100 % of switch-His engage the
conserved cation set {1,6,27,39,60,77,85,93,98,107,123}, none the divergent set {26,126,133}. So #2 is maxed on
the provable core across the whole portfolio and **no #2-motivated swap exists**. (Caveat: single-protomer
msa-empty co-fold — TNF is a trimer — so footprints are approximate, but the residue-identity conclusion is
robust to pose.)
