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
