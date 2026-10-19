# Assign design method: field text (v27, 18 proteins)

## Method name

RFdiffusion + ProteinMPNN + AF2; His pH-switch confirmed by alchemical FE

(73 characters, under the portal's 80-char limit. Alternatives that also fit:
"RFdiffusion + ProteinMPNN + AF2ig; pH-switch validated by alchemical FE" (71);
"De-novo pH-switch TNF binders; switch confirmed by alchemical free energy" (73);
"RFdiffusion + ProteinMPNN + AF2; FE-validated multi-His pH switch" (65))

## Description

I designed these as de-novo minibinders to the conserved receptor-binding groove of TNF-α, then built
in a pH switch so they hold at pH 7.4 and let go at pH 6.0. Backbones come from beta-pairing
RFdiffusion, sequences from ProteinMPNN with no cysteines, and I keep only the ones AF2 initial-guess
believes as binders. I install the switch as two or more histidines, each placed next to a different
conserved lysine on the target, so that charging them at low pH pushes the binder off. The part I lean
on hardest is that I confirm the switch direction by direct alchemical free energy rather than by
design alone, and I run that same measurement against mouse TNF to check cross-reactivity. I select on
confirmed switch direction plus lineage diversity, because switch magnitude and ON-state affinity are
not reliably rankable in silico on this target. The 18 proteins are ordered so the designs whose switch
fires in both human and mouse come first.

## Configuration

1. Generation. Beta-pairing RFdiffusion at the conserved Region-I edge of the receptor-binding groove,
   two-protomer target from 1TNF, with sampling noise set to zero so the binder makes real contact.
2. Sequence design. ProteinMPNN, zero cysteines.
3. Binding filter. AF2 initial-guess against a full-atom two-protomer target, keep pae_interaction
   below 10, with a known de-novo binder as a positive control to calibrate the threshold.
4. Amplification. Partial diffusion of the confirmed binders to grow each cross-reactive shape.
5. Switch install. Two or more spatially separated histidines, each adjacent to a distinct conserved
   target lysine, so the protonated histidine repels the lysine at pH 6.
6. Switch direction, confirmed four ways. (a) Two-state explicit-solvent MD occupancy in OpenMM with a
   His-to-Gln isostere as a falsifiable causal control; (b) PROPKA-Wyman and pKAI linkage; (c) GFN2-xTB
   on the histidine-lysine pair; (d) direct alchemical free energy of His0 to His+ in the bound complex
   versus the free binder (positive means release at low pH), OpenMM amber14SB and TIP3P, PME, hydrogen
   mass repartitioning at 4 fs, 11 lambda windows by thermodynamic integration, box-matched legs, run
   at 6 then 12 then 24 independent seeds.
7. Cross-species. Co-fold each binder on mature mouse-TNF (Boltz-2), gate each mouse pose on its
   binder-side paratope and switch-histidine placement against the human designed pose, then run the
   same alchemical free-energy sign on the mouse complex (216 production legs across the set).
8. Cross-reactivity footprint. Conserved epitope, 89 to 95% human-to-mouse interface identity, with the
   cationic partners present in both species.
9. Developability and novelty. Zero cysteines, an orthogonal multi-signal expression screen, structural
   novelty by foldseek TM and mmseqs, and final novelty on the Proteinbase portal (pass at 3/4 or more).
10. Selection. Diversity across distinct lineages; two anti-switches found and removed by the
    free-energy audit; 18 designs submitted, ordered with the five dual-species releasers first.
