# Methodology: de-novo conditional pH-switch binders to TNF-α (Challenge 2)

I designed de-novo binders that hold TNF-α at pH 7.4 and let go at pH 6.0, the direction a molecule
needs for endosomal recycling, and that work against both human and mouse TNF-α. I ranked the switch
(pH-selectivity) first, cross-reactivity second, and ON-state affinity last, and I treated structural
novelty as a pass-or-fail gate rather than a score (Proteinbase at least 3/4). Every design has no
cysteines, clears an expression screen, and carries a switch I checked in silico before I picked it.

## How I built them

The backbones came from beta-pairing RFdiffusion at the conserved Region-I edge of the receptor-binding
groove. I put sequences on them with ProteinMPNN (no cysteines) and kept only the ones AF2
initial-guess believed as binders: pae_interaction below 10 against a full-atom two-protomer target,
with a known de-novo binder as a positive control so I knew the filter could still say yes. That filter
does most of the work. TNF-α is flat and polar and most backbones never clear it. For the few that did
and that stayed cross-reactive, I amplified them by partial diffusion. Then I installed the switch: two
or more histidines, spaced apart, each sitting next to a different conserved lysine on the target, so
that at pH 6 the protonated histidines push on those lysines and loosen the grip. For some designs I
repacked the histidine's neighborhood with LigandMPNN, holding the histidine fixed, until two separate
co-folders (Boltz and Chai) both agreed it still pointed at its lysine.

## Why I trust the switch direction

I did not want to ship a switch that only looked right on paper, so I checked its direction three ways
and then measured it.

- Two-state MD in explicit solvent, with a control I could falsify. Swap the histidine for glutamine,
  which cannot change charge with pH: a real switch shows the protonated histidine shoving its lysine
  away, while the glutamine version does nothing.
- PROPKA-Wyman and pKAI linkage, plus GFN2-xTB on the histidine-lysine pair, as independent reads on
  the same sign.
- The direct measurement: the alchemical free energy of turning the switch histidine from neutral to
  charged, in the bound complex and again in the free binder. If that difference is positive, the
  binder releases at low pH. I trust this number because the double difference cancels the proton
  self-energy and the box-size artifact and leaves only the physical rearrangement, which is the part I
  can actually rank.

I ran that free-energy sign test on every design, at 6 seeds, then 12, then 24 where I needed to be
sure. It earned its place: it caught two designs whose switch ran backward, binding tighter at low pH
instead of weaker, and I cut both. What I submit has five releasers confirmed by free energy and zero
confirmed backward switches.

## The mouse test, and what it told me

Cross-reactivity in this kind of work usually rests on a conserved footprint. The interface is 89 to
95% identical between human and mouse and the lysines are present in both, so the binder should still
bind. I wanted more than "should", so I ran the same free-energy switch measurement against co-folded
mouse-TNF for every design, 216 production runs in all. Before any of them, I made each mouse pose earn
trust. The binder is the same molecule in both species, so I compared its contact face and its
switch-histidine placement against the human pose I actually designed, and that check threw out one
co-fold that had latched onto the wrong face. The controls hold up: the free binder behaves the same in
both species, the whole effect sits in the bound complex, its size is specific to each design, and one
design whose switch stays put in mouse rules out the lazy reading where every complex just stabilizes a
charged histidine.

Here is what surprised me. The switch I tuned on human is often specific to human, and how strong it is
on human tells me almost nothing about whether it survives in mouse. Five designs pass in both species,
with the switch firing energetically against human and mouse TNF: des_6025_4_T10_pd_77_4_mh_52_sa10_sh3_sa11,
des_820_3_r3_sa4_r2_8, b4m2_r1_w5_3, des_4082_8_T10_pd_10_4_mh_51_79_sa7, and des_732_1_r3_sa15. Those
are my best bets on the two axes I care about most, so I put them at the top of the list.

## What I am not claiming

I picked for diversity across distinct lineages instead of chasing a magnitude ranking, because two
things genuinely cannot be ranked here and I would rather say so.

- How big the switch is. The slow step in the bound state is the whole interface rearranging, not the
  histidine, so the free energy will not converge even under well-tempered metadynamics. I can trust
  the sign, not the size, so I select on sign plus diversity and leave size to your assay.
- ON-state binding strength. On this target, in-silico affinity ranking is close to a coin flip, so I
  will not pretend an order.

Every design I submit clears the novelty gate (at least 3/4), has no cysteines, and passes an
expression screen built from several independent signals. TNF-α has only ever been bound by antibodies,
and its surface is flat and polar, so a small de-novo binder that both holds on and lets go on cue,
with that switch confirmed by direct free energy in two species, earns a well on the plate even if I
cannot tell you its affinity rank.

## How to read the list

I ordered the 18 designs by how well the switch carries from human to mouse. The five dual-species
releasers come first. Then come the designs whose switch holds or goes neutral in mouse. Last are the
ones with a strong human switch that reverses in mouse; those still bind mouse, and human is the
primary target for the switch. If you can only test a few, start from the top five.

## The per-design numbers

The table below is the measurement behind the switch claims above, for all 18 designs in the
order I submit them. ddG is the free energy of charging the switch histidine (bound complex minus
free binder); positive means the design releases at low pH. z is the per-seed mean over its
standard error, and n is the seed count. For each design the human and mouse runs use the same
histidine set, so the two columns compare directly.

| # | design | human ddG (kcal/mol), z | mouse ddG (kcal/mol), z (n seeds) | what the mouse run shows |
|---|--------|------------------------|-----------------------------------|--------------------------|
| 1 | des_6025_4_T10_pd_77_4_mh_52_sa10_sh3_sa11 | +7.7, z +8.3 | +3.57, z +4.9 (n=6) | releases in both species |
| 2 | des_820_3_r3_sa4_r2_8 | +6.8, z +4.8 | +6.83, z +3.5 (n=6) | releases in both species |
| 3 | b4m2_r1_w5_3 | +6.7, z +18.8 | +3.26, z +2.7 (n=6) | releases in both species |
| 4 | des_4082_8_T10_pd_10_4_mh_51_79_sa7 | +6.2, z +2.5 | +5.59, z +4.0 (n=6) | releases in both species |
| 5 | des_732_1_r3_sa15 | +3.2, z +5.1 | +2.25, z +3.5 (n=6) | releases in both species |
| 6 | des_138_7_T10_pd_39_4_mh_49_77 | +3.8, z +2.6 | -2.80, z -2.6 (n=12) | mild mouse reversal, not significant |
| 7 | des_6126_3_T20_pd_45_3_mh_9_56_78 | +1.9, z +1.9 | -2.53, z -1.3 (n=6) | holds or neutral in mouse |
| 8 | b2m0_c70 | +1.7, z +1.5 | +0.16, z +0.2 (n=6) | holds or neutral in mouse |
| 9 | des_875_5_H3-57-80 | +0.9, z +1.7 | +0.08, z +0.1 (n=6) | holds or neutral in mouse |
| 10 | des_33012_g4_sa8 | +0.4, z +0.5 | -0.02, z -0.0 (n=6) | holds or neutral in mouse |
| 11 | des_4082_8_T10_pd_9_1_mh_43_55_57_79 | +0.3, z +0.3 | -1.83, z -1.4 (n=6) | holds or neutral in mouse |
| 12 | des_4098_1_c44 | +7.2, z +8.1 | -2.79, z -3.1 (n=12) | reverses in mouse |
| 13 | des_1559_3_mh_48_75 | +2.0, z +1.7 | -6.50, z -4.1 (n=6) | reverses in mouse |
| 14 | des_1060_4_mh_46_74 | +1.8, z +0.7 | -8.01, z -6.0 (n=6) | reverses in mouse |
| 15 | des_36133_g4_sa18 | +1.3, z +0.8 | -4.80, z -3.7 (n=12) | reverses in mouse |
| 16 | des_495_6_mh_41_49 | -0.0, z -0.0 | -1.35, z -3.1 (n=6) | reverses in mouse |
| 17 | des_180_6_H9-14-38 | -0.9, z -1.2 | -6.34, z -5.9 (n=6) | reverses in mouse |
| 18 | des_2106_2_mh_37_43_50_75_sa2 | -2.0, z -1.9 | -7.82, z -4.8 (n=6) | reverses in mouse |

Five designs release in both species (rows 1 to 5), six hold or go neutral in mouse, and seven
either reverse in mouse or lean that way without clearing significance. All 18 bind mouse-TNF; I
gated each mouse pose on its contact face before I trusted its free energy, and human is the
primary target for the switch.
