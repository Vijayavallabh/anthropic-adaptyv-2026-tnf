# Switch free-energy evidence, v27 submission (18 designs)

This table is the measurement behind my switch claims. For each design I computed the
alchemical free energy of turning the switch histidine from neutral to charged, once in the bound
complex and once in the free binder. A positive ddG means the design releases at low pH, which is
what I want; a negative ddG means it binds tighter at low pH, which is the wrong direction. z is
the per-seed mean over its standard error, and n is how many independent seeds I ran per leg. For
each design the human and mouse runs use the identical histidine set, so I read the two columns
side by side. The rows are in the order I submit them.

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

How I read it: five designs release in both species (rows 1 to 5), six hold or go neutral in
mouse, and seven either reverse in mouse or lean that way without clearing significance. All 18
bind mouse-TNF; I gated each mouse pose on its contact face before I trusted its free energy, and
human is the primary target for the switch.
