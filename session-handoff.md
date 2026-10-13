# Session handoff — TNF Challenge 2

Updated 2026-10-12/13 (SESSION 13). **Deadline ~2026-10-19 11:30 UTC** (user reconfirmed 7d 02h 59m left on
2026-10-12 ~08:25Z — the Stop-hook /goal still quotes a stale "1d 17h"; IGNORE that, real runway ~7 days).

## ⚖ v21-vs-v19 HEAD-TO-HEAD (SESSION 13, identical gauntlet on both — HONEST verdict)

Ran the SAME in-silico gauntlet on v19's 20 and v21's 17. **v21 is NOT "better in every metric":**
- **Binding floor:** TIE (both all-pass af2ig pae_int 4–6).
- **Cross-react:** TIE (both ~90–97%).
- **Expression:** ~TIE (v21 better floor — it cut des_84_3 at sol 0.452).
- **Switch DIRECTION (#1):** TIE — v19 raw multi-His switchers are ALREADY 13/15 MD His-causal with strong gaps; v21 also ~all His-causal. (H43 is MD-strong +1.85 despite weak static Wyman; des_84_3 cut for expression/novelty, NOT switch.)
- **Novelty (hard gate):** **v19 WINS** — all 20 portal-confirmed ≥3/4 vs v21's 13 pending re-check.
- **Fold-robustness (Boltz∩Chai both-fold) + diversity (15 lineages):** **v21 WINS** (confidence + hedge, not judged magnitude).

**So the v19↔v21 choice is a TRADEOFF, not a domination.** v19 = fully-gated, strong switches, comparable everything. v21 = same measurable axes + more diversity + fold-robust confidence, BUT needs the portal re-check on 13 members. Recommend v21 only IF the re-check clears; otherwise v19 is the safe strong pick. Full data: dossier `COMPARISON_v21_vs_v19`.

## ★ RECOMMENDATION (SESSION 13) — PRIORITY-DEPENDENT: v21 leads #1, v22 = pragmatic, v19 = zero-work

The three options differ on re-check burden and on #1-axis depth, NOT on a provable domination (switch magnitude
is unrankable in-silico; binding magnitude is wet-lab-decided; novelty is portal-gated). Pick by priority:

- **MAX #1-axis strength → v21** (`runs/push/submission_v21.csv`). On the top judged axis (expected pH
  fold-change, proxied by the both-fold-confirmed COUNT of separated in-window His) v21 LEADS: 8 members carry
  2–3 both-fold-geometry-confirmed cation-repelling His (3 three-cat + 5 two-cat), vs v22's 2 and v19's 0. Cost:
  **13 portal re-checks.** Also wins diversity (15 lineages) + fold-robustness.
- **PRAGMATIC strict-improvement over confirmed v19 → v22** (`runs/push/submission_v22.csv`, md5
  **`abaec06790766f25282133d929083187`**, 20 designs, validator PASS, LOCAL-ONLY). = v19 with ONLY its **2
  MD-FAILING switchers** replaced: b1m9_A12H (d_repel −1.307 wrong-sign) + des_138_7_pd_24_4 (−0.097 non-causal) →
  des_36133_g4_sa18 (3-cat, MD gap +2.30) + des_732_1_r3_sa15 (3-cat). **≥ v19 on every axis, STRICTLY > on
  switch-DIRECTION (fixes 2 dead switchers) + diversity; TIE on binding/cross-react/expression; only 2 portal
  re-checks** (the other 18 already portal-confirmed ≥3/4). NOTE: does NOT dominate v21 on #1 — its 18
  v19-derived members are MD-direction-OK but NOT both-fold geometry-confirmed. (Earlier "strictly better on
  binding" was an overclaim — pae 4.5 vs 4.9/6.0 is noise; corrected to TIE.) Dossier
  `docs/SUBMISSION_v22_candidate.md`.
- **ZERO-WORK fully-gated fallback → v19** (all 20 portal-confirmed ≥3/4). Recommended if no portal re-check is done.

## v21 — max-diversity option (SESSION 13), USER-PORTAL-GATED (13 re-checks; ties v19 on measurable axes)

- **`runs/push/submission_v21.csv`** (md5 **`f7f246f41a3a062b1922f0f980a9928a`**) — **17 LEANER DIVERSE designs**,
  validator **PASS**, LOCAL-ONLY. = v20b minus same-fold cousins + weakest-switch members, + **7 NEW-lineage
  refills** (des_732_1_r3_sa15 [3-cat] / des_1169_2_r3_sa7 / des_243_3_r3_sa2 / **des_820_3_r3_sa4_r2_8** [2-cat,
  TM 0.575] + fresh-campaign trio des_36133_g4_sa18 [3-cat] / des_34125_g4_sa15 [2-cat] / des_33012_g4_sa8 [2-cat],
  all MD-direction-confirmed His-causal); sh3_sa11 replaces sa10 (novelty-safety, struct-TM
  0.913 vs 0.959). **Selection principle:** switch magnitude is NOT reliably rankable in-silico (pKa err 2-3u >
  1.4u window; two-engine Wyman rank disagrees ρ=0.56) → pick many foldseek-DISTINCT direction-confirmed switchers,
  each strong on binding / cross-react ≥90% / structural novelty / solubility, NOT a magnitude rank. Dossier
  `docs/SUBMISSION_v21_candidate.md` (UPLOAD_CHECKLIST + FINAL_COMPOSITION are authoritative).
  **EXPANSION COMPLETE:** regen round-3 + round-2 cousins → des_820; then a FRESH de-novo campaign (1200 bpair
  backbones → 49 binders → 10 ultra-novel switch-piped) → round-4 added 3 (des_36133/des_34125/des_33012). Pattern:
  ALL 1-cation switchers fail the MD direction gate; the ≥2-cation class survives → only ≥2-cation kept.
  **UPLOAD HARD-GATE: only 4 members are byte-identical to portal-confirmed v19 (c44, H9-14-38, b14, H43); the
  other 13 (5 shell-adapt saN + sh3_sa11 + 7 new-lineage) have sequences NOT in v19 → EACH needs USER portal
  novelty re-check before upload.** (Corrects an earlier "10 byte-identical" claim — shell-adapt saN were never
  portal-checked.)

## Gated-safe default (fully portal-novelty-gated)

- **`runs/push/submission_v19.csv`** (md5 **`7756055f0b20bcd6994fefb644271db7`**) — 20 designs, validator
  **PASS**, user portal-checked **all 20 ≥3/4**, 20/20 switchers. **Recommended upload if no portal re-check.**
- **Gated fallback:** `submission_v18.csv` (md5 **`524722158467273c0e9f57f5402650ed`**, 19/20 switchers,
  4×4/4-novel). Floor: `submission_v16.csv` (`a524d77f`).
- **Upload is USER-GATED** — do not auto-upload. CSVs under `runs/` (gitignored). Account `j_v_v_07`.
- Portal CSV format: `name,sequence,molecule_class`; de-novo minis `single_chain`; unique names AND seqs.

## SESSION 12 (2026-10-12) — v20 fold-robust shell-adaptation campaign (the live improvement)

**v20 = v19 with 9 fold-FRAGILE switchers replaced by shell-adapted fold-ROBUST variants.** Fold-robust = the
switch-His engages a cross-conserved TNF cation in BOTH folders (Boltz∩Chai), not one. Portfolio fold-robust
**8/20 → 17/20** (round-2 adds b2m0_c70 sa17 + des_875_5 sa1-r2 → **19/20**; only des_138_7_pd_39_4 stays
fragile — Chai never agrees in 2 rounds). `runs/push/submission_v20.csv` md5
**`6f4b61d0a55de2e4e0db205970ec6821`**, validator PASS. Dossier `docs/SUBMISSION_v20_candidate.md`;
picks `runs/foldrobust/v20_final_picks.tsv`.

- **Mechanism:** `scripts/shelladapt_shell.py` — LigandMPNN shell-redesign, switch-His FIXED, repack binder
  residues within 10Å of target OR 8Å of a switch-His → locks the His rotamer so both folders agree.
  Per-member campaign = `shelladapt_shell.py` → LigandMPNN (protein_mpnn, 16 variants/member, omit-C) → Boltz +
  Chai fold (`run_boltz_tnf.sh` / `run_chai.sh`, target `data/targets/human_tnf.seq`) → Boltz∩Chai per-His
  triple-confirm (`scripts/triple_confirm_mut.py`) → `cif2pdb_af2ig.py` + `run_af2ig.py` for binding re-confirm.
  sa2 = 176 variants / 11 members; round-2 = 96 / 3 stubborn. (Per-run driver shell scripts are session scratch,
  regenerable from these durable scripts + the dossier.)
- **Binding re-confirmed per design (REQUIRED — shell-adapt rewrote 60-75/80 residues).** af2ig initial-guess
  on each NEW sequence (binder=chainA + 152aa protomer; `scripts/cif2pdb_af2ig.py` converts the Boltz fold).
  **All 9 adopted PASS** (rmsd<2/pae_int<10/plddt_b>80); 8/9 match-or-beat parent. The most-SOLUBLE variants
  sa5/sa3 FAILED binding (pae 23/22) → rejected: af2ig is seq-specific, the binding filter is load-bearing.
  Scores `runs/sa2_af2ig{,2}/out.sc`. des_6025 sa10 binding CONFIRMED (4.78/92/0.45) — never actually run before.
- **Joint re-selection** (fold-robust AND max-solubility over 16 variants/member): all 9 Protein-Sol ≥0.59.
  Flagship by geometry = des_2106_2 sa2 (3 His, only bc=3 of 176) — BUT MD-weak (2 His crowd K93). Cleanest =
  des_6025 sa10 + des_4082_8 ×2 (2 His on 2 SEPARATED distinct cations, MD-selective TRUE).
- **Cross-react:** every adopted winner 89-95% epitope id human↔mouse, all engaged cations K/R in both species
  (`runs/foldrobust/crossreact_interface_conservation.md`). Mechanistic, fold-independent.
- **Novelty:** sequence 0/20 >70% id to PDB (all de novo). Backbone UNCHANGED + winners fold to parent backbone
  (af2ig rmsd<2) → structural novelty ~ parent's portal ≥3/4. **REQUIRED: USER portal re-check on the 9 (soon 11)
  reseqenced swaps before any v20 upload.** Local single-seq structural pre-screen staged (`tmp/sa2/struct_novelty.sh`,
  esmfold→foldseek vs `runs/novelty_struct/pdbdb`).
- **HONEST SCOPE:** both-fold = a STRICTER switch-GEOMETRY filter (51% single-fold → 28% both-fold across 176;
  higher-confidence #1), NOT in-silico binding proof (TNF docking stays wall-limited; #3 affinity wet-lab-decided).
  Two-state MD (3-rep, a noisy CONFIRMER not a hard cutter) corroborates the separated-cation designs and flags
  the crowded-His ones — vindicates separated-His. Firming up: higher-rep MD + parent MD on final flagships.

**RUNNING at handoff (all 8 GPUs):** two-state MD on the 9 adopted (`runs/sa2_md`, GPUs 0-3) + campaign-2 funnel
(1200 fresh bpair-RFd backbones → MPNN→af2ig, GPUs 4-7, `runs/bpair/bpair_fresh/funnel`) for new-lineage binders
(novelty/diversity). **Staged when GPUs free:** round-2 winners' af2ig (→19/20), structural-novelty pre-screen.

**NEXT (7-day runway → bold, thorough):** (1) finish MD + round-2 af2ig → finalize v20/v20b. (2) Funnel passers
→ new fold-robust switchers from diverse lineages (novelty-4/4 + expression hedge). (3) Push switch SHARPNESS:
designs with His on 3 SEPARATED distinct conserved cations (des_2106_2 showed 2-on-K93 is anti-cooperative).
(4) Assemble the optimal 20 = max fold-robust/MD-selective + novelty-diverse + cross-react + soluble. Keep
v18/v19 gated; upload USER-GATED.

## What v16 is (the LONG SHOT win)

v16 = v14 with **4 genuinely NEW-shape de-novo switchers** swapped in for the **4 weakest v14 members**
(2 wrong-sign + 1 noisy + 1 flat, all MD-NONRESPOND). Portfolio two-state-MD trigger **13/20 → 17/20
CONFIRMED**; +2 new lineages; diversity preserved.

Swap-IN (new lineages des_138_7, des_4082_8 — ≤2 each):
- `des_138_7_T10_pd_24_4_mh_49_58_77`  dSwitch +1.11, x-react 0.87, pLDDT 0.87, novelty 3/4
- `des_138_7_T10_pd_39_4_mh_49_77`     dSwitch +1.09, x-react 0.90, pLDDT 0.90, novelty 3/4
- `des_4082_8_T10_pd_10_4_mh_51_79`    dSwitch +0.95, x-react 0.90, pLDDT 0.81, novelty 3/4
- `des_4082_8_T10_pd_9_1_mh_43_55_57_79` dSwitch +0.46, x-react 0.91, pLDDT 0.86, novelty 3/4

Swap-OUT: `b1m7_2A57H_K183_r1_w5_3` (−1.08 wrong-sign), `des_59_1_H51` (−0.62 wrong-sign),
`b1m12_A41H_vsR27` (+0.07, 2·SEM 1.70 noisy), `des_3050_6_c19` (+0.25 flat). Full dossier:
`docs/SUBMISSION_v16_candidate.md`; the complete 20-design v14 MD trigger matrix: `docs/TRIGGER_MATRIX_v14.md`.

## The long-shot pipeline (reproducible; all in scripts/)

User directive "go with the long shot" = genuinely NEW shapes, not partial-diffusion cousins of the 4
already-confirmed lineages. The winning chain:
1. **de-novo generation** — `scripts/run_bpair_rfd.sh <subdir> <ndes> <gpu> [startnum]` (β-pairing RFd1,
   fold-conditioned to the conserved Region-I edge strand B114-120). Scaffold-guided → structural TM high
   (novelty caps at 3/4, seq-novel) but cross-reactive.
2. **binding funnel** — `scripts/denovo_stage2.sh "<glob>" <base> "<gpus>" [nseq]` (ProteinMPNN chain A
   0-Cys → full-atom af2in → af2ig). De-novo binding yield ~0.5% → few shapes.
3. **cross-react** — `XRBASE=<base> scripts/v15_xreact.sh [TOPN] "<gpus>"` (Boltz human+mouse, min-iptm≥0.70).
4. **AMPLIFY the cross-reactive new shapes** — `PD_BASE=<base2> PD_GPUS="..." scripts/partial_denovo.sh [ndes]`
   (partial-diffuse the proven cross-reactive de-novo shapes; HIGH yield ~25-30% → backbone-diverse cousins
   of the NEW lineages). **This is the multiplier that made the long shot work.**
   ⚠ renumber seed chain B to 1-152 first (af2pred poses number it 81-232; contig `B1-152` else finds nothing).
5. funnel the cousins (step 2) → cross-react (step 3) → **multi-His** `MHBASE=<base> MHPOSEGLOB=<poses>
   scripts/v15_multihis.sh [CAP] "<gpus>"` (bpair_multihis_build.py separated-His + Boltz + mh_screen) →
   **two-state MD** `scripts/run_md_release.sh <manifest> <out> "<gpus>" 3000 4` + `md_verdict.py`.
6. expression: `esmfold_monomer.py` (pLDDT) + `fold_monomer_boltz.sh` + `solubility_predict.py`; novelty
   proxy: foldseek TMalign vs `runs/novelty_struct/pdbdb/pdb` (<0.80=novel) + mmseqs seq-id; portal = USER.

**KEY LESSON (reconfirmed):** static Wyman linkage ≠ MD trigger. The lead parent `des_138_7_mh_47_53_80`
was MD-NONRESPOND (wrong-sign −0.78) despite static −1.05; only the backbone-diversified cousins MD-confirmed.
Always gate on the two-state MD, never static linkage alone.

## v17 candidate — RESOLVED, USER-PORTAL-GATED (SESSION 11 continued)

The v17 push played out: **des_6025 yielded ONE MD-confirmed new-lineage switcher**
(`des_6025_4_T20_pd_64_1_mh_34_50_52_74`); **des_6126 and des_8xxx were dead ends** (des_8xxx 0/17 strong).
That one switcher was **fully re-validated on the exact His-variant sequence**:
- #1 switch: two-state MD d_repel **+0.372 CONFIRMED** (2·SEM ≈0.27), isostere ctrl +0.308, selective_release true.
- #2 cross-react: **re-folded** the His-variant (the mouse Boltz output had been cleaned up) — human iptm
  **0.962**, mouse iptm **0.938**, min **0.938 PASS**. (The earlier "0.887" was stale/parent — actual is stronger.)
- expression: solub 0.87 PASS. Seq-novelty: **NO PDB hit (novel)**.
- **One open risk:** foldseek structural max-TM **0.945** (above v16's 0.86–0.93) → likely 3/4 but
  **portal-UNVERIFIED**. Novelty is a hard per-design gate; only the portal (USER) is authoritative.

→ `runs/push/submission_v17.csv` (md5 **`9357b6905ed81179da5137a16dd9dd5a`**, validator PASS) =
v16 with `b4m2_A57H` → des_6025 (18/20 MD-confirmed, 3 de-novo lineages). Dossier:
`docs/SUBMISSION_v17_candidate.md`.

## BOLD thrust IN FLIGHT (SESSION 11 cont.) — Genie3 at the CRO-validated Region-I GROOVE

Re-examined the "bold headroom exhausted" claim (goal-hook pushback) and found a genuinely untried,
evidence-grounded lever: **the Anthropic paper's CRO-measured TNF campaign got 8% (12/150) and EVERY binder
came from Genie3 + PXDesign (RFd3+BoltzGen=0); Genie3 gave ALL 3 human/mouse cross-reactive binders — but I
only ever ran Genie3 at K90 (wrong site). Genie3 at the validated Region-I GROOVE is UNTRIED**
(the `runs/regionI*` dirs are RFd3, which walls at 0). See [[lesson-tnf-generators-genie3-pxdesign-not-rfd3]].

STAGED + ready to fire the instant a GPU window opens (do NOT clobber the foreign job):
- **Generation:** `scripts/gen_genie3_regI.sh "<gpus>" [N=600] [seedbase] [tag]` — shards `g3batch.py` on the
  `tnf_regIbench`/`tnf_regI` dataset (hotspots B69/B70/B92 + C109/C110 = groove V74/L75/I97 + W114/Y115),
  ~26 s/design, per-GPU triton cache. Output: `runs/genie3_regI/shard_g*/tnf_regI/pdbs/*.pdb` (binder=A).
- **Target for funnel:** `data/targets/1TNF_regI_BC.pdb` (304 CA, two protomers B+C). `build_fullatom_af2in_tnf.py`
  supports a 2-chain target via `--target-chains "B,C"` (maps to B,C). Funnel = MPNN(chain A, 0-Cys) →
  build_fullatom(2-chain) → af2ig — BUT af2ig is near-chance for TNF (AP 0.16), so OVERSAMPLE + pose-gate at
  the groove, don't hard-floor (see [[reference-anthropic-tnf-dataset-regionI-recipe]]). Then switch → two-state
  MD → novelty, as v16/v17.
- **LIVE (2026-10-11 ~20:15Z):** GPU window opened (mmbu released 4-7); **Genie3 generation RUNNING on GPUs
  4,5,6,7** (PID in `$TMPDIR/tmp/gen_regI_main.log`, N=800, ~182/800 so far, ~70min left). mmbu still
  holds 0-3 — do NOT touch. Gen-completion waiter running (re-invokes to run the funnel).
- **fbb builder works + Genie3 BREAKS the groove-placement floor:** `scripts/build_genie3_fullbb_bc.py`
  (PULCHRA at `.../recon_venv/bin/pulchra`) reconstructs binder full-backbone + joint-Kabsch grafts the real
  2-protomer target; 69/69 built, and the binder **docks the groove (min CA 4.9 Å, 13/93 binder res <8 Å)** —
  PXDesign couldn't (0/320). So Genie3 CAN place at the two-protomer groove.
- **Funnel built + validated:** `scripts/funnel_genie3_regI.sh "<gpus>" [nseq] [tag]` = fullbb_bc → ProteinMPNN
  (chain A, 0-Cys, oversample) → `build_fullatom_af2in_tnf.py --target data/targets/1TNF_regI_BC.pdb
  --target-chains "B,C"` → af2ig (WIDE filter: af2ig is weak for TNF). Then ipSAE+pose ranking (Boltz h+m +
  tools/ipsae.py + groove-hotspot pose-gate) does the REAL ranking, calibrated on the 12 CRO binders
  (`_research_downloads/tnfa_claude_validated_binders.fasta`). Then switch (bpair_multihis_build.py near
  conserved K65/K90/K98) → two-state MD → novelty.
- **Research (2 agents) findings adopted:** (1) Genie3-at-groove is THE method that broke TNF + gave mouse
  x-react — confirms the thrust. (2) Rank by **ipSAE_min ensemble** (ESMFold2+Protenix-v2, 5 seeds, z-avg) +
  Protenix∧af2ig AND-gate + pose-gate; NOT raw ipTM/af2ig. (3) **MEMORY CORRECTION:** the "0.16" is the
  Spearman ρ(co-fold score, AFFINITY) among binders (weakly +), NOT a TNF AP, and NOT anti-correlation; af2ig/
  ipSAE DO separate binder/non-binder above chance (ipSAE 1.4× af2ig AP, Overath). So my bpair af2ig pae<10
  passes are valid positives. (4) pH-switch: separated-His-on-conserved-cations is the FRONTIER, not superseded;
  His uniquely matched to 6.0-7.4; Proton-PottsMPNN de-prioritized (window 4.0-5.8). Optional CPU upgrade:
  JustHISpKa (His-specific burial-aware pKa, JCIM 2025) as a 3rd pKa class — logged, not chased (Amber-heavy).
- **DERISKED POSITIVE (2026-10-11 ~20:30Z):** (a) Genie3 PLACES at the two-protomer groove — on a 20-backbone
  sample, **50% cover ≥2 groove hotspots at <6 Å** (PXDesign floored 0/320). (b) The groove has **3 CONSERVED
  cations within reach** (dataset→mature offset +5, verified by residue identity): **K98 (3.8 Å), K65 (4.4 Å),
  K112 (6.6 Å)** — all human=mouse conserved. So a separated multi-His switch repelling K65+K98+K112 is
  installable AND cross-reactive on the same conserved cations = the Baker SOTA form; all 3 scored axes align
  at the groove. (Avoid R138/dataset-R133, divergent.)
- **Pipeline fully staged + validated:** `build_genie3_fullbb_bc.py` (2-protomer fbb, 69/69 built),
  `funnel_genie3_regI.sh` (fullbb→MPNN→af2ig, syntax OK), `rank_genie3_regI.sh` (2-protomer Boltz h+m →
  ipSAE_min + groove-pose-gate + cross-react, python compiles; ipSAE col=f[5] type=f[4]),
  `runs/genie3_regI_calib.csv` (12 CRO binders + 12 scrambles for ranker calibration). SWITCH-STAGE TODO: tweak
  `bpair_multihis_build.py` to scan chain C cations too (currently chain B only, line ~43 `cc=="B"`) — the
  groove's K65 is on chain C.
- **Next (GPU-gated, waiter re-invokes on gen completion ~20min):** funnel on all 800 → calibrate ranker on the
  12 → rank Genie3 binders by ipSAE+pose → multi-His (K65/K98/K112) → two-state MD → novelty → v18 hedge.
- **MD-STAGE 2-PROTOMER FIX (identified, apply when survivors exist):** `md_release_tnf.py` picks binder/target
  by CHAIN SIZE (2 largest), so a 3-chain groove complex (binder A 88 + B 152 + C 152) mis-identifies a
  protomer as the binder. FIX before MD: merge the two target protomers into ONE 304-res target chain
  (relabel C→B, renumber C residues +152) → standard 2-chain complex (target 304 + binder 88) → run
  md_release_tnf.py unmodified (correct: the MD restrains the whole target rigid; nearest-Lys then searches
  both protomers for the conserved anchor K65/K98/K112). Don't build this until the funnel yields binders.
- **LIVE (2026-10-11 ~20:45Z):** funnel RUNNING on 4-7 (PID 667950, fullbb→MPNN→af2ig, ~2h); ranker
  CALIBRATION RUNNING on 6,7 (PID 687999, 12 CRO binders + 12 scrambles, decisive ranker validation, ~18min).
  Waiters set. mmbu holds 0-3. Rational order: WAIT for the calibration + funnel gates before building more
  downstream (don't speculatively plumb past an unconfirmed gate).
- **CALIBRATION RESULT — Boltz ipSAE ranker is UNUSABLE for Genie3-type (rigorous, via 12 CRO binders):**
  with corrected parsing, Boltz-msa-empty ipSAE recovers only 3/12 known binders (all PXDesign-type:
  rank23 ipSAE 0.447/iptm 0.69, rank9, rank19) and MISSES all 8 Genie3-L83 binders (ipSAE 0.0, incl the 2nM
  cross-reactive ones) — Boltz blind-fold can't place the Genie3 fold-type. So `rank_genie3_regI.sh` (Boltz
  ipSAE) is DROPPED for selection. Use instead: **af2ig (initial-guess, NOT blind-dock — the funnel) +
  groove-placement + the MD-confirmable switch + mechanism-based cross-react** (conserved cations). Do NOT
  set up Protenix-v2 ensemble (heavy; even the Anthropic ensemble was only "above chance" for TNF). The
  calib run also showed Boltz DOES dock+score PXDesign-type, so the method isn't broken — it's fold-specific.
- **BUG FIXED — fbb name collision (was losing 600/800 backbones):** the 4 generation shards each emit
  `tnf_regI_0..199.pdb`, so basenames collided across shards (last-writer-wins → only 200 distinct fbb).
  Fixed `build_genie3_fullbb_bc.py` to prefix the shard (`g4_tnf_regI_0`); now 800 distinct. (Same class as
  the denovo_stage2 collision.) Also added SHARD_I/SHARD_N env for parallel fullbb workers + skip-existing.
- **LIVE (2026-10-11 ~21:06Z):** mmbu FINISHED (all 8 GPUs free). Funnel RE-running clean on ALL 8 GPUs
  (parallel 16-worker fullbb → MPNN → af2ig, 800 backbones × 4 seq, ~1h). Completion waiter alive. Process-kill
  note: `pkill -f` self-matches the shell (exit 144) — kill funnel-pipeline procs by PID or /proc scan;
  af2ig children are `uplifting-af2ig/bin/python` (not matched by 'funnel'/'build_genie3' patterns).
- **THRUST RESULT — rigorous NO-GO (2026-10-11 ~21:40Z):** Genie3 at the groove PLACES binders (~50% cover
  ≥2 hotspots <6 Å; PXDesign floored 0/320) but **af2ig FLOORS: 0/3200 at pae<16, distribution min 23.2 /
  median 27.5** — the flat-polar wall (binders fold great, af2ig never believes the interface). Ruled out the
  2-protomer dilution confound: single-protomer re-score of the top-16 gave pae 26-27, EVEN WORSE. Boltz ipSAE
  calibration also missed all 8 Genie3-L83 known binders. **No in-silico tool (af2ig/Boltz) can confirm or rank
  Genie3-type TNF binders — Anthropic needed WET LAB.** So a Genie3-groove design can't be shown to bind and
  would be a blind bet displacing an af2ig-confirmed bpair design. **v16/v17 stand as the submission.** The bold
  headroom for improving the in-silico-confirmable portfolio is exhausted (every generator walls at af2ig for
  de-novo TNF binding). See [[lesson-genie3-groove-floors-af2ig-like-all-generators]]. Pipeline committed
  (gen/fullbb/funnel/select/rank/switch_screen/merge_protomers) for reuse if a better binding oracle appears.
- **RECONSIDERED (bold hedge, 2026-10-11 ~21:50Z):** af2ig is near-chance for TNF AND my v16/v17 rest on a
  pipeline never CRO-validated for TNF, while Genie3+groove is the ONLY CRO-validated cross-reactive recipe.
  So betting all 20 slots on af2ig-pass (bpair) isn't obviously optimal — a DIVERSIFICATION hedge (1-2
  Genie3-groove designs from the validated recipe) is defensible. Key enabler: **169/800 Genie3-groove designs
  place at the groove (≥3 hotspots) AND reach ALL 3 conserved cations K65/K98/K112 within switch-range** — a
  BETTER cross-reactive multi-His switch setup than bpair/cb1 (which reach 1-2). So a Genie3-groove hedge can
  carry separated His on all 3 conserved cations = strongest cross-reactive switch, at the validated epitope;
  the ONLY unconfirmed axis is binding (af2ig floors, but af2ig is near-chance + this is THE recipe). Top-16
  selected → `runs/genie3_regI/hedge_top.txt`. **PREPPED (2026-10-11 ~22:05Z):** `scripts/prep_genie3_hedge_md.py`
  builds the MD switch-probe. Pivoted 2-protomer-merge → **SINGLE-protomer target**: merging B+C into one chain
  gave OpenMM internal-terminus template errors (fused C-term/N-term across protomers); TNF is a HOMOTRIMER so
  all 3 conserved cations (mature K65/K98/K112 = dataset 60/93/107) exist on EACH protomer — keeping binder A +
  the one protomer the His best cover loses nothing for the switch test and is topologically clean. Result: **12
  variants, 12/12 MD-loadable (CPU dry-run through build()), 11/12 install separated His on ALL 3 conserved
  cations** (manifest `$TMPDIR/tmp/hedge_md/manifest.tsv`; regenerable from the committed script +
  hedge_top.txt + runs/genie3_regI/af2in). GPU-gated behind the bpair MD. NEXT (when GPUs free):
  `bash scripts/run_md_release.sh $TMPDIR/tmp/hedge_md/manifest.tsv <out> "0..7" 3000 4` → `md_verdict.py`
  → any CONFIRMED = an MD-switch-confirmed, groove-placed, all-3-conserved-cation, novel Genie3 design. Present
  as a USER-decided diversification hedge — binding stays a recipe-bet (af2ig floors; af2ig near-chance for TNF),
  labeled honestly; NOT submitted over confirmed designs unilaterally. A positive switch result strengthens the
  #1-axis pitch; a negative result closes the hedge.
- **LIVE (2026-10-11 ~22:35Z, ~13.4h to deadline): BOTH MD tracks running in PARALLEL (no conflict).**
  Track A (bpair harvest, 13 variants) on GPUs 0-4, 8/13 verdicts done, ~30min to finish — waiter PID
  1590848 (fires on n>=13) re-invokes me → run `harvest_v18.sh` (Boltz xreact on freed 0-4) → assemble v18.
  Hedge switch-probe (12 variants) on the freed GPUs 5,6,7, ~3-4h — waiter (bg job) re-invokes me on 12/12 →
  `md_verdict.py $TMPDIR/tmp/hedge_md_out`. **SEQ-NOVELTY PRE-COMPUTED (mmseqs vs PDB, CPU):** ALL 13
  Track A candidates incl BOTH des_6126 variants = FULLY seq-novel (0 PDB hits, <<30% id; POSCTRL TNF hit 100%
  = valid). So the seq axis of portal novelty is satisfied for every v18 candidate; only the STRUCTURAL axis
  (foldseek TM, GPU single-seq fold) remains — des_6025 variants ≈ des_6025 (TM 0.945, 3/4), des_6126 TBD.
  **harvest_v18 solubility step FIXED** (was wrong `-o` invocation + wrong column `solubility_score`; was
  silently dropping the expression signal). md_verdict parse (CONFIRMED=$NF, dSwitch=f[5]) verified correct.
- **v18 PRODUCED + VALIDATED (2026-10-12 ~00:20Z).** Track A done: 5/13 CONFIRMED (des_6126 1.36 [4th
  lineage, 4/4 reps positive], des_6025 mh_52_74 0.93, +3 weaker des_6025). ALL 5 harvested: xr(min h/m)
  0.91–0.95, solub 0.82–0.92, 0-Cys, seq-novel. **v18 = v16 with 2 swaps** (des_3050_6→des_6126,
  b4m2→des_6025 mh_52_74), KEPT des_818_4 (4/4-novel, avoid a des_6025 cousin → "diverse not cousins").
  19/20 MD-confirmed, 4 de-novo lineages, all ≥3/4 novel (foldseek TM des_6126 0.948 / des_6025 0.870, both
  helical→3/4). `runs/push/submission_v18.csv` md5 **524722158467273c0e9f57f5402650ed**, validator PASS.
  Doc: docs/SUBMISSION_v18_candidate.md. Bugs fixed in harvest/assemble this session: solub -o/column,
  _af2pred seq lookup, v18_candidates needs a sequence column (augment from confirmed_seqs.csv).
  **Upload USER-GATED; v16 safe default, v17 one-check, v18 three-checks (portal-novelty the 2 new members).**
- **Genie3-groove switch-HEDGE: rigorous NO-GO (2026-10-12 ~02:20Z).** Full 12-variant two-state MD: 4
  CONFIRMED / 7 NONRESPOND / 1 crashed — BUT the 4 "confirmed" are pose-INSTABILITY artifacts: binder_rmsd
  4.5–7.1 Å, iface_frac 0.44–0.66 (the binder DRIFTS + sheds 35–55% interface over 3 ns; vs v18 des_6126
  stable rmsd ~2.0 / iface 0.84). So the His-Lys "switch" (noisy 2·SEM 1.0–1.3) is on a DETACHING pose, not a
  bound-state repulsion. The af2ig binding wall shows up as MD instability. Genie3-groove can't be confirmed to
  bind OR to switch-while-bound → NOT submitted. v18 (bpair, af2ig pae<10 + MD-stable + clean switch) stands.
  Method lesson: gate switch-MD verdicts on binder_rmsd/iface, not just His-Lys dSwitch.
- **v16/v17 stand regardless;** the two live improvement tracks are the bpair switch-harvest (v18 core, MD
  running) + the Genie3-groove switch-hedge (diversification, queued behind it).

## Open items / decision (USER)

- **Upload decision (USER-gated, nothing uploaded):**
  1. **Default (zero extra work): upload v19** (`7756055f`) — fully portal-verified, 20/20 ≥3/4, 20/20 switchers.
  1b. **STRONGEST (recommended): upload v21** (`f7f246f4`) — 17 LEANER DIVERSE designs (brutal-review cut of cousins/weak-switch + 7 NEW-lineage refills: des_820 [2-cat, TM 0.575] + the fresh-campaign trio des_36133 [3-cat] / des_34125 [2-cat] / des_33012 [2-cat], all MD-direction-confirmed). LOCAL-ONLY. REQUIRES USER portal novelty re-check on the 13 non-byte-identical members (only 4 — c44/H9-14-38/b14/H43 — are byte-identical to portal-confirmed v19; the 5 shell-adapt saN + sh3_sa11 + 7 new-lineage were never portal-checked). Expansion COMPLETE. Dossier docs/SUBMISSION_v21_candidate.md (UPLOAD_CHECKLIST authoritative).
  2. **Stronger on #1, needs portal re-checks: v20** (`6f4b61d0`) — 17/20 (→19/20 w/ round-2) fold-robust,
     binding re-confirmed per design. REQUIRES USER portal novelty re-check on the 9 (soon 11) reseqenced swaps
     (sequence changed ~60-75/80 residues; backbone unchanged → expect ≥3/4, but hard-gated so must be verified).
- **v18/v19 are the safe gated fallbacks.** The fold-robust v20 is the bold #1-axis improvement in flight; the
  honest claim is higher-CONFIDENCE switch geometry across more members, corroborated (not proven) by MD.

## Key invariants (see AGENTS.md + docs/METHODS.md)
- pH-switch = ≥2–3 SEPARATED His on DISTINCT conserved target cations; prove sign per His (PROPKA dpKa<0 AND
  xTB dE>0) AND the two-state MD trigger (dSwitch>2·SEM, >0, dIso>0) on the exact submitted sequence. k_off is
  wet-lab-only; d_iface≈0 in 3ns is EXPECTED (trigger, not release).
- Novelty is a hard gate (≥3/4). Scaffold-guided de-novo = seq-novel + structural-capped → 3/4 (passes).
  Portal-check before locking; USER-gated.
- GPUs: 8× H100 local, run CUDA 0-7 directly; NEVER clobber a foreign job (verify nvidia-smi; VLLM/scGPT
  cycle on 0,1,6,7 — kill only by explicit PID). Per-GPU Triton cache isolation for parallel RFdiffusion.
- Commit source/docs/harness only; `runs/ results/ envs/ third_party/ tools/ weights/` local-only.
