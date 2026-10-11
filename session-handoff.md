# Session handoff — TNF Challenge 2

Updated 2026-10-11 (SESSION 11). Deadline **2026-10-11 23:59 AoE = 2026-10-12 11:59 UTC**.

## Current best submission (validated, portal-novelty-verified)

- **`runs/push/submission_v16.csv`** (md5 **`a524d77fa7a63ef33b28346ce3a7d87d`**) — 20 designs, validator
  **PASS**, all swap-ins portal-novelty **3/4**. **Recommended upload.**
- **Fallbacks (fully portal-verified):** `submission_v14.csv` (md5 `7994cf46f0a55921a986cf73cee5e13d`),
  `submission_v13.csv` (`4c682b80…`), `submission_v12.csv` (`e34ca96d…`).
- **Upload is USER-GATED** — do not auto-upload. CSVs under `runs/` (gitignored). Account `j_v_v_07`.
- Portal CSV format: `name,sequence,molecule_class`; de-novo minis `single_chain`; unique names AND seqs.

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
  4,5,6,7** (PID in `$CLAUDE_JOB_DIR/tmp/gen_regI_main.log`, N=800, ~182/800 so far, ~70min left). mmbu still
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
- **v16/v17 stand regardless;** this thrust could only have ADDED a Genie3 hedge — it floored at binding.

## Open items / decision (USER)

- **Upload decision (USER-gated, nothing uploaded):**
  1. **Default (zero extra work): upload v16** (`a524d77f`) — fully portal-verified, 20/20 ≥3/4 novelty.
  2. **If one portal novelty check can be spared: check des_6025.** ≥3/4 → **upload v17** (`9357b690`,
     stronger on #1: 18/20 confirmed + a 3rd de-novo lineage). <3/4 → stay on v16.
- De-novo long shot is **played out** — strong yield captured in v16 (+1 marginal in v17). Diminishing returns.

## Key invariants (see AGENTS.md + docs/METHODS.md)
- pH-switch = ≥2–3 SEPARATED His on DISTINCT conserved target cations; prove sign per His (PROPKA dpKa<0 AND
  xTB dE>0) AND the two-state MD trigger (dSwitch>2·SEM, >0, dIso>0) on the exact submitted sequence. k_off is
  wet-lab-only; d_iface≈0 in 3ns is EXPECTED (trigger, not release).
- Novelty is a hard gate (≥3/4). Scaffold-guided de-novo = seq-novel + structural-capped → 3/4 (passes).
  Portal-check before locking; USER-gated.
- GPUs: 8× H100 local, run CUDA 0-7 directly; NEVER clobber a foreign job (verify nvidia-smi; VLLM/scGPT
  cycle on 0,1,6,7 — kill only by explicit PID). Per-GPU Triton cache isolation for parallel RFdiffusion.
- Commit source/docs/harness only; `runs/ results/ envs/ third_party/ tools/ weights/` local-only.
