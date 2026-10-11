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
