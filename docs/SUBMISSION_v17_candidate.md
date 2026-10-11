# Submission v17 candidate — 3rd de-novo new-lineage switcher (USER-PORTAL-GATED)

Updated 2026-10-11 (SESSION 11, continued). Deadline **2026-10-11 23:59 AoE = 2026-10-12 11:59 UTC**.

## TL;DR

v17 = **v16 with one swap**: the weakest redundant MD-NONRESPOND member `b4m2_A57H` (+0.22, below
threshold) is replaced by a **3rd de-novo new-lineage MD-CONFIRMED switcher**,
`des_6025_4_T20_pd_64_1_mh_34_50_52_74`. Portfolio two-state-MD trigger **17/20 → 18/20 CONFIRMED**;
new de-novo lineages **2 → 3** (`des_138_7`, `des_4082_8`, `des_6025`). The `b4m2` lineage keeps its
CONFIRMED member `b4m2_r1_w5_3`, so no lineage is lost.

- **v17 CSV:** `runs/push/submission_v17.csv` (md5 **`9357b6905ed81179da5137a16dd9dd5a`**), validator **PASS**
  (20 designs, 0 errors, 0-Cys, unique names+seqs, all `single_chain`).
- **v16 CSV (unchanged, the SAFE default):** `runs/push/submission_v16.csv` (md5 `a524d77fa7a63ef33b28346ce3a7d87d`).

## The swap-in: des_6025_4_T20_pd_64_1_mh_34_50_52_74

`SAAVAAQQAEFNKLLAERVEAVAKELRAKGVEVHVEVSENGGVALIARSHDHAAALEAVDLMKKMTAEVIEELHAAAAAA` (80 aa, 0 Cys)

Fully re-validated **on the exact His-variant sequence** (not inherited from the parent — the
re-validate-the-mutated-artifact invariant):

| Axis | Metric | Value | Verdict |
|---|---|---|---|
| **#1 pH-switch** | two-state MD d_repel (HIP−HIE), 4 reps × 3 ns | **+0.372** (2·SEM ≈ 0.27) | **CONFIRMED** (marginal) |
|  | isostere control (HIP−GLN) | +0.308 | positive (ablation holds) |
|  | `selective_release` | `true` | — |
|  | His→Lys targets | His 34/50/52/74 → Lys 93/107/60/6 | separated multi-His |
| **#2 cross-react** | human iptm (Boltz-2, re-fold) | **0.962** | — |
|  | mouse iptm (Boltz-2, re-fold) | **0.938** | — |
|  | min(h,m) | **0.938** | **PASS** (≥0.70) |
| **expression** | Protein-Sol scaled solubility | 0.87 | PASS (≥0.50) |
| **novelty seq-id** | mmseqs vs PDB seqres | **NO PDB hit (novel)** | passes seq criterion |
| **novelty struct** | foldseek TMalign max-TM vs PDB | **0.945** | high — see risk below |

Wyman linkage (static, human co-fold): **−1.114** (strongest of all 15 des_6025 His-variants);
per-His dpKa·linkage [(34,−0.82),(50,−0.60),(52,−0.36),(74,−1.16)].

## The one open risk: portal structural-novelty (USER-GATED)

The local foldseek structural max-TM is **0.945**, *higher* than the v16 de-novo swap-ins' 0.86–0.93
band. Per the SESSION-11 novelty calibration, scaffold-guided de-novo shapes score portal **3/4** because
they are **seq-novel** (no PDB hit — confirmed here) and that criterion plus two others carry them past the
single structural-TM criterion. So des_6025 is **likely 3/4**, but at a higher structural TM it is a
**higher-risk, portal-UNVERIFIED 3/4**, unlike the v16 members which are all user-portal-verified. Novelty
is a **hard per-design gate (DQ if <3/4)**, and only the portal (account `j_v_v_07`) is authoritative.

## Recommendation (decision is the user's)

1. **Zero extra work → upload v16** (`a524d77f`): fully portal-verified, 20/20 at ≥3/4 novelty, 17/20
   MD-confirmed switch. The safe, complete submission.
2. **If one portal novelty check can be spared → check des_6025.**
   - **≥3/4 → upload v17** (`9357b690`): strictly better on the #1 axis (18/20 confirmed switch), +1
     new de-novo lineage, strong cross-react (0.938). The b4m2 lineage is retained via `b4m2_r1_w5_3`.
   - **<3/4 → stay on v16.**

**Upload is USER-GATED either way** — do not auto-upload. CSVs under `runs/` are gitignored (local-only).

## Reproduction

- Re-fold (cross-react): `scripts/run_boltz_tnf.sh <1-row csv> <target_txt> <tag> <gpu>` vs
  `runs/partial_denovo2/funnel/xreact/human_tnf.txt` (human) and `runs/buried2/mouse_tnf_mature.txt` (mouse).
- Two-state MD verdict: `runs/partial_denovo2/funnel/md/des_6025_4_T20_pd_64_1_mh_34_50_52_74.json`
  (`d_repel` +0.372, `selective_release` true).
- Seq novelty: `scripts/novelty_check.sh <csv>` (fail-loud, positive control).
- Struct novelty: `tools/foldseek_dl/foldseek/bin/foldseek easy-search <binder.pdb>
  runs/novelty_struct/pdbdb/pdb ... --alignment-type 1 --format-output "...,alntmscore,..."`.
