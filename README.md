# anthropic-adaptyv-2026-tnf

De-novo conditional (pH-switch) binders for **TNF-α** — the Anthropic × Adaptyv 2026 Protein Design
Competition, **Challenge 2**. The designs bind TNF-α at **pH 7.4 and release at pH 6.0** (the
endosomal-recycling direction), with **human↔mouse cross-reactivity**, and the switch is proven in
silico before any wet-lab test.

Built on the reusable conditional-binder toolkit (generalized from the Challenge-1 EGFR project); the
TNF-specific target, epitope, and switch live in `configs/` and `docs/`.

## Objective priority

1. **pH-selectivity** — bind ON at pH 7.4, release OFF at pH 6.0 (highest)
2. **cross-reactivity** — one sequence binds human and mouse TNF-α
3. **affinity** in the ON state

Novelty (distance from known proteins) is a hard gate, not a score — final novelty is the Proteinbase portal.

## The submission

**The pick is priority-dependent** — options differ on re-check burden and on #1-axis depth, not on a provable
domination (switch magnitude is unrankable in-silico, binding magnitude is wet-lab-decided, novelty is portal-gated).
**User portal uploads (2026-10-13) revealed the portal also runs an internal cousin-dedup** ("Similar designs in
this upload"): on v22 it flagged 4 cousins (des_4098_1 ×3, des_2106_2 ×2) → only 16 unique (hits v19/v21 too). The
dedup-clean **v23 and v24 were then both uploaded — ALL 20 PASSED each** → both fully portal-confirmed, 0 re-checks.

**The v27 ladder (2026-10-16 brutal whole-portfolio stress test)** — re-opened the full v26 roster adversarially
("best 20 even if it means rejecting all selected"). A decisive direct-FE sign extension (3 remaining 6-seed
negative-leaners → 12 seeds; the one that firmed → 24) found a **second anti-switch hiding in v26**:
`des_34125_g4_sa15` = **−3.96 kcal, z−4.46, 20/24 seeds negative**, firming monotonically with sampling
(−0.68→−3.26→−3.96) to a magnitude larger than both known anti-switches (b14 −3.88, des_84_3 −2.16). **An external
review (2026-10-17) rightly flagged that this cut overrode the analyzer's own `NEG_LEAN_UNCONFIRMED` label (20/24
misses its ≥23/24 "all-but-one-seed" rule). A falsifiable control extended *both* reference anti-switches to the
same 24-seed depth and settled it: b14 lands 22/24 and des_84_3 19/24, so both *also* fail the ≥23/24 rule — it is a
6-seed heuristic that breaks at high n, where the z-test and distribution-free sign-test govern. By those,
des_34125 (p=7.7e-4) is more one-signed than the already-removed des_84_3 (p=3.3e-3); the cut is vindicated and
v27=18 stands.** All three options below are **#1-pure (0 confirmed anti-switch, 5 confirmed releasers, validator PASS)**; they trade count
against novelty-certainty at one fixed core. A clean 20 (all portal-confirmed AND all releasing) is not
data-supported — two anti-switches had to be removed and the built pool has no novelty-safe non-cousin releaser left.

`runs/push/submission_v27.csv` — **the recommended #1-pure primary** (md5 `f5147388ceb064481315a33f20a85bc7`,
**18 designs**, validator PASS). = **v26 minus `des_34125_g4_sa15`** → **2 anti-switches removed total**
(des_84_3 + des_34125). Every member is portal-confirmed (strict subset of the portal-passed v25), so **zero
re-check, zero dedup, zero novelty risk**. 5 confirmed releasers, 0 confirmed anti-switch. Dossier `docs/SWITCH_FE_findings.md`.

`runs/push/submission_v27d.csv` — **diversity alternative, 19 designs** (md5 `17b45312840f89597a0869f2c49f054d`,
validator PASS). = v27 + `des_243_3_r3_sa2` (NEW lineage, neutral switch, best-in-pool novelty TM 0.808 — more
novel than the already-portal-passed des_732 at 0.842, so strongly inferred to pass). Adds one lineage of wet-lab
hedge at the cost of one portal novelty-check.

`runs/push/submission_v27dd.csv` — **the full-20 option** (md5 `25bf9322eccaa1c5ef7196b825aeb3e5`, 20 designs,
validator PASS). = v27d + `des_1169_2_r3_sa7` (NEW lineage, neutral switch, strong binder, novelty TM 0.934 —
marginal, comparable to the already-passed sh3_sa11 at 0.913). Reaches the requested 20; the soft spot is
des_1169's novelty (two portal-checks needed, des_1169 the riskier). Upload user-gated.

`runs/push/submission_v26.csv` — **superseded by v27** (md5 `73119a41fee4d726a39263229ad3f9b0`,
**19 designs**, validator PASS). = **v25 minus `des_84_3_mh_49_56_76`**, a strict subset of portal-confirmed v25
(every design already portal-PASSED, no re-check, no dedup risk). A whole-portfolio direct-FE **sign** audit
(all 20 v25 members, bound+free, 6→12 seeds — `scripts/fe/analyze_sign.py`) found des_84_3 is a **confirmed
anti-switch** (paired ΔΔG −2.16, z=−4.6, 11/12 seeds negative) — the one v25 member proven to switch the WRONG
direction on the #1 axis, same class as b14 (~half the magnitude, higher-seed confidence). No clean replacement
exists (des_84_3 is the only des_84-lineage design with FE coverage; portal-confirmed non-v25 designs are all
negative-leaning, dedup-cousins, or lack an FE sign), so the #1-correct move is to drop it. v26 holds **5
FE-confirmed releasers and 0 confirmed anti-switch**. Dossier `docs/SWITCH_FE_findings.md`.

`runs/push/submission_v25.csv` — **the 20-design max-diversity alternative** (md5 `41d5e12b4d89a353c44c7c3eb1289064`,
20 designs, validator PASS, **portal-parsed 2026-10-15: ALL 20 PASSED, 0 needing fixes, no cousin/dedup flags**). = v24
with ONE direct-free-energy swap `des_4098_1_b14 → des_4098_1_c44` (same lineage): alchemical ΔΔG (His⁰→His⁺, bound
complex vs free binder) shows **b14 is a reproducible anti-switch** (all 6 seeds negative, z=−6.3) while **c44 is a
confirmed strong releaser** (+8.0 kcal/mol, z=+6.5; c44 = b14 + His63/His66, the release mechanism). Clears every axis —
novelty (c44 portal-confirmed ≥3/4, byte-identical in v19), binding (co-fold ipTM 0.974), cross-react (94% footprint
conserved), dedup (portal 20-unique). **Trade vs v26:** v25 keeps one extra lineage of wet-lab diversity hedge, at the
cost of one slot (des_84_3) now proven to switch the wrong way. Supersedes v24. Dossiers
`docs/SUBMISSION_v25_candidate.md` + `docs/SWITCH_FE_findings.md`.

`runs/push/submission_v24.csv` — **conservative alternative** (md5 `34f3c096bf63dd9ef70d05308c6c3f1d`, 20
designs, validator PASS, **all 20 portal-PASSED 2026-10-13**). v23 with 3 same-lineage swaps to v21's solid ≥2-cat
both-fold saN variants → **8 both-fold members = the #1-axis ordinal ceiling, dedup-clean 20-unique, 0 re-checks.** The
exact earlier-portal-blessed set — pick it if you distrust the in-silico FE swap. Supersedes v21. Dossier `docs/SUBMISSION_v24_candidate.md`.

`runs/push/submission_v23.csv` — **confirmed alternative** (md5 `65bf89bb38f5498ad90522542655eae2`, 20 designs,
validator PASS, **all 20 portal-PASSED**). 5 both-fold, 19 distinct lineages, 0 re-checks; fewer shell-adapt members
than v24 (pick it if you prefer less shell-adaptation over #1-depth). Dossier `docs/SUBMISSION_v23_candidate.md`.

**Not in-silico-provable:** "dominate every metric" — binding magnitude (#3) is wet-lab-decided (Anthropic's own
AP≈0.16 ≈ near-chance in-silico ranking) and switch magnitude is thermodynamically unrankable (confirmed by direct
alchemical FE *and* enhanced-sampling metadynamics — the bound-leg slow mode is interface reorganization, not the
His CV). The provable axes
are maximized: novelty 20/20 confirmed (v23); #1 both-fold ordinal = 8 (v24, the ceiling — a 4-His campaign
confirmed 3-cat/member is the max); diversity 19 lineages. **Switch direction is now audited portfolio-wide by
direct FE** (6→12 seeds, every member): v26 is all-releaser/neutral with **0 confirmed anti-switch**; v25 carries
exactly one (des_84_3, −2.2, z−4.6), which is why v26 drops it.

`runs/push/submission_v22.csv` / `submission_v21.csv` — **SUPERSEDED** (by v23 and v24 respectively); kept as records.

`runs/push/submission_v19.csv` — **20 designs**, validator PASS, **user portal-checked all 20 ≥3/4**, 20/20
carrying a pH-release switch. The **zero-work fallback** — but it carries 4 cousins (→ ~16 unique under the
portal's internal dedup), so it is **dominated by the fully-confirmed v24/v23**; use only as a last resort
(further fallback `submission_v18.csv`).

`runs/push/submission_v20.csv` — the **fold-robust predecessor to v21** (dossier `docs/SUBMISSION_v20_candidate.md`):
v19 with 9 fold-fragile switchers replaced by **shell-adapted fold-ROBUST** variants whose switch-His engages a
cross-conserved TNF cation in **both** independent folders (Boltz∩Chai), lifting the portfolio from 8/20 to
17/20 (→19/20) fold-robust. Binding re-confirmed per design (AF2 initial-guess). Superseded by v21's leaner
diverse cut, but its fold-robust members carry into v21. Every design carries:

- a **pH-release switch** (interface His adjacent to a conserved target Lys → His⁺ repels Lys⁺ at pH 6),
  validated by explicit two-state MD (+ Gln-isostere causal control), PROPKA–Wyman and pKAI linkage, and
  cross-folder (Boltz∩Chai) His–cation agreement for the v20 fold-robust set;
- **cross-reactivity** by conserved-epitope mechanism (89–95% human↔mouse interface identity, cations K/R in both);
- a confident binder (AF2 initial-guess), **0 Cys**, expression-screened, and sequence-novel; **structural
  novelty (the HARD portal ≥3/4 gate) must be re-checked on every reseqenced design before upload** — v20's 9
  shell-adapted swaps, and in v21 the 10 members whose sequence is not byte-identical to portal-confirmed v19.

Upload to the Proteinbase portal is **user-gated**. CSVs under `runs/` are local-only (gitignored).

## Pipeline (what actually ran)

Most submitted designs come from the **β-pairing RFdiffusion long-shot** at the conserved Region-I edge →
ProteinMPNN (0-Cys) → **AF2 initial-guess binding** (pae_int<10, the load-bearing filter) → partial-diffusion
amplification of the cross-reactive new shapes → install **≥2 SEPARATED His, each repelling a DISTINCT
conserved target cation** (K6/K93/K60…) → **two-state MD** switch validation (+ Gln-isostere control) → PROPKA
+ pKAI linkage → conserved-epitope cross-reactivity → novelty (foldseek TM + mmseqs; portal is authoritative).
**v20 adds LigandMPNN shell-adaptation** (`scripts/shelladapt_shell.py`): repack the His microenvironment with
the His FIXED so BOTH independent folders (Boltz∩Chai) agree the His engages its conserved cation (fold-robust),
re-confirming binding per design. `scripts/README.md` maps every stage to its scripts.

## Environment & compute

**uv** only for the core env; heavy design models install into isolated per-tool uv venvs under `envs/`
(see `docs/TOOLS.md`). Compute is local multi-GPU (run CUDA directly, one process per GPU). All scripts
resolve their repo root dynamically — no hardcoded paths.

## Quickstart

```bash
./init.sh                 # uv sync + compile + pytest + validator self-check (the gate)
uv run pytest -q          # tests (schema, length, numbering, liabilities)
uv run python scripts/validate_submission.py <designs.csv> --track 3   # submission gate
```

## Layout

- `AGENTS.md` / `CLAUDE.md` — the working contract and invariants (start here)
- `progress.md` — current status and next step
- `configs/` — the target/epitope/switch contract (`target.yaml`), filter thresholds (`filters.yaml`), tool matrix (`uplifting.yaml`)
- `scripts/` — the pipeline implementation (`scripts/README.md` is the manifest: active TNF pipeline vs reusable utils vs legacy EGFR code)
- `src/adaptyv_tnf/` — the stage skeleton and the submission contract/validator
- `docs/PLAN_TNF.md` — the campaign plan and outcome · `docs/METHODS.md` — the validation playbook · `docs/TOOLS.md` — the tool stack
- `runs/submission/` — the final portfolio + methods (local-only; gitignored)
