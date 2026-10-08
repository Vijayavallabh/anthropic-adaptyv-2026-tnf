# Tool stack and installation

**Environment policy:** everything is driven by **uv**. The *core* env
(`uv sync`) is CPU-only and light (target prep, numbering, liability scan,
submission build and validation, tests). Heavy GPU models each install into
their **own isolated uv venv**.

## Primary provider: inference-optimized model kits

H100-tuned, inference-optimized kits are the primary tool provider. Each kit is
self-contained (`run.sh`, a stock manifest, an `h100.env`), Apache-2.0, cloned
to `third_party/` (gitignored), with weights downloaded on demand into
`weights/` (SHA-256 verified).

- Modes: `off` | `exact` | `fast` (default) | `big` (`--n_gpu P` splits across GPUs).
- The stage-to-kit mapping lives in `configs/uplifting.yaml`.

### Install a kit (uv)

```bash
scripts/setup_uplifting.sh <kit>            # venv + the exact install command, no weights
scripts/setup_uplifting.sh <kit> --weights  # also downloads weights (several GB)
```

Under the hood: `uv venv --seed --managed-python --python 3.11 envs/uplifting-<kit>`,
then the kit's `run.sh install`, then `run.sh check`. Persist the JIT cache with
`export MODEL_OPT_JIT_ROOT=$PWD/.jit`. Python is pinned per kit in its stock
manifest (mostly 3.11, some 3.12).

### Kits per pipeline stage

| Stage | Primary kits | Notes |
|---|---|---|
| 1 generate | `colabdesign` (BindCraft), `boltzgen`, `rfdiffusion1` | opt: `rfdiffusion3`, `pxdesign`, `genie3`, `complexa`; one process per GPU |
| 2 sequence | `proteinmpnn`, `caliby` (soluble) | `proteinmpnn` needs an explicit `--mode`; `caliby` is py3.12 |
| 3 fold/score | `af2ig` (AF2 initial-guess), `boltz2`, `chai1` | multi-GPU via `--n_gpu`. Real protein-protein affinity = Rosetta InterfaceAnalyzer `dG_separated` (PyRosetta) + PRODIGY, NOT a co-folder's ligand affinity head |
| 6 rank | `esmc`, `e1` | length-normalized LM scoring |

## Conditional-switch scoring (the decisive engine)

The trigger is real physics and the folders are blind to it (see
`docs/METHODS.md` section 3). For a pH switch, two orthogonal models must agree
on the exact submitted sequence:

- **Thermodynamic pKa linkage** (`scripts/obj1_switch.py`, PROPKA3 from the core
  env): integrate the proton-uptake difference across the pH window from the
  complex/target/binder pKas, restricted to the sensor residues, with buried
  charge artifacts clamped. The same script also *installs* the switch (scan
  binder positions near the sensor, place a carboxylate, FastRelax, keep best).
  Needs PyRosetta (academic license) plus PROPKA3.
- **Explicit-state MD occupancy** (`scripts/md_pH_switch.py`, OpenMM on GPU):
  build the sensor in each protonation state, run both, and confirm the switch
  only when the ON-state interaction occupancy exceeds the OFF-state by a clear
  margin. Sharded across GPUs by `scripts/run_md_campaign.py`.

A design is confirmed only when both legs agree. A positive static score alone
proves direction, not selectivity. Independent extra legs where available: a PB-class
pKa shift (pKAI), a burial-aware pKa (JustHISpKa), a fair quantum switch energy
(GFN2-xTB with an implicit solvent, for cation-pi where fixed-charge MD is
biased), and a two-state MM-GBSA binding delta.

Dead ends: Rosetta `-pH_mode` does not titrate on ref2015; constant-pH MD is a
setup NO-GO; no structure-based predictor beats a pKa-6.5 null for a lone
histidine, so trust cross-method sign-consensus.

## Standalone pip fallbacks (if not using the kits)

| Package | Stage | Install | Notes |
|---|---|---|---|
| `boltz` | 3,4 | `uv pip install boltz` | co-fold ipTM, MIT. Its affinity head is small-molecule only, out-of-distribution for protein mini-binders |
| `chai_lab` | 3 | `uv pip install chai_lab` | Apache-2.0, no gating |
| `ligandmpnn` | 2 | `uv pip install ligandmpnn` | soluble-MPNN weights bundled |
| `fair-esm` | 6 | `uv pip install fair-esm` | length-normalized pseudolikelihood |
| `foldseek` + `mmseqs2` | 6 | bioconda / prebuilt binaries (non-uv) | de-novo novelty gate |

## GPU strategy

Local multi-GPU, run CUDA directly (no SSH). Throughput, not model size, is the
bottleneck: one process per GPU (the sweep scripts shard across GPUs), or
`--mode big --n_gpu P` for large inputs. Queue heavy sweeps into low-utilisation
windows, and cap per-process CPU threads (OMP_NUM_THREADS near cores/n_gpu) so
concurrent jobs do not oversubscribe the box.

## A fold is not a switch

AlphaFold, Boltz, and Chai are trigger-agnostic, so a good fold or ipTM is not
evidence of the switch. Confirm selectivity only via the two-state engine above,
on the exact submitted sequence. When engineering a stronger switch, re-validate
the mutated artifact end to end: a high linkage score that fails MD, or that
breaks binding, is an artifact, not a switch.
