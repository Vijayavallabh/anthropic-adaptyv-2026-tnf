"""Stage 3 - fold + interface scoring (high-throughput filter).

1. Cheap first filter with AF2 initial-guess (dl_binder_design): interface pAE, ipTM,
   pLDDT, binder RMSD. Use ipTM >= 0.8 and pae_interaction < 7 as SANITY GATES ONLY -
   the 2024/25 Adaptyv rounds showed iPAE/ipTM/pLDDT were near-random for predicting
   binding, so never rank on them alone.
2. Re-score survivors with Boltz-2 (affinity head -> objective 3) and cross-check the
   top set with Chai-1 (Apache-2.0, no gating). AF3 only if weights already granted.

Parallelise across the 8 H100 (one model per GPU).

CAUTION: AlphaFold/Boltz/Chai are pH-AGNOSTIC. A good fold metric is NOT evidence of the
pH switch or cross-reactivity - those are decided in stage 4.

Primary tools: Anthropic uplifting kits af2ig / boltz2 / chai1 (--mode fast, --n_gpu for multi-GPU).
See docs/TOOLS.md and configs/uplifting.yaml.
"""
from __future__ import annotations

# Sanity gates only (not selection criteria).
IPTM_MIN = 0.80
PAE_INTERACTION_MAX = 7.0


def fold_score(design_seq: str, target_pdb, *, predictor: str = "af2_initial_guess") -> dict:
    """Return interface metrics for one design against the target."""
    raise NotImplementedError(
        "stage3.fold_score needs the folding stack (af2_initial_guess / boltz / chai_lab). "
        "Install via `./init.sh --tools`."
    )


def main() -> None:  # pragma: no cover
    print(__doc__)


if __name__ == "__main__":
    main()
