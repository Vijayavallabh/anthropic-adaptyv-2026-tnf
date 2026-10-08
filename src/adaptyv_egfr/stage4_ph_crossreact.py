"""Stage 4 - pH-selectivity + cross-reactivity specialization (THE decisive stage).

This stage decides the two highest-priority objectives.

pH-selectivity (objective 1, highest):
  * Repack/minimise each pose at pH 6.5 and pH 7.4 with EXPLICIT His protonation toggling.
  * Score interface ddG at each pH (PyRosetta fa_elec dddG_elec, or pHDock).
  * RANK by pH-selectivity margin = ddG(7.4) - ddG(6.5): binding ON at 6.5, OFF at 7.4.
  * Tumor-ON SIGN: >=2-3 interface His forming salt bridges to Domain III Asp/Glu ONLY
    when protonated at 6.5. Do NOT copy endosomal-recycling (bind-OFF-at-low-pH) His motifs.
  * Two-tier pKa QC: PROPKA3 fast triage -> PypKa + JCIM-2025 His predictor on the shortlist;
    confirm each interface His titrates ~6.3-7.0.

cross-reactivity (objective 2):
  * Fold/score every finalist against BOTH human (6ARU chain A) and the mouse EGFR Domain III
    model; keep only dual binders whose energetic contacts land on the conserved anchor set.

A small PILOT (protonated vs neutral scoring on a positive control) must validate the switch
before scaling - fold models cannot confirm pH dependence.

CPU-bound: parallelise over cores while stages 1-3 use the GPUs.
Requires PyRosetta + PypKa (see docs/TOOLS.md). PROPKA3 (`propka`) is in the core env.
"""
from __future__ import annotations

PH_LOW = 6.5   # tumor microenvironment - binding ON
PH_HIGH = 7.4  # normal tissue - binding OFF
HIS_PKA_TARGET = (6.3, 7.0)
MIN_INTERFACE_HIS = 2


def ph_selectivity_margin(pose, *, ph_low: float = PH_LOW, ph_high: float = PH_HIGH) -> float:
    """ddG(pH_high) - ddG(pH_low); larger = sharper tumor-ON switch."""
    raise NotImplementedError(
        "stage4.ph_selectivity_margin needs PyRosetta (fa_elec / pHDock). See docs/TOOLS.md."
    )


def dual_species_gate(design_seq: str) -> bool:
    """True iff the design is predicted to bind BOTH human and mouse EGFR Domain III."""
    raise NotImplementedError("stage4.dual_species_gate needs the folding stack + mouse model.")


def main() -> None:  # pragma: no cover
    print(__doc__)


if __name__ == "__main__":
    main()
