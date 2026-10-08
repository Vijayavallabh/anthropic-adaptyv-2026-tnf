"""Stage 2 - sequence design (inverse folding).

For RFdiffusion / BoltzGen backbones, design ~8 sequences/backbone with
LigandMPNN / ProteinMPNN using the SolubleMPNN weights on the solvent-exposed binder
face (maximise expressibility, ~95% achievable historically). Use fixed-position / bias
flags to PLANT the interface histidines that pair with EGFR Domain III Asp/Glu residues
(the pH switch, stage 4). BindCraft already emits sequences (optionally re-surface).

Trivial on a single H100; batch thousands.

Primary tools: Anthropic uplifting kits proteinmpnn / caliby (SolubleCaliby); fallback ligandmpnn.
See docs/TOOLS.md and configs/uplifting.yaml.
"""
from __future__ import annotations

from pathlib import Path

DEFAULT_SEQS_PER_BACKBONE = 8
SOLUBLE_WEIGHTS = "soluble_mpnn"  # SolubleMPNN weights shipped with ligandmpnn


def design(backbone_pdb: Path, *, n_seqs: int = DEFAULT_SEQS_PER_BACKBONE,
           fixed_positions: list[int] | None = None, his_positions: list[int] | None = None) -> list[str]:
    """Return designed sequences for one backbone (planting His at `his_positions`)."""
    raise NotImplementedError(
        "stage2.design needs ligandmpnn (torch). Install via `./init.sh --tools`."
    )


def main() -> None:  # pragma: no cover
    print(__doc__)


if __name__ == "__main__":
    main()
