"""Stage 1 - de novo backbone / binder generation.

Generate diverse backbones against the CONSERVED Domain III hotspots on 6ARU chain A
(the cross-reactive, pH-sensor face). Throughput, not model sharding, is the bottleneck:
run one generator process per GPU across all 8x H100.

Generators (open-source; installed in isolated envs by `init.sh --tools`, see docs/TOOLS.md):
  * BindCraft   - primary AF2-hallucination one-shot generator (best historical de-novo EGFR yield)
  * BoltzGen    - newest fully-open (MIT) universal all-atom generator (orthogonal diversity)
  * RFdiffusion - battle-tested targeted backbone diffusion (v1) + ProteinMPNN + AF2 initial-guess
  * RFantibody  - nanobody / scFv / Fab formats

Hotspots come from data/epitope/epitope_map.json (conserved anchors 408/432/433/436/464/489),
biasing His placement toward the Domain III Asp/Glu partners for the pH switch (stage 4).

Primary tools: Anthropic uplifting kits colabdesign / boltzgen / rfdiffusion1
(scripts/setup_uplifting.sh; configs/uplifting.yaml). NOT in the core uv env. See docs/TOOLS.md.
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGET_PDB = REPO / "data" / "targets" / "6ARU.pdb"
TARGET_CHAIN = "A"
GENERATORS = ("bindcraft", "boltzgen", "rfdiffusion", "rfantibody")


def run(generator: str, *, n_designs: int, out_dir: Path, gpu: int, hotspots: list[int]) -> None:
    """Dispatch one generator on one GPU. Implemented per-tool in the isolated envs."""
    raise NotImplementedError(
        f"stage1.run('{generator}') needs the generator tool stack. "
        "Install with `./init.sh --tools` and wire the per-tool driver here (see docs/TOOLS.md)."
    )


def main() -> None:  # pragma: no cover - orchestration entry
    print(__doc__)
    print(f"target: {TARGET_PDB} chain {TARGET_CHAIN}; generators: {GENERATORS}")


if __name__ == "__main__":
    main()
