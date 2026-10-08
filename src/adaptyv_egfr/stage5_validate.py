"""Stage 5 - antibody validation & developability.

For nanobody / scFv / Fab designs:
  * number with anarcism (core env) / riot_na; confirm a valid Ig V-domain, CDRs in
    canonical positions, intact frameworks, and VL chain_type ('K'/'L') matching the
    declared fab_kappa/fab_lambda (this satisfies the ANARCI hard filter).
  * graft designed CDRs onto known humanized frameworks:
        VHH  -> h-NbBcII10FGLA
        scFv/Fab -> IGHV3-23/DP-47 with IGKV1-39 (kappa) or an IGLV (lambda)
  * run the CDR liability scan below; optionally TAP + BioPhi OASis (see docs/TOOLS.md).

`scan_liabilities` is implemented here (pure regex, core env) and is reused by the ranking
stage to penalise designs with developability liabilities in the CDRs.
"""
from __future__ import annotations

import re

# Sequence-liability motifs (screen especially within CDRs).
_NGLYC = re.compile(r"N[^P][ST]")   # N-linked glycosylation sequon
_DEAMIDATION = re.compile(r"N[GS]")  # NG / NS deamidation hotspots
_ISOMERIZATION = re.compile(r"D[GST]")  # DG / DS / DT isomerization hotspots


def scan_liabilities(seq: str, cdr_spans: list[tuple[int, int]] | None = None) -> dict:
    """Flag common developability liabilities in a single chain (1-based positions).

    cdr_spans (optional) restricts reporting to CDR regions for the ``in_cdr`` summary.
    """
    seq = (seq or "").upper()

    def positions(rx):
        return [m.start() + 1 for m in rx.finditer(seq)]

    cys = [i + 1 for i, c in enumerate(seq) if c == "C"]
    hits = {
        "n_glycosylation": positions(_NGLYC),
        "deamidation_NG_NS": positions(_DEAMIDATION),
        "isomerization_DG_DS_DT": positions(_ISOMERIZATION),
        "cysteines": cys,
        "unpaired_cys": (len(cys) % 2 == 1),
        "met_oxidation": [i + 1 for i, c in enumerate(seq) if c == "M"],
        "trp_oxidation": [i + 1 for i, c in enumerate(seq) if c == "W"],
    }
    if cdr_spans:
        def in_cdr(pos):
            return any(lo <= pos <= hi for lo, hi in cdr_spans)
        hits["in_cdr"] = {
            k: [p for p in v if in_cdr(p)]
            for k, v in hits.items()
            if isinstance(v, list)
        }
    return hits


def liability_count(seq: str) -> int:
    """Total liability motif count (for ranking penalties)."""
    h = scan_liabilities(seq)
    return (
        len(h["n_glycosylation"])
        + len(h["deamidation_NG_NS"])
        + len(h["isomerization_DG_DS_DT"])
        + int(h["unpaired_cys"])
    )


def main() -> None:  # pragma: no cover
    print(__doc__)


if __name__ == "__main__":
    main()
