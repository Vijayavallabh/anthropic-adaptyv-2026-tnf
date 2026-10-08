"""Stage 6 - novelty gate, LLM shortlist, rank & CSV export.

  * Enforce the de-novo hard filter: sequence (MMseqs2 / %identity-to-PDB) AND structural
    (Foldseek TMscore) diversity vs known proteins; deduplicate.
  * Shortlist with a Claude-generated multi-metric ranking policy over
    ipTM/ipSAE/pH-margin/normalized-ESM/interface proxies (arXiv 2608.20755), mirroring how
    the organizers use Claude to select Tracks 2/3 designs.
  * Rank by the objective priority: pH-margin > cross-reactivity > affinity.
  * Cap per track (Track 1: 40; Tracks 2/3: 20) and export the CSV through the validated
    builder in ``submission.py`` (the definition-of-done gate must pass before upload).

The novelty gate and LLM shortlist need external tools (Foldseek/MMseqs2 + Claude);
``export_submission`` below is fully implemented and reuses the core validator.
"""
from __future__ import annotations

from pathlib import Path

from .submission import TRACK_CAPS, Report, validate_rows, write_submission

# Ranking key: objective priority order (see __init__.OBJECTIVES).
RANK_KEYS = ("ph_margin", "cross_reactive_score", "affinity_score")


def novelty_gate(sequences: list[str]) -> list[bool]:
    """True per sequence if sufficiently novel vs known proteins (Foldseek + MMseqs2)."""
    raise NotImplementedError("stage6.novelty_gate needs foldseek + mmseqs2. See docs/TOOLS.md.")


def llm_shortlist(candidates: list[dict], *, policy: str | None = None) -> list[dict]:
    """Claude multi-metric shortlist (arXiv 2608.20755). Needs ANTHROPIC_API_KEY."""
    raise NotImplementedError("stage6.llm_shortlist needs the anthropic SDK + ANTHROPIC_API_KEY.")


def export_submission(ranked_rows: list[dict], out_path, *, track: int = 3, extra_columns=None) -> Report:
    """Cap to the track limit, write the ranked CSV, and validate it. Returns the Report.

    ``ranked_rows`` must already be sorted best-first and contain name/sequence/molecule_class.
    Raises ValueError if the resulting CSV fails the submission contract.
    """
    cap = TRACK_CAPS.get(track)
    rows = list(ranked_rows)[:cap] if cap else list(ranked_rows)
    write_submission(rows, out_path, extra_columns=extra_columns)
    rep = validate_rows(rows, track=track)
    if not rep.ok:
        raise ValueError(
            f"exported CSV fails the submission contract ({rep.n_errors} errors); "
            f"run: uv run python scripts/validate_submission.py {out_path} --track {track}"
        )
    return rep


def main() -> None:  # pragma: no cover
    print(__doc__)


if __name__ == "__main__":
    main()
