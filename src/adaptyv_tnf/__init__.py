"""adaptyv_tnf - pipeline + submission tooling for the Anthropic x Adaptyv 2026
TNF-alpha conditional-binder challenge (Challenge 2).

Objective priority (this order IS the scoring priority):
    1. pH-selective binding        bind @ pH 7.4, release @ pH 6.0       (highest)
    2. human/mouse cross-reactive  the same sequence binds both
    3. affinity to TNF-alpha in the ON (pH 7.4) state

See README.md and docs/METHODS.md for the method and pipeline.
"""
from __future__ import annotations

__version__ = "0.1.0"

# Objective keys in priority order (see stage4/stage6 ranking).
OBJECTIVES = ("ph_selective", "cross_reactive", "human_affinity")


def main() -> None:
    """Entry point for the ``adaptyv-tnf`` console script."""
    import textwrap

    print(
        textwrap.dedent(
            f"""
            adaptyv-tnf v{__version__}
            Anthropic x Adaptyv 2026 - TNF-alpha conditional-binder challenge (Challenge 2)

            Objective priority : pH-selectivity > cross-reactivity > affinity
            Start here         : README.md
            Validate a CSV     : uv run python scripts/validate_submission.py <file.csv> --track 3
            """
        ).strip()
    )
