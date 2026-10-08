"""adaptyv_egfr - pipeline + submission tooling for the Anthropic x Adaptyv 2026
EGFR conditional-binder challenge (Challenge 1).

Objective priority (this order IS the scoring priority):
    1. pH-selective binding        bind @ pH 6.5, no binding @ pH 7.4   (highest)
    2. human/mouse cross-reactive  the same sequence binds both
    3. high affinity to human EGFR extracellular region

See AGENTS.md for the harness contract and feature_list.json for pipeline state.
"""
from __future__ import annotations

__version__ = "0.1.0"

# Objective keys in priority order (see stage4/stage6 ranking).
OBJECTIVES = ("ph_selective", "cross_reactive", "human_affinity")


def main() -> None:
    """Entry point for the ``adaptyv-egfr`` console script."""
    import textwrap

    print(
        textwrap.dedent(
            f"""
            adaptyv-egfr v{__version__}
            Anthropic x Adaptyv 2026 - EGFR conditional-binder challenge (Challenge 1)

            Objective priority : pH-selectivity > cross-reactivity > affinity
            Start here         : AGENTS.md
            Pipeline state     : feature_list.json  (stages 0-6)
            Validate a CSV     : uv run python scripts/validate_submission.py <file.csv> --track 3
            """
        ).strip()
    )
