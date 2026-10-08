#!/usr/bin/env python3
"""WAVE 19 re-ranker: re-order a submission CSV (name,sequence,molecule_class) best-first
by the stated objective priority (pH-selectivity > cross-reactivity > affinity), force the
E66Q knockout control to the last row, and flag (never silently fix) the case where the
knockout's parent does not land in the top tier.

Objective priority (report ranking section, "Rank by MD selectivity and keep the knockout
in slot twenty"):
  1. pH-selectivity: MD d_occ descending, then occ_HIE ascending (lower neutral-pH leak first).
  2. tie-break: cross-reactivity (smaller xreact_delta_pae), then binding (lower af2
     pae_interaction, then more negative MM-GBSA dd_bind).
  3. the knockout control is always the last row.
  4. the knockout's parent (egfr_phsw_22, or the design the knockout is a one-residue mutant
     of) must be in the top tier; if it is not, that is a WARNING to stderr, not a reorder.

Metrics source: results/metrics_final_w19.csv for the new designs, falling back to
results/metrics_20.csv. Column names differ between the two files, so both schemas are read.
Designs with no d_occ are ranked by d_occ alone and sort to the bottom of the real designs;
every field is None-safe. This script only reads metrics and writes a new ranked CSV; it runs
no GPU job, touches no runs/ path, and never overwrites results/designs_20.csv.
"""
import argparse
import csv
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# metrics source precedence: new-design table first, then the shipped 20-row table
METRICS_PRIMARY = os.path.join(REPO, "results", "metrics_final_w19.csv")
METRICS_FALLBACK = os.path.join(REPO, "results", "metrics_20.csv")

# number of real designs that count as the top tier for the parent-placement guarantee
DEFAULT_TOP_TIER = 6

INF = float("inf")

# expression-risk tie-break + lead guard, from developability_scan.py expr_risk. Best-of scoring tests
# the top few designs first, so an unexpressible charge-blob at the lead wastes the best slot; prefer
# CLEAN leads and warn (do not silently override selectivity) if a SEVERE design sits in the top tier.
RISK_RANK = {"CLEAN": 0, "MODERATE": 1, "SEVERE": 2}
DEV_CANON = os.path.join(REPO, "results", "developability.csv")


def ff(x):
    """Parse a float, returning None for blank or non-numeric cells."""
    try:
        if x is None:
            return None
        x = x.strip() if isinstance(x, str) else x
        if x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def _first(row, *keys):
    """First non-empty value among the given column names (handles both file schemas)."""
    for k in keys:
        if k in row and row[k] not in (None, ""):
            return row[k]
    return None


def normalize(row):
    """Map one metrics row (either schema) to the fields the ranker uses."""
    return {
        "d_occ": ff(_first(row, "d_occ")),
        "occ_HIE": ff(_first(row, "occ_HIE")),
        "xreact_delta_pae": ff(_first(row, "xreact_delta_pae")),
        "af2_pae": ff(_first(row, "af2_pae", "af2ig_iface_pae", "pae_interaction")),
        "mmgbsa_dd_bind": ff(_first(row, "mmgbsa_dd_bind", "dd_bind_kcal")),
        "tier": (_first(row, "md_tier", "tier") or ""),
    }


def load_metrics(path):
    """Build a name -> normalized-metrics lookup, keyed on both new_name and name so a
    submission name (egfr_phsw_NN) finds its row in either schema."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            norm = normalize(row)
            for key_col in ("new_name", "name"):
                k = (row.get(key_col) or "").strip()
                if k:
                    out.setdefault(k, norm)
    return out


def _trailing_index(name):
    """Trailing integer of a name like egfr_phsw_32 -> 32, else None."""
    tail = ""
    for ch in reversed(name or ""):
        if ch.isdigit():
            tail = ch + tail
        else:
            break
    return int(tail) if tail else None


def _hamming(a, b):
    """Residue differences between two equal-length strings, or None if not comparable."""
    if a is None or b is None or len(a) != len(b):
        return None
    return sum(1 for x, y in zip(a, b) if x != y)


def _is_control_tier(tier):
    return "control" in (tier or "").lower()


def rank_rows(submission_rows, metrics_by_name, dev_by_name=None, top_tier_size=DEFAULT_TOP_TIER):
    """Core re-ranking logic, file-free so it can be unit-tested in memory.

    submission_rows: list of dicts, each with at least 'name' and 'sequence'.
    metrics_by_name: dict name -> normalized metrics (keys: d_occ, occ_HIE,
                     xreact_delta_pae, af2_pae, mmgbsa_dd_bind, tier).
    Returns (ordered_rows, warnings) where ordered_rows is the re-ordered list and
    warnings is a list of strings (the caller decides where to print them).
    """
    warnings = []
    dev_by_name = dev_by_name or {}
    if not submission_rows:
        return [], warnings

    def met(row):
        return metrics_by_name.get((row.get("name") or "").strip(), {})

    # knockout detection: tier marked control, OR the maximum trailing index. Either flags it.
    indices = [_trailing_index(r.get("name")) for r in submission_rows]
    max_idx = max((i for i in indices if i is not None), default=None)

    def is_knockout(row):
        if _is_control_tier(met(row).get("tier")):
            return True
        idx = _trailing_index(row.get("name"))
        return idx is not None and idx == max_idx

    knockouts = [r for r in submission_rows if is_knockout(r)]
    designs = [r for r in submission_rows if not is_knockout(r)]

    def sort_key(row):
        m = met(row)
        d_occ = m.get("d_occ")
        occ_hie = m.get("occ_HIE")
        xreact = m.get("xreact_delta_pae")
        pae = m.get("af2_pae")
        ddbind = m.get("mmgbsa_dd_bind")
        return (
            (-d_occ) if d_occ is not None else INF,   # d_occ descending, missing sorts last
            occ_hie if occ_hie is not None else INF,   # occ_HIE ascending, missing sorts last
            xreact if xreact is not None else INF,     # cross-reactivity, smaller first
            pae if pae is not None else INF,           # interface PAE, lower first
            ddbind if ddbind is not None else INF,     # MM-GBSA, more negative first
            RISK_RANK.get(dev_by_name.get((row.get("name") or "").strip()), 0),  # final tie-break: CLEAN leads
        )

    designs.sort(key=sort_key)   # stable: equal keys keep input order
    ordered = designs + knockouts

    # expression-risk lead guard: flag (never auto-override selectivity) SEVERE designs in the top tier
    severe_lead = [(d.get("name") or "?").strip() for d in designs[:top_tier_size]
                   if dev_by_name.get((d.get("name") or "").strip()) == "SEVERE"]
    if severe_lead:
        warnings.append(
            "WARNING: expression-risk SEVERE design(s) in the top %d by selectivity: %s. These are "
            "charge-blob or low-complexity sequences at risk of not expressing in the cell-free assay; "
            "best-of scoring tests the top few first, so weigh whether to lead with them."
            % (top_tier_size, ", ".join(severe_lead)))

    # parent-placement guarantee for each knockout (normally exactly one)
    for ko in knockouts:
        ko_seq = ko.get("sequence")
        cands = [d for d in designs if _hamming(d.get("sequence"), ko_seq) == 1]
        parent = None
        if cands:
            # prefer the literal parent name if present among the one-residue mutants
            parent = next(
                (c for c in cands if (c.get("name") or "").strip() == "egfr_phsw_22"),
                cands[0],
            )
        else:
            parent = next(
                (d for d in designs if (d.get("name") or "").strip() == "egfr_phsw_22"),
                None,
            )
        ko_name = (ko.get("name") or "?").strip()
        if parent is None:
            warnings.append(
                "WARNING: knockout %s has no identifiable parent (no one-residue-mutant match "
                "and no egfr_phsw_22 in the set); parent-placement guarantee not checked." % ko_name
            )
            continue
        pos = designs.index(parent)   # 0-indexed rank among real designs
        if pos >= top_tier_size:
            warnings.append(
                "WARNING: knockout %s parent %s ranked %d of %d real designs, outside the top "
                "%d tier; a knockout whose parent is not tested is dead weight. Not reordered "
                "automatically -- fix the parent's selectivity inputs or drop the knockout."
                % (ko_name, (parent.get("name") or "?").strip(), pos + 1, len(designs), top_tier_size)
            )
    return ordered, warnings


def resolve_metrics_sources(explicit):
    """Ordered, de-duplicated list of existing metrics files, earlier = higher precedence.
    metrics_20 (the full 20-row table) is ALWAYS appended as a fallback so the 8 locked designs,
    which are absent from metrics_final_w19, still contribute their d_occ to the ranking."""
    want = [explicit if explicit else METRICS_PRIMARY, METRICS_FALLBACK]
    seen, out = set(), []
    for p in want:
        if p and os.path.exists(p) and os.path.abspath(p) not in seen:
            seen.add(os.path.abspath(p)); out.append(p)
    if not out:
        sys.exit("no metrics source found (looked for %s then %s)" % (METRICS_PRIMARY, METRICS_FALLBACK))
    return out


def load_metrics_merged(paths):
    """Merge several metrics files into one name -> metrics map; earlier paths win per design."""
    merged = {}
    for p in reversed(paths):   # load lowest-precedence first; earlier paths overwrite
        merged.update(load_metrics(p))
    return merged


def load_dev(path):
    """name -> expr_risk ('CLEAN'/'MODERATE'/'SEVERE') from a developability_scan.py CSV, or {}."""
    out = {}
    if path and os.path.exists(path):
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                k = (row.get("name") or "").strip()
                if k and row.get("expr_risk"):
                    out[k] = row["expr_risk"].strip().upper()
    return out


def run_cli(submission_csv, metrics_path=None, out_path=None, top_tier_size=DEFAULT_TOP_TIER, dev_path=None):
    sources = resolve_metrics_sources(metrics_path)
    metrics_by_name = load_metrics_merged(sources)
    # expression-risk: prefer an explicit --dev, then developability_<submission stem>.csv, then the canonical file
    stem = os.path.splitext(os.path.basename(submission_csv))[0]
    dev_by_name = (load_dev(dev_path)
                   or load_dev(os.path.join(REPO, "results", "developability_%s.csv" % stem))
                   or load_dev(DEV_CANON))

    with open(submission_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not fieldnames:
        sys.exit("submission CSV has no header: %s" % submission_csv)

    if out_path is None:
        base, ext = os.path.splitext(submission_csv)
        out_path = base + "_ranked" + (ext or ".csv")
    # never clobber the canonical design file
    if os.path.abspath(out_path) == os.path.abspath(os.path.join(REPO, "results", "designs_20.csv")):
        sys.exit("refusing to overwrite results/designs_20.csv; choose another --out path")

    ordered, warnings = rank_rows(rows, metrics_by_name, dev_by_name=dev_by_name, top_tier_size=top_tier_size)
    for w in warnings:
        print(w, file=sys.stderr)

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(ordered)

    print("re-ranked %d rows (metrics: %s) -> %s"
          % (len(ordered), ", ".join(os.path.relpath(p, REPO) for p in sources), out_path))
    return out_path


def _selftest():
    """Synthetic in-memory checks: ordering, occ_HIE tie-break, knockout-last, parent warning."""
    import io

    # four real designs + one knockout. d19 has the highest d_occ but is the knockout (control
    # tier) so it must still land last. d_tieA/d_tieB share d_occ and split on occ_HIE.
    sub = [
        {"name": "egfr_phsw_02", "sequence": "AAAA", "molecule_class": "single_chain"},
        {"name": "egfr_phsw_tieA", "sequence": "BBBB", "molecule_class": "single_chain"},
        {"name": "egfr_phsw_tieB", "sequence": "CCCC", "molecule_class": "single_chain"},
        {"name": "egfr_phsw_22", "sequence": "DDDE", "molecule_class": "single_chain"},
        {"name": "egfr_phsw_99", "sequence": "DDDQ", "molecule_class": "single_chain"},  # knockout: E->Q of _22
    ]
    met = {
        "egfr_phsw_02": {"d_occ": 0.30, "occ_HIE": None, "xreact_delta_pae": None,
                         "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "A_locked"},
        "egfr_phsw_tieA": {"d_occ": 0.80, "occ_HIE": 0.30, "xreact_delta_pae": None,
                           "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "B"},
        "egfr_phsw_tieB": {"d_occ": 0.80, "occ_HIE": 0.05, "xreact_delta_pae": None,
                           "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "B"},
        "egfr_phsw_22": {"d_occ": 0.95, "occ_HIE": 0.04, "xreact_delta_pae": None,
                         "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "B"},
        "egfr_phsw_99": {"d_occ": 0.99, "occ_HIE": 0.00, "xreact_delta_pae": None,
                         "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "C_control"},
    }
    ordered, warnings = rank_rows(sub, met, top_tier_size=6)
    names = [r["name"] for r in ordered]
    assert names == ["egfr_phsw_22", "egfr_phsw_tieB", "egfr_phsw_tieA",
                     "egfr_phsw_02", "egfr_phsw_99"], names
    assert names[-1] == "egfr_phsw_99", "knockout must be last"
    # tieB before tieA: same d_occ, lower occ_HIE wins
    assert names.index("egfr_phsw_tieB") < names.index("egfr_phsw_tieA")
    assert warnings == [], warnings   # parent _22 is rank 1, inside top tier

    # parent _22 is rank 1, so even a top_tier_size of 1 leaves it inside -> no warning
    _, warn2 = rank_rows(sub, met, top_tier_size=1)
    assert warn2 == [], warn2
    # now push the parent's d_occ down so it falls behind the ties -> warning fires
    _, warn3 = rank_rows(
        [sub[1], sub[2], sub[3], sub[4]],   # drop _02; order _22 behind the two ties by d_occ
        {**met, "egfr_phsw_22": {**met["egfr_phsw_22"], "d_occ": 0.10}},
        top_tier_size=2,
    )
    assert any("parent" in w for w in warn3), warn3

    # knockout detected by max index alone when no metrics row exists for it
    sub2 = [
        {"name": "egfr_phsw_21", "sequence": "AAAA"},
        {"name": "egfr_phsw_40", "sequence": "AAAQ"},   # max index, no metrics -> knockout
    ]
    met2 = {"egfr_phsw_21": {"d_occ": 0.5, "occ_HIE": 0.1, "xreact_delta_pae": None,
                             "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "B"}}
    ordered2, _ = rank_rows(sub2, met2)
    assert [r["name"] for r in ordered2][-1] == "egfr_phsw_40"

    # None-safe: a design missing d_occ sorts below one that has it (the _88 row is the
    # knockout by max index, so the two compared designs are both real designs)
    sub3 = [
        {"name": "egfr_phsw_03", "sequence": "AAAA"},   # no d_occ
        {"name": "egfr_phsw_04", "sequence": "BBBB"},   # has d_occ
        {"name": "egfr_phsw_88", "sequence": "CCCC"},   # knockout (max index)
    ]
    met3 = {"egfr_phsw_04": {"d_occ": 0.2, "occ_HIE": None, "xreact_delta_pae": None,
                             "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "A"}}
    ordered3, _ = rank_rows(sub3, met3)
    assert [r["name"] for r in ordered3] == ["egfr_phsw_04", "egfr_phsw_03", "egfr_phsw_88"]

    # expr_risk: a SEVERE design in the top tier triggers a warning, selectivity order unchanged
    _, warn4 = rank_rows(sub, met, dev_by_name={"egfr_phsw_22": "SEVERE"}, top_tier_size=6)
    assert any("SEVERE" in w for w in warn4), warn4
    # expr_risk final tie-break: with all objective axes equal, CLEAN sorts before MODERATE
    subT = [{"name": "egfr_phsw_61", "sequence": "AAAA"},
            {"name": "egfr_phsw_62", "sequence": "BBBB"},
            {"name": "egfr_phsw_90", "sequence": "CCCQ"}]   # max index -> knockout, absorbs the last slot
    metT = {"egfr_phsw_61": {"d_occ": 0.8, "occ_HIE": 0.1, "xreact_delta_pae": None, "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "B"},
            "egfr_phsw_62": {"d_occ": 0.8, "occ_HIE": 0.1, "xreact_delta_pae": None, "af2_pae": None, "mmgbsa_dd_bind": None, "tier": "B"}}
    orderedT, _ = rank_rows(subT, metT, dev_by_name={"egfr_phsw_62": "MODERATE", "egfr_phsw_61": "CLEAN"}, top_tier_size=6)
    assert [r["name"] for r in orderedT] == ["egfr_phsw_61", "egfr_phsw_62", "egfr_phsw_90"], [r["name"] for r in orderedT]

    _ = io  # keep import referenced for anyone extending the test with CSV fixtures
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser(description="Re-rank a submission CSV best-first by objective priority.")
    ap.add_argument("submission_csv", nargs="?", help="input CSV (name,sequence,molecule_class)")
    ap.add_argument("--metrics", help="metrics source CSV (default: metrics_final_w19.csv, else metrics_20.csv)")
    ap.add_argument("--out", help="output ranked CSV (default: <input>_ranked.csv)")
    ap.add_argument("--top-tier-size", type=int, default=DEFAULT_TOP_TIER,
                    help="number of lead designs that count as the top tier (default %d)" % DEFAULT_TOP_TIER)
    ap.add_argument("--dev", help="developability CSV with expr_risk (default: developability_<stem>.csv, else developability.csv)")
    ap.add_argument("--selftest", action="store_true", help="run in-memory unit checks and exit")
    args = ap.parse_args()

    if args.selftest:
        _selftest()
        return
    if not args.submission_csv:
        ap.error("submission_csv is required (or pass --selftest)")
    run_cli(args.submission_csv, metrics_path=args.metrics, out_path=args.out,
            top_tier_size=args.top_tier_size, dev_path=args.dev)


if __name__ == "__main__":
    main()
