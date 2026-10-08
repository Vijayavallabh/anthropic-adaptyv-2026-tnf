#!/usr/bin/env python3
"""WAVE 19 portfolio selection: order a submission CSV to maximize best-of (max-of-K) success.

Adaptyv scores best-of: standing is the single best measured design and the lab synthesizes only the
top few (T ~ 3-8) in rank order, so the right objective for the submitted ORDER is not "the K highest
point estimates" but P(at least one of the tested top-T expresses AND binds AND switches), under
CORRELATED failure: designs from the same partial-diffusion lineage / fold share a failure mode, so a
second cousin adds almost no new coverage. That is a submodular max-coverage problem; greedily covering
decorrelated groups formalizes the "diversity not cousins" rule and pushes a high-upside design from a
fresh group above a slightly-stronger cousin of an already-covered group.

Model (transparent, adjustable priors; the RELATIVE ordering and the correlation structure carry the
decision, not the absolute probabilities):
  p_success(i) = p_express * p_bind * p_switch
    p_express : from developability expr_risk (CLEAN/MODERATE/SEVERE) -- cell-free expression is a real
                attrition filter, so an unexpressible charge-blob contributes little even if it would switch.
    p_bind    : from the binding consensus (co-fold iPTM; falls back to AF2-ig interface PAE).
    p_switch  : from the MD selectivity (d_occ, penalized by neutral-pH leak occ_HIE).
  group(i)   = partial-diffusion lineage (dominant shared fold/expression failure mode).
  Full within-group failure correlation: a group's coverage = the best p_success of its SELECTED members,
  so cousins beyond the first add ~0. P(>=1 success in a set) = 1 - prod_over_groups(1 - coverage_g).
  The knockout control (p_switch ~ 0 by design) is forced to the last slot.

Read-only; writes results/<stem>_portfolio.csv (the recommended order) and prints the coverage curve and
a diff vs the input order. No GPU, no runs/ access.
"""
import argparse, csv, os, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METRICS_PRIMARY = os.path.join(REPO, "results", "metrics_final_w19.csv")
METRICS_FALLBACK = os.path.join(REPO, "results", "metrics_20.csv")

# --- adjustable priors (documented heuristics, not calibrated truths) ---
P_EXPRESS = {"CLEAN": 0.90, "MODERATE": 0.70, "SEVERE": 0.45, "": 0.75}  # by developability expr_risk


def ff(x):
    try:
        if x is None or (isinstance(x, str) and x.strip() == ""):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def clip(x, lo, hi):
    return max(lo, min(hi, x))


def p_bind(bind_consensus, af2_pae):
    """Binding prior: prefer the co-fold consensus iptm; else map AF2-ig interface PAE (lower=better)."""
    if bind_consensus is not None and bind_consensus > 0:
        return clip(0.45 + 0.45 * bind_consensus, 0.45, 0.92)
    if af2_pae is not None:
        return clip(0.92 - 0.045 * af2_pae, 0.30, 0.92)   # pae 6 -> 0.65, pae 12 -> 0.38
    return 0.55


def p_switch(d_occ, occ_HIE):
    """Switch prior: high d_occ is a real pH switch; neutral-pH leak (occ_HIE) erodes it."""
    if d_occ is None:
        return 0.15
    leak = 0.0 if occ_HIE is None else 0.6 * max(0.0, occ_HIE - 0.2)
    return clip(d_occ - leak, 0.08, 0.95)


def _first(row, *keys):
    for k in keys:
        if k in row and row[k] not in (None, ""):
            return row[k]
    return None


def load_metrics():
    """name -> {d_occ, occ_HIE, bind_consensus, af2_pae, lineage, tier}, merging w19 over the 20-row table."""
    out = {}
    for path in (METRICS_FALLBACK, METRICS_PRIMARY):   # primary loaded second so it wins
        if not os.path.exists(path):
            continue
        for row in csv.DictReader(open(path)):
            chai = ff(_first(row, "chai_iptm_human", "chai_human"))
            boltz = ff(_first(row, "boltz2_iptm", "boltz"))
            bc = ff(_first(row, "bind_consensus"))
            if bc is None and (chai is not None or boltz is not None):
                bc = min([v for v in (chai, boltz) if v is not None])
            rec = {"d_occ": ff(_first(row, "d_occ")), "occ_HIE": ff(_first(row, "occ_HIE")),
                   "bind_consensus": bc, "af2_pae": ff(_first(row, "af2_pae", "af2ig_iface_pae", "pae_interaction")),
                   "lineage": (_first(row, "lineage", "lin") or ""), "tier": (_first(row, "md_tier", "tier") or "")}
            for key_col in ("new_name", "name"):
                k = (row.get(key_col) or "").strip()
                if k:
                    out[k] = rec
    return out


def load_dev(path):
    """name -> expr_risk from one developability_scan.py CSV (the one matching the submission)."""
    out = {}
    if path and os.path.exists(path):
        for row in csv.DictReader(open(path)):
            k = (row.get("name") or "").strip()
            if k and row.get("expr_risk"):
                out[k] = row["expr_risk"].strip().upper()
    return out


def trailing_index(name):
    tail = ""
    for ch in reversed(name or ""):
        if ch.isdigit():
            tail = ch + tail
        else:
            break
    return int(tail) if tail else None


def coverage(cov):
    """P(>=1 success) given a dict group -> best p_success among selected members of that group."""
    prod = 1.0
    for p in cov.values():
        prod *= (1.0 - p)
    return 1.0 - prod


def portfolio_order(rows, metrics, dev, group_of=None):
    """Greedy submodular max-coverage ordering. Returns (ordered_rows, p_success_by_name, group_by_name)."""
    def grp(name):
        if group_of and name in group_of:
            return group_of[name]
        return (metrics.get(name, {}).get("lineage") or name)  # lineage = shared fold/failure group

    psucc, group = {}, {}
    knockout = None
    idxs = [trailing_index(r.get("name")) for r in rows]
    max_idx = max((i for i in idxs if i is not None), default=None)
    for r in rows:
        name = (r.get("name") or "").strip()
        m = metrics.get(name, {})
        risk = dev.get(name, "")
        is_ko = "control" in (m.get("tier") or "").lower() or (trailing_index(name) == max_idx and max_idx is not None
                                                               and p_switch(m.get("d_occ"), m.get("occ_HIE")) < 0.2)
        ps = P_EXPRESS.get(risk, P_EXPRESS[""]) * p_bind(m.get("bind_consensus"), m.get("af2_pae")) * \
            p_switch(m.get("d_occ"), m.get("occ_HIE"))
        psucc[name] = round(ps, 4)
        group[name] = grp(name)
        if is_ko:
            knockout = r

    pool = [r for r in rows if r is not knockout]
    selected, cov = [], {}
    while pool:
        best, best_gain = None, -1.0
        for r in pool:
            name = (r.get("name") or "").strip()
            g, p = group[name], psucc[name]
            trial = dict(cov)
            trial[g] = max(trial.get(g, 0.0), p)
            gain = coverage(trial) - coverage(cov)
            # tie-break: higher raw p_success first (stable, deterministic)
            if gain > best_gain + 1e-12 or (abs(gain - best_gain) <= 1e-12 and best is not None
                                            and p > psucc[(best.get("name") or "").strip()]):
                best, best_gain = r, gain
        name = (best.get("name") or "").strip()
        cov[group[name]] = max(cov.get(group[name], 0.0), psucc[name])
        selected.append(best)
        pool.remove(best)
    if knockout is not None:
        selected.append(knockout)
    return selected, psucc, group


def main():
    ap = argparse.ArgumentParser(description="Order a submission by best-of portfolio coverage.")
    ap.add_argument("submission_csv", nargs="?")
    ap.add_argument("--out")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        _selftest(); return
    if not a.submission_csv:
        ap.error("submission_csv required (or --selftest)")

    stem = os.path.splitext(os.path.basename(a.submission_csv))[0]
    # resolve the developability file robustly: try the exact stem, then the stem with ordering/variant
    # suffixes stripped (so re-running on <stem>_FINAL / _portfolio / _ranked still finds the base dev
    # file instead of silently falling back to a mismatched developability.csv -> flat expression prior).
    dev = {}
    stems = [stem]
    for suf in ("_FINAL", "_portfolioOrder", "_portfolio", "_ranked", "_evt", "_kdpp"):
        if stem.endswith(suf):
            stems.append(stem[: -len(suf)])
    for s in stems:
        dev = load_dev(os.path.join(REPO, "results", "developability_%s.csv" % s))
        if dev:
            break
    dev = dev or load_dev(os.path.join(REPO, "results", "developability.csv"))
    metrics = load_metrics()
    with open(a.submission_csv, newline="") as fh:
        reader = csv.DictReader(fh); fields = reader.fieldnames; rows = list(reader)
    ordered, psucc, group = portfolio_order(rows, metrics, dev)

    out = a.out or (os.path.splitext(a.submission_csv)[0] + "_portfolio.csv")
    if os.path.abspath(out) == os.path.abspath(os.path.join(REPO, "results", "designs_20.csv")):
        sys.exit("refusing to overwrite results/designs_20.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(ordered)

    # coverage curve vs T, and the ordering
    cov = {}
    print(f"{'rank':>4}  {'design':<14}{'group':<16}{'p_succ':>7}  P(>=1 success in top-T)")
    for t, r in enumerate(ordered, 1):
        name = (r.get("name") or "").strip()
        cov[group[name]] = max(cov.get(group[name], 0.0), psucc[name])
        mark = "  <= tested band" if t in (3, 5, 8) else ""
        print(f"{t:>4}  {name:<14}{group[name]:<16}{psucc[name]:>7.3f}  {coverage(cov):.3f}{mark}")
    print(f"\nwrote {os.path.relpath(out, REPO)}  (knockout forced last)")


def _selftest():
    # two strong cousins in group g1, one slightly-weaker design in a fresh group g2, plus a knockout.
    rows = [{"name": "egfr_phsw_21", "sequence": "AAAA"},
            {"name": "egfr_phsw_22", "sequence": "BBBB"},   # same group as 21
            {"name": "egfr_phsw_40", "sequence": "CCCC"},   # fresh group, slightly weaker
            {"name": "egfr_phsw_99", "sequence": "DDDQ"}]   # knockout
    metrics = {
        "egfr_phsw_21": {"d_occ": 0.95, "occ_HIE": 0.03, "bind_consensus": 0.85, "af2_pae": 6.0, "lineage": "L1", "tier": "B"},
        "egfr_phsw_22": {"d_occ": 0.93, "occ_HIE": 0.05, "bind_consensus": 0.85, "af2_pae": 6.0, "lineage": "L1", "tier": "B"},
        "egfr_phsw_40": {"d_occ": 0.80, "occ_HIE": 0.10, "bind_consensus": 0.80, "af2_pae": 6.5, "lineage": "L2", "tier": "B"},
        "egfr_phsw_99": {"d_occ": 0.05, "occ_HIE": 0.9, "bind_consensus": 0.80, "af2_pae": 6.0, "lineage": "L1", "tier": "C_control"}}
    ordered, psucc, group = portfolio_order(rows, metrics, {})
    names = [r["name"] for r in ordered]
    # best of group L1 first, THEN the fresh group L2 (coverage), THEN the L1 cousin, knockout last
    assert names[0] == "egfr_phsw_21", names
    assert names[1] == "egfr_phsw_40", ("fresh group should outrank the cousin for coverage", names)
    assert names[2] == "egfr_phsw_22", names
    assert names[-1] == "egfr_phsw_99", ("knockout last", names)
    # coverage is monotone non-decreasing
    cov, series = {}, []
    for r in ordered:
        n = r["name"]; cov[group[n]] = max(cov.get(group[n], 0.0), psucc[n]); series.append(coverage(cov))
    assert all(series[i] <= series[i + 1] + 1e-9 for i in range(len(series) - 1)), series
    # adding the cousin (rank 3) adds ~0 coverage vs rank 2 (same group already covered)
    assert series[2] - series[1] < 1e-6, ("cousin should add ~0 coverage", series)
    print("selftest ok")


if __name__ == "__main__":
    main()
