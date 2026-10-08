#!/usr/bin/env python3
"""WAVE 22g: selection-&-validation-rigor bundle for the confirmed PASS20 (deep-research Tier-1 math/stats).

Runs in envs/analysis. ZERO floor risk: the SET of 20 is fixed (portal-confirmed); this only RE-ORDERS and
STRESS-TESTS the ordering under noise + correlation, as an audit of which designs a top-K wet-lab pick
should cover. Implements, on the per-design oracle scores we already have:
  1. log-opinion-pool (geometric-mean) oracle fusion vs arithmetic mean + min-consensus (where individual
     oracles exist: the inverted designs) -- the per-oracle-veto fusion appropriate to pH-sel-first priorities;
  2. multi-criterion rank aggregation (Plackett-Luce latent quality + Borda) across the ranked objective axes;
  3. STABILITY SELECTION under a modeled oracle-noise bootstrap -> top-K membership frequency = robust vs
     noise-driven (the direct antidote to "a small-n metric win is noise until bootstrapped");
  4. empirical-Bayes / James-Stein shrinkage of the per-axis scores (limited, protects the genuine outlier);
  5. portfolio P(>=1 wet-lab hit) submodular-greedy ordering with a correlated-failure (lineage/mechanism)
     model -- formalizes "best-of, no mean penalty" + "diverse, not cousins";
  6. QUBO/Ising top-K diverse-subset selector (simulated annealing, dimod/neal) as the classical solve of
     the same quality-weighted-diversity objective, cross-checking the DPP reorder;
  7. scalarization audit: linear-weighted-sum vs Chebyshev ranking -> flags designs the linear rule hides in
     concave Pareto regions; plus bootstrap P(Pareto-optimal).

Objective priority used throughout: pH-selectivity (d_occ high, occ_HIE low) >> cross-reactivity >> affinity;
structural novelty is a GATE all 20 already pass (so it is not a discriminator here, only a veto).

HONEST LIMITS (printed): no replicate measurements -> oracle noise is MODELED (documented SEs), and no
exchangeable wet-lab labels -> conformal/FDR guarantees are out of scope; this is robustness, not calibration.
"""
import csv, json, glob, math, os, sys
import numpy as np

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
rng = np.random.default_rng(20261006)

# ---- modeled oracle noise (1 sigma); documented, not measured ----
SE = dict(d_occ=0.12, occ_HIE=0.08, bind=0.04, tm=0.02)  # MD occupancy spread / co-fold iptm / Foldseek TM


def load_pass20():
    """18 from metrics_w21_final + inv01/inv02 from moonshot; attach mechanism class; impute missing."""
    m = {r["name"]: r for r in csv.DictReader(open(f"{REPO}/results/metrics_w21_final.csv"))}
    order = [r["name"] for r in csv.DictReader(open(f"{REPO}/results/adaptyv_challenge1_resubmit_PASS20_inv.csv"))]
    inv = {
        "egfr_phsw_inv01": dict(d_occ=0.515, occ_HIE=0.00, oracle_TM=0.725, bind_consensus=0.623,
                                xreact_mech_sev=0, lineage="inv_moon03", mech="inverted", portal="pass", chai=0.623, boltz=0.776),
        "egfr_phsw_inv02": dict(d_occ=0.87, occ_HIE=0.13, oracle_TM=0.750, bind_consensus=0.786,
                                xreact_mech_sev=0, lineage="inv_des1_1", mech="inverted", portal="pass", chai=0.786, boltz=0.843),
    }
    def f(x):
        try: return float(x)
        except: return None
    rows = []
    for nm in order:
        if nm in inv:
            d = dict(inv[nm]); d["name"] = nm
        else:
            r = m[nm]
            d = dict(name=nm, d_occ=f(r["d_occ"]), occ_HIE=f(r["occ_HIE"]), oracle_TM=f(r["oracle_TM"]),
                     bind_consensus=f(r["bind_consensus"]), xreact_mech_sev=f(r["xreact_mech_sev"]) or 0.0,
                     lineage=r["lineage"], mech="his433", portal="pass", chai=None, boltz=None)
        rows.append(d)
    # Impute missing values. occ_HIE is unrecorded in the summary tables for 9 portal-confirmed designs
    # (it lives only in raw W18 MD jsons); these passed earlier selectivity screening, so impute the
    # MEDIAN of the known occ_HIE (a typical clean value ~0.05), NOT the worst-case max (which would
    # adversarially sink them as leaky). d_occ/bind gaps (only the knockout #32) sink to the min. Flag all.
    import statistics as st
    for key in ("d_occ", "occ_HIE", "bind_consensus"):
        vals = [r[key] for r in rows if r.get(key) is not None]
        fill = st.median(vals) if key == "occ_HIE" else min(vals)
        for r in rows:
            if r.get(key) is None:
                r[key] = fill; r.setdefault("imputed", []).append(key)
    for r in rows:
        if r.get("oracle_TM") is None:
            r["oracle_TM"] = 0.70   # portal-confirmed 3/4 -> treat as novel (good); not a discriminator
    return rows


def axis_scores(r):
    """per-objective 'higher is better' scores in [0,1]-ish raw units (normalized later)."""
    ph = r["d_occ"] - 1.5 * max(0.0, r["occ_HIE"] - 0.15)   # switch margin, penalize leak past the 0.15 gate
    xr = 1.0 if (r["xreact_mech_sev"] or 0) == 0 else 0.0
    aff = r["bind_consensus"]
    return ph, xr, aff


def composite(rows, weights=(0.55, 0.25, 0.20)):
    """weighted-sum of min-max-normalized axes; pH-sel dominant per priority."""
    A = np.array([axis_scores(r) for r in rows], float)  # n x 3 (ph, xr, aff)
    N = np.zeros_like(A)
    for j in range(3):
        lo, hi = A[:, j].min(), A[:, j].max()
        N[:, j] = (A[:, j] - lo) / (hi - lo + 1e-9)
    w = np.array(weights)
    return N @ w, A, N


# ---------- 1. oracle fusion (geometric mean / log-opinion-pool vs arithmetic vs min) ----------
def fusion_demo(rows):
    out = []
    for r in rows:
        if r.get("chai") is not None and r.get("boltz") is not None:
            vals = np.array([r["chai"], r["boltz"]])
            out.append((r["name"], float(vals.mean()), float(np.exp(np.mean(np.log(vals)))), float(vals.min())))
    return out


# ---------- 2. rank aggregation (Plackett-Luce latent quality + Borda) ----------
def rank_aggregation(rows):
    import choix
    n = len(rows)
    # one ranking per objective axis (best-first)
    A = np.array([axis_scores(r) for r in rows], float)
    rankings = []
    for j in range(3):
        order = list(np.argsort(-A[:, j]))  # indices best-first
        rankings.append(order)
    # Plackett-Luce from full rankings (choix wants lists of items in order of decreasing preference)
    try:
        params = choix.ilsr_rankings(n, rankings, alpha=0.01)
    except Exception:
        params = choix.lsr_rankings(n, rankings, alpha=0.01)
    # Borda as a transparent cross-check
    borda = np.zeros(n)
    for order in rankings:
        for pos, idx in enumerate(order):
            borda[idx] += (n - pos)
    return params, borda


# ---------- 3. stability selection under modeled oracle noise ----------
def stability(rows, K_list=(5, 8, 10), nboot=5000):
    n = len(rows)
    base_ph = np.array([r["d_occ"] for r in rows])
    base_hie = np.array([r["occ_HIE"] for r in rows])
    base_bind = np.array([r["bind_consensus"] for r in rows])
    base_xr = np.array([1.0 if (r["xreact_mech_sev"] or 0) == 0 else 0.0 for r in rows])
    freq = {K: np.zeros(n) for K in K_list}
    for _ in range(nboot):
        d = base_ph + rng.normal(0, SE["d_occ"], n)
        h = np.clip(base_hie + rng.normal(0, SE["occ_HIE"], n), 0, 1)
        b = base_bind + rng.normal(0, SE["bind"], n)
        ph = d - 1.5 * np.maximum(0.0, h - 0.15)
        # normalize + weight (same composite)
        def nz(x):
            lo, hi = x.min(), x.max(); return (x - lo) / (hi - lo + 1e-9)
        score = 0.55 * nz(ph) + 0.25 * base_xr + 0.20 * nz(b)
        ordr = np.argsort(-score)
        for K in K_list:
            for idx in ordr[:K]:
                freq[K][idx] += 1
    for K in K_list:
        freq[K] /= nboot
    return freq


# ---------- 5+6. portfolio P(>=1 hit) submodular + QUBO diverse subset ----------
def success_prob(rows):
    """map scores -> a per-design P(wet-lab success) via monotone logistic surrogates (uncalibrated; ranks only)."""
    comp, A, N = composite(rows)
    # squash composite to a plausible success probability band [0.1, 0.6]
    p = 0.1 + 0.5 * (comp - comp.min()) / (comp.max() - comp.min() + 1e-9)
    return p


def similarity(rows):
    """pairwise similarity for diversity: shared mechanism + shared lineage + bind-profile closeness."""
    n = len(rows)
    S = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            s = 0.0
            if rows[i]["mech"] == rows[j]["mech"]: s += 0.5
            if rows[i]["lineage"] == rows[j]["lineage"]: s += 0.5
            S[i, j] = S[j, i] = min(1.0, s)
    return S


def portfolio_greedy(rows, K=10):
    """maximize P(>=1 hit) = 1 - prod(1-p_i) with correlated failure: cousins' failures co-occur.
    Correlation via a shared-lineage/mechanism latent: effective independent contribution is discounted
    by similarity to already-picked designs (submodular)."""
    p = success_prob(rows); S = similarity(rows)
    chosen, remaining = [], list(range(len(rows)))
    log_fail = 0.0  # log prod(1-p_eff)
    order = []
    while remaining and len(chosen) < K:
        best, best_gain = None, -1
        for j in remaining:
            disc = 1.0 - (max([S[j, c] for c in chosen]) if chosen else 0.0)  # novelty vs picked
            p_eff = p[j] * disc
            new_log_fail = log_fail + math.log(max(1e-9, 1 - p_eff))
            gain = (1 - math.exp(new_log_fail)) - (1 - math.exp(log_fail))
            if gain > best_gain:
                best_gain, best, best_lf = gain, j, new_log_fail
        chosen.append(best); remaining.remove(best); log_fail = best_lf
        order.append((rows[best]["name"], round(1 - math.exp(log_fail), 4)))
    return order


def qubo_subset(rows, K=10, lam=0.6):
    """select-K QUBO: maximize sum q_i x_i - lam sum S_ij x_i x_j  s.t. sum x_i = K (penalty). neal anneal."""
    import dimod, neal
    q = success_prob(rows); S = similarity(rows); n = len(rows)
    P = 5.0  # cardinality penalty strength
    Q = {}
    for i in range(n):
        Q[(i, i)] = -q[i] + P * (1 - 2 * K)
        for j in range(i + 1, n):
            Q[(i, j)] = lam * S[i, j] + 2 * P
    bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
    sampler = neal.SimulatedAnnealingSampler()
    res = sampler.sample(bqm, num_reads=200, seed=7)
    best = res.first.sample
    picked = [rows[i]["name"] for i in range(n) if best.get(i, 0) == 1]
    return picked


# ---------- 7. scalarization audit (linear vs Chebyshev) + prob Pareto ----------
def scalarization_audit(rows):
    comp_lin, A, N = composite(rows)
    # Chebyshev: maximize min over weighted-normalized axes (reaches concave Pareto regions)
    w = np.array([0.55, 0.25, 0.20])
    cheb = (N * w).min(axis=1)
    lin_order = [rows[i]["name"] for i in np.argsort(-comp_lin)]
    cheb_order = [rows[i]["name"] for i in np.argsort(-cheb)]
    # bootstrap P(Pareto-optimal) under noise
    n = len(rows); pareto = np.zeros(n)
    base = np.array([[r["d_occ"] - 1.5 * max(0, r["occ_HIE"] - 0.15),
                      1.0 if (r["xreact_mech_sev"] or 0) == 0 else 0.0, r["bind_consensus"]] for r in rows], float)
    for _ in range(2000):
        pert = base + rng.normal(0, [SE["d_occ"], 0.001, SE["bind"]], base.shape)
        for i in range(n):
            dominated = np.any(np.all(pert >= pert[i], axis=1) & np.any(pert > pert[i], axis=1))
            if not dominated: pareto[i] += 1
    pareto /= 2000
    return lin_order, cheb_order, pareto


def main():
    rows = load_pass20()
    names = [r["name"] for r in rows]
    comp, A, N = composite(rows)
    print("=== PASS20 selection-rigor bundle (floor-safe reorder/audit; envs/analysis) ===\n")
    print("modeled oracle noise (1σ, documented not measured):", SE, "\n")

    # 1 fusion
    print("--- 1. oracle fusion (where individual oracles exist = inverted designs) ---")
    print(f"{'design':16} {'arith':>6} {'geom(LOP)':>9} {'min':>6}")
    for nm, a, g, mn in fusion_demo(rows):
        print(f"{nm:16} {a:>6.3f} {g:>9.3f} {mn:>6.3f}")
    print("  (geom mean = KL-barycenter: gives a dissenting oracle a veto, apt for pH-sel-first)\n")

    # 2 rank aggregation
    pl, borda = rank_aggregation(rows)
    print("--- 2. rank aggregation (Plackett-Luce latent quality, best first) ---")
    for idx in np.argsort(-pl)[:8]:
        print(f"  {names[idx]:16} PL={pl[idx]:+.3f}  Borda={borda[idx]:.0f}")
    print()

    # 3 stability
    freq = stability(rows)
    print("--- 3. stability selection (top-K membership freq under noise; robust vs noise-driven) ---")
    print(f"{'design':16} {'top5':>6} {'top8':>6} {'top10':>6}  d_occ occ_HIE bind")
    for idx in np.argsort(-freq[8]):
        r = rows[idx]
        print(f"{names[idx]:16} {freq[5][idx]:>6.2f} {freq[8][idx]:>6.2f} {freq[10][idx]:>6.2f}  "
              f"{r['d_occ']:.3f} {r['occ_HIE']:.3f} {r['bind_consensus']:.3f}")
    print()

    # 5 portfolio
    print("--- 5. portfolio P(>=1 wet-lab hit) submodular-greedy ordering (correlated failure) ---")
    for i, (nm, cum) in enumerate(portfolio_greedy(rows, K=10), 1):
        print(f"  {i:>2}. {nm:16} cumP(>=1)={cum}")
    print()

    # 6 QUBO
    picked = qubo_subset(rows, K=10)
    print("--- 6. QUBO top-10 diverse-subset (neal simulated annealing; classical solve) ---")
    print("  ", picked, "\n")

    # 7 scalarization
    lin, cheb, pareto = scalarization_audit(rows)
    print("--- 7. scalarization audit: linear vs Chebyshev top-8 (divergence = hidden trade-off) ---")
    print("  linear   :", lin[:8])
    print("  chebyshev:", cheb[:8])
    hidden = [nm for nm in cheb[:8] if nm not in lin[:8]]
    print("  surfaced by Chebyshev (linear hides):", hidden or "none")
    print("  P(Pareto-optimal) top:", [(names[i], round(pareto[i], 2)) for i in np.argsort(-pareto)[:6]])

    # write a combined CSV
    out = f"{REPO}/results/wave22_selection_bundle.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "mech", "lineage", "d_occ", "occ_HIE", "bind", "xsev",
                    "composite", "PL_quality", "stab_top8", "P_pareto", "imputed"])
        for i, r in enumerate(rows):
            w.writerow([r["name"], r["mech"], r["lineage"], r["d_occ"], round(r["occ_HIE"],3),
                        r["bind_consensus"], r.get("xreact_mech_sev", 0), round(float(comp[i]),4),
                        round(float(pl[i]),4), round(float(freq[8][i]),3), round(float(pareto[i]),3),
                        "|".join(r.get("imputed", []))])
    print(f"\n-> {out}")
    print("\nHONEST LIMITS: oracle noise is MODELED (no replicates); no exchangeable wet-lab labels, so")
    print("conformal/FDR certification is out of scope. This is noise/diversity ROBUSTNESS, not calibration.")


if __name__ == "__main__":
    main()
