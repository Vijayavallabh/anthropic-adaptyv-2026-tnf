#!/usr/bin/env python3
"""WAVE 22e: diversity-maximizing reorder of the confirmed PASS20 set via DPP-MAP (greedy).

The Challenge-1 CSV is submitted best-first and a downstream process picks a SUBSET for wet-lab, with
best-of (no mean penalty) standing. So the right objective is to maximize expected DISTINCT successes --
quality-weighted DIVERSITY over the top of the list -- not sum-of-scores. The current order buries the two
inverted-mechanism designs (the most distinct entries, a different switch class) at slots 19-20, so a top-K
pick would test 18 His433-face cousins and exclude both new mechanisms.

Model the 20 as a determinantal point process: L = diag(q) . S . diag(q), with q a gentle quality factor
and S a PSD cosine-similarity kernel over features [k-mer sequence composition (+) mechanism class (+)
lineage]. Greedy DPP-MAP (Chen et al. 2018, log-det marginal gain) yields a best-first order that is
high-quality AND spreads sequence/fold/mechanism/lineage -- directly the "diverse portfolio, not cousins"
axis. ZERO floor risk: the SET of 20 is fixed (all portal-confirmed), this only reorders.
"""
import argparse, csv, json, math, os, re
import numpy as np

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
PASS20 = f"{REPO}/results/adaptyv_challenge1_resubmit_PASS20_inv.csv"
TABLE = f"{REPO}/docs/DESIGN_TABLE_W21.md"


def parse_table():
    """name -> {lineage, d_occ, occ_HIE, bind} for the His433-face designs."""
    out = {}
    for ln in open(TABLE):
        if not re.match(r"\| *\d+ *\| *egfr_phsw_", ln):
            continue
        c = [x.strip() for x in ln.strip().strip("|").split("|")]
        # rank name tier len lineage d_occ occ_HIE novelty oracleTM bind mousesev expr parent
        name, lineage = c[1], c[4]
        def f(x):
            try: return float(x)
            except: return None
        out[name] = dict(lineage=lineage, d_occ=f(c[5]), occ_HIE=f(c[6]), bind=f(c[9]))
    return out


def kmer_vec(seq, k=3):
    from collections import Counter
    cnt = Counter(seq[i:i+k] for i in range(len(seq)-k+1))
    return cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w-seq", type=float, default=1.0)
    ap.add_argument("--w-mech", type=float, default=1.1, help="mechanism-class feature weight")
    ap.add_argument("--w-lin", type=float, default=0.8, help="lineage feature weight")
    ap.add_argument("--qmin", type=float, default=0.6, help="floor of the gentle quality factor")
    ap.add_argument("--out", default=f"{REPO}/results/adaptyv_challenge1_resubmit_PASS20_dpp.csv")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(PASS20)))
    tbl = parse_table()
    # inverted-mechanism metrics (this session) for inv01/inv02
    inv = {
        "egfr_phsw_inv01": dict(lineage="inv_moon03", mech="inverted", d_occ=0.515, occ_HIE=0.00, bind=0.623),
        "egfr_phsw_inv02": dict(lineage="inv_des1_1", mech="inverted", d_occ=0.87,  occ_HIE=0.13, bind=0.786),
    }
    meta = []
    for r in rows:
        nm = r["name"]
        if nm in inv:
            m = dict(inv[nm])
        else:
            t = tbl.get(nm, {})
            m = dict(lineage=t.get("lineage") or nm, mech="his433", d_occ=t.get("d_occ"), bind=t.get("bind"))
        m["name"] = nm; m["seq"] = r["sequence"]
        meta.append(m)
    n = len(meta)

    # ---- features -> PSD cosine kernel S ----
    # k-mer vocabulary
    vocab = sorted({km for m in meta for km in kmer_vec(m["seq"])})
    vidx = {km: i for i, km in enumerate(vocab)}
    mechs = sorted({m["mech"] for m in meta}); midx = {x: i for i, x in enumerate(mechs)}
    lins = sorted({m["lineage"] for m in meta}); lidx = {x: i for i, x in enumerate(lins)}
    D = len(vocab) + len(mechs) + len(lins)
    F = np.zeros((n, D))
    for i, m in enumerate(meta):
        kc = kmer_vec(m["seq"]); tot = sum(kc.values()) or 1
        for km, c in kc.items():
            F[i, vidx[km]] = a.w_seq * c / tot
        F[i, len(vocab) + midx[m["mech"]]] = a.w_mech
        F[i, len(vocab) + len(mechs) + lidx[m["lineage"]]] = a.w_lin
    # cosine normalize -> S = F F^T is PSD with unit diagonal
    F = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    S = F @ F.T

    # ---- gentle quality factor q in [qmin, 1] from validated metrics ----
    def norm(vals):
        v = np.array([x if x is not None else np.nan for x in vals], float)
        lo, hi = np.nanmin(v), np.nanmax(v)
        v = np.where(np.isnan(v), lo, v)
        return (v - lo) / (hi - lo + 1e-12)
    qd = norm([m.get("d_occ") for m in meta])      # switch strength (obj-1)
    qb = norm([m.get("bind") for m in meta])        # co-fold bind
    qual = 0.5 * qd + 0.5 * qb
    q = a.qmin + (1 - a.qmin) * qual

    L = (q[:, None] * q[None, :]) * S
    L += 1e-9 * np.eye(n)  # numerical PSD

    # ---- greedy DPP-MAP ordering (incremental log-det marginal gain) ----
    order, remaining = [], list(range(n))
    while remaining:
        best, bestgain = None, -1e18
        for j in remaining:
            sub = order + [j]
            sign, logdet = np.linalg.slogdet(L[np.ix_(sub, sub)])
            gain = logdet if sign > 0 else -1e18
            if gain > bestgain:
                bestgain, best = gain, j
        order.append(best); remaining.remove(best)

    dpp_names = [meta[i]["name"] for i in order]
    incumbent = [m["name"] for m in meta]

    def mech_cov(names, K):
        top = names[:K]
        return sum(1 for nm in top if nm.startswith("egfr_phsw_inv"))
    print("=== mechanism coverage (# inverted-mechanism designs in top-K) ===")
    print(f"{'K':>3} {'incumbent':>10} {'DPP':>5}")
    for K in (3, 5, 8, 10):
        print(f"{K:>3} {mech_cov(incumbent,K):>10} {mech_cov(dpp_names,K):>5}")
    print("\n=== DPP order ===")
    for rk, i in enumerate(order, 1):
        m = meta[i]
        print(f"{rk:>2}. {m['name']:16} mech={m['mech']:9} lineage={m['lineage']:14} "
              f"d_occ={m.get('d_occ')} bind={m.get('bind')} q={q[i]:.2f}")

    # ---- write reordered CSV (same 20, byte-identical rows, reordered) ----
    byname = {r["name"]: r for r in rows}
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["name", "sequence", "molecule_class"])
        for nm in dpp_names:
            r = byname[nm]; w.writerow([r["name"], r["sequence"], r["molecule_class"]])
    # integrity: same set
    assert set(dpp_names) == set(incumbent) and len(dpp_names) == 20
    print(f"\nwrote {a.out} (same 20, reordered; set identical to PASS20)")


if __name__ == "__main__":
    main()
