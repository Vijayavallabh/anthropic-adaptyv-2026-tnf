#!/usr/bin/env python3
"""WAVE 19 linkage/cooperativity: quantify switch MAGNITUDE and STEEPNESS, not just direction.

The PROPKA/pKAI ddg_pH score is the Wyman-Tanford linkage integral and only tells us the switch points
the right way. Best-of scoring and wet-lab detectability reward a STEEP, LARGE pH response. This module
extends the same linkage math to the predicted tumor/blood affinity fold-change and an effective
cooperativity, and shows quantitatively why a single histidine caps selectivity near 5-6x while multiple
coupled ionizable contacts (the bidentate His433+His370 route) give a sharper, larger switch.

Thermodynamics (Wyman linkage), proton uptake on binding Dn_H(pH) = sum over binding-linked titratable
groups of [theta(pH; pKa_complex) - theta(pH; pKa_free)], theta = protonated fraction 1/(1+10^(pH-pKa)):
  log10[ Kd(blood)/Kd(tumor) ] = integral_{tumor}^{blood}? -- computed as I = integral_{pH_tumor}^{pH_blood} Dn_H dpH
  fold_change (Kd_blood / Kd_tumor, tumor-ON when >1) = 10**I
  ddg_pH (kcal/mol, the Wyman score) = RT*ln10*I   (RT = 0.616 kcal/mol at 310 K)
This reproduces the stored ddg_pH for a single-His design (validated in the self-test), so it is the same
quantity extended to a fold-change and a steepness (max |dDn_H/dpH|, which rises when groups titrate
together). Pure stdlib, read-only, no GPU.
"""
import argparse, csv, os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RT = 0.616        # kcal/mol at 310 K (body/tumor temperature)
LN10 = 2.302585
PH_TUMOR, PH_BLOOD = 6.5, 7.4


def theta_protonated(pH, pka):
    """Fraction of a titratable base (His) in the protonated, positively charged form at pH."""
    return 1.0 / (1.0 + 10 ** (pH - pka))


def delta_nH(pH, groups):
    """Protons taken up on binding at this pH, summed over binding-linked titratable groups.
    Each group is (pka_free, pka_complex): binding that raises a His pKa holds its proton to higher pH."""
    return sum(theta_protonated(pH, c) - theta_protonated(pH, f) for (f, c) in groups)


def _trapz(ys, xs):
    return sum((xs[i + 1] - xs[i]) * (ys[i + 1] + ys[i]) / 2.0 for i in range(len(xs) - 1))


def linkage(groups, pH_lo=PH_TUMOR, pH_hi=PH_BLOOD, n=181):
    """Return the linkage observables for a set of binding-linked titratable groups."""
    xs = [pH_lo + (pH_hi - pH_lo) * i / (n - 1) for i in range(n)]
    ys = [delta_nH(x, groups) for x in xs]
    integral = _trapz(ys, xs)                  # integral of Dn_H over [tumor, blood]
    slopes = [abs((ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])) for i in range(len(xs) - 1)]
    return {
        "integral_dnH": round(integral, 4),
        "ddg_pH_kcal": round(RT * LN10 * integral, 3),       # the Wyman score (tumor-ON if > 0)
        "fold_change": round(10 ** integral, 2),             # Kd(blood)/Kd(tumor)
        "peak_dnH": round(max(ys), 3),                        # how many protons the switch moves (>1 => multi-group)
        "max_steepness": round(max(slopes), 3),               # d Dn_H / d pH; higher => sharper (cooperative) switch
    }


def demo():
    """Single His vs two coupled His (the bidentate His433+His370 route): show the ceiling and the gain."""
    print("Linkage / cooperativity: single-His ceiling vs a multi-ionizable (bidentate) switch")
    print("(pKa_free -> pKa_complex per group; window pH %.1f tumor / %.1f blood)\n" % (PH_TUMOR, PH_BLOOD))
    cases = [
        ("single His433 (6.27->8.45, observed egfr_phsw_22 parent)", [(6.27, 8.45)]),
        ("single His, pKa centered in window (6.3->7.3)", [(6.3, 7.3)]),
        ("single His, strong shift (6.2->9.0)", [(6.2, 9.0)]),
        ("bidentate His433+His370 (6.27->8.45 and 6.3->7.6)", [(6.27, 8.45), (6.3, 7.6)]),
        ("bidentate, both centered + coupled (6.4->7.4 and 6.4->7.4)", [(6.4, 7.4), (6.4, 7.4)]),
    ]
    print(f"{'case':<56}{'dDn_H':>7}{'ddg_pH':>8}{'fold':>7}{'peak':>6}{'steep':>7}")
    for label, groups in cases:
        r = linkage(groups)
        print(f"{label:<56}{r['integral_dnH']:>7.3f}{r['ddg_pH_kcal']:>8.3f}{r['fold_change']:>6.1f}x{r['peak_dnH']:>6.2f}{r['max_steepness']:>7.2f}")
    print("\nTakeaway: one His moves at most ~1 proton, capping the fold-change in the 0.9-pH tumor/blood")
    print("window at ~5-6x; a second coupled ionizable contact that also titrates in the window moves the")
    print("integral (and so the fold-change) toward 2 protons and sharpens the transition -- the quantitative")
    print("rationale for the obj1_switch_network.py bidentate bet (His433 + His370).")


def apply_designs(path="results/pkai_w18.csv", top=15):
    """Rank designs by predicted single-His switch magnitude from stored His433 free/complex pKa."""
    if not os.path.exists(path):
        print("no %s" % path); return
    rows = []
    for r in csv.DictReader(open(path)):
        try:
            f, c = float(r["his433_pka_free"]), float(r["his433_pka_complex"])
        except (TypeError, ValueError, KeyError):
            continue
        if c <= f:
            continue   # not tumor-ON (complex pKa must exceed free for a low-pH-ON His switch)
        res = linkage([(f, c)])
        rows.append((r["name"], f, c, res))
    rows.sort(key=lambda x: -x[3]["fold_change"])
    print(f"\nPredicted single-His switch magnitude (His433 free->complex, from {os.path.basename(path)}):")
    print(f"{'design':<34}{'free':>6}{'cplx':>6}{'ddg_pH':>8}{'fold':>7}")
    for name, f, c, res in rows[:top]:
        print(f"{name:<34}{f:>6.2f}{c:>6.2f}{res['ddg_pH_kcal']:>8.3f}{res['fold_change']:>6.1f}x")


def _selftest():
    # single His433 6.27 -> 8.45 must reproduce the stored pKAI ddg_pH (~0.97) and a ~5x fold-change
    r = linkage([(6.27, 8.45)])
    assert 0.80 <= r["ddg_pH_kcal"] <= 1.10, r          # stored ddg_pH_pkai for this design is 0.972
    assert 3.5 <= r["fold_change"] <= 6.5, r             # single-His ceiling regime
    assert r["peak_dnH"] <= 1.01, r                       # one His moves at most ~1 proton
    # a second coupled group raises the integral, fold-change, and peak proton count
    r2 = linkage([(6.27, 8.45), (6.3, 7.6)])
    assert r2["fold_change"] > r["fold_change"], (r, r2)
    assert r2["peak_dnH"] > r["peak_dnH"], (r, r2)
    assert r2["integral_dnH"] > r["integral_dnH"], (r, r2)
    # direction: no pKa shift on binding => no switch
    r0 = linkage([(6.3, 6.3)])
    assert abs(r0["integral_dnH"]) < 1e-6 and abs(r0["fold_change"] - 1.0) < 1e-6, r0
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser(description="Wyman-linkage switch magnitude and steepness.")
    ap.add_argument("--demo", action="store_true", help="single-His ceiling vs bidentate comparison")
    ap.add_argument("--designs", nargs="?", const="results/pkai_w18.csv",
                    help="rank designs by predicted single-His fold-change from a pKAI CSV")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        _selftest(); return
    if a.demo:
        demo()
    if a.designs:
        apply_designs(a.designs)
    if not (a.demo or a.designs):
        demo()


if __name__ == "__main__":
    main()
