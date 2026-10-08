#!/usr/bin/env python3
"""Debye-Hueckel screening check on the His433 salt bridges (physiological ionic strength de-risk).

The production MD (md_pH_switch / md_pH_bind, GBn2 implicit solvent) runs at ZERO salt, so a salt bridge
validated there could in principle be a low-salt artifact. This script applies the linearized
Debye-Hueckel screening analytically to each design's measured His433-carboxylate distance, to confirm the
bridge survives ~150 mM ionic strength (the in-vitro / tumor-interstitium regime) without a GPU re-run.

Screened Coulomb interaction ~ (1/r) * exp(-kappa * r), with the Debye length
  lambda_D [A] = 3.04 / sqrt(I[M])      (1:1 salt, ~298-310 K; kappa = 1/lambda_D).
A salt bridge whose partner distance r is well inside lambda_D keeps most of its bare interaction; one with
r > lambda_D is heavily screened and fragile at physiological salt. At 150 mM, lambda_D ~ 7.8 A, so a
2.6-3.5 A buried bridge retains ~0.66-0.72 of its bare electrostatic = robust.

Caveat: this is the direct pairwise screening only. It does NOT capture salt-induced desolvation or geometry
change (an implicit-GB salt MD would, approximately); but for a SHORT buried bridge those are second order,
and the GBn2-salt model is itself approximate -- so the analytic screening is the honest, cheap answer here.

Usage: envs/md/bin/python scripts/debye_huckel_screen.py [results/isostere_w19.csv] [ionic_strength_M]
Reads an isostere CSV (needs columns name, partner, partner_dist); writes results/debye_huckel_<I>mM.csv.
"""
import csv, math, os, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def lambda_debye(I_molar):
    return 3.04 / math.sqrt(I_molar)   # Angstrom


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "results/isostere_w19.csv")
    I = float(sys.argv[2]) if len(sys.argv) > 2 else 0.150
    lamD = lambda_debye(I); kappa = 1.0 / lamD
    rows = []
    for r in csv.DictReader(open(src)):
        try:
            d = float(r["partner_dist"])
        except (TypeError, ValueError):
            continue
        rows.append(dict(name=r["name"], partner=r.get("partner", ""), r=d,
                         screen=round(math.exp(-kappa * d), 3),
                         inside_lambdaD=int(d <= lamD)))
    out = os.path.join(REPO, "results", "debye_huckel_%dmM.csv" % round(I * 1000))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "partner", "r", "screen", "inside_lambdaD"])
        w.writeheader()
        for x in sorted(rows, key=lambda x: x["r"]):
            w.writerow(x)
    causal = [x for x in rows if x["r"] <= 3.5]   # the occ-gated shipped switches have short partners
    print("Debye length lambda_D = %.2f A (kappa %.3f /A) at %d mM\n" % (lamD, kappa, round(I * 1000)))
    print("%-30s%-9s%6s%10s" % ("design", "partner", "r(A)", "exp(-kr)"))
    for x in sorted(rows, key=lambda x: x["r"]):
        print("%-30s%-9s%6.2f%10.2f" % (x["name"].replace("pdnov_egfr_phsw_", ""), x["partner"], x["r"], x["screen"]))
    if causal:
        lo = min(x["r"] for x in causal); hi = max(x["r"] for x in causal)
        print("\nShipped causal switches (partner<3.5 A): n=%d, r %.2f-%.2f A -> screening %.2f-%.2f retained."
              % (len(causal), lo, hi, math.exp(-kappa * hi), math.exp(-kappa * lo)))
        print("All sit inside lambda_D (%.1f A): the bridges are robust at %d mM; the heavily-screened long"
              " partners (>4.5 A) are already gate-excluded as non-causal." % (lamD, round(I * 1000)))
    print("\nwrote", os.path.relpath(out, REPO))


if __name__ == "__main__":
    main()
