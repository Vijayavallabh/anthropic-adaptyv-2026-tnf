#!/usr/bin/env python3
"""WAVE 25 probe: per-binder-His Wyman linkage -- does a design achieve Hill n>1 (multiple binder His each
protonation-gated for binding)?

For each binder His we compute its pKa in the COMPLEX vs the APO binder (PROPKA3), then its Wyman contribution
  contrib = RT*ln10 * integral_{6.5}^{7.4} [ f(pH; pKa_complex) - f(pH; pKa_apo) ] dpH,   f=1/(1+10^(pH-pKa))
A His that shifts its pKa UP on binding (pKa_complex > pKa_apo) holds more protons when bound -> releases them
on unbinding -> contributes POSITIVELY to tumor-ON selectivity. The COLLECTIVE binder-His linkage = sum over
binder His; n_contributing = # His with contrib > +0.05 kcal/mol. Hill n>1 needs >=2 contributing His.
This is the fair continuum/empirical-pKa route (sums pKa SHIFTS, not a single fixed-charge salt-bridge energy).
"""
import argparse, json, math, os, re, subprocess, tempfile, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
PROPKA = f"{REPO}/envs/pyrosetta/bin/propka3"
RT_LN10 = 1.419   # kcal/mol at 310 K


class _Sel(Select):
    def __init__(self, chains): self.chains = chains
    def accept_chain(self, ch): return ch.id in self.chains


def propka_his_pkas(pdb, workdir, chain):
    """run propka3; return {resid: pKa} for all HIS in `chain`."""
    subprocess.run([PROPKA, os.path.basename(pdb)], capture_output=True, text=True, cwd=workdir)
    pka_file = os.path.join(workdir, os.path.basename(pdb).replace(".pdb", ".pka"))
    if not os.path.exists(pka_file):
        raise RuntimeError("propka3 produced no .pka")
    out = {}
    for ln in open(pka_file):
        m = re.match(r"\s*HIS\s+(\d+)\s+(\w)\s+([-\d.]+)", ln)
        if m and m.group(2) == chain:
            out[int(m.group(1))] = float(m.group(3))
    return out


def f(pH, pKa):
    return 1.0 / (1.0 + 10 ** (pH - pKa))


def wyman_contrib(pKa_cplx, pKa_apo, lo=6.5, hi=7.4, n=19):
    grid = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    vals = [f(pH, pKa_cplx) - f(pH, pKa_apo) for pH in grid]
    integ = sum((vals[i] + vals[i + 1]) / 2 * (grid[i + 1] - grid[i]) for i in range(n - 1))
    return RT_LN10 * integ


def run(pdb, binder_chain, target_chain, name):
    model = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
        io = PDBIO(); io.set_structure(model)
        cplx = os.path.join(wd, "cplx.pdb"); io.save(cplx, _Sel([binder_chain, target_chain]))
        apo = os.path.join(wd, "apo.pdb"); io.save(apo, _Sel([binder_chain]))
        pk_c = propka_his_pkas(cplx, wd, binder_chain)
        pk_a = propka_his_pkas(apo, wd, binder_chain)
    per = []
    for resid, pkc in sorted(pk_c.items()):
        pka = pk_a.get(resid)
        if pka is None:
            continue
        c = wyman_contrib(pkc, pka)
        per.append(dict(his=resid, pKa_complex=round(pkc, 2), pKa_apo=round(pka, 2),
                        dpKa=round(pkc - pka, 2), contrib=round(c, 3)))
    collective = round(sum(p["contrib"] for p in per), 3)
    n_contrib = sum(1 for p in per if p["contrib"] > 0.05)
    return dict(name=name, n_binder_his=len(per), n_contributing=n_contrib,
                collective_linkage=collective, per_his=per)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--binder-chain", default="A"); ap.add_argument("--target-chain", default="B")
    ap.add_argument("--name", default=None)
    a = ap.parse_args()
    res = run(a.pdb, a.binder_chain, a.target_chain, a.name or os.path.basename(a.pdb))
    print(json.dumps(res))


if __name__ == "__main__":
    main()
