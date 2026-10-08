#!/usr/bin/env python3
"""WAVE 22i: His pKa-window + constitutive-bridge-shift check (deep-research chemistry Tier-1).

For an inverted-switch design, the binder His must (a) titrate IN the tumor-pH window and (b) NOT have its
pKa driven far up by binding (a large complex-vs-free up-shift is the constitutive-bridge signature that
fools static scores yet binds at both pH). We compute the binder switch-His pKa with PROPKA3 in the COMPLEX
and in the APO binder, and report:
  pKa(complex), pKa(free), shift = pKa(complex) - pKa(free)
  in_window  := pKa(complex) in [6.3, 6.8]  (sharp transition parked between tumor 6.5 and blood 7.4)
  shift_ok   := shift <= +0.5                (bridge does not over-elevate the pKa -> keeps an OFF state)
PROPKA is empirical and single-structure (the harshest proxy); trust the shift SIGN/magnitude qualitatively.
"""
import argparse, os, re, subprocess, sys, tempfile, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROPKA = f"{REPO}/envs/pyrosetta/bin/propka3"


class _Sel(Select):
    def __init__(self, chains): self.chains = chains
    def accept_chain(self, ch): return ch.id in self.chains


def propka_his_pka(pdb, workdir, his_chain, his_resid):
    """run propka3, return the pKa of the given His (chain+resid)."""
    r = subprocess.run([PROPKA, os.path.basename(pdb)], capture_output=True, text=True, cwd=workdir)
    pka_file = os.path.join(workdir, os.path.basename(pdb).replace(".pdb", ".pka"))
    if not os.path.exists(pka_file):
        raise RuntimeError(f"propka3 produced no .pka (rc={r.returncode}): {r.stderr[-200:]}")
    # SUMMARY section: lines like "HIS  41 A   <pKa>  ..."
    best = None
    for ln in open(pka_file):
        m = re.match(r"\s*HIS\s+(\d+)\s+(\w)\s+([-\d.]+)", ln)
        if m and int(m.group(1)) == his_resid and m.group(2) == his_chain:
            best = float(m.group(3))
    if best is None:
        raise RuntimeError(f"HIS {his_chain}{his_resid} pKa not found in {pka_file}")
    return best


def run_design(pdb, binder_chain, target_chain, his_resid, name):
    model = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
        cplx = os.path.join(wd, "cplx.pdb"); io = PDBIO(); io.set_structure(model)
        io.save(cplx, _Sel([binder_chain, target_chain]))
        apo = os.path.join(wd, "apo.pdb"); io.save(apo, _Sel([binder_chain]))
        pka_c = propka_his_pka(cplx, wd, binder_chain, his_resid)
        pka_f = propka_his_pka(apo, wd, binder_chain, his_resid)
    shift = pka_c - pka_f
    return dict(name=name, pKa_complex=round(pka_c, 2), pKa_free=round(pka_f, 2), shift=round(shift, 2),
                in_window=(6.3 <= pka_c <= 6.8), shift_ok=(shift <= 0.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--binder-chain", default="A"); ap.add_argument("--target-chain", default="B")
    ap.add_argument("--his-resid", type=int, required=True); ap.add_argument("--name", default=None)
    a = ap.parse_args()
    res = run_design(a.pdb, a.binder_chain, a.target_chain, a.his_resid, a.name or os.path.basename(a.pdb))
    import json; print(json.dumps(res))


if __name__ == "__main__":
    main()
