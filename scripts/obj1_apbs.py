#!/usr/bin/env python3
"""WAVE 22h: continuum Poisson-Boltzmann pH-linkage leg (deep-research physics flagship, orthogonal).

An INDEPENDENT selectivity sign-vote from first-principles continuum electrostatics (APBS nonlinear PB),
orthogonal to PROPKA (empirical), pKAI/JustHISpKa (ML), and fixed-charge MD d_occ. It alone computes
explicit Born desolvation + full charge-charge coupling of the His+ / Asp- pair on binding.

Binding-coupled electrostatic linkage (the continuum analog of the Wyman linkage we already integrate):
    dG_prot(complex) = G_elec(complex, His=HIP) - G_elec(complex, His=HIE)
    dG_prot(binder)  = G_elec(binder,  His=HIP) - G_elec(binder,  His=HIE)
    DDG_link = dG_prot(complex) - dG_prot(binder)     # target cancels (its state is unchanged)
DDG_link < 0  <=>  protonating the switch His is electrostatically CHEAPER in the complex than free,
i.e. the target (conserved Asp460) stabilizes His+, raising its pKa on binding = tumor-ON switch.
We trust the SIGN (continuum PB pKa error is ~1.1 pH units for buried/salt-bridged His); never the magnitude.

Robustness: 4 APBS runs/design (cplx-HIP, cplx-HIE, bind-HIP, bind-HIE); His state forced by pre-naming
the residue HIP/HIE and pdb2pqr30 --noopt --ff=AMBER (keeps the explicit state); fail-loud on any APBS/pdb2pqr
nonzero exit or unparseable energy. Run in the repo base python; uses envs/apbs binaries.
"""
import argparse, os, re, subprocess, sys, tempfile, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB2PQR = f"{REPO}/envs/apbs/bin/pdb2pqr30"
APBS = f"{REPO}/envs/apbs/bin/apbs"


class _Sel(Select):
    def __init__(self, chains): self.chains = chains
    def accept_chain(self, ch): return ch.id in self.chains


def write_subset(model, chains, his_chain, his_resid, his_state, out_pdb):
    """Write the given chains with the switch His renamed to his_state (HIP/HIE)."""
    # rename the switch His in-place (restore after)
    target_res = None
    for ch in model:
        if ch.id == his_chain:
            for r in ch:
                if r.id[1] == his_resid and r.resname in ("HIS", "HIP", "HIE", "HID"):
                    target_res = r; break
    if target_res is None:
        raise RuntimeError(f"switch His {his_chain}{his_resid} not found")
    orig = target_res.resname
    target_res.resname = his_state
    io = PDBIO(); io.set_structure(model); io.save(out_pdb, _Sel(chains))
    target_res.resname = orig


def run_apbs(pdb_in, workdir, tag):
    """pdb2pqr30 (AMBER, noopt, keep explicit His state) -> APBS nonlinear PB -> total elec energy (kJ/mol)."""
    pqr = os.path.join(workdir, f"{tag}.pqr")
    apbs_in = os.path.join(workdir, f"{tag}.in")
    r = subprocess.run([PDB2PQR, "--ff=AMBER", "--noopt", "--keep-chain",
                        f"--apbs-input={apbs_in}", pdb_in, pqr],
                       capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0 or not os.path.exists(pqr):
        raise RuntimeError(f"pdb2pqr failed ({tag}) rc={r.returncode}: {r.stderr[-300:]}")
    # force nonlinear PB (npbe) in the generated input
    txt = open(apbs_in).read().replace("lpbe", "npbe")
    open(apbs_in, "w").write(txt)
    r = subprocess.run([APBS, os.path.basename(apbs_in)], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        raise RuntimeError(f"apbs failed ({tag}) rc={r.returncode}: {r.stderr[-300:]}")
    # parse "Global net ELEC energy" / "totEnergy" (kJ/mol)
    m = re.findall(r"Global net (?:ELEC )?energy\s*=\s*([-\d.Ee+]+)", r.stdout)
    if not m:
        m = re.findall(r"totEnergy\s+([-\d.Ee+]+)", r.stdout)
    if not m:
        raise RuntimeError(f"no APBS energy parsed ({tag}); tail: {r.stdout[-300:]}")
    return float(m[-1])


def linkage(pdb, binder_chain, target_chain, his_resid, workdir):
    model = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    E = {}
    for subset_tag, chains in (("cplx", [binder_chain, target_chain]), ("bind", [binder_chain])):
        for state in ("HIP", "HIE"):
            tag = f"{subset_tag}_{state}"
            sub_pdb = os.path.join(workdir, f"{tag}.pdb")
            write_subset(model, chains, binder_chain, his_resid, state, sub_pdb)
            E[tag] = run_apbs(sub_pdb, workdir, tag)
    dG_prot_cplx = E["cplx_HIP"] - E["cplx_HIE"]
    dG_prot_bind = E["bind_HIP"] - E["bind_HIE"]
    ddg_link = dG_prot_cplx - dG_prot_bind
    return dict(E=E, dG_prot_cplx=round(dG_prot_cplx, 2), dG_prot_bind=round(dG_prot_bind, 2),
                ddg_link_kJ=round(ddg_link, 2), tumor_on=(ddg_link < 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--target-chain", default="B")
    ap.add_argument("--his-resid", type=int, required=True, help="switch His residue id in the binder chain")
    ap.add_argument("--name", default=None)
    a = ap.parse_args()
    name = a.name or os.path.basename(a.pdb)
    with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
        res = linkage(a.pdb, a.binder_chain, a.target_chain, a.his_resid, wd)
    import json
    res["name"] = name
    print(json.dumps(res, indent=2))
    print(f"\n[APBS continuum linkage] {name}: DDG_link = {res['ddg_link_kJ']} kJ/mol -> "
          f"{'TUMOR-ON (sign agrees with switch)' if res['tumor_on'] else 'NOT tumor-on (sign DISAGREES)'}")


if __name__ == "__main__":
    main()
