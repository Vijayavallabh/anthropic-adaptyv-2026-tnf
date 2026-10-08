#!/usr/bin/env python3
"""WAVE 28: fair two-state QM for the carboxylate-DYAD pH switch on a FOLDED pose. Protonate the binder-Glu
as GLH (neutral, tumor/ON) vs GLU (deprotonated -1, blood/OFF) via pdbfixer; the target Asp460 stays ASP
(-1). Cap each sidechain carboxyl to a small acid/acetate and compute GFN2-xTB E_int(binder-Glu, Asp460) in
each state.

  dE_switch = E_int(ON: GluH + Asp-) - E_int(OFF: Glu- + Asp-)   (want < 0: protonation turns binding ON)
  E_int(OFF)                                                      (want > 0 = repulsive = genuinely off)

This is the make-or-break selectivity test for the repulsion-gated dyad, fair to the chemistry (xtb captures
the charged-H-bond ON and the like-charge-repulsion OFF). Reuses qm_switch_fair's extract/xtb helpers."""
import argparse, json, os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qm_switch_fair as Q

H2KCAL = Q.H2KCAL


def protonate_glu(pdb, glu_chain_idx, glu_resid, protonated):
    """add H with the binder-Glu forced to GLH (protonated) or GLU (deprotonated); Asp stays default (-1)."""
    from pdbfixer import PDBFixer
    from openmm.app import Modeller, ForceField
    from openmm import unit
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    residues = list(modeller.topology.residues())
    variants = [None] * len(residues)
    for i, r in enumerate(residues):
        if r.chain.index == glu_chain_idx and r.name in ("GLU", "GLH") and str(r.id) == str(glu_resid):
            variants[i] = "GLH" if protonated else "GLU"
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    pos = np.array(modeller.positions.value_in_unit(unit.angstrom))
    return modeller.topology, pos


def e_int(pdb, glu_ci, glu_resid, asp_ci, asp_resid, protonated, wd, solvent):
    topo, pos = protonate_glu(pdb, glu_ci, glu_resid, protonated)
    glu = Q.extract_cluster(topo, pos, glu_ci, glu_resid, "GLU")
    asp = Q.extract_cluster(topo, pos, asp_ci, asp_resid, "ASP")
    qg = 0 if protonated else -1
    qa = -1
    Q.write_xyz(os.path.join(wd, "cplx.xyz"), glu + asp)
    Q.write_xyz(os.path.join(wd, "g.xyz"), glu)
    Q.write_xyz(os.path.join(wd, "a.xyz"), asp)
    Ec = Q.xtb_energy("cplx.xyz", qg + qa, wd); Eg = Q.xtb_energy("g.xyz", qg, wd); Ea = Q.xtb_energy("a.xyz", qa, wd)
    return (Ec - Eg - Ea) * H2KCAL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--glu", type=int, required=True, help="binder Glu resid")
    ap.add_argument("--asp", type=int, required=True, help="target Asp460 resid (af2pred numbering)")
    ap.add_argument("--binder-chain", default="A"); ap.add_argument("--target-chain", default="B")
    ap.add_argument("--opt", action="store_true"); ap.add_argument("--solvent", default="chcl3")
    ap.add_argument("--name", default=None)
    a = ap.parse_args()
    Q.OPT = a.opt; Q.SOLVENT = a.solvent
    from pdbfixer import PDBFixer
    chain_ids = [c.id for c in PDBFixer(filename=a.pdb).topology.chains()]
    gci, aci = chain_ids.index(a.binder_chain), chain_ids.index(a.target_chain)
    with tempfile.TemporaryDirectory(dir=f"{Q.REPO}/runs") as wd:
        e_on = e_int(a.pdb, gci, a.glu, aci, a.asp, True, wd, a.solvent)    # GluH (tumor, ON)
    with tempfile.TemporaryDirectory(dir=f"{Q.REPO}/runs") as wd:
        e_off = e_int(a.pdb, gci, a.glu, aci, a.asp, False, wd, a.solvent)  # Glu- (blood, OFF)
    res = dict(name=a.name or os.path.basename(a.pdb), glu=a.glu, asp=a.asp, solvent=a.solvent, opt=a.opt,
               Eint_ON=round(e_on, 2), Eint_OFF=round(e_off, 2), dE_switch=round(e_on - e_off, 2),
               off_is_repulsive=bool(e_off > 0))
    print(json.dumps(res))


if __name__ == "__main__":
    main()
