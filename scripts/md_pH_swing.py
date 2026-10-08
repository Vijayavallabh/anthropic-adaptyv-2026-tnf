#!/usr/bin/env python3
"""WAVE 29 (ALLOSTERIC PivotLid swing): 2-state MD of a His-gated rigid-body LID swing.

A PivotLid binder presents its EGFR epitope on a small LID lobe, hinged to a BASE lobe and held DOCKED by a
SURFACE His...Asp staple. Tumor pH (HIP, His+) -> staple holds the lid DOCKED -> binds; blood pH (HIE, His0)
-> staple weakens -> the whole lid PEELS off -> paratope undocks -> OFF. Unlike md_pH_switch (which pins ALL
CAs and reads one His-carboxyl distance), here we:
  * RESTRAIN the BASE-lobe + TARGET CAs (k_base), FREE the LID CAs (weak k_lid only, to stop spurious GBn2
    fraying), so the rigid-body swing is physically possible;
  * read COLLECTIVE lid observables each frame (base-superposed): lid displacement vs the docked ref, the
    lid-epitope<->target contact occupancy, and the LID-INTERNAL rmsd (does the freed lid stay FOLDED?).

Make-or-break (the validation handle itself): (1) lid_internal_rmsd stays LOW in BOTH HIP and HIE (freed lid
stays folded, not disordered -> the handle is usable); (2) HIE peels the lid MORE than HIP beyond replica
noise: d_lid_disp = disp(HIE)-disp(HIP) > 0 AND d_epitope_occ = occ(HIP)-occ(HIE) > 0. ECC (--charge-scale
0.8) must preserve the SIGN. Clones md_pH_switch.py build/ECC machinery."""
import argparse, json, os, sys
import numpy as np
from openmm import app, unit, Platform, LangevinMiddleIntegrator, CustomExternalForce
from pdbfixer import PDBFixer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import md_pH_switch as M   # reuse detect_chains, scale_ionizable_charges, ECC machinery


def log(*a): print(*a, file=sys.stderr, flush=True)


def parse_resids(spec):
    """'45-70,80,85-88' -> set of ints."""
    out = set()
    for tok in str(spec).split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok:
            a, b = tok.split("-"); out |= set(range(int(a), int(b) + 1))
        else:
            out.add(int(tok))
    return out


def kabsch(P, Q):
    """rotation+translation aligning P onto Q (both Nx3); returns P aligned."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (P - Pc) @ R.T + Qc


def build(pdb, his_variant, lid_resids, switch_his_pos, k_base, k_lid, charge_scale):
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    target, binder = M.detect_chains(modeller.topology)
    # the staple His is on the binder (lid) at switch_his_pos
    def find_staple_his(top):
        tgt, bnd = M.detect_chains(top)
        his = [r for r in bnd.residues() if r.name.startswith("HI") and str(r.id) == str(switch_his_pos)]
        if not his:
            raise RuntimeError("staple His %s not found on binder chain" % switch_his_pos)
        return tgt, bnd, his[0]
    target, binder, his = find_staple_his(modeller.topology)
    residues = list(modeller.topology.residues())
    variants = [None] * len(residues); variants[residues.index(his)] = his_variant
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    target, binder, his = find_staple_his(modeller.topology)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds,
                             hydrogenMass=4.0 * unit.amu)
    if charge_scale != 1.0:
        M.scale_ionizable_charges(system, modeller.topology, charge_scale)
    # selective CA restraints: base+target at k_base, lid at k_lid
    binder_resids = [r.id for r in binder.residues()]
    lid_set = set(str(x) for x in lid_resids)
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    pos = modeller.positions
    kb = (k_base * unit.kilocalories_per_mole / unit.angstrom**2).value_in_unit(unit.kilojoules_per_mole/unit.nanometer**2)
    kl = (k_lid * unit.kilocalories_per_mole / unit.angstrom**2).value_in_unit(unit.kilojoules_per_mole/unit.nanometer**2)
    lid_ca, base_ca = [], []
    for atom in modeller.topology.atoms():
        if atom.name != "CA":
            continue
        x, y, z = pos[atom.index].value_in_unit(unit.nanometer)
        res = atom.residue
        on_lid = (res.chain.id == binder.id and str(res.id) in lid_set)
        force.addParticle(atom.index, [kl if on_lid else kb, x, y, z])
        (lid_ca if on_lid else base_ca).append(atom.index)
    system.addForce(force)
    return ff, modeller, system, target, binder, his, lid_ca, base_ca


def run_state(pdb, variant, a, device, replica):
    ff, modeller, system, target, binder, his, lid_ca, base_ca = build(
        pdb, variant, parse_resids(a.lid_resids), a.switch_his_pos, a.restraint_k, a.lid_restraint_k,
        getattr(a, "charge_scale", 1.0))
    if not lid_ca:
        raise RuntimeError("no lid CAs matched --lid-resids on the binder chain")
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(1000 * replica + (1 if variant == "HIP" else 2))
    sim = app.Simulation(modeller.topology, system, integ, Platform.getPlatformByName("OpenCL"),
                         {"OpenCLDeviceIndex": str(device), "Precision": "mixed"})
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 1234 + replica)
    sim.step(int(a.equil_ps / 0.004))

    # references (post-min docked pose)
    p0 = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))
    lid_ca = np.array(lid_ca); base_ca = np.array(base_ca)
    ref_lid = p0[lid_ca].copy(); ref_base = p0[base_ca].copy()
    # epitope lid atoms = lid residues with any atom within 4.5 A of target (the docked contacts)
    tgt_atoms = np.array([at.index for at in target.atoms()])
    lid_res_atoms = {}
    for r in binder.residues():
        if str(r.id) in set(str(x) for x in parse_resids(a.lid_resids)):
            lid_res_atoms[r.id] = [at.index for at in r.atoms()]
    tgt_xyz0 = p0[tgt_atoms]
    epi_atoms = []
    for rid, idxs in lid_res_atoms.items():
        dmin = np.min(np.linalg.norm(p0[idxs][:, None, :] - tgt_xyz0[None, :, :], axis=-1))
        if dmin < 4.5:
            epi_atoms += idxs
    epi_atoms = np.array(epi_atoms) if epi_atoms else None

    nsteps = int(a.prod_ps / 0.004); stride = int(a.interval_ps / 0.004)
    disp, intern, occ = [], [], []
    for _ in range(nsteps // stride):
        sim.step(stride)
        p = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))
        # base-superpose the whole frame onto the reference base, then lid displacement vs docked ref
        p_al = kabsch(p[base_ca], ref_base)  # aligned base (sanity); use transform on lid
        # compute rotation/translation from base and apply to lid
        Pc, Qc = p[base_ca].mean(0), ref_base.mean(0)
        H = (p[base_ca] - Pc).T @ (ref_base - Qc)
        U, _, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        lid_now = (p[lid_ca] - Pc) @ R.T + Qc
        disp.append(float(np.sqrt(np.mean(np.sum((lid_now - ref_lid) ** 2, axis=1)))))
        # lid-internal rmsd: superpose lid on its own ref (folding integrity)
        lid_fit = kabsch(p[lid_ca], ref_lid)
        intern.append(float(np.sqrt(np.mean(np.sum((lid_fit - ref_lid) ** 2, axis=1)))))
        if epi_atoms is not None:
            dmin = np.min(np.linalg.norm(p[epi_atoms][:, None, :] - p[tgt_atoms][None, :, :], axis=-1))
            occ.append(1.0 if dmin < 4.0 else 0.0)
    return {
        "variant": variant, "replica": replica, "n_frames": len(disp),
        "lid_displacement": round(float(np.mean(disp)), 2) if disp else None,
        "lid_internal_rmsd": round(float(np.mean(intern)), 2) if intern else None,
        "epitope_occ": round(float(np.mean(occ)), 3) if occ else None,
        "n_epitope_atoms": int(len(epi_atoms)) if epi_atoms is not None else 0,
        "n_lid_ca": int(len(lid_ca)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--lid-resids", required=True, help="binder-chain resids forming the LID, e.g. '55-80'")
    ap.add_argument("--switch-his-pos", type=int, required=True, help="the staple His residue id on the lid")
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--equil-ps", type=float, default=200.0); ap.add_argument("--prod-ps", type=float, default=2000.0)
    ap.add_argument("--interval-ps", type=float, default=10.0); ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--restraint-k", type=float, default=2.0, help="base+target CA restraint k (kcal/mol/A^2)")
    ap.add_argument("--lid-restraint-k", type=float, default=0.1, help="weak lid CA restraint (anti-fray)")
    ap.add_argument("--charge-scale", type=float, default=1.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    name = os.path.splitext(os.path.basename(a.pdb))[0]
    log("=== MD pH-swing:", name, "device", a.device, "lid", a.lid_resids, "His", a.switch_his_pos)
    res = {"name": name, "pdb": os.path.abspath(a.pdb), "lid_resids": a.lid_resids,
           "switch_his_pos": a.switch_his_pos, "charge_scale": a.charge_scale, "states": []}
    for rep in range(a.replicas):
        for v in ("HIP", "HIE"):
            try:
                r = run_state(a.pdb, v, a, a.device, rep); res["states"].append(r)
                log("  rep%d %s disp=%s intern=%s occ=%s" % (rep, v, r["lid_displacement"], r["lid_internal_rmsd"], r["epitope_occ"]))
            except Exception as e:
                log("  rep%d %s FAILED: %r" % (rep, v, e)); res["states"].append({"variant": v, "replica": rep, "error": repr(e)})

    def agg(metric, v):
        vals = [s[metric] for s in res["states"] if s.get("variant") == v and s.get(metric) is not None]
        return float(np.mean(vals)) if vals else None
    disp_hip, disp_hie = agg("lid_displacement", "HIP"), agg("lid_displacement", "HIE")
    occ_hip, occ_hie = agg("epitope_occ", "HIP"), agg("epitope_occ", "HIE")
    int_hip, int_hie = agg("lid_internal_rmsd", "HIP"), agg("lid_internal_rmsd", "HIE")
    res["disp_HIP"], res["disp_HIE"] = disp_hip, disp_hie
    res["occ_HIP"], res["occ_HIE"] = occ_hip, occ_hie
    res["internal_HIP"], res["internal_HIE"] = int_hip, int_hie
    res["d_lid_disp"] = round(disp_hie - disp_hip, 2) if (disp_hip is not None and disp_hie is not None) else None
    res["d_epitope_occ"] = round(occ_hip - occ_hie, 3) if (occ_hip is not None and occ_hie is not None) else None
    # make-or-break flags
    res["lid_stays_folded"] = (int_hip is not None and int_hie is not None and int_hip < 3.0 and int_hie < 3.0)
    res["peels_at_HIE"] = (res["d_lid_disp"] is not None and res["d_lid_disp"] > 0 and
                           res["d_epitope_occ"] is not None and res["d_epitope_occ"] > 0)
    res["swing_GO"] = bool(res["lid_stays_folded"] and res["peels_at_HIE"])
    outpath = os.path.join(a.out, name + ".json")
    json.dump(res, open(outpath, "w"), indent=2)
    log("WROTE", outpath, "| d_lid_disp=%s d_epitope_occ=%s lid_folded=%s peels=%s GO=%s" % (
        res["d_lid_disp"], res["d_epitope_occ"], res["lid_stays_folded"], res["peels_at_HIE"], res["swing_GO"]))
    print(outpath)


if __name__ == "__main__":
    main()
