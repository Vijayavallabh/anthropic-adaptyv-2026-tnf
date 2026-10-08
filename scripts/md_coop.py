#!/usr/bin/env python3
"""WAVE 23 downstream: explicit-protonation 2-state MD for the COOPERATIVE two-His switch.

Sets BOTH binder His to the same protonation state -- HIP (tumor pH 6.5) vs HIE (blood 7.4) -- and
measures BOTH salt bridges (His_a .. Asp460 ordinal 127, His_b .. Asp344 ordinal 35) each frame. A
cooperative tumor-ON switch requires BOTH bridges to form together in the all-protonated state and break
together in the all-neutral state. Reports per-bridge occ_HIP / occ_HIE / d_occ AND the key COOPERATIVE
metric occ_BOTH (fraction of frames with BOTH bridges formed) in each state -> d_occ_both = occ_BOTH(HIP)
- occ_BOTH(HIE), the two-contact switch signal that a single-His design cannot provide.

Reuses md_pH_switch's OpenMM/forcefield machinery (detect_chains, CA restraints, integrator, HMR).
"""
import argparse, glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import md_pH_switch as M
app, unit = M.app, M.unit

ASP_A_ORD = 127  # Asp460
ASP_B_ORD = 35   # Asp344


def build_two_his(pdb, variant, his_a, his_b):
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    target, binder = M.detect_chains(modeller.topology)
    # set both binder His (by residue id) to the variant
    residues = list(modeller.topology.residues())
    variants = [None] * len(residues)
    nset = 0
    for i, r in enumerate(residues):
        if r.chain.index == binder.index and r.name.startswith("HI") and str(r.id) in (str(his_a), str(his_b)):
            variants[i] = variant; nset += 1
    if nset != 2:
        M.log("  WARN: set %d His variants (expected 2) for %s/%s" % (nset, his_a, his_b))
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds,
                             hydrogenMass=4.0 * unit.amu)
    target, binder = M.detect_chains(modeller.topology)
    return ff, modeller, system, target, binder


def his_bridge_atoms(topology, binder, target, his_id, asp_ord):
    """return (his-N indices, asp-carboxyl-O indices) for the His at binder residue his_id and the
    target Asp at ordinal asp_ord."""
    his_N, asp_O = [], []
    for r in binder.residues():
        if r.name.startswith("HI") and str(r.id) == str(his_id):
            his_N = [a.index for a in r.atoms() if a.name in ("ND1", "NE2")]
    tgt_res = list(target.residues())
    if asp_ord - 1 < len(tgt_res):
        ar = tgt_res[asp_ord - 1]
        if ar.name in ("ASP", "GLU"):
            asp_O = [a.index for a in ar.atoms() if a.name in ("OD1", "OD2", "OE1", "OE2")]
    return his_N, asp_O


def run_state(pdb, variant, his_a, his_b, device, equil_ps, prod_ps, interval_ps, replica, restraint_k):
    ff, modeller, system, target, binder = build_two_his(pdb, variant, his_a, his_b)
    nA = his_bridge_atoms(modeller.topology, binder, target, his_a, ASP_A_ORD)
    nB = his_bridge_atoms(modeller.topology, binder, target, his_b, ASP_B_ORD)
    M.add_ca_restraints(system, modeller, restraint_k)
    integ = M.LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(1000 * replica + (1 if variant == "HIP" else 2))
    plat = M.Platform.getPlatformByName("OpenCL")
    sim = app.Simulation(modeller.topology, system, integ, plat,
                         {"OpenCLDeviceIndex": str(device), "Precision": "mixed"})
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 1234 + replica)
    sim.step(int(equil_ps / 0.004))
    nsteps, stride = int(prod_ps / 0.004), int(interval_ps / 0.004)
    da, db = [], []
    for _ in range(nsteps // stride):
        sim.step(stride)
        p = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))
        for (hN, aO), store in ((nA, da), (nB, db)):
            if hN and aO:
                store.append(float(np.min(np.linalg.norm(p[aO][:, None, :] - p[hN][None, :, :], axis=-1))))
    da, db = np.array(da), np.array(db)
    both = (da < 4.0) & (db < 4.0) if len(da) and len(db) else np.array([])
    return dict(variant=variant, replica=replica, n_frames=len(da),
                occ_a=round(float(np.mean(da < 4.0)), 3) if len(da) else None,
                occ_b=round(float(np.mean(db < 4.0)), 3) if len(db) else None,
                occ_both=round(float(np.mean(both)), 3) if len(both) else None,
                mindist_a=round(float(np.mean(da)), 2) if len(da) else None,
                mindist_b=round(float(np.mean(db)), 2) if len(db) else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--his-a", type=int, required=True); ap.add_argument("--his-b", type=int, required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--equil-ps", type=float, default=200.0); ap.add_argument("--prod-ps", type=float, default=2000.0)
    ap.add_argument("--interval-ps", type=float, default=10.0); ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--restraint-k", type=float, default=2.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    name = os.path.splitext(os.path.basename(a.pdb))[0]
    M.log("=== MD coop 2-His:", name, "His", a.his_a, a.his_b, "device", a.device)
    states = []
    for rep in range(a.replicas):
        for variant in ("HIP", "HIE"):
            states.append(run_state(a.pdb, variant, a.his_a, a.his_b, a.device,
                                    a.equil_ps, a.prod_ps, a.interval_ps, rep, a.restraint_k))
    def avg(key, var):
        xs = [s[key] for s in states if s["variant"] == var and s[key] is not None]
        return round(float(np.mean(xs)), 3) if xs else None
    res = dict(name=name, pdb=os.path.abspath(a.pdb), his_a=a.his_a, his_b=a.his_b, states=states,
               occ_a_HIP=avg("occ_a", "HIP"), occ_a_HIE=avg("occ_a", "HIE"),
               occ_b_HIP=avg("occ_b", "HIP"), occ_b_HIE=avg("occ_b", "HIE"),
               occ_both_HIP=avg("occ_both", "HIP"), occ_both_HIE=avg("occ_both", "HIE"))
    res["d_occ_a"] = round((res["occ_a_HIP"] or 0) - (res["occ_a_HIE"] or 0), 3)
    res["d_occ_b"] = round((res["occ_b_HIP"] or 0) - (res["occ_b_HIE"] or 0), 3)
    res["d_occ_both"] = round((res["occ_both_HIP"] or 0) - (res["occ_both_HIE"] or 0), 3)
    # cooperative clean switch: both bridges form together at tumor, break together at blood
    res["coop_switch"] = (res["occ_both_HIP"] or 0) > 0.3 and (res["occ_both_HIE"] or 0) < 0.15
    json.dump(res, open(os.path.join(a.out, name + ".json"), "w"), indent=1)
    M.log("  occ_both HIP=%s HIE=%s  d_occ_both=%s  bridges d_occ_a=%s d_occ_b=%s  coop_switch=%s" % (
        res["occ_both_HIP"], res["occ_both_HIE"], res["d_occ_both"], res["d_occ_a"], res["d_occ_b"], res["coop_switch"]))


if __name__ == "__main__":
    main()
