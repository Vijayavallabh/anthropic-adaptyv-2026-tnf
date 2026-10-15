#!/usr/bin/env python
"""Build a solvated OpenMM base system for the pH-switch alchemical ddG campaign.

Switch-His are built as HIP (protonated) so the Hd1 atom + its bond/angle/constraint
terms exist natively (the alchemy driver scales HIE<->HIP via NonbondedForce offsets on
the His atoms). Non-switch His (binder + TNF) are neutral HIE. Disulfides -> CYX.

Legs:
  bound = binder (chain A) + TNF (chain B); free = binder only (--free).
For artifact cancellation the free leg is solvated into the SAME box vectors as its bound
leg (--box-from <bound>_box.txt), so the net-charge PME finite-size term cancels in
ddG = dG_bound - dG_free.

Outputs: <out>.pdb, <out>_system.xml, <out>_state.xml, <out>_box.txt, <out>_meta.json
Run with envs/md/bin/python; the CIF is read only as coordinates (never imported).
"""
import sys, argparse, time, json
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import (PDBxFile, PDBFile, Modeller, ForceField, PME, HBonds,
                        Simulation, StateDataReporter)
from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat, Platform,
                    XmlSerializer, Vec3)
from openmm.unit import nanometer, molar, kelvin, picosecond, femtosecond, bar, amu

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("out")
    ap.add_argument("--switch-his", default="", help="comma binder His resSeqs built as HIP")
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--free", action="store_true")
    ap.add_argument("--pad", type=float, default=1.2)
    ap.add_argument("--box-from", default="", help="read box vectors (nm, 3x3) from this file")
    ap.add_argument("--equil-ps", type=float, default=200.0)
    ap.add_argument("--device", default="0")
    args = ap.parse_args()
    prot = set(int(x) for x in args.switch_his.split(",") if x.strip())

    # --- load CIF, optionally drop TNF, write a clean PDB for PDBFixer ---
    cif = PDBxFile(args.cif)
    m0 = Modeller(cif.topology, cif.positions)
    if args.free:
        m0.delete([r for ch in m0.topology.chains() if ch.id != args.binder_chain
                   for r in ch.residues()])
    tmp = args.out + "_raw.pdb"
    with open(tmp, "w") as f:
        PDBFile.writeFile(m0.topology, m0.positions, f)

    fixer = PDBFixer(filename=tmp)
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        ch = chains[key[0]]; n = len(list(ch.residues()))
        if key[1] == 0 or key[1] == n:         # drop terminal gaps (no loop-building at ends)
            del fixer.missingResidues[key]
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()                     # adds OXT etc.

    model = Modeller(fixer.topology, fixer.positions)

    # --- disulfides: SG-SG < 2.5 A -> add bond + mark CYX ---
    pos = np.array(model.positions.value_in_unit(nanometer))
    sg = [(r, next(a for a in r.atoms() if a.name == "SG"))
          for ch in model.topology.chains() for r in ch.residues()
          if r.name == "CYS" and any(a.name == "SG" for a in r.atoms())]
    cyx = set()
    for i in range(len(sg)):
        for j in range(i+1, len(sg)):
            d = np.linalg.norm(pos[sg[i][1].index] - pos[sg[j][1].index])
            if d < 0.25:
                model.topology.addBond(sg[i][1], sg[j][1]); cyx.add(sg[i][0]); cyx.add(sg[j][0])

    his = [(ch.id, r.id, r.name) for ch in model.topology.chains() for r in ch.residues()
           if r.name in ("HIS","HID","HIE","HIP")]
    print(f"[build] free={args.free} chains={[c.id for c in model.topology.chains()]} "
          f"His={his} disulfides={len(cyx)//2}")

    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    variants = []
    for ch in model.topology.chains():
        for r in ch.residues():
            if r in cyx: variants.append("CYX")
            elif r.name in ("HIS","HID","HIE","HIP"):
                if ch.id == args.binder_chain and int(r.id) in prot: variants.append("HIP")
                else: variants.append("HIE")
            else: variants.append(None)
    model.addHydrogens(ff, pH=7.0, variants=variants)
    print(f"[build] HIP switch-His applied = {sorted(prot)}  (count={sum(1 for v in variants if v=='HIP')})")

    if args.box_from:
        V = np.loadtxt(args.box_from)
        bv = [Vec3(*V[0]), Vec3(*V[1]), Vec3(*V[2])] * nanometer
        model.addSolvent(ff, model="tip3p", boxVectors=bv, ionicStrength=0.15*molar, neutralize=True)
    else:
        model.addSolvent(ff, model="tip3p", padding=args.pad*nanometer, ionicStrength=0.15*molar, neutralize=True)
    natom = model.topology.getNumAtoms()
    bv = model.topology.getPeriodicBoxVectors().value_in_unit(nanometer)
    np.savetxt(args.out + "_box.txt", np.array(bv))
    print(f"[build] solvated atoms = {natom}  box(nm) = {[round(bv[k][k],2) for k in range(3)]}")

    system = ff.createSystem(model.topology, nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds, rigidWater=True, hydrogenMass=4*amu)
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin, 25))
    integ = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 4*femtosecond)
    plat = Platform.getPlatformByName("CUDA")
    sim = Simulation(model.topology, system, integ, plat, {"DeviceIndex": args.device, "Precision": "mixed"})
    sim.context.setPositions(model.positions)
    t0 = time.time(); sim.minimizeEnergy(maxIterations=5000)
    print(f"[min] {time.time()-t0:.0f}s  E={sim.context.getState(getEnergy=True).getPotentialEnergy()}")
    sim.context.setVelocitiesToTemperature(300*kelvin)
    nsteps = int(args.equil_ps*1000/4)
    sim.reporters.append(StateDataReporter(sys.stdout, max(nsteps//4,1), step=True,
                         potentialEnergy=True, temperature=True, density=True, speed=True))
    t0 = time.time(); sim.step(nsteps); wall = time.time()-t0
    print(f"[equil] {args.equil_ps:.0f}ps in {wall:.0f}s => {(args.equil_ps/1000)/(wall/86400):.1f} ns/day  ({natom} atoms)")

    st = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    with open(args.out+".pdb","w") as f: PDBFile.writeFile(model.topology, st.getPositions(), f)
    with open(args.out+"_system.xml","w") as f: f.write(XmlSerializer.serialize(system))
    with open(args.out+"_state.xml","w") as f: f.write(XmlSerializer.serialize(st))
    json.dump({"natom":natom,"switch_his":sorted(prot),"free":args.free,
               "box_nm":[[bv[i][j] for j in range(3)] for i in range(3)]},
              open(args.out+"_meta.json","w"), indent=0)
    print(f"[save] {args.out}.{{pdb,_system.xml,_state.xml,_box.txt,_meta.json}}")

if __name__ == "__main__":
    main()
