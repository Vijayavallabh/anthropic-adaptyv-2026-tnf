#!/usr/bin/env python
"""One leg of the pH-switch alchemical ddG by thermodynamic integration (TI).

Perturbs switch-His HIE(lam=0)->HIP(lam=1) via NonbondedForce per-particle offsets on global
'lam'; OpenMM computes dU/dlam natively (addEnergyParameterDerivative). Independent lam windows
(sequential config carryover), <dU/dlam> per window, integrate -> dG(0->1) for this leg.
ddG_switch = dG_bound - dG_free (combine step). 1-GPU-efficient; parallelize legs across GPUs.
Run with envs/md2 + OPENMM_CUDA_COMPILER=cuda-12.3 nvcc.
"""
import argparse, json, numpy as np
from openmm import (XmlSerializer, NonbondedForce, MonteCarloBarostat, LangevinMiddleIntegrator,
                    Platform)
from openmm.app import PDBFile, ForceField, Simulation
from openmm.unit import kelvin, picosecond, femtosecond, kilojoule_per_mole

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sys", required=True)
    ap.add_argument("--switch-his", required=True)
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--nwin", type=int, default=11)
    ap.add_argument("--equil-ns", type=float, default=0.5)
    ap.add_argument("--prod-ns", type=float, default=2.0)
    ap.add_argument("--sample-ps", type=float, default=1.0, help="dU/dlam sampling stride")
    ap.add_argument("--device", default="0")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    prot = set(int(x) for x in args.switch_his.split(",") if x.strip())

    system = XmlSerializer.deserialize(open(args.sys+"_system.xml").read())
    pdb = PDBFile(args.sys+".pdb")
    st = XmlSerializer.deserialize(open(args.sys+"_state.xml").read())
    top = pdb.topology
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i), MonteCarloBarostat):
            system.removeForce(i)
    nbf = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]

    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    qH = {v: {a.name: a.parameters["charge"] for a in ff._templates[v].atoms} for v in ("HIE", "HIP")}
    nbf.addGlobalParameter("lam", 0.0)
    # Only atoms whose partial charge actually changes on HIE->HIP (the imidazole ring + the appearing
    # HD1 proton) are perturbed. Backbone / CB / terminal atoms (e.g. a C-terminal OXT or N-terminal
    # H1/H2/H3) carry an identical charge in both protonation states, so they are LEFT UNTOUCHED. This
    # keeps a terminal switch-His correct — its FF terminal-variant charges differ from the internal
    # HIE/HIP template, so overwriting them with the template value corrupts the leg (and trips the net-dq
    # assertion) — and is byte-identical for internal His (those atoms were previously no-op overwrites
    # with a zero lambda-offset). Verified byte-identical on des_820; unblocks C-terminal His (des_875).
    CHG = {n for n in set(qH["HIE"]) | set(qH["HIP"])
           if abs(qH["HIP"].get(n, 0.0) - qH["HIE"].get(n, 0.0)) > 1e-6}
    tot_dq = 0.0; nmod = 0; applied = []
    for res in top.residues():
        if res.chain.id == args.binder_chain and res.name in ("HIS","HID","HIE","HIP") and int(res.id) in prot:
            for a in res.atoms():
                if a.name not in CHG:
                    continue
                q_cur, sig, eps = nbf.getParticleParameters(a.index)
                eps_v = eps.value_in_unit(eps.unit)
                q_hie = qH["HIE"].get(a.name, 0.0); q_hip = qH["HIP"].get(a.name, 0.0)
                if a.name == "HD1":
                    nbf.setParticleParameters(a.index, 0.0, sig, 0.0)
                    nbf.addParticleParameterOffset("lam", a.index, q_hip, 0.0, eps_v)
                else:
                    nbf.setParticleParameters(a.index, q_hie, sig, eps)
                    nbf.addParticleParameterOffset("lam", a.index, q_hip - q_hie, 0.0, 0.0)
                tot_dq += (q_hip - q_hie); nmod += 1
            applied.append(int(res.id))
    assert abs(tot_dq - len(applied)) < 0.05, f"net dq {tot_dq:+.3f} != +{len(applied)}"
    print(f"[TI-FD] leg={args.sys.split('/')[-1]} His={applied} atoms={nmod} dq={tot_dq:+.2f} seed={args.seed}")

    integ = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 4*femtosecond)
    integ.setRandomNumberSeed(args.seed)
    plat = Platform.getPlatformByName("CUDA")
    sim = Simulation(top, system, integ, plat, {"DeviceIndex": args.device, "Precision": "mixed"})
    sim.context.setPositions(st.getPositions())
    sim.context.setPeriodicBoxVectors(*st.getPeriodicBoxVectors())
    sim.context.setVelocities(st.getVelocities())

    lambdas = np.linspace(0.0, 1.0, args.nwin)
    eq = int(args.equil_ns*1e6/4); prod = int(args.prod_ns*1e6/4); stride = max(int(args.sample_ps*1000/4),1)
    dl = 0.01  # finite-difference step for dU/dlam
    def potE():
        return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    means = []; sems = []
    for l in lambdas:
        lp = min(float(l)+dl, 1.0); lm = max(float(l)-dl, 0.0); hh = lp - lm
        sim.context.setParameter("lam", float(l))
        sim.context.setVelocitiesToTemperature(300*kelvin, args.seed + int(l*1000))
        sim.step(eq)
        s = []
        nsamp = max(prod//stride, 1)
        for _ in range(nsamp):
            sim.step(stride)
            sim.context.setParameter("lam", lp); Up = potE()
            sim.context.setParameter("lam", lm); Um = potE()
            sim.context.setParameter("lam", float(l))
            s.append((Up - Um)/hh)  # kJ/mol  (central diff; one-sided at endpoints)
        s = np.array(s)
        # block SEM (5 blocks)
        nb = min(5, len(s)); blk = np.array([b.mean() for b in np.array_split(s, nb)])
        means.append(float(s.mean())); sems.append(float(blk.std(ddof=1)/np.sqrt(nb)) if nb > 1 else float("nan"))
        print(f"   lam={l:.2f}  <dU/dlam>={s.mean():8.1f} +/- {sems[-1]:5.1f} kJ/mol  (n={len(s)})")

    means = np.array(means); sems = np.array(sems); lam = lambdas
    # trapezoid weights (used for both the integral and error propagation)
    w = np.zeros_like(lam)
    w[1:-1] = (lam[2:]-lam[:-2])/2; w[0] = (lam[1]-lam[0])/2; w[-1] = (lam[-1]-lam[-2])/2
    dG_kJ = float(np.sum(w*means))
    err_kJ = float(np.sqrt(np.sum((w*np.nan_to_num(sems))**2)))
    dG = dG_kJ/4.184; err = err_kJ/4.184
    out = {"leg": args.sys.split("/")[-1], "dG_kcal": dG, "err_kcal": err, "switch_his": applied,
           "nwin": args.nwin, "prod_ns": args.prod_ns, "equil_ns": args.equil_ns, "seed": args.seed,
           "total_dq": float(tot_dq), "dudl_kJ": [round(x,1) for x in means.tolist()],
           "lambdas": [round(x,3) for x in lam.tolist()]}
    json.dump(out, open(args.out+".json","w"), indent=2)
    print(f"[RESULT] {out['leg']} dG(HIE->HIP)= {dG:+.2f} +/- {err:.2f} kcal/mol")

if __name__ == "__main__":
    main()
