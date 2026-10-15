#!/usr/bin/env python
"""One leg of the pH-switch alchemical ddG by TI, with WELL-TEMPERED METADYNAMICS on the
His(+)<->carboxylate salt-bridge coordination CV to break bound-leg basin-trapping.

Why: plain TI gives a reproducible SIGN for ddG_switch (the bound-free double-difference cancels the
correlated slow-reorganization noise) but per-leg seed spread of 4-8 kcal/mol makes MAGNITUDE
unrankable. The slow DOF is the His(+) salt-bridge state. We bias a coordination-number CV so the
dynamics crosses salt-bridge basins, then reweight <dU/dlam> to the UNBIASED ensemble.

Key identity: the metadynamics bias V(s) depends only on the CV s, NOT on lam. In the finite-difference
dU/dlam = [U(lam+d) - U(lam-d)]/(2d) the bias cancels exactly, so each per-config dU/dlam is unbiased;
we recover the unbiased ensemble average by importance weights w_i = exp(+V(s_i)/kT). For well-tempered
metadynamics at (quasi-)convergence, V(s) = -(1 - 1/gamma) * F(s), with F from Metadynamics.getFreeEnergy()
and gamma = biasFactor. A FRESH bias is grown per lam window (F(s) is lam-dependent).

Reuses cached equilibrated systems in runs/fe/sys/. Native OpenMM metadynamics (no PLUMED / no MPI).
Run with envs/md2 + OPENMM_CUDA_COMPILER=cuda-12.3 nvcc. One leg per GPU process.
"""
import argparse, json, math, numpy as np
from openmm import (XmlSerializer, NonbondedForce, MonteCarloBarostat, LangevinMiddleIntegrator,
                    Platform, CustomBondForce)
from openmm.app import PDBFile, ForceField, Simulation
from openmm.app.metadynamics import Metadynamics, BiasVariable
from openmm.unit import kelvin, picosecond, femtosecond, kilojoule_per_mole, nanometer

KT = 0.00831446262 * 300.0  # kJ/mol at 300 K
CARBOXYL_O = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "GLH": ("OE1", "OE2"), "ASH": ("OD1", "OD2")}

def load_ff_charges():
    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    return {v: {a.name: a.parameters["charge"] for a in ff._templates[v].atoms} for v in ("HIE", "HIP")}

def apply_alchemy(nbf, top, prot, binder_chain, qH):
    """HIE(lam=0) -> HIP(lam=1) via NonbondedForce per-particle offsets on global 'lam'. Returns applied His ids."""
    nbf.addGlobalParameter("lam", 0.0)
    tot_dq = 0.0; applied = []
    for res in top.residues():
        if res.chain.id == binder_chain and res.name in ("HIS","HID","HIE","HIP") and int(res.id) in prot:
            for a in res.atoms():
                q_cur, sig, eps = nbf.getParticleParameters(a.index)
                eps_v = eps.value_in_unit(eps.unit)
                q_hie = qH["HIE"].get(a.name, 0.0)
                q_hip = qH["HIP"].get(a.name, q_cur.value_in_unit(q_cur.unit))
                if a.name == "HD1":
                    nbf.setParticleParameters(a.index, 0.0, sig, 0.0)
                    nbf.addParticleParameterOffset("lam", a.index, q_hip, 0.0, eps_v)
                else:
                    nbf.setParticleParameters(a.index, q_hie, sig, eps)
                    nbf.addParticleParameterOffset("lam", a.index, q_hip - q_hie, 0.0, 0.0)
                tot_dq += (q_hip - q_hie)
            applied.append(int(res.id))
    assert abs(tot_dq - len(applied)) < 0.05, f"net dq {tot_dq} != +{len(applied)}"
    return applied

def build_cv_pairs(top, positions_nm, prot, binder_chain, cutoff_nm=0.80):
    """His(+) ring N (ND1/NE2 of switch-His on binder) x carboxylate O (any chain) within cutoff at start."""
    hisN = []; carbO = []
    for res in top.residues():
        if res.chain.id == binder_chain and res.name in ("HIS","HID","HIE","HIP") and int(res.id) in prot:
            for a in res.atoms():
                if a.name in ("ND1", "NE2"): hisN.append(a.index)
        onames = CARBOXYL_O.get(res.name)
        for a in res.atoms():
            if (onames and a.name in onames) or a.name == "OXT":
                carbO.append(a.index)
    P = np.asarray(positions_nm)
    pairs = []
    for i in hisN:
        for j in carbO:
            if np.linalg.norm(P[i] - P[j]) < cutoff_nm:
                pairs.append((i, j))
    return pairs, len(hisN), len(carbO)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sys", required=True)
    ap.add_argument("--switch-his", required=True)
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--nwin", type=int, default=11)
    ap.add_argument("--equil-ns", type=float, default=0.5, help="per-window bias-fill before sampling")
    ap.add_argument("--prod-ns", type=float, default=1.0, help="per-window production (sampled)")
    ap.add_argument("--sample-ps", type=float, default=0.5)
    ap.add_argument("--meta-height", type=float, default=1.0, help="Gaussian height kJ/mol")
    ap.add_argument("--meta-freq", type=int, default=500, help="deposition stride (steps)")
    ap.add_argument("--bias-factor", type=float, default=8.0)
    ap.add_argument("--bias-width", type=float, default=0.3)
    ap.add_argument("--cv-max", type=float, default=8.0)
    ap.add_argument("--r0", type=float, default=0.33, help="coordination switch radius nm")
    ap.add_argument("--device", default="0")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    prot = set(int(x) for x in args.switch_his.split(",") if x.strip())
    qH = load_ff_charges()

    pdb = PDBFile(args.sys + ".pdb")
    top = pdb.topology
    st0 = XmlSerializer.deserialize(open(args.sys + "_state.xml").read())
    pos0 = st0.getPositions(asNumpy=True).value_in_unit(nanometer)
    box0 = st0.getPeriodicBoxVectors()
    vel0 = st0.getVelocities()
    pairs, nN, nO = build_cv_pairs(top, pos0, prot, args.binder_chain, cutoff_nm=0.80)
    use_meta = len(pairs) > 0
    plat = Platform.getPlatformByName("CUDA")
    gamma = args.bias_factor
    wt_coef = (1.0 - 1.0 / gamma)  # V(s) = -wt_coef * F(s)

    def make_sim(lam, start_pos):
        """Fresh system + fresh metad bias (bias is lam-specific). Positions carry over across windows
        (sequential relaxation, matching the validated plain-TI protocol); velocities reinit each window."""
        system = XmlSerializer.deserialize(open(args.sys + "_system.xml").read())
        for i in reversed(range(system.getNumForces())):
            if isinstance(system.getForce(i), MonteCarloBarostat):
                system.removeForce(i)
        nbf = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        applied = apply_alchemy(nbf, top, prot, args.binder_chain, qH)
        meta = None
        if use_meta:
            cvf = CustomBondForce("1.0/(1.0+(r/%f)^6)" % args.r0)
            for (i, j) in pairs:
                cvf.addBond(i, j, [])
            bv = BiasVariable(cvf, 0.0, args.cv_max, args.bias_width, False, gridWidth=121)
            meta = Metadynamics(system, [bv], 300 * kelvin, gamma, args.meta_height * kilojoule_per_mole,
                                args.meta_freq)
        integ = LangevinMiddleIntegrator(300 * kelvin, 1.0 / picosecond, 4 * femtosecond)
        integ.setRandomNumberSeed(args.seed + int(lam * 997) + 1)
        sim = Simulation(top, system, integ, plat, {"DeviceIndex": args.device, "Precision": "mixed"})
        sim.context.setPositions(start_pos)
        sim.context.setPeriodicBoxVectors(*box0)
        sim.context.setParameter("lam", float(lam))
        sim.context.setVelocitiesToTemperature(300 * kelvin, args.seed + int(lam * 1000) + 7)
        return sim, meta, applied

    def potE(sim):
        return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)

    lambdas = np.linspace(0.0, 1.0, args.nwin)
    eq = int(args.equil_ns * 1e6 / 4); prod = int(args.prod_ns * 1e6 / 4)
    stride = max(int(args.sample_ps * 1000 / 4), 1)
    applied_ref = None
    means = []; sems = []; cvdiag = []
    start_pos = pos0
    for l in lambdas:
        lp = min(float(l) + 0.01, 1.0); lm = max(float(l) - 0.01, 0.0); hh = lp - lm
        sim, meta, applied = make_sim(l, start_pos)
        applied_ref = applied
        # bias-fill
        if meta is not None: meta.step(sim, eq)
        else: sim.step(eq)
        dudl = []; cvs = []
        nsamp = max(prod // stride, 1)
        for _ in range(nsamp):
            if meta is not None: meta.step(sim, stride)
            else: sim.step(stride)
            if meta is not None:
                cvs.append(float(meta.getCollectiveVariables(sim)[0]))
            sim.context.setParameter("lam", lp); Up = potE(sim)
            sim.context.setParameter("lam", lm); Um = potE(sim)
            sim.context.setParameter("lam", float(l))
            dudl.append((Up - Um) / hh)
        dudl = np.array(dudl)
        if meta is not None and len(cvs):
            F = meta.getFreeEnergy().value_in_unit(kilojoule_per_mole)  # grid kJ/mol
            F = np.nan_to_num(F - np.nanmin(F), nan=0.0, posinf=0.0, neginf=0.0)
            cvs = np.array(cvs)
            idx = np.clip(((cvs - 0.0) / (args.cv_max - 0.0) * (len(F) - 1)).round().astype(int), 0, len(F) - 1)
            Vbias = -wt_coef * F[idx]                 # V(s_i) kJ/mol
            w = np.exp((Vbias - Vbias.max()) / KT)     # importance weights (shifted for stability)
            w /= w.sum()
            m = float(np.sum(w * dudl))
            # weighted block SEM (5 blocks)
            nb = min(5, len(dudl)); bm = []
            for bi in range(nb):
                sl = slice(bi * len(dudl) // nb, (bi + 1) * len(dudl) // nb)
                ws = w[sl]
                bm.append(float(np.sum(ws * dudl[sl]) / ws.sum()) if ws.sum() > 0 else float("nan"))
            bm = np.array(bm)
            sem = float(np.nanstd(bm, ddof=1) / math.sqrt(nb)) if nb > 1 else float("nan")
            cvdiag.append({"lam": round(float(l), 2), "cv_mean": round(float(cvs.mean()), 2),
                           "cv_min": round(float(cvs.min()), 2), "cv_max": round(float(cvs.max()), 2),
                           "cv_sd": round(float(cvs.std()), 2), "n_eff": round(float(1.0 / np.sum(w**2)), 1)})
        else:
            m = float(dudl.mean())
            nb = min(5, len(dudl)); blk = np.array([b.mean() for b in np.array_split(dudl, nb)])
            sem = float(blk.std(ddof=1) / math.sqrt(nb)) if nb > 1 else float("nan")
        means.append(m); sems.append(sem)
        extra = cvdiag[-1] if cvdiag else {}
        print(f"   lam={l:.2f} <dU/dl>={m:8.1f} +/- {sem:5.1f} kJ/mol (n={len(dudl)}) cv={extra.get('cv_mean','-')} neff={extra.get('n_eff','-')}", flush=True)
        start_pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)  # carryover
        del sim

    means = np.array(means); sems = np.array(sems); lam = lambdas
    w = np.zeros_like(lam)
    w[1:-1] = (lam[2:] - lam[:-2]) / 2; w[0] = (lam[1] - lam[0]) / 2; w[-1] = (lam[-1] - lam[-2]) / 2
    dG_kJ = float(np.sum(w * means)); err_kJ = float(np.sqrt(np.sum((w * np.nan_to_num(sems)) ** 2)))
    out = {"leg": args.sys.split("/")[-1], "dG_kcal": dG_kJ / 4.184, "err_kcal": err_kJ / 4.184,
           "switch_his": applied_ref, "nwin": args.nwin, "prod_ns": args.prod_ns, "equil_ns": args.equil_ns,
           "seed": args.seed, "method": "meta-TI" if use_meta else "plain-TI", "n_cv_pairs": len(pairs),
           "n_hisN": nN, "n_carbO": nO, "dudl_kJ": [round(x, 1) for x in means.tolist()],
           "cv_diag": cvdiag}
    json.dump(out, open(args.out + ".json", "w"), indent=2)
    print(f"[RESULT] {out['leg']} ({out['method']}, {len(pairs)} CV pairs) dG(HIE->HIP)= {out['dG_kcal']:+.2f} +/- {out['err_kcal']:.2f} kcal/mol", flush=True)

if __name__ == "__main__":
    main()
