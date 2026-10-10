#!/usr/bin/env python3
"""BOLD BET D: steered-MD (SMD) binder<->target SEPARATION work in BOTH His protonation states = direct
FUNCTIONAL pH-release evidence (beyond the thermodynamic PROPKA/xTB linkage proxy). Pull the binder COM away
from the (restrained) target at constant velocity in implicit solvent; accumulate Jarzynski pulling work W.
Release signal: W(HIE, pH7.4 neutral) - W(HIP, pH6.0 protonated) > 0  (neutral binds tighter, protonated
lets go more easily). GLN-isostere run = His-removed binding baseline (His-specificity). Reuses
md_release_tnf building blocks (same implicit-solvent system, target-restrained). envs/md2/bin/python.
Usage: md_pull_tnf.py --pdb X --his-pos N[,M..] --out dir [--variants HIE,HIP,GLN] [--replicas 3]
       [--pull-nm 1.5 --pull-ps 500 --equil-ps 200 --interval-ps 1 --k 2000 --device G]"""
import sys, os, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import md_release_tnf as M
from openmm import CustomCentroidBondForce, LangevinMiddleIntegrator, unit
from openmm import app

def run_pull(pdb, variant, his_list, a, device, replica):
    mut = (variant == "GLN")
    ff, modeller, system, target, binder, sws = M.build(pdb, variant, his_list, a.binder_chain, mut)
    bh = [at.index for at in binder.atoms() if at.element is not None and at.element.symbol != "H"]
    th = [at.index for at in target.atoms() if at.element is not None and at.element.symbol != "H"]
    M.add_ca_restraints(system, modeller, target, a.restraint_k)          # target rigid, binder free
    pull = CustomCentroidBondForce(2, "0.5*kpull*(distance(g1,g2)-d0)^2")  # 'kpull' not 'k' (k is used by add_ca_restraints)
    pull.addGlobalParameter("kpull", a.k); pull.addGlobalParameter("d0", 0.0)  # kpull kJ/mol/nm^2, d0 nm
    g1 = pull.addGroup(bh, [1.0]*len(bh)); g2 = pull.addGroup(th, [1.0]*len(th))  # geometric centroids
    pull.addBond([g1, g2], []); system.addForce(pull)
    integ = LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond, 4.0*unit.femtosecond)
    integ.setRandomNumberSeed(9000 + replica)
    plat, props, pname = M.pick_platform(device)
    sim = app.Simulation(modeller.topology, system, integ, plat, props)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    def com_dist():
        p = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer))
        return float(np.linalg.norm(p[bh].mean(0) - p[th].mean(0)))
    d0 = com_dist(); sim.context.setParameter("d0", d0)
    sim.context.setVelocitiesToTemperature(310*unit.kelvin, 9000 + replica)
    sim.step(int(a.equil_ps/0.004))
    d0 = com_dist(); sim.context.setParameter("d0", d0)                    # recenter after equil
    n_str = int(a.pull_ps/a.interval_ps); steps_per = max(1, int(a.interval_ps/0.004))
    dd = a.pull_nm/n_str; W = 0.0; traj = []
    for _ in range(n_str):
        sim.step(steps_per)
        dc = com_dist()
        W += a.k*(d0 - dc)*dd                                             # kJ/mol (Jarzynski pulling work)
        traj.append(dc); d0 += dd; sim.context.setParameter("d0", d0)
    return {"variant": variant, "replica": replica, "W_kcal": round(W/4.184, 2),
            "d_start": round(traj[0], 2), "d_end": round(traj[-1], 2), "platform": pname}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True); ap.add_argument("--his-pos", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--variants", default="HIE,HIP,GLN"); ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--pull-nm", type=float, default=1.5); ap.add_argument("--pull-ps", type=float, default=500.0)
    ap.add_argument("--equil-ps", type=float, default=200.0); ap.add_argument("--interval-ps", type=float, default=1.0)
    ap.add_argument("--k", type=float, default=2000.0); ap.add_argument("--restraint-k", type=float, default=5.0)
    ap.add_argument("--device", type=int, default=0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    his_list = [int(x) for x in str(a.his_pos).split(",")]
    variants = a.variants.split(","); name = os.path.splitext(os.path.basename(a.pdb))[0]
    runs = []
    for v in variants:
        for r in range(a.replicas):
            try:
                res = run_pull(a.pdb, v, his_list, a, a.device, r); runs.append(res)
                M.log(f"[pull] {name} {v} rep{r}: W={res['W_kcal']} kcal d {res['d_start']}->{res['d_end']}")
            except Exception as e:
                M.log(f"[pull] {name} {v} rep{r} FAIL: {str(e)[:160]}")
    def meanW(v):
        ws = [x["W_kcal"] for x in runs if x["variant"] == v]; return (round(float(np.mean(ws)),2), round(float(np.std(ws)),2), len(ws)) if ws else (None,None,0)
    wHIE, sHIE, nHIE = meanW("HIE"); wHIP, sHIP, nHIP = meanW("HIP"); wGLN, sGLN, nGLN = meanW("GLN")
    dW = round(wHIE - wHIP, 2) if (wHIE is not None and wHIP is not None) else None   # >0 = release at low pH
    out = {"name": name, "his_pos": his_list, "pdb": os.path.abspath(a.pdb),
           "W_HIE": wHIE, "sd_HIE": sHIE, "W_HIP": wHIP, "sd_HIP": sHIP, "W_GLN": wGLN, "sd_GLN": sGLN,
           "dW_release": dW, "replicas": a.replicas, "runs": runs,
           "release_sign": (dW is not None and dW > 0)}
    json.dump(out, open(f"{a.out}/{name}_pull.json", "w"), indent=1)
    M.log(f"[pull] {name}: W_HIE={wHIE}+/-{sHIE} W_HIP={wHIP}+/-{sHIP} W_GLN={wGLN}  dW_release={dW} "
          f"({'RELEASE sign OK' if out['release_sign'] else 'no release sign'})")
    print(f"PULL_DONE {name} dW_release={dW}")

if __name__ == "__main__":
    main()
