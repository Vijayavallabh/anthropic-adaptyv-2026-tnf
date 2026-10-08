#!/usr/bin/env python3
"""WAVE 24: two-pH MD for a CATION-pi pH switch (the ECC-robust alternative to the His+...Asp- salt bridge).

Measures, per frame, the stacking geometry between each binder His imidazole ring and its target conserved-Phe
ring (ring-centroid distance + inter-plane angle). A cation-pi contact = the His+ ring sitting over/near the
aromatic face. We run the binder His in HIP (protonated, tumor pH 6.5) vs HIE (neutral, blood pH 7.4) and
report occupancy in each state; d_occ = occ(HIP) - occ(HIE) is the per-contact switch, and for a DUAL design
occ_both / d_occ_both is the cooperative two-contact signal.

Crucially this script supports --charge-scale (ECC): scaling ionizable-group charges by 0.8 is the stricter
robustness test that COLLAPSED our His+...Asp- salt bridge (d_occ 0.87 -> 0.003). A cation-pi partner is a
NEUTRAL aromatic, so there is no counter-charge to scale; the contact should survive ECC far better. The
head-to-head that decides "far better on SELECTIVITY" is: ECC-corrected d_occ(cation-pi) vs ECC-corrected
d_occ(salt bridge), both measured here on the same footing.

Reuses md_pH_switch's OpenMM/forcefield machinery (chains, CA restraints, integrator, HMR, ECC scaler).
"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
import md_pH_switch as M
app, unit = M.app, M.unit

HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
PHE_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
TRP_RING = ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")
CC_CUT = 5.5   # A; His-ring..Phe-ring centroid distance counted as a cation-pi contact (incl. offset stacks)


def build_his_variants(pdb, his_ids, variant, binder_idx, charge_scale=1.0):
    """set the listed binder His (by residue id) to HIP or HIE; optional ECC charge scaling."""
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    residues = list(modeller.topology.residues())
    variants = [None] * len(residues)
    want = {str(h) for h in his_ids}
    nset = 0
    for i, r in enumerate(residues):
        if r.chain.index == binder_idx and r.name.startswith("HI") and str(r.id) in want:
            variants[i] = variant; nset += 1
    if nset != len(his_ids):
        M.log("  WARN: set %d/%d His variants" % (nset, len(his_ids)))
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds,
                             hydrogenMass=4.0 * unit.amu)
    if charge_scale != 1.0:
        M.scale_ionizable_charges(system, modeller.topology, charge_scale)
    return ff, modeller, system


def ring_idx(topology, chain_index, resid, names):
    for r in topology.residues():
        if r.chain.index == chain_index and str(r.id) == str(resid):
            return [a.index for a in r.atoms() if a.name in names]
    return []


def centroid(p, idx):
    return p[idx].mean(0) if idx else None


def run_state(pdb, pairs, variant, binder_idx, target_idx, device, equil_ps, prod_ps,
              interval_ps, replica, restraint_k, charge_scale):
    his_ids = [hp[0] for hp in pairs]
    ff, modeller, system = build_his_variants(pdb, his_ids, variant, binder_idx, charge_scale)
    M.add_ca_restraints(system, modeller, restraint_k)
    # ring atom indices for each (His, Phe) pair
    idxs = []
    for his_id, phe_id, phe_names in pairs:
        hi = ring_idx(modeller.topology, binder_idx, his_id, HIS_RING)
        pi = ring_idx(modeller.topology, target_idx, phe_id, phe_names)
        idxs.append((hi, pi))
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
    dists = [[] for _ in pairs]
    for _ in range(nsteps // stride):
        sim.step(stride)
        p = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))
        for k, (hi, pi) in enumerate(idxs):
            if hi and pi:
                dists[k].append(float(np.linalg.norm(centroid(p, hi) - centroid(p, pi))))
    out = []
    both_mask = None
    for k, d in enumerate(dists):
        d = np.array(d)
        occ = float(np.mean(d < CC_CUT)) if len(d) else None
        out.append(dict(occ=round(occ, 3) if occ is not None else None,
                        mindist=round(float(d.min()), 2) if len(d) else None,
                        meandist=round(float(d.mean()), 2) if len(d) else None))
        m = (d < CC_CUT) if len(d) else None
        both_mask = m if both_mask is None else (both_mask & m if m is not None else both_mask)
    occ_both = round(float(np.mean(both_mask)), 3) if (both_mask is not None and len(pairs) > 1) else None
    return dict(variant=variant, replica=replica, contacts=out, occ_both=occ_both)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--pairs", required=True,
                    help="comma list of hisResid:pheResid[:TRP] in the af2pred numbering, e.g. 37:206,54:174")
    ap.add_argument("--binder-chain", default="A"); ap.add_argument("--target-chain", default="B")
    ap.add_argument("--out", required=True); ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--equil-ps", type=float, default=200.0); ap.add_argument("--prod-ps", type=float, default=2000.0)
    ap.add_argument("--interval-ps", type=float, default=10.0); ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--restraint-k", type=float, default=2.0)
    ap.add_argument("--charge-scale", type=float, default=1.0, help="ECC: scale ionizable charges (0.8 = test)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    name = os.path.splitext(os.path.basename(a.pdb))[0]
    # resolve chain indices (order of appearance)
    from pdbfixer import PDBFixer
    fx = PDBFixer(filename=a.pdb)
    chain_ids = [c.id for c in fx.topology.chains()]
    binder_idx = chain_ids.index(a.binder_chain); target_idx = chain_ids.index(a.target_chain)
    pairs = []
    for tok in a.pairs.split(","):
        parts = tok.split(":")
        his_id, phe_id = int(parts[0]), int(parts[1])
        names = TRP_RING if (len(parts) > 2 and parts[2].upper() == "TRP") else PHE_RING
        pairs.append((his_id, phe_id, names))
    M.log("=== MD cation-pi:", name, "pairs", a.pairs, "charge_scale", a.charge_scale, "device", a.device)
    states = []
    for rep in range(a.replicas):
        for variant in ("HIP", "HIE"):
            states.append(run_state(a.pdb, pairs, variant, binder_idx, target_idx, a.device,
                                     a.equil_ps, a.prod_ps, a.interval_ps, rep, a.restraint_k, a.charge_scale))

    def avg_contact(k, field, var):
        xs = [s["contacts"][k][field] for s in states if s["variant"] == var and s["contacts"][k][field] is not None]
        return round(float(np.mean(xs)), 3) if xs else None

    def avg_both(var):
        xs = [s["occ_both"] for s in states if s["variant"] == var and s["occ_both"] is not None]
        return round(float(np.mean(xs)), 3) if xs else None

    res = dict(name=name, pdb=os.path.abspath(a.pdb), pairs=a.pairs, charge_scale=a.charge_scale, states=states,
               contacts=[], n_contacts=len(pairs))
    for k, (his_id, phe_id, _) in enumerate(pairs):
        hip, hie = avg_contact(k, "occ", "HIP"), avg_contact(k, "occ", "HIE")
        res["contacts"].append(dict(his=his_id, phe=phe_id, occ_HIP=hip, occ_HIE=hie,
                                    d_occ=round((hip or 0) - (hie or 0), 3)))
    if len(pairs) > 1:
        bh, be = avg_both("HIP"), avg_both("HIE")
        res["occ_both_HIP"] = bh; res["occ_both_HIE"] = be
        res["d_occ_both"] = round((bh or 0) - (be or 0), 3)
        res["coop_switch"] = (bh or 0) > 0.3 and (be or 0) < 0.15
    # a cation-pi switch is clean if every contact is protonation-gated (d_occ>0.2, occ_HIE<0.3)
    res["switch_clean"] = all(c["d_occ"] > 0.2 and (c["occ_HIE"] or 0) < 0.3 for c in res["contacts"])
    json.dump(res, open(os.path.join(a.out, name + ".json"), "w"), indent=1)
    msg = " | ".join("His%s-Phe%s occ %.2f/%.2f d=%.2f" % (c["his"], c["phe"], c["occ_HIP"] or 0,
                     c["occ_HIE"] or 0, c["d_occ"]) for c in res["contacts"])
    M.log("  [scale=%.2f] %s  both=%s  clean=%s" % (a.charge_scale, msg, res.get("d_occ_both"), res["switch_clean"]))


if __name__ == "__main__":
    main()
