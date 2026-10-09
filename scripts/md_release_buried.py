#!/usr/bin/env python3
"""Explicit two-state MD selectivity gate for the TNF pH-RELEASE switch (Challenge-2), WITH CONTROLS.

Mechanism: an interface binder-His sits next to a CONSERVED target Lys (K90/K98).
  * HIE = neutral imidazole (0)    -> plasma pH 7.4  -> binder HELD      (ON  state)
  * HIP = protonated imidazole (+1)-> endosome pH 6.0 -> His+ repels Lys+ -> RELEASE (OFF state)

A static PROPKA linkage proves DIRECTION, not SELECTIVITY (a constitutive binder, held by other
contacts regardless of His state, scores linkage<0 yet binds in both). Only an explicit two-state
occupancy delta catches it.  BUT a few-ns implicit-solvent MD cannot show full dissociation and the
delta can be noise.  So this gate is only trustworthy WITH CONTROLS, and must not change a submission
until the controls separate:
  * CAUSAL DUD (isostere): mutate the switch His -> Gln (neutral amide, CANNOT titrate). A real switch
    shows ON(HIE) ~= GLN (both neutral, binder held) and OFF(HIP) deviating. If HIP deviates no more
    than GLN does from HIE, the "signal" is pose noise, not a switch.
  * NOISE FLOOR: replicate spread of each state (HIE rep-to-rep) = test-retest noise; the HIE->HIP
    signal must exceed it.
Primary signal is the DIRECT His-Lys distance (Coulomb repulsion responds in <1 ns = fast, robust):
  d_repel = his_lys_mindist(HIP) - his_lys_mindist(HIE)   (>0 = His+/Lys+ repulsion at low pH)
Secondary (slower): interface persistence, local (near-switch) and global:
  d_iface_release = iface_frac(HIE) - iface_frac(HIP)       (>0 = holds at 7.4, releases at 6.0)
All 11 switches are the SAME chemistry (His+ vs Lys+) so fixed-charge GBn2 occupancy is a fair
within-set ranking.  TARGET CA restrained (rigid receptor); binder FREE to release.  Only the switch
residue differs between states; all other titratables held at pH-7 default (common mode).

Usage: md_release_tnf.py --pdb complex.pdb --his-pos 19 --out dir [--binder-chain A]
       [--variants HIP,HIE,GLN] [--equil-ps 100 --prod-ps 5000 --interval-ps 10 --replicas 4]
"""
import argparse, json, os, sys
import numpy as np
from openmm import app, unit, Platform, LangevinMiddleIntegrator, CustomExternalForce
from pdbfixer import PDBFixer

IFACE_CONTACT_CUT = 4.5    # A, heavy-atom residue contact defining the initial interface
IFACE_KEEP_CUT    = 5.0    # A, a contact counts as "retained" below this
LOCAL_RADIUS      = 10.0   # A, binder residues within this of the switch CA define the LOCAL interface
IFACE_THRESH      = 0.08   # min d_iface_release to call a (secondary) release signal

def log(*a): print(*a, file=sys.stderr, flush=True)

def pick_platform(device):
    for name in ("OpenCL", "CUDA"):   # OpenCL FIRST: envs/md nvrtc drifted to 13.4 > driver 580 (CUDA13.0) -> CUDA PTX fails
        try: plat = Platform.getPlatformByName(name)
        except Exception: continue
        names = set(plat.getPropertyNames()); props = {}
        if "DeviceIndex" in names: props["DeviceIndex"] = str(device)
        if "Precision" in names: props["Precision"] = "mixed"
        return plat, props, name
    return Platform.getPlatformByName("CPU"), {}, "CPU"

def detect_chains(topology):
    sizes = sorted(((c, sum(1 for _ in c.residues())) for c in topology.chains()), key=lambda x: -x[1])
    return sizes[0][0], (sizes[1][0] if len(sizes) > 1 else None)

def res_by_id(chain, pos):
    c = [r for r in chain.residues() if str(r.id) == str(pos)]
    if not c: raise RuntimeError(f"no residue at binder id {pos}; ids {[r.id for r in chain.residues()][:5]}...")
    return c[0]

def switch_polar_atoms(res):
    """Charged/polar sidechain termini for the switch residue: His ND1/NE2, Gln OE1/NE2, else CB."""
    nm = {a.name: a.index for a in res.atoms()}
    if res.name.startswith("HI"): keys = ("ND1", "NE2")
    elif res.name == "GLN": keys = ("OE1", "NE2")
    else: keys = ("CB",)
    got = [nm[k] for k in keys if k in nm]
    return got or [nm.get("CA")]

def nearest_lys_nz(positions, sw_atoms, target_chain):
    # BURIED-SWITCH variant: find nearest buried CATION (Lys NZ or Arg NH1/NH2/NE) in the SEARCHED chain.
    pos = np.array(positions.value_in_unit(unit.angstrom))
    best, best_d = None, 1e9
    for r in target_chain.residues():
        if r.name == "LYS":   catn = [a.index for a in r.atoms() if a.name == "NZ"]
        elif r.name == "ARG": catn = [a.index for a in r.atoms() if a.name in ("NH1","NH2","NE")]
        else: continue
        if not catn: continue
        d = float(np.min(np.linalg.norm(pos[catn][:, None] - pos[sw_atoms][None, :], axis=-1)))
        if d < best_d: best_d, best = d, (r, catn)
    return best, best_d

def build(pdb, his_variant, his_list, binder_chain="A", mutate_gln=False):
    fixer = PDBFixer(filename=pdb)
    # Identify the binder chain the SAME way the His-Lys distance is measured (detect_chains = size-based:
    # 2nd-largest chain). Mutating by the --binder-chain arg instead would, on any mismatch, no-op the
    # isostere so the GLN "dud" still contains the His and the causal control passes vacuously (d_iso~0).
    _t, _b = detect_chains(fixer.topology)
    binder_cid = _b.id if _b is not None else binder_chain
    if mutate_gln:
        fixer.applyMutations([f"HIS-{h}-GLN" for h in his_list], binder_cid)  # causal dud: neutral isostere
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    target, binder = detect_chains(modeller.topology)
    residues = list(modeller.topology.residues())
    variants = [None] * len(residues)
    if not mutate_gln:
        for h in his_list:
            variants[residues.index(res_by_id(binder, h))] = his_variant
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    target, binder = detect_chains(modeller.topology)
    sws = [res_by_id(binder, h) for h in his_list]
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds, hydrogenMass=4.0 * unit.amu)
    return ff, modeller, system, target, binder, sws

def add_ca_restraints(system, modeller, chain, k_kcal):
    k = k_kcal * unit.kilocalories_per_mole / unit.angstrom**2
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"): force.addPerParticleParameter(p)
    pos = modeller.positions; idx = {a.index for a in chain.atoms()}; n = 0
    for atom in modeller.topology.atoms():
        if atom.name == "CA" and atom.index in idx:
            x, y, z = pos[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, [x, y, z]); n += 1
    system.addForce(force); return n

def heavy_by_res(chain):
    out = {}
    for r in chain.residues():
        hs = [a.index for a in r.atoms() if a.element is not None and a.element.symbol != "H"]
        if hs: out[r] = np.array(hs)
    return out

def initial_interface(pos, binder, target, sw_ca):
    """interface residue pairs (binder,target) in heavy-atom contact at t0; flag 'local' = binder
    residue CA within LOCAL_RADIUS of the switch CA (where the switch acts)."""
    bh, th = heavy_by_res(binder), heavy_by_res(target)
    bca = {r: next((a.index for a in r.atoms() if a.name == "CA"), None) for r in bh}
    tca = {r: next((a.index for a in r.atoms() if a.name == "CA"), None) for r in th}
    pairs, local = [], []
    for rb, ib in bh.items():
        cb = bca[rb]
        if cb is None: continue
        isloc = (sw_ca is not None and np.linalg.norm(pos[cb] - pos[sw_ca]) < LOCAL_RADIUS)
        for rt, it in th.items():
            ct = tca[rt]
            if ct is None or np.linalg.norm(pos[cb] - pos[ct]) > 14.0: continue
            d = np.min(np.linalg.norm(pos[ib][:, None] - pos[it][None, :], axis=-1))
            if d < IFACE_CONTACT_CUT: pairs.append((ib, it)); local.append(isloc)
    return pairs, np.array(local, dtype=bool)

def pair_mindists(pos, pairs):
    return np.array([np.min(np.linalg.norm(pos[ib][:, None] - pos[it][None, :], axis=-1)) for ib, it in pairs])

def run_state(pdb, variant, his_list, args, device, replica):
    mut = (variant == "GLN")
    ff, modeller, system, target, binder, sws = build(pdb, variant, his_list, args.binder_chain, mut)
    pos0 = np.array(modeller.positions.value_in_unit(unit.angstrom))
    # per switch His: its polar atoms + its nearest target Lys NZ (the His-Lys pair that repels)
    groups = []; d0s = []; lys_ids = []
    for sw in sws:
        sa = switch_polar_atoms(sw)
        lys, d0 = nearest_lys_nz(modeller.positions, sa, binder)
        if lys is None: raise RuntimeError("no Lys in target chain")
        groups.append((sa, lys[1])); d0s.append(round(d0, 2)); lys_ids.append(str(lys[0].id))
    sw_ca = next((a.index for a in sws[0].atoms() if a.name == "CA"), None)
    pairs, local = initial_interface(pos0, binder, target, sw_ca)
    n_init = len(pairs); n_loc = int(local.sum())
    binder_ca = [a.index for a in binder.atoms() if a.name == "CA"]
    add_ca_restraints(system, modeller, target, args.restraint_k)   # TARGET rigid, binder FREE
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(7000 + replica)                       # PAIRED across variants within a replica
    plat, props, pname = pick_platform(device)
    sim = app.Simulation(modeller.topology, system, integ, plat, props)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    ref = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))[binder_ca]
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 7000 + replica)
    sim.step(int(args.equil_ps / 0.004))
    nsteps = int(args.prod_ps / 0.004); stride = int(args.interval_ps / 0.004)
    hl, fg, fl, rm = [], [], [], []
    for _ in range(nsteps // stride):
        sim.step(stride)
        p = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))
        # mean over the switch His of each His's min distance to ITS nearest target Lys NZ
        hl.append(float(np.mean([np.min(np.linalg.norm(p[nz][:, None] - p[sa][None, :], axis=-1)) for sa, nz in groups])))
        if n_init:
            d = pair_mindists(p, pairs); keep = d < IFACE_KEEP_CUT
            fg.append(float(np.mean(keep)))
            if n_loc: fl.append(float(np.mean(keep[local])))
        rm.append(float(np.sqrt(np.mean(np.sum((p[binder_ca] - ref) ** 2, axis=1)))))
    return {"variant": variant, "replica": replica, "platform": pname,
            "lys_resid": ",".join(lys_ids), "lys_init_dist": d0s,
            "n_iface_init": n_init, "n_local_init": n_loc,
            "his_lys_mindist": round(float(np.mean(hl)), 2),
            "iface_frac": round(float(np.mean(fg)), 3) if fg else None,
            "local_iface_frac": round(float(np.mean(fl)), 3) if fl else None,
            "binder_rmsd": round(float(np.mean(rm)), 2), "n_frames": len(hl)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True); ap.add_argument("--his-pos", required=True, help="switch His resnum(s), comma-separated for cooperative multi-His")
    ap.add_argument("--out", required=True); ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--variants", default="HIP,HIE,GLN")
    ap.add_argument("--equil-ps", type=float, default=100.0); ap.add_argument("--prod-ps", type=float, default=5000.0)
    ap.add_argument("--interval-ps", type=float, default=10.0); ap.add_argument("--replicas", type=int, default=4)
    ap.add_argument("--restraint-k", type=float, default=2.0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    name = os.path.splitext(os.path.basename(args.pdb))[0]
    his_list = [int(x) for x in str(args.his_pos).split(",")]   # 1+ switch His positions (comma list)
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]
    log("=== MD release gate:", name, "his", his_list, "variants", variants, "dev", args.device)
    res = {"name": name, "pdb": os.path.abspath(args.pdb), "his_pos": his_list,
           "variants": variants, "replicas": args.replicas, "prod_ps": args.prod_ps, "states": []}
    for rep in range(args.replicas):
        for variant in variants:
            try:
                r = run_state(args.pdb, variant, his_list, args, args.device, rep)
                res["states"].append(r)
                log("  rep%d %s hisLys=%s iface=%s local=%s rmsd=%s lys=%s d0=%s" %
                    (rep, variant, r["his_lys_mindist"], r["iface_frac"], r["local_iface_frac"],
                     r["binder_rmsd"], r["lys_resid"], r["lys_init_dist"]))
            except Exception as e:
                res["states"].append({"variant": variant, "replica": rep, "error": repr(e)})
                log("  rep%d %s FAILED: %r" % (rep, variant, e))
    def vals(metric, v):
        return [s[metric] for s in res["states"] if s.get("variant") == v and s.get(metric) is not None]
    def mean(metric, v):
        xs = vals(metric, v); return float(np.mean(xs)) if xs else None
    def sem(metric, v):
        xs = vals(metric, v); return float(np.std(xs) / max(1, len(xs) ** 0.5)) if len(xs) > 1 else None
    agg = {}
    for v in variants:
        agg[v] = {m: mean(m, v) for m in ("his_lys_mindist", "iface_frac", "local_iface_frac", "binder_rmsd")}
        agg[v]["sem_his_lys"] = sem("his_lys_mindist", v)
        agg[v]["sem_iface"] = sem("iface_frac", v)
    res["agg"] = agg
    def d(metric, a, b):
        return (round(agg[a][metric] - agg[b][metric], 3)
                if agg.get(a, {}).get(metric) is not None and agg.get(b, {}).get(metric) is not None else None)
    if "HIP" in agg and "HIE" in agg:
        res["d_repel"] = d("his_lys_mindist", "HIP", "HIE")            # >0 = repulsion at low pH (PRIMARY)
        res["d_iface_release"] = d("iface_frac", "HIE", "HIP")          # >0 = releases at low pH
        res["d_local_release"] = d("local_iface_frac", "HIE", "HIP")
        res["d_rmsd"] = d("binder_rmsd", "HIP", "HIE")
        if "GLN" in agg:   # causal dud: neutral isostere must NOT release vs HIE (both neutral)
            res["dud_repel"] = d("his_lys_mindist", "GLN", "HIE")
            res["dud_iface"] = d("iface_frac", "HIE", "GLN")
        noise = agg["HIE"].get("sem_his_lys")
        res["selective_release"] = bool(res["d_repel"] is not None and res["d_repel"] > 0
                                        and (noise is None or res["d_repel"] > 2 * noise))
    outp = os.path.join(args.out, name + ".json")
    json.dump(res, open(outp, "w"), indent=2)
    log("WROTE %s | d_repel=%s d_iface=%s d_local=%s dud_repel=%s selective=%s" %
        (outp, res.get("d_repel"), res.get("d_iface_release"), res.get("d_local_release"),
         res.get("dud_repel"), res.get("selective_release")))
    print(outp)

if __name__ == "__main__":
    main()
