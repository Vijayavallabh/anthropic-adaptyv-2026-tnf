#!/usr/bin/env python3
"""WAVE 26: measure how BURIED a His+...carboxylate salt bridge is in the bound complex -- the filter that
separates a buried (strong + ECC-robust) switch from an ECC-fragile surface one.

For the bridge atoms (binder-His ND1/NE2 + partner-acid OD1/OD2/OE1/OE2) we compute solvent-accessible surface
area (Shrake-Rupley) in the FULL complex and in the isolated binder+target; burial = 1 - SASA(complex)/SASA(free
reference of the same atoms in isolation). A BURIED bridge has low complex-SASA on the bridge atoms AND high
burial on complex formation. Reports per-atom and summed SASA + burial fraction. Threshold suggestion: bridge
complex-SASA < ~5 A^2 (near-fully desolvated) = buried; > ~20 A^2 = surface (fragile). Works on af2pred-style
complexes (binder chain A, target chain B); auto-finds the His's nearest acid if --acid not given.
"""
import argparse, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, Select
from Bio.PDB.SASA import ShrakeRupley

HISN = ("ND1", "NE2")
ACIDO = ("OD1", "OD2", "OE1", "OE2")


def sasa_of_atoms(structure, atom_keys):
    """atom_keys = set of (chain_id, resid, atomname). return {key: sasa} after SR compute on `structure`."""
    ShrakeRupley().compute(structure, level="A")
    out = {}
    model = structure[0]
    for ch in model:
        for r in ch:
            for a in r:
                k = (ch.id, r.id[1], a.name)
                if k in atom_keys:
                    out[k] = float(a.sasa)
    return out


class _Keep(Select):
    def __init__(self, chains): self.chains = chains
    def accept_chain(self, ch): return ch.id in self.chains


def run(pdb, binder_chain, target_chain, his_resid, acid_resid, his_atoms=HISN):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("x", pdb)[0]
    binder = model[binder_chain]; target = model[target_chain]
    his = [r for r in binder if r.id[1] == his_resid and r.resname.startswith("HI")]
    if not his:
        raise RuntimeError(f"His {his_resid} not in binder chain {binder_chain}")
    his = his[0]
    hn = {a.name: a.coord for a in his if a.name in his_atoms}
    # find bridging acid if not given: nearest target ASP/GLU carboxylate O to the His N's
    if acid_resid is None:
        best = None
        for r in target:
            if r.resname in ("ASP", "GLU"):
                for a in r:
                    if a.name in ACIDO:
                        d = min(np.linalg.norm(a.coord - c) for c in hn.values())
                        if best is None or d < best[0]:
                            best = (d, r.id[1])
        if best is None:
            raise RuntimeError("no target ASP/GLU carboxylate found")
        acid_resid = best[1]
    acid = [r for r in target if r.id[1] == acid_resid][0]
    his_keys = {(binder_chain, his_resid, n) for n in his_atoms}
    acid_keys = {(target_chain, acid_resid, n) for n in ACIDO if n in [a.name for a in acid]}
    bridge_keys = his_keys | acid_keys
    # min bridge distance (is it actually formed?)
    ao = {a.name: a.coord for a in acid if a.name in ACIDO}
    dmin = min(np.linalg.norm(hn[n] - ao[o]) for n in hn for o in ao)

    # SASA in the full complex
    sasa_cplx = sasa_of_atoms(parser.get_structure("c", pdb), bridge_keys)
    # SASA of the isolated binder (His exposure when target removed) + isolated target (acid exposure)
    import tempfile, os
    from Bio.PDB import PDBIO
    with tempfile.TemporaryDirectory(dir="/home/prachh/v/anthropic-adaptyv-2026-egfr/runs") as wd:
        io = PDBIO(); io.set_structure(model)
        bp = os.path.join(wd, "b.pdb"); io.save(bp, _Keep([binder_chain]))
        tp = os.path.join(wd, "t.pdb"); io.save(tp, _Keep([target_chain]))
        sasa_b = sasa_of_atoms(parser.get_structure("b", bp), his_keys)
        sasa_t = sasa_of_atoms(parser.get_structure("t", tp), acid_keys)
    sasa_free = {**sasa_b, **sasa_t}
    cplx_sum = round(sum(sasa_cplx.get(k, 0.0) for k in bridge_keys), 1)
    free_sum = round(sum(sasa_free.get(k, 0.0) for k in bridge_keys), 1)
    burial = round(1.0 - cplx_sum / free_sum, 3) if free_sum > 1e-6 else None
    return dict(his=his_resid, acid=acid_resid, acid_resname=acid.resname, bridge_dist=round(float(dmin), 2),
                bridge_sasa_complex=cplx_sum, bridge_sasa_free=free_sum, burial_fraction=burial,
                buried=(cplx_sum < 5.0), surface=(cplx_sum > 20.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--binder-chain", default="A"); ap.add_argument("--target-chain", default="B")
    ap.add_argument("--his", type=int, required=True)
    ap.add_argument("--acid", type=int, default=None, help="acid resid in target (auto = nearest)")
    ap.add_argument("--mode", default="full", choices=["full", "interion"],
                    help="full = both His N's; interion = NE2+carboxylate only (asymmetric burial, allows distal ND1 exposed)")
    ap.add_argument("--name", default=None)
    a = ap.parse_args()
    his_atoms = ("NE2",) if a.mode == "interion" else HISN
    res = run(a.pdb, a.binder_chain, a.target_chain, a.his, a.acid, his_atoms)
    res["name"] = a.name
    print(json.dumps(res))


if __name__ == "__main__":
    main()
