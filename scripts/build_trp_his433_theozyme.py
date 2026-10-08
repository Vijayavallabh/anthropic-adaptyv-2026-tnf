#!/usr/bin/env python3
"""WAVE 27 MOONSHOT: binder-Trp ... EGFR-His433 CATION-pi switch -- the Trp-strength, ECC-robust, clean
version of the PROVEN His433 pH switch.

Rationale (derived this session):
 - Buried His+...Asp- salt bridge = NO-GO: fair GFN2-xTB shows the neutral-His OFF-state H-bonds the
   desolvated carboxylate at -11.6..-12.0 kcal/mol (n=4) -> CONSTITUTIVE (confirms MD occ_HIE 0.78).
 - cation-pi OFF-state is only -2.8 kcal/mol (neutral His has nothing to grip on a pure aromatic) -> CLEAN
   + ECC-robust (no counter-charge for fixed-charge FF / ECC to over-stabilize).
 - Trp is the STRONGEST cation-pi partner (biggest quadrupole; model-cluster on-state -6.6 vs Phe -5.0).
   EGFR has NO conserved Trp to target, BUT His433 (the conserved, proven, pKa~6.3 sensor on the ONLY
   proven-bindable epitope) can be the CATION, and we place the Trp on OUR binder -> free choice of the
   strongest partner. Tumor pH -> His433+ ... Trp cation-pi -> binds; blood pH -> His433 neutral -> weak -> off.

Construction: parallel-stack a binder Trp (chain C res 1) ~3.7 A over the His433 (A409 mature) imidazole
face from the solvent side, CB pointing outward (so a binder backbone can hold it clash-free), searching
TILT x DISTANCE x in-plane OFFSET x SPIN and maximizing CB-outward among antigen-clash-free poses. Trp
template geometry is lifted from a real indole rotamer (A386). Output: theozyme PDB + RFd3 spec that FIXES
the Trp ring atoms and HOTSPOTS the His433 ring (A409), with NO h-bond term (pi contact, not a salt bridge).
Adapts scripts/build_cationpi_theozyme.py (roles swapped: aromatic on binder, cation on target)."""
import os
import json, sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = f"{REPO}/specs/rfd3/theozyme.pdb"   # chain A = DIII target (mature numbering 310-500); A409 = His433
OUTDIR = f"{REPO}/specs/rfd3"

STACK_D = 3.7     # A, Trp indole centroid above the His433 ring plane (cation-pi stacking distance)
HIS433 = 409      # mature numbering of the conserved sensor His433
TRP_TEMPLATE = 386  # a real Trp rotamer in chain A -> lift its internal indole geometry

TRP_ORDER = ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
TRP_ELEM = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C", "CD1": "C", "CD2": "C",
            "NE1": "N", "CE2": "C", "CE3": "C", "CZ2": "C", "CZ3": "C", "CH2": "C"}
INDOLE = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")  # 9 ring atoms
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def ring_frame(coords):
    c = coords.mean(0)
    _, _, vt = np.linalg.svd(coords - c)
    return c, unit(vt[2])


def _axis_rot(axis, deg):
    axis = unit(axis); th = np.radians(deg)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def align_R(a_vec, b_vec):
    v = np.cross(a_vec, b_vec); s = np.linalg.norm(v); c = np.dot(a_vec, b_vec)
    if s < 1e-8:
        return np.eye(3) if c > 0 else _axis_rot(unit(np.cross(a_vec, [1, 0, 0]) + 1e-6), 180)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def build_trp_cationpi(struct, his_author=HIS433, template_author=TRP_TEMPLATE):
    """Place a binder TRP (chain C res 1) with its indole ring parallel-stacked over the His433 imidazole
    face from the solvent side; search poses, keep antigen-clash-free, maximize CB-outward."""
    chA = struct[0]["A"]
    his = chA[his_author]
    assert his.resname == "HIS", f"A{his_author} is {his.resname}, not HIS (expected His433)"
    his_ring = np.array([his[a].coord for a in HIS_RING if a in his])
    c_his, n_his = ring_frame(his_ring)
    centroid = np.mean([a.coord for r in chA for a in r], axis=0)
    if np.dot(n_his, c_his - centroid) < 0:              # His433 normal -> outward (solvent)
        n_his = -n_his
    # displaceable flexible sidechains (reposition on binding) excluded from clash set; His433 = intended contact
    FLEX = {"ASN", "GLN", "GLU", "ASP", "LYS", "ARG", "SER", "THR", "CYS", "MET"}
    BB_CB = {"N", "CA", "C", "O", "CB"}
    antigen = np.array([a.coord for r in chA for a in r
                        if r.id[1] != his_author and not (r.resname in FLEX and a.name not in BB_CB)])
    near = antigen[np.linalg.norm(antigen - c_his, axis=1) < 10.0]
    e = unit(c_his - (near.mean(0) if len(near) else centroid))     # local solvent-escape direction
    e_perp = e - np.dot(e, n_his) * n_his
    e_perp = unit(e_perp) if np.linalg.norm(e_perp) > 1e-6 else unit(np.cross(n_his, [1.0, 0, 0]))
    tilt_axis = unit(np.cross(n_his, e))

    tmpl = chA[template_author]
    assert tmpl.resname == "TRP", f"A{template_author} is {tmpl.resname}, not TRP"
    T0 = {n: tmpl[n].coord.copy() for n in TRP_ORDER}
    tmpl_cen, tmpl_norm = ring_frame(np.array([T0[a] for a in INDOLE]))
    T0 = {n: T0[n] - tmpl_cen for n in TRP_ORDER}         # indole centroid -> origin

    CLASH_MIN = 2.6
    best_ok, best_any = None, None
    for tilt in (0, 20, 35, 50, 65):
        a_norm = _axis_rot(tilt_axis, tilt) @ n_his
        R0 = align_R(tmpl_norm, a_norm)
        for D in (3.5, 3.8, 4.1, 4.5, 4.9):
            for off in (0.0, 1.0, 2.0, 3.0):
                center = c_his + D * n_his + off * e_perp
                cc = float(np.linalg.norm(center - c_his))
                for deg in range(0, 360, 6):
                    R = _axis_rot(a_norm, deg) @ R0
                    placed = {n: R @ T0[n] + center for n in TRP_ORDER}
                    sc = np.array([placed[n] for n in ("CB",) + INDOLE])
                    clash = float(np.linalg.norm(antigen[None, :, :] - sc[:, None, :], axis=2).min())
                    cb_out = float(np.dot(unit(placed["CB"] - center), e))
                    off_face = float(np.linalg.norm((center - c_his) - np.dot(center - c_his, n_his) * n_his))
                    if clash >= CLASH_MIN and 3.3 <= cc <= 5.2 and off_face <= 3.5:
                        if best_ok is None or cb_out > best_ok[0]:
                            best_ok = (cb_out, placed, clash, cb_out)
                    if best_any is None or clash > best_any[2]:
                        best_any = (cb_out, placed, clash, cb_out)
    if best_ok is not None:
        _, placed, clash, cb_out = best_ok
    else:
        _, placed, clash, cb_out = best_any
        print("  WARN: no clash-free cation-pi pose (best clash %.2f) -> least-clashing" % clash)

    ind_cen, ind_norm = ring_frame(np.array([placed[a] for a in INDOLE]))
    cc = float(np.linalg.norm(ind_cen - c_his))
    ang = float(np.degrees(np.arccos(min(1.0, abs(np.dot(ind_norm, n_his))))))
    off_axis = float(np.linalg.norm((ind_cen - c_his) - np.dot(ind_cen - c_his, n_his) * n_his))
    diag = dict(his433=his_author, ring_cc=round(cc, 2), plane_angle=round(ang, 1),
                off_face=round(off_axis, 2), clash=round(float(clash), 2), cb_outward=round(float(cb_out), 2))
    return placed, diag


def write_theozyme(struct, placed, out_pdb):
    model = struct[0]
    if "C" in [c.id for c in model]:
        model.detach_child("C")
    ch = Chain("C")
    res = Residue((" ", 1, " "), "TRP", "")
    for i, n in enumerate(TRP_ORDER):
        res.add(Atom(n, placed[n].astype(float), 0.0, 1.0, " ", (" " + n).ljust(4), i + 1, TRP_ELEM[n]))
    ch.add(res)
    model.add(ch)
    io = PDBIO(); io.set_structure(struct); io.save(out_pdb)


def write_spec(name, out_pdb, out_json, binder="25-50"):
    """cation-pi spec: FIX the binder-Trp indole ring; HOTSPOT the His433 (A409) ring + the conserved
    neighbors Q408/F412; NO h-bond term (pi contact). infer_ori from hotspots."""
    spec = {name: {
        "dialect": 2,
        "input": out_pdb,
        "contig": f"{binder},C1-1,{binder},/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL",
                               "C1": "CB,CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2"},
        "select_hotspots": {"A409": "CG,ND1,CD2,CE1,NE2", "A408": "CB,CG", "A412": "CG,CD1,CZ"},
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
    }}
    json.dump(spec, open(out_json, "w"), indent=2)


def main():
    struct = PDBParser(QUIET=True).get_structure("x", SRC)
    placed, diag = build_trp_cationpi(struct)
    print("Trp...His433 cation-pi theozyme:", json.dumps(diag))
    out_pdb = f"{OUTDIR}/theozyme_trp_his433.pdb"
    write_theozyme(struct, placed, out_pdb)
    out_json = f"{OUTDIR}/spec_trp_his433_cationpi.json"
    write_spec("egfr_trp_his433", out_pdb, out_json)
    print("wrote", out_pdb, "and", out_json)


if __name__ == "__main__":
    main()
