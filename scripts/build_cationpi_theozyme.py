#!/usr/bin/env python3
"""WAVE 24 MOONSHOT: a CATION-pi pH switch -- an ECC-robust alternative to the His+...Asp- salt bridge.

The stricter robustness test (electronic-continuum-correction, charge-scaling d_occ) showed our inverted
His+...Asp- salt-bridge switch COLLAPSES under 0.8 charge scaling: it leans on the ~3-4 kcal/mol
over-stabilization that fixed-charge force fields carry for ionic pairs. Doubling that chemistry (the
cooperative two-His...two-Asp design) inherits the same fragility. So we change the PHYSICS of the switch.

A cation-pi contact (protonated His imidazolium stacked on an aromatic ring) is strengthened by protonation
(His+...Trp -5.4 vs His0...Trp -4.0 kcal/mol aqueous; JCTC 2024 PMC11325542) yet the aromatic partner is
NEUTRAL -- there is NO counter-charge for the fixed-charge FF to over-stabilize, and a neutral aromatic does
NOT raise the His pKa (so it structurally avoids the constitutive-bridge / too-high-pKa trap). It is therefore
robust to the exact ECC scaling that killed the salt bridge, and it is a switch chemistry never used in de
novo pH-conditional binder design.

Target partner: EGFR Domain III PHE412 (mature numbering) = Phe436 (P00533 precursor) -- a CONSERVED
(human==mouse) epitope anchor, ring SASA ~40 A^2 (face accessible), only 6.8 A from the sensor His433 on the
PROVEN-BINDABLE His433 face (the one epitope af2ig confirms binders dock -- the lesson that killed the pure
conserved-Asp inverted de novo). So: cross-reactive by construction (conserved Phe), ECC-robust (neutral
aromatic), and on the one epitope that actually binds.

Output: an RFd3 theozyme (chain A = DIII target 310-500; chain C res 1 = a binder His whose imidazole ring is
placed PARALLEL-stacked ~3.7 A over the Phe412 ring face, on the solvent side) + the cation-pi RFd3 spec
(fix the His ring atoms; hotspot Phe412 ring + the His433 face A409; NO h-bond donor/acceptor -- this is a
pi contact, not a salt bridge). Reuses the clash-gated rotamer search pattern from build_his_asp_theozyme.
"""
import os
import json, sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = f"{REPO}/specs/rfd3/theozyme.pdb"   # chain A = DIII target (mature numbering 310-500)
OUTDIR = f"{REPO}/specs/rfd3"

STACK_D = 3.7     # A, His ring centroid above the Phe ring plane (cation-pi stacking distance)
LAT_OFF = 1.2     # A, lateral offset -> parallel-DISPLACED stack (energetically favoured over sandwich)
ELEM = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C",
        "ND1": "N", "CD2": "C", "CE1": "C", "NE2": "N"}
ORDER = ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"]
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
PHE_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def ring_frame(coords):
    """centroid + unit normal (smallest-variance axis via SVD) of a set of ring-atom coords."""
    c = coords.mean(0)
    _, _, vt = np.linalg.svd(coords - c)
    return c, unit(vt[2])


def _axis_rot(axis, deg):
    axis = unit(axis); th = np.radians(deg)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def build_his_cationpi(struct, phe_author, template_his_author=394):
    """Place a HIS (chain C res 1) so its imidazole ring engages the Phe ring face at cation-pi distance
    FROM THE SOLVENT SIDE, with CB pointing outward so a binder backbone can hold it without clashing. A
    perfectly flat parallel stack forces the His backbone in-plane (into the antigen), so we search over
    TILT (ring normal tilted from the Phe normal toward solvent) x in-plane OFFSET x stacking DISTANCE x
    SPIN, require clash-free to the antigen (Phe412 excluded -- it is the intended contact), and maximize
    CB-outward while keeping the ring over the Phe face at ~3.4-4.8 A. Returns atomname->coord + diag."""
    chA = struct[0]["A"]
    phe = chA[phe_author]
    assert phe.resname in ("PHE", "TYR", "TRP"), f"A{phe_author} is {phe.resname}, not aromatic"
    phe_ring = np.array([phe[a].coord for a in PHE_RING if a in phe])
    c_phe, n_phe = ring_frame(phe_ring)
    centroid = np.mean([a.coord for r in chA for a in r], axis=0)
    if np.dot(n_phe, c_phe - centroid) < 0:          # Phe normal -> outward (solvent)
        n_phe = -n_phe
    # Clash model: the target BACKBONE (+CB) is rigid, but flexible surface-polar sidechains reposition on
    # binding, so their atoms beyond CB are DISPLACEABLE and excluded from the clash set (the Phe partner is
    # excluded too -- it is the intended contact). This is why a static crystal Gln411 rotamer should not
    # veto a stack on the Phe412 face. af2ig (fully flexible) is the real binding test downstream.
    FLEX = {"ASN", "GLN", "GLU", "ASP", "LYS", "ARG", "SER", "THR", "CYS", "MET"}
    BB_CB = {"N", "CA", "C", "O", "CB"}
    antigen = np.array([a.coord for r in chA for a in r
                        if r.id[1] != phe_author and not (r.resname in FLEX and a.name not in BB_CB)])
    # local solvent-escape direction: from the nearby antigen bulk out to the Phe ring
    near = antigen[np.linalg.norm(antigen - c_phe, axis=1) < 10.0]
    e = unit(c_phe - (near.mean(0) if len(near) else centroid))
    e_perp = e - np.dot(e, n_phe) * n_phe             # in-plane escape component
    e_perp = unit(e_perp) if np.linalg.norm(e_perp) > 1e-6 else unit(np.cross(n_phe, [1.0, 0, 0]))
    tilt_axis = unit(np.cross(n_phe, e))              # tilt the ring toward solvent about this

    this = chA[template_his_author]                   # real His internal geometry
    T0 = {n: this[n].coord.copy() for n in ORDER}
    tmpl_cen, tmpl_norm = ring_frame(np.array([T0[a] for a in HIS_RING]))
    T0 = {n: T0[n] - tmpl_cen for n in ORDER}         # ring centroid -> origin

    def align_R(a_vec, b_vec):                        # rotation mapping a_vec onto b_vec
        v = np.cross(a_vec, b_vec); s = np.linalg.norm(v); c = np.dot(a_vec, b_vec)
        if s < 1e-8:
            return np.eye(3) if c > 0 else _axis_rot(unit(np.cross(a_vec, [1, 0, 0]) + 1e-6), 180)
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))

    CLASH_MIN = 2.6
    best_ok, best_any = None, None
    # cation-pi tolerates an OFFSET/edge stack over the exposed part of a semi-buried ring; search wide and
    # score clash-free poses by (CB-outward) while keeping the His+ ring near the Phe pi face.
    for tilt in (0, 20, 35, 50, 65):                  # ring normal tilted from Phe normal toward solvent
        a_norm = _axis_rot(tilt_axis, tilt) @ n_phe
        R0 = align_R(tmpl_norm, a_norm)               # template ring normal -> a_norm
        for D in (3.6, 3.9, 4.2, 4.6, 5.0):
            for off in (0.0, 1.0, 2.0, 3.0):
                center = c_phe + D * n_phe + off * e_perp
                cc = float(np.linalg.norm(center - c_phe))
                off_face = float(np.linalg.norm((center - c_phe) - np.dot(center - c_phe, n_phe) * n_phe))
                for deg in range(0, 360, 6):
                    R = _axis_rot(a_norm, deg) @ R0
                    placed = {n: R @ T0[n] + center for n in ORDER}
                    sc = np.array([placed[n] for n in ("CB",) + HIS_RING])
                    clash = float(np.linalg.norm(antigen[None, :, :] - sc[:, None, :], axis=2).min())
                    cb_out = float(np.dot(unit(placed["CB"] - center), e))   # CB toward solvent
                    # the His+ cation sits over/near the pi face: require ring-centroid within cation-pi
                    # range AND not drifted too far off the face edge
                    if clash >= CLASH_MIN and 3.4 <= cc <= 5.2 and off_face <= 3.5:
                        if best_ok is None or cb_out > best_ok[0]:
                            best_ok = (cb_out, placed, clash, cb_out)
                    if best_any is None or clash > best_any[2]:
                        best_any = (cb_out, placed, clash, cb_out)
    if best_ok is not None:
        _, placed, clash, cb_out = best_ok
    else:
        _, placed, clash, cb_out = best_any
        print("  WARN: no clash-free cation-pi pose (best clash %.2f < %.1f) -> least-clashing" % (clash, CLASH_MIN))

    his_cen, his_norm = ring_frame(np.array([placed[a] for a in HIS_RING]))
    cc = float(np.linalg.norm(his_cen - c_phe))
    ang = float(np.degrees(np.arccos(min(1.0, abs(np.dot(his_norm, n_phe))))))
    # lateral offset of the His centroid off the Phe-normal axis (how far off the ring face)
    off_axis = float(np.linalg.norm((his_cen - c_phe) - np.dot(his_cen - c_phe, n_phe) * n_phe))
    diag = dict(phe=phe_author, ring_cc=round(cc, 2), plane_angle=round(ang, 1),
                off_face=round(off_axis, 2), clash=round(float(clash), 2), cb_outward=round(float(cb_out), 2))
    return placed, diag


def write_theozyme(struct, placed, out_pdb):
    model = struct[0]
    if "C" in [c.id for c in model]:
        model.detach_child("C")
    ch = Chain("C")
    res = Residue((" ", 1, " "), "HIS", "")
    for i, n in enumerate(ORDER):
        res.add(Atom(n, placed[n].astype(float), 0.0, 1.0, " ", (" " + n).ljust(4), i + 1, ELEM[n]))
    ch.add(res)
    model.add(ch)
    io = PDBIO(); io.set_structure(struct); io.save(out_pdb)


def write_spec(name, out_pdb, phe_author, out_json, binder="25-50"):
    """cation-pi spec: fix the His ring; hotspot the Phe ring + the His433 face; NO h-bond terms."""
    spec = {name: {
        "dialect": 2,
        "input": out_pdb,
        "contig": f"{binder},C1-1,{binder},/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL", "C1": "CB,CG,ND1,CD2,CE1,NE2"},
        # cation-pi is a pi contact, NOT a salt bridge or h-bond -> no donor/acceptor constraint.
        "select_hotspots": {f"A{phe_author}": "CG,CD1,CD2,CE1,CE2,CZ", "A409": "ND1,NE2"},
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
    }}
    json.dump(spec, open(out_json, "w"), indent=2)


def write_theozyme_dual(struct, placed_a, placed_b, out_pdb):
    """chain C: res 1 = His for Phe_a, res 2 = His for Phe_b (both cation-pi stacked)."""
    model = struct[0]
    if "C" in [c.id for c in model]:
        model.detach_child("C")
    ch = Chain("C")
    for resid, placed in ((1, placed_a), (2, placed_b)):
        res = Residue((" ", resid, " "), "HIS", "")
        for i, n in enumerate(ORDER):
            res.add(Atom(n, placed[n].astype(float), 0.0, 1.0, " ", (" " + n).ljust(4), i + 1, ELEM[n]))
        ch.add(res)
    model.add(ch)
    io = PDBIO(); io.set_structure(struct); io.save(out_pdb)


def write_spec_dual(name, out_pdb, phe_a, phe_b, out_json):
    """DUAL cation-pi: binder presents two His, one stacked on each conserved Phe; hotspot both Phe rings
    + the His433 face. Pure pi contacts -> no donor/acceptor. Larger binder spans both sites (~9 A apart)."""
    spec = {name: {
        "dialect": 2,
        "input": out_pdb,
        "contig": "25-45,C1-1,20-40,C2-2,25-45,/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL",
                               "C1": "CB,CG,ND1,CD2,CE1,NE2", "C2": "CB,CG,ND1,CD2,CE1,NE2"},
        "select_hotspots": {f"A{phe_a}": "CG,CD1,CD2,CE1,CE2,CZ",
                            f"A{phe_b}": "CG,CD1,CD2,CE1,CE2,CZ", "A409": "ND1,NE2"},
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
    }}
    json.dump(spec, open(out_json, "w"), indent=2)


def _report(tag, diag):
    print(f"  {tag}: ring_cc={diag['ring_cc']}A plane_angle={diag['plane_angle']}deg off_face={diag['off_face']}A "
          f"clash={diag['clash']}A CB_out={diag['cb_outward']}")
    warn = []
    if not (3.4 <= diag["ring_cc"] <= 4.8): warn.append(f"ring_cc {diag['ring_cc']} off 3.4-4.8")
    if diag["off_face"] > 2.8: warn.append(f"off_face {diag['off_face']}>2.8")
    if diag["clash"] < 2.2: warn.append(f"clash {diag['clash']}<2.2")
    if diag["cb_outward"] < 0.2: warn.append(f"CB_out {diag['cb_outward']}<0.2")
    print("     WARNING: " + "; ".join(warn) if warn else "     geometry OK")
    return not warn


def main():
    struct = PDBParser(QUIET=True).get_structure("t", SRC)
    PHE_A, PHE_B = 412, 380   # mature: Phe436 + Phe404 (precursor); BOTH conserved human==mouse; 9.1A apart; His433 face

    # ---- DUAL cation-pi (primary: robust AND strong enough to beat the single-His selectivity ceiling) ----
    placed_a, diag_a = build_his_cationpi(struct, PHE_A)
    placed_b, diag_b = build_his_cationpi(struct, PHE_B)
    hca, hna = ring_frame(np.array([placed_a[a] for a in HIS_RING]))
    hcb, hnb = ring_frame(np.array([placed_b[a] for a in HIS_RING]))
    his_his = float(np.linalg.norm(hca - hcb))
    dpdb = f"{OUTDIR}/theozyme_cationpi_dual.pdb"; djson = f"{OUTDIR}/egfr_cationpi_dual.json"
    write_theozyme_dual(struct, placed_a, placed_b, dpdb)
    write_spec_dual("egfr_cationpi_dual", dpdb, PHE_A, PHE_B, djson)
    print("=== DUAL cation-pi (Phe412 + Phe380, both conserved) ===")
    ok_a = _report("His_a->Phe412", diag_a)
    ok_b = _report("His_b->Phe380", diag_b)
    print(f"  His_a .. His_b ring-centroid dist = {his_his:.1f}A  (want >7A so the two His+ stay independent)")
    print(f"  -> {dpdb}\n  -> {djson}")
    if his_his < 7.0:
        print("  WARNING: two His closer than 7A -> like-charge coupling risk")

    # ---- SINGLE cation-pi on Phe412 (fallback / ablation control) ----
    spdb = f"{OUTDIR}/theozyme_cationpi_phe412.pdb"; sjson = f"{OUTDIR}/egfr_cationpi_phe412.json"
    write_theozyme(struct, placed_a, spdb)
    write_spec("egfr_cationpi_phe412", spdb, PHE_A, sjson)
    print("=== SINGLE cation-pi (Phe412; fallback/ablation) ===")
    _report("His->Phe412", diag_a)
    print(f"  -> {spdb}\n  -> {sjson}")


if __name__ == "__main__":
    main()
