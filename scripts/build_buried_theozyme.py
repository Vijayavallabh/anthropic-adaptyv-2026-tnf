#!/usr/bin/env python3
"""WAVE 26 MOONSHOT: BURIED His+...Asp460 salt-bridge pH switch -- the one escape from the
strength-vs-desolvation-vs-ECC trade-off. A buried (desolvated) salt bridge is STRONG (full Coulomb -> net
tumor-on, overcomes the His desolvation penalty) AND ECC-robust (a deep well keeps the bridge formed under
0.8 charge-scaling, where surface inv01/inv02 collapsed). Must titrate in-window (not constitutive).

Two theozyme strategies (from the design-exploration workflow), both adapting build_his_asp_theozyme on the
conserved Asp460 (A436) on the proven-bindable His433 face:

  S4 (PRIMARY) groove-cap: the TARGET groove around Asp460 (walls Phe412/Gln408/His409) already supplies the
    pocket (Asp460 relSASA 0.19); the binder caps the mouth. Place the His climbing the groove bore along a
    LOCAL mouth-normal, CB toward the mouth (backbone attaches from the mouth), and seat the ring DEEPEST
    between the conserved walls (maximize target contacts around the ring = burial proxy).

  S1 (HEDGE) deep-pocket engulf: pre-encode concavity. Approach the carboxylate along a surface-TANGENT so the
    His sits BESIDE it against the target wall (CB inward = pocket-wall residue), plus 2 ALA rim anchors on a
    6.5-7.5A shell (CB inward), DELIBERATELY leaving the His-distal azimuth open (asymmetric solvent channel so
    the His stays titratable). Build-time pocket-QC: carboxylate solvent-hemisphere occlusion > 0.5.
"""
import os
import json, sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_his_asp_theozyme import unit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = f"{REPO}/specs/rfd3/theozyme.pdb"
OUTDIR = f"{REPO}/specs/rfd3"
ELEM = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C",
        "ND1": "N", "CD2": "C", "CE1": "C", "NE2": "N"}
ORDER = ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"]
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
HBOND = 3.0   # raised from 2.8 (slightly off-linear) to bound the pKa up-shift (anti-constitutive)
CLASH_MIN = 2.5


def _rot_axis(axis, deg):
    axis = unit(axis); th = np.radians(deg)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def _align(a, b):
    v = np.cross(a, b); s = np.linalg.norm(v); c = np.dot(a, b)
    if s < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def _his_template(chA, template_his_author=394):
    this = chA[template_his_author]
    T = {n: this[n].coord.copy() for n in ORDER}
    ne2 = T["NE2"].copy()
    for n in ORDER:
        T[n] = T[n] - ne2
    ring_cen = np.mean([T[a] for a in HIS_RING], axis=0)
    return T, unit(ring_cen)


def build_his_buried(struct, asp_author, approach, template_his_author=394):
    """Place a His NE2-H-bonding the Asp460 carboxylate with the ring seated DEEP against the target (buried).
    approach='groove' (S4: climb the local mouth-normal, CB toward mouth) or 'tangent' (S1: beside the
    carboxylate along a surface tangent, CB inward). Objective among clash-free rotamers = MAXIMIZE target
    heavy-atom contacts around the His ring (burial proxy)."""
    chA = struct[0]["A"]
    asp = chA[asp_author]
    od1, od2 = asp["OD1"].coord, asp["OD2"].coord
    o_mid = (od1 + od2) / 2.0
    centroid = np.mean([a.coord for r in chA for a in r], axis=0)
    u_out = unit(o_mid - centroid)
    o_acc = od1 if np.dot(od1 - centroid, u_out) >= np.dot(od2 - centroid, u_out) else od2

    if approach == "groove":
        # local groove-mouth normal from the groove-floor centroid (Asp460 mid, Phe412 ring, Gln408 CG) -> solvent
        phe = chA[412]; gln = chA[408]
        phe_cen = np.mean([phe[a].coord for a in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ") if a in phe], axis=0)
        floor = np.mean([o_mid, phe_cen, gln["CG"].coord], axis=0)
        u_app = unit(floor - centroid)        # mouth normal (outward through the groove)
        donor_pos = o_acc + HBOND * u_app
        cb_sign = +1   # CB toward the mouth (backbone attaches from the mouth)
    else:  # tangent (S1)
        # a surface-tangent direction (perp to u_out) so the His sits BESIDE the carboxylate against the wall
        ref = np.array([0, 0, 1.0]) if abs(u_out[2]) < 0.9 else np.array([1.0, 0, 0])
        u_app = unit(np.cross(u_out, ref))
        donor_pos = o_acc + 2.8 * u_app
        cb_sign = -1   # CB inward (pocket-wall residue)

    T, tmpl_dir = _his_template(chA, template_his_author)
    R0 = _align(tmpl_dir, u_app)              # map template NE2->ring onto the approach dir
    antigen = np.array([a.coord for r in chA for a in r if r.id[1] != asp_author])

    best = None; best_any = None
    for deg in range(0, 360, 6):
        R = _rot_axis(u_app, deg) @ R0
        placed = {n: R @ T[n] + donor_pos for n in ORDER}
        sc = np.array([placed[n] for n in ("CB",) + HIS_RING])
        clash = float(np.linalg.norm(antigen[None, :, :] - sc[:, None, :], axis=2).min())
        cb_out = float(np.dot(placed["CB"] - donor_pos, u_out))   # + = toward solvent
        ring = np.array([placed[n] for n in HIS_RING])
        contacts = int((np.linalg.norm(antigen[None, :, :] - ring[:, None, :], axis=2) < 6.0).any(axis=0).sum())
        ne2_o = float(np.linalg.norm(placed["NE2"] - o_acc))
        ok_cb = (cb_out >= 0) if cb_sign > 0 else (cb_out <= 0.5)
        if clash >= CLASH_MIN and ok_cb and ne2_o <= HBOND + 0.6:
            if best is None or contacts > best[0]:
                best = (contacts, placed, clash, cb_out, ne2_o)
        if best_any is None or clash > best_any[2]:
            best_any = (contacts, placed, clash, cb_out, ne2_o)
    if best is None:
        best = best_any
        print(f"  WARN[{approach}]: no rotamer met gates -> least-clashing (clash {best[2]:.2f})")
    contacts, placed, clash, cb_out, ne2_o = best
    diag = dict(approach=approach, asp=asp_author, ne2_O=round(ne2_o, 2), clash=round(clash, 2),
                cb_out=round(cb_out, 2), ring_target_contacts=contacts)
    return placed, o_mid, u_out, diag


def build_rim_anchors(struct, asp_author, o_mid, u_out, his_placed, n_anchors=2):
    """S1: ALA rim anchors on a 6.5-7.5A shell around the carboxylate, solvent hemisphere, CB inward, leaving
    the His-distal azimuth open (asymmetric). Returns list of atom-dicts (N,CA,C,O,CB) + QC occlusion frac."""
    chA = struct[0]["A"]
    antigen = np.array([a.coord for r in chA for a in r if r.id[1] != asp_author])
    his_ring_cen = np.mean([his_placed[n] for n in HIS_RING], axis=0)
    # basis in the solvent hemisphere plane
    e1 = unit(his_ring_cen - o_mid - np.dot(his_ring_cen - o_mid, u_out) * u_out)
    if np.linalg.norm(e1) < 1e-6:
        e1 = unit(np.cross(u_out, [1.0, 0, 0]))
    e2 = unit(np.cross(u_out, e1))
    his_az = np.arctan2(np.dot(his_ring_cen - o_mid, e2), np.dot(his_ring_cen - o_mid, e1))
    anchors = []
    # place anchors ~120-240 deg away from the His azimuth (leave His edge open)
    for k in range(n_anchors):
        az = his_az + np.radians(120 + 120 * k / max(n_anchors, 1))
        for R in (6.5, 7.0, 7.5):
            cb_pos = o_mid + 0.5 * u_out * R + R * (np.cos(az) * e1 + np.sin(az) * e2)
            # CB inward toward o_mid; build a crude ALA backbone outside it
            inward = unit(o_mid - cb_pos)
            ca = cb_pos - 1.5 * inward
            n = ca - 1.46 * unit(np.cross(inward, e2)) - 0.5 * inward
            c = ca + 1.52 * unit(np.cross(inward, e1)) - 0.5 * inward
            o = c + 1.23 * unit(c - ca)
            atoms = {"N": n, "CA": ca, "C": c, "O": o, "CB": cb_pos}
            sc = np.array(list(atoms.values()))
            if float(np.linalg.norm(antigen[None, :, :] - sc[:, None, :], axis=2).min()) >= CLASH_MIN:
                anchors.append(atoms); break
    # QC: occlusion of the carboxylate solvent hemisphere by motif atoms within 8A
    motif = [his_placed[n] for n in ORDER] + [a for anc in anchors for a in anc.values()]
    dirs = []
    for _ in range(200):
        v = np.random.randn(3); v = unit(v)
        if np.dot(v, u_out) > 0:
            dirs.append(v)
    occ = 0
    for v in dirs:
        p = o_mid + 4.0 * v
        if any(np.linalg.norm(p - m) < 4.0 for m in motif):
            occ += 1
    return anchors, round(occ / max(len(dirs), 1), 2)


def write_chainC(struct, residues, out_pdb):
    """residues = list of (resid, resname, atomdict)."""
    model = struct[0]
    if "C" in [c.id for c in model]:
        model.detach_child("C")
    ch = Chain("C")
    for resid, resname, atoms in residues:
        res = Residue((" ", resid, " "), resname, "")
        names = ORDER if resname == "HIS" else ["N", "CA", "C", "O", "CB"]
        for i, n in enumerate(names):
            res.add(Atom(n, atoms[n].astype(float), 0.0, 1.0, " ", (" " + n).ljust(4), i + 1,
                         ELEM.get(n, n[0])))
        ch.add(res)
    model.add(ch)
    io = PDBIO(); io.set_structure(struct); io.save(out_pdb)


def main():
    ASP = 436  # Asp460 (conserved human=mouse) on the His433 face
    # ---- S4 groove-cap (PRIMARY) ----
    struct = PDBParser(QUIET=True).get_structure("t", SRC)
    placed, o_mid, u_out, diag = build_his_buried(struct, ASP, "groove")
    s4_pdb = f"{OUTDIR}/theozyme_buried_s4.pdb"; s4_json = f"{OUTDIR}/egfr_buried_s4.json"
    write_chainC(struct, [(1, "HIS", placed)], s4_pdb)
    json.dump({"egfr_buried_s4": {
        "dialect": 2, "input": s4_pdb, "contig": "40-55,C1-1,45-60,/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL", "C1": "CB,CG,ND1,CD2,CE1,NE2"},
        "select_hbond_donor": {"C1": "ND1,NE2"}, "select_hbond_acceptor": {f"A{ASP}": "OD1,OD2"},
        "select_hotspots": {f"A{ASP}": "OD1,OD2", "A412": "CG,CE1,CZ", "A409": "CB,CE1,NE2", "A408": "CB,CG"},
        "infer_ori_strategy": "hotspots", "is_non_loopy": True}}, open(s4_json, "w"), indent=2)
    print(f"S4 groove: NE2..O={diag['ne2_O']} clash={diag['clash']} cb_out={diag['cb_out']} "
          f"ring_target_contacts={diag['ring_target_contacts']} -> {s4_pdb.split('/')[-1]}, {s4_json.split('/')[-1]}")

    # ---- S1 deep-pocket engulf (HEDGE, 2 rim anchors) ----
    struct2 = PDBParser(QUIET=True).get_structure("t", SRC)
    placed2, o_mid2, u_out2, diag2 = build_his_buried(struct2, ASP, "tangent")
    anchors, occ = build_rim_anchors(struct2, ASP, o_mid2, u_out2, placed2, n_anchors=2)
    residues = [(1, "HIS", placed2)] + [(2 + i, "ALA", a) for i, a in enumerate(anchors)]
    s1_pdb = f"{OUTDIR}/theozyme_buried_s1.pdb"; s1_json = f"{OUTDIR}/egfr_buried_s1.json"
    write_chainC(struct2, residues, s1_pdb)
    fixed = {"A310-500": "ALL", "C1": "CB,CG,ND1,CD2,CE1,NE2"}
    for i in range(len(anchors)):
        fixed[f"C{2+i}"] = "CA,CB"
    anchor_contig = "".join(f"8-16,C{2+i}-{2+i}," for i in range(len(anchors)))
    json.dump({"egfr_buried_s1": {
        "dialect": 2, "input": s1_pdb,
        "contig": f"18-40,C1-1,{anchor_contig}18-40,/0,A310-500",
        "select_fixed_atoms": fixed, "select_hbond_donor": {"C1": "ND1,NE2"},
        "select_hbond_acceptor": {f"A{ASP}": "OD1,OD2"},
        "select_hotspots": {f"A{ASP}": "OD1,OD2", "A409": "ND1,NE2", "A412": "CG,CD1,CD2,CE1,CE2,CZ"},
        "infer_ori_strategy": "hotspots", "is_non_loopy": True}}, open(s1_json, "w"), indent=2)
    print(f"S1 pocket: NE2..O={diag2['ne2_O']} clash={diag2['clash']} cb_out={diag2['cb_out']} "
          f"anchors={len(anchors)} occlusion={occ} -> {s1_pdb.split('/')[-1]}, {s1_json.split('/')[-1]}")
    if occ < 0.5:
        print(f"  WARN S1: occlusion {occ} < 0.5 (pocket may be shallow)")


if __name__ == "__main__":
    main()
