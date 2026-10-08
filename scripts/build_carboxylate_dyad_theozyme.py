#!/usr/bin/env python3
"""WAVE 28 MOONSHOT: BURIED carboxylate-DYAD pH switch -- binder-Glu + conserved target-Asp460, the
REPULSION-GATED escape from the stickiness that killed the buried salt bridge.

Premise CONFIRMED (fair GFN2-xTB, scripts/qm_carboxylate_switch_premise.py): a buried Glu/Asp pair gives an
ON state (tumor: protonated GluH...Asp- H-bond, -6.6 kcal/mol @ chcl3, -16.7 opt) and an OFF state (blood:
Glu-...Asp- like-charge REPULSION, +31.8 @ chcl3) -- the off-state is strongly REPULSIVE and burial AMPLIFIES
it (+9.7 surface -> +31.8 buried), the exact inverse of the buried salt bridge's sticky neutral-His off-state
(-11.9). dE_switch -38 kcal/mol. The adjacent Asp460 up-shifts the binder-Glu pKa into the tumor/blood window
(carboxylate dyad = shared proton); burial pushes the binder AWAY at blood pH instead of holding it.

Reuses the WAVE 26 S4 groove-cap (proven-bindable Asp460/A436 groove, walls Phe412/Gln408/His409): place a
binder GLU climbing the groove mouth-normal, carboxyl H-bonding Asp460 at DYAD distance (~2.6A O...O), CB
toward the mouth, seated DEEP (maximize target contacts around the carboxyl = burial proxy). RFd3 spec fixes
the Glu carboxyl, donor=Glu carboxyl, acceptor=Asp460 carboxyl. Adapts build_buried_theozyme.build_his_buried
(binder residue His->Glu, NE2-H-bond -> carboxyl-carboxyl dyad)."""
import json, sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
sys.path.insert(0, "/home/prachh/v/anthropic-adaptyv-2026-egfr/scripts")
from build_his_asp_theozyme import unit
from build_buried_theozyme import _rot_axis, _align

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
SRC = f"{REPO}/specs/rfd3/theozyme.pdb"
OUTDIR = f"{REPO}/specs/rfd3"
ORDER_GLU = ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"]
GLU_ELEM = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C", "CD": "C", "OE1": "O", "OE2": "O"}
GLU_CARB = ("CD", "OE1", "OE2")
DYAD = 2.6     # carboxyl O...O H-bond distance (shorter than His NE2..O; a shared-proton dyad)
CLASH_MIN = 2.5


def _glu_template(chA, template_author=431):
    this = chA[template_author]
    T = {n: this[n].coord.copy() for n in ORDER_GLU}
    carb_mid = np.mean([T[a] for a in ("OE1", "OE2")], axis=0)
    for n in ORDER_GLU:
        T[n] = T[n] - carb_mid                 # carboxyl O-midpoint -> origin
    point = unit(np.mean([T[a] for a in ("OE1", "OE2")], axis=0) - T["CG"])  # CG -> carboxyl direction
    return T, point


def build_glu_dyad(struct, asp_author, template_author=431):
    """Place a binder GLU whose carboxyl H-bonds the Asp460 carboxyl at dyad distance, buried in the groove
    (climb the mouth-normal, CB toward the mouth), maximizing target heavy-atom contacts around the carboxyl."""
    chA = struct[0]["A"]
    asp = chA[asp_author]
    od1, od2 = asp["OD1"].coord, asp["OD2"].coord
    o_mid = (od1 + od2) / 2.0
    centroid = np.mean([a.coord for r in chA for a in r], axis=0)
    u_out = unit(o_mid - centroid)
    o_acc = od1 if np.dot(od1 - centroid, u_out) >= np.dot(od2 - centroid, u_out) else od2
    # groove mouth-normal (floor = Asp460 mid, Phe412 ring, Gln408 CG -> solvent)
    phe = chA[412]; gln = chA[408]
    phe_cen = np.mean([phe[a].coord for a in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ") if a in phe], axis=0)
    floor = np.mean([o_mid, phe_cen, gln["CG"].coord], axis=0)
    u_app = unit(floor - centroid)
    T, point = _glu_template(chA, template_author)
    R0 = _align(point, -u_app)                 # carboxyl points BACK toward the Asp
    antigen = np.array([a.coord for r in chA for a in r if r.id[1] != asp_author])
    asp_os = [asp[o].coord for o in ("OD1", "OD2")]

    best = None; best_any = None
    # sweep the carboxyl-midpoint offset along the approach AND the spin about it; gate on the dyad O..O
    for d in (2.2, 2.5, 2.8, 3.1, 3.4, 3.7):
        carb_pos = o_acc + d * u_app
        for deg in range(0, 360, 4):
            R = _rot_axis(u_app, deg) @ R0
            placed = {n: R @ T[n] + carb_pos for n in ORDER_GLU}
            sc = np.array([placed[n] for n in ("CB", "CG", "CD", "OE1", "OE2")])
            clash = float(np.linalg.norm(antigen[None, :, :] - sc[:, None, :], axis=2).min())
            cb_out = float(np.dot(placed["CB"] - carb_pos, u_out))
            oo = min(np.linalg.norm(placed[o] - a) for o in ("OE1", "OE2") for a in asp_os)  # dyad O..O
            carb = np.array([placed[n] for n in GLU_CARB])
            contacts = int((np.linalg.norm(antigen[None, :, :] - carb[:, None, :], axis=2) < 6.0).any(axis=0).sum())
            if clash >= CLASH_MIN and cb_out >= 0 and 2.4 <= oo <= 3.0:
                if best is None or contacts > best[0]:
                    best = (contacts, placed, clash, cb_out, oo)
            if best_any is None or clash > best_any[2]:
                best_any = (contacts, placed, clash, cb_out, oo)
    if best is None:
        best = best_any
        print(f"  WARN: no rotamer met gates -> least-clashing (clash {best[2]:.2f}, O..O {best[4]:.2f})")
    contacts, placed, clash, cb_out, oo = best
    diag = dict(asp=asp_author, dyad_OO=round(float(oo), 2), clash=round(clash, 2),
                cb_out=round(cb_out, 2), carboxyl_target_contacts=contacts)
    return placed, diag


def write_theozyme(struct, placed, out_pdb):
    model = struct[0]
    if "C" in [c.id for c in model]:
        model.detach_child("C")
    ch = Chain("C")
    res = Residue((" ", 1, " "), "GLU", "")
    for i, n in enumerate(ORDER_GLU):
        res.add(Atom(n, placed[n].astype(float), 0.0, 1.0, " ", (" " + n).ljust(4), i + 1, GLU_ELEM[n]))
    ch.add(res); model.add(ch)
    io = PDBIO(); io.set_structure(struct); io.save(out_pdb)


def main():
    ASP = 436  # Asp460 (conserved human=mouse)
    struct = PDBParser(QUIET=True).get_structure("t", SRC)
    placed, diag = build_glu_dyad(struct, ASP)
    pdb = f"{OUTDIR}/theozyme_dyad_glu.pdb"; spec = f"{OUTDIR}/egfr_dyad_glu.json"
    write_theozyme(struct, placed, pdb)
    json.dump({"egfr_dyad_glu": {
        "dialect": 2, "input": pdb, "contig": "40-55,C1-1,45-60,/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL", "C1": "CB,CG,CD,OE1,OE2"},
        "select_hbond_donor": {"C1": "OE1,OE2"}, "select_hbond_acceptor": {f"A{ASP}": "OD1,OD2"},
        "select_hotspots": {f"A{ASP}": "OD1,OD2", "A412": "CG,CE1,CZ", "A409": "CB,CE1,NE2", "A408": "CB,CG"},
        "infer_ori_strategy": "hotspots", "is_non_loopy": True}}, open(spec, "w"), indent=2)
    print(f"dyad Glu: O..O={diag['dyad_OO']} clash={diag['clash']} cb_out={diag['cb_out']} "
          f"carboxyl_target_contacts={diag['carboxyl_target_contacts']} -> {pdb.split('/')[-1]}, {spec.split('/')[-1]}")


if __name__ == "__main__":
    main()
