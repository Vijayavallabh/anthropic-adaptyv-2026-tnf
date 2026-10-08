#!/usr/bin/env python3
"""WAVE 23 MOONSHOT: a COOPERATIVE (Hill n>1) tumor-ON pH switch to break the single-His ~5x ceiling.

A single His switch is capped at ~4.5x fold-selectivity over pH 6.5->7.4 (Henderson-Hasselbalch). TWO
protonation-gated contacts that are BOTH required for binding give an effective Hill n~2 -> ~25x, and n~3
-> ~250x (Baker 2025 found sharpness empirically correlates with >=2 protonated-His contacts). The dual-His
NO-GO was two His bridging the SAME region (like-charge repulsion -> pKa SPLITTING -> anti-cooperative
broadening). We avoid that by placing two binder His on SPATIALLY SEPARATED, conserved antigen acids:
  His_a (chain C res 1)  ->  Asp460 (A436)   [proven-bindable His433 face, 10.8A from A409]
  His_b (chain C res 2)  ->  Asp344 (A344)   [12.0A from A409]
Asp460..Asp344 = 15.1A apart: far enough that the two His+ do not electrostatically repel (no splitting),
close enough for one 70-130aa binder to span. BOTH acids are conserved human=mouse (D/D) -> the switch is
cross-reactive by construction AND tumor-ON (both His+ bridges form at 6.5, both break at 7.4).

Output: a two-His RFd3 theozyme PDB + spec (two motif His, two hotspots, dock at both Asp). Reuses the
hardened clash-gated rotamer placement from build_his_asp_theozyme.build_his_for_asp.
"""
import os
import json, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_his_asp_theozyme import build_his_for_asp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = f"{REPO}/specs/rfd3/theozyme.pdb"
OUTDIR = f"{REPO}/specs/rfd3"
ELEM = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C",
        "ND1": "N", "CD2": "C", "CE1": "C", "NE2": "N"}
ORDER = ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"]


def write_two_his_theozyme(struct, placed_a, placed_b, out_pdb):
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


def write_coop_spec(out_pdb, asp_a, asp_b, out_json):
    """two-His motif: binder spans both separated conserved Asp; dock at both + the His433 face."""
    spec = {"egfr_coop_dualhis": {
        "dialect": 2,
        "input": out_pdb,
        # larger binder that presents His_a (C1) and His_b (C2) ~15A apart to the two Asp
        "contig": f"25-45,C1-1,20-40,C2-2,25-45,/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL",
                               "C1": "CB,CG,ND1,CD2,CE1,NE2", "C2": "CB,CG,ND1,CD2,CE1,NE2"},
        "select_hbond_donor": {"C1": "ND1,NE2", "C2": "ND1,NE2"},   # both protonated His donate
        "select_hbond_acceptor": {f"A{asp_a}": "OD1,OD2", f"A{asp_b}": "OD1,OD2"},
        "select_hotspots": {"A409": "ND1,NE2", f"A{asp_a}": "OD1,OD2", f"A{asp_b}": "OD1,OD2"},
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
    }}
    json.dump(spec, open(out_json, "w"), indent=2)


def main():
    struct = PDBParser(QUIET=True).get_structure("t", SRC)
    ASP_A, ASP_B = 436, 344   # Asp460 (face, central) + Asp344 (15.1A away), both conserved human=mouse
    placed_a, diag_a = build_his_for_asp(struct, ASP_A)
    placed_b, diag_b = build_his_for_asp(struct, ASP_B)
    # verify the two His do not clash each other (separated sites should be safe)
    sc_a = np.array([placed_a[n] for n in ("CB", "CG", "ND1", "CD2", "CE1", "NE2")])
    sc_b = np.array([placed_b[n] for n in ("CB", "CG", "ND1", "CD2", "CE1", "NE2")])
    his_his = float(np.linalg.norm(sc_a[:, None, :] - sc_b[None, :, :], axis=2).min())
    ne2_ne2 = float(np.linalg.norm(placed_a["NE2"] - placed_b["NE2"]))
    out_pdb = f"{OUTDIR}/theozyme_coop_dualhis.pdb"
    out_json = f"{OUTDIR}/egfr_coop_dualhis.json"
    write_two_his_theozyme(struct, placed_a, placed_b, out_pdb)
    write_coop_spec(out_pdb, ASP_A, ASP_B, out_json)
    print(f"His_a->Asp460(A436): NE2..O={diag_a['ne2_O']}A clash={diag_a['clash']}A CB_out={diag_a['cb_outward']}")
    print(f"His_b->Asp344(A344): NE2..O={diag_b['ne2_O']}A clash={diag_b['clash']}A CB_out={diag_b['cb_outward']}")
    print(f"His_a .. His_b sidechain min-dist = {his_his:.1f}A  (NE2..NE2 = {ne2_ne2:.1f}A; want >8A, no repulsion)")
    print(f"-> {out_pdb}\n-> {out_json}")
    if his_his < 6.0:
        print("  WARNING: the two His are closer than 6A -> like-charge repulsion risk (anti-cooperativity)")
    for tag, d in (("His_a", diag_a), ("His_b", diag_b)):
        if d["ne2_O"] > 3.3 or d["clash"] < 2.0:
            print(f"  WARNING {tag}: geometry marginal (NE2..O {d['ne2_O']}, clash {d['clash']})")


if __name__ == "__main__":
    main()
