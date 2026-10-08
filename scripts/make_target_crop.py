#!/usr/bin/env python3
"""Crop 6ARU chain A to the EGFR Domain III epitope region for BoltzGen target conditioning.

Produces a clean single-chain PDB (chain A only, standard amino acids, glycans/other
chains stripped) keeping the ORIGINAL 6ARU author numbering, plus a JSON mapping that
records, for each hotspot / avoid residue, its 1-based SEQUENTIAL ordinal among the
residues present in the crop. BoltzGen's spec `binding` / `not_binding` / `res_index`
select residues by sequential position within the parsed structure (parse_range over the
chain's residue array), NOT by author number -- so the spec must use these ordinals.

Usage:
  python scripts/make_target_crop.py [LO] [HI]     # 6ARU author resSeq window (default 310 500)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "data" / "targets" / "6ARU.pdb"
OUT_PDB = REPO / "data" / "targets" / "6ARU_domainIII.pdb"
OUT_MAP = REPO / "data" / "targets" / "6ARU_domainIII_map.json"

CHAIN = "A"
# 6ARU chain A author resSeq -> P00533 for reference (from data/epitope/structure_epitope.json)
# Conserved cross-reactive anchors (human==mouse) -> BINDING hotspots:
CONSERVED_6ARU = {384: 408, 408: 432, 409: 433, 412: 436, 440: 464, 465: 489, 472: 496}
# Cetuximab contacts that DIFFER human vs mouse -> NOT_BINDING (steer the interface away):
DIVERGENT_6ARU = {353: 377, 418: 442, 443: 467, 467: 491, 468: 492, 473: 497}
PH_SENSOR_6ARU = 409  # His433
ACIDIC_PARTNER_6ARU = 472  # Glu496

STANDARD3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def main() -> None:
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 310
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else 500

    kept_lines: list[str] = []
    order: list[tuple[int, str]] = []  # unique (resSeq, icode) in file order
    resname_by_key: dict[tuple[int, str], str] = {}

    for line in SRC.read_text().splitlines():
        rec = line[:6].strip()
        if rec != "ATOM" and not (rec == "HETATM" and line[17:20].strip() == "MSE"):
            continue
        if line[21] != CHAIN:
            continue
        resn = line[17:20].strip()
        if resn == "MSE":
            resn_out = "MET"
        elif resn in STANDARD3:
            resn_out = resn
        else:
            continue
        # altloc: keep blank or 'A' only
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if not (lo <= resseq <= hi):
            continue
        icode = line[26]
        key = (resseq, icode)
        if key not in resname_by_key:
            order.append(key)
            resname_by_key[key] = resn_out
        # normalize: force record ATOM, resname (MSE->MET), altloc blank, chain A
        newline = "ATOM  " + line[6:16] + " " + resn_out.ljust(3) + line[20:]
        kept_lines.append(newline)

    if not order:
        raise SystemExit(f"No residues kept in {CHAIN}:{lo}-{hi}; check inputs.")

    ordinal = {key: i + 1 for i, key in enumerate(order)}  # 1-based sequential

    def ords(author_set) -> dict[int, int]:
        out = {}
        for a in author_set:
            k = (a, " ")
            if k in ordinal:
                out[a] = ordinal[k]
        return out

    conserved_ord = ords(CONSERVED_6ARU)
    divergent_ord = ords(DIVERGENT_6ARU)

    # Write PDB
    with OUT_PDB.open("w") as fh:
        fh.write(f"REMARK  EGFR Domain III crop of 6ARU chain A, author resSeq {lo}-{hi}\n")
        fh.write("REMARK  standard residues only; glycans/other chains stripped; author numbering kept\n")
        for ln in kept_lines:
            fh.write(ln.rstrip() + "\n")
        fh.write("TER\nEND\n")

    n = len(order)
    first_auth, last_auth = order[0][0], order[-1][0]
    result = {
        "source": "6ARU.pdb chain A",
        "crop_range_6aru": [lo, hi],
        "n_residues_present": n,
        "first_author_resseq": first_auth,
        "last_author_resseq": last_auth,
        "numbering_note": "PDB keeps 6ARU author numbering; ordinals below are 1-based sequential positions among present residues (what BoltzGen spec binding/not_binding/res_index use).",
        "ph_sensor_his433": {"6aru": PH_SENSOR_6ARU, "ordinal": ordinal.get((PH_SENSOR_6ARU, " "))},
        "acidic_partner_glu496": {"6aru": ACIDIC_PARTNER_6ARU, "ordinal": ordinal.get((ACIDIC_PARTNER_6ARU, " "))},
        "conserved_binding": {str(a): {"p00533": CONSERVED_6ARU[a], "ordinal": o} for a, o in conserved_ord.items()},
        "divergent_not_binding": {str(a): {"p00533": DIVERGENT_6ARU[a], "ordinal": o} for a, o in divergent_ord.items()},
        "binding_ordinals_csv": ",".join(str(o) for o in sorted(conserved_ord.values())),
        "not_binding_ordinals_csv": ",".join(str(o) for o in sorted(divergent_ord.values())),
    }
    OUT_MAP.write_text(json.dumps(result, indent=2))

    print(f"Wrote {OUT_PDB} ({n} residues, author {first_auth}-{last_auth})")
    print(f"Wrote {OUT_MAP}")
    print("BINDING (conserved) ordinals:", result["binding_ordinals_csv"])
    print("NOT_BINDING (divergent) ordinals:", result["not_binding_ordinals_csv"])
    print("His433 ordinal:", result["ph_sensor_his433"]["ordinal"],
          "| Glu496 ordinal:", result["acidic_partner_glu496"]["ordinal"])


if __name__ == "__main__":
    main()
