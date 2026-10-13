#!/usr/bin/env python
"""Generic shell-adaptation input builder.

Given a binder-target complex pose (Boltz/Chai CIF or af2ig PDB), define the LigandMPNN
redesign set = chain-A (binder) residues in the His shell + interface, with the switch-His
positions FIXED. Writes a clean PDB (chain A = binder 1..N, chain B = target) for LigandMPNN
and prints a JSON with the redesigned-residues string.

Shell rule (mechanistic): redesign binder residue i if it is
  (a) within CONTACT_TARGET A (heavy-atom) of the target, OR
  (b) within CONTACT_HIS A of any switch-His heavy atom,
then REMOVE the switch-His positions themselves (they stay fixed as the pH trigger).
This repacks the environment that locks each His rotamer across folders (fold-robustness)
while keeping the far scaffold and the His identities fixed.

Usage: shelladapt_shell.py NAME STRUCT HIS_CSV OUT_PDB [CONTACT_TARGET=10.0] [CONTACT_HIS=8.0]
"""
import sys, json
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

name, struct, his_csv, out_pdb = sys.argv[1:5]
CT = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
CH = float(sys.argv[6]) if len(sys.argv) > 6 else 8.0
his_pos = [int(x) for x in his_csv.split(",") if x.strip()]

parser = MMCIFParser(QUIET=True) if struct.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
model = parser.get_structure("x", struct)[0]

chains = list(model)
def _npoly(c):
    return len([r for r in c if r.id[0] == ' '])
# binder = chain A if present (and non-empty) else the shortest POLYMER chain (never a het/solvent chain)
by_id = {c.id: c for c in chains}
polymer = [c for c in chains if _npoly(c) > 0]
if not polymer:
    print(json.dumps({"error": "no polymer chains in input (need a binder-target complex)"}))
    sys.exit(1)
binder = by_id.get("A") or min(polymer, key=_npoly)
targets = [c for c in chains if c is not binder]

HIS_NAMES = {"HIS", "HIE", "HID", "HIP", "HSD", "HSE", "HSP"}  # accept protonation-state names too
bres = [r for r in binder if r.id[0] == ' ']
# verify His present at switch positions
resmap = {r.id[1]: r for r in bres}
bad = [p for p in his_pos if resmap.get(p) is None or resmap[p].resname not in HIS_NAMES]
if bad:
    print(json.dumps({"error": f"HIS not found at {bad} in binder (got "
                      + ",".join(f"{p}:{resmap[p].resname if resmap.get(p) else 'NA'}" for p in his_pos) + ")"}))
    sys.exit(1)

def atoms(res):
    # heavy-atom coords only: the shell rule (CONTACT_TARGET/CONTACT_HIS) is a heavy-atom distance,
    # so explicit H (if a protonated input is ever passed) must not shrink the measured distances.
    hv = [a.coord for a in res if a.element != "H"]
    return np.array(hv) if hv else np.empty((0, 3))

tlist = [atoms(r) for c in targets for r in c if r.id[0] == ' ']
if not tlist:
    print(json.dumps({"error": "no target polymer residues found (need a binder-target complex)"}))
    sys.exit(1)
tatoms = np.concatenate(tlist)
# His-shell distance is optional: with no switch-His (target-contact-only run) hatoms is empty -> skip d_h
hatoms = np.concatenate([atoms(resmap[p]) for p in his_pos]) if his_pos else None

redesign = []
for r in bres:
    i = r.id[1]
    if i in his_pos:
        continue
    ra = atoms(r)
    d_t = np.min(np.linalg.norm(ra[:, None, :] - tatoms[None, :, :], axis=-1))
    d_h = (np.min(np.linalg.norm(ra[:, None, :] - hatoms[None, :, :], axis=-1))
           if hatoms is not None else np.inf)
    if d_t < CT or d_h < CH:
        redesign.append(i)

# write clean PDB: binder -> chain A, ALL target protomers kept -> chains B,C,D...
# (full trimer context for LigandMPNN; the redesign set above was computed against all protomers)
class Clean(Select):
    def accept_residue(self, res):
        return res.id[0] == ' '

# relabel binder -> 'A', targets -> 'B','C',... COLLISION-SAFE: assign unique temp ids first so a
# non-canonical input order (e.g. A,C,B) cannot hit Biopython's "id already used for a sibling".
used = "BCDEFGHIJK"
if len(targets) > len(used):
    print(json.dumps({"error": f"too many target chains ({len(targets)}) to relabel (max {len(used)})"}))
    sys.exit(1)
for k, c in enumerate([binder] + targets):
    c.id = f"__T{k}__"
binder.id = "A"
for ti, c in enumerate(targets):
    c.id = used[ti]
io = PDBIO()
io.set_structure(model)
io.save(out_pdb, Clean())

print(json.dumps({
    "name": name, "pdb": out_pdb,
    "his": his_pos, "n_binder": len(bres),
    "redesigned": " ".join(f"A{i}" for i in redesign),
    "n_redesign": len(redesign),
}))
