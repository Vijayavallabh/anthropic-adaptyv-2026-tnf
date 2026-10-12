#!/usr/bin/env python
"""Convert a Boltz/Chai binder-target complex CIF to an af2ig-ready PDB:
binder -> chain A (1..N), first target protomer -> chain B, full-atom (keeps target sidechains).
Usage: cif2pdb.py IN.cif OUT.pdb [BINDER_CHAIN=A]"""
import sys
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

inp, outp = sys.argv[1], sys.argv[2]
bchain = sys.argv[3] if len(sys.argv) > 3 else "A"
P = MMCIFParser(QUIET=True) if inp.lower().endswith(".cif") else PDBParser(QUIET=True)
model = P.get_structure("x", inp)[0]
chains = [c for c in model]
by = {c.id: c for c in chains}
binder = by.get(bchain, min(chains, key=lambda c: len([r for r in c if r.id[0] == ' '])))
targets = [c for c in chains if c is not binder]
# keep only first target protomer (af2ig construct = binder + 1 protomer, 152aa)
target = max(targets, key=lambda c: len([r for r in c if r.id[0] == ' ']))
binder.id = "A"; target.id = "B"
keep = {"A", "B"}

class Sel(Select):
    def accept_chain(self, c): return c.id in keep
    def accept_residue(self, r): return r.id[0] == ' '

io = PDBIO(); io.set_structure(model); io.save(outp, Sel())
nb = len([r for r in binder if r.id[0] == ' '])
nt = len([r for r in target if r.id[0] == ' '])
print(f"{outp} binderA={nb} targetB={nt}")
