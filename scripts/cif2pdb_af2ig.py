#!/usr/bin/env python
"""Convert a Boltz/Chai binder-target complex CIF to an af2ig-ready PDB:
binder -> chain A, the largest target protomer -> chain B, full-atom (keeps target sidechains).
Fail-loud: errors out (exit 2) rather than silently emitting a binder-less/empty construct if the
binder cannot be identified, if there is no target chain, or if the chosen "binder" is target-sized.
Usage: cif2pdb_af2ig.py IN.cif OUT.pdb [BINDER_CHAIN=A] [MAX_BINDER_LEN=130]
NOTE: residue numbers are PRESERVED (not renumbered); binder is relabeled chain A, target chain B."""
import sys
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

def die(msg):
    sys.stderr.write(f"cif2pdb_af2ig: ERROR: {msg}\n"); sys.exit(2)

inp, outp = sys.argv[1], sys.argv[2]
bchain = sys.argv[3] if len(sys.argv) > 3 else "A"
MAXB = int(sys.argv[4]) if len(sys.argv) > 4 else 130  # a binder longer than this is almost certainly a target
P = MMCIFParser(QUIET=True) if inp.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
model = P.get_structure("x", inp)[0]

def npoly(c):  # count standard polymer residues (blank hetflag)
    return len([r for r in c if r.id[0] == ' '])

chains = [c for c in model if npoly(c) > 0]
if not chains:
    die("no polymer chains in input")
by = {c.id: c for c in chains}
if len(sys.argv) > 3:
    # an explicit binder chain was requested -> it MUST exist (fail-loud, per this script's contract);
    # silently falling back to the smallest chain would emit a construct built around the wrong chain.
    if bchain not in by:
        die(f"requested binder chain {bchain!r} not found (polymer chains: {sorted(by)})")
    binder = by[bchain]
else:
    binder = by.get("A") or min(chains, key=npoly)
targets = [c for c in chains if c is not binder]
if not targets:
    die("only one polymer chain present — need a binder-target complex")
target = max(targets, key=npoly)
nb, nt = npoly(binder), npoly(target)
if nb > MAXB and nt <= nb:
    die(f"identified binder (chain {binder.id}, {nb} res) is target-sized (>{MAXB}) and not shorter than the "
        f"target (chain {target.id}, {nt} res) — chain A is probably NOT the binder; pass the binder chain id explicitly")

# collision-safe: detach every chain except binder+target BEFORE relabeling, then rename
keep_ids = {binder.id, target.id}
for c in list(model):
    if c.id not in keep_ids:
        model.detach_child(c.id)
# rename via temporary ids to avoid an A/B swap collision
binder.id = "__BINDER__"; target.id = "__TARGET__"
binder.id = "A"; target.id = "B"

class Sel(Select):
    def accept_residue(self, r):
        return r.id[0] == ' '

io = PDBIO(); io.set_structure(model); io.save(outp, Sel())
print(f"{outp} binderA={nb} targetB={nt}")
