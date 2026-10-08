#!/usr/bin/env python3
"""Convert a Boltz cif complex to pdb (preserve chain ids A=binder,B=target)."""
import sys
from Bio.PDB import MMCIFParser, PDBIO
s = MMCIFParser(QUIET=True).get_structure("x", sys.argv[1])
io = PDBIO(); io.set_structure(s); io.save(sys.argv[2])
print("wrote", sys.argv[2])
