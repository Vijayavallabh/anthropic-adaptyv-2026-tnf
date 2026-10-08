#!/usr/bin/env python3
"""Adaptyv protein-qc mirror: shape complementarity (SC) + interface_dG (Rosetta InterfaceAnalyzer)
for binder(chain A):target(chain B) complexes. Gates: SC>0.50/0.62, interface_dG<-10/-15.
Usage: interface_qc.py <complex1.pdb> [complex2.pdb ...]   (run in envs/pyrosetta)
"""
import sys
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
print(f"{'complex':52} {'SC':>5} {'dG_sep':>8} {'dSASA':>8} {'nres_int':>8}")
print(f"{'':52} {'(>.5)':>5} {'(<-10)':>8}")
for pdb in sys.argv[1:]:
    try:
        pose = pose_from_pdb(pdb)
        iam = InterfaceAnalyzerMover("A_B")
        iam.set_compute_interface_sc(True)
        iam.set_compute_packstat(False)
        iam.set_calc_dSASA(True)
        iam.apply(pose)
        d = iam.get_all_data()
        sc = d.sc_value
        dG = iam.get_interface_dG()
        dsasa = iam.get_interface_delta_sasa()
        nres = d.interface_nres if hasattr(d, "interface_nres") else 0
        import os
        nm = os.path.basename(pdb).replace("_af2pred.pdb", "").replace(".pdb", "")
        flag = "" if (sc > 0.5 and dG < -10) else "  <-- below a gate"
        print(f"{nm[:52]:52} {sc:5.2f} {dG:8.1f} {dsasa:8.0f} {str(nres):>8}{flag}")
    except Exception as e:
        print(f"{pdb}: ERROR {e}")
