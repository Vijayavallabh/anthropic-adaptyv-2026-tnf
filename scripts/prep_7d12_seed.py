#!/usr/bin/env python3
"""Prep the 7D12-nanobody-seeded partial-diffusion input.

7D12 (PDB 4KRL chain B, a VHH) binds EGFR Domain III (4KRL chain A) at pH 6.0 and DIRECTLY contacts
His433 (= mature residue 409, 3.6 A) -- an ideal, validated pH-region binding pose to seed from.
4KRL and our canonical target data/targets/6ARU_domainIII.pdb both use MATURE EGFR numbering, so we
Kabsch-superpose 4KRL chain A onto the 6ARU DIII target by matched residue number, apply the same
transform to 7D12, and emit a complex in OUR target frame:
    chain A = 7D12 binder (renumbered 1..Lb)
    chain B = 6ARU DIII target (as-is; His433 = mature 409)
Then scripts/canon_renumber.py renumbers B -> 1..191 (His433 at ordinal 100) for run_rfd_partial.sh.
Usage: prep_7d12_seed.py [4KRL.pdb] [6ARU_domainIII.pdb] [out.pdb]
"""
import sys
import numpy as np

KRL = sys.argv[1] if len(sys.argv) > 1 else "runs/7d12/4KRL.pdb"
TGT = sys.argv[2] if len(sys.argv) > 2 else "data/targets/6ARU_domainIII.pdb"
OUT = sys.argv[3] if len(sys.argv) > 3 else "runs/7d12/7d12_on_6aru.pdb"


def parse(path, chain, ca_only=False):
    """return list of (resnum, resname, atomname, (x,y,z), raw_line)"""
    out = []
    for l in open(path):
        if l[:4] != "ATOM":
            continue
        if l[21] != chain:
            continue
        an = l[12:16].strip()
        if ca_only and an != "CA":
            continue
        try:
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError:
            continue
        out.append((int(l[22:26]), l[17:20].strip(), an, xyz, l))
    return out


def kabsch(P, Q):
    """R,t so that R@P.T + t ~ Q (P,Q are Nx3). Classic Kabsch."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


def main():
    krl_egfr_ca = {r[0]: r[3] for r in parse(KRL, "A", ca_only=True)}   # 4KRL EGFR (mature num)
    tgt_ca = {r[0]: r[3] for r in parse(TGT, "A", ca_only=True)}        # 6ARU DIII (mature num)
    common = sorted(set(krl_egfr_ca) & set(tgt_ca))
    if len(common) < 20:
        sys.exit(f"too few matched EGFR residues ({len(common)}) between 4KRL and 6ARU -- numbering mismatch")
    P = np.array([krl_egfr_ca[i] for i in common])
    Q = np.array([tgt_ca[i] for i in common])
    R, t = kabsch(P, Q)
    rmsd = np.sqrt(((R @ P.T).T + t - Q) ** 2).sum(1).mean() ** 0.5
    print(f"superposed 4KRL EGFR -> 6ARU DIII on {len(common)} matched CA, RMSD {rmsd:.2f} A")

    # transform 7D12 (4KRL chain B) into the 6ARU frame, write as binder chain A renumbered 1..Lb
    vhh = parse(KRL, "B")
    resmap, nxt = {}, 0
    lines = []
    for rn, rname, an, xyz, raw in vhh:
        if rn not in resmap:
            nxt += 1; resmap[rn] = nxt
        x, y, z = (R @ np.array(xyz)) + t
        newnum = resmap[rn]
        lines.append(f"ATOM  {len(lines)+1:5d} {raw[12:16]} {raw[17:20]} A{newnum:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{raw[54:66] if len(raw) > 66 else '  1.00  0.00'}          "
                     f"{raw[76:78] if len(raw) > 78 else ' C'}\n")
    nbind = len(resmap)
    # target chain B = 6ARU DIII, as-is (mature numbering preserved; canon_renumber will 1..191 it)
    for r in parse(TGT, "A"):
        raw = r[4]
        lines.append("ATOM  " + f"{len(lines)+1:5d}" + raw[11:21] + "B" + raw[22:])
    with open(OUT, "w") as fh:
        fh.writelines(lines); fh.write("TER\nEND\n")
    print(f"wrote {OUT}: binder(7D12) chain A len {nbind} + target(6ARU DIII) chain B")


if __name__ == "__main__":
    main()
