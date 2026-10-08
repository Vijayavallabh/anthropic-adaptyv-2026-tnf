#!/usr/bin/env python3
"""Build AF2-ig inputs for the TNF-alpha TWO-protomer target (full-atom).

RFdiffusion/MPNN output here = target chains A+B (backbone-only) + binder chain C.
AF2 initial-guess needs the target with sidechains (backbone-only target -> plddt_target ~40,
0% yield), and expects binder = chain A. RFd keeps the A+B target rigid (one rigid-body move of
the whole system), so we Kabsch-superpose the real full-atom 1TNF_AB (chains A,B) onto each
design's target backbone (matched per chain by residue order) and emit:
  chain A = binder (designed, backbone; renumbered 1..Lb)
  chain B = full-atom target protomer 1 (real A, in the design frame)
  chain C = full-atom target protomer 2 (real B, in the design frame)

Usage:
  build_fullatom_af2in_tnf.py --mpnn 'runs/.../backbones/*.pdb' \
     --target data/targets/1TNF_AB.pdb --out runs/.../af2in [--binder-chain C]
"""
import argparse, glob, os
import numpy as np


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cP, cQ


def xyz_of(ln):
    return np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])


def put_xyz(ln, xyz):
    return "%s%8.3f%8.3f%8.3f%s" % (ln[:30], xyz[0], xyz[1], xyz[2], ln[54:].rstrip("\n"))


def read_chain(path, chain, backbone_only=False):
    """ordered (resnum, line) list + {resnum: CA xyz} for one chain's ATOM records."""
    lines, ca = [], {}
    for ln in open(path):
        if ln.startswith("ATOM") and ln[21] == chain:
            an = ln[12:16].strip(); rn = int(ln[22:26])
            if backbone_only and an not in ("N", "CA", "C", "O"):
                continue
            lines.append((rn, ln))
            if an == "CA":
                ca[rn] = xyz_of(ln)
    return lines, ca


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpnn", required=True, help="glob of design complexes (target A+B + binder)")
    ap.add_argument("--target", default="data/targets/1TNF_AB.pdb")
    ap.add_argument("--target-chains", default="A,B")
    ap.add_argument("--binder-chain", default="C")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-rmsd", type=float, default=2.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tchains = a.target_chains.split(",")

    # real full-atom target, per chain
    real = {c: read_chain(a.target, c) for c in tchains}

    n_ok = n_skip = 0; worst = 0.0
    for comp in sorted(glob.glob(a.mpnn)):
        # design-frame target CA (per chain) + binder lines
        dca = {c: read_chain(comp, c)[1] for c in tchains}
        binder = read_chain(comp, a.binder_chain)[0]
        if not binder or any(len(dca[c]) == 0 for c in tchains):
            n_skip += 1; continue
        # matched CA arrays, per chain by residue order
        P, Q = [], []
        ok = True
        for c in tchains:
            rl = sorted(real[c][1]); dl = sorted(dca[c])
            m = min(len(rl), len(dl))
            if m < 40:
                ok = False; break
            # positional pairing (real/design differ only by a constant resnum offset);
            # if a gap breaks that constant offset, skip rather than mis-align silently.
            off = rl[0] - dl[0]
            if any(rl[i] - dl[i] != off for i in range(m)):
                ok = False; break
            P += [real[c][1][rl[i]] for i in range(m)]
            Q += [dca[c][dl[i]] for i in range(m)]
        if not ok:
            n_skip += 1; continue
        P = np.array(P); Q = np.array(Q)
        R, cP, cQ = kabsch(P, Q)
        rmsd = float(np.sqrt(np.mean(np.sum(((P - cP) @ R.T + cQ - Q) ** 2, axis=1))))
        worst = max(worst, rmsd)
        if rmsd > a.max_rmsd:
            n_skip += 1; continue

        out = []; serial = 1
        # binder -> chain A, renumber 1..Lb
        resmap = {}; nr = 0
        for rn, ln in binder:
            if rn not in resmap:
                nr += 1; resmap[rn] = nr
            out.append("%s%5d %sA%4d %s\n" % (ln[:6], serial, ln[12:21], resmap[rn], ln[27:].rstrip("\n")))
            serial += 1
        out.append("TER\n")
        # target protomers -> chains B, C (full-atom, moved into the design frame)
        for c, outc in zip(tchains, ["B", "C"]):
            resmap = {}; nr = 0
            for rn, ln in real[c][0]:
                if rn not in resmap:
                    nr += 1; resmap[rn] = nr
                nl = put_xyz(ln, (xyz_of(ln) - cP) @ R.T + cQ)
                out.append("%s%5d %s%s%4d %s\n" % (nl[:6], serial, nl[12:21], outc, resmap[rn], nl[27:].rstrip("\n")))
                serial += 1
            out.append("TER\n")
        out.append("END\n")
        open(os.path.join(a.out, os.path.basename(comp)), "w").write("".join(out))
        n_ok += 1
    print(f"wrote {n_ok} af2ig inputs (full-atom 2-protomer target) -> {a.out} | skipped {n_skip} | worst RMSD {worst:.2f} A")


if __name__ == "__main__":
    main()
