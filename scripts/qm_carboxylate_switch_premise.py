#!/usr/bin/env python3
"""WAVE 28 premise: the REPULSION-GATED pH switch. Pair the conserved target Asp460 carboxylate with a
BURIED binder carboxylate (Glu/Asp) whose pKa is up-shifted to ~6.5 by the desolvated interface pocket.

  blood pH 7.4: binder-Glu deprotonated (-)  -> Glu(-)...Asp(-) LIKE-CHARGE REPULSION -> strongly OFF
  tumor pH 6.5: binder-Glu protonated (COOH) -> GluH...Asp(-) H-bond                 -> ON

This inverts the stickiness that killed the buried salt bridge: there the NEUTRAL His still H-bonded the
desolvated carboxylate (off-state -12, constitutive). Here the OFF state is two like charges held in a low
dielectric -> they REPEL (off-state should be POSITIVE/repulsive, the opposite of sticky), and the ON state
is a neutral-carboxyl H-bond (less ECC-sensitive than a +/- salt bridge).

Model clusters (acetate = Asp/Glu sidechain), fixed salt-bridge geometry (O...O target), GFN2-xTB single
point E_int at two dielectrics (water eps80 surface, chcl3 eps4.8 buried). We report:
  dE_switch = E_int(ON: GluH + Asp-) - E_int(OFF: Glu- + Asp-)   (want < 0: protonation turns binding ON)
  E_int(OFF)  (want POSITIVE = repulsive = genuinely off, NOT sticky)
Compare to the buried His...Asp salt bridge (off -11.9, sticky). Also --opt the ON state for the real bound E."""
import argparse, json, os, subprocess, tempfile
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
XTB = f"{REPO}/envs/xtb/bin/xtb"
H2KCAL = 627.509474

ACETATE = ("CC(=O)[O-]", -1)      # Asp/Glu sidechain, deprotonated
ACETIC = ("CC(=O)O", 0)           # Glu sidechain, protonated (COOH)


def embed(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(m)
    return m


def coords(m):
    c = m.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])


def carboxyl(m):
    """return (C_idx, [O_idxs]) of the carboxyl group."""
    patt = Chem.MolFromSmarts("[CX3](=O)[OX1,OX2]")
    match = m.GetSubstructMatch(patt)
    c = match[0]
    os_ = [a.GetIdx() for a in m.GetAtomWithIdx(c).GetNeighbors() if a.GetSymbol() == "O"]
    return c, os_


def place_pair(m1, m2, oo=2.65):
    """orient both carboxyls to face along x, separate so nearest O...O ~= oo."""
    x1, x2 = coords(m1), coords(m2)
    c1, o1 = carboxyl(m1); c2, o2 = carboxyl(m2)
    def orient(x, c, os_, sign):
        omid = x[os_].mean(0)
        d = omid - x[c]; d = d / (np.linalg.norm(d) + 1e-9)
        # rotate so carboxyl points along sign*x-axis
        tgt = np.array([sign, 0.0, 0.0])
        v = np.cross(d, tgt); s = np.linalg.norm(v); cth = np.dot(d, tgt)
        if s < 1e-8:
            R = np.eye(3) if cth > 0 else np.diag([1.0, -1, -1])
        else:
            vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            R = np.eye(3) + vx + vx @ vx * ((1 - cth) / (s * s))
        xr = (x - x[c]) @ R.T
        return xr
    x1 = orient(x1, c1, o1, +1.0)    # frag1 carboxyl points +x
    x2 = orient(x2, c2, o2, -1.0)    # frag2 carboxyl points -x (facing frag1)
    # separate along x: place frag2 so its nearest O is ~oo from frag1's nearest O
    o1max = x1[o1][np.argmax(x1[o1][:, 0])]       # frag1 O with largest x
    o2min = x2[o2][np.argmin(x2[o2][:, 0])]       # frag2 O with smallest x (after facing)
    shift = (o1max[0] - o2min[0]) + oo
    x2 = x2 + np.array([shift, 0, 0])
    return x1, x2


def write_xyz(path, syms, xyz):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n\n")
        for s, p in zip(syms, xyz):
            f.write(f"{s} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}\n")


def xtb_e(xyz_name, chrg, wd, solvent, opt=False):
    cmd = [XTB, xyz_name, "--gfn", "2", "--chrg", str(chrg), "--opt" if opt else "--sp"]
    if solvent != "vacuum":
        cmd[4:4] = ["--alpb", solvent]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
    e = None
    for ln in r.stdout.splitlines():
        if "TOTAL ENERGY" in ln:
            e = float(ln.split()[-3])
    if e is None:
        raise RuntimeError(f"xtb failed chrg={chrg}: {r.stdout[-300:]} {r.stderr[-150:]}")
    return e


def e_int(frag1_smi, q1, frag2_smi, q2, solvent, oo, opt=False):
    m1, m2 = embed(frag1_smi), embed(frag2_smi)
    x1, x2 = place_pair(m1, m2, oo)
    s1 = [a.GetSymbol() for a in m1.GetAtoms()]; s2 = [a.GetSymbol() for a in m2.GetAtoms()]
    with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
        write_xyz(os.path.join(wd, "c.xyz"), s1 + s2, np.vstack([x1, x2]))
        write_xyz(os.path.join(wd, "a.xyz"), s1, x1)
        write_xyz(os.path.join(wd, "b.xyz"), s2, x2)
        Ec = xtb_e("c.xyz", q1 + q2, wd, solvent, opt)
        Ea = xtb_e("a.xyz", q1, wd, solvent)
        Eb = xtb_e("b.xyz", q2, wd, solvent)
    return (Ec - Ea - Eb) * H2KCAL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oo", type=float, default=2.65, help="target carboxyl O...O distance (A)")
    ap.add_argument("--out", default="runs/wave28_premise/carboxylate_switch.json")
    a = ap.parse_args()
    rows = []
    for solvent in ("water", "chcl3"):
        # ON = protonated binder-Glu (acetic acid, 0) + target Asp (acetate, -1); H-bond
        on_sp = e_int(ACETIC[0], ACETIC[1], ACETATE[0], ACETATE[1], solvent, a.oo, opt=False)
        on_opt = e_int(ACETIC[0], ACETIC[1], ACETATE[0], ACETATE[1], solvent, a.oo, opt=True)
        # OFF = deprotonated binder-Glu (acetate, -1) + target Asp (acetate, -1); like-charge repulsion
        off_sp = e_int(ACETATE[0], ACETATE[1], ACETATE[0], ACETATE[1], solvent, a.oo, opt=False)
        row = dict(solvent=solvent, oo=a.oo, Eint_ON_sp=round(on_sp, 2), Eint_ON_opt=round(on_opt, 2),
                   Eint_OFF_sp=round(off_sp, 2), dE_switch_sp=round(on_sp - off_sp, 2),
                   off_is_repulsive=bool(off_sp > 0))
        rows.append(row); print(json.dumps(row), flush=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(rows, open(a.out, "w"), indent=1)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
