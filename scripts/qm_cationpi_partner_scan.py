#!/usr/bin/env python3
"""WAVE 27 premise check: which aromatic partner gives the BEST pH-gated cation-pi switch?

Builds idealized capped model clusters (His sidechain = 4-methylimidazole(ium); partner = toluene[Phe],
3-methylindole[Trp], p-cresol[Tyr], methylguanidinium-free) in a canonical stacked cation-pi geometry
(imidazolium ring centroid ~3.5 A above the partner ring centroid, faces parallel), then GFN2-xTB --opt
E_int in the HIP (protonated, on) and neutral (off) states, at two dielectrics (water eps80 = surface,
chcl3 eps4.8 = buried protein interior).

The switch GAP = E_int(HIP) - E_int(neutral) is the selectivity; the OFF = E_int(neutral) is the
constitutivity leak (must stay weak, > ~-4 kcal/mol, to stay unoccupied at blood pH). A good switch wants
a LARGE gap AND a WEAK off. This fairly ranks Trp vs Phe vs Tyr as the cation-pi partner, correcting the
fixed-charge MD bias (user directive) and extending the WAVE24 Phe-only result."""
import argparse, json, os, subprocess, tempfile
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XTB = f"{REPO}/envs/xtb/bin/xtb"
H2KCAL = 627.509474

# capped sidechain models: methyl stands in for CB (His) / the attachment carbon
HIS_HIP = ("Cc1c[nH+]c[nH]1", 1)   # 4-methylimidazolium (both N protonated), +1
HIS_NEU = ("Cc1cnc[nH]1", 0)       # 4-methylimidazole (neutral, one N-H), 0
PARTNERS = {
    "Phe_toluene":     ("Cc1ccccc1", 0),
    "Trp_3Me_indole":  ("Cc1c[nH]c2ccccc12", 0),
    "Tyr_p_cresol":    ("Cc1ccc(O)cc1", 0),
}


def embed(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(m)
    return m


def aromatic_ring_atoms(m, prefer_six=True):
    ri = m.GetRingInfo()
    rings = [r for r in ri.AtomRings() if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    if not rings:
        raise RuntimeError("no aromatic ring")
    rings.sort(key=lambda r: (len(r) != (6 if prefer_six else 5), -len(r)))
    return list(rings[0])


def coords(m):
    c = m.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])


def elements(m):
    return [a.GetSymbol() for a in m.GetAtoms()]


def ring_frame(xyz, ring):
    cen = xyz[ring].mean(0)
    # best-fit plane normal via SVD of centered ring coords
    u, s, vt = np.linalg.svd(xyz[ring] - cen)
    n = vt[2]
    return cen, n / np.linalg.norm(n)


def place_stacked(his_xyz, his_ring, par_xyz, par_ring, sep=3.5, offset=1.0):
    """translate+rotate HIS so its ring is parallel to and `sep` A above the partner ring centroid,
    with a small lateral `offset` (slipped-stacked, avoids perfect eclipse)."""
    pc, pn = ring_frame(par_xyz, par_ring)
    hc, hn = ring_frame(his_xyz, his_ring)
    # rotate His so hn -> pn
    v = np.cross(hn, pn); s = np.linalg.norm(v); c = np.dot(hn, pn)
    if s < 1e-6:
        R = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))
    hx = (his_xyz - hc) @ R.T
    # lateral offset direction in the partner plane
    ref = np.array([1.0, 0, 0]); ref = ref - np.dot(ref, pn) * pn
    ref = ref / (np.linalg.norm(ref) + 1e-9)
    target = pc + sep * pn + offset * ref
    return hx + target


def write_xyz(path, els, xyz):
    with open(path, "w") as f:
        f.write(f"{len(els)}\n\n")
        for e, p in zip(els, xyz):
            f.write(f"{e} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}\n")


def xtb_e(xyz_name, chrg, wd, solvent):
    cmd = [XTB, xyz_name, "--gfn", "2", "--chrg", str(chrg), "--opt"]
    if solvent != "vacuum":
        cmd[4:4] = ["--alpb", solvent]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
    e = None
    for ln in r.stdout.splitlines():
        if "TOTAL ENERGY" in ln:
            e = float(ln.split()[-3])
    if e is None:
        raise RuntimeError(f"xtb failed chrg={chrg}: {r.stdout[-300:]} {r.stderr[-200:]}")
    return e


def e_int(his_smi, his_q, par_smi, par_q, solvent, wd, prefer_six=True):
    hm, pm = embed(his_smi), embed(par_smi)
    hx, px = coords(hm), coords(pm)
    hr = aromatic_ring_atoms(hm, prefer_six=False)  # imidazole 5-ring
    pr = aromatic_ring_atoms(pm, prefer_six=prefer_six)
    hx2 = place_stacked(hx, hr, px, pr)
    els = elements(hm) + elements(pm)
    allx = np.vstack([hx2, px])
    write_xyz(os.path.join(wd, "cplx.xyz"), els, allx)
    write_xyz(os.path.join(wd, "his.xyz"), elements(hm), hx2)
    write_xyz(os.path.join(wd, "par.xyz"), elements(pm), px)
    Ec = xtb_e("cplx.xyz", his_q + par_q, wd, solvent)
    Eh = xtb_e("his.xyz", his_q, wd, solvent)
    Ep = xtb_e("par.xyz", par_q, wd, solvent)
    return (Ec - Eh - Ep) * H2KCAL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solvents", default="water,chcl3")
    ap.add_argument("--out", default="runs/wave27_premise/partner_scan.json")
    a = ap.parse_args()
    rows = []
    for solvent in a.solvents.split(","):
        for pname, (psmi, pq) in PARTNERS.items():
            prefer_six = "indole" not in pname  # Trp: use the electron-rich 5-ring (pyrrole) for cation-pi
            with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
                e_on = e_int(HIS_HIP[0], HIS_HIP[1], psmi, pq, solvent, wd, prefer_six)
            with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
                e_off = e_int(HIS_NEU[0], HIS_NEU[1], psmi, pq, solvent, wd, prefer_six)
            row = dict(partner=pname, solvent=solvent, Eint_on=round(e_on, 2), Eint_off=round(e_off, 2),
                       gap=round(e_on - e_off, 2))
            rows.append(row)
            print(json.dumps(row), flush=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(rows, open(a.out, "w"), indent=1)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
