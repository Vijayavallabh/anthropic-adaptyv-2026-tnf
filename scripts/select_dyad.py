#!/usr/bin/env python3
"""WAVE 28 downstream: from af2ig-folded carboxylate-dyad binders, SELECT designs where a binder GLU carboxyl
forms a buried dyad with the conserved target Asp460, on the folded pose, AND that pass the af2ig binding gate.

af2pred convention: chain A = binder (1..N), chain B = target renumbered CONTIGUOUSLY after the binder.
Native target r -> af2pred (r + offset), offset = (first target-chain residue id) - 310. Asp460 = mature 436
-> af2pred (436 + offset). For each binder Glu we measure min carboxyl O...O to Asp460; a design is a dyad
switcher if a Glu is within OO_SELECT. Emits a manifest + the (gluResid:aspResid) pair for the downstream
two-state QM/APBS/PROPKA. Inverts select_trp_his433.py."""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

ASP460_MATURE = 436
NATIVE_START = 310
ASP_O = ("OD1", "OD2")
GLU_O = ("OE1", "OE2")
OO_SELECT = 3.6


def load_binders_from_sc(af2ig_dir):
    keep = {}
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "out.sc")):
        lines = open(f).read().splitlines()
        if not lines:
            continue
        h = lines[0].split()
        try:
            ip, ipl, ir, idesc = (h.index("pae_interaction"), h.index("plddt_binder"),
                                  h.index("binder_aligned_rmsd"), h.index("description"))
        except ValueError:
            continue
        for ln in lines[1:]:
            p = ln.split()
            if len(p) <= max(ip, ipl, ir, idesc):
                continue
            try:
                pae, plddt, rms = float(p[ip]), float(p[ipl]), float(p[ir])
            except ValueError:
                continue
            if pae < 10 and plddt > 80 and rms < 2:
                keep[p[idesc]] = (round(pae, 2), round(plddt, 1), round(rms, 2))
    return keep


def find_pdb(af2ig_dir, desc):
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", desc + "*.pdb")):
        return f
    hits = glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", "*" + desc + "*af2pred.pdb"))
    return hits[0] if hits else None


def analyze(pdb, binder_chain="A", target_chain="B"):
    model = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    chains = {c.id: c for c in model}
    if binder_chain not in chains or target_chain not in chains:
        return None
    tgt_res = [r for r in chains[target_chain] if r.id[0] == " "]
    offset = tgt_res[0].id[1] - NATIVE_START
    rid = ASP460_MATURE + offset
    hit = [r for r in tgt_res if r.id[1] == rid and r.resname == "ASP"]
    if not hit:
        return None
    asp = hit[0]; ao = [asp[o].coord for o in ASP_O if o in asp]
    if not ao:
        return None
    best = None
    for gr in chains[binder_chain]:
        if gr.resname != "GLU":
            continue
        go = [gr[o].coord for o in GLU_O if o in gr]
        if not go:
            continue
        oo = float(min(np.linalg.norm(x - y) for x in go for y in ao))
        if best is None or oo < best[0]:
            best = (round(oo, 2), gr.id[1], rid)
    return dict(offset=offset, best=best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2ig", default="runs/wave28_dyad/af2ig")
    ap.add_argument("--out", default="runs/wave28_dyad/selected")
    ap.add_argument("--oo-select", type=float, default=OO_SELECT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    binders = load_binders_from_sc(a.af2ig)
    print(f"af2ig binders passing gate (pae<10 & plddt>80 & rmsd<2): {len(binders)}")
    rows = []
    for desc, (pae, plddt, rms) in binders.items():
        pdb = find_pdb(a.af2ig, desc)
        if not pdb:
            continue
        info = analyze(pdb)
        if not info or info["best"] is None:
            continue
        oo, glu, asp = info["best"]
        if oo > a.oo_select:
            continue
        rows.append(dict(desc=desc, pdb=os.path.abspath(pdb), pae=pae, plddt=plddt, rmsd=rms,
                         glu=glu, asp460=asp, dyad_oo=oo))
    rows.sort(key=lambda r: (r["dyad_oo"], r["pae"]))
    with open(os.path.join(a.out, "dyad_selected.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["desc", "pae", "plddt", "rmsd", "glu", "asp460", "dyad_oo", "pdb"])
        w.writeheader(); w.writerows(rows)
    json.dump(rows, open(os.path.join(a.out, "dyad_selected.json"), "w"), indent=1)
    print(f"carboxylate-dyad switchers selected: {len(rows)}")
    for r in rows[:15]:
        print(f"  {r['desc'][:44]:44s} pae={r['pae']:5.2f} Glu{r['glu']}-Asp{r['asp460']} O..O={r['dyad_oo']}")
    print(f"-> {a.out}/dyad_selected.csv")


if __name__ == "__main__":
    main()
