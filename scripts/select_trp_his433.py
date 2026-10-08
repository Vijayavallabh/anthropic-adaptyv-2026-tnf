#!/usr/bin/env python3
"""WAVE 27 downstream: from af2ig-folded Trp...His433 binders, SELECT designs where a binder TRP indole
cation-pi-stacks the conserved target His433, on the real (flexible) folded pose, AND that pass the af2ig
binding gate.

af2pred convention: chain A = binder (1..N), chain B = target, renumbered CONTIGUOUSLY after the binder.
Native target residue r -> af2pred (r + offset), offset = (first target-chain residue id) - 310. His433 =
mature 409 -> af2pred (409 + offset). For each binder Trp we measure indole-centroid distance + inter-plane
angle to His433; a design is a cation-pi switcher if a Trp is within CC_SELECT of the His433 ring. Emits a
manifest + the (trpResid:his433Resid) pair in af2pred numbering that qm_switch_fair/md need. Inverts
moonshot_prep_cationpi.py (binder aromatic Trp, target cation His433)."""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

HIS433_MATURE = 409
NATIVE_START = 310
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
INDOLE = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")
CC_SELECT = 6.0


def ring_frame(res, names):
    cs = np.array([res[a].coord for a in names if a in res])
    if len(cs) < 3:
        return None, None
    c = cs.mean(0); _, _, vt = np.linalg.svd(cs - c)
    return c, vt[2]


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
    rid = HIS433_MATURE + offset
    hit = [r for r in tgt_res if r.id[1] == rid and r.resname.startswith("HI")]
    if not hit:
        return None
    hc, hn = ring_frame(hit[0], HIS_RING)
    if hc is None:
        return None
    best = None
    for tr in chains[binder_chain]:
        if tr.resname != "TRP":
            continue
        tc, tn = ring_frame(tr, INDOLE)
        if tc is None:
            continue
        cc = float(np.linalg.norm(tc - hc))
        ang = float(np.degrees(np.arccos(min(1.0, abs(np.dot(tn, hn))))))
        if best is None or cc < best[0]:
            best = (round(cc, 2), round(ang, 1), tr.id[1], rid)
    return dict(offset=offset, best=best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2ig", default="runs/wave27_trp_his433/af2ig")
    ap.add_argument("--out", default="runs/wave27_trp_his433/selected")
    ap.add_argument("--cc-select", type=float, default=CC_SELECT)
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
        cc, ang, trp, his = info["best"]
        if cc > a.cc_select:
            continue
        rows.append(dict(desc=desc, pdb=os.path.abspath(pdb), pae=pae, plddt=plddt, rmsd=rms,
                         trp=trp, his433=his, ring_cc=cc, plane_angle=ang))
    rows.sort(key=lambda r: (r["ring_cc"], r["pae"]))
    with open(os.path.join(a.out, "trp_his433_selected.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["desc", "pae", "plddt", "rmsd", "trp", "his433",
                                           "ring_cc", "plane_angle", "pdb"])
        w.writeheader(); w.writerows(rows)
    json.dump(rows, open(os.path.join(a.out, "trp_his433_selected.json"), "w"), indent=1)
    print(f"Trp...His433 cation-pi switchers selected: {len(rows)}")
    for r in rows[:15]:
        print(f"  {r['desc'][:44]:44s} pae={r['pae']:5.2f} Trp{r['trp']}-His{r['his433']} "
              f"cc={r['ring_cc']} ang={r['plane_angle']}")
    print(f"-> {a.out}/trp_his433_selected.csv")


if __name__ == "__main__":
    main()
