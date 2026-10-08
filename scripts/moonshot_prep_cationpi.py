#!/usr/bin/env python3
"""WAVE 24 downstream: from af2ig-folded cation-pi binders, SELECT designs where a binder His imidazole
cation-pi-stacks a conserved target Phe (Phe412 and/or Phe380), on the real (flexible) folded pose.

af2pred convention (learned from the salt-bridge run): chain A = binder (res 1..N), chain B = target,
renumbered CONTIGUOUSLY after the binder. So native target residue r -> af2pred (r + offset), with
offset = (first target-chain residue id) - 310. We map Phe412/Phe380/His433(409) that way, verify the
residue name, and measure, for every binder His, the ring-centroid distance + inter-plane angle to each
conserved Phe. A design is a cation-pi switcher if >=1 His is within CC_SELECT of a conserved Phe face;
DUAL if two distinct His engage the two distinct Phe. Emits a manifest + the --pairs string md_cationpi needs
(hisResid:pheResid in af2pred numbering), restricted to binders that pass the af2ig binding gate.
"""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

NATIVE_PHES = {"Phe412": 412, "Phe380": 380}   # mature; both conserved human==mouse; His433 face
NATIVE_START = 310
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
PHE_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
CC_SELECT = 6.0     # A; generous ring-centroid cutoff for SELECTION (MD then refines with CC_CUT 5.5)


def ring_frame(res, names):
    cs = np.array([res[a].coord for a in names if a in res])
    if len(cs) < 3:
        return None, None
    c = cs.mean(0); _, _, vt = np.linalg.svd(cs - c)
    return c, vt[2]


def load_binders_from_sc(af2ig_dir):
    """return {design_base: (pae, plddt, rmsd)} for rows passing the binding gate."""
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
    tgt = chains[target_chain]
    tgt_res = [r for r in tgt if r.id[0] == " "]
    offset = tgt_res[0].id[1] - NATIVE_START     # af2pred = native + offset
    phe_map = {}
    for nm, nat in NATIVE_PHES.items():
        rid = nat + offset
        hit = [r for r in tgt_res if r.id[1] == rid and r.resname == "PHE"]
        if hit:
            phe_map[nm] = (rid, ring_frame(hit[0], PHE_RING))
    if not phe_map:
        return None
    his = [r for r in chains[binder_chain] if r.resname in ("HIS", "HIP", "HIE", "HID")]
    best = {}   # phe_name -> (ring_cc, angle, his_resid)
    for nm, (rid, (pc, pn)) in phe_map.items():
        if pc is None:
            continue
        for hr in his:
            hc, hn = ring_frame(hr, HIS_RING)
            if hc is None:
                continue
            cc = float(np.linalg.norm(hc - pc))
            ang = float(np.degrees(np.arccos(min(1.0, abs(np.dot(hn, pn))))))
            if nm not in best or cc < best[nm][0]:
                best[nm] = (round(cc, 2), round(ang, 1), hr.id[1], rid)
    return dict(offset=offset, best=best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2ig", default="runs/wave24_cationpi/af2ig")
    ap.add_argument("--out", default="runs/wave24_cationpi/selected")
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
        if not info:
            continue
        engaged = {nm: v for nm, v in info["best"].items() if v[0] <= a.cc_select}
        if not engaged:
            continue
        his_used = {v[2] for v in engaged.values()}
        dual = len(engaged) >= 2 and len(his_used) >= 2
        pairs = ",".join(f"{v[2]}:{v[3]}" for v in engaged.values())
        rows.append(dict(desc=desc, pdb=os.path.abspath(pdb), pae=pae, plddt=plddt, rmsd=rms,
                         n_cationpi=len(engaged), dual=dual, pairs=pairs,
                         phe412_cc=info["best"].get("Phe412", [None])[0],
                         phe412_ang=info["best"].get("Phe412", [None, None])[1] if "Phe412" in info["best"] else None,
                         phe380_cc=info["best"].get("Phe380", [None])[0],
                         phe380_ang=info["best"].get("Phe380", [None, None])[1] if "Phe380" in info["best"] else None))
    rows.sort(key=lambda r: (not r["dual"], r["pae"]))
    with open(os.path.join(a.out, "cationpi_selected.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["desc", "pae", "plddt", "rmsd", "n_cationpi", "dual", "pairs",
                                           "phe412_cc", "phe412_ang", "phe380_cc", "phe380_ang", "pdb"])
        w.writeheader(); w.writerows(rows)
    json.dump(rows, open(os.path.join(a.out, "cationpi_selected.json"), "w"), indent=1)
    ndual = sum(1 for r in rows if r["dual"])
    print(f"cation-pi switchers selected: {len(rows)}  (DUAL={ndual}, single={len(rows) - ndual})")
    for r in rows[:12]:
        print(f"  {r['desc'][:42]:42s} pae={r['pae']:5.2f} npi={r['n_cationpi']} dual={int(r['dual'])} "
              f"Phe412 cc={r['phe412_cc']} Phe380 cc={r['phe380_cc']}  pairs={r['pairs']}")
    print(f"-> {a.out}/cationpi_selected.csv")


if __name__ == "__main__":
    main()
