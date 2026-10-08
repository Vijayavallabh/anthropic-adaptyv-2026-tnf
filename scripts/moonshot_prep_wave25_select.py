#!/usr/bin/env python3
"""WAVE 25 selection: from af2ig-folded designs, find MULTI-His switchers -- designs where >=2 binder His
each H-bond a NEUTRAL target acceptor (backbone C=O, Ser/Thr-OH, Asn/Gln amide-O), spread out (no pKa
splitting). These are the candidate Hill n>1 collective switches (premise verified: such contacts raise the
His pKa -- QM + APBS -- yet have a neutral partner -> ECC-robust). Emits per design the His/acceptor pairs for
collective-linkage + QM scoring, restricted to af2ig binders.
"""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

HBOND = 3.5           # His N .. acceptor O
HIS_SEP = 8.0         # min His-His sidechain separation to avoid like-charge pKa splitting
NEUT_SC = {"SER": ["OG"], "THR": ["OG1"], "ASN": ["OD1"], "GLN": ["OE1"]}  # neutral sidechain acceptors
HISN = ("ND1", "NE2")


def load_binders(af2ig_dir):
    keep = {}
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "out.sc")):
        lines = open(f).read().splitlines()
        if not lines:
            continue
        h = lines[0].split()
        try:
            ip, ipl, ir, idc = (h.index("pae_interaction"), h.index("plddt_binder"),
                                h.index("binder_aligned_rmsd"), h.index("description"))
        except ValueError:
            continue
        for ln in lines[1:]:
            p = ln.split()
            if len(p) <= max(ip, ipl, ir, idc):
                continue
            try:
                pae, plddt, rms = float(p[ip]), float(p[ipl]), float(p[ir])
            except ValueError:
                continue
            if pae < 10 and plddt > 80 and rms < 2:
                keep[p[idc]] = (round(pae, 2), round(plddt, 1), round(rms, 2))
    return keep


def find_pdb(af2ig_dir, desc):
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", desc + "*.pdb")):
        return f
    hits = glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", "*" + desc + "*af2pred.pdb"))
    return hits[0] if hits else None


def neutral_acceptors(tgt_res):
    """list of (coord, kind, resid) for neutral acceptors: backbone C=O of every res + neutral sidechain O."""
    acc = []
    for r in tgt_res:
        if "O" in r:
            acc.append((r["O"].coord, "bbCO", r.id[1]))
        for a in NEUT_SC.get(r.resname, []):
            if a in r:
                acc.append((r[a].coord, r.resname, r.id[1]))
    return acc


def ring_centroid(his):
    cs = [his[a].coord for a in ("CG", "ND1", "CD2", "CE1", "NE2") if a in his]
    return np.mean(cs, axis=0) if len(cs) >= 3 else None


def analyze(pdb, binder_chain="A", target_chain="B"):
    model = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    chains = {c.id: [r for r in c if r.id[0] == " "] for c in model}
    if binder_chain not in chains or target_chain not in chains:
        return None
    acc = neutral_acceptors(chains[target_chain])
    if not acc:
        return None
    acc_xyz = np.array([a[0] for a in acc])
    # carboxylate O (to flag salt-bridge-dominated His we DON'T want)
    carbox = []
    for r in chains[target_chain]:
        if r.resname in ("ASP", "GLU"):
            for a in ("OD1", "OD2", "OE1", "OE2"):
                if a in r:
                    carbox.append(r[a].coord)
    carbox = np.array(carbox) if carbox else np.zeros((0, 3))
    switch_his = []
    for h in chains[binder_chain]:
        if not h.resname.startswith("HI"):
            continue
        hn = [h[a].coord for a in HISN if a in h]
        if not hn:
            continue
        d = np.min([np.linalg.norm(acc_xyz - x, axis=1).min() for x in hn])
        if d > HBOND:
            continue
        # which acceptor
        flat = [(np.linalg.norm(a[0] - x), a) for a in acc for x in hn]
        dmin, a_best = min(flat, key=lambda t: t[0])
        # is a carboxylate closer? (salt-bridge-dominated -> not the robust neutral contact)
        dcarb = float(np.linalg.norm(carbox - np.mean(hn, axis=0), axis=1).min()) if len(carbox) else 99.0
        switch_his.append(dict(his=h.id[1], acc_kind=a_best[1], acc_resid=a_best[2],
                               d=round(float(dmin), 2), d_carbox=round(dcarb, 2),
                               cen=ring_centroid(h)))
    # keep distinct, spatially separated His (no pKa splitting)
    kept = []
    for s in sorted(switch_his, key=lambda s: s["d"]):
        if s["cen"] is None:
            continue
        if all(np.linalg.norm(s["cen"] - k["cen"]) > HIS_SEP for k in kept):
            kept.append(s)
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2ig", default="runs/wave25_multihis/af2ig")
    ap.add_argument("--out", default="runs/wave25_multihis/selected")
    ap.add_argument("--min-his", type=int, default=2, help="min spread-out neutral-acceptor His to keep")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    binders = load_binders(a.af2ig)
    print(f"af2ig binders passing gate: {len(binders)}")
    rows = []
    for desc, (pae, plddt, rms) in binders.items():
        pdb = find_pdb(a.af2ig, desc)
        if not pdb:
            continue
        kept = analyze(pdb)
        if not kept or len(kept) < a.min_his:
            continue
        n_robust = sum(1 for s in kept if s["d_carbox"] > s["d"] + 0.5)  # neutral contact closer than any carboxylate
        rows.append(dict(desc=desc, pdb=os.path.abspath(pdb), pae=pae, plddt=plddt, rmsd=rms,
                         n_his=len(kept), n_robust=n_robust,
                         his=";".join(str(s["his"]) for s in kept),
                         acceptors=";".join(f"{s['his']}:{s['acc_kind']}{s['acc_resid']}@{s['d']}" for s in kept)))
    rows.sort(key=lambda r: (-r["n_robust"], -r["n_his"], r["pae"]))
    with open(os.path.join(a.out, "wave25_selected.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["desc", "pae", "plddt", "rmsd", "n_his", "n_robust", "his", "acceptors", "pdb"])
        w.writeheader(); w.writerows(rows)
    json.dump(rows, open(os.path.join(a.out, "wave25_selected.json"), "w"), indent=1)
    print(f"multi-His switchers (>= {a.min_his} spread-out neutral-acceptor His): {len(rows)}")
    for r in rows[:12]:
        print(f"  {r['desc'][:40]:40s} pae={r['pae']:5.2f} n_his={r['n_his']} n_robust={r['n_robust']}  {r['acceptors']}")
    print(f"-> {a.out}/wave25_selected.csv")


if __name__ == "__main__":
    main()
