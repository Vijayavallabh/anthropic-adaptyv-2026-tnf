#!/usr/bin/env python3
"""WAVE 29 PivotLid topology selector: among af2ig binders (pae<10), find designs that are BILOBED with the
epitope+staple on a SMALL lid lobe spatially separable from a base lobe that does NOT contact the target.

For each binder: spatially 2-cluster the binder CAs (KMeans); the lid = cluster containing the staple His
(the LOCKED motif His from staging, fixed_residues_multi.json; falls back to the His nearest the conserved
target Asp460 with a warning). Require (a) lid small (<=LID_MAX res), (b) lid holds most epitope contacts,
(c) base contacts the target in <=BASE_EPI residues (base is non-binding), (d) the two lobes are spatially
distinct (<=INTER_MAX inter-lobe CA contacts = a hinge, not a merged domain), and (e) the His...Asp460 staple
is actually FORMED (<=STAPLE_MAX A). Emits lid residue ids per candidate for md_pH_swing. A PivotLid that can
PEEL requires this topology; a monolithic binder (one lobe, distributed epitope, merged lobes, or a broken
staple) is rejected."""
import argparse, csv, glob, json, os, re, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

ASP460_MATURE = 436
NATIVE_START = 310
LID_MAX = 40
BASE_EPI = 3
INTER_MAX = 20     # inter-lobe CA pairs <8A: a thin hinge has few, a merged domain many
STAPLE_MAX = 4.0   # His-N...Asp460-O: the designed pH staple must actually be formed


def load_binders(af2ig_dir):
    keep = {}
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "out.sc")):
        L = open(f).read().splitlines()
        if not L:
            continue
        h = L[0].split()
        try:
            ip, ipl, ir, idc = (h.index("pae_interaction"), h.index("plddt_binder"),
                                h.index("binder_aligned_rmsd"), h.index("description"))
        except ValueError:
            continue
        for ln in L[1:]:
            p = ln.split()
            if len(p) <= max(ip, ipl, ir, idc):
                continue
            try:
                pae, pl, rm = float(p[ip]), float(p[ipl]), float(p[ir])
            except ValueError:
                continue
            if pae < 10 and pl > 80 and rm < 2:
                keep[p[idc]] = round(pae, 2)
    return keep


def _core(name):
    """Strip the leading gpuN_ shard prefix and any trailing _<seq>_af2pred so a staged name and an af2pred
    desc reduce to the same core key (the af2ig shard gpu differs from the gen gpu in the staged name)."""
    n = re.sub(r"^gpu\d+_", "", name)
    n = re.sub(r"_\d+_af2pred$", "", n)  # strip MPNN seq idx + af2pred (desc only; no-op on staged names)
    return n


def load_locked(staged_dir):
    """core-name -> locked staple His resnum (int), from fixed_residues_multi.json ({staged_pdb: 'A<num>'})."""
    fx = os.path.join(staged_dir, "fixed_residues_multi.json")
    out = {}
    if not os.path.exists(fx):
        return out
    for pdb, spec in json.load(open(fx)).items():
        m = re.match(r"[A-Za-z](\d+)", str(spec))
        if m:
            out[_core(os.path.basename(pdb).replace(".pdb", ""))] = int(m.group(1))
    return out


def find_pdb(af2ig_dir, desc):
    for f in glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", desc + "*.pdb")):
        return f
    hits = glob.glob(os.path.join(af2ig_dir, "gpu*", "pdbs", "*" + desc + "*af2pred.pdb"))
    return hits[0] if hits else None


def kmeans2(X, iters=50):
    c = X[[0, len(X) - 1]].astype(float)  # seed with the two termini (lid is usually terminal)
    lab = np.zeros(len(X), dtype=int)
    for _ in range(iters):
        d = np.linalg.norm(X[:, None, :] - c[None, :, :], axis=2)
        lab = d.argmin(1)
        newc = np.array([X[lab == k].mean(0) if (lab == k).any() else c[k] for k in range(2)])
        if np.allclose(newc, c):
            break
        c = newc
    return lab


def analyze(pdb, locked_his=None, binder_chain="A", target_chain="B"):
    m = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    chains = {c.id: c for c in m}
    if binder_chain not in chains or target_chain not in chains:
        return None
    A = [r for r in chains[binder_chain] if r.id[0] == " "]
    B = [r for r in chains[target_chain] if r.id[0] == " "]
    off = B[0].id[1] - NATIVE_START
    asp = [r for r in B if r.id[1] == ASP460_MATURE + off and r.resname == "ASP"]
    if not asp:
        return None
    ao = [asp[0][o].coord for o in ("OD1", "OD2") if o in asp[0]]
    if not ao:
        return None

    def his_asp_d(r):
        hn = [r[n].coord for n in ("ND1", "NE2") if n in r]
        return min((np.linalg.norm(x - y) for x in hn for y in ao), default=None)

    # nearest-Asp His (reference), and the LOCKED motif His (preferred staple)
    nearest = None; nd = 1e9
    for r in A:
        if r.resname.startswith("HI"):
            d = his_asp_d(r)
            if d is not None and d < nd:
                nd = d; nearest = r
    staple = None; warn = ""
    if locked_his is not None:
        cand = [r for r in A if r.id[1] == locked_his and r.resname.startswith("HI")]
        if cand:
            staple = cand[0]
            if nearest is not None and nearest.id[1] != staple.id[1]:
                warn = f"locked His{staple.id[1]} != nearest-Asp His{nearest.id[1]}"
        else:
            warn = f"locked pos {locked_his} not a His in af2pred; using nearest-Asp"
    if staple is None:
        staple = nearest
    if staple is None:
        return None
    staple_d = his_asp_d(staple)

    ca = np.array([r["CA"].coord for r in A if "CA" in r])
    ca_ids = np.array([r.id[1] for r in A if "CA" in r])
    lab = kmeans2(ca)
    # lid cluster = the one containing the staple His
    sidx = int(np.where(ca_ids == staple.id[1])[0][0])
    mask = (lab == lab[sidx])
    lid_ids = ca_ids[mask].tolist()
    base_ids = ca_ids[~mask].tolist()
    lid_ca, base_ca = ca[mask], ca[~mask]
    tgt = np.array([a.coord for r in B for a in r])
    epi = [r.id[1] for r in A if min(np.linalg.norm(a.coord - tgt, axis=1).min() for a in r) < 4.5]
    epi_in_lid = sum(1 for e in epi if e in lid_ids)
    base_epi = sum(1 for e in epi if e in base_ids)
    # inter-lobe CA contacts (hinge-like if few): count CA pairs <8A across lobes
    inter = int((np.linalg.norm(lid_ca[:, None, :] - base_ca[None, :, :], axis=2) < 8.0).sum())
    bilobed = (len(lid_ids) <= LID_MAX and len(base_ids) >= 20 and epi_in_lid >= max(3, 0.7 * len(epi))
               and base_epi <= BASE_EPI and inter <= INTER_MAX and staple_d <= STAPLE_MAX)
    return dict(staple=int(staple.id[1]), staple_d=round(staple_d, 2), lid_n=len(lid_ids),
                base_n=len(base_ids), epi_n=len(epi), epi_in_lid=epi_in_lid, base_epi=base_epi,
                inter_lobe=inter, lid_ids=sorted(lid_ids), bilobed=bool(bilobed), warn=warn)


def _rng(ids):
    """compact contiguous-range string for md_pH_swing --lid-resids."""
    ids = sorted(ids); out = []; s = p = ids[0]
    for x in ids[1:]:
        if x == p + 1:
            p = x
        else:
            out.append(f"{s}-{p}" if s != p else f"{s}"); s = p = x
    out.append(f"{s}-{p}" if s != p else f"{s}")
    return ",".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2ig", default="runs/wave29_pivotlid2/af2ig")
    ap.add_argument("--staged", default="runs/wave29_pivotlid2/staged")
    ap.add_argument("--out", default="runs/wave29_pivotlid2/selected")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    locked = load_locked(a.staged)
    print(f"locked staple-His map entries: {len(locked)}")
    binders = load_binders(a.af2ig)
    print(f"af2ig binders (pae<10): {len(binders)}")
    rows = []
    for desc, pae in sorted(binders.items(), key=lambda kv: kv[1]):
        pdb = find_pdb(a.af2ig, desc)
        if not pdb:
            continue
        try:
            info = analyze(pdb, locked_his=locked.get(_core(desc)))
        except Exception as e:
            print(f"  FAIL {desc[:40]}: {repr(e)[:50]}"); continue
        if not info:
            continue
        info.update(desc=desc, pae=pae, pdb=os.path.abspath(pdb),
                    lid_resids=_rng(info["lid_ids"]) if info["lid_ids"] else "")
        rows.append(info)
        tag = "BILOBED" if info["bilobed"] else ""
        wtag = f" !{info['warn']}" if info["warn"] else ""
        print(f"  {desc[:38]:38s} pae={pae:5.2f} staple His{info['staple']}(d{info['staple_d']}) "
              f"lid={info['lid_n']} base={info['base_n']} epi_lid={info['epi_in_lid']}/{info['epi_n']} "
              f"base_epi={info['base_epi']} inter={info['inter_lobe']} {tag}{wtag}")
    bilobed = [r for r in rows if r["bilobed"]]
    for r in rows:
        r.pop("lid_ids", None)
    with open(os.path.join(a.out, "pivotlid_selected.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["desc", "pae", "staple", "staple_d", "lid_n", "base_n", "epi_n",
                                           "epi_in_lid", "base_epi", "inter_lobe", "bilobed", "warn",
                                           "lid_resids", "pdb"])
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.out, "pivotlid_bilobed.json"), "w") as fh:
        json.dump(bilobed, fh, indent=1)
    print(f"\nBILOBED PivotLid candidates (small lid + non-binding base + hinge + formed staple): "
          f"{len(bilobed)} / {len(rows)}")
    for r in bilobed[:8]:
        print(f"  {r['desc'][:40]:40s} pae={r['pae']} staple=His{r['staple']} lid={r['lid_resids']}")
    print(f"-> {a.out}/pivotlid_selected.csv")


if __name__ == "__main__":
    main()
