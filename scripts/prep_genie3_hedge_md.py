#!/usr/bin/env python3
"""Prep the Genie3-groove switch-HEDGE for two-state release MD.

Context: Genie3 at the CRO-validated Region-I groove is the ONLY recipe that produced cross-reactive TNF
binders in the Anthropic campaign, but af2ig CANNOT confirm its binding (the flat-polar wall; see
lesson-genie3-groove-floors-af2ig-like-all-generators). This probe does NOT claim binding confirmation.
It asks the one untested question that governs the hedge's #1-axis (pH-selectivity) value: does a
multi-His switch on separated conserved groove cations, installed on the groove-PLACED Genie3 pose, give
a CONFIRMED two-state release signal? A positive result strengthens the (honest, user-decided) hedge pitch
on the top-scored axis; a negative result closes it.

The groove-placed pose is the af2in complex (chain A = MPNN-designed binder with real resnames but
BACKBONE-ONLY atoms; chains B,C = the two full-atom TNF protomers). md_release's detect_chains keeps the
2 LARGEST chains, so we must hand it exactly TWO: binder A + ONE TNF protomer as target. TNF is a
homotrimer, so all three conserved Region-I cations (mature K65/K98/K112 = dataset 60/93/107) exist on
EITHER protomer; we keep whichever protomer (B or C) the binder's His anchors best cover, relabel it "B",
and drop the other. This sidesteps all two-protomer merge/topology problems (fused internal termini).
Then:
  - find binder CA near that protomer's cations (CA proxy: the pose has no CB), one anchor per DISTINCT
    conserved cation, up to 3 mutually SEPARATED (>=SEP CA-CA) -> the Baker separated-multi-His form;
  - install HIS resname at those anchors on the backbone binder (md_release's PDBFixer rebuilds every
    sidechain incl His; build() then sets HIE/HIP/GLN). Target protomer stays full-atom.
Emits per-variant PDBs + md_release manifest (design<TAB>his_pos<TAB>anchor<TAB>pdb_abspath).
Usage: envs/md2/bin/python scripts/prep_genie3_hedge_md.py <hedge_top.txt> <af2in_dir> <outdir> [DMAX=10] [SEP=8] [top_designs=12]
"""
import os, sys, glob, itertools
import numpy as np
R = "/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
TOPF, AF2IN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DMAX  = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
SEP   = float(sys.argv[5]) if len(sys.argv) > 5 else 8.0
TOPN  = int(sys.argv[6]) if len(sys.argv) > 6 else 12
os.makedirs(OUT, exist_ok=True)
CONSERVED = {60, 93, 107}   # dataset numbering (= mature K65/K98/K112); conserved human+mouse Region-I cations

def parse_atoms(pdb):
    out = []
    for l in open(pdb):
        if l[:4] != "ATOM": continue
        out.append((l[21], int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l.rstrip("\n")))
    return out

def cations_of(atoms, chain):
    cats = {}
    for ch, rn, res, at, xyz, _ in atoms:
        if ch == chain and ((res == "LYS" and at == "NZ") or (res == "ARG" and at in ("CZ", "NH1", "NH2"))):
            cats.setdefault(rn, []).append(xyz)
    return cats

def best_combo_for(A_ca, cats, seqlen):
    """anchors on ONE protomer: binder pos -> nearest distinct cation <DMAX; pick a separated combo of
    up to 3 preferring conserved cations then tight dist. returns (his_pos_sorted, anchor_descs, ncons) or None."""
    pos_cat = {}
    for rn, (xyz, res) in A_ca.items():
        if res in ("GLY", "PRO", "HIS") or rn <= 1 or rn >= seqlen: continue
        cand = [(min(np.linalg.norm(xyz - c) for c in cs), cat) for cat, cs in cats.items()]
        if not cand: continue
        d, cat = min(cand)
        if d <= DMAX: pos_cat[rn] = (cat, round(float(d), 1), xyz)
    by_cat = {}
    for pos, (cat, dist, xyz) in pos_cat.items():
        if cat not in by_cat or dist < by_cat[cat][1]: by_cat[cat] = (pos, dist, xyz, cat)
    anchors = sorted(by_cat.values(), key=lambda x: x[1])
    if len(anchors) < 2: return None
    def sep_ok(c): return all(np.linalg.norm(a[2] - b[2]) >= SEP for a, b in itertools.combinations(c, 2))
    for k in (3, 2):
        cand = [c for c in itertools.combinations(anchors, k) if sep_ok(c)]
        if cand:
            best = max(cand, key=lambda c: (sum(1 for a in c if a[3] in CONSERVED), -np.mean([a[1] for a in c])))
            his = sorted(a[0] for a in best)
            desc = ",".join(f"{a[0]}->{a[3]}@{a[1]}{'*' if a[3] in CONSERVED else ''}" for a in sorted(best, key=lambda x: x[1]))
            return his, desc, sum(1 for a in best if a[3] in CONSERVED)
    return None

designs = [l.strip() for l in open(TOPF) if l.strip()][:TOPN]
manifest = []
print(f"{'design':20s} {'prot':>4} {'bbmin':>6} {'anchors(pos->cat@dist,*=conserved)':42s} {'his':>12} {'cons'}")
for bb in designs:
    cands = sorted(glob.glob(f"{AF2IN}/{bb}_*.pdb"))
    if not cands: print(f"{bb}: no af2in pose, skip"); continue
    pdb = ([c for c in cands if c.endswith("_1.pdb")] or cands)[0]
    dname = os.path.basename(pdb)[:-4]
    atoms = parse_atoms(pdb)
    A_ca = {rn: (xyz, res) for ch, rn, res, at, xyz, _ in atoms if ch == "A" and at == "CA"}
    if not A_ca: print(f"{dname}: no chain-A CA, skip"); continue
    seqlen = max(A_ca)
    # pick the protomer whose conserved-cation coverage is best
    choice = None
    for prot in ("B", "C"):
        cats = cations_of(atoms, prot)
        if not cats: continue
        bc = best_combo_for(A_ca, cats, seqlen)
        if bc and (choice is None or (bc[2], len(bc[0])) > (choice[2][2], len(choice[2][0]))):
            choice = (prot, cats, bc)
    if choice is None: print(f"{dname}: no 2+ cation anchors on either protomer, skip"); continue
    prot, cats, (his_pos, desc, ncons) = choice
    # groove-placement sanity: binder CA -> kept-protomer heavy min
    tgt_xyz = np.array([xyz for ch, rn, res, at, xyz, _ in atoms if ch == prot])
    a_xyz = np.array([v[0] for v in A_ca.values()])
    bbmin = float(np.min(np.linalg.norm(a_xyz[:, None] - tgt_xyz[None], axis=-1)))
    # write 2-chain pose: binder A (backbone, HIS at anchors) + kept protomer relabeled "B" (full-atom)
    outpdb = os.path.abspath(f"{OUT}/{dname}_mh_{'_'.join(map(str, his_pos))}.pdb")
    with open(outpdb, "w") as fh:
        for ch, rn, res, at, xyz, l in atoms:
            if ch != "A": continue
            rn3 = "HIS" if rn in his_pos else res
            fh.write(l[:17] + f"{rn3:>3s}" + l[20:] + "\n")
        fh.write("TER\n")
        for ch, rn, res, at, xyz, l in atoms:
            if ch != prot: continue
            fh.write(l[:21] + "B" + l[22:] + "\n")   # relabel chain to B, keep original resSeq
        fh.write("TER\nEND\n")
    manifest.append((f"{dname}_mh_{'_'.join(map(str, his_pos))}", ",".join(map(str, his_pos)), str(his_pos[0]), outpdb))
    print(f"{dname:20s} {prot:>4} {bbmin:6.1f} {desc:42s} {'/'.join(map(str,his_pos)):>12} {ncons}/3")

with open(f"{OUT}/manifest.tsv", "w") as fh:
    for did, hp, anc, pdb in manifest:
        fh.write(f"{did}\t{hp}\t{anc}\t{pdb}\n")
print(f"\nwrote {len(manifest)} hedge switch-probe variants + manifest -> {OUT}/manifest.tsv")
print("NOTE: binding is af2ig-UNCONFIRMABLE for Genie3-TNF (wet-lab-only); this probes the SWITCH axis only.")
