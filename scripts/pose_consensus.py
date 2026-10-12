#!/usr/bin/env python3
"""Cross-folder pose-consensus binding oracle for TNF (bold lever: individual co-fold SCORES are near-chance
for TNF, but if two independent seq-only co-folders [Boltz-2, Chai-1] dock a binder to the SAME epitope in
the SAME orientation, that pose AGREEMENT is a stronger real-interface signal than any single score).

For a design: parse the Boltz and Chai complexes, superpose on the TARGET (TNF) Cα, then measure
  - binder_rmsd : Cα RMSD of the BINDER after target-superposition (low = same placement)
  - iface_jacc  : Jaccard overlap of the TARGET residues contacted (<5 Å) by the binder in each fold
A real, designed interface should drive both folders to the same pose (low rmsd, high jaccard); a scramble
(no real interface) should dock randomly (high rmsd, low jaccard).

Usage: pose_consensus.py <boltz_cif> <chai_cif> <binder_seq> [name]
Prints one JSON line.  Env: run with a Biopython python (envs/uplifting-esmc/bin/python -I)."""
import sys, json
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
from Bio.SeqUtils import seq1

def chains_seqs(model):
    out = []
    for ch in model:
        res = [r for r in ch if r.id[0] == " " and "CA" in r]
        if not res:
            continue
        s = "".join(seq1(r.resname) if r.resname != "UNK" else "X" for r in res)
        out.append((ch.id, res, s))
    return out

def pick_binder_target(model, bseq):
    cs = chains_seqs(model)
    if len(cs) < 2:
        return None, None
    # binder = chain whose seq best matches bseq (by identity fraction over min length)
    def idf(s):
        n = min(len(s), len(bseq));  return sum(1 for i in range(n) if s[i] == bseq[i]) / n if n else 0
    binder = max(cs, key=lambda c: idf(c[2]))
    target = max([c for c in cs if c[0] != binder[0]], key=lambda c: len(c[2]))
    return binder, target

def ca(res_list):
    return np.array([r["CA"].coord for r in res_list])

def iface_target_resids(binder_res, target_res, cut=5.0):
    ba = np.array([a.coord for r in binder_res for a in r])
    hit = set()
    for i, r in enumerate(target_res):
        ta = np.array([a.coord for a in r])
        d = np.sqrt(((ta[:, None, :] - ba[None, :, :]) ** 2).sum(-1))
        if d.min() < cut:
            hit.add(i)  # positional index in target_res (same target seq across folds)
    return hit

def main():
    bcif, ccif, bseq = sys.argv[1], sys.argv[2], sys.argv[3]
    name = sys.argv[4] if len(sys.argv) > 4 else "x"
    P = MMCIFParser(QUIET=True)
    try:
        bm = P.get_structure("b", bcif)[0]; cm = P.get_structure("c", ccif)[0]
    except Exception as e:
        print(json.dumps(dict(name=name, err=f"parse:{e}"))); return
    bB, bT = pick_binder_target(bm, bseq)
    cB, cT = pick_binder_target(cm, bseq)
    if not (bB and cB):
        print(json.dumps(dict(name=name, err="chains"))); return
    # align target residues by sequence position (both are TNF; use the shorter common prefix of matched CA)
    nT = min(len(bT[1]), len(cT[1])); nB = min(len(bB[1]), len(cB[1]))
    tb, tc = ca(bT[1][:nT]), ca(cT[1][:nT])
    # superpose target: rotation that maps chai target onto boltz target
    sup = Superimposer()
    sup.set_atoms([a for r in bT[1][:nT] for a in (r["CA"],)], [a for r in cT[1][:nT] for a in (r["CA"],)])
    rot, tran = sup.rotran
    # apply to chai binder CA, compare to boltz binder CA
    cbCA = ca(cB[1][:nB]) @ rot + tran
    bbCA = ca(bB[1][:nB])
    binder_rmsd = float(np.sqrt(((cbCA - bbCA) ** 2).sum(1).mean()))
    # interface overlap (target positional indices)
    ib = iface_target_resids(bB[1], bT[1][:nT])
    ic = iface_target_resids(cB[1], cT[1][:nT])
    jacc = len(ib & ic) / len(ib | ic) if (ib | ic) else 0.0
    print(json.dumps(dict(name=name, binder_rmsd=round(binder_rmsd, 2), iface_jacc=round(jacc, 3),
                          n_iface_b=len(ib), n_iface_c=len(ic), target_sup_rmsd=round(sup.rms, 2))))

if __name__ == "__main__":
    main()
