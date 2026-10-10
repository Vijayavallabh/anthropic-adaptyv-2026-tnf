#!/usr/bin/env python
# Build joint multi-state complexes for cross-react MPNN:
#   chain A = binder (human-docked), B = human target, C = binder (mouse-docked, same numbering), D = mouse target.
# Emit redesigned-residue + symmetry strings: redesign ONLY binder positions within CUT A of either target,
# tie A_i<->C_i, keep the switch His + core fixed. Preserves validated fold/expression; optimizes the
# interface jointly for both orthologs.
import sys, json, glob, os
import gemmi, numpy as np

R = "/home/prachh/v/anthropic-adaptyv-2026-tnf"
CUT = 6.0
OUT = f"{R}/runs/jxr"; os.makedirs(f"{OUT}/pdb", exist_ok=True)
# (design, his_pos, htag, mtag) — switch His position + which Boltz run holds the human/mouse folds.
# des59_H48 = the hiscan cross-react winner (min ipSAE 0.577); joint-optimize its interface further.
DESIGNS = [("des59_H48", 48, "hs_h", "hs_m"),
           ("b1m9_A12H", 12, "rf_h", "rf_m"),
           ("b1m7_A66H", 66, "rf_h", "rf_m"),
           ("b4m2_A57H", 57, "rf_h", "rf_m")]

def load(cif):
    st = gemmi.read_structure(cif); st.setup_entities(); m = st[0]
    chs = list(m)
    # Boltz: chain 0 = binder, chain 1 = target
    return chs[0], chs[1]

def res_atoms(res):
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res])

def one(res):
    t = gemmi.find_tabulated_residue(res.name)
    return t.one_letter_code.upper() if t else "X"

meta = {}
for d, hp, htag, mtag in DESIGNS:
    hc = glob.glob(f"{R}/runs/boltz_{htag}*/out/by_seed/{d}/s0/{d}_model_0.cif")
    mc = glob.glob(f"{R}/runs/boltz_{mtag}*/out/by_seed/{d}/s0/{d}_model_0.cif")
    if not hc or not mc:
        print(f"SKIP {d}: missing h={bool(hc)} m={bool(mc)}"); continue
    hb, ht = load(hc[0]); mb, mt = load(mc[0])
    bseq = "".join(one(r) for r in hb)
    N = len(list(hb))
    if hp > N or bseq[hp-1] != "H":
        print(f"WARN {d}: pos {hp} is '{bseq[hp-1] if hp<=N else '?'}' not H (N={N})")
    # interface binder positions: within CUT of any target atom in EITHER complex
    def iface(binder, tgt):
        tA = np.concatenate([res_atoms(r) for r in tgt])
        s = set()
        for i, r in enumerate(binder, 1):
            bA = res_atoms(r)
            dmin = np.sqrt(((bA[:, None, :] - tA[None, :, :])**2).sum(-1)).min()
            if dmin < CUT: s.add(i)
        return s
    ifh = iface(hb, ht); ifm = iface(mb, mt)
    iface_pos = sorted((ifh | ifm) - {hp})   # redesign interface, NOT the switch His
    # write joint PDB A/B/C/D
    st = gemmi.Structure(); st.add_model(gemmi.Model("1")); M = st[0]
    def addchain(name, src, renum=True):
        ch = gemmi.Chain(name)
        for i, r in enumerate(src, 1):
            nr = gemmi.Residue(); nr.name = r.name
            nr.seqid = gemmi.SeqId(i if renum else r.seqid.num, ' ')
            for a in r:
                na = gemmi.Atom(); na.name = a.name; na.element = a.element
                na.pos = a.pos; na.b_iso = a.b_iso; na.occ = 1.0
                nr.add_atom(na)
            ch.add_residue(nr)
        M.add_chain(ch)
    addchain("A", hb); addchain("B", ht); addchain("C", mb); addchain("D", mt)
    pdb = f"{OUT}/pdb/{d}.pdb"; st.write_pdb(pdb)
    redes = " ".join(f"A{i}" for i in iface_pos) + " " + " ".join(f"C{i}" for i in iface_pos)
    symm = "|".join(f"A{i},C{i}" for i in iface_pos)
    wts = "|".join("0.5,0.5" for _ in iface_pos)
    meta[d] = dict(pdb=pdb, his=hp, N=N, n_iface=len(iface_pos),
                   iface=iface_pos, redesigned=redes, symmetry=symm, weights=wts)
    print(f"{d}: N={N} his={hp} iface={len(iface_pos)} (h={len(ifh)} m={len(ifm)}) -> {pdb}")

json.dump(meta, open(f"{OUT}/meta.json", "w"), indent=1)
print(f"\nwrote {OUT}/meta.json  ({len(meta)} designs)")
