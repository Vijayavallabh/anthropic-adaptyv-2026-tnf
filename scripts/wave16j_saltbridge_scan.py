#!/usr/bin/env python3
"""WAVE 16j / deep-research M4 prerequisite: confirm the designed His433 salt bridge exists.

Standalone mechanism-geometry check (and the prerequisite for a full isostere ablation): in each
design's complex, find the closest binder Asp/Glu carboxylate-O to a target-His sidechain-N. A
real binder-acidic <-> His433 salt bridge is <=~4.0 A (O-N). Flags designs whose nearest acidic
partner is marginal (>4 A) or contacts a His OTHER than His433 (target ordinal 100 == His100 in
the canon/crop numbering used by obj1_switch). Reads the SAME complex PDBs M1 used
(runs/ens_propka/*.json 'pdb' field). Pure stdlib, CPU-seconds. Validation only; frozen set.
"""
import json, glob, math, os
from collections import defaultdict

HIS433_ORD = 100   # target ordinal of His433 (obj1_switch HIS433_NUM)

def heavy(pdb):
    ch = defaultdict(lambda: defaultdict(list))
    for l in open(pdb):
        if l[:6] in ("ATOM  ", "HETATM"):
            an = l[12:16].strip(); rn = l[17:20].strip(); c = l[21]
            try:
                num = int(l[22:26]); xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                continue
            ch[c][(rn, num)].append((an, xyz))
    return ch

def main():
    rows = []
    for jf in sorted(glob.glob(os.path.join(os.path.dirname(__file__), "..", "runs/ens_propka/egfr_phsw_*.json"))):
        d = json.load(open(jf)); name = d["name"]; pdb = d["pdb"]
        if not os.path.exists(pdb):
            rows.append((name, "PDB-MISSING", None, None, None)); continue
        ch = heavy(pdb)
        chains = sorted(ch, key=lambda c: -sum(len(v) for v in ch[c].values()))
        tc, bc = chains[0], chains[1]
        his_N = {k: [xyz for an, xyz in ats if an in ("ND1", "NE2")]
                 for k, ats in ch[tc].items() if k[0] == "HIS"}
        best = (None, None, 999.0)
        for k, ats in ch[bc].items():
            if k[0] in ("ASP", "GLU"):
                os_ = [xyz for an, xyz in ats if an in ("OD1", "OD2", "OE1", "OE2")]
                for hk, hns in his_N.items():
                    for o in os_:
                        for n in hns:
                            dd = math.dist(o, n)
                            if dd < best[2]:
                                best = (k, hk, dd)
        bres, hres, dmin = best
        is_his433 = (hres is not None and hres[1] == HIS433_ORD)
        rows.append((name,
                     f"{bres[0]}{bres[1]}{bc}" if bres else "NO-ACIDIC",
                     f"{hres[0]}{hres[1]}{tc}" if hres else None,
                     round(dmin, 2) if bres else None, is_his433))
    print(f"{'design':15s}{'binder partner':16s}{'target His':>12}{'O-N (A)':>9}  flag")
    nbridge = nhis433 = 0
    for n, b, h, dd, ishis in rows:
        bridge = dd is not None and dd <= 4.0
        nbridge += bridge; nhis433 += (ishis and bridge)
        flag = ""
        if not bridge: flag = "MARGINAL/none (>4A)"
        elif not ishis: flag = f"contacts {h} (NOT His433)"
        print(f"{n:15s}{b:16s}{str(h):>12}{str(dd):>9}  {flag}")
    print(f"\nsalt bridge <=4.0A: {nbridge}/{len(rows)} | to His433 specifically: {nhis433}/{len(rows)}")

if __name__ == "__main__":
    main()
