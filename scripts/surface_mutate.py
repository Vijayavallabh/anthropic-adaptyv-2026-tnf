#!/usr/bin/env python
"""WAVE 16f surface-developability fix: reduce a design's solvent-exposed hydrophobic patch
(the structure-based SAP outlier axis, scripts/sap_scan.py) WITHOUT changing function, then
RE-VALIDATE and retain only strict non-regressions. Generalizes scripts/cys_clean.py (WAVE 16d,
Cys-only) to an arbitrary set of position->residue point mutations.

WHY: a large solvent-exposed, spatially-clustered hydrophobic patch is a wet-lab aggregation /
low-expression liability (Chennamsetty 2009 SAP; Lauer 2012 developability index). Expression
GATES every scored objective. sap_scan.py flags designs whose largest hydrophobic patch, as a
%-of-surface, exceeds the de-novo-class p95; this tool removes the most-exposed, NON-interface,
far-from-switch patch residues by a conservative polar substitution (default Gln).

METHOD (per design, zero data leakage — submitted sequence + fixed target only):
  1. Pick targets with sap_scan + a per-residue relSASA / interface / distance-to-His433 scan: the
     MOST solvent-exposed hydrophobics that are NOT at the binding interface and FAR (>~12 A) from
     the His433 switch (so binding and the switch are untouched). Prefer aromatics (Tyr/Phe/Trp).
  2. BUILD (this script): mutate those positions -> Gln (or --to), keeping backbone + CB so
     af2-initial-guess rebuilds the sidechain. Keyed off SEQUENTIAL CA position (numbering-safe),
     with fail-loud verification that each source residue is the expected type.
  3. RE-FOLD: run_af2ig.py --mode fast on the mutant(s) + the unmutated original (same run = fair
     control). Binding preserved iff pae_interaction / binder_aligned_rmsd / plddt_binder within noise.
  4. RE-VALIDATE: SAP patch (must drop below class p95); His433 switch pKAI AND PROPKA (canon-renumber
     so His433==B100 first); cross-reactivity mechanism severity (still 0); PRODIGY affinity; novelty.
  5. RETAIN the mutant ONLY if it is 0-Cys, the patch is reduced, AND binding + switch + x-react +
     novelty are all non-regressed; else keep the original (the patch is tolerated / flagged).

First applied: egfr_phsw_15 v3 (seqpos 30,48,51 Leu/Leu/Tyr -> Gln) — SAP patch 19.9%->8.5%
(below pool mean), switch pKAI 2.73->3.22 / PROPKA 2.12->2.03, binding pae 6.68->6.33 /
PRODIGY -10.7->-11.7, xreact sev 0, de novo 60.7% — retained (RESULTS_LOG WAVE 16f).

Usage (build step):
  python scripts/surface_mutate.py --struct <complex.pdb> --binder A --pos 30 48 51 \
      --expect L L Y --to Q --name egfr_phsw_15_v3 --out runs/sap15_fix_in
then: envs/md/bin/python scripts/run_af2ig.py --pdbdir runs/sap15_fix_in --out ... --gpus 0-2 --mode fast
"""
from __future__ import annotations
import argparse, os

KEEP = {"N", "CA", "C", "O", "CB"}  # af2-initial-guess rebuilds the rest from the sequence
THREE = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS',
         'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
ONE = {v: k for k, v in THREE.items()}


def binder_pos_index(struct: str, binder: str) -> dict[int, tuple[int, str]]:
    """sequential CA index (1..N) -> (pdb_resnum, 1-letter) for the binder chain."""
    out, i = {}, 0
    for ln in open(struct):
        if ln.startswith("ATOM") and ln[21] == binder and ln[12:16].strip() == "CA":
            i += 1
            out[i] = (int(ln[22:26]), ONE.get(ln[17:20].strip(), "X"))
    return out


def build_mutant(struct: str, binder: str, pos: list[int], expect: list[str], to: list[str], dst: str) -> None:
    idx = binder_pos_index(struct, binder)
    tos = to * len(pos) if len(to) == 1 else to   # one value broadcasts; else per-position
    if len(tos) != len(pos):
        raise SystemExit(f"[FATAL] --to must be 1 value or {len(pos)} values (got {len(to)})")
    for t in tos:
        if t not in THREE:
            raise SystemExit(f"[FATAL] --to {t} is not a standard 1-letter residue")
    resnum_to_new = {}
    for p, e, t in zip(pos, expect, tos):
        if p not in idx:
            raise SystemExit(f"[FATAL] binder seqpos {p} out of range (len {len(idx)})")
        rn, aa = idx[p]
        if aa != e:
            raise SystemExit(f"[FATAL] binder seqpos {p} is {aa}, expected {e} (numbering?) — abort")
        resnum_to_new[rn] = THREE[t]
    with open(dst, "w") as o:
        for ln in open(struct):
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    rn = int(ln[22:26])
                except ValueError:
                    rn = None
                if ln[21] == binder and rn in resnum_to_new:
                    if ln[12:16].strip() not in KEEP:
                        continue
                    ln = ln[:17] + f"{resnum_to_new[rn]:>3}" + ln[20:]
                o.write(ln)
            else:
                o.write(ln)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a surface-hydrophobic-reduced mutant complex for af2ig re-fold")
    ap.add_argument("--struct", required=True, help="engineered binder:target complex PDB (submit==folded seq)")
    ap.add_argument("--binder", default="A")
    ap.add_argument("--pos", type=int, nargs="+", required=True, help="binder SEQUENTIAL positions (1..N, CA order)")
    ap.add_argument("--expect", nargs="+", required=True, help="expected source 1-letter per --pos (fail-loud check)")
    ap.add_argument("--to", nargs="+", default=["Q"], help="replacement 1-letter: one value (all) or one per --pos (default Q = Gln)")
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", required=True, help="output dir (also receives the unmutated original control)")
    a = ap.parse_args()
    if len(a.pos) != len(a.expect):
        ap.error("--pos and --expect must have the same length")
    os.makedirs(a.out, exist_ok=True)
    build_mutant(a.struct, a.binder, a.pos, a.expect, a.to, f"{a.out}/{a.name}.pdb")
    with open(f"{a.out}/{a.name.rsplit('_',1)[0]}_orig.pdb", "w") as o:
        o.write(open(a.struct).read())
    tos = a.to * len(a.pos) if len(a.to) == 1 else a.to
    print(f"wrote {a.name}.pdb (+ orig control) -> {a.out}  [{'/'.join(f'{e}{p}{t}' for p,e,t in zip(a.pos,a.expect,tos))}]")


if __name__ == "__main__":
    main()
