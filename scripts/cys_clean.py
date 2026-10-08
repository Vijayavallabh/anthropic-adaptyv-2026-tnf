#!/usr/bin/env python
"""WAVE 16c developability cleaning: remove UNPAIRED (free) cysteines from a submitted binder
without changing function, then RE-VALIDATE and retain only strict non-regressions.

WHY: a free Cys is the dominant expression/aggregation liability for de-novo mini-binders
(Adaptyv's own developability guidance; Bennett 2023, Cao 2022). Expression GATES every scored
objective (a design that doesn't express scores zero on pH-selectivity / cross-reactivity /
affinity). The project already PREFERS 0-Cys designs (the WAVE-12 Option-G+ swap-ins were chosen
0-Cys); this applies the same standard to the retained floor designs that still carried free Cys.

METHOD (per design, zero data leakage — uses only the submitted sequence + the fixed 6ARU target):
  1. BUILD: take the design's engineered binder:target complex; mutate the free Cys -> Ser (surface)
     or Ala (buried), keeping backbone + CB so af2-initial-guess rebuilds the sidechain. Freeze
     nothing else. (this script)
  2. RE-FOLD: run_af2ig.py --mode fast on the mutant + the original (same run = fair control);
     binding is PRESERVED iff pae_interaction / binder_aligned_rmsd / plddt_binder are within noise.
  3. RE-VALIDATE #1 (switch): obj1_pkai.score_one + obj1_switch.ddg_ph_from_pdb on the mutant
     (renumber the af2pred chain B to canonical so His433 == B100 first; offset = min(chain-B resnum)-1).
     Require His433 still driven (pKAI dpKa > 0.3) and the PROPKA His433 term ~unchanged.
  4. RE-VALIDATE #2 (cross-reactivity): xreact_mechanism.py on the mutant; require core-4A severity
     still 0 (same divergent contacts, all conservative).
  5. RETAIN the mutant ONLY if it is 0-Cys AND binding + switch + x-react are all non-regressed;
     otherwise keep the original (the free Cys is tolerated). Preserve the pre-clean submission as a
     checkpoint (results/submission_track3_v2_pre_cysclean.csv).

Evidence written to results/wave16c_cysclean.csv (orig vs mutant for each metric + verdict).

Usage (build step):
  python scripts/cys_clean.py --struct runs/eng_foldin_pd3/<design>.pdb --binder A \
      --cys 62 75 --repl S --name egfr_phsw_04 --out runs/cys_clean_in04
then: envs/md/bin/python scripts/run_af2ig.py --pdbdir runs/cys_clean_in04 --out ... --gpus 0-2 --mode fast
"""
from __future__ import annotations
import argparse, os, shutil

KEEP = {"N", "CA", "C", "O", "CB"}  # af2-initial-guess rebuilds the rest from sequence
THREE = {"S": "SER", "A": "ALA", "V": "VAL", "T": "THR", "L": "LEU", "M": "MET", "I": "ILE"}


def build_mutant(struct: str, binder: str, cys_pos: list[int], repl: str, dst: str) -> None:
    """Write a complex with the binder's `cys_pos` residues mutated to `repl` (1-letter)."""
    newresn = THREE[repl]
    # verify the targets really are CYS (fail loud on a numbering mismatch)
    is_cys = {}
    for ln in open(struct):
        if ln.startswith("ATOM") and ln[21] == binder and ln[12:16].strip() == "CA":
            is_cys[int(ln[22:26])] = ln[17:20].strip()
    for p in cys_pos:
        if is_cys.get(p) != "CYS":
            raise SystemExit(f"[FATAL] binder {binder}{p} is {is_cys.get(p)}, not CYS (numbering?) — abort")
    with open(dst, "w") as o:
        for ln in open(struct):
            if ln.startswith(("ATOM", "HETATM")):
                try:                              # match the verify loop's int() parse (handles e.g. negative resnums);
                    rn = int(ln[22:26])           # using .isdigit() here would silently skip a verified Cys at a non-digit resnum
                except ValueError:
                    rn = None
                if ln[21] == binder and rn in cys_pos:
                    if ln[12:16].strip() not in KEEP:
                        continue
                    ln = ln[:17] + f"{newresn:>3}" + ln[20:]
                o.write(ln)
            else:
                o.write(ln)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a free-Cys-removed mutant complex for af2ig re-fold")
    ap.add_argument("--struct", required=True, help="engineered binder:target complex PDB (submit==folded seq)")
    ap.add_argument("--binder", default="A")
    ap.add_argument("--cys", type=int, nargs="+", required=True, help="binder residue numbers of the free Cys")
    ap.add_argument("--repl", default="S", choices=list(THREE), help="replacement AA (S surface / A buried)")
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", required=True, help="output dir (also receives the unmutated original control)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tag = "".join(f"C{p}{a.repl}" for p in a.cys)
    build_mutant(a.struct, a.binder, a.cys, a.repl, f"{a.out}/{a.name}_{tag}.pdb")
    shutil.copy(a.struct, f"{a.out}/{a.name}_orig.pdb")
    print(f"wrote {a.name}_{tag}.pdb + {a.name}_orig.pdb -> {a.out}")


if __name__ == "__main__":
    main()
