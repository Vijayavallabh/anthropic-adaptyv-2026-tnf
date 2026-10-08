#!/usr/bin/env python3
"""Prepare af2ig inputs from BoltzGen inverse-folded CIFs (ORIGINAL sequences).

ProteinMPNN re-sequencing destroys AF2-confident binding (sequence-specific), so to find more
binders we fold MORE of BoltzGen's OWN sequences (which bind at ~2-4%). Converts each selected
design's complex CIF to PDB (gemmi) and strips the placeholder (0,0,0) side-chain atoms so af2ig
gets binder backbone+sequence + full target. Selection = His433-contacting, novel, cross-reactive-
ish, reasonable length; excludes designs already folded.

Usage: python scripts/prep_af2ig_from_cif.py --out runs/boltz_more_foldin --max-len 100 [--short-only]
"""
from __future__ import annotations
import argparse, csv, glob, os
from pathlib import Path
import gemmi

REPO = Path(__file__).resolve().parents[1]


def num(r, k, d=0.0):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return d


def strip_origin(pdb_text):
    out = []
    for l in pdb_text.splitlines():
        if l[:6] in ("ATOM  ", "HETATM"):
            try:
                x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
            except ValueError:
                continue
            if abs(x) < 1e-3 and abs(y) < 1e-3 and abs(z) < 1e-3:
                continue
        out.append(l)
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyzed", default=str(REPO / "results" / "designs_analyzed.csv"))
    ap.add_argument("--out", default=str(REPO / "runs" / "boltz_more_foldin"))
    ap.add_argument("--min-len", type=int, default=55)
    ap.add_argument("--max-len", type=int, default=100)
    ap.add_argument("--short-only", action="store_true", help="cap length at 82 (binders cluster here)")
    ap.add_argument("--min-cons", type=int, default=2)
    ap.add_argument("--max-div", type=int, default=2)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--exclude-dirs", nargs="*", default=[], help="skip designs already converted here")
    a = ap.parse_args()
    if a.short_only:
        a.max_len = min(a.max_len, 82)

    folded = set()
    for r in csv.DictReader(open(REPO / "results" / "filtered_designs.csv")):
        folded.add(r["design"])
    for dpath in a.exclude_dirs:
        for p in glob.glob(os.path.join(dpath, "*.pdb")):
            folded.add(Path(p).stem)

    rows = list(csv.DictReader(open(a.analyzed)))
    seen_seq = set()
    sel = []
    for r in rows:
        if str(r.get("his433_contact")) != "True":
            continue
        if num(r, "ident_to_ubiquitin", 1) >= 0.6:
            continue
        L = num(r, "binder_len")
        if not (a.min_len <= L <= a.max_len):
            continue
        if num(r, "n_conserved_contacts") < a.min_cons or num(r, "n_divergent_contacts") > a.max_div:
            continue
        if num(r, "frac_C") > 0.1:
            continue
        stem = f"{r['source_spec']}_{r['source_seed']}_{r['design_id']}"
        if stem in folded:
            continue
        seq = r.get("binder_seq", "")
        if seq in seen_seq:
            continue
        seen_seq.add(seq)
        sel.append((stem, r["cif"]))
    sel.sort(key=lambda x: x[0])
    if a.limit:
        sel = sel[:a.limit]

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ok = fail = 0
    for stem, cif in sel:
        try:
            s = gemmi.read_structure(cif)
            s.setup_entities()
            txt = strip_origin(s.make_pdb_string())
            (out / f"{stem}.pdb").write_text(txt)
            ok += 1
        except Exception:
            fail += 1
    print(f"selected {len(sel)} | converted {ok} (fail {fail}) -> {out}")


if __name__ == "__main__":
    main()
