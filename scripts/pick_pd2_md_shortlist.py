#!/usr/bin/env python3
"""Pick a diverse shortlist of pd2 (diverse-PD) candidates to MD-validate, so they can dual-confirm
and compete fairly in the family-capped selection. Top-K per baseline family by ddg_pH (pae<max-pae),
canonicalize each complex (binder=A/target=B 1-191) into cand_final2_canon (the dir the MD list reads),
and append to runs/md_final2_list.txt. MD is additive (run_md_campaign skips designs already scored)."""
import csv, os, re, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CANON = REPO/"scripts"/"canon_renumber.py"
PY = sys.executable
K = int(sys.argv[1]) if len(sys.argv) > 1 else 5
MAXPAE = float(sys.argv[2]) if len(sys.argv) > 2 else 11.0


def fam(d):
    m = re.search(r"egfr_dom3_([A-Za-z]+(?:_[sl])?)", d)
    return m.group(1) if m else "other"


def ff(d, k, dv=-9.0):
    try: return float(d[k])
    except (KeyError, ValueError, TypeError): return dv


obj1 = {r["design"]: r for r in csv.DictReader(open(REPO/"results"/"rfd_obj1_switch_pd2.csv")) if not r.get("error")}
pae = {r["design"]: ff(r, "pae_interaction", 99) for r in csv.DictReader(open(REPO/"results"/"rfd_af2_passers_pd2.csv"))}

rows = []
for d, o in obj1.items():
    p = pae.get(d, 99)
    dd = max(ff(o, "ddg_pH"), ff(o, "orig_ddg_pH"))
    if p < MAXPAE and dd > 0:
        rows.append((fam(d), d, dd, p))

by_fam = {}
for f, d, dd, p in sorted(rows, key=lambda x: -x[2]):
    by_fam.setdefault(f, [])
    if len(by_fam[f]) < K:
        by_fam[f].append((d, dd, p))

canon_dir = REPO/"runs"/"cand_final2_canon"; canon_dir.mkdir(parents=True, exist_ok=True)
md_list = REPO/"runs"/"md_final2_list.txt"
existing = set()
if md_list.exists():
    existing = {ln.split("\t")[0] for ln in md_list.read_text().splitlines() if ln.strip()}
added = []
with open(md_list, "a") as fh:
    for f, picks in by_fam.items():
        for d, dd, p in picks:
            if d in existing: continue
            raw = REPO/"runs"/"rfd_passers_pd2"/f"{d}.pdb"
            cpdb = canon_dir/f"{d}.pdb"
            subprocess.run([PY, str(CANON), str(raw), str(cpdb)], check=True, capture_output=True)
            fh.write(f"{d}\t{cpdb}\n")
            added.append((f, d, dd, p))

print(f"shortlist: {len(added)} new pd2 designs across families", {f: len(v) for f, v in by_fam.items()})
for f, d, dd, p in sorted(added, key=lambda x: (x[0], -x[2])):
    print(f"  {f:10s} ddg_pH={dd:5.2f} pae={p:5.1f} {d}")
# paths-only MD list for run_md_campaign
paths = [ln.split("\t")[1] for ln in md_list.read_text().splitlines() if ln.strip()]
(REPO/"runs"/"md_pd2_paths.txt").write_text("\n".join(paths) + "\n")
print(f"MD paths written: {len(paths)} -> runs/md_pd2_paths.txt")
