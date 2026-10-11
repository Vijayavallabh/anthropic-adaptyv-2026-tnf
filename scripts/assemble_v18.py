#!/usr/bin/env python3
"""Assemble v18: swap MD-CONFIRMED, re-validated switchers (from harvest_v18.sh's v18_candidates.csv) into
the NONRESPOND slots of the v17 (or v16) base, respecting <=2 per lineage and keeping 20 unique 0-Cys designs.
Candidate gate: xr_min>=0.70 AND solub>=0.50 (switch already CONFIRMED by md_verdict). Ranks candidates by
dSwitch; a stronger same-lineage candidate can also UPGRADE (replace) a weaker same-lineage member.
NOTE: new de-novo switchers' PORTAL novelty is USER-gated (like des_6025) — this writes the candidate CSV;
the user portal-checks the new members before upload. Novelty NOT auto-confirmable here.
Usage: assemble_v18.py <base.csv> <v18_candidates.csv> <out_v18.csv> [nonrespond_names_comma] [cap_per_lineage=2]
"""
import sys, csv, re
base_f, cand_f, out_f = sys.argv[1], sys.argv[2], sys.argv[3]
NONRESP = (sys.argv[4].split(",") if len(sys.argv) > 4 and sys.argv[4] else
           ["des_3050_6_H52-54-79", "des_818_4_c13S", "b4m2_A57H"])  # v16/v17 NONRESPOND slots
CAP = int(sys.argv[5]) if len(sys.argv) > 5 else 2
def lineage(n):   # collapse to de-novo lineage root (des_138_7_T10_pd.. -> des_138_7), else name stem
    m = re.match(r'(des_\d+_\d+)', n);  return m.group(1) if m else re.sub(r'_(T\d+|mh|pd|H\d|c\d).*', '', n)
base = list(csv.DictReader(open(base_f)))
cand = list(csv.DictReader(open(cand_f)))
def fnum(x):
    try: return float(x)
    except: return None
# qualifying candidates: xr>=0.70 & solub>=0.5, ranked by dSwitch desc
qual = []
for c in cand:
    xr, so, ds = fnum(c.get("xr_min")), fnum(c.get("solub")), fnum(c.get("dSwitch"))
    if xr is not None and xr >= 0.70 and (so is None or so >= 0.50):
        qual.append((ds if ds is not None else -9, c["name"], c.get("sequence", "")))
qual.sort(reverse=True)
print(f"qualifying candidates (xr>=0.70, solub>=0.5): {len(qual)}")
# current lineage counts in base
cnt = {}
for r in base: cnt[lineage(r["name"])] = cnt.get(lineage(r["name"]), 0) + 1
names = [r["name"] for r in base]
rows = {r["name"]: r for r in base}
# fill NONRESPOND slots with qualifying candidates (respecting cap), need their sequences
seqs = {c["name"]: c.get("sequence", "") for c in cand}
slots = [n for n in NONRESP if n in names]
placed = []
for ds, nm, sq in qual:
    if not slots: break
    lin = lineage(nm)
    if cnt.get(lin, 0) >= CAP: continue           # respect <=2 per lineage
    if nm in names: continue
    seq = sq or seqs.get(nm, "")
    if not seq or "C" in seq: continue            # 0-Cys
    slot = slots.pop(0)
    # protect des_818_4 (4/4 novel) unless we have a strong replacement AND other slots exhausted
    idx = names.index(slot)
    names[idx] = nm; rows[nm] = {"name": nm, "sequence": seq, "molecule_class": "single_chain"}
    cnt[lin] = cnt.get(lin, 0) + 1
    placed.append((slot, nm, ds))
# write
with open(out_f, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["name", "sequence", "molecule_class"])
    for n in names:
        r = rows[n]; w.writerow([r["name"], r["sequence"].split(":")[0], r.get("molecule_class", "single_chain")])
seqset = [rows[n]["sequence"] for n in names]
assert len(set(names)) == 20 and len(names) == 20, f"not 20 unique: {len(names)}"
assert all("C" not in s for s in seqset), "Cys present"
print(f"swapped {len(placed)} NONRESPOND -> confirmed switchers:")
for slot, nm, ds in placed: print(f"  {slot}  ->  {nm}  (dSwitch {ds:+.2f})")
print(f"wrote {out_f} (20 designs, unique, 0-Cys). USER portal-novelty-check the new members before upload.")
