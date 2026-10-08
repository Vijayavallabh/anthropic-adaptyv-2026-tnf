#!/usr/bin/env python3
"""WAVE 19: from the MD screen (runs/md_w18_screen) of the 157 unscreened switch-passers, build an
EXTENDED md-selective table = the existing Tier B (md_selective_w18) + new floor-passing hits that ALSO
clear the SAME upstream gates the current set cleared (cross-react delta_pae < 3.0, pKAI-driven switch,
0 free Cys). Every leg except MD occupancy was already computed for all 157 in metrics_w18_full.csv;
the screen supplies the missing selectivity leg. Writes results/md_selective_w19.csv and does NOT touch
the w18 file, so the validated 20 stays as the fallback. The MD-selectivity floor (d_occ>=0.2 AND
occ_HIP>occ_HIE) is applied here AND re-asserted fail-loud in wave19_finalize.
"""
import csv, glob, json, os

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"; os.chdir(REPO)
SCREEN = "runs/md_w18_screen"; FLOOR_DOCC = 0.2


def ff(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def tb(x): return str(x).strip().lower() in ("true", "1", "yes")


def stem(n):
    return os.path.basename(n)[:-5].replace("__best", "") if n.endswith(".json") else n.replace("__best", "")


def tier(docc, hie):
    if docc >= 0.5 and hie < 0.2: return "STRONG"
    if docc >= 0.4 and hie < 0.35: return "GOOD"
    return "FLOOR"


def gated(f):
    """Same upstream gate the current Tier B cleared. MD selectivity is applied separately via the floor;
    the pH switch itself is also confirmed by MD here (these were all static-pH switch-passers)."""
    xr, cys = ff(f.get("xreact_delta_pae")), str(f.get("free_cys"))
    return (xr is not None and xr < 3.0 and tb(f.get("pkai_driven")) and cys in ("0", "0.0"))


full = {r["name"]: r for r in csv.DictReader(open("results/metrics_w18_full.csv"))}
af2 = {r["design"]: r for r in csv.DictReader(open("results/rfd_af2_passers_w18.csv"))}

new, rej = [], {"no_floor": 0, "no_gate": 0, "no_meta": 0}
for j in sorted(glob.glob(SCREEN + "/*.json")):
    try: d = json.load(open(j))
    except Exception: continue
    name = stem(d.get("name", os.path.basename(j)))
    hip, hie = ff(d.get("occ_HIP")), ff(d.get("occ_HIE"))
    if hip is None or hie is None: continue
    docc = round(hip - hie, 4)
    if not (docc >= FLOOR_DOCC and hip > hie):
        rej["no_floor"] += 1; continue
    f = full.get(name)
    if not f:
        rej["no_meta"] += 1; continue
    if not gated(f):
        rej["no_gate"] += 1; continue
    a = af2.get(name, {})
    new.append(dict(name=name, hip=hip, hie=hie, docc=docc, tier=tier(docc, hie),
                    fold=f.get("cluster") or name, lin=f.get("lineage", "?"),
                    pae=ff(a.get("pae_interaction")), tm=ff(f.get("struct_tm")),
                    xreact=ff(f.get("xreact_delta_pae")), novelty=f.get("novelty_level", "")))

cols_out = ["name", "hip", "hie", "docc", "tier", "fold", "lin", "pae", "tm"]
existing = list(csv.DictReader(open("results/md_selective_w18.csv")))
have = {r["name"] for r in existing}
merged = existing[:]
for r in new:
    if r["name"] in have: continue
    merged.append({k: r[k] for k in cols_out})
with open("results/md_selective_w19.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols_out, extrasaction="ignore"); w.writeheader(); w.writerows(merged)

new.sort(key=lambda r: -r["docc"])
screened = len(glob.glob(SCREEN + "/*.json"))
print(f"screened={screened}/157 | new floor+gate-passing hits={len(new)} | "
      f"rejected: floor={rej['no_floor']} gate={rej['no_gate']} no_meta={rej['no_meta']}")
print(f"merged selective total={len(merged)} (existing {len(existing)} + {len(merged)-len(existing)} new)")
print(f"\n{'name':34}{'lineage':16}{'docc':>6}{'hie':>6}{'pae':>6}{'xreact':>7}{'nov':>4}  tier")
for r in new:
    print(f"{r['name']:34}{r['lin']:16}{r['docc']:6.2f}{r['hie']:6.2f}{(r['pae'] or 0):6.1f}"
          f"{(r['xreact'] or 0):7.2f}{str(r['novelty']):>4}  {r['tier']}")
