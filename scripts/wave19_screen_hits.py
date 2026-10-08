#!/usr/bin/env python3
"""WAVE 19: rank the MD pH-switch screen (runs/md_w18_screen) of the 157 unscreened switch-passers.

Joins MD occupancy (the selectivity gate) with lineage (switch_passers_w18) and real af2ig binding
(rfd_af2_passers_w18: pae_interaction, plddt_binder). Applies the locked-8 floor gate and ranks by
selectivity, flagging designs from lineages NOT already in the submitted Tier B (the diversity axis).
Re-runnable while the screen is still filling in. Writes results/md_screen_hits_w19.csv.
"""
import csv, glob, json, os

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"; os.chdir(REPO)
SCREEN = "runs/md_w18_screen"
FLOOR_DOCC = 0.2
# lineages already represented in the submitted Tier B (egfr_phsw_21..31)
REPRESENTED = {"egfr_phsw_11", "egfr_phsw_06", "egfr_phsw_03", "egfr_phsw_08", "egfr_phsw_09",
               "egfr_phsw_12", "egfr_phsw_16", "egfr_phsw_07", "egfr_phsw_10"}


def ff(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d


def stem(n):
    return os.path.basename(n)[:-5].replace("__best", "") if n.endswith(".json") else n.replace("__best", "")


lin = {r["design"]: r.get("lineage", "?") for r in csv.DictReader(open("results/switch_passers_w18.csv"))}
af2 = {r["design"]: r for r in csv.DictReader(open("results/rfd_af2_passers_w18.csv"))}

rows = []
for j in sorted(glob.glob(SCREEN + "/*.json")):
    try:
        d = json.load(open(j))
    except Exception:
        continue
    name = stem(d.get("name", os.path.basename(j)))
    hip, hie = ff(d.get("occ_HIP")), ff(d.get("occ_HIE"))
    if hip is None or hie is None:
        continue
    docc = round(hip - hie, 4)
    a = af2.get(name, {})
    rows.append(dict(
        name=name, lineage=lin.get(name, "?"),
        occ_HIP=round(hip, 3), occ_HIE=round(hie, 3), d_occ=docc,
        pae_interaction=ff(a.get("pae_interaction")), plddt_binder=ff(a.get("plddt_binder")),
        passes_floor=(docc >= FLOOR_DOCC and hip > hie),
        clean=(hie < 0.2), strong=(docc >= 0.5),
        new_lineage=(lin.get(name, "?") not in REPRESENTED),
    ))

# rank: floor first, then selectivity (d_occ) with a low-leak (occ_HIE) tiebreak
rows.sort(key=lambda r: (not r["passes_floor"], -r["d_occ"], r["occ_HIE"]))

cols = ["name", "lineage", "d_occ", "occ_HIP", "occ_HIE", "clean", "strong", "new_lineage",
        "passes_floor", "pae_interaction", "plddt_binder"]
with open("results/md_screen_hits_w19.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)

done = len(rows)
flo = [r for r in rows if r["passes_floor"]]
cleanstrong = [r for r in flo if r["clean"] and r["strong"]]
newhits = [r for r in flo if r["new_lineage"]]
print(f"screened={done}/157 | floor-pass={len(flo)} | clean+strong={len(cleanstrong)} | from NEW lineages={len(newhits)}")
print(f"\n{'name':34}{'lin':16}{'d_occ':>6}{'HIE':>6}{'pae':>6}  flags")
for r in flo[:30]:
    flags = " ".join(f for f, on in [("CLEAN", r["clean"]), ("STRONG", r["strong"]), ("NEW", r["new_lineage"])] if on)
    print(f"{r['name']:34}{r['lineage']:16}{r['d_occ']:6.2f}{r['occ_HIE']:6.2f}"
          f"{(r['pae_interaction'] or 0):6.1f}  {flags}")
