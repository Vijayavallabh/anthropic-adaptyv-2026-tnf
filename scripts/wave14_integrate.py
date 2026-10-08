#!/usr/bin/env python3
"""WAVE 14 comprehensive integration consumer (hardened per 2026-10-02 /code-review).

Fuses the new pH-aware MM-GBSA binding axis (md_pH_bind) into the G+ manifest and adds:
  * 4-leg pH-selectivity sign-consensus (PROPKA Wyman + pKAI + MD-occupancy + MD-binding-dG),
    with a pKAI-reverse HARD VETO (blood-ON designs are rejected, cf. live #20; mirrors
    select_final_v2.py:287-288) and an SD-aware dd_bind verdict (a reverse binding-switch is only
    flagged when the per-replica CI excludes 0 -- a bare sign on a noisy single-traj MM-GBSA point
    estimate would hand the newest/noisiest axis veto power via sub-kcal jitter).
  * pH-window fold-change from the Wyman legs only (exp(ddG_pH/RT), RT=0.616 kcal/mol @310 K).
    We do NOT exp() the MM-GBSA interface dG (it is NOT a pH-linkage free energy -> astronomically
    inflated); dd_bind is reported as its value+sign only.
  * N-glycan epitope-proximity to the DIII sequon Asn SIDECHAIN (ND2/CG/OD1/CB = glycan-attachment
    point), with a robust heavy-atom filter (digit-leading PDB-v2 H names like 1HB are hydrogens).
Legs track PRESENCE: a missing upstream cell is None (not a measured failing leg); quad_present
reports how many legs were actually measured. Outputs an ADDITIVE wave14 manifest (never overwrites
wave13) + a decision summary.

Target-fact note (single-source): the DIII crop numbering (His433 = target chain ordinal 100) is
defined in data/targets/6ARU_domainIII_map.json; we assert ordinal-100 is a His (fail-loud) rather
than trust a hard-coded residue count.
"""
import csv, json, os, glob, math, sys
import numpy as np

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
MAN = f"{REPO}/results/submission_manifest_v2_ALT_optGplus_wave13.csv"
SCREEN = f"{REPO}/runs/md_bind_gplus/out"
DEEP = f"{REPO}/runs/md_bind_gplus/deep"
CPLX = f"{REPO}/runs/md_bind_gplus/in"
OUT_CSV = f"{REPO}/results/wave14_phbind.csv"
OUT_MAN = f"{REPO}/results/submission_manifest_v2_wave14.csv"
RT = 0.616  # kcal/mol at 310 K

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
       'HIS': 'H', 'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
       'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
ASN_SIDECHAIN = {"CB", "CG", "OD1", "ND2"}  # glycan projects from ND2 (sidechain), not backbone


def warn(*a):
    print("WARN:", *a, file=sys.stderr, flush=True)


def present_float(x):
    """Return float(x) if the cell is a parseable number, else None (missing/blank/unparseable)."""
    if x is None:
        return None
    s = str(x).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def is_hydrogen(line, name):
    elem = line[76:78].strip() if len(line) >= 78 else ""
    if elem:
        return elem == "H"
    return name.lstrip("0123456789")[:1] == "H"   # 1HB / 2HD1 / HG etc. are all hydrogens


def parse_heavy(pdb):
    """chain -> ordered {(resi,resn): [(atomname,x,y,z)]}, hydrogens excluded robustly."""
    ch = {}
    for l in open(pdb):
        if not l.startswith("ATOM"):
            continue
        name = l[12:16].strip()
        if is_hydrogen(l, name):
            continue
        ch.setdefault(l[21], {}).setdefault((l[22:26].strip(), l[17:20].strip()), []).append(
            (name, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return ch


def select_chains(ch):
    """target = largest chain, binder = 2nd-largest (matches md_pH_switch.detect_chains)."""
    order = sorted(ch, key=lambda c: -len(ch[c]))
    target = order[0]
    binder = order[1] if len(order) > 1 else None
    return target, binder


def glycan_prox(pdb):
    ch = parse_heavy(pdb)
    tc, bc = select_chains(ch)
    tres = list(ch[tc].values())
    treskeys = list(ch[tc].keys())
    # fail-loud single-source sanity: His433 == target ordinal 100 (6ARU_domainIII_map.json)
    if len(treskeys) < 100 or not treskeys[99][1].startswith("HI"):
        warn(f"{os.path.basename(pdb)}: target ordinal-100 is not His "
             f"({treskeys[99][1] if len(treskeys) >= 100 else 'n/a'}); numbering desync — skipping glycan")
        return None, []
    tseq = ''.join(AA3.get(rn, 'X') for (_, rn) in treskeys)
    sequ = [i for i in range(len(tseq) - 2) if tseq[i] == 'N' and tseq[i + 1] != 'P' and tseq[i + 2] in 'ST']
    if bc is None:
        return None, []
    batoms = np.array([[a[1], a[2], a[3]] for res in ch[bc].values() for a in res])
    if len(batoms) == 0:
        return None, []
    best = 99.0
    contacted = []
    for si in sequ:
        sc = np.array([[a[1], a[2], a[3]] for a in tres[si] if a[0] in ASN_SIDECHAIN])
        if len(sc) == 0:   # glycine-like / missing sidechain atoms -> skip this sequon
            continue
        d = float(np.min(np.linalg.norm(batoms[:, None, :] - sc[None, :, :], axis=-1)))
        if d < 4.5:
            contacted.append(si + 1)
        best = min(best, d)
    return round(best, 1), contacted


def per_rep_ddbind(d):
    """per-replica dd_bind = dG_bind(HIP,rep) - dG_bind(HIE,rep); returns (mean, sd, n, has_pair)."""
    byrep = {}
    for s in d.get("states", []):
        if s.get("dG_bind_kcal") is None:
            continue
        byrep.setdefault(s["replica"], {})[s["variant"]] = s["dG_bind_kcal"]
    vals = [vv["HIP"] - vv["HIE"] for vv in byrep.values() if "HIP" in vv and "HIE" in vv]
    if not vals:
        return None, None, 0, d.get("has_sensor_pair", True)
    m = sum(vals) / len(vals)
    sd = (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5
    return m, sd, len(vals), d.get("has_sensor_pair", True)


def fold(ddg):
    if ddg is None:
        return ""
    try:
        return round(math.exp(ddg / RT), 1)
    except OverflowError:
        return ""


def main():
    rows = list(csv.DictReader(open(MAN)))
    # load screen (2-rep) first, then overlay deep (4-rep) -- deep wins for the verified designs
    dd = {}
    for d_dir, label in ((SCREEN, "screen"), (DEEP, "deep")):
        for j in glob.glob(f"{d_dir}/*.json"):
            try:
                d = json.load(open(j))
            except (json.JSONDecodeError, OSError) as e:
                warn(f"could not parse binding json {j}: {e!r}")   # fail-loud, not silent skip
                continue
            dd[os.path.splitext(os.path.basename(j))[0]] = d

    gp = {}
    for pdb in glob.glob(f"{CPLX}/*.pdb"):
        nm = os.path.basename(pdb)[:-4]
        gp[nm] = glycan_prox(pdb)   # exceptions here are real bugs -> let them surface (no mask)

    out = []
    for r in rows:
        nm = r["name"]
        d = dd.get(nm)
        rec = dict(r)
        m_ddb, sd_ddb, n_ddb, has_pair = per_rep_ddbind(d) if d else (None, None, 0, None)
        rec["dG_bind_HIP"] = (d or {}).get("dG_bind_HIP", "")
        rec["dG_bind_HIE"] = (d or {}).get("dG_bind_HIE", "")
        rec["dd_bind_kcal"] = round(m_ddb, 2) if m_ddb is not None else ""
        rec["dd_bind_sd"] = round(sd_ddb, 2) if sd_ddb is not None else ""
        rec["dd_bind_nrep"] = n_ddb
        rec["mdbind_d_occ"] = (d or {}).get("d_occupancy", "")
        rec["has_sensor_pair"] = has_pair if has_pair is not None else ""

        # --- four legs, presence-tracked (None = not measured, not a failing leg) ---
        v_propka = present_float(r.get("ddg_pH"))
        v_pkai = present_float(r.get("pkai_his433_dpKa"))
        v_mdocc = present_float(r.get("md_docc"))
        leg_propka = None if v_propka is None else int(v_propka > 0.0)
        leg_pkai = None if v_pkai is None else int(v_pkai > 0.0)
        leg_mdocc = None if v_mdocc is None else int(v_mdocc > 0.1)   # 0.1 = engine's own switch margin
        if m_ddb is None or (has_pair is False):
            leg_ddb = None   # no data, or no mechanistic sensor pair -> not a valid binding leg
        else:
            leg_ddb = int(m_ddb < 0.0)
        legs = [leg_propka, leg_pkai, leg_mdocc, leg_ddb]
        rec["leg_propka"], rec["leg_pkai"], rec["leg_mdocc"] = leg_propka, leg_pkai, leg_mdocc
        rec["leg_ddbind"] = leg_ddb if leg_ddb is not None else ""
        rec["quad_score"] = sum(1 for x in legs if x == 1)
        rec["quad_present"] = sum(1 for x in legs if x is not None)

        # dd_bind robustness
        robust_reverse = (m_ddb is not None and sd_ddb is not None and m_ddb > 0.0 and (m_ddb - sd_ddb) > 0.0)
        ambiguous_ddb = (m_ddb is not None and m_ddb >= 0.0 and not robust_reverse)

        # --- verdict (priority order) ---
        if v_pkai is not None and v_pkai < 0.0:
            rec["ph_verdict"] = f"VETO: pKAI-reverse (blood-ON, dpKa={v_pkai})"
        elif has_pair is False:
            rec["ph_verdict"] = "no His433 sensor pair (dd_bind not a switch readout)"
        elif m_ddb is None:
            rec["ph_verdict"] = "dd_bind pending"
        elif robust_reverse:
            rec["ph_verdict"] = f"RECONSIDER: robust reverse binding-switch (dd_bind={round(m_ddb,2)}±{round(sd_ddb,2)})"
        elif ambiguous_ddb:
            rec["ph_verdict"] = f"ambiguous binding-dG (dd_bind={round(m_ddb,2)}±{round(sd_ddb,2)}, CI spans 0)"
        elif rec["quad_score"] == rec["quad_present"] == 4:
            rec["ph_verdict"] = "quad-confirmed (4/4)"
        else:
            rec["ph_verdict"] = f"{rec['quad_score']}/{rec['quad_present']} legs"

        # pH-window fold-change (Wyman legs only; sign/order-of-magnitude, not calibrated magnitude)
        rec["fold_propka"] = fold(present_float(r.get("ddg_pH")))
        rec["fold_pkai"] = fold(present_float(r.get("pkai_ddg_pH")))
        g = gp.get(nm, (None, []))
        rec["glycan_min_A"] = g[0] if g[0] is not None else ""
        rec["glycan_contacted"] = ";".join(str(x) for x in g[1]) if g[1] else ""
        out.append(rec)

    base = list(rows[0].keys())
    newc = ["dG_bind_HIP", "dG_bind_HIE", "dd_bind_kcal", "dd_bind_sd", "dd_bind_nrep", "mdbind_d_occ",
            "has_sensor_pair", "leg_propka", "leg_pkai", "leg_mdocc", "leg_ddbind", "quad_score",
            "quad_present", "fold_propka", "fold_pkai", "glycan_min_A", "glycan_contacted", "ph_verdict"]
    cols = base + [c for c in newc if c not in base]
    with open(OUT_MAN, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for rec in out:
            w.writerow({k: rec.get(k, "") for k in cols})
    small = ["name", "source", "dd_bind_kcal", "dd_bind_sd", "dd_bind_nrep", "mdbind_d_occ",
             "ddg_pH", "pkai_his433_dpKa", "md_docc", "leg_propka", "leg_pkai", "leg_mdocc", "leg_ddbind",
             "quad_score", "quad_present", "fold_propka", "fold_pkai", "glycan_min_A", "ph_verdict"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=small)
        w.writeheader()
        for rec in out:
            w.writerow({k: rec.get(k, "") for k in small})

    ndone = sum(1 for rec in out if rec["dd_bind_kcal"] != "")
    nquad = sum(1 for rec in out if rec["ph_verdict"] == "quad-confirmed (4/4)")
    print(f"dd_bind available {ndone}/{len(out)} | quad-confirmed (4/4): {nquad}/{len(out)}\n")
    print("%-14s %8s %5s | %4s%4s%4s%4s | %-5s | %6s %6s | %5s %s" % (
        "name", "dd_bind", "nrep", "prp", "pki", "occ", "ddb", "quad", "f_prp", "f_pki", "glyA", "verdict"))
    for rec in sorted(out, key=lambda x: (-(x["quad_score"]), present_float(x.get("dd_bind_kcal")) if present_float(x.get("dd_bind_kcal")) is not None else 99)):
        dds = f"{rec['dd_bind_kcal']}±{rec['dd_bind_sd']}" if rec["dd_bind_kcal"] != "" else "-"
        print("%-14s %8s %5s | %4s%4s%4s%4s | %s/%s | %6s %6s | %5s %s" % (
            rec["name"], dds, rec["dd_bind_nrep"],
            rec["leg_propka"], rec["leg_pkai"], rec["leg_mdocc"], rec["leg_ddbind"],
            rec["quad_score"], rec["quad_present"], rec["fold_propka"], rec["fold_pkai"],
            rec["glycan_min_A"], rec["ph_verdict"]))
    flags = [rec["name"] for rec in out if rec["ph_verdict"].startswith(("RECONSIDER", "VETO"))]
    amb = [rec["name"] for rec in out if rec["ph_verdict"].startswith("ambiguous")]
    print(f"\nRECONSIDER/VETO: {flags or 'none'}")
    print(f"ambiguous binding-dG (keep, note uncertainty): {amb or 'none'}")
    gl = sorted([rec for rec in out if rec["glycan_min_A"] != ""], key=lambda x: present_float(x["glycan_min_A"]))
    print("glycan closest: " + ", ".join(f"{rec['name']}={rec['glycan_min_A']}A" for rec in gl[:3]))
    print(f"wrote {OUT_CSV}\nwrote {OUT_MAN}")


if __name__ == "__main__":
    main()
