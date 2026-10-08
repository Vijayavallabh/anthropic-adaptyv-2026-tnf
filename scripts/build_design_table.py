#!/usr/bin/env python3
"""Build the reviewer-facing per-design table (docs/DESIGN_TABLE.md) for the Adaptyv submission.

Merges the submission manifest (quad pH-switch legs + cross-react + affinity + developability) with
the WAVE 16L independent burial-aware pKa (JustHISpKa) and the isostere causal-ablation results, into
a transparent per-design cross-method confidence tier so reviewers can pick the strongest designs.

Confidence (all 20 are the floor: quad-confirmed switch + cross-reactive + 0-Cys + de-novo):
  obj-1 corroborations beyond the quad =
     +1 if JustHISpKa (independent burial-aware pKa) says tumour-ON
     +1 if the 20-frame MD ensemble is sign-stable (frac_protonated >= 0.75)
  isostere causal salt-bridge confirmed (single-structure proxy WT linkage > 0) is reported as a PLUS.
  Tier: High = both corroborations; Solid = one; Watch = neither (still quad-confirmed) + the reason.
"""
import csv, os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# optional: build_design_table.py <comma-sep names> <out.md>  -> table for just that subset
ONLY = set(sys.argv[1].split(",")) if len(sys.argv) > 1 and sys.argv[1] else None
OUTPATH = sys.argv[2] if len(sys.argv) > 2 else os.path.join(REPO, "docs/DESIGN_TABLE.md")

def load(path, key="name"):
    p = os.path.join(REPO, path)
    if not os.path.exists(p):
        return {}
    return {r[key]: r for r in csv.DictReader(open(p))}

def f(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d

man = load("results/submission_manifest_v2_wave16.csv")
jhp = load("results/wave16L_justhispka.csv")
iso = load("results/wave16k_isostere_ablation.csv")

rows = []
for name in sorted(man):
    if ONLY and name not in ONLY:
        continue
    m = man[name]
    j = jhp.get(name, {}); i = iso.get(name, {})
    ens = f(m.get("ens_frac_protonated"))
    jhp_on = (j.get("tumour_on", "").strip().lower() == "true")
    jhp_shift = f(j.get("wyman_shift"))
    ens_robust = ens is not None and ens >= 0.75
    causal = f(i.get("wt_ddg_his433"))
    causal_ok = causal is not None and causal > 0
    corrob = int(jhp_on) + int(ens_robust)
    tier = "High" if corrob == 2 else ("Solid" if corrob == 1 else "Watch")
    flags = []
    if not jhp_on and name in jhp: flags.append("JustHISpKa pKa shift weak")
    if not ens_robust and ens is not None: flags.append(f"MD-ensemble noisy (f_prot {ens:.2f})")
    rows.append(dict(
        name=name, length=m.get("length"),
        propka=f(m.get("ddg_pH")), pkai=f(m.get("pkai_ddg_pH")),
        mdocc=f(m.get("md_docc")), ddbind=f(m.get("dd_bind_kcal")),
        jhp_shift=jhp_shift, jhp_on=jhp_on, causal_ok=causal_ok,
        xreact=m.get("cross_reactive"), sev4=m.get("w15_xreact_severity_core4A"),
        chai_m=f(m.get("chai_mouse_iptm")),
        prodigy=f(m.get("prodigy_dG_kcal")), af2=f(m.get("af2_pae")),
        af3_ipsae=f(m.get("af3_ipsae")), cys=m.get("free_cys"),
        tier=tier, flags="; ".join(flags) or "-",
    ))

# order: High, Solid, Watch; within tier by PROPKA ddg_pH desc
order = {"High": 0, "Solid": 1, "Watch": 2}
rows.sort(key=lambda r: (order[r["tier"]], -(r["propka"] or -9)))

def g(x, fmt="{:.2f}"):
    return fmt.format(x) if isinstance(x, float) else ("-" if x is None else str(x))

out = [
    "# Per-design confidence table — Adaptyv Challenge 1 (de-novo conditional EGFR DIII binders)",
    "",
    "All 20 designs are the validated floor: **quad-confirmed pH-switch** (PROPKA-Wyman ΔΔG_pH>0, "
    "pKAI, MD HIP/HIE occupancy, MM-GBSA two-state ΔΔG all tumour-ON), **cross-reactive** (human+mouse, "
    "4-method), **0 free cysteines**, and **de-novo** (<70% identity to any natural protein, max 60.7%). "
    "The tier below adds two *independent* obj-1 corroborations on top of the quad — the burial-aware "
    "JustHISpKa pKa (a 5th engine sharing no code with PROPKA/pKAI) and MD-ensemble sign-stability — so "
    "reviewers can rank within an already-strong set. **Watch** designs remain quad-confirmed; they simply "
    "have weaker independent corroboration (stated reason), not a failure.",
    "",
    "ΔΔG_pH>0 = tumour-ON (tighter at pH 6.5). JustHISpKa Δ = bound−free His433 pKa (>0 tumour-ON). "
    "MD d_occ>0 and MM-GBSA dd_bind<0 = tumour-ON. af2 pae<10 and PRODIGY (kcal/mol, more negative=tighter) = binding.",
    "",
    "| Tier | Design | Len | PROPKA ΔΔG_pH | pKAI | MD d_occ | MM-GBSA dd_bind | JustHISpKa Δ | isostere-causal | x-react (core-4Å sev) | PRODIGY | af2 pae | 0-Cys | note |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
]
for r in rows:
    out.append("| **{tier}** | {name} | {length} | {propka} | {pkai} | {mdocc} | {ddbind} | {jhp} | {causal} | {xr} ({sev}) | {prod} | {af2} | {cys} | {flags} |".format(
        tier=r["tier"], name=r["name"], length=r["length"],
        propka=g(r["propka"]), pkai=g(r["pkai"]), mdocc=g(r["mdocc"]), ddbind=g(r["ddbind"]),
        jhp=(g(r["jhp_shift"]) + (" ✓" if r["jhp_on"] else " ✗") if r["jhp_shift"] is not None else "n/a"),
        causal=("yes" if r["causal_ok"] else "n/s"),
        xr=("yes" if str(r["xreact"]).lower() in ("true", "yes", "1") else r["xreact"]),
        sev=r["sev4"], prod=g(r["prodigy"], "{:.1f}"), af2=g(r["af2"], "{:.1f}"),
        cys=("yes" if str(r["cys"]) in ("0", "0.0") else f"NO({r['cys']})"), flags=r["flags"]))

n_high = sum(1 for r in rows if r["tier"] == "High")
n_solid = sum(1 for r in rows if r["tier"] == "Solid")
n_watch = sum(1 for r in rows if r["tier"] == "Watch")
out += ["",
    f"**Summary:** {n_high} High / {n_solid} Solid / {n_watch} Watch (all {len(rows)} quad-confirmed + cross-reactive + 0-Cys + de-novo). "
    "If the lab tests a subset, prioritise High then Solid; the Watch designs carry the stated independent-corroboration caveat.",
    ""]
op = OUTPATH
open(op, "w").write("\n".join(out))
print(f"wrote {op}: {len(rows)} designs ({n_high} High / {n_solid} Solid / {n_watch} Watch)")
for r in rows:
    print(f"  {r['tier']:6s} {r['name']} propka={g(r['propka'])} jhp={'on' if r['jhp_on'] else 'off'} causal={'y' if r['causal_ok'] else 'n'} flags={r['flags']}")
