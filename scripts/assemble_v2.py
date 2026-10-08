#!/usr/bin/env python3
"""Assemble the final diverse portfolio (data-driven, no hardcoded metrics).
  cb1-family floor (floor_core8.csv) + cb1-epitope de-novo CONFIRMED folds (confirmed.tsv per batch),
  cross-react joined by FULL design name from xreact.tsv (missing xr => fails the cross-react gate).
Caps the cb1 lineage (diversity: cap cousins), keeps one design per distinct de-novo backbone, ranks
by objective priority (pH-sel > cross-react > affinity, + expression), fills <= MAX.
Usage: assemble_v2.py --out <prefix> [--cb1cap 6] [--max 20] --confirmed f1.tsv f2.tsv --xreact x.tsv
Gates (documented Adaptyv + switch): dual_sw(d_switch>0 & ddg<0), xr_clean(xr>-0.10, missing=FAIL),
  expr_ok(cys==0 & instab<40 & gravy<0.4); expr_accept(cys==0 & instab<50 & gravy<0.5).
"""
import csv, re, argparse, os
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=f"{R}/runs/submission/final_portfolio")
ap.add_argument("--cb1cap", type=int, default=6)
ap.add_argument("--max", type=int, default=20)
ap.add_argument("--confirmed", nargs="*", default=[f"{R}/runs/cb1ep/final/confirmed.tsv", f"{R}/runs/cb1ep2/final/confirmed.tsv"])
ap.add_argument("--xreact", default=f"{R}/runs/cb1ep/final/xreact.tsv")
ap.add_argument("--clusters", default=f"{R}/runs/submission/diversity/clusters.tsv",
                help="name<TAB>cluster from fold_diversity.py; cb1ep designs dedup by structural CLUSTER (true shape), not backbone id")
ap.add_argument("--portal", default=f"{R}/runs/submission/portal_novelty.tsv",
                help="design_id<TAB>novelty from the Proteinbase portal; designs with score<3 are DROPPED (authoritative novelty gate). Designs absent = unchecked (kept, flagged).")
ap.add_argument("--extra", default="", help="comma-separated design_ids to force-include as deliberate hedges/diversity shots (bypass ipSAE/cluster dedup; still require switch + portal>=3)")
A = ap.parse_args()
portal = {}
if os.path.exists(A.portal):
    import csv as _csv
    for r in _csv.DictReader(open(A.portal), delimiter="\t"):
        try: portal[r["design_id"]] = int(r["novelty"])
        except (ValueError, KeyError): pass
# worst-case-over-species ipSAE (sharper interface metric; adopted cross-react #2): token -> ipsae_wc
ipsae_wc = {}
_ipf = f"{R}/runs/submission/ipsae_all.tsv"
if os.path.exists(_ipf):
    import csv as _csv2
    for r in _csv2.DictReader(open(_ipf), delimiter="\t"):
        try: ipsae_wc[r["token"]] = float(r["ipsae_wc"])
        except (ValueError, KeyError): pass
# CONTACT gate (switch_contact_gate.py): a d_switch is a REAL protonation switch only if the switch-His
# is in CONTACT with the anchor Lys in the ON state (HIE<8A) AND d_switch>2*max(sem incl GLN isostere).
# Added after an adversarial review caught a non-contact N-terminal drift artifact (A1M, HIE 13.3A) that
# passed the distance-ladder-only gate. real_switch[name] -> bool; hie_contact[name] -> float.
real_switch, hie_contact = {}, {}
_scf = f"{R}/runs/submission/switch_contact.tsv"
if not os.path.exists(_scf):
    # FAIL LOUD: the contact gate drives real_switch_gate; without it every cb1ep design silently
    # fails the gate and the whole de-novo portfolio is dropped. Refuse rather than ship a cb1-only CSV.
    raise SystemExit(f"ERROR: contact-gate table missing: {_scf}\n"
                     f"  Run `python scripts/switch_contact_gate.py` first (required: the real_switch gate "
                     f"needs it; without it all cb1ep designs are dropped).")
import csv as _csv3
for r in _csv3.DictReader(open(_scf), delimiter="\t"):
    real_switch[r["name"]] = (r.get("real_switch") == "True")
    try: hie_contact[r["name"]] = float(r["hie_contact"])
    except (ValueError, KeyError): pass

def fnum(x, d=None):
    try: return float(x)
    except (ValueError, TypeError): return d

rows = []
# ---- cb1-family floor ----
for r in csv.DictReader(open(f"{R}/runs/submission/floor_core8.csv")):
    rows.append(dict(id=r["design_id"], lineage="cb1", fold="cb1", cluster="cb1", seq=r["sequence"], length=int(r["length"]),
        d_switch=fnum(r["d_switch_MD"], 0), ddg=fnum(r["ddg_pH_PROPKA"]), xr=fnum(r.get("xreact_delta")),
        iptm_h=fnum(r.get("iptm_human")), instab=fnum(r["instability"], 999), gravy=fnum(r["gravy"], 9),
        cys=int(r.get("n_cys", 0)), tm=fnum(r.get("novelty_TM"), 0.99), switch=r.get("switch_type", "")))
# ---- cb1-epitope CONFIRMED folds (one per distinct backbone across batches) ----
xr = {}
if os.path.exists(A.xreact):
    for r in csv.DictReader(open(A.xreact), delimiter="\t"): xr[r["name"]] = fnum(r.get("xr"))
clu = {}
if os.path.exists(A.clusters):
    for r in csv.DictReader(open(A.clusters), delimiter="\t"): clu[r["name"]] = "C" + r["cluster"]
for tsv in A.confirmed:
    if not os.path.exists(tsv): continue
    tag = "2" if "/cb1ep2/" in tsv else "1"   # batch tag: same cb1_epitope_N_model_M index is a DIFFERENT structure across RFd3 seeds
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        bb = r["backbone"]
        short = re.search(r"cb1_epitope_(\d+)_model_(\d+)", bb)
        fold = f"cb1ep{tag}:b{short.group(1)}m{short.group(2)}" if short else f"cb1ep{tag}:" + bb
        # fold-diversity dedup key = STRUCTURAL CLUSTER (true shape from fold_diversity.py), not backbone id
        cluster = clu.get(r["name"], f"{tag}:{bb}")
        swtok = re.search(r"A\d+[A-Z]2H_vsK\d+", r["name"])
        rows.append(dict(id=f"{fold.split(':')[1]}_{tag}A{r['his_pos']}H_{r['anchor']}", lineage="cb1ep", fold=fold,
            cluster=cluster, seq=r["sequence"], length=int(r["length"]), d_switch=fnum(r["d_switch"], 0), ddg=fnum(r["ddg_pH"]),
            xr=xr.get(r["name"]), iptm_h=None, instab=fnum(r["instab"], 999), gravy=fnum(r["gravy"], 9),
            ipsae=ipsae_wc.get(swtok.group(0)) if swtok else None,
            real_sw=real_switch.get(r["name"]), hie=hie_contact.get(r["name"]),
            cys=int(r["cys"]), tm=None, switch="denovo-His", _name=r["name"]))
# ---- gates ----
for r in rows:
    r["dual_sw"] = (r["d_switch"] or 0) > 0 and (r["ddg"] is not None and r["ddg"] < 0)
    # CONTACT-gated switch: cb1ep must pass the contact gate (real contact switch, not a drift artifact);
    # cb1-family switches are contact-confirmed separately (HIE 3.5-6 A, runs/md/gate) -> gate True.
    r["real_switch_gate"] = True if r["lineage"] == "cb1" else bool(r.get("real_sw"))
    r["xr_clean"] = r["xr"] is not None and r["xr"] > -0.10
    r["expr_ok"] = r["cys"] == 0 and (r["instab"] or 999) < 40 and (r["gravy"] or 9) < 0.4
    r["expr_accept"] = r["cys"] == 0 and (r["instab"] or 999) < 50 and (r["gravy"] or 9) < 0.5
    r["portal"] = portal.get(r["id"])   # authoritative Proteinbase novelty score (None = unchecked)
    r["novel"] = (r["portal"] >= 3) if r["portal"] is not None else ((r["lineage"] == "cb1ep") or ((r["tm"] or 0.99) < 0.80))
    r.setdefault("ipsae", None)
    # cross-react ROBUST (sharper than the ipTM delta): worst-case-over-species ipSAE>0.5 (cb1-family: no ipSAE, rely on matched-Boltz xr_clean)
    r["xreact_robust"] = (r["ipsae"] is None) or (r["ipsae"] > 0.5)
    # composite: pH-sel (d_switch, -ddg) dominant, then cross-react, then expression
    r["score"] = (2.0*(r["d_switch"] or 0) + 0.8*(-(r["ddg"] or 0)) + 1.0*((r["xr"] or -0.1)+0.1)*5
                  + (1.0 if r["expr_ok"] else 0.3 if r["expr_accept"] else 0))

# ---- print full gated table ----
print(f"{'id':26} {'fold':14} {'len':>3} {'dSw':>5} {'ddg':>6} {'xr':>6} {'instab':>6} {'grv':>5} flags")
for r in sorted(rows, key=lambda r: (r["lineage"], -r["score"])):
    fl = [f for f, ok in (("dual?", r["dual_sw"]), ("contact?", r["real_switch_gate"]), ("xr?", r["xr_clean"]), ("ipSAE?", r["xreact_robust"]), ("expr", r["expr_ok"] or r["expr_accept"]), ("nov?", r["novel"])) if not ok]
    xs = f"{r['xr']:+.2f}" if r["xr"] is not None else "  NA"
    print(f"{r['id'][:26]:26} {r['fold']:14} {r['length']:>3} {(r['d_switch'] or 0):5.2f} {(r['ddg'] or 0):6.2f} {xs:>6} {(r['instab'] or 0):6.1f} {(r['gravy'] or 0):5.2f} {','.join(fl) if fl else 'CLEAN'}")

# ---- select ----
cb1 = [r for r in rows if r["lineage"] == "cb1" and r["dual_sw"] and r["novel"] and (r["expr_accept"] or r["d_switch"] > 2.0)]
cb1 = sorted(cb1, key=lambda r: -r["score"])[:A.cb1cap]
# cb1ep: switch + cross-react-ROBUST (ipSAE worst-case>0.5) + PORTAL-validated(>=3) + expr_accept(or strong switch), then ONE best design per STRUCTURAL CLUSTER
cb1ep = [r for r in rows if r["lineage"] == "cb1ep" and r["dual_sw"] and r["real_switch_gate"] and r["xr_clean"] and r["xreact_robust"]
         and (r["portal"] is not None and r["portal"] >= 3) and (r["expr_accept"] or r["d_switch"] > 2.0)]
cb1ep = sorted(cb1ep, key=lambda r: -r["score"])
seen_clu, cb1ep_sel = set(), []
for r in cb1ep:
    if r["cluster"] in seen_clu: continue   # dedup by true shape, not backbone id
    seen_clu.add(r["cluster"]); cb1ep_sel.append(r)
# force-include deliberate hedges / diversity risk-shots (--extra): switch + portal gated, ipSAE/dedup bypassed
extra_ids = set(x.strip() for x in A.extra.split(",") if x.strip())
sel_ids = {r["id"] for r in cb1 + cb1ep_sel}
extras = [r for r in rows if r["id"] in extra_ids and r["id"] not in sel_ids
          and r["dual_sw"] and r["real_switch_gate"] and (r["portal"] is not None and r["portal"] >= 3)]
for r in extras: r["hedge"] = True
final = (cb1 + cb1ep_sel + extras)[:A.max]
final = sorted(final, key=lambda r: (r["lineage"], -r["score"]))

def nov_str(r): return f"{r['portal']}/4" if r.get("portal") is not None else (str(r["tm"]) if r["tm"] is not None else "denovo")
cols = ["design_id","fold","sequence","length","d_switch_MD","hie_contact_A","ddg_pH_PROPKA","xreact_delta","ipSAE_wc","instability","gravy","n_cys","novelty_portal","switch_type"]
with open(A.out + ".csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(cols)
    for r in final:
        w.writerow([r["id"], r["fold"], r["seq"], r["length"], round(r["d_switch"] or 0, 3),
                    round(r["hie"], 2) if r.get("hie") is not None else "",
                    round(r["ddg"], 3) if r["ddg"] is not None else "", r["xr"] if r["xr"] is not None else "",
                    r.get("ipsae") if r.get("ipsae") is not None else "", r["instab"], r["gravy"], r["cys"], nov_str(r), r["switch"]])
with open(A.out + ".fasta", "w") as f:
    for r in final:
        f.write(f">{r['id']} fold={r['fold']} dSw={round(r['d_switch'] or 0,2)} ddg={round(r['ddg'],2) if r['ddg'] is not None else 'NA'} xr={r['xr'] if r['xr'] is not None else 'NA'} nov={nov_str(r)}\n{r['seq']}\n")
clusters = sorted(set(r["cluster"] for r in final))
print(f"\nSELECTED {len(final)} designs | cb1={sum(1 for r in final if r['lineage']=='cb1')} cb1ep={sum(1 for r in final if r['lineage']=='cb1ep')} | {len(clusters)} distinct SHAPES (TM-clusters)")
print("shapes:", clusters)
print(f"wrote {A.out}.{{csv,fasta}}")
