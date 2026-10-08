#!/usr/bin/env python3
"""WAVE 18 final merge + gate. Fuse per design: His433 switch (PROPKA ddg_pH) + pKAI (independent pKa)
+ cross-react (af2ig mouse-vs-human-repack delta) + novelty (foldseek TM) + binding (af2ig pae) +
lineage/cluster. Gate to DUAL-confirmed switch (PROPKA>0 AND pKAI-driven) AND cross-reactive
(delta<3.0), score by objective priority (pH-switch > cross-react > binding; novelty secondary), and
write a scored pool for diversity-aware final selection (select_diverse.py)."""
import csv, os, re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(REPO)

def norm(n):
    n = re.sub(r'_(mouse|humanrepack)_af2pred$', '', n)
    n = re.sub(r'_af2pred$', '', n)
    n = re.sub(r'__best$', '', n)
    return n

def load(path, key):
    return {norm(r[key]): r for r in csv.DictReader(open(path))} if os.path.exists(path) else {}

def ff(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def read_sc(path):
    """Rosetta out.sc -> {base_name: pae_interaction}."""
    out, hdr = {}, None
    if not os.path.exists(path): return out
    for l in open(path):
        if not l.startswith("SCORE:"): continue
        p = l.split()[1:]
        if hdr is None and "description" in p: hdr = p; continue
        if hdr is None: continue
        d = dict(zip(hdr, p))
        out[norm(d["description"])] = ff(d.get("pae_interaction"))
    return out

def tb(x): return str(x).strip().lower() in ("true", "1", "yes")

sw   = load("results/rfd_obj1_switch_w18.csv", "design")
pk   = load("results/pkai_w18.csv", "name")
nov  = load("results/novelty_w18.csv", "name")
bind = load("results/rfd_af2_passers_w18.csv", "design")
sl   = load("results/shortlist_w18.csv", "name")
mouse = read_sc("runs/af2ig_xreact_w18_mouse/out.sc")
human = read_sc("runs/af2ig_xreact_w18_human/out.sc")

rows = []
for base, s in sw.items():
    p, nv, b, sh = pk.get(base, {}), nov.get(base, {}), bind.get(base, {}), sl.get(base, {})
    mp, hp = mouse.get(base), human.get(base)
    delta = (mp - hp) if (mp is not None and hp is not None) else None
    ddg, orig = ff(s.get("ddg_pH")), ff(s.get("orig_ddg_pH"))
    tm = ff(nv.get("struct_tm"))
    rows.append({
        "name": base, "sequence": s.get("binder_seq", ""), "length": s.get("binder_len", ""),
        "free_cys": (s.get("binder_seq", "") or "").count("C"),
        "ddg_pH": ddg, "orig_ddg_pH": orig, "his433_dpKa_propka": s.get("his433_dpKa", ""),
        "pkai_his433_dpka": p.get("his433_dpka_pkai", ""), "pkai_driven": tb(p.get("his433_driven_pkai")),
        "pkai_ddg_pH": p.get("ddg_pH_pkai", ""),
        "pae_human": (hp if hp is not None else ff(b.get("pae_interaction"))), "pae_mouse": mp,
        "xreact_delta_pae": (round(delta, 3) if delta is not None else None),
        "struct_tm": tm, "novelty_level": nv.get("level", ""),
        "lineage": sh.get("_lin", ""), "cluster": sh.get("_clu", ""),
        "native_switch": (orig is not None and orig > 0),
    })

def gated(r):
    return (r["ddg_pH"] is not None and r["ddg_pH"] > 0 and r["pkai_driven"]
            and r["xreact_delta_pae"] is not None and r["xreact_delta_pae"] < 3.0
            and r["free_cys"] == 0)   # 0-Cys developability gate, matching the 8 locked

def score(r):
    d = r["xreact_delta_pae"]
    tier = (3 if (r["native_switch"] and d is not None and d < 1.5) else
            2 if (d is not None and d < 1.5) else
            1 if (d is not None and d < 3.0) else 0)
    nov = 1.0 - (r["struct_tm"] if r["struct_tm"] is not None else 1.0)   # lower TM -> higher
    pae = r["pae_mouse"] if r["pae_mouse"] is not None else 12.0
    return round(tier + nov + 0.05 * max(0.0, 12.0 - pae), 4)

for r in rows: r["score"] = score(r)
val = sorted((r for r in rows if gated(r)), key=lambda r: -r["score"])

cols = list(rows[0].keys()) if rows else ["name"]
def wr(path, data):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(data)
wr("results/metrics_w18_full.csv", rows)
wr("results/validated_w18.csv", val)

dd = [r["xreact_delta_pae"] for r in rows if r["xreact_delta_pae"] is not None]
print(f"total switch designs: {len(rows)}")
print(f"PROPKA ddg_pH>0: {sum(1 for r in rows if r['ddg_pH'] and r['ddg_pH']>0)} | native (orig>0): {sum(1 for r in rows if r['native_switch'])}")
print(f"pKAI-driven: {sum(1 for r in rows if r['pkai_driven'])}")
print(f"xreact delta computed: {len(dd)} | delta<1.5: {sum(1 for d in dd if d<1.5)} | delta<3.0: {sum(1 for d in dd if d<3.0)}")
print(f"GATED (dual-switch AND xreact<3.0): {len(val)} | native: {sum(1 for r in val if r['native_switch'])} | delta<1.5: {sum(1 for r in val if r['xreact_delta_pae']<1.5)}")
print(f"  spans {len(set(r['lineage'] for r in val))} lineages, {len(set(r['cluster'] for r in val))} folds")
print("wrote results/metrics_w18_full.csv + results/validated_w18.csv")
