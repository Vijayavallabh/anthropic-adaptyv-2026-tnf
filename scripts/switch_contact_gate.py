#!/usr/bin/env python3
"""Contact gate for the pH-switch validation (added after an adversarial review caught a
non-contact N-terminal drift artifact that passed the distance-ladder gate).

A d_switch (= d_HIP - d_HIE of the His-anchor-Lys distance) is only a REAL protonation switch
if the switch-His is actually IN CONTACT with the anchor Lys in the ON state. A His that drifts
13-17 A from the Lys gives a large |d_switch| that is pure drift noise between two non-contact
distances, not a charge-repulsion switch. We therefore require, jointly:
  (a) CONTACT:  ON-state (HIE) His-Lys min-distance < CONTACT_A  (the His can actually reach the Lys)
  (b) LADDER :  d_switch > 2 * max(SEM over HIP, HIE, GLN)       (signal exceeds the noise floor,
                                                                  INCLUDING the isostere control's SEM)
  (c) CAUSAL :  d_iso = d_HIP - d_GLN > 0                        (Gln isostere cannot protonate => no release)
Writes name<TAB>his_pos<TAB>hie_contact<TAB>d_switch<TAB>two_sem_all<TAB>d_iso<TAB>contact_ok<TAB>real_switch.
Usage: switch_contact_gate.py --mddirs runs/cb1ep/md runs/cb1ep2/md --confirmed ...tsv ... --out switch_contact.tsv
"""
import json, glob, os, csv, argparse
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONTACT_A = 8.0   # His sidechain (~4 A) + Lys sidechain (~5 A) reach; 13+ A = non-contact drift
ap = argparse.ArgumentParser()
ap.add_argument("--mddirs", nargs="+", default=[f"{R}/runs/cb1ep/md", f"{R}/runs/cb1ep2/md"])
ap.add_argument("--confirmed", nargs="+", default=[f"{R}/runs/cb1ep/final/confirmed.tsv", f"{R}/runs/cb1ep2/final/confirmed.tsv"])
ap.add_argument("--out", default=f"{R}/runs/submission/switch_contact.tsv")
A = ap.parse_args()

def md_json(name):
    for d in A.mddirs:
        p = os.path.join(d, name + ".json")
        if os.path.exists(p): return p
    hits = []
    for d in A.mddirs: hits += glob.glob(os.path.join(d, f"*{name}*.json"))
    return hits[0] if hits else None

def metric(agg, variant, key):
    return (agg.get(variant) or {}).get(key)

rows = []
for tsv in A.confirmed:
    if not os.path.exists(tsv): continue
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        name = r["name"]
        j = md_json(name)
        if not j:
            rows.append((name, r.get("his_pos", "?"), None, None, None, None, False, False)); continue
        d = json.load(open(j)); agg = d.get("agg", {})
        hip = metric(agg, "HIP", "his_lys_mindist"); hie = metric(agg, "HIE", "his_lys_mindist"); gln = metric(agg, "GLN", "his_lys_mindist")
        sems = [metric(agg, v, "sem_his_lys") for v in ("HIP", "HIE", "GLN")]
        two_sem = 2 * max([s for s in sems if s is not None] or [0])
        dsw = (hip - hie) if None not in (hip, hie) else None
        diso = (hip - gln) if None not in (hip, gln) else None
        contact_ok = hie is not None and hie < CONTACT_A
        real = bool(contact_ok and dsw is not None and dsw > two_sem and dsw > 0 and (diso is None or diso > 0))
        rows.append((name, d.get("his_pos"), hie, dsw, two_sem, diso, contact_ok, real))

os.makedirs(os.path.dirname(A.out), exist_ok=True)
with open(A.out, "w") as f:
    f.write("name\this_pos\thie_contact\td_switch\ttwo_sem_all\td_iso\tcontact_ok\treal_switch\n")
    for name, hp, hie, dsw, ts, diso, cok, real in rows:
        def s(x): return f"{x:.3f}" if isinstance(x, float) else ("" if x is None else str(x))
        f.write(f"{name}\t{s(hp)}\t{s(hie)}\t{s(dsw)}\t{s(ts)}\t{s(diso)}\t{cok}\t{real}\n")

n_real = sum(1 for r in rows if r[7]); n_art = sum(1 for r in rows if r[6] is False and r[2] is not None)
print(f"wrote {A.out}  ({len(rows)} designs; {n_real} REAL contact switches; {n_art} non-contact artifacts)")
print(f"{'name (short)':34} {'his':>4} {'HIE':>6} {'dSw':>6} {'2sem':>6} {'dIso':>6} contact real")
for name, hp, hie, dsw, ts, diso, cok, real in sorted(rows, key=lambda r: (r[6] is not False, r[3] or 0), reverse=True):
    short = name.split("__")[-1].replace("_af2pred", "")
    f6 = lambda x: f"{x:6.2f}" if isinstance(x, float) else "    NA"
    flag = "" if real else ("  <-- NON-CONTACT ARTIFACT" if cok is False and hie is not None else "  <-- weak/failed ladder")
    print(f"{short:34} {str(hp):>4} {f6(hie)} {f6(dsw)} {f6(ts)} {f6(diso)}   {str(cok):5} {str(real):5}{flag}")
