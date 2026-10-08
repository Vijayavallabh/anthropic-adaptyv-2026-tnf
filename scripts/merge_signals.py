#!/usr/bin/env python3
"""Merge every selection signal per candidate variant into one master table for final portfolio choice.

Signals (objective priority: pH-selectivity > cross-react > affinity; novelty + diversity as judging axes):
  * MD gate (runs/md/gate/*.json): d_switch (HIP-HIE), d_iso (HIP-GLN), verdict  [SELECTIVITY, primary]
  * static ddg (switch_table_rescored.tsv): release linkage, more negative = stronger            [SELECTIVITY]
  * novelty TM (novelty_struct/pool_tm.tsv): max structural TM vs PDB, <0.80 = novel              [NOVELTY]
  * binding (switch_table_rescored.tsv): af2ig pae_interaction, plddt, rmsd                        [AFFINITY]
  * lineage + anchor (manifest_pool.tsv): for diversity cap + conserved-Lys cross-react
Prints a ranked master table; does NOT auto-select (final portfolio is a judged call).
Usage: merge_signals.py
"""
import glob, json, os, re

def lin(v):
    if v.startswith("cb_1"): return "cb1"
    if v.startswith("cb_7"): return "cb7"
    m = re.match(r'tnf_K90_tnf_K90_(\d+)_', v); return "dn" + m.group(1) if m else "?"

# static binding + ddg
stat = {}
for ln in open("runs/pool/switch_table_rescored.tsv").read().splitlines()[1:]:
    pae, plddt, rmsd, old, new, v = ln.split("\t")
    stat[v] = dict(pae=float(pae), plddt=float(plddt), rmsd=float(rmsd), ddg=float(new))
# novelty
tm = {}
for ln in open("runs/novelty_struct/pool_tm.tsv"):
    v, t, l = ln.rstrip("\n").split("\t"); tm[v] = float(t)
# his + anchor
man = {}
for ln in open("runs/md/manifest_pool.tsv"):
    v, hp, la, pdb = ln.rstrip("\n").split("\t"); man[v] = (int(hp), la)
# MD verdict
md = {}
for f in glob.glob("runs/md/gate/*.json"):
    d = json.load(open(f)); agg = d.get("agg", {})
    def gg(va, m): return agg.get(va, {}).get(m)
    hip, hie, gln = gg("HIP", "his_lys_mindist"), gg("HIE", "his_lys_mindist"), gg("GLN", "his_lys_mindist")
    sem = max([x for x in (gg("HIP", "sem_his_lys"), gg("HIE", "sem_his_lys"), gg("GLN", "sem_his_lys")) if x is not None] or [0.0])
    dsw = (hip - hie) if None not in (hip, hie) else None
    diso = (hip - gln) if None not in (hip, gln) else None
    vd = "NA"
    if dsw is not None:
        vd = "CONFIRMED" if (dsw > 2 * sem and dsw > 0 and (diso is None or diso > 0)) else \
             ("NONRESPOND" if dsw <= 2 * sem else "MARGINAL")
    md[d["name"]] = dict(dsw=dsw, diso=diso, sem=sem, vd=vd)

rows = []
for v in stat:
    s = stat[v]; m = md.get(v, {})
    rows.append(dict(v=v, lin=lin(v), his=man.get(v, ('?', '?'))[0], anc=man.get(v, ('?', '?'))[1],
                     pae=s["pae"], plddt=s["plddt"], rmsd=s["rmsd"], ddg=s["ddg"], tm=tm.get(v),
                     dsw=m.get("dsw"), diso=m.get("diso"), vd=m.get("vd", "pending")))
# rank: CONFIRMED first, then by static ddg (release), novelty tiebreak
order = {"CONFIRMED": 0, "MARGINAL": 1, "pending": 2, "NA": 3, "NONRESPOND": 4}
rows.sort(key=lambda r: (order.get(r["vd"], 5), r["ddg"], r["tm"] if r["tm"] else 9))
print(f"{'variant':50} {'lin':4} {'anc':10} {'pae':>5} {'ddg':>6} {'TM':>5} {'dSw':>5} {'dIso':>5} verdict")
for r in rows:
    def fm(x, w=5, p=2): return (f"{x:>{w}.{p}f}" if isinstance(x, (int, float)) else f"{'--':>{w}}")
    nov = "*" if (r["tm"] and r["tm"] < 0.80) else " "
    print(f"{r['v'][:50]:50} {r['lin']:4} {r['anc']:10} {fm(r['pae'])} {fm(r['ddg'],6)} {fm(r['tm'])}{nov}{fm(r['dsw'])} {fm(r['diso'])} {r['vd']}")
print("\n* = structurally novel (TM<0.80). verdict from MD ladder (pending = MD not finished).")

if __name__ == "__main__":
    pass
