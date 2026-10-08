#!/usr/bin/env python3
"""Finalize a cb1-epitope MD batch: read two-state MD verdict JSONs, keep CONFIRMED switches, and
attach all per-design metrics (sequence, PROPKA ddg_pH, developability) into one TSV.
Boltz cross-react is added separately (GPU) via boltz on the emitted name,sequence CSV.
Usage: finalize_cb1ep_batch.py <md_dir> <out_tsv> [keep=CONFIRMED|CONFIRMED,MARGINAL]
Verdict logic mirrors md_verdict.py exactly (d_switch=HIP-HIE, d_iso=HIP-GLN, 2*sem floor).
"""
import sys, os, re, glob, json, subprocess
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R = "/home/prachh/v/anthropic-adaptyv-2026-tnf"
aa3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
md_dir, out_tsv = sys.argv[1], sys.argv[2]
keep = set((sys.argv[3] if len(sys.argv) > 3 else "CONFIRMED").split(","))

def verdict(d):
    agg = d.get("agg", {})
    def g(v, m): return agg.get(v, {}).get(m)
    hip, hie, gln = g("HIP","his_lys_mindist"), g("HIE","his_lys_mindist"), g("GLN","his_lys_mindist")
    sem = max([x for x in (g("HIP","sem_his_lys"), g("HIE","sem_his_lys")) if x is not None] or [0.0])
    if None in (hip, hie): return None, None, None, "NA"
    dsw, diso = hip - hie, (hip - gln if gln is not None else None)
    if dsw > 2*sem and dsw > 0 and (diso is None or diso > 0): v = "CONFIRMED"
    elif dsw <= 2*sem: v = "NONRESPOND"
    else: v = "MARGINAL"
    return dsw, diso, sem, v

def chainA_seq(pdb):
    seq = {}
    for ln in open(pdb):
        if ln[:4] == "ATOM" and ln[21] == "A" and ln[12:16].strip() == "CA":
            seq[int(ln[22:26])] = aa3.get(ln[17:20].strip(), "X")
    return "".join(seq[k] for k in sorted(seq))

def propka_ddg(pdb):
    try:
        out = subprocess.run(["uv","run","python","scripts/ph_release_score.py",pdb], cwd=R,
                             capture_output=True, text=True, timeout=600).stdout
        m = re.search(r"ddg_pH \(sum.*?=\s*(-?[0-9.]+)", out)
        return float(m.group(1)) if m else None
    except Exception: return None

rows = []
for f in sorted(glob.glob(os.path.join(md_dir, "*.json"))):
    d = json.load(open(f))
    dsw, diso, sem, v = verdict(d)
    if v not in keep: continue
    name, pdb = d["name"], d["pdb"]
    if not os.path.exists(pdb):
        alt = glob.glob(os.path.join(md_dir, "..", "switch_af2ig", "**", "pdbs", name + ".pdb"), recursive=True)
        pdb = alt[0] if alt else pdb
    seq = chainA_seq(pdb)
    bb = (re.search(r"cb1_epitope_\d+_model_\d+", name) or [None]) and re.search(r"cb1_epitope_\d+_model_\d+", name).group(0)
    hm = re.search(r"_A(\d+)[A-Z]2H", name); anc = re.search(r"_vsK(\d+)", name)
    pa = ProteinAnalysis(seq.replace("X",""))
    rows.append(dict(name=name, backbone=bb, his_pos=hm.group(1) if hm else "",
                     anchor="K"+anc.group(1) if anc else "", length=len(seq), sequence=seq,
                     d_switch=round(dsw,3), d_iso=round(diso,3) if diso is not None else "",
                     sem=round(sem,3), verdict=v, ddg_pH=propka_ddg(pdb),
                     gravy=round(pa.gravy(),3), instab=round(pa.instability_index(),1), cys=seq.count("C")))
cols = ["name","backbone","his_pos","anchor","length","d_switch","d_iso","sem","verdict","ddg_pH","gravy","instab","cys","sequence"]
with open(out_tsv, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: -r["d_switch"]):
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
nrel = sum(1 for r in rows if (r["ddg_pH"] or 0) < 0)
print(f"wrote {len(rows)} {'/'.join(keep)} designs -> {out_tsv}  ({nrel} with ddg_pH<0, {sum(1 for r in rows if r['cys']==0)} 0-Cys)")
