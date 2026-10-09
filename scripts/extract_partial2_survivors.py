#!/usr/bin/env python3
"""Extract strong-binder survivors from the partial2 af2ig run, map each back to its
ProteinMPNN sequence, dedupe vs the current portfolio + among themselves, and emit a
ranked candidate CSV (name,sequence,pae_i,plddt_b,rmsd,seed_fold,instab,gravy,expr).
Gate: pae_interaction<10 AND plddt_binder>80 AND binder_aligned_rmsd<2 (af2ig trust gate).
Usage: extract_partial2_survivors.py [pae_cut=10] [plddt_cut=80] [rmsd_cut=2.0]"""
import sys, os, glob, csv, re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
SC=f"{R}/runs/partial2/af2ig/out.sc"
PAE=float(sys.argv[1]) if len(sys.argv)>1 else 10.0
PLD=float(sys.argv[2]) if len(sys.argv)>2 else 80.0
RMS=float(sys.argv[3]) if len(sys.argv)>3 else 2.0

# portfolio seqs to dedupe against
port=set()
for r in csv.DictReader(open(f"{R}/runs/push/upload_final20.csv")):
    port.add(r["sequence"])

# index MPNN fastas: {seed}_pd_{batch} -> {seqidx: binder_seq}
mpnn={}
for fa in glob.glob(f"{R}/runs/partial2/mpnn/gpu*/seqs/*.fa"):
    key=os.path.basename(fa)[:-3]  # e.g. b1m12_A91H_pd_0
    d={}
    hdr=None
    for ln in open(fa):
        ln=ln.rstrip("\n")
        if ln.startswith(">"):
            m=re.search(r"id=(\d+)", ln); hdr=int(m.group(1)) if m else None
        elif ln and hdr is not None:
            d[hdr]=ln.split(":")[0]  # binder = before colon
    mpnn[key]=d

def seq_for(desc):
    # desc = {seed}_pd_{batch}_{seqidx}[_af2pred]
    d=desc[:-8] if desc.endswith("_af2pred") else desc
    m=re.match(r"(.+_pd_\d+)_(\d+)$", d)
    if not m: return None
    key, idx = m.group(1), int(m.group(2))
    return mpnn.get(key,{}).get(idx)

rows=[]
if not os.path.exists(SC):
    print("NO out.sc yet"); sys.exit(1)
seen=set()
for ln in open(SC):
    if not ln.startswith("SCORE:"): continue
    p=ln.split()
    try:
        rmsd=float(p[1]); pae_i=float(p[3]); plddt_b=float(p[5]); desc=p[-1]
    except: continue
    if not (pae_i<PAE and plddt_b>PLD and rmsd<RMS): continue
    s=seq_for(desc)
    if not s or "X" in s or s in seen or s in port: continue
    seen.add(s)
    d2=desc[:-8] if desc.endswith("_af2pred") else desc
    seed_fold=re.match(r"(.+)_pd_\d+_\d+$", d2).group(1)
    pa=ProteinAnalysis(s); ins=round(pa.instability_index(),1); gr=round(pa.gravy(),2)
    expr="PASS" if ins<40 and gr<0.4 and s.count("C")==0 else "FLAG"
    rows.append(dict(name=desc,sequence=s,pae_i=round(pae_i,2),plddt_b=round(plddt_b,1),
                     rmsd=round(rmsd,2),seed_fold=seed_fold,instab=ins,gravy=gr,expr=expr))
rows.sort(key=lambda r:r["pae_i"])
out=f"{R}/runs/partial2/survivors.csv"
with open(out,"w") as fh:
    w=csv.DictWriter(fh,fieldnames=["name","sequence","pae_i","plddt_b","rmsd","seed_fold","instab","gravy","expr"])
    w.writeheader()
    for r in rows: w.writerow(r)
import collections
byfold=collections.Counter(r["seed_fold"] for r in rows)
print(f"strong survivors (pae<{PAE} plddt>{PLD} rmsd<{RMS}), unique + novel-vs-portfolio: {len(rows)}")
print(f"by seed-fold: {dict(byfold)}")
print(f"expr-PASS among survivors: {sum(1 for r in rows if r['expr']=='PASS')}/{len(rows)}")
print(f"-> {out}")
for r in rows[:15]:
    print(f"  {r['name']:28} pae_i={r['pae_i']:5} plddt={r['plddt_b']:5} rmsd={r['rmsd']:4} {r['expr']:4} {r['seed_fold']}")
