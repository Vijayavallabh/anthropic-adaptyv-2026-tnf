#!/usr/bin/env python3
"""Compute ipSAE_min (Anthropic/Overath-validated binder metric) from Boltz PAE for a set of designs.
ipSAE_min = min(ipSAE(A->B), ipSAE(B->A)); Overath F1-optimal cutoff ~0.61 (monomer target).
Usage: batch_ipsae.py <design_list.csv (design_id col)> <boltz_glob_dirs colon-sep> [tag]
Finds each design's Boltz human complex (cif + pae npz), runs tools/ipsae.py, parses the two asym rows."""
import sys, os, csv, glob, subprocess, re
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
listcsv=sys.argv[1]
globdirs=sys.argv[2].split(":") if len(sys.argv)>2 else ["runs/boltz_*_h/out/by_seed"]
ids=[r.get("design_id") or r.get("name") for r in csv.DictReader(open(listcsv))]
def find_boltz(did):
    for gd in globdirs:
        for cif in glob.glob(f"{R}/{gd}/{did}/s0/{did}_model_0.cif"):
            pae=os.path.join(os.path.dirname(cif),f"pae_{did}_model_0.npz")
            if os.path.exists(pae): return cif,pae
    return None,None
def ipsae_min(did):
    cif,pae=find_boltz(did)
    if not cif: return None,"no-boltz-pae"
    out=cif.replace(".cif","_10_10.txt")
    if not os.path.exists(out):
        subprocess.run([PY,f"{R}/tools/ipsae.py",pae,cif,"10","10"],capture_output=True,cwd=os.path.dirname(cif))
    if not os.path.exists(out): return None,"ipsae-failed"
    asym=[]
    for ln in open(out):
        f=ln.split()
        if len(f)>6 and f[4]=="asym":
            try: asym.append(float(f[5]))  # ipSAE column
            except: pass
    return (min(asym) if asym else None),("ok" if asym else "no-asym")
rows=[]
for did in ids:
    v,st=ipsae_min(did); rows.append((did,v,st))
rows.sort(key=lambda r:-(r[1] or -1))
print(f"{'design':24}{'ipSAE_min':>10}  (>=0.61 pass; TNF ranking compresses, use relative)")
for did,v,st in rows:
    print(f"{did:24}{(f'{v:.3f}' if v is not None else st):>10}")
import statistics
good=[v for _,v,_ in rows if v is not None]
if good: print(f"\nn={len(good)} | median {statistics.median(good):.3f} | >=0.61: {sum(1 for v in good if v>=0.61)} | >=0.5: {sum(1 for v in good if v>=0.5)}")
