#!/usr/bin/env python3
"""WAVE-6 stage A: score human ipSAE_min for all wave6 variants (Boltz wv6h0-2), rank, select binders.
Calibrates against a few v10 members. third_party/LigandMPNN/.venv/bin/python -I."""
import os, subprocess, glob, csv, json
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
def ipsae(did, tags):
    for tag in tags:
        b=f"{R}/runs/boltz_{tag}/out/by_seed/{did}/s0"; cif=f"{b}/{did}_model_0.cif"
        if not os.path.exists(cif): continue
        out=cif.replace(".cif","_10_10.txt")
        if not os.path.exists(out):
            npz=glob.glob(f"{b}/pae_{did}_model_0.npz")
            if npz: subprocess.run([PY,f"{R}/tools/ipsae.py",os.path.basename(npz[0]),f"{did}_model_0.cif","10","10"],capture_output=True,cwd=b)
        if not os.path.exists(out): continue
        a=[float(l.split()[5]) for l in open(out) if len(l.split())>6 and l.split()[4]=="asym"]
        if a: return min(a)
    return None
def iptm(did, tags):
    for tag in tags:
        b=f"{R}/runs/boltz_{tag}/out/by_seed/{did}/s0"
        cj=glob.glob(f"{b}/confidence_{did}_model_0.json")
        if cj:
            try: return json.load(open(cj[0])).get("iptm")
            except: pass
    return None
rows=list(csv.DictReader(open(f"{R}/runs/wave6/variants.csv")))
HS=["wv6h0","wv6h1","wv6h2"]
out=[]
for r in rows:
    d=r["name"]; h=ipsae(d,HS); it=iptm(d,HS)
    out.append(dict(variant=d,ipsae_h=h,iptm_h=it))
scored=[o for o in out if o["ipsae_h"] is not None]
scored.sort(key=lambda o:-o["ipsae_h"])
json.dump(scored, open(f"{R}/runs/wave6/scored_h.json","w"), indent=1)
print(f"scored {len(scored)}/{len(rows)} human")
print(f"{'variant':34s} {'ipSAE_h':>8} {'iptm_h':>7}")
for o in scored[:40]:
    print(f"{o['variant']:34s} {o['ipsae_h']:8.3f} {str(round(o['iptm_h'],3)) if o['iptm_h'] else '-':>7}")
import numpy as np
v=np.array([o["ipsae_h"] for o in scored])
print(f"\nipSAE_h dist: max={v.max():.3f} p75={np.percentile(v,75):.3f} median={np.median(v):.3f} min={v.min():.3f}")
# binders = ipSAE_h in the strong band; pick top-30 or >0.55 (calibrate vs v10 after)
binders=[o["variant"] for o in scored if o["ipsae_h"]>=0.55][:36]
open(f"{R}/runs/wave6/binders_h.txt","w").write("\n".join(binders)+"\n")
print(f"\nselected {len(binders)} human-binders (ipSAE_h>=0.55, cap 36) -> binders_h.txt")
