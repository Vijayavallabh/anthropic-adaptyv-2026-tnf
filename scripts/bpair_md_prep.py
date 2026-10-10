#!/usr/bin/env python3
"""Prep brutal two-state switch MD for bpair His-switch candidates.
For each switch candidate (scored.json, switch_cand=True), convert its human Boltz fold cif->pdb and
measure His(ND1/NE2)->nearest target Lys NZ distance (the pair md_release_tnf.py tracks; it uses Lys NZ
only). Select up to 8 candidates that (a) keep cross-react (highest min-ipSAE) and (b) have a Lys within
CONTACT range (<7.5A) so the switch is measurable. Writes PDBs + runs/bpair_mature/md_manifest.json."""
import json, os
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
MD=f"{R}/runs/bpair_mature/md"; os.makedirs(MD, exist_ok=True)
CONTACT=7.5; NSEL=8

scored=json.load(open("runs/bpair_mature/scored.json"))
cands=[o for o in scored if o.get("switch_cand")]
cands.sort(key=lambda o:-o["min"])
def cifpath(v):
    for t in ("bpmh0","bpmh1","bpmh2","bpmh3"):
        p=f"runs/boltz_{t}/out/by_seed/{v}/s0/{v}_model_0.cif"
        if os.path.exists(p): return p
    return None
parser=MMCIFParser(QUIET=True); io=PDBIO()
info=[]
for o in cands:
    v=o["variant"]
    if "_H" not in v: continue
    hpos=int(v.split("_H")[-1])
    c=cifpath(v)
    if not c: continue
    s=parser.get_structure("x",c); model=s[0]
    his_at=[]; lysnz=[]; argcat=[]
    for ch in model:
        for res in ch:
            rn=res.id[1]
            if ch.id=="A" and rn==hpos and res.resname in ("HIS","HID","HIE","HIP"):
                his_at=[a.coord for a in res if a.name in ("ND1","NE2")]
            if ch.id=="B" and res.resname=="LYS":
                for a in res:
                    if a.name=="NZ": lysnz.append((a.coord,rn))
            if ch.id=="B" and res.resname=="ARG":
                for a in res:
                    if a.name in ("CZ","NH1","NH2"): argcat.append((a.coord,rn))
    if not his_at: continue
    dl=[(min(np.linalg.norm(h-nz) for h in his_at),rn) for nz,rn in lysnz]
    da=[(min(np.linalg.norm(h-cz) for h in his_at),rn) for cz,rn in argcat]
    dlys,klys=min(dl) if dl else (999,None)
    darg,karg=min(da) if da else (999,None)
    pdb=f"{MD}/{v}.pdb"; io.set_structure(s); io.save(pdb)
    info.append(dict(variant=v,his_pos=hpos,min_ipsae=o["min"],wt_min=o["wt_min"],
                     his_lysNZ=round(float(dlys),2),lys_res=klys,
                     his_argCZ=round(float(darg),2),arg_res=karg,pdb=pdb))

# select: Lys within contact range, highest cross-react first
contactable=[x for x in info if x["his_lysNZ"]<CONTACT]
contactable.sort(key=lambda x:-x["min_ipsae"])
sel=contactable[:NSEL]
for i,x in enumerate(sel): x["gpu"]=i
json.dump(sel, open("runs/bpair_mature/md_manifest.json","w"), indent=1)
# also dump ALL contactable (ranked) so later batches can MD ranks beyond the first NSEL
json.dump(contactable, open("runs/bpair_mature/md_manifest_all.json","w"), indent=1)

print(f"{'variant':20}{'his':>4}{'min':>7}{'hisLysNZ':>9}{'lys':>5}{'hisArgCZ':>9}{'arg':>5}  sel")
for x in sorted(info,key=lambda x:-x["min_ipsae"]):
    s="GPU%d"%x["gpu"] if x in sel else ("" if x["his_lysNZ"]<CONTACT else "noLys")
    print(f"{x['variant']:20}{x['his_pos']:>4}{x['min_ipsae']:>7.3f}{x['his_lysNZ']:>9.2f}{str(x['lys_res']):>5}{x['his_argCZ']:>9.2f}{str(x['arg_res']):>5}  {s}")
print(f"\nselected {len(sel)} for brutal MD (Lys within {CONTACT}A): "+", ".join(x['variant'] for x in sel))
