#!/usr/bin/env python3
"""WAVE-2 MD prep: from runs/bpair_mature2/scored.json take the DUAL-STRONG His variants, convert each one's
human Boltz fold (bpm2h*) cif->pdb, measure His(ND1/NE2)->nearest target Lys NZ (md_release tracks Lys), and
select up to NSEL with a Lys within CONTACT (<7.5A). Writes runs/bpair_mature2/md_manifest.json + PDBs.
Run with .venv/bin/python -I."""
import json, os, glob
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
MD=f"{R}/runs/bpair_mature2/md"; os.makedirs(MD, exist_ok=True)
CONTACT=7.5; NSEL=8

scored=json.load(open("runs/bpair_mature2/scored.json"))
cands=[o for o in scored if o.get("dual_strong")]
cands.sort(key=lambda o:-o["min"])
def cifpath(v):
    for t in ("bpm2h0","bpm2h1","bpm2h2","bpm2h3"):
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
    his_at=[]; lysnz=[]
    for ch in model:
        for res in ch:
            rn=res.id[1]
            if ch.id=="A" and rn==hpos and res.resname in ("HIS","HID","HIE","HIP"):
                his_at=[a.coord for a in res if a.name in ("ND1","NE2")]
            if ch.id=="B" and res.resname=="LYS":
                for a in res:
                    if a.name=="NZ": lysnz.append((a.coord,rn))
    if not his_at or not lysnz: continue
    dlys,klys=min((min(np.linalg.norm(h-nz) for h in his_at),rn) for nz,rn in lysnz)
    pdb=f"{MD}/{v}.pdb"; io.set_structure(s); io.save(pdb)
    info.append(dict(variant=v,his_pos=hpos,min_ipsae=o["min"],dHM=o["dHM"],his_lysNZ=round(float(dlys),2),lys_res=klys,pdb=pdb))

contactable=[x for x in info if x["his_lysNZ"]<CONTACT]
contactable.sort(key=lambda x:-x["min_ipsae"])
sel=contactable[:NSEL]
for i,x in enumerate(sel): x["gpu"]=i
json.dump(sel, open("runs/bpair_mature2/md_manifest.json","w"), indent=1)
json.dump(contactable, open("runs/bpair_mature2/md_manifest_all.json","w"), indent=1)
print(f"{'variant':18}{'his':>4}{'min':>7}{'dHM':>6}{'hisLysNZ':>9}{'lys':>5}  sel")
for x in sorted(info,key=lambda x:-x["min_ipsae"]):
    s="GPU%d"%x["gpu"] if x in sel else ("" if x["his_lysNZ"]<CONTACT else "noLys")
    print(f"{x['variant']:18}{x['his_pos']:>4}{x['min_ipsae']:>7.3f}{x['dHM']:>6.3f}{x['his_lysNZ']:>9.2f}{str(x['lys_res']):>5}  {s}")
print(f"\nselected {len(sel)} dual-strong for brutal MD (Lys<{CONTACT}A): "+", ".join(x['variant'] for x in sel))
