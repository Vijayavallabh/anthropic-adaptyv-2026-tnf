#!/usr/bin/env python3
"""Parametric multi-His linkage screen (generalizes wave10_screen): for every variant in a manifest,
read its Boltz human complex, compute aggregate Wyman linkage (collective_his_linkage) + human iptm, rank
most-negative first. Usage:
  envs/uplifting-esmc/bin/python scripts/mh_screen.py <manifest.tsv> <boltz_tag_glob> <out.tsv>
  e.g. ... runs/rebuild/multihis/manifest.tsv 'runs/boltz_rbmh*' runs/rebuild/multihis/screen.tsv
manifest lines: name<TAB>design<TAB>his_csv[<TAB>ncat]
"""
import os, glob, json, subprocess, tempfile, sys
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB import MMCIFParser, PDBIO
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PY=f"{R}/envs/uplifting-esmc/bin/python"
MAN, TAGGLOB, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
man=[l.split("\t") for l in open(MAN).read().splitlines() if l.strip()]

def find_cif(name):
    g=glob.glob(f"{R}/{TAGGLOB}/out/by_seed/{name}/s0/{name}_model_0.cif")
    return g[0] if g else None
def iptm(name):
    g=glob.glob(f"{R}/{TAGGLOB}/out/by_seed/{name}/s0/confidence_{name}_model_0.json")
    return round(json.load(open(g[0])).get("iptm",0),3) if g else None

def work(rec):
    name,design=rec[0],rec[1]; ncat=int(rec[3]) if len(rec)>3 else (rec[2].count(",")+1)
    cif=find_cif(name)
    if not cif: return None
    with tempfile.TemporaryDirectory() as wd:
        pdb=f"{wd}/{name}.pdb"
        try:
            m=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
            io=PDBIO(); io.set_structure(m); io.save(pdb)
            r=subprocess.run([PY,f"{R}/scripts/collective_his_linkage.py","--pdb",pdb,
                              "--binder-chain","A","--target-chain","B","--name",name],
                             capture_output=True,text=True,timeout=300)
            d=json.loads(r.stdout.strip().splitlines()[-1])
        except Exception as e:
            return (name,design,ncat,None,None,None)
    return (name,design,ncat,d.get("collective_linkage"),
            [(p["his"],round(p["dpKa"],2)) for p in d.get("per_his",[])], iptm(name))

out=[]
with ProcessPoolExecutor(max_workers=16) as ex:
    for res in ex.map(work, man):
        if res and res[3] is not None: out.append(res)
out.sort(key=lambda x:x[3])
with open(OUT,"w") as fh:
    fh.write("name\tdesign\tnHis\tlinkage\tiptm_h\tper_his\n")
    for name,design,ncat,link,per,ipt in out:
        fh.write(f"{name}\t{design}\t{ncat}\t{link}\t{ipt}\t{per}\n")
print(f"{'name':30s}{'nHis':>5}{'linkage':>9}{'iptm_h':>8}  per_his")
for name,design,ncat,link,per,ipt in out[:30]:
    print(f"{name:30s}{ncat:>5}{link:>9.3f}{(ipt or 0):>8.3f}  {per}")
strong=[o for o in out if o[3] is not None and o[3]<=-0.7 and o[2]>=2 and (o[5] or 0)>=0.7]
print(f"\nMULTI-His strong (<= -0.7) + cross-binding (iptm_h>=0.7): {len(strong)}")
for o in strong: print(f"  {o[0]} ({o[1]}) nHis={o[2]} link={o[3]} iptm_h={o[5]}")
print(f"total scored: {len(out)}")
