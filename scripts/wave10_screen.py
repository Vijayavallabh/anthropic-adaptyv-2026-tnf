#!/usr/bin/env python3
"""WAVE-10 screen: for every multi-His variant, read its Boltz human complex, compute the aggregate
Wyman linkage (collective_his_linkage) and human iptm, and rank. Most-negative linkage = strongest
correct-direction switch. Multi-His variants with strong linkage (comparable to the c19 flagship -1.274
or the v11 fold-robust members) are the 2nd-flagship candidates -> go to Chai + xTB both-fold triple-confirm.
Run: envs/uplifting-esmc/bin/python scripts/wave10_screen.py
"""
import os, glob, json, csv, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB import MMCIFParser, PDBIO
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PY=f"{R}/envs/uplifting-esmc/bin/python"
man=[l.split("\t") for l in open(f"{R}/runs/wave10/manifest.tsv").read().splitlines() if l.strip()]

def find_cif(name):
    for g in range(8):
        c=glob.glob(f"{R}/runs/boltz_w10h{g}/out/by_seed/{name}/s0/{name}_model_0.cif")
        if c: return c[0]
    c=glob.glob(f"{R}/runs/boltz_w10h*/out/**/{name}_model_0.cif",recursive=True)
    return c[0] if c else None

def iptm(name):
    for g in range(8):
        j=glob.glob(f"{R}/runs/boltz_w10h{g}/out/by_seed/{name}/s0/confidence_{name}_model_0.json")
        if j: return round(json.load(open(j[0])).get("iptm",0),3)
    return None

def work(rec):
    name,design,hpos,ncat=rec[0],rec[1],rec[2],int(rec[3])
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
            return (name,design,ncat,None,None,None,str(e)[:40])
    return (name,design,ncat,d.get("collective_linkage"),
            [(p["his"],round(p["dpKa"],2)) for p in d.get("per_his",[])], iptm(name), "")

recs=[r for r in man]
out=[]
with ProcessPoolExecutor(max_workers=16) as ex:
    for res in ex.map(work, recs):
        if res: out.append(res)
out=[o for o in out if o[3] is not None]
out.sort(key=lambda x:(x[3] if x[3] is not None else 0))  # most negative first
with open(f"{R}/runs/wave10/screen.tsv","w") as fh:
    fh.write("name\tdesign\tnHis\tlinkage\tiptm_h\tper_his\n")
    for name,design,ncat,link,per,ipt,err in out:
        fh.write(f"{name}\t{design}\t{ncat}\t{link}\t{ipt}\t{per}\n")
print(f"{'name':28s}{'nHis':>5}{'linkage':>9}{'iptm_h':>8}  per_his(dpKa)")
for name,design,ncat,link,per,ipt,err in out[:30]:
    print(f"{name:28s}{ncat:>5}{link:>9.3f}{(ipt if ipt else 0):>8.3f}  {per}")
strong=[o for o in out if o[3] is not None and o[3]<=-0.7 and o[2]>=2 and (o[5] or 0)>=0.5]
print(f"\nMULTI-His strong-linkage (<= -0.7) cross-binding candidates: {len(strong)}")
for o in strong: print(f"  {o[0]} ({o[1]}) nHis={o[2]} linkage={o[3]} iptm_h={o[5]}")
print(f"\ntotal scored: {len(out)}  ; flagship bar = c19 -1.274 / v11 fold-robust members -0.76..-1.27")
