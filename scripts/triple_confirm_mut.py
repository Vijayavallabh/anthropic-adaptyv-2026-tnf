#!/usr/bin/env python3
"""Standalone multi-His triple-confirm for ad-hoc mutants (His-headroom test). Same gate/tools as
wave4b_triple.py: per His within 9A of a target cation, confirm if PROPKA dpKa<0 AND GFN2-xTB dE_switch>0;
design CONFIRMED if >=2 His confirmed. Reads a Boltz human-complex cif. .venv/bin/python -I.
Usage: triple_confirm_mut.py <name> <cif> <binder_seq>"""
import json, os, copy, subprocess, sys
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
VENV=f"{R}/.venv/bin/python"; MD2=f"{R}/envs/md2/bin/python"
OUT="runs/hhead/triple_pdb"; os.makedirs(OUT, exist_ok=True)
P=MMCIFParser(QUIET=True); io=PDBIO()
name, cif, seq = sys.argv[1], sys.argv[2], sys.argv[3]
his=[i+1 for i,ch in enumerate(seq) if ch=="H"]
model=P.get_structure("x",cif)[0]
chs=sorted(model,key=lambda ch:sum(1 for _ in ch)); binder,target=chs[0],chs[-1]
st=Structure.Structure("x"); md=Model.Model(0); st.add(md)
nb=Chain.Chain("A"); i=1
for res in binder:
    r=copy.deepcopy(res); r.id=(' ',i,' '); nb.add(r); i+=1
nt=Chain.Chain("B"); j=1; cat=[]
for res in target:
    r=copy.deepcopy(res); r.id=(' ',j,' '); nt.add(r)
    if r.resname=="LYS":
        for a in r:
            if a.name=="NZ": cat.append((j,a.coord,"LYS"))
    elif r.resname=="ARG":
        for a in r:
            if a.name in("NH1","NH2","NE"): cat.append((j,a.coord,"ARG"))
    j+=1
md.add(nb); md.add(nt)
pdb=f"{OUT}/{name}.pdb"; io.set_structure(st); io.save(pdb)
bres={r.id[1]:r for r in nb}; hp={}
for h in his:
    r=bres.get(h)
    if r is None or r.resname not in("HIS","HID","HIE","HIP"): continue
    ha=[a.coord for a in r if a.name in("ND1","NE2")]
    if not ha or not cat: continue
    d,rn,kind=min((float(np.linalg.norm(x-cc)),rn,kind) for x in ha for rn,cc,kind in cat)
    if d<9.0: hp[h]=(rn,kind,round(d,2))
dpka={}
r=subprocess.run([VENV,"scripts/collective_his_linkage.py","--pdb",pdb,"--binder-chain","A","--target-chain","B"],capture_output=True,text=True)
try:
    jj=json.loads(r.stdout.strip().splitlines()[-1])
    for ph in jj["per_his"]: dpka[ph["his"]]=ph["dpKa"]
except Exception as e: print(f"[propka err] {str(e)[:80]}")
conf=0; per=[]
print(f"{'His':>4}{'part':>7}{'dist':>6}{'dpKa':>8}{'dE_xtb':>9}  conf")
for h,(rn,kind,d) in sorted(hp.items()):
    dp=dpka.get(h)
    r2=subprocess.run([MD2,"scripts/qm_switch_fair.py","--pdb",pdb,"--his",str(h),"--partner",str(rn),"--partner-kind",kind],capture_output=True,text=True)
    dE=None
    try: dE=json.loads(r2.stdout.strip().splitlines()[-1])["dE_switch"]
    except Exception: pass
    okp=(dp is not None and dp<0); okx=(dE is not None and dE>0); ok=okp and okx
    if ok: conf+=1
    per.append(dict(his=h,partner=rn,kind=kind,dist=d,dpKa=dp,dE=dE,ok=ok))
    sp=f"{dp:.2f}" if dp is not None else "-"; se=f"{dE:.3f}" if dE is not None else "-"
    print(f"{h:>4}{str(rn)+kind[0]:>7}{d:>6}{sp:>8}{se:>9}  {'Y' if ok else '.'}")
print(f"-> {name}: {conf} confirmed His ({'TRIPLE-OK' if conf>=2 else 'weak'}); His positions={his}")
json.dump(dict(name=name,confirmed_his=conf,per_his=per,his=his),open(f"runs/hhead/{name}_triple.json","w"),indent=1)
