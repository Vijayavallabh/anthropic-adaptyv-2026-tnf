#!/usr/bin/env python3
"""WAVE-3B multi-His triple-confirm (interface-adapted designs). Same as wave3_triple but reads
runs/wave3b/switch_ranked.json (switch=True), derives His FROM SEQUENCE (runs/wave3b/variants.csv),
and uses the wv3bh Boltz folds. His switch-confirmed if PROPKA dpKa<0 AND xTB dE>0; design CONFIRMED if
>=2 His confirmed. Writes runs/wave3b/triple.json. Run with .venv/bin/python -I."""
import json, os, copy, subprocess, csv
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
VENV=f"{R}/.venv/bin/python"; MD2=f"{R}/envs/md2/bin/python"
OUT="runs/wave3b/triple_pdb"; os.makedirs(OUT, exist_ok=True)
P=MMCIFParser(QUIET=True); io=PDBIO()
ranked=json.load(open("runs/wave3b/switch_ranked.json"))
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/wave3b/variants.csv"))}
cands=[z for z in ranked if z["switch"]]
print(f"[wave3b triple] {len(cands)} interface-adapted switch candidates")
def cifpath(v):
    for k in range(4):
        p=f"runs/boltz_wv3bh{k}/out/by_seed/{v}/s0/{v}_model_0.cif"
        if os.path.exists(p): return p
    return None
def convert(v, his):
    c=cifpath(v)
    if not c: return None,None
    model=P.get_structure("x",c)[0]
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
    pdb=f"{OUT}/{v}.pdb"; io.set_structure(st); io.save(pdb)
    bres={r.id[1]:r for r in nb}; hp={}
    for h in his:
        r=bres.get(h)
        if r is None or r.resname not in("HIS","HID","HIE","HIP"): continue
        ha=[a.coord for a in r if a.name in("ND1","NE2")]
        if not ha or not cat: continue
        d,rn,kind=min((float(np.linalg.norm(x-cc)),rn,kind) for x in ha for rn,cc,kind in cat)
        if d<9.0: hp[h]=(rn,kind,round(d,2))   # only confirm His that are actually near a cation
    return pdb, hp
out=[]
print(f"{'variant':24}{'His':>4}{'part':>6}{'dpKa':>8}{'dE_xtb':>9}  conf")
for z in cands:
    v=z["variant"]; s=seqs.get(v,""); his=[i+1 for i,ch in enumerate(s) if ch=="H"]
    pdb,hp=convert(v,his)
    if not pdb or not hp: print(f"{v}: convert/no-near-His"); continue
    dpka={}
    r=subprocess.run([VENV,"scripts/collective_his_linkage.py","--pdb",pdb,"--binder-chain","A","--target-chain","B"],capture_output=True,text=True)
    try:
        j=json.loads(r.stdout.strip().splitlines()[-1])
        for ph in j["per_his"]: dpka[ph["his"]]=ph["dpKa"]
    except Exception as e: print(f"  [propka err {v}] {str(e)[:50]}")
    conf=0; per=[]
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
        print(f"{v:24}{h:>4}{str(rn)+kind[0]:>6}{sp:>8}{se:>9}  {'Y' if ok else '.'}")
    dc=conf>=2
    out.append(dict(variant=v,min=z["min"],dHM=z["dHM"],nhis=z["nhis"],ndistinct_cat=z["ndistinct_cat_h"],
                    confirmed_his=conf,design_confirmed=bool(dc),per_his=per))
    print(f"  -> {v}: {conf} confirmed  {'TRIPLE-OK' if dc else 'weak'}")
json.dump(out, open("runs/wave3b/triple.json","w"), indent=1)
ok=[o for o in out if o["design_confirmed"]]
print(f"\n=== TRACK-B TRIPLE-CONFIRMED multi-His switches: {len(ok)} ===")
for o in sorted(ok,key=lambda z:(-z["confirmed_his"],-z["min"])):
    print(f"  {o['variant']:24} min={o['min']:.3f} dHM={o['dHM']:.3f} his={o['nhis']} confirmed={o['confirmed_his']}")
