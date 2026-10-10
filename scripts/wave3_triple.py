#!/usr/bin/env python3
"""WAVE-3 multi-His triple-confirm. For each multi-His SWITCH candidate (runs/wave3/switch_ranked.json,
switch=True): convert its Boltz HUMAN fold to a clean complex PDB (binder=chain A renum 1.., target=chain B
renum 1..), then for EACH engineered His run PROPKA Wyman linkage (collective_his_linkage -> dpKa<0 = pKa
drops on binding = release at low pH) and GFN2-xTB (qm_switch_fair vs its nearest target cation -> dE>0 =
protonation weakens the contact). His switch-confirmed if dpKa<0 AND dE>0. Design CONFIRMED if >=2 His
confirmed (Baker: >=2 His-cation contacts drive sensitivity). Writes runs/wave3/triple.json.
Run with .venv/bin/python -I."""
import json, os, copy, subprocess
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
VENV=f"{R}/.venv/bin/python"; MD2=f"{R}/envs/md2/bin/python"
OUT="runs/wave3/triple_pdb"; os.makedirs(OUT, exist_ok=True)
P=MMCIFParser(QUIET=True); io=PDBIO()
ranked=json.load(open("runs/wave3/switch_ranked.json"))
mani=json.load(open("runs/wave3/manifest.json"))
cands=[z for z in ranked if z["switch"]]
print(f"[wave3 triple] {len(cands)} multi-His switch candidates")

def cifpath(v):
    for k in range(4):
        p=f"runs/boltz_wv3h{k}/out/by_seed/{v}/s0/{v}_model_0.cif"
        if os.path.exists(p): return p
    return None
def convert(v):
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
    # His (new pos==design pos) -> nearest target cation (new resnum, kind)
    his_partner={}
    bres={r.id[1]:r for r in nb}
    for hp in mani.get(v,{}).get("his",[]):
        r=bres.get(hp)
        if r is None: continue
        ha=[a.coord for a in r if a.name in("ND1","NE2")]
        if not ha or not cat: continue
        d,rn,kind=min((float(np.linalg.norm(h-cc)),rn,kind) for h in ha for rn,cc,kind in cat)
        his_partner[hp]=(rn,kind,round(d,2))
    return pdb, his_partner

out=[]
print(f"{'variant':22}{'His':>4}{'part':>6}{'dpKa':>8}{'dE_xtb':>9}  conf")
for z in cands:
    v=z["variant"]; pdb,hp=convert(v)
    if not pdb or not hp: print(f"{v}: convert fail"); continue
    # PROPKA once
    dpka={}
    r=subprocess.run([VENV,"scripts/collective_his_linkage.py","--pdb",pdb,"--binder-chain","A","--target-chain","B"],
                     capture_output=True,text=True)
    try:
        j=json.loads(r.stdout.strip().splitlines()[-1])
        for ph in j["per_his"]: dpka[ph["his"]]=ph["dpKa"]
    except Exception as e:
        print(f"  [propka err {v}] {str(e)[:60]} :: {r.stdout[-120:]} {r.stderr[-120:]}")
    confirmed=0; perhis=[]
    for h,(rn,kind,d) in sorted(hp.items()):
        dp=dpka.get(h)
        r2=subprocess.run([MD2,"scripts/qm_switch_fair.py","--pdb",pdb,"--his",str(h),"--partner",str(rn),"--partner-kind",kind],
                          capture_output=True,text=True)
        dE=None
        try: dE=json.loads(r2.stdout.strip().splitlines()[-1])["dE_switch"]
        except Exception: pass
        okp=(dp is not None and dp<0); okx=(dE is not None and dE>0); ok=okp and okx
        if ok: confirmed+=1
        perhis.append(dict(his=h,partner=rn,kind=kind,dist=d,dpKa=dp,dE=dE,ok=ok))
        sp=f"{dp:.2f}" if dp is not None else "-"; se=f"{dE:.3f}" if dE is not None else "-"
        print(f"{v:22}{h:>4}{str(rn)+kind[0]:>6}{sp:>8}{se:>9}  {'Y' if ok else '.'}")
    design_conf=confirmed>=2
    out.append(dict(variant=v,min=z["min"],dHM=z["dHM"],nhis=z["nhis"],ndistinct_cat=z["ndistinct_cat_h"],
                    confirmed_his=confirmed,design_confirmed=bool(design_conf),per_his=perhis))
    print(f"  -> {v}: {confirmed} His confirmed  {'TRIPLE-OK' if design_conf else 'weak'}")
json.dump(out, open("runs/wave3/triple.json","w"), indent=1)
ok=[o for o in out if o["design_confirmed"]]
print(f"\n=== TRIPLE-CONFIRMED multi-His switches (>=2 His dpKa<0 AND dE>0): {len(ok)} ===")
for o in sorted(ok,key=lambda z:(-z["confirmed_his"],-z["min"])):
    print(f"  {o['variant']:22} min={o['min']:.3f} dHM={o['dHM']:.3f} his={o['nhis']} confirmed={o['confirmed_his']}")
