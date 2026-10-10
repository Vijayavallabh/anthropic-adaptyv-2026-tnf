#!/usr/bin/env python3
"""WAVE-9 dual-metric screen: for top binders, per-His fold-robust (Boltz n Chai) + both-fold Wyman linkage
+ cross-species iptm. Flag new flagships (>=3 per-His robust OR both-fold linkage both < -0.9).
third_party/LigandMPNN/.venv/bin/python -I."""
import glob, json, os, subprocess, copy, csv
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
PYB=f"{R}/third_party/LigandMPNN/.venv/bin/python"
P=MMCIFParser(QUIET=True); io=PDBIO()
def bo_h(n):
    for k in(0,1,2):
        f=f"runs/boltz_wv9h{k}/out/by_seed/{n}/s0/{n}_model_0.cif"
        if os.path.exists(f): return f
def ch_h(n):
    f=f"runs/chai_wv9ch/out/fast/{n}/seed_0/pred.model_idx_0.cif"
    return f if os.path.exists(f) else None
def iptm(n,kind):
    for b in glob.glob(f"runs/boltz_wv9{kind}*/out/by_seed/{n}/s0"):
        cj=glob.glob(f"{b}/confidence_{n}_model_0.json")
        if cj:
            try: return round(json.load(open(cj[0]))["iptm"],3)
            except: pass
def okset(p):
    try:
        j=json.load(open(p)); return set((r["his"],r.get("partner")) for r in j["per_his"] if r.get("ok"))
    except: return None
def cif2pdb(n,cif,tag):
    m=P.get_structure("x",cif)[0]; chs=sorted(m,key=lambda c:sum(1 for _ in c)); b,t=chs[0],chs[-1]
    st=Structure.Structure("x"); md=Model.Model(0); st.add(md)
    nb=Chain.Chain("A"); i=1
    for r in b:
        if r.id[0]!=' ': continue
        rr=copy.deepcopy(r); rr.id=(' ',i,' '); nb.add(rr); i+=1
    nt=Chain.Chain("B"); j=1
    for r in t:
        if r.id[0]!=' ': continue
        rr=copy.deepcopy(r); rr.id=(' ',j,' '); nt.add(rr); j+=1
    md.add(nb); md.add(nt); os.makedirs("runs/wave9/lk",exist_ok=True)
    pdb=f"runs/wave9/lk/{n}_{tag}.pdb"; io.set_structure(st); io.save(pdb); return pdb
def linkage(pdb):
    try:
        r=subprocess.run([PYB,"scripts/collective_his_linkage.py","--pdb",pdb,"--binder-chain","A","--target-chain","B"],capture_output=True,text=True,timeout=300)
        return json.loads(r.stdout.strip().splitlines()[-1]).get("collective_linkage")
    except: return None
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/wave9/binders_top.csv"))}
res=[]
for n,seq in seqs.items():
    bc=bo_h(n); cc=ch_h(n)
    if not bc or not cc: continue
    subprocess.run([PYB,"-I","scripts/triple_confirm_mut.py",n+"__B",bc,seq],capture_output=True,timeout=2400)
    subprocess.run([PYB,"-I","scripts/triple_confirm_mut.py",n+"__C",cc,seq],capture_output=True,timeout=2400)
    b=okset(f"runs/hhead/{n}__B_triple.json"); c=okset(f"runs/hhead/{n}__C_triple.json")
    both=(b&c) if (b is not None and c is not None) else set()
    cats=sorted(set(p for (h,p) in both if p))
    lb=linkage(cif2pdb(n,bc,"B")); lc=linkage(cif2pdb(n,cc,"C"))
    r=dict(variant=n,n_robust=len(both),cations=cats,lk_b=lb,lk_c=lc,iptm_h=iptm(n,'h'),iptm_m=iptm(n,'m'))
    flag=(len(both)>=3) or (lb is not None and lc is not None and lb<-0.9 and lc<-0.9)
    r["FLAGSHIP"]=flag
    res.append(r)
    print(f"{n}: robust={len(both)} cats={cats} lk={lb}/{lc} iptm={r['iptm_h']}/{r['iptm_m']} {'<<<FLAGSHIP' if flag else ''}",flush=True)
res.sort(key=lambda x:(-x["n_robust"], (x["lk_b"] or 0)+(x["lk_c"] or 0)))
json.dump(res, open("runs/wave9/flagships.json","w"), indent=1)
print("\n=== WAVE9 (vs c19: 4 robust, -1.27/-1.12) ===",flush=True)
for r in res[:12]:
    print(f"{r['variant']}: robust={r['n_robust']} cats={r['cations']} lk={r['lk_b']}/{r['lk_c']} iptm={r['iptm_h']}/{r['iptm_m']} {'FLAGSHIP' if r['FLAGSHIP'] else ''}",flush=True)
print(f"\nNEW flagships: {[r['variant'] for r in res if r['FLAGSHIP']]}",flush=True)
print("WAVE9_CONFIRM_DONE",flush=True)
