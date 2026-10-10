#!/usr/bin/env python3
"""WAVE-6 stage C: fold-robust confirm (Boltz∩Chai human) + mouse binding + Wyman linkage for top binders.
Ranks 3-His separated-cation flagship candidates. third_party/LigandMPNN/.venv/bin/python -I."""
import os, glob, csv, json, subprocess
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
PYB=f"{R}/third_party/LigandMPNN/.venv/bin/python"
def boltz_h(n):
    for k in (0,1,2):
        f=f"runs/boltz_wv6h{k}/out/by_seed/{n}/s0/{n}_model_0.cif"
        if os.path.exists(f): return f
    return None
def chai_h(n):
    f=f"runs/chai_wv6ch/out/fast/{n}/seed_0/pred.model_idx_0.cif"
    return f if os.path.exists(f) else None
def mouse_ipsae(n):
    for k in (0,1,2):
        b=f"runs/boltz_wv6m{k}/out/by_seed/{n}/s0"; cif=f"{b}/{n}_model_0.cif"
        if not os.path.exists(cif): continue
        out=cif.replace(".cif","_10_10.txt")
        if not os.path.exists(out):
            npz=glob.glob(f"{b}/pae_{n}_model_0.npz")
            if npz: subprocess.run([PYB,f"{R}/tools/ipsae.py",os.path.basename(npz[0]),f"{n}_model_0.cif","10","10"],capture_output=True,cwd=b)
        if os.path.exists(out):
            a=[float(l.split()[5]) for l in open(out) if len(l.split())>6 and l.split()[4]=="asym"]
            if a: return min(a)
    return None
def okset(p):
    try:
        j=json.load(open(p)); return set((r["his"],r.get("partner")) for r in j["per_his"] if r.get("ok"))
    except: return None
def linkage(pdb):
    try:
        r=subprocess.run([PYB,"scripts/collective_his_linkage.py","--pdb",pdb,"--binder-chain","A","--target-chain","B"],capture_output=True,text=True,timeout=600)
        return json.loads(r.stdout.strip().splitlines()[-1]).get("collective_linkage")
    except: return None
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open("runs/wave6/binders_top.csv"))}
sc={o["variant"]:o for o in json.load(open("runs/wave6/scored_h.json"))}
res=[]
for n,seq in seqs.items():
    bc=boltz_h(n); cc=chai_h(n)
    if not bc or not cc: 
        print(f"{n}: missing fold (boltz={bool(bc)} chai={bool(cc)})",flush=True); continue
    subprocess.run([PYB,"-I","scripts/triple_confirm_mut.py",n+"__B",bc,seq],capture_output=True,text=True,timeout=2400)
    subprocess.run([PYB,"-I","scripts/triple_confirm_mut.py",n+"__C",cc,seq],capture_output=True,text=True,timeout=2400)
    b=okset(f"runs/hhead/{n}__B_triple.json"); c=okset(f"runs/hhead/{n}__C_triple.json")
    both=(b&c) if (b is not None and c is not None) else set()
    cats=sorted(set(p for (h,p) in both if p))
    mip=mouse_ipsae(n); hip=sc.get(n,{}).get("ipsae_h")
    pdb=f"runs/hhead/triple_pdb/{n}__B.pdb"
    link=linkage(pdb) if os.path.exists(pdb) else None
    res.append(dict(variant=n,n_robust=len(both),cations=cats,robust=sorted(list(both)),
                    ipsae_h=hip,ipsae_m=mip,dHM=abs(hip-mip) if (hip and mip) else None,link=link))
    print(f"{n:30s} robust={len(both)} cats={cats} ipsae_h={hip} ipsae_m={mip} link={link}",flush=True)
res.sort(key=lambda x:(-x["n_robust"], -(x["link"] or 0)))
json.dump(res, open("runs/wave6/flagship_candidates.json","w"), indent=1)
print("\n=== WAVE6 FLAGSHIP RANKING ===",flush=True)
for r in res:
    m=" <<<FLAGSHIP(>=3)" if r["n_robust"]>=3 else (" (2)" if r["n_robust"]==2 else "")
    print(f"{r['variant']:30s} robust={r['n_robust']} cats={r['cations']} link={r['link']} ipsae_h={r['ipsae_h']} ipsae_m={r['ipsae_m']}{m}",flush=True)
print("WAVE6_CONFIRM_DONE",flush=True)
