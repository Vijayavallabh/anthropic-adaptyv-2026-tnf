#!/usr/bin/env python3
"""Rebuild-gauntlet screen: for each candidate, gather the cheap gates BEFORE spending MD.
  - both-fold switch SIGN: human Boltz linkage (from mh_screen.tsv) AND Chai linkage (computed here)
    must BOTH be negative (correct release direction) -> guards the des_2375_4 anti-switch failure mode.
  - cross-reactivity: min(human_iptm, mouse_iptm) >= 0.70.
  - EXPRESSION (user hard co-gate: switch AND makeability): ESMFold monomer pLDDT, Boltz-2 monomer
    pLDDT, cross-folder TM agreement (indexed-CA Kabsch, same seq), Protein-Sol solubility, total His load.
Only designs passing sign + cross-react + expression advance to the two-state MD His-causality gate.
Usage: envs/uplifting-esmc/bin/python scripts/gauntlet_screen.py
"""
import os, glob, json, csv, subprocess, tempfile
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB import MMCIFParser, PDBIO
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PY=f"{R}/envs/uplifting-esmc/bin/python"
NAMES=[l.strip() for l in open(f"{R}/runs/rebuild/gauntlet.csv").read().splitlines()[1:] if l.strip()]
NAMES=[n.split(",")[0] for n in NAMES]
SEQ={r["name"]:r["sequence"].split(":")[0] for r in csv.DictReader(open(f"{R}/runs/rebuild/gauntlet.csv"))}

# human linkage + iptm from the rebuild screen
HUM={}
for r in csv.DictReader(open(f"{R}/runs/rebuild/multihis/screen.tsv"),delimiter="\t"):
    HUM[r["name"]]=(float(r["linkage"]), float(r["iptm_h"] or 0))

def chai_cif(n):
    g=glob.glob(f"{R}/runs/chai_rbg_ch/out/fast/{n}/seed_0/pred.model_idx_0.cif"); return g[0] if g else None
def mouse_iptm(n):
    g=glob.glob(f"{R}/runs/boltz_rbgm*/out/by_seed/{n}/s0/confidence_{n}_model_0.json")
    return round(json.load(open(g[0])).get("iptm",0),3) if g else None
def boltz_mono(n):
    g=glob.glob(f"{R}/runs/rebuild/gaunt/boltz_mono/by_seed/{n}/s0/confidence_{n}_model_0.json")
    if not g: return None,None
    d=json.load(open(g[0])); cif=g[0].replace("confidence_","").replace(".json",".cif")
    return round(d.get("complex_plddt",d.get("plddt",0)),3), cif
def esm_pdb(n):
    p=f"{R}/runs/rebuild/gaunt/esmf/{n}.pdb"; return p if os.path.exists(p) else None
def esm_plddt(pdb):
    bf=[float(l[60:66]) for l in open(pdb) if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    if not bf: return None
    m=sum(bf)/len(bf); return round(m/100 if m>1.5 else m,3)

def ca_coords_pdb(pdb):
    out={}
    for l in open(pdb):
        if l[:4]=="ATOM" and l[12:16].strip()=="CA":
            out[int(l[22:26])]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return out
def ca_coords_cif(cif):
    m=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
    out={}
    for ch in m:
        for r in ch:
            if r.id[0]!=' ': continue
            if 'CA' in r: out[r.id[1]]=r['CA'].coord
    return out
def tmscore(a,b):
    """indexed-CA Kabsch TM-score; both monomers share sequence so residue index = correspondence."""
    keys=sorted(set(a)&set(b))
    if len(keys)<10: return None
    P=np.array([a[k] for k in keys]); Q=np.array([b[k] for k in keys])
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(Pc.T@Qc)
    d=np.sign(np.linalg.det(V@Wt)); D=np.diag([1,1,d])
    U=V@D@Wt
    Pr=Pc@U
    di=np.linalg.norm(Pr-Qc,axis=1)
    L=len(keys); d0=1.24*(max(L,19)-15)**(1/3)-1.8
    return round(float(np.mean(1/(1+(di/d0)**2))),3)

def chai_link(n):
    cif=chai_cif(n)
    if not cif: return None
    with tempfile.TemporaryDirectory() as wd:
        pdb=f"{wd}/{n}.pdb"
        try:
            m=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
            io=PDBIO(); io.set_structure(m); io.save(pdb)
            r=subprocess.run([PY,f"{R}/scripts/collective_his_linkage.py","--pdb",pdb,
                              "--binder-chain","A","--target-chain","B","--name",n],
                             capture_output=True,text=True,timeout=300)
            return round(json.loads(r.stdout.strip().splitlines()[-1]).get("collective_linkage"),3)
        except Exception as e:
            return None

_SOL={}
for r in csv.DictReader(open(f"{R}/runs/rebuild/gaunt/solub.csv")):
    try: _SOL[r["name"]]=round(float(r.get("solubility_score") or r.get("scaled_sol")),3)
    except Exception: pass
def solub(name): return _SOL.get(name)

def work(n):
    hl,hi=HUM.get(n,(None,None))
    cl=chai_link(n); mi=mouse_iptm(n)
    ep=esm_pdb(n); esm=esm_plddt(ep) if ep else None
    bp,bcif=boltz_mono(n)
    tm=None
    if ep and bcif and os.path.exists(bcif):
        try: tm=tmscore(ca_coords_pdb(ep), ca_coords_cif(bcif))
        except Exception: tm=None
    sol=solub(n)
    totH=SEQ[n].count("H")
    return dict(name=n,totH=totH,hum_link=hl,chai_link=cl,hum_iptm=hi,mouse_iptm=mi,
                esm_plddt=esm,boltz_mono_plddt=bp,cross_tm=tm,solub=sol)

res=[]
with ProcessPoolExecutor(max_workers=14) as ex:
    for r in ex.map(work, NAMES): res.append(r)

# gates
def gate(r):
    sign = (r["hum_link"] is not None and r["hum_link"]<0 and r["chai_link"] is not None and r["chai_link"]<0)
    xr   = (min(r["hum_iptm"] or 0, r["mouse_iptm"] or 0) >= 0.70)
    expr = ((r["esm_plddt"] or 0)>=0.70 and (r["boltz_mono_plddt"] or 0)>=0.70
            and (r["cross_tm"] or 0)>=0.60 and (r["solub"] or 0)>=0.50)
    return sign,xr,expr

print(f"{'name':28s}{'H':>3}{'hLink':>7}{'cLink':>7}{'hIp':>6}{'mIp':>6}{'esm':>6}{'bz':>6}{'TM':>6}{'sol':>6}  SIGN XR EXPR")
res.sort(key=lambda r:(r["hum_link"] or 0))
surv=[]
for r in res:
    s,x,e=gate(r)
    flag="".join(["S" if s else "-","X" if x else "-","E" if e else "-"])
    print(f"{r['name']:28s}{r['totH']:>3}{(r['hum_link'] or 0):>7.2f}{(r['chai_link'] or 0):>7.2f}"
          f"{(r['hum_iptm'] or 0):>6.2f}{(r['mouse_iptm'] or 0):>6.2f}{(r['esm_plddt'] or 0):>6.2f}"
          f"{(r['boltz_mono_plddt'] or 0):>6.2f}{(r['cross_tm'] or 0):>6.2f}{(r['solub'] or 0):>6.2f}   {flag}")
    if s and x and e: surv.append(r["name"])
json.dump(res,open(f"{R}/runs/rebuild/gaunt/screen.json","w"),indent=1)
print(f"\nSURVIVORS (sign+xreact+express -> MD): {len(surv)}")
for n in surv: print(f"  {n}")
open(f"{R}/runs/rebuild/gaunt/survivors.txt","w").write("\n".join(surv)+"\n")
print("wrote runs/rebuild/gaunt/screen.json + survivors.txt")
