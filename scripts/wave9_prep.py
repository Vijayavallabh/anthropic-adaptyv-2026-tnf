import os, glob, json, copy
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Structure, Model, Chain
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
OUT="runs/wave9"; os.makedirs(f"{OUT}/staged", exist_ok=True)
CONS=[1,6,27,39,60,77,85,93,98,107,123]; HIS_BIAS=1.5; CUT_HIS=9.0; CUT_IFACE=10.0
io=PDBIO()
# donor complexes: c19 (have pdb), c44 + b14 (from boltz cif)
donors={
 "des_3050_6_c19":"runs/hhead/triple_pdb/des_3050_6_c19_boltzchk.pdb",
 "des_4098_1_c44":"runs/boltz_wv4bm2/out/by_seed/des_4098_1_c44/s0/des_4098_1_c44_model_0.cif",
 "des_4098_1_b14":"runs/boltz_wv3bh3/out/by_seed/des_4098_1_b14/s0/des_4098_1_b14_model_0.cif",
}
mmp=MMCIFParser(QUIET=True); pp=PDBParser(QUIET=True)
def load_AB(path):
    m=(mmp if path.endswith(".cif") else pp).get_structure("x",path)[0]
    chs=sorted(m,key=lambda c:sum(1 for _ in c)); b,t=chs[0],chs[-1]
    st=Structure.Structure("x"); md=Model.Model(0); st.add(md)
    nb=Chain.Chain("A"); i=1
    for r in b:
        if r.id[0]!=' ': continue
        rr=copy.deepcopy(r); rr.id=(' ',i,' '); nb.add(rr); i+=1
    nt=Chain.Chain("B"); j=1
    for r in t:
        if r.id[0]!=' ': continue
        rr=copy.deepcopy(r); rr.id=(' ',j,' '); nt.add(rr); j+=1
    md.add(nb); md.add(nt); return st,nb,nt
bias={}; jobs=[]
for name,path in donors.items():
    if not os.path.exists(path): print(f"{name}: no donor"); continue
    st,nb,nt=load_AB(path)
    # cations (conserved) charged atoms
    cats=[]
    for r in nt:
        if r.id[1] in CONS and r.resname in("LYS","ARG"):
            for a in r:
                if a.name in("NZ","NH1","NH2","NE"): cats.append(a.coord)
    cats=np.array(cats)
    tatoms=np.array([a.coord for r in nt for a in r])
    iface=[]; hissites=[]
    for r in nb:
        cb=[a.coord for a in r if a.name in("CB","CA")]
        if not cb: continue
        cb=cb[0]
        if tatoms.size and np.linalg.norm(tatoms-cb,axis=1).min()<CUT_IFACE:
            iface.append(r.id[1])
            if cats.size and np.linalg.norm(cats-cb,axis=1).min()<CUT_HIS: hissites.append(r.id[1])
    if len(iface)<4: print(f"{name}: too few iface"); continue
    sp=f"{OUT}/staged/{name}.pdb"; io.set_structure(st); io.save(sp)
    bias[os.path.abspath(sp)]={f"A{p}":{"H":HIS_BIAS} for p in hissites}
    jobs.append((name, os.path.abspath(sp), " ".join(f"A{p}" for p in sorted(iface))))
    print(f"{name}: iface={len(iface)} His-biased={len(hissites)}")
json.dump(bias, open(f"{OUT}/bias.json","w"))
with open(f"{OUT}/jobs.tsv","w") as f:
    for n,p,red in jobs: f.write(f"{n}\t{p}\t{red}\n")
print(f"wrote {len(jobs)} jobs")
