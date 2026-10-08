#!/usr/bin/env python3
"""Rank Region-I (inter-subunit groove) RFd3->MPNN->af2ig designs by POSE-AT-EPITOPE, not pae.
Calibration (reference-anthropic-tnf-dataset-regionI-recipe): af2ig pae is pessimistic/unreliable
for TNF affinity ranking, but fold self-consistency (rmsd<2, plddt_binder>80) and whether the binder
docks at the validated epitope ARE trustworthy. Epitope = V74/L75/I97 (protomer1) + W114/Y115 (protomer2),
which in the renumbered af2ig output (crystal = output+5) is chain B {69,70,92} + chain C {109,110}.
Usage: rank_regionI_pose.py <af2ig_dir> [out_csv]"""
import sys, os, glob, csv, numpy as np
AF=sys.argv[1] if len(sys.argv)>1 else "runs/regionI_af2ig"
OUTCSV=sys.argv[2] if len(sys.argv)>2 else f"{AF}/regionI_ranked.csv"
# af2ig MERGES the 2 target protomers into one chain (B), renumbered after the binder. So locate the
# epitope by POSITION WITHIN each 152-res protomer (robust to renumber + variable binder length):
#   protomer1 (first 152 target res): V74->pos69, L75->pos70, I97->pos92   (af2in protomer numbering)
#   protomer2 (next 152 target res):  W114->pos109, Y115->pos110
EPI_P1=[69,70,92]; EPI_P2=[109,110]; PROTLEN=152
EPI_EXPECT={69:'V',70:'L',92:'I',109:'W',110:'Y'}   # identity asserts (catch mapping drift)
RMSD_MAX=2.0; PLDDT_MIN=80.0; CONTACT=5.0; NEAR=8.0
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def parse_sc(path):
    rows={}
    if not os.path.exists(path): return rows
    with open(path) as fh:
        hdr=None
        for ln in fh:
            t=ln.split()
            if not t: continue
            if t[0]=="SCORE:" and hdr is None and "description" in ln: hdr=t[1:]; continue
            if t[0]=="SCORE:" and hdr:
                v=t[1:]
                d=dict(zip(hdr,v)); rows[d["description"]]=d
    return rows

# merge all shard scorefiles (or the merged one)
sc={}
for p in [f"{AF}/out.sc"]+glob.glob(f"{AF}/gpu*/out.sc"):
    sc.update(parse_sc(p))
print(f"scored designs: {len(sc)}")

# index output pdbs by description
pdbmap={}
for p in glob.glob(f"{AF}/gpu*/pdbs/*_af2pred.pdb")+glob.glob(f"{AF}/**/*_af2pred.pdb",recursive=True):
    base=os.path.basename(p).replace(".pdb","")
    # description in out.sc usually lacks the "gpuN__" prefix and the "_af2pred" suffix handling varies
    pdbmap[base]=p
    pdbmap[base.split("__",1)[-1]]=p

def load_chains(pdb):
    # returns {chain: [(resi,resname,[xyz...]) in file order]}
    ch={}; cur={}
    order=[]
    for ln in open(pdb):
        if not ln.startswith("ATOM"): continue
        c=ln[21]; resi=int(ln[22:26]); resn=ln[17:20].strip(); el=ln[76:78].strip() or ln[12:16].strip()[0]
        if el=='H': continue
        x,y,z=float(ln[30:38]),float(ln[38:46]),float(ln[46:54])
        key=(c,resi)
        if key not in cur:
            cur[key]=(resi,resn,[]); ch.setdefault(c,[]).append(cur[key])
        cur[key][2].append((x,y,z))
    return ch

def pose_score(pdb):
    ch=load_chains(pdb)
    if 'A' not in ch: return None
    A=np.array([xyz for _,_,atoms in ch['A'] for xyz in atoms])
    # target = all non-binder residues, in file order, split into 152-res protomers
    tgt=[]
    for c,reslist in ch.items():
        if c=='A': continue
        tgt.extend(reslist)
    if len(tgt)<PROTLEN+max(EPI_P2): return (0,0,None)   # incomplete target
    p1=tgt[:PROTLEN]; p2=tgt[PROTLEN:2*PROTLEN]
    picks=[(p1,EPI_P1),(p2,EPI_P2)]
    nc=0; nnear=0; mind=[]; bad=[]
    for prot,poss in picks:
        for pos in poss:
            if pos-1>=len(prot): continue
            resi,resn,atoms=prot[pos-1]
            if AA3.get(resn,'?')!=EPI_EXPECT[pos]: bad.append(f"{pos}:{resn}")
            E=np.array(atoms)
            d=np.sqrt(((A[:,None,:]-E[None,:,:])**2).sum(-1)).min()
            mind.append(d)
            if d<CONTACT: nc+=1
            if d<NEAR: nnear+=1
    pose_score.bad=bad
    return nc,nnear,(round(float(np.mean(mind)),2) if mind else None)
pose_score.bad=[]

rows=[]
for desc,d in sc.items():
    try:
        rmsd=float(d["binder_aligned_rmsd"]); plb=float(d["plddt_binder"]); pae=float(d["pae_interaction"])
    except: continue
    fold_ok = rmsd<RMSD_MAX and plb>PLDDT_MIN
    pdb=pdbmap.get(desc) or pdbmap.get("gpu2__"+desc)
    if pdb is None:
        cand=[v for k,v in pdbmap.items() if desc in k]
        pdb=cand[0] if cand else None
    ps=pose_score(pdb) if (pdb and fold_ok) else None
    nc,nnear,mind=(ps if ps else (0,0,None))
    rows.append(dict(design=desc,rmsd=round(rmsd,2),plddt_binder=round(plb,1),pae_int=round(pae,1),
                     fold_ok=fold_ok,epi_contacts=nc,epi_near=nnear,epi_mindist=mind,pdb=pdb))

# rank: fold-quality gate, then docks-at-epitope (contacts, then near), then plddt_binder
ranked=sorted([r for r in rows if r["fold_ok"]],
              key=lambda r:(-r["epi_contacts"],-r["epi_near"],-r["plddt_binder"]))
with open(OUTCSV,"w") as fh:
    w=csv.DictWriter(fh,fieldnames=["design","rmsd","plddt_binder","pae_int","epi_contacts","epi_near","epi_mindist","pdb"])
    w.writeheader()
    for r in ranked: w.writerow({k:r[k] for k in w.fieldnames})
nfold=sum(1 for r in rows if r["fold_ok"])
ndock=sum(1 for r in ranked if r["epi_contacts"]>=3)
print(f"fold-self-consistent (rmsd<2 & plddt_binder>80): {nfold}/{len(rows)}")
print(f"  of those, docking >=3/5 epitope residues: {ndock}")
print(f"\nTOP 15 by pose-at-epitope:")
print(f"{'design':48} {'rmsd':>5} {'plddtB':>6} {'pae':>5} {'epiC':>4} {'near':>4} {'mind':>5}")
for r in ranked[:15]:
    print(f"{r['design'][:48]:48} {r['rmsd']:>5} {r['plddt_binder']:>6} {r['pae_int']:>5} {r['epi_contacts']:>4} {r['epi_near']:>4} {str(r['epi_mindist']):>5}")
print(f"\nwrote {OUTCSV}")
