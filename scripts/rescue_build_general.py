#!/usr/bin/env python3
"""General interface-frozen rescue builder: freeze (human iface ∪ mouse iface ∪ ALL His) so SolubleMPNN
re-sequences only the solvent-facing surface, preserving both-species binding AND the His switch.
Input CSV cols: name,seq,human_cif,mouse_cif. Writes <out>/scaffolds/<name>.pdb + <out>/manifest.tsv.
Usage: rescue_build_general.py <manifest_in.csv> <out_dir>"""
import csv, os, sys
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, PDBIO, Structure, Model, Chain
from Bio.SeqUtils import seq1
OUT=sys.argv[2]; os.makedirs(f"{OUT}/scaffolds", exist_ok=True)
CUT=float(os.environ.get("IFACE_CUT","7.0"))
def load(p):
    par=MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)
    return par.get_structure("x",p)[0]
def binder_and_iface(model, seq):
    chains=list(model)
    def cs(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
    best=None; bf=-1
    for ch in chains:
        s=cs(ch)
        if not s: continue
        ov=sum(1 for a,b in zip(s,seq) if a==b); f=ov/min(len(s),len(seq))
        if f>bf: bf=f; best=ch
    if bf<0.5: return None,None,bf
    bres=[r for r in best if r.id[0]==" "]
    tat=[a for ch in chains if ch.id!=best.id for r in ch if r.id[0]==" " for a in r if a.element!="H"]
    ns=NeighborSearch(tat)
    idx=set(i for i,r in enumerate(bres) if any(ns.search(a.coord,CUT) for a in r if a.element!="H"))
    return best,idx,bf
man=open(f"{OUT}/manifest.tsv","w")
man.write("design_id\tbackbone_pdb\tbinder_chain\tfixed_residues\this_pos\tfinal_seq\n")
for row in csv.DictReader(open(sys.argv[1])):
    did=row["name"]; seq=row["seq"]; hc=row["human_cif"]; mc=row.get("mouse_cif") or hc
    if not hc or not os.path.exists(hc): print(f"{did}: no human cif"); continue
    hm=load(hc); bch,hidx,hf=binder_and_iface(hm,seq)
    if bch is None: print(f"{did}: no binder match (hf={hf:.2f})"); continue
    midx=set()
    if mc and os.path.exists(mc):
        _,mi,mf=binder_and_iface(load(mc),seq); midx=mi or set()
    else: mf=-1
    his=set(i for i,c in enumerate(seq) if c=="H")
    fixed=sorted((hidx|midx)|his)
    bres=[r for r in bch if r.id[0]==" "]
    news=Structure.Structure("s"); nm=Model.Model(0); news.add(nm)
    ca=Chain.Chain("A")
    for i,res in enumerate(bres):
        rc=res.copy(); rc.id=(" ",i+1," ")
        try: ca.add(rc)
        except Exception: pass
    nm.add(ca)
    cb=Chain.Chain("B"); j=1
    for ch in [c for c in list(hm) if c.id!=bch.id]:
        for res in ch:
            if res.id[0]!=" ": continue
            rc=res.copy(); rc.id=(" ",j," "); j+=1
            try: cb.add(rc)
            except Exception: pass
    nm.add(cb)
    sp=f"{OUT}/scaffolds/{did}.pdb"; io=PDBIO(); io.set_structure(news); io.save(sp)
    fixed_str=" ".join(f"A{i+1}" for i in fixed)
    his_str=",".join(str(i+1) for i in sorted(his))
    man.write(f"{did}\t{sp}\tA\t{fixed_str}\t{his_str}\t{seq}\n")
    print(f"{did}: L={len(bres)} h_if={len(hidx)} m_if={len(midx)} union={len(hidx|midx)} fixed(+His)={len(fixed)} free={len(bres)-len(fixed)}")
man.close(); print(f"wrote {OUT}/manifest.tsv")
