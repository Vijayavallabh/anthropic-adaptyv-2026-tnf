#!/usr/bin/env python3
"""Union-freeze rescue builder: freeze (human interface ∪ mouse interface ∪ His) so soluble
redesign preserves BOTH species' binding. Uses Boltz human cif as scaffold (binder=final seq w/ His).
Writes scaffold PDB (chain A=binder renum 1..L, B=target) + manifest for soluble-MPNN."""
import csv, os
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, PDBIO, Structure, Model, Chain
from Bio.SeqUtils import seq1
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
OUT=f"{R}/runs/push/rescue"; os.makedirs(f"{OUT}/scaffolds", exist_ok=True)
XH=f"{R}/runs/boltz_xrall_h/out/by_seed"; XM=f"{R}/runs/boltz_xrall_m/out/by_seed"
def xr(tok, base):  # find cif by token
    import glob
    h=glob.glob(f"{base}/*{tok}*/s0/*model_0.cif"); return h[0] if h else None
# design -> (final_seq_key, human_cif, mouse_cif)
final={r["design_id"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/submission/final_portfolio.csv"))}
conf ={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/submission/confirmed_all.csv"))}
A20V_name="gpu0__spec_cb1ep_cb1_epitope_0_model_13_1__A20V2H_vsK197_af2pred"
DES={
 "b4m2": (final["b4m2_2A87H_K101"], xr("A87E2H_vsK101",XH), xr("A87E2H_vsK101",XM)),
 "b1m12":(final["b1m12_2A67H_K112"], xr("A67V2H_vsK112",XH), xr("A67V2H_vsK112",XM)),
 "b1m7": (final["b1m7_2A57H_K183"], xr("A57K2H_vsK183",XH), xr("A57K2H_vsK183",XM)),
 "b2m0": (final["b2m0_1A67H_K84"], xr("A67W2H_vsK84",XH), xr("A67W2H_vsK84",XM)),
 "b1m9": (final["b1m9_1A55H_K93"], xr("A55R2H_vsK93",XH), xr("A55R2H_vsK93",XM)),
 "A20V": (conf[A20V_name], xr("A20V2H_vsK197",XH), xr("A20V2H_vsK197",XM)),
 "TNFpH01":(final["TNFpH01"], f"{R}/runs/boltz_core8_h/out/by_seed/TNFpH01/s0/TNFpH01_model_0.cif",
                              f"{R}/runs/boltz_core8_m/out/by_seed/TNFpH01/s0/TNFpH01_model_0.cif"),
}
CUT=float(os.environ.get("IFACE_CUT","7.0"))
HUMAN_ONLY=os.environ.get("HUMAN_ONLY","0")=="1"
ONLY=set(x for x in os.environ.get("ONLY","").split(",") if x)
if ONLY: DES={k:v for k,v in DES.items() if k in ONLY}
def load(path):
    par=MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return par.get_structure("x",path)[0]
def binder_and_iface(model, seq):
    chains=list(model)
    def cs(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
    best=None;bf=-1
    for ch in chains:
        s=cs(ch)
        if not s: continue
        ov=sum(1 for a,b in zip(s,seq) if a==b); f=ov/min(len(s),len(seq))
        if f>bf: bf=f;best=ch
    if bf<0.5: return None,None,bf
    bres=[r for r in best if r.id[0]==" "]
    tat=[a for ch in chains if ch.id!=best.id for r in ch if r.id[0]==" " for a in r if a.element!="H"]
    ns=NeighborSearch(tat)
    idx=set(i for i,r in enumerate(bres) if any(ns.search(a.coord,CUT) for a in r if a.element!="H"))
    return best,idx,bf
man=open(f"{OUT}/manifest.tsv","w")
man.write("design_id\tbackbone_pdb\tbinder_chain\tfixed_residues\this_pos\tfinal_seq\n")
for did,(seq,hc,mc) in DES.items():
    if not hc or not mc: print(f"{did}: MISSING complex h={bool(hc)} m={bool(mc)}"); continue
    hm=load(hc); bch,hidx,hf=binder_and_iface(hm,seq)
    if bch is None: print(f"{did}: no binder in human (f={hf:.2f})"); continue
    mm=load(mc); _,midx,mf=binder_and_iface(mm,seq)
    if midx is None: midx=set();
    his=set(i for i,c in enumerate(seq) if c=="H")
    union=hidx if HUMAN_ONLY else (hidx|midx)
    fixed=sorted(union|his)   # include His (scaffold already has His -> fixing keeps it)
    # write scaffold pdb: binder chain A renum 1..L, target chain B
    bres=[r for r in bch if r.id[0]==" "]
    news=Structure.Structure("s"); nm=Model.Model(0); news.add(nm)
    ca=Chain.Chain("A")
    for i,res in enumerate(bres):
        rc=res.copy(); rc.id=(" ",i+1," ");
        try: ca.add(rc)
        except Exception: pass
    nm.add(ca)
    tchains=[ch for ch in list(hm) if ch.id!=bch.id]
    cb=Chain.Chain("B"); j=1
    for ch in tchains:
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
    print(f"{did}: L={len(bres)} h_if={len(hidx)} m_if={len(midx)} union={len(union)} fixed(+His)={len(fixed)} free={len(bres)-len(fixed)} (hf={hf:.2f} mf={mf:.2f})")
man.close(); print(f"\nwrote {OUT}/manifest.tsv")
