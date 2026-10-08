#!/usr/bin/env python3
"""Baker His-cation contact metric (bioRxiv 2025.09.29.678932): working pH-switches contain >=2
contacts between (binder) protonated-His and surrounding cationic residues (Lys/Arg/His) on either
chain. Interface His...target-cation = like-charge repulsion/release; binder His...binder-cation =
buried-core destabilization. Count both; Baker floor = 2 (best designs 8-11).
Geometry: His imidazole N (ND1/NE2) within CUT A of a cationic charge-group atom
(Lys NZ; Arg NE/NH1/NH2/CZ; His ND1/NE2). Chain id'd by seq match to the design.
Usage: his_cation_contacts.py <names_seqs.csv(name,sequence[,...])> [out.tsv]"""
import sys, os, csv, glob
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.SeqUtils import seq1
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
SEARCH=[f"{R}/runs/boltz_push_h/out/by_seed",f"{R}/runs/boltz_fin_h/out/by_seed",
        f"{R}/runs/boltz_coop_h/out/by_seed",f"{R}/runs/boltz_rescue_h/out/by_seed"]
CUTS=[5.0,6.0]
CAT_ATOMS={"LYS":{"NZ"},"ARG":{"NE","NH1","NH2","CZ"},"HIS":{"ND1","NE2"}}
HIS_N={"ND1","NE2"}
def find_cif(name):
    for base in SEARCH:
        p=f"{base}/{name}/s0/{name}_model_0.cif"
        if os.path.exists(p): return p
    g=glob.glob(f"{R}/runs/boltz_*_h/out/by_seed/{name}/s0/{name}_model_0.cif")
    return g[0] if g else None
def chainseq(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
def idfrac(a,b):
    n=min(len(a),len(b));
    return (sum(1 for i in range(n) if a[i]==b[i])/n) if n else 0
def analyze(name,seq):
    cif=find_cif(name)
    if not cif: return None
    model=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
    chains=[c for c in model if chainseq(c)]
    # binder chain = best seq-id match to design seq
    bch=max(chains,key=lambda c: idfrac(chainseq(c),seq))
    if idfrac(chainseq(bch),seq)<0.5: return None
    tchs=[c for c in chains if c.id!=bch.id]
    his=[r for r in bch if r.id[0]==" " and r.resname=="HIS"]
    # cationic residues on binder (excluding the switching His itself) and target
    def cats(ch): return [r for r in ch if r.id[0]==" " and r.resname in CAT_ATOMS]
    bcat=cats(bch); tcat=[r for ch in tchs for r in cats(ch)]
    allcat=[(r,"bind") for r in bcat]+[(r,"tgt") for r in tcat]
    res={c:{"tot":0,"tgt":0,"bind":0} for c in CUTS}
    perhis={}
    for h in his:
        hN=[a for a in h if a.name in HIS_N]
        hid=h.id[1]; perhis[hid]={c:[] for c in CUTS}
        for cr,side in allcat:
            if cr is h: continue
            cN=[a for a in cr if a.name in CAT_ATOMS[cr.resname]]
            if not cN or not hN: continue
            d=min((a-b) for a in hN for b in cN)
            for c in CUTS:
                if d<=c:
                    res[c]["tot"]+=1; res[c][side]+=1
                    perhis[hid][c].append((cr.resname,cr.id[1],side,round(d,1)))
    return dict(name=name,nhis=len(his),cif=os.path.relpath(cif,R),res=res,perhis=perhis)
rows=[]
inp=sys.argv[1]
for r in csv.DictReader(open(inp)):
    name=r.get("name") or r.get("design_id"); seq=r["sequence"]
    a=analyze(name,seq)
    if a is None: print(f"{name:20} NO COMPLEX / no chain match"); continue
    rows.append(a)
    c5=a["res"][5.0]; c6=a["res"][6.0]
    flag="BAKER-PASS" if c6["tot"]>=2 else ("weak" if c6["tot"]==1 else "NONE")
    det=[]
    for hid,d in a["perhis"].items():
        parts=d[6.0]
        if parts: det.append(f"H{hid}->"+",".join(f"{rn}{ri}{'(tgt)' if s=='tgt' else ''}@{dd}" for rn,ri,s,dd in parts))
    print(f"{name:20} nHis={a['nhis']} | contacts@5A={c5['tot']}(tgt{c5['tgt']}) @6A={c6['tot']}(tgt{c6['tgt']}) | {flag}")
    for x in det: print(f"    {x}")
out=sys.argv[2] if len(sys.argv)>2 else None
if out:
    with open(out,"w") as fh:
        fh.write("design_id\tn_his_binder\tcontacts_5A\tcontacts_5A_tgt\tcontacts_6A\tcontacts_6A_tgt\tbaker_class\n")
        for a in rows:
            c5=a["res"][5.0]; c6=a["res"][6.0]
            cl="PASS" if c6["tot"]>=2 else ("weak" if c6["tot"]==1 else "none")
            fh.write(f"{a['name']}\t{a['nhis']}\t{c5['tot']}\t{c5['tgt']}\t{c6['tot']}\t{c6['tgt']}\t{cl}\n")
    print(f"\nwrote {out}")
