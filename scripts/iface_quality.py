#!/usr/bin/env python3
"""Geometric interface-quality / affinity proxy (EGFR post-mortem correlates: contacts, buried SASA,
H-bonds) on binder(A):target(B) Boltz complexes. No PyRosetta. Usage: iface_quality.py <dir-of-cifs>"""
import sys, glob, os
from Bio.PDB import MMCIFParser, NeighborSearch, ShrakeRupley
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
DIR=sys.argv[1] if len(sys.argv)>1 else f"{R}/runs/boltz_push_h/out/by_seed"
sr=ShrakeRupley()
def analyze(cif):
    p=MMCIFParser(QUIET=True); m=p.get_structure("x",cif)[0]
    chains={c.id:c for c in m}
    A=chains.get("A"); B=chains.get("B")
    if A is None or B is None:
        cs=list(m); A,B=cs[0],cs[1]
    ares=[r for r in A if r.id[0]==" "]; bres=[r for r in B if r.id[0]==" "]
    batoms=[a for r in bres for a in r if a.element!="H"]
    ns=NeighborSearch(batoms)
    # interface contacts: binder atom within 4.5A of target atom
    contacts=0; iface_res=set(); hbonds=0
    NO=set(("N","O"))
    for r in ares:
        hit=False
        for a in r:
            if a.element=="H": continue
            near=ns.search(a.coord,4.5)
            if near:
                contacts+=len(near); hit=True
                if a.element in NO:
                    for b in near:
                        if b.element in NO and (a-b)<3.5: hbonds+=1
        if hit: iface_res.add(r.id[1])
    # buried SASA = SASA(A alone)+SASA(B alone)-SASA(complex) over binder+target
    try:
        sr.compute(m,level="R")
        sasa_cplx=sum(getattr(r,"sasa",0) for ch in (A,B) for r in ch if r.id[0]==" ")
        from Bio.PDB import Structure,Model,Chain
        def iso(ch):
            s=Structure.Structure("i");mm=Model.Model(0);s.add(mm)
            c=Chain.Chain(ch.id)
            for r in ch:
                if r.id[0]==" ":
                    try:c.add(r.copy())
                    except:pass
            mm.add(c); sr.compute(s,level="R")
            return sum(getattr(r,"sasa",0) for r in c)
        bsa=(iso(A)+iso(B)-sasa_cplx)
    except Exception as e:
        bsa=float("nan")
    return len(iface_res),contacts,hbonds,round(bsa,0)
print(f"{'design':18} {'iface_res':>9} {'contacts':>8} {'hbonds':>6} {'buried_SASA':>11}")
rows=[]
for d in sorted(glob.glob(f"{DIR}/*")):
    cif=glob.glob(f"{d}/s0/*_model_0.cif")
    if not cif: continue
    nm=os.path.basename(d)
    ir,ct,hb,bsa=analyze(cif[0]); rows.append((nm,ir,ct,hb,bsa))
for nm,ir,ct,hb,bsa in sorted(rows,key=lambda x:-x[4] if x[4]==x[4] else 0):
    print(f"{nm[:18]:18} {ir:>9} {ct:>8} {hb:>6} {bsa:>11.0f}")
