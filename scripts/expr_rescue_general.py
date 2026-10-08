#!/usr/bin/env python3
"""General minimal-edit expression rescue for the remaining FLAG designs.
Freeze interface (binder res <5A of target) + ALL His (switch/structural). Greedy fewest
surface substitutions from a neutral hydrophilic panel to reach instability<TARGET, guarding
against new N-glyc sequons / KKK|RRR runs / homopolymer runs. Preserves fold + switch region."""
import gemmi, numpy as np, re, json, os, sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
OUT=f"{R}/runs/push/exprrescue"; os.makedirs(OUT,exist_ok=True)
PANEL=list("QNST"); NMAX=10; TARGET=38.0
JOBS={
 "b1m9_1A55H_K93":"runs/boltz_push_h/out/by_seed/b1m9_1A55H_K93/s0/b1m9_1A55H_K93_model_0.cif",
 "b2m0_c20":"runs/boltz_coop_h/out/by_seed/b2m0_c20/s0/b2m0_c20_model_0.cif",
 "b2m0_c70":"runs/boltz_coop_h/out/by_seed/b2m0_c70/s0/b2m0_c70_model_0.cif",
 "b1m7_2A57H_K183_sh5":"runs/boltz_sharp2_6_h/out/by_seed/b1m7_2A57H_K183_sh5/s0/b1m7_2A57H_K183_sh5_model_0.cif",
}
def penalty(sq):
    pa=ProteinAnalysis(sq); ins=pa.instability_index(); gr=pa.gravy()
    return max(0,ins-TARGET)/10.0+max(0,gr-0.4)*2.0, ins, gr
def liability(sq,pos,aa):
    s=sq[:pos]+aa+sq[pos+1:]; lo=max(0,pos-2); hi=min(len(s),pos+3)
    if re.search(r"N[^P][ST]",s[max(0,pos-2):min(len(s),pos+3)]): return True   # N-glyc sequon
    if re.search(r"[KR]{3,}",s[lo:hi]): return True
    if re.search(r"(.)\1\1",s[lo:hi]): return True                              # homopolymer run
    return False
out_rows=[]
for did,cif in JOBS.items():
    cif=f"{R}/{cif}"
    s=gemmi.read_structure(cif); s.setup_entities(); m=s[0]
    A=[r for r in m['A']]; B=[r for r in m['B']]
    seq="".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in A)
    his={i+1 for i,c in enumerate(seq) if c=='H'}
    Bxyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for r in B for a in r if a.element.name!='H'])
    # representative side-chain atom per binder residue (CB, or CA for Gly) for burial proxy
    def rep(r):
        for nm in ("CB","CA"):
            a=r.find_atom(nm,'*')
            if a: return np.array([a.pos.x,a.pos.y,a.pos.z])
        return None
    reps=[rep(r) for r in A]
    iface=set(); buried=set()
    for i,r in enumerate(A):
        axyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in r if a.element.name!='H'])
        if axyz.size and np.sqrt(((axyz[:,None,:]-Bxyz[None,:,:])**2).sum(-1)).min()<5.0: iface.add(i+1)
        if reps[i] is not None:
            nb=sum(1 for j,q in enumerate(reps) if j!=i and q is not None and np.linalg.norm(reps[i]-q)<10.0)
            if nb>=10: buried.add(i+1)          # core residue -> freeze (not interface but fold-critical)
    frozen=iface|his|buried
    # never edit Pro/Gly (backbone-special: caps, turns)
    free=[i for i in range(len(seq)) if (i+1) not in frozen and seq[i] not in "PG"]
    cur=seq; edits=[]
    for _ in range(NMAX):
        p0,_,_=penalty(cur)
        if p0<=0: break
        best=None; bestp=p0
        for pos in free:
            if pos in [e[0] for e in edits]: continue
            for aa in PANEL:
                if aa==cur[pos] or liability(cur,pos,aa): continue
                p,_,_=penalty(cur[:pos]+aa+cur[pos+1:])
                if p<bestp-1e-9: bestp=p; best=(pos,aa,cur[:pos]+aa+cur[pos+1:])
        if best is None: break
        edits.append((best[0],best[1])); cur=best[2]
    pa0,pa1=ProteinAnalysis(seq),ProteinAnalysis(cur)
    ok=pa1.instability_index()<40 and pa1.gravy()<0.4 and cur.count('C')==0
    assert all(seq[p]==cur[p] for p in range(len(seq)) if (p+1) in frozen),"FROZEN VIOLATED"
    editstr=",".join(f"{seq[p]}{p+1}{aa}" for p,aa in edits)
    print(f"{did:30} iface={len(iface):2} buried={len(buried):2} His={sorted(his)} edits={len(edits)} [{editstr}]")
    print(f"{'':30} instab {pa0.instability_index():.1f}->{pa1.instability_index():.1f}  gravy {pa0.gravy():.2f}->{pa1.gravy():.2f}  {'PASS' if ok else 'STILL-FAIL'}")
    out_rows.append((did,seq,cur,len(edits),round(pa1.instability_index(),1),round(pa1.gravy(),2),ok,editstr,sorted(his)))
with open(f"{OUT}/rescued.csv","w") as fh:
    fh.write("name,sequence\n")
    for did,orig,new,ne,ins,gr,ok,es,his in out_rows:
        if ok and ne>0: fh.write(f"{did}_rx,{new}\n")
json.dump([dict(name=d+"_rx",orig=o,rescued=n,n_edit=ne,instab=i,gravy=g,ok=ok,edits=es,his=h) for d,o,n,ne,i,g,ok,es,h in out_rows],
          open(f"{OUT}/rescued.json","w"),indent=2)
print(f"\nwrote {OUT}/rescued.csv ({sum(1 for r in out_rows if r[6] and r[3]>0)} PASS) + rescued.json")
