#!/usr/bin/env python3
"""Orthogonal pKa SIGN-CONSENSUS for the pH switch: pKAI (ML engine, PROPKA-independent) version of
switch_sharpness. Per binder His: pKa complex (binder+target) vs apo (binder only); dpKa<0 and
peak Dn_H>0 => bind-7.4/release-6.0 (agrees with PROPKA). Guards against single-engine pKa artifacts
(lesson: ph-switch triple-confirm). Usage: pkai_sharpness.py <names_seqs.csv> [out.tsv]"""
import sys, os, csv, glob, subprocess, json, tempfile
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
from Bio.SeqUtils import seq1
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PKAIPY=f"{R}/envs/pka/bin/python"
SEARCH=[f"{R}/runs/boltz_push_h/out/by_seed",f"{R}/runs/boltz_fin_h/out/by_seed",f"{R}/runs/boltz_coop_h/out/by_seed"]
if os.environ.get("BOLTZ_SEARCH"): SEARCH=[os.environ["BOLTZ_SEARCH"]]+SEARCH
def find_cif(n):
    for b in SEARCH:
        p=f"{b}/{n}/s0/{n}_model_0.cif"
        if os.path.exists(p): return p
    g=glob.glob(f"{R}/runs/boltz_*_h/out/by_seed/{n}/s0/{n}_model_0.cif"); return g[0] if g else None
def cseq(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
def idf(a,b):
    n=min(len(a),len(b)); return (sum(1 for i in range(n) if a[i]==b[i])/n) if n else 0
def f(pH,pKa): return 1.0/(1.0+10**(pH-pKa))
def relabel(bres,tres):
    st=Structure.Structure("s"); mo=Model.Model(0); st.add(mo); a=Chain.Chain("A")
    for i,r in enumerate(bres):
        rc=r.copy(); rc.id=(" ",i+1," "); a.add(rc)
    mo.add(a)
    if tres is not None:
        b=Chain.Chain("B")
        for j,r in enumerate(tres):
            rc=r.copy(); rc.id=(" ",j+1," "); b.add(rc)
        mo.add(b)
    return st
def pkai_his(pdb,chain):
    code=("import sys,json;from pkai.pKAI import pKAI;r=pKAI(sys.argv[1]);"
          "print('PKAI_JSON:'+json.dumps([[x[1],x[3]] for x in r if x[2]=='HIS' and x[0]==sys.argv[2]]))")
    out=subprocess.run([PKAIPY,"-c",code,pdb,chain],capture_output=True,text=True)
    for ln in out.stdout.splitlines():
        if ln.startswith("PKAI_JSON:"):
            return {int(k):float(v) for k,v in json.loads(ln[10:])}
    return {}
def analyze(name,seq):
    cif=find_cif(name)
    if not cif: return None
    model=MMCIFParser(QUIET=True).get_structure("x",cif)[0]
    chains=[c for c in model if cseq(c)]
    bch=max(chains,key=lambda c: idf(cseq(c),seq))
    if idf(cseq(bch),seq)<0.5: return None
    tch=max([c for c in chains if c.id!=bch.id],key=lambda c:len(cseq(c)))
    bres=[r for r in bch if r.id[0]==" "]; tres=[r for r in tch if r.id[0]==" "]
    with tempfile.TemporaryDirectory(dir=f"{R}/runs") as wd:
        io=PDBIO()
        io.set_structure(relabel(bres,tres)); cp=f"{wd}/c.pdb"; io.save(cp)
        io.set_structure(relabel(bres,None)); ap=f"{wd}/a.pdb"; io.save(ap)
        pc=pkai_his(cp,"A"); pa=pkai_his(ap,"A")
    per=[(h,pa[h],pc[h]) for h in sorted(set(pc)&set(pa))]
    grid=[6.0+0.05*i for i in range(29)]
    dnh=[sum(f(pH,a)-f(pH,c) for _,a,c in per) for pH in grid]
    integ=sum((dnh[i]+dnh[i+1])/2*0.05 for i in range(len(grid)-1))
    peak=max(dnh) if dnh else 0
    return dict(name=name,per=per,peak=peak,fc=10**integ)
rows=[]
for r in csv.DictReader(open(sys.argv[1])):
    name=r.get("name") or r.get("design_id"); a=analyze(name,r["sequence"])
    if a: rows.append(a)
rows.sort(key=lambda a:-a["peak"])
print(f"{'design':22} {'pKAI_peak':>9} {'pKAI_FC':>8}  per-His(apo->cplx)")
for a in rows:
    ph=" ".join(f"H{h}:{pa:.1f}->{pc:.1f}" for h,pa,pc in a["per"])
    tag="SHARP" if a["peak"]>=1.5 else ("ok" if a["peak"]>=1.0 else "weak")
    print(f"{a['name']:22} {a['peak']:9.2f} {a['fc']:7.0f}x [{tag}] {ph}")
if len(sys.argv)>2:
    with open(sys.argv[2],"w") as fh:
        fh.write("design_id\tpkai_peak_dnH\tpkai_fold_change\n")
        for a in rows: fh.write(f"{a['name']}\t{a['peak']:.2f}\t{a['fc']:.0f}\n")
    print(f"\nwrote {sys.argv[2]}")
