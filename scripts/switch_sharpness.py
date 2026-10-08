#!/usr/bin/env python3
"""pH-switch sharpness (Wyman/Hill) for bind-7.4 / release-6.0, the Challenge-2 direction.
For each binder His: pKa in complex (binder+target) vs apo (binder only) via PROPKA3.
  Dn_H(pH) = sum_His [ f(pH;pKa_apo) - f(pH;pKa_cplx) ],  f = protonated frac = 1/(1+10^(pH-pKa))
  FC = 10^( integral_{6.0}^{7.4} Dn_H dpH )   = Kd(6.0)/Kd(7.4) release factor  (>1 = releases at low pH, wanted)
  peak Dn_H over [6.0,7.4] = effective Hill n = SHARPNESS (physics: Dn_H=1 ~25x, 2 ~630x, 3 ~16000x).
Anti-cooperative His (clustered +charges) split their pKas -> peak stays ~1 even with 2 His (the obstacle).
Usage: switch_sharpness.py <names_seqs.csv> [out.tsv]"""
import sys, os, csv, glob, re, subprocess, tempfile
from Bio.PDB import MMCIFParser, PDBIO, Select, Structure, Model, Chain
from Bio.SeqUtils import seq1
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PROPKA=f"{R}/envs/pyrosetta/bin/propka3"
SEARCH=[f"{R}/runs/boltz_push_h/out/by_seed",f"{R}/runs/boltz_fin_h/out/by_seed",f"{R}/runs/boltz_coop_h/out/by_seed"]
def find_cif(n):
    for b in SEARCH:
        p=f"{b}/{n}/s0/{n}_model_0.cif"
        if os.path.exists(p): return p
    g=glob.glob(f"{R}/runs/boltz_*_h/out/by_seed/{n}/s0/{n}_model_0.cif"); return g[0] if g else None
def cseq(ch): return "".join(seq1(r.resname) if r.id[0]==" " else "" for r in ch)
def idf(a,b):
    n=min(len(a),len(b)); return (sum(1 for i in range(n) if a[i]==b[i])/n) if n else 0
def f(pH,pKa): return 1.0/(1.0+10**(pH-pKa))
def his_pkas(pdb,wd,chain):
    subprocess.run([PROPKA,os.path.basename(pdb)],capture_output=True,text=True,cwd=wd)
    pk=os.path.join(wd,os.path.basename(pdb).replace(".pdb",".pka"))
    if not os.path.exists(pk): return {}
    out={}
    for ln in open(pk):
        m=re.match(r"\s*HIS\s+(\d+)\s+(\w)\s+([-\d.]+)",ln)
        if m and m.group(2)==chain: out[int(m.group(1))]=float(m.group(3))
    return out
def relabel(binder_res,target_res):
    st=Structure.Structure("s"); mo=Model.Model(0); st.add(mo)
    a=Chain.Chain("A")
    for i,r in enumerate(binder_res):
        rc=r.copy(); rc.id=(" ",i+1," "); a.add(rc)
    mo.add(a)
    if target_res is not None:
        b=Chain.Chain("B")
        for j,r in enumerate(target_res):
            rc=r.copy(); rc.id=(" ",j+1," "); b.add(rc)
        mo.add(b)
    return st
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
        io.set_structure(relabel(bres,tres)); cpdb=os.path.join(wd,"c.pdb"); io.save(cpdb)
        io.set_structure(relabel(bres,None)); apdb=os.path.join(wd,"a.pdb"); io.save(apdb)
        pc=his_pkas(cpdb,wd,"A"); pa=his_pkas(apdb,wd,"A")
    per=[]
    for resid in sorted(set(pc)&set(pa)):
        per.append(dict(his=resid,apo=round(pa[resid],2),cplx=round(pc[resid],2),dpKa=round(pc[resid]-pa[resid],2)))
    # Dn_H(pH) over window
    grid=[6.0+0.05*i for i in range(29)]  # 6.0..7.4
    dnh=[sum(f(pH,p["apo"])-f(pH,p["cplx"]) for p in per) for pH in grid]
    integ=sum((dnh[i]+dnh[i+1])/2*(grid[i+1]-grid[i]) for i in range(len(grid)-1))
    peak=max(dnh) if dnh else 0
    fc=10**integ
    n_switch=sum(1 for p in per if p["dpKa"]<-0.1 and p["apo"]>5.3)  # titrates in-window, right direction
    return dict(name=name,nhis=len(per),fc=fc,peak=peak,integ=integ,n_switch=n_switch,per=per)
rows=[]
for r in csv.DictReader(open(sys.argv[1])):
    name=r.get("name") or r.get("design_id"); seq=r["sequence"]
    a=analyze(name,seq)
    if a is None: print(f"{name:20} NO COMPLEX/match"); continue
    rows.append(a)
rows.sort(key=lambda a:(-a["peak"],-a["fc"]))
print(f"{'design':20} {'nHis':>4} {'n_sw':>4} {'peakDnH':>7} {'FC(rel@6.0)':>11}  per-His(apo->cplx)")
for a in rows:
    ph=" ".join(f"H{p['his']}:{p['apo']}->{p['cplx']}" for p in a["per"])
    tag="SHARP" if a["peak"]>=1.5 else ("ok" if a["peak"]>=1.0 else "weak")
    print(f"{a['name']:20} {a['nhis']:4} {a['n_switch']:4} {a['peak']:7.2f} {a['fc']:10.0f}x  [{tag}] {ph}")
if len(sys.argv)>2:
    with open(sys.argv[2],"w") as fh:
        fh.write("design_id\tn_his\tn_switch\tpeak_dnH\tfold_change\tlog10FC\n")
        for a in rows: fh.write(f"{a['name']}\t{a['nhis']}\t{a['n_switch']}\t{a['peak']:.2f}\t{a['fc']:.0f}\t{a['integ']:.2f}\n")
    print(f"\nwrote {sys.argv[2]}")
