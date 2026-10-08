#!/usr/bin/env python3
"""Release-at-low-pH switch score via PROPKA Wyman linkage (propka3 only; no PyRosetta).
Splits a binder(chain A):target(chain B) complex, runs propka3 on complex/binder/target,
and computes, per His group, the proton-linkage over [pH_lo, pH_hi]:
  n_H(pH)=1/(1+10^(pH-pKa));  dn_H = n_H(complex) - n_H(binder) - n_H(target)
  ddg_pH = RT*ln10 * integral_{lo}^{hi} dn_H dpH    (RT*ln10 = 1.419 kcal/mol @310K)
Convention (same as Challenge-1 obj1_switch): ddg_pH = dG_bind(hi) - dG_bind(lo).
  ddg_pH > 0 => tighter at LOW pH;  ddg_pH < 0 => WEAKER at low pH = RELEASE (what Challenge-2 wants).
His groups are keyed on (chain, resnum): a binder His and a target His that share a residue
number (common when the target is NOT renumbered above the binder) must not collide.
Each His's contribution is n_H(complex) - n_H(free), where "free" is the free binder for a
binder His and the free target for a target His. Also prints each His's pKa in complex vs free
(suppression in complex => release-competent).
Usage: ph_release_score.py <complex.pdb> [binder_chain=A] [target_chain=B] [pH_lo=6.0] [pH_hi=7.4]
"""
import sys,os,subprocess,tempfile,numpy as np
PROPKA=os.path.join(os.path.dirname(__file__),"..","envs","pyrosetta","bin","propka3")
RT_LN10=1.419
def run_propka(pdb,wd):
    b=os.path.basename(pdb)
    r=subprocess.run([os.path.abspath(PROPKA),b],cwd=wd,capture_output=True,timeout=600)
    pka=os.path.join(wd,b[:-4]+".pka")
    if not os.path.exists(pka):
        sys.exit(f"PROPKA-FAIL on {b}: rc={r.returncode} stderr={r.stderr.decode()[-300:]}")
    return pka
def parse_his(pkafile):
    # return {(chain,resnum): pKa} for HIS groups  (format: HIS <resnum> <chain> <pKa> <model-pKa>)
    out={}
    if not os.path.exists(pkafile): return out
    inblock=False
    for ln in open(pkafile):
        if ln.startswith("SUMMARY OF THIS PREDICTION"): inblock=True; continue
        if inblock:
            p=ln.split()
            if len(p)>=4 and p[0]=="HIS":
                try: out[(p[2],int(p[1]))]=float(p[3])
                except: pass
            if ln.strip().startswith("-") or ln.strip()=="":
                if out: break
    return out
def split_write(pdb,wd,bch,tch):
    A=[];B=[]
    for ln in open(pdb):
        if ln[:4]=="ATOM":
            c=ln[21]
            if c==bch: A.append(ln)
            elif c==tch: B.append(ln)
    open(os.path.join(wd,"cx.pdb"),'w').write("".join(A+B)+"END\n")
    open(os.path.join(wd,"bd.pdb"),'w').write("".join(A)+"END\n")
    open(os.path.join(wd,"tg.pdb"),'w').write("".join(B)+"END\n")
def nH(pka,pH): return 1.0/(1.0+10**(pH-pka))
def main():
    pdb=sys.argv[1]
    bch=sys.argv[2] if len(sys.argv)>2 else 'A'
    tch=sys.argv[3] if len(sys.argv)>3 else 'B'
    lo=float(sys.argv[4]) if len(sys.argv)>4 else 6.0
    hi=float(sys.argv[5]) if len(sys.argv)>5 else 7.4
    with tempfile.TemporaryDirectory() as wd:
        split_write(pdb,wd,bch,tch)
        hcx=parse_his(run_propka("cx.pdb",wd))
        hbd=parse_his(run_propka("bd.pdb",wd))
        htg=parse_his(run_propka("tg.pdb",wd))
        free={**hbd,**htg}   # (chain,resnum)->free pKa; binder/target chains distinct => no clash
        pHs=np.linspace(lo,hi,15)
        print(f"His pKa (complex vs free) — suppression in complex => release-competent:")
        tot=0.0
        for key in sorted(hcx):
            ch,rn=key
            pc=hcx[key]; pf=free.get(key)
            if pf is None:
                print(f"  His {ch}{rn}: complex={pc} free=None (missing, skip)"); continue
            dn=np.array([nH(pc,p)-nH(pf,p) for p in pHs])
            contrib=RT_LN10*np.trapezoid(dn,pHs)
            tot+=contrib
            side="binder" if ch==bch else "target"
            print(f"  His {ch}{rn} ({side}): pKa complex={pc:.2f} free={pf:.2f}  dpKa={pc-pf:+.2f}  contrib_ddg_pH={contrib:+.3f} kcal/mol")
        print(f"\nddg_pH (sum over His, [{lo},{hi}]) = {tot:+.3f} kcal/mol")
        print(f"  <0 => binds WEAKER at low pH = RELEASE (Challenge-2 target); magnitude ~ switch strength")
if __name__=="__main__": main()
