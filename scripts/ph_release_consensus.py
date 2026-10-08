#!/usr/bin/env python3
"""Multi-class pKa SIGN-CONSENSUS for the pH-release switch (hardens objective #1).

PROPKA3 benchmarks as the WEAKEST His pKa predictor (JCIM 2023, His R^2~0.05), so a ddg_pH<0 from
PROPKA alone is thin evidence. This adds an ORTHOGONAL class -- pKAI (GNN trained to Poisson-Boltzmann
PypKa pKas) -- and reports whether both methods agree on the SIGN of the Wyman release linkage:
  n_H(pH)=1/(1+10^(pH-pKa));  dn_H = n_H(complex) - n_H(free)   [free=binder for A-His, target for B-His]
  ddg_pH = RT*ln10 * integral_{6.0}^{7.4} sum_His dn_H dpH       (RT*ln10=1.419 @310K)
  ddg_pH < 0  => weaker binding at low pH = RELEASE (what Challenge-2 wants).
Together with the two-state MD (geometric class) this is a 3-class consensus on the decisive objective.
Usage: ph_release_consensus.py <name_TAB_pdb.tsv> <out.tsv> [--model pKAI|pKAI+]
PROPKA leg runs via scripts/ph_release_score.py (envs/pyrosetta); pKAI leg via envs/pka.
"""
import sys, os, csv, math, re, subprocess, contextlib, argparse
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
TMP = REPO/"runs/tmp/consensus"; TMP.mkdir(parents=True, exist_ok=True)
RT_LN10 = 1.419
HIS_VARIANTS = {"HIP","HIE","HID","HSP","HSE","HSD"}
PH_LO, PH_HI = 6.0, 7.4

def extract_chain(pdb, chain, out):
    with open(out,"w") as fh:
        for l in open(pdb):
            if l.startswith("ATOM") and l[21]==chain: fh.write(l)
        fh.write("TER\nEND\n")

def normalize_his(pdb):
    txt=open(pdb).read()
    if not any(v in txt for v in HIS_VARIANTS): return pdb
    out=str(TMP/(Path(pdb).stem+"_hn.pdb"))
    with open(out,"w") as fh:
        for l in txt.splitlines(keepends=True):
            if l.startswith(("ATOM","HETATM")) and l[17:20].strip() in HIS_VARIANTS: l=l[:17]+"HIS"+l[20:]
            fh.write(l)
    return out

def run_pkai(pdb, model):
    from pkai.pKAI import pKAI
    with open(os.devnull,"w") as dn, contextlib.redirect_stdout(dn), contextlib.redirect_stderr(dn):
        res = pKAI(normalize_his(pdb), model_name=model, device="cpu")
    return {(str(c),int(rn)):float(pk) for (c,rn,resn,pk) in res if resn=="HIS"}  # His groups only

def nH(pka,pH): return 1.0/(1.0+10**(pH-pka))

def pkai_ddg(pdb, bchain, tchain, model):
    tgt=str(TMP/"t.pdb"); bnd=str(TMP/"b.pdb")
    extract_chain(pdb,tchain,tgt); extract_chain(pdb,bchain,bnd)
    cx=run_pkai(pdb,model); fb=run_pkai(bnd,model); ft=run_pkai(tgt,model)
    groups=[]  # (pka_complex, pka_free) per His, chain-aware free reference
    for (c,rn),pk in cx.items():
        free = fb.get((c,rn)) if c==bchain else ft.get((c,rn))
        if free is not None: groups.append((pk,free))
    if not groups: return None, 0
    pHs=[PH_LO+(PH_HI-PH_LO)*i/40 for i in range(41)]
    integ=sum(0.5*(sum(nH(c,pHs[k])-nH(f,pHs[k]) for c,f in groups)
                   +sum(nH(c,pHs[k+1])-nH(f,pHs[k+1]) for c,f in groups))*(pHs[k+1]-pHs[k])
              for k in range(len(pHs)-1))
    return RT_LN10*integ, len(groups)

def propka_ddg(pdb, bchain, tchain):
    try:
        out=subprocess.run(["uv","run","python","scripts/ph_release_score.py",pdb,bchain,tchain,
                            str(PH_LO),str(PH_HI)], cwd=REPO, capture_output=True, text=True, timeout=600).stdout
        m=re.search(r"ddg_pH \(sum.*?=\s*(-?[0-9.]+)", out)
        return float(m.group(1)) if m else None
    except Exception: return None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("tsv"); ap.add_argument("out")
    ap.add_argument("--model",default="pKAI",choices=["pKAI","pKAI+"])
    ap.add_argument("--binder",default="A"); ap.add_argument("--target",default="B")
    a=ap.parse_args()
    rows=list(csv.DictReader(open(a.tsv),delimiter="\t"))
    res=[]; agree=both=0
    for r in rows:
        name,pdb=r["name"],r["pdb"]
        dp=propka_ddg(pdb,a.binder,a.target)
        dk,ng=pkai_ddg(pdb,a.binder,a.target,a.model)
        sa = (dp is not None and dk is not None) and ((dp<0)==(dk<0))
        rb = (dp is not None and dk is not None) and dp<0 and dk<0
        agree += int(sa); both += int(rb)
        res.append(dict(name=name, ddg_propka=round(dp,3) if dp is not None else "",
                        ddg_pkai=round(dk,3) if dk is not None else "", n_his=ng,
                        sign_agree=sa, release_both=rb))
    with open(a.out,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["name","ddg_propka","ddg_pkai","n_his","sign_agree","release_both"])
        w.writeheader(); [w.writerow(x) for x in res]
    print(f"{'name':44} {'PROPKA':>8} {'pKAI':>8} {'sign?':>6} {'rel2?':>6}")
    for x in res: print(f"{x['name'][:44]:44} {str(x['ddg_propka']):>8} {str(x['ddg_pkai']):>8} {str(x['sign_agree']):>6} {str(x['release_both']):>6}")
    print(f"\nsign-agreement PROPKA vs pKAI: {agree}/{len(rows)} | release-confirmed by BOTH: {both}/{len(rows)}")
    print(f"wrote {a.out}")

if __name__=="__main__": main()
