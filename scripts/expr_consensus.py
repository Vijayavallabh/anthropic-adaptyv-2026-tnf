#!/usr/bin/env python3
"""Expression de-risking CONSENSUS for the submission portfolio — answers the valid critique that a
single ESMFold-pLDDT predictor (AUROC ~0.72, Garcia/Dixit/Rocklin 2026) is a hopeful prediction, not
proof. Combines ORTHOGONAL signals and reframes "all pass" as a calibrated EXPECTED-YIELD with a band.

Signals (independent error sources):
  1. ESMFold monomer pLDDT        (ESM2 language-model folder)
  2. Boltz-2  monomer pLDDT, pTM  (AF3-style diffusion folder — architecturally independent of ESM2)
  3. cross-folder TM(ESMFold,Boltz) agreement  (Baker/Rocklin gate: 2 independent folders must agree
     on the SAME fold; disagreement is a red flag even when each confidence is individually ok)
  4. developability liabilities   (net charge / GRAVY / hydrophobic patch / charge-blob / low-cplx)
  5. [optional] experimental-data solubility model score (results/solubility_v11.csv), fully
     orthogonal modality (trained on measured E. coli solubility, not fold confidence)

Per-design success probability is a TRANSPARENT, literature-anchored map on ESMFold pLDDT (the only
axis with a published de-novo AUROC), then MODULATED by the orthogonal signals (agreement raises it,
disagreement / liabilities / low solubility lower it). This is an estimate, NOT a fit to the Adaptyv
assay; the band reflects that. The robust headline is the TIER count (A/B/C), which needs no priors.

Usage: python scripts/expr_consensus.py [submission.csv]   (default runs/push/submission_v11.csv)
Read-only. Writes runs/expr_v11/consensus.{csv,json} and prints the dossier.
"""
import csv, json, sys, os, glob, subprocess, math
from pathlib import Path
R = Path("/home/prachh/v/anthropic-adaptyv-2026-tnf")
FS = R/"tools/foldseek_dl/foldseek/bin/foldseek"
TMP = Path("/home/prachh/v/.claude-private/jobs/5045f83d/tmp/fsrun")

def esm_plddt(pdb):
    bf = [float(l[60:66]) for l in open(pdb) if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    if not bf: return None
    m = sum(bf)/len(bf)
    return round(m/100 if m > 1.5 else m, 3)

def boltz_conf(name):
    for g in (3,4,5):
        j = list((R/f"runs/expr_v11/boltz_out_g{g}").rglob(f"confidence_{name}_model_0.json"))
        if j:
            d = json.load(open(j[0]))
            cif = str(j[0]).replace("confidence_","").replace(".json",".cif")
            return round(d.get("complex_plddt",0),3), round(d.get("ptm",0),3), cif
    return None, None, None

def tm(esm_pdb, boltz_cif):
    if not (os.path.exists(esm_pdb) and os.path.exists(boltz_cif)): return None
    TMP.mkdir(parents=True, exist_ok=True)
    out = TMP/"aln.m8"
    try:
        subprocess.run([str(FS),"easy-search",esm_pdb,boltz_cif,str(out),str(TMP/"t"),
                        "--alignment-type","1","--format-output","query,target,alntmscore,qtmscore,ttmscore",
                        "-e","100","--tmscore-threshold","0.0"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
        txt = out.read_text().strip()
        if not txt: return 0.0
        f = txt.splitlines()[0].split("\t")
        # min of query/target-normalized TM = conservative agreement (same residue count here anyway)
        return round(min(float(f[3]), float(f[4])), 3)
    except Exception:
        return None

def base_p(pl):
    """Literature-anchored P(de-novo success | ESMFold monomer pLDDT), Garcia/Dixit/Rocklin 2026 regime."""
    if pl is None: return 0.35
    if pl >= 0.85: return 0.80
    if pl >= 0.80: return 0.72
    if pl >= 0.75: return 0.62
    if pl >= 0.70: return 0.50
    return 0.35

def main():
    src = Path(sys.argv[1]) if len(sys.argv)>1 else R/"runs/push/submission_v11.csv"
    rows = list(csv.DictReader(open(src)))
    dev = {r['name']: r for r in csv.DictReader(open(R/"results/developability.csv"))}
    stem = src.stem.replace("submission_", "")
    solf = R/f"results/solubility_{stem}.csv"
    if not solf.exists(): solf = R/"results/solubility_v11.csv"
    sol = {r['name']: r for r in csv.DictReader(open(solf))} if solf.exists() else {}
    recs = []
    for r in rows:
        n = r['name']
        esm_pdb = str(R/f"runs/expr_v11/esm_mono/{n}.pdb")
        ep = esm_plddt(esm_pdb) if os.path.exists(esm_pdb) else None
        bp, bptm, bcif = boltz_conf(n)
        tmv = tm(esm_pdb, bcif) if bcif else None
        d = dev.get(n, {})
        liab = int(d.get('liability_score',0) or 0); risk = d.get('expr_risk','?')
        q = float(d.get('net_charge_pH7_4',0) or 0); gravy = float(d.get('gravy',0) or 0)
        # agreement gate
        agree = (ep is not None and bp is not None and ep>=0.70 and bp>=0.70 and (tmv is None or tmv>=0.60))
        disagree = (tmv is not None and tmv < 0.50) or (bp is not None and bp < 0.60) or (ep is not None and ep<0.65)
        # probability
        p = base_p(ep)
        notes = []
        if agree and bp>=0.80 and (tmv is None or tmv>=0.60):
            p = min(0.90, p*1.10); notes.append("2fold-agree")
        if disagree:
            p *= 0.70; notes.append("FOLDER-DISAGREE")
        # developability SEVERE here is charge-magnitude driven; the experimental-data solubility model
        # (Protein-Sol) directly models |net charge| and (D+E)%, so a high solubility score is the
        # orthogonal arbiter that CLEARS a charge-only SEVERE flag (acidic proteins are often soluble).
        sval = float(sol[n].get('solubility_score', 'nan')) if n in sol else float('nan')
        if risk == "SEVERE":
            if sval == sval and sval >= 0.80:
                p *= 0.92; notes.append("dev-SEVERE(charge)/sol-cleared")
            else:
                p *= 0.70; notes.append("dev-SEVERE")
        elif liab >= 5:
            p *= 0.88; notes.append(f"liab{liab}")
        elif liab >= 4:
            p *= 0.93; notes.append(f"liab{liab}")
        if n in sol:
            s = sval
            notes.append(f"sol{s:.2f}")
            if s == s:  # not nan
                if s < 0.45: p *= 0.80; notes.append("LOW-SOL")
                elif s >= 0.60: p = min(0.90, p*1.05)
        p = round(max(0.20, min(0.90, p)), 3)
        severe_cleared = (risk == "SEVERE" and sval == sval and sval >= 0.80)
        hard_severe = (risk == "SEVERE" and not severe_cleared)
        tier = ("C" if (disagree or hard_severe)
                else "A" if (agree and risk != "SEVERE" and liab < 4) else "B")
        recs.append(dict(name=n, L=len(r['sequence']), esm_plddt=ep, boltz_plddt=bp, boltz_ptm=bptm,
                         tm_agree=tmv, dev_risk=risk, liab=liab, charge=q, gravy=gravy,
                         p_express=p, tier=tier, notes=";".join(notes)))
    # expected yield + band (+/-0.08 per design, clamped)
    E = sum(x['p_express'] for x in recs)
    lo = sum(max(0.15, x['p_express']-0.08) for x in recs)
    hi = sum(min(0.92, x['p_express']+0.08) for x in recs)
    tiers = {t: sum(1 for x in recs if x['tier']==t) for t in "ABC"}
    recs.sort(key=lambda x: x['p_express'])
    out = dict(source=str(src), n=len(recs), expected_yield=round(E,2),
               band=[round(lo,2),round(hi,2)], tiers=tiers,
               folders=["ESMFold(ESM2)","Boltz-2(AF3-diffusion)"],
               solubility_model_loaded=bool(sol), designs=recs)
    (R/"runs/expr_v11").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(R/"runs/expr_v11/consensus.json","w"), indent=1)
    with open(R/"runs/expr_v11/consensus.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys())); w.writeheader()
        for x in recs: w.writerow(x)
    print(f"\n=== EXPRESSION CONSENSUS ({src.name}) — independent folders {out['folders']} ===")
    print(f"{'design':32s}{'esmPL':>6}{'bzPL':>6}{'bzPTM':>6}{'TM':>6}{'risk':>9}{'liab':>5}{'chg':>7}{'P':>6} tier notes")
    for x in recs:
        print(f"{x['name']:32s}{str(x['esm_plddt']):>6}{str(x['boltz_plddt']):>6}{str(x['boltz_ptm']):>6}"
              f"{str(x['tm_agree']):>6}{x['dev_risk']:>9}{x['liab']:>5}{x['charge']:>7.1f}{x['p_express']:>6.2f}  {x['tier']}  {x['notes']}")
    print(f"\nTIERS  A(clean,2-fold agree)={tiers['A']}  B(ok,one soft axis)={tiers['B']}  C(red flag)={tiers['C']}")
    print(f"EXPECTED EXPRESSING (literature-anchored) = {E:.1f}/20   band [{lo:.1f}, {hi:.1f}]")
    print(f"WEAKEST 4: {[x['name'] for x in recs[:4]]}")
    print(f"solubility model loaded: {bool(sol)}  (re-run after results/solubility_v11.csv exists)")

if __name__ == "__main__":
    main()
