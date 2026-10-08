#!/usr/bin/env python3
"""Developability / wet-lab-liability scan for the final submission sequences.

Flags the sequence liabilities that most depress experimental expression/solubility
and bias de-novo-binder wet-lab outcomes (the in-silico -> wet-lab gap):
  - free cysteines (odd count => unpaired thiol; any free Cys => scrambling/aggregation)
  - N-glycosylation sequons  N-X-[S/T], X != P
  - deamidation hotspots     N[G/S/T]  (Asn, fast)
  - isomerization hotspots   D[G/S/T/H]  (Asp)
  - Met (oxidation) / Trp (oxidation) counts
  - net charge @ pH 7.4, isoelectric-ish sign, GRAVY (Kyte-Doolittle)
  - fraction hydrophobic (AILMFWV) + longest hydrophobic run

Pure stdlib. Reads a submission CSV (name,sequence,...), writes results/developability.csv
and prints a flagged summary. Does NOT modify the submission.

Usage: python scripts/developability_scan.py [results/submission_track3_v2.csv]
"""
import csv, sys, re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
KD = {  # Kyte-Doolittle hydropathy
    'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
    'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
    'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
HYDRO = set("AILMFWV")
# simple side-chain pKa (Nterm/Cterm added once)
PKA = {'D':3.65,'E':4.25,'C':8.3,'Y':10.07,'H':6.0,'K':10.53,'R':12.48}

def net_charge(seq, pH=7.4):
    def frac_pos(pka): return 1.0/(1.0+10**(pH-pka))
    def frac_neg(pka): return 1.0/(1.0+10**(pka-pH))
    q = frac_pos(8.0) - frac_neg(3.1)  # termini
    for aa in seq:
        if aa in ('K','R'): q += frac_pos(PKA[aa])
        elif aa == 'H': q += frac_pos(PKA['H'])
        elif aa in ('D','E'): q -= frac_neg(PKA[aa])
        elif aa in ('C','Y'): q -= frac_neg(PKA[aa])
    return round(q, 2)

ACIDIC = set("DE")
BASIC = set("KR")

def longest_run(seq, charset):
    best = cur = 0
    for aa in seq:
        cur = cur + 1 if aa in charset else 0
        best = max(best, cur)
    return best

def longest_hydro_run(seq):
    return longest_run(seq, HYDRO)

def aa_entropy(seq):
    """Shannon entropy (bits) of the amino-acid composition; low entropy = low-complexity sequence."""
    import math
    from collections import Counter
    n = len(seq)
    if not n:
        return 0.0
    return round(-sum((c / n) * math.log2(c / n) for c in Counter(seq).values()), 2)


# Manning counterion condensation: counterions condense onto a polyelectrolyte when the Manning
# parameter xi = l_B / b exceeds 1, where l_B ~ 7.1 A is the Bjerrum length in water at ~body temperature
# and b is the axial spacing between like charges. In a run of consecutive same-sign charged residues
# every position carries a like charge, so b ~ the per-residue rise (~1.5 A on a helix) and xi ~ 4.7 >> 1:
# such runs are short polyelectrolyte segments above the condensation threshold, which lowers the
# effective charge and degrades solubility / cell-free expression. This is the physical basis for the
# charge-run and charged-fraction liabilities; polyelectrolyte_frac measures how much of the chain sits
# in such same-sign condensation runs. (A 1D rod approximation for a folded protein, used qualitatively.)
L_BJERRUM = 7.1   # Angstrom, Bjerrum length in water near physiological temperature
HELIX_RISE = 1.5  # Angstrom per residue, axial rise (these designs are helical bundles)

def polyelectrolyte_frac(seq):
    """Fraction of residues inside a maximal same-sign charged run of length >= 3 (Manning xi > 1 regime)."""
    def sign(a):
        return "ac" if a in ACIDIC else ("ba" if a in BASIC else None)
    n, j, total = len(seq), 0, 0
    while j < n:
        if sign(seq[j]) is None:
            j += 1; continue
        k = j
        while k < n and sign(seq[k]) == sign(seq[j]):
            k += 1
        if k - j >= 3:
            total += k - j
        j = k
    return round(total / max(n, 1), 3)

def scan(name, seq):
    seq = (seq or "").strip().upper()
    n = len(seq)
    cys = seq.count('C')   # NB: TOTAL cysteines, not free thiols -- without a 3D structure we cannot know
                           # disulfide pairing, so this is a CONSERVATIVE proxy (every Cys treated as a risk).
                           # The manifest column is named free_cys for historical reasons; it is total Cys.
    glyc = len(re.findall(r'(?=(N[^P][ST]))', seq))   # lookahead => count overlapping sequons too
    deamid = len(re.findall(r'(?=(N[GST]))', seq))
    isom = len(re.findall(r'(?=(D[GSTH]))', seq))
    gravy = round(sum(KD.get(a,0) for a in seq)/max(n,1), 3)      # max(n,1): never divide by zero on an empty cell
    hyd_frac = round(sum(1 for a in seq if a in HYDRO)/max(n,1), 3)
    q = net_charge(seq)
    # charge-blob / low-complexity liabilities: these de-novo binders are charge-rich, and long
    # same-sign charge runs plus low sequence complexity are the dominant cell-free expression and
    # solubility risks here (more than the hydrophobic runs the standard panel flags).
    charged_frac = round(sum(1 for a in seq if a in "DEKR")/max(n,1), 3)
    de_run = longest_run(seq, ACIDIC); kr_run = longest_run(seq, BASIC)
    max_charge_run = max(de_run, kr_run)
    ent = aa_entropy(seq)
    poly = polyelectrolyte_frac(seq)   # fraction in same-sign Manning-condensation runs (xi>1)
    if (charged_frac > 0.55 and max_charge_run >= 5) or ent < 3.0 or abs(q) > 15 or poly > 0.35:
        expr_risk = "SEVERE"
    elif charged_frac > 0.50 or max_charge_run >= 4 or poly > 0.20:
        expr_risk = "MODERATE"
    else:
        expr_risk = "CLEAN"
    rec = dict(name=name, length=n, free_cys=cys, cys_odd=int(cys % 2 == 1),
               nglyc=glyc, deamid_NG=deamid, isom_DG=isom,
               met=seq.count('M'), trp=seq.count('W'),
               net_charge_pH7_4=q, gravy=gravy,
               hydrophobic_frac=hyd_frac, max_hydro_run=longest_hydro_run(seq),
               charged_frac=charged_frac, max_charge_run=max_charge_run,
               aa_entropy=ent, poly_frac=poly, expr_risk=expr_risk)
    flags = []
    if cys >= 2: flags.append(f"Cys{cys}")
    elif cys == 1: flags.append("Cys1(unpaired)")
    if glyc: flags.append(f"Nglyc{glyc}")
    if gravy > 0: flags.append(f"GRAVY+{gravy}")
    if hyd_frac > 0.45: flags.append(f"hyd{hyd_frac}")
    if abs(q) > 12: flags.append(f"charge{q}")
    if rec['max_hydro_run'] >= 5: flags.append(f"hrun{rec['max_hydro_run']}")
    if max_charge_run >= 5: flags.append(f"chgrun{max_charge_run}")
    if ent < 3.0: flags.append(f"lowcplx{ent}")
    if charged_frac > 0.55: flags.append(f"chg{charged_frac}")
    if poly > 0.20: flags.append(f"poly{poly}")
    rec['flags'] = ";".join(flags)
    rec['liability_score'] = (2*cys + 2*glyc + (2 if gravy>0 else 0)
                              + (1 if hyd_frac>0.45 else 0) + rec['max_hydro_run']//5
                              + {"SEVERE": 2, "MODERATE": 1, "CLEAN": 0}[expr_risk])
    return rec

def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO/"results/submission_track3_v2.csv"
    rows = list(csv.DictReader(open(src)))
    recs = [scan(r['name'], r['sequence']) for r in rows]
    cols = ['name','length','free_cys','cys_odd','nglyc','deamid_NG','isom_DG','met','trp',
            'net_charge_pH7_4','gravy','hydrophobic_frac','max_hydro_run',
            'charged_frac','max_charge_run','aa_entropy','poly_frac','expr_risk','liability_score','flags']
    out = REPO/"results"/("developability.csv" if "submission" in src.name else f"developability_{src.stem}.csv")
    with open(out,'w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in recs: w.writerow(r)
    print(f"wrote {out}\n")
    print(f"{'name':<13}{'len':>4}{'Cys':>4}{'glyc':>5}{'chg%':>6}{'crun':>5}{'ent':>6}{'poly%':>6}{'risk':>9}{'LIAB':>5}  flags")
    for r in sorted(recs, key=lambda x:-x['liability_score']):
        print(f"{r['name']:<13}{r['length']:>4}{r['free_cys']:>4}{r['nglyc']:>5}"
              f"{r['charged_frac']*100:>5.0f}%{r['max_charge_run']:>5}{r['aa_entropy']:>6}{r['poly_frac']*100:>5.0f}%{r['expr_risk']:>9}"
              f"{r['liability_score']:>5}  {r['flags']}")
    flagged = [r['name'] for r in recs if r['liability_score'] >= 4]
    severe = [r['name'] for r in recs if r['expr_risk'] == "SEVERE"]
    print(f"\nHIGH-LIABILITY (score>=4): {flagged or 'none'}")
    print(f"EXPRESSION-RISK SEVERE (charge-blob / low-complexity): {severe or 'none'}")

if __name__ == "__main__":
    main()
