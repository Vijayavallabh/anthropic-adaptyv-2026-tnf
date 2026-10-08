#!/usr/bin/env python3
"""WAVE 16L / deep-research M2: INDEPENDENT burial-aware His433 pKa via JustHISpKa (NRC, JCIM 2025).

An independent 5th pKa leg, orthogonal to PROPKA/pKAI (different engine: Amber-FF + Solvated-
Interaction-Energy + rotamer sampling + a His-specific buried correction; MUE 0.40 on buried His).
Resolves M1's ~1.5-unit uncertainty caveat on the engineered, buried, salt-bridged His433.

Per design: prep the complex (pdb4amber -> tleap ff14SB -> ambpdb mol2 -> push_renum, with the
5->4 col renum fix for this pdb4amber version), then run justHisPKa on His433 (= target chain-B
crop-ordinal 100) with the binder Asp/Glu salt-bridge partner FORCED to its salt-bridge
(deprotonated) state to resolve the His<->Asp/Glu coupling. Bound pKa vs the unbound-DIII
reference (His433 = author His409 = SID 100; JustHISpKa pKa 6.088) = an independent Wyman shift.
A positive shift (binding raises His433 pKa) = tumour-ON, same sign as the quad legs.

CPU only (AmberTools + justHisPKa 6 threads). Validation; frozen set unchanged.
Usage: PATH=.../ambertools/bin:$PATH python scripts/wave16L_justhispka.py [--only egfr_phsw_18]
"""
import os, sys, csv, json, glob, subprocess, tempfile, shutil, math
from collections import defaultdict

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
AMBER = "/home/prachh/miniconda3/envs/ambertools"
JHP = os.path.join(REPO, "third_party/justhispka/justhispka-1.0.2")
EXE = os.path.join(JHP, "justHisPKa"); LIC = os.path.join(JHP, "justHisPKa.lic")
PUSH = os.path.join(JHP, "prep.dir/bin/push_renum.sh")
FREE_PKA = 6.088   # JustHISpKa unbound-DIII His433 (author His409 = SID 100)
ENV = dict(os.environ, AMBERHOME=AMBER, PATH=f"{AMBER}/bin:" + os.environ.get("PATH", ""))
HIS_ORD = 100      # His433 crop ordinal in the complex

def sh(cmd, cwd, **kw):
    return subprocess.run(cmd, cwd=cwd, env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, **kw)

def heavy(pdb):
    ch = defaultdict(lambda: defaultdict(list))
    for l in open(pdb):
        if l[:6] in ("ATOM  ", "HETATM"):
            an = l[12:16].strip(); rn = l[17:20].strip(); c = l[21]
            try:
                num = int(l[22:26]); xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                continue
            ch[c][(rn, num)].append((an, xyz))
    return ch

def partner(pdb):
    """binder (2nd-largest chain) ASP/GLU with sidechain-O nearest His(ord100) in target -> (resn,num,chain,dist)."""
    ch = heavy(pdb)
    order = sorted(ch, key=lambda c: -sum(len(v) for v in ch[c].values()))
    tc, bc = order[0], order[1]
    # His433 = target residue at ordinal 100
    tkeys = list(ch[tc].keys())
    his = tkeys[HIS_ORD - 1] if len(tkeys) >= HIS_ORD and tkeys[HIS_ORD - 1][0] == "HIS" else \
          next((k for k in tkeys if k[0] == "HIS" and k[1] == HIS_ORD), None)
    if his is None:
        return None
    hn = [xyz for an, xyz in ch[tc][his] if an in ("ND1", "NE2")]
    best = (None, 1e9)
    for k, ats in ch[bc].items():
        if k[0] in ("ASP", "GLU"):
            for an, o in ats:
                if an in ("OD1", "OD2", "OE1", "OE2"):
                    for n in hn:
                        d = math.dist(o, n)
                        if d < best[1]:
                            best = ((k[0], k[1], bc), d)
    return (best[0][0], best[0][1], best[0][2], round(best[1], 2)) if best[0] else None

def prep(pdb, wd):
    shutil.copyfile(pdb, os.path.join(wd, "in.pdb"))
    sh(["pdb4amber", "-i", "in.pdb", "-o", "p2.pdb", "--dry", "--nohyd"], wd)
    # 5-col ("name chain orignum newname newseq") -> 4-col ("name orignum chain newseq")
    rn = os.path.join(wd, "p2_renum.txt")
    if not os.path.exists(rn):
        return None, None, "no renum"
    with open(rn) as fh, open(os.path.join(wd, "p2fix_renum.txt"), "w") as out:
        for line in fh:
            t = line.split()
            if len(t) >= 5:
                out.write(f"{t[0]} {t[2]} {t[1]} {t[4]}\n")
    tleap = "source leaprc.protein.ff14SB\nm = loadpdb p2.pdb\nsaveamberparm m p3.prmtop p3.inpcrd\nquit\n"
    open(os.path.join(wd, "leap.in"), "w").write(tleap)
    r = sh(["tleap", "-f", "leap.in"], wd)
    if not os.path.exists(os.path.join(wd, "p3.prmtop")):
        return None, None, "tleap:" + r.stdout[-200:]
    r = sh(["bash", "-c", "ambpdb -p p3.prmtop -c p3.inpcrd -mol2 -sybyl > p4.mol2"], wd)
    r = sh(["bash", PUSH, "p2.pdb", "p2fix_renum.txt", "p4.mol2"], wd)
    open(os.path.join(wd, "prep.mol2"), "w").write(r.stdout if r.stdout and r.stdout.startswith("@") else "")
    # push_renum prints mol2 to stdout; capture properly
    with open(os.path.join(wd, "prep.mol2"), "w") as out:
        rr = subprocess.run(["bash", PUSH, "p2.pdb", "p2fix_renum.txt", "p4.mol2"], cwd=wd, env=ENV,
                            stdout=out, stderr=subprocess.PIPE, text=True)
    if os.path.getsize(os.path.join(wd, "prep.mol2")) == 0:
        return None, None, "push_renum:" + rr.stderr[-150:]
    return os.path.join(wd, "prep.mol2"), os.path.join(wd, "p3.prmtop"), ""

def sub_sid_map(mol2):
    """SUBSTRUCTURE records -> {(chain, resname3, orignum): sid}."""
    m = {}; ins = False
    for l in open(mol2):
        if l.startswith("@<TRIPOS>SUBSTRUCTURE"):
            ins = True; continue
        if ins:
            if l.startswith("@"):
                break
            t = l.split()
            if len(t) >= 6:
                m[(t[5], t[1])] = t[0]   # (chain, resname+orignum e.g. "ASP81") -> SID
    return m

def run_jhp(mol2, prmtop, his_res, force_opts):
    cmd = [EXE, "-b", f"MOL={mol2}", f"PRM={prmtop}", f"RES={his_res}", f"LICENSE={LIC}"] + force_opts
    r = subprocess.run(cmd, env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = r.stdout.strip()
    for line in reversed(out.splitlines()):      # -b brief mode = bare float on the last line
        s = line.strip()
        try:
            return float(s), out
        except ValueError:
            pass
        for tok in s.split():
            if tok.startswith("PKA="):
                try:
                    return float(tok.split("=")[1]), out
                except ValueError:
                    pass
    return None, out

def main():
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]
    rows = []
    jsons = sorted(glob.glob(os.path.join(REPO, "runs/ens_propka/egfr_phsw_*.json")))
    outdir = os.path.join(REPO, "runs/wave16L_jhp"); os.makedirs(outdir, exist_ok=True)
    for jf in jsons:
        d = json.load(open(jf)); name = d["name"]; pdb = d["pdb"]
        if only and name != only:
            continue
        rec = {"name": name, "free_pka": FREE_PKA, "err": ""}
        if not os.path.exists(pdb):
            rec["err"] = "no pdb"; rows.append(rec); print(name, "no pdb", flush=True); continue
        p = partner(pdb)
        if not p:
            rec["err"] = "no partner"; rows.append(rec); print(name, "no partner", flush=True); continue
        prn, pnum, pch, pdist = p
        rec.update(partner=f"{prn}{pnum}{pch}", partner_dist=pdist)
        wd = tempfile.mkdtemp(prefix=f"jhp_{name}_", dir=outdir)
        try:
            mol2, prmtop, err = prep(pdb, wd)
            if not mol2:
                rec["err"] = "prep:" + err; rows.append(rec); print(name, "PREP FAIL", err, flush=True); continue
            smap = sub_sid_map(mol2)
            psid = smap.get((pch, f"{prn}{pnum}"))
            his_res = f"B:HIS{HIS_ORD}"
            force_key = {"ASP": "FORCE_ASP", "GLU": "FORCE_GLU"}[prn]
            # no-force first; if His couples ("multiple titratable groups"), force the partner deprotonated (salt-bridge state)
            pka, out = run_jhp(mol2, prmtop, his_res, [])
            if pka is None and "multiple" in out.lower() and psid:
                pka, out = run_jhp(mol2, prmtop, his_res, [f"{force_key}={psid}"])
                rec["forced"] = f"{force_key}={psid}"
            rec.update(partner_sid=psid, bound_pka=pka, jhp_out=out[:120])
            if pka is not None:
                rec["wyman_shift"] = round(pka - FREE_PKA, 2)
                rec["tumour_on"] = pka > FREE_PKA
            print(f"{name}: partner {rec['partner']} sid={psid} | bound His433 pKa={pka} "
                  f"| shift vs free {rec.get('wyman_shift')} | {out[:80]}", flush=True)
        except Exception as e:
            rec["err"] = "exc:" + repr(e)[:80]; print(name, "EXC", rec["err"], flush=True)
        finally:
            rows.append(rec)
    cols = ["name", "partner", "partner_dist", "partner_sid", "free_pka", "bound_pka",
            "wyman_shift", "tumour_on", "jhp_out", "err"]
    out = os.path.join(REPO, "results/wave16L_justhispka.csv")
    if not only:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
            for r in rows: w.writerow({k: r.get(k, "") for k in cols})
        ok = [r for r in rows if r.get("bound_pka") is not None]
        print(f"\nDONE {len(ok)}/{len(rows)} -> {out}")
        if ok:
            pos = sum(1 for r in ok if r.get("tumour_on"))
            sh_ = [r["wyman_shift"] for r in ok]
            print(f"tumour-ON (bound pKa > free {FREE_PKA}): {pos}/{len(ok)} | Wyman shift min/med/max: "
                  f"{min(sh_)}/{sorted(sh_)[len(sh_)//2]}/{max(sh_)}")

if __name__ == "__main__":
    main()
