#!/usr/bin/env python3
"""Control-gated verdict for the two-state MD release gate, scored on the His-Lys distance LADDER.

For a real pH-release switch the His-Lys min distance should order:
    HIP (His+, pH6)  >  HIE (His0, pH7.4)  >=  GLN (neutral isostere dud)
i.e. protonation PUSHES the His away from the conserved Lys (Coulomb repulsion), while the neutral
isostere does not (often it even H-bonds closer).  Signals:
    d_switch = his_lys(HIP) - his_lys(HIE)      pH response (>0 = repulsion at low pH)   [PRIMARY]
    d_iso    = his_lys(HIP) - his_lys(GLN)      causal: charge-on vs neutral-isostere    [CAUSAL]
    noise    = max SEM(HIP,HIE) over replicas   test-retest floor
Verdict:
    CONFIRMED  : d_switch > 2*noise AND d_switch > 0 AND d_iso > 0
    NONRESPOND : d_switch <= 2*noise (flat ladder -> MD cannot confirm a switch)
Secondary (slow, corroborative only): d_local = local_iface(HIE) - local_iface(HIP) (>0 = local release).
MD is a CONFIRMER/flagger, not a sole cutter: a NONRESPOND design is not auto-dropped, but is
deprioritised vs CONFIRMED designs of equal static/novelty quality.

Usage: md_verdict.py <dir-or-glob of *.json>
"""
import sys, glob, json, os

def g(agg, v, m):
    return agg.get(v, {}).get(m) if agg else None

def main():
    pat = sys.argv[1] if len(sys.argv) > 1 else "runs/md/gate"
    files = sorted(glob.glob(os.path.join(pat, "*.json")) if os.path.isdir(pat) else glob.glob(pat))
    rows = []
    for f in files:
        d = json.load(open(f)); agg = d.get("agg", {})
        hip, hie, gln = g(agg, "HIP", "his_lys_mindist"), g(agg, "HIE", "his_lys_mindist"), g(agg, "GLN", "his_lys_mindist")
        sem = max([x for x in (g(agg, "HIP", "sem_his_lys"), g(agg, "HIE", "sem_his_lys")) if x is not None] or [0.0])
        d0 = next((s.get("lys_init_dist") for s in d.get("states", []) if "error" not in s), None)
        d_switch = (hip - hie) if None not in (hip, hie) else None
        d_iso = (hip - gln) if None not in (hip, gln) else None
        d_local = d.get("d_local_release")
        verdict = "NA"
        if d_switch is not None:
            if d_switch > 2 * sem and d_switch > 0 and (d_iso is None or d_iso > 0):
                verdict = "CONFIRMED"
            elif d_switch <= 2 * sem:
                verdict = "NONRESPOND"
            else:
                verdict = "MARGINAL"
        rows.append(dict(name=d["name"], d0=d0, hip=hip, hie=hie, gln=gln,
                         d_switch=d_switch, d_iso=d_iso, sem=sem, d_local=d_local, verdict=verdict))
    rows.sort(key=lambda r: (r["verdict"] != "CONFIRMED", -(r["d_switch"] or -9)))
    print(f"{'design':50} {'d0':>5} {'HIP':>5} {'HIE':>5} {'GLN':>5} {'dSwitch':>8} {'dIso':>6} {'2sem':>5} {'dLocal':>6}  verdict")
    for r in rows:
        def fm(x, w=5, p=2): return (f"{x:>{w}.{p}f}" if isinstance(x, (int, float)) else f"{'--':>{w}}")
        print(f"{r['name'][:50]:50} {fm(r['d0'])} {fm(r['hip'])} {fm(r['hie'])} {fm(r['gln'])} "
              f"{fm(r['d_switch'],8)} {fm(r['d_iso'],6)} {fm(2*r['sem'],5)} {fm(r['d_local'],6)}  {r['verdict']}")
    nconf = sum(1 for r in rows if r["verdict"] == "CONFIRMED")
    nnr = sum(1 for r in rows if r["verdict"] == "NONRESPOND")
    print(f"\nCONFIRMED {nconf} | MARGINAL {sum(1 for r in rows if r['verdict']=='MARGINAL')} | NONRESPOND {nnr} | total {len(rows)}")

if __name__ == "__main__":
    main()
