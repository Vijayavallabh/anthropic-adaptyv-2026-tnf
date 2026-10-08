#!/usr/bin/env python3
"""WAVE 16j / deep-research M1: His433 in-window proton-linkage MAGNITUDE + ensemble uncertainty.

The quad gate confirms the SIGN of the pH switch; it does not prove the switch has real
MAGNITUDE inside the tumour-vs-blood window [6.5, 7.4]. A buried, salt-bridged His433 can sit
where PROPKA carries a ~2 pH-unit one-sided error, so a positive sign can coexist with a
negligible in-window proton linkage (the "sign-right, magnitude-dead" failure the deep-research
report flags from four converging lenses).

This re-uses the project's own Wyman definition (obj1_switch.ddg_ph_from_pdb, His433 term):
    S_His433 = RT*ln10 * integral_{6.5}^{7.4} [ theta(pH; pKa_complex) - theta(pH; pKa_free) ] dpH
    theta(pH,pKa) = 1/(1+10^(pH-pKa));  RT*ln10 = 1.419 kcal/mol @310K
but instead of a single structure it propagates the 20-frame MD ENSEMBLE of the COMPLEX His433
pKa (runs/ens_propka/egfr_phsw_NN.json) -> S_His433 mean +/- sd, sign-stability across frames,
and a +/-1.5 pKa worst-case sensitivity. pKa_free = the UNBOUND Domain-III His433 pKa (one PROPKA
run on the reference crop); the binder is designed to RAISE His433's pKa on binding, so
pKa_complex > pKa_free and S_His433 > 0 = tumour-ON.

Pure stdlib (+ propka3 for the single free-reference run). Writes results/wave16j_switch_physicality.csv.
Does NOT modify the frozen submission; this is a validation annotation only.
"""
import csv, json, glob, os, subprocess, tempfile, shutil, statistics, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RT_LN10 = 1.419
GRID = [6.5 + 0.05 * k for k in range(19)]           # 6.5 .. 7.4, matches obj1_switch
HIS_CLAMP_HI = 8.5                                   # obj1_switch PKA_CLAMP HIS upper (desolvation artifact guard)

def theta(pH, pka):
    return 1.0 / (1.0 + 10 ** (pH - pka))

def trapz(vals):
    return sum((vals[i] + vals[i+1]) / 2 * (GRID[i+1] - GRID[i]) for i in range(len(vals) - 1))

def S_link(pka_cx, pka_free):
    """His433 Wyman in-window linkage (kcal/mol); >0 = binds tighter at pH 6.5 (tumour-ON)."""
    return RT_LN10 * trapz([theta(p, pka_cx) - theta(p, pka_free) for p in GRID])

def free_ref_pKa():
    """Unbound Domain-III His433 pKa via one PROPKA run on the reference crop (author resseq 433)."""
    mp = json.load(open(REPO / "data/targets/6ARU_domainIII_map.json"))
    his_resseq = int(mp.get("ph_sensor_his433", 433)) if str(mp.get("ph_sensor_his433", "433")).isdigit() else 433
    pk = (sorted(glob.glob(str(REPO / "envs/*/bin/propka3"))) or [None])[0]
    if not pk:
        sys.exit("propka3 not found under envs/*/bin")
    wd = tempfile.mkdtemp(prefix="m1free_")
    try:
        shutil.copyfile(REPO / "data/targets/6ARU_domainIII.pdb", os.path.join(wd, "diii.pdb"))
        subprocess.run([pk, "diii.pdb"], cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pka = os.path.join(wd, "diii.pka")
        his = {}
        inss = False
        for line in open(pka):
            if "SUMMARY OF THIS PREDICTION" in line:
                inss = True; continue
            if inss:
                t = line.split()
                if len(t) >= 4 and t[0] == "HIS":
                    try: his[int(t[1])] = float(t[3])
                    except ValueError: pass
                elif line.strip().startswith("-") and his:
                    break
        if his_resseq in his:
            return his[his_resseq], his_resseq, his
        # fallback: the His whose pKa is closest to a titratable value if exact resseq absent
        return (his.get(his_resseq) or statistics.median(list(his.values()))), his_resseq, his
    finally:
        shutil.rmtree(wd, ignore_errors=True)

def main():
    free, his_resseq, allhis = free_ref_pKa()
    print(f"free-target (unbound DIII) His433 pKa @resseq {his_resseq} = {free:.2f}  (all His in crop: {allhis})")
    free_titr = theta(6.5, free) - theta(7.4, free)   # how much the FREE His433 titrates in-window (max 0.476 @pKa6.95)
    rows = []
    for jf in sorted(glob.glob(str(REPO / "runs/ens_propka/egfr_phsw_*.json"))):
        d = json.load(open(jf))
        name = d["name"]
        frames = [min(p, HIS_CLAMP_HI) for p in d["his433_pka_frames"]]   # clamp complex pKa as the pipeline does
        raw = d["his433_pka_frames"]
        S = [S_link(p, free) for p in frames]
        m, sd = statistics.mean(S), (statistics.pstdev(S) if len(S) > 1 else 0.0)
        fpos = sum(1 for s in S if s > 0) / len(S)
        # +/-1.5 worst-case: push complex DOWN and free UP (both toward killing the switch)
        S_worst = statistics.mean([S_link(max(p - 1.5, 3.0), min(free + 1.5, 8.5)) for p in frames])
        cx_mean = statistics.mean(frames)
        cx_self = theta(6.5, cx_mean) - theta(7.4, cx_mean)   # complex self-titration in-window (want ~0 = saturated)
        # verdict: real in-window switch needs a sign-stable, non-negligible magnitude
        dead = (m < 0.15) or (fpos < 0.9) or (m - sd <= 0)
        rows.append(dict(name=name, S_his433_mean=round(m, 3), S_his433_sd=round(sd, 3),
                         S_minus_sd=round(m - sd, 3), frac_frames_pos=round(fpos, 2),
                         S_worst_pm1p5=round(S_worst, 3), cx_pka_mean=round(cx_mean, 2),
                         cx_pka_raw_mean=round(statistics.mean(raw), 2), cx_self_titr_inwin=round(cx_self, 3),
                         verdict=("WEAK/dead-in-window" if dead else "robust-in-window")))
    out = REPO / "results/wave16j_switch_physicality.csv"
    cols = list(rows[0].keys())
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow(r)
    print(f"\nfree His433 in-window titration theta(6.5)-theta(7.4) = {free_titr:.3f}  (max 0.476 @pKa 6.95)\n")
    print(f"{'design':15s}{'S_His433':>9}{'±sd':>7}{'m-sd':>7}{'frac>0':>7}{'S±1.5':>7}{'cx_pKa':>7}{'verdict':>22}")
    for r in sorted(rows, key=lambda x: x["S_his433_mean"]):
        print(f"{r['name']:15s}{r['S_his433_mean']:>9}{r['S_his433_sd']:>7}{r['S_minus_sd']:>7}"
              f"{r['frac_frames_pos']:>7}{r['S_worst_pm1p5']:>7}{r['cx_pka_mean']:>7}{r['verdict']:>22}")
    Sm = [r["S_his433_mean"] for r in rows]
    weak = [r["name"] for r in rows if r["verdict"].startswith("WEAK")]
    print(f"\nN={len(rows)} | S_His433 mean: min={min(Sm):.3f} median={statistics.median(Sm):.3f} max={max(Sm):.3f}")
    print(f"sign-stable robust-in-window: {len(rows)-len(weak)}/{len(rows)} | WEAK/dead: {weak or 'none'}")
    print(f"wrote {out}")

if __name__ == "__main__":
    main()
