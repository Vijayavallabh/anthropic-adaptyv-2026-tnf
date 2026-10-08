#!/usr/bin/env python3
"""Ensemble-PROPKA: conformation-averaged His433 pKa — the research-endorsed, in-env substitute for
constant-pH MD (which has no pip/uv GPU engine here). For each complex we run short Cα-restrained
implicit-solvent MD with His433 built as HIP (the pH-6.5 bound state), dump N snapshots, run PROPKA3
on each, and report the His433 pKa distribution (mean ± sd). A robust tumor-ON switch shows the
ensemble His433 pKa elevated in the complex with a tight spread — removing PROPKA's single-structure
blind spot without the CpHMD install risk. Reuses scripts/md_pH_switch.build (Amber14+GBn2, HMR 4fs).

Usage: envs/md/bin/python scripts/ensemble_propka.py --pdb-list L.tsv --out runs/ens_propka \
         [--device 0] [--nframes 20] [--equil-ps 100] [--prod-ps 1000] [--restraint-k 2.0]
  L.tsv: a header row then <name>\\t...\\t<pdb> (final20_pdbs.tsv works), or plain "<name>\\t<pdb>".
Needs propka3 on PATH (export PATH=$PWD/envs/pyrosetta/bin:$PATH). Writes <out>/<name>.json.
"""
import sys, os, json, subprocess, tempfile, argparse, io, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from openmm import unit, LangevinMiddleIntegrator, Platform
from openmm import app
from md_pH_switch import build, add_ca_restraints

PROPKA = os.environ.get("PROPKA3", "propka3")
HISVAR = {"HIP", "HIE", "HID", "HSP", "HSE", "HSD"}

def write_frame(topology, positions, path):
    """heavy-atom PDB with His variants renamed HIS (PROPKA needs standard residue names)."""
    buf = io.StringIO()
    app.PDBFile.writeFile(topology, positions, buf)
    with open(path, "w") as fh:
        for l in buf.getvalue().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                if l[76:78].strip() == "H":
                    continue
                rn = l[17:20].strip()
                if rn in HISVAR:
                    l = l[:17] + "HIS" + l[20:]
                fh.write(l + "\n")
        fh.write("END\n")

def propka_his433(frame_pdb, workdir):
    """run propka3 in workdir; parse His 100 chain B pKa from the .pka SUMMARY block."""
    try:
        subprocess.run([PROPKA, os.path.abspath(frame_pdb)], cwd=workdir,
                       check=True, capture_output=True, timeout=300)
    except Exception:
        return None
    pka = Path(workdir) / f"{Path(frame_pdb).stem}.pka"
    if not pka.exists():
        return None
    insum = False
    for line in open(pka):
        if "SUMMARY OF THIS PREDICTION" in line:
            insum = True; continue
        if insum:
            if line.strip().startswith("-") or "Free energy" in line or "Protein charge" in line:
                break
            p = line.split()
            if len(p) >= 4 and p[0] == "HIS" and p[1] == "100" and p[2] == "B":
                try: return float(p[3])
                except ValueError: pass
    return None

def run_design(name, pdb, a):
    ff, modeller, system, target, binder, his = build(pdb, "HIP")
    add_ca_restraints(system, modeller, a.restraint_k)
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(777)
    sim = app.Simulation(modeller.topology, system, integ,
                         Platform.getPlatformByName("OpenCL"),
                         {"OpenCLDeviceIndex": str(a.device), "Precision": "mixed"})
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 1234)
    sim.step(int(a.equil_ps / 0.004))
    stride = max(1, int(a.prod_ps / 0.004 / a.nframes))
    pkas = []
    Path(a.out).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(a.out)) as wd:
        for i in range(a.nframes):
            sim.step(stride)
            st = sim.context.getState(getPositions=True)
            fp = os.path.join(wd, f"f{i}.pdb")
            write_frame(sim.topology, st.getPositions(), fp)
            pk = propka_his433(fp, wd)
            if pk is not None:
                pkas.append(pk)
    arr = np.array(pkas) if pkas else np.array([])
    return {"name": name, "pdb": pdb, "his433_pka_frames": [round(float(x), 2) for x in pkas],
            "his433_pka_mean": round(float(arr.mean()), 2) if len(arr) else None,
            "his433_pka_sd": round(float(arr.std()), 2) if len(arr) else None,
            "frac_protonation_favoring": round(float((arr > 6.5).mean()), 2) if len(arr) else None,
            "n_frames": len(pkas), "nframes_req": a.nframes}

def load_list(path):
    rows = []
    with open(path) as fh:
        head = fh.readline()
        cols = head.rstrip("\n").split("\t")
        if "pdb" in cols and "name" in cols:
            ni, pi = cols.index("name"), cols.index("pdb")
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) > max(ni, pi) and c[pi]:
                    rows.append((c[ni], c[pi]))
        else:  # plain "name\tpdb", first line was data
            for line in [head] + fh.readlines():
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2 and c[1]:
                    rows.append((c[0], c[1]))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--nframes", type=int, default=20)
    ap.add_argument("--equil-ps", type=float, default=100.0)
    ap.add_argument("--prod-ps", type=float, default=1000.0)
    ap.add_argument("--restraint-k", type=float, default=2.0)
    a = ap.parse_args()
    if a.nframes < 1: ap.error("--nframes must be >= 1")   # guard: stride divides by nframes (ZeroDivisionError at 0)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in out.glob("*.json")}
    for name, pdb in load_list(a.pdb_list):
        if name in done:
            print(f"skip {name} (done)", flush=True); continue
        if not os.path.exists(pdb):
            print(f"MISSING {name}: {pdb}", flush=True); continue
        try:
            r = run_design(name, pdb, a)
            json.dump(r, open(out / f"{name}.json", "w"), indent=1)
            print(f"{name}: His433 pKa {r['his433_pka_mean']}±{r['his433_pka_sd']} "
                  f"(n={r['n_frames']}, frac>6.5={r['frac_protonation_favoring']})", flush=True)
        except Exception as e:
            print(f"ERROR {name}: {str(e)[:160]}", flush=True)

if __name__ == "__main__":
    main()
