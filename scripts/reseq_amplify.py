#!/usr/bin/env python3
"""Amplify good-fold binder backbones: re-sequence the binder (chain A) with soluble ProteinMPNN
(LigandMPNN), keeping the target (chain B) fixed, to turn each proven backbone into a family of
candidate sequences. BoltzGen's own inverse folding is the weak link (~2-4% AF2 bind rate);
ProteinMPNN sequences fold/bind far better. Collects each designed-backbone PDB (designed sequence
on the backbone + fixed target) into one dir for scripts/run_af2ig.py.
"""
from __future__ import annotations
import argparse, glob, os, shutil, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
L = REPO / "third_party" / "LigandMPNN"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True, help="file of design names (one per line)")
    ap.add_argument("--n", type=int, default=48, help="top-N backbones from the list")
    ap.add_argument("--nseq", type=int, default=8)
    ap.add_argument("--srcdir", default=str(REPO / "runs" / "obj1"), help="where <name>__best.pdb live")
    ap.add_argument("--outdir", default=str(REPO / "runs" / "reseq_foldin"))
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--temperature", default="0.2")
    ap.add_argument("--variant", default="soluble", choices=["soluble", "vanilla"])
    a = ap.parse_args()

    names = [l.strip() for l in open(a.list) if l.strip()][:a.n]
    out = Path(a.outdir).resolve()        # ABSOLUTE: run.py runs with cwd=LigandMPNN, so a relative --out_folder lands there
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_mpnn"
    work.mkdir(exist_ok=True)
    ckpt = str(L / "model_params" / ("solublempnn_v_48_020.pt" if a.variant == "soluble" else "proteinmpnn_v_48_020.pt"))
    mtype = "soluble_mpnn" if a.variant == "soluble" else "protein_mpnn"
    ckflag = "--checkpoint_soluble_mpnn" if a.variant == "soluble" else "--checkpoint_protein_mpnn"
    py = str(L / ".venv" / "bin" / "python")
    env = {k: v for k, v in os.environ.items() if k not in ("VIRTUAL_ENV", "CONDA_PREFIX")}
    env["CUDA_VISIBLE_DEVICES"] = a.gpu

    collected = 0
    srcdir = os.path.abspath(a.srcdir)
    for nm in names:
        pdb = os.path.join(srcdir, nm + "__best.pdb")
        if not os.path.exists(pdb):
            print(f"skip {nm}: no {pdb}", flush=True)
            continue
        wd = work / nm
        cmd = [py, "run.py", "--model_type", mtype, ckflag, ckpt,
               "--pdb_path", pdb, "--chains_to_design", "A", "--out_folder", str(wd),
               "--batch_size", str(a.nseq), "--number_of_batches", "1", "--temperature", a.temperature]
        r = subprocess.run(cmd, cwd=str(L), env=env, capture_output=True, text=True)
        bbs = glob.glob(str(wd / "backbones" / "*.pdb"))
        for bb in bbs:
            i = Path(bb).stem.rsplit("_", 1)[-1]
            shutil.copyfile(bb, out / f"{nm}__mpnn{i}.pdb")
            collected += 1
        msg = f"{nm}: {len(bbs)} seqs"
        if not bbs:
            msg += f"  [rc={r.returncode}] {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ''}"
        print(msg, flush=True)
    print(f"\ncollected {collected} designed PDBs -> {out}  (next: run_af2ig.py --pdbdir {out})")


if __name__ == "__main__":
    main()
