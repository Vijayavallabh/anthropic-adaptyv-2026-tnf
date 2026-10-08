#!/usr/bin/env python3
"""RFdiffusion backbones -> soluble ProteinMPNN (design binder) -> AF2-ig-ready complexes.

RFdiffusion outputs target=chain A (real EGFR DIII seq), binder=chain B (poly-GLY backbone).
This driver, sharded across GPUs:
  1. stages RFd des_*.pdb under unique names (dirs collide on des_<i>.pdb),
  2. runs LigandMPNN (soluble_mpnn) designing ONLY chain B, target A fixed, N seqs each,
  3. reorders each designed complex to binder=chain A / target=chain B (the AF2-ig convention),
     renumbering each chain 1.. and writing to the output dir for scripts/run_af2ig.py.

Usage:
  python scripts/rfd_to_af2ig.py --rfd-dirs 'runs/rfd/w1_*' --out runs/rfd_af2in \
      --work runs/rfd_mpnn --gpus 1-7 --nseq 2 [--limit 0]
Then: python scripts/run_af2ig.py --pdbdir runs/rfd_af2in --out runs/rfd_af2ig --gpus 0-2 --mode fast
"""
from __future__ import annotations
import argparse, glob, json, os, shutil, subprocess, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
L = REPO / "third_party" / "LigandMPNN"


def parse_gpus(s):
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-"); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def reorder(in_pdb, out_pdb, binder_chain="B", target_chain="A"):
    """Write binder (old chain B) first as chain A, then target (old A) as chain B.
    Each chain renumbered from 1; atom serials reset. Returns binder length or None."""
    groups = {"binder": [], "target": []}
    for ln in open(in_pdb):
        if ln.startswith(("ATOM", "HETATM")):
            ch = ln[21]
            if ch == binder_chain:
                groups["binder"].append(ln)
            elif ch == target_chain:
                groups["target"].append(ln)
    if not groups["binder"] or not groups["target"]:
        return None
    out = []
    serial = 1
    for grp, newch in (("binder", "A"), ("target", "B")):
        resmap = {}
        newres = 0
        for ln in groups[grp]:
            old_res = ln[22:27]  # resseq + icode
            if old_res not in resmap:
                newres += 1
                resmap[old_res] = newres
            rn = resmap[old_res]
            nl = "%s%5d %s%s%4d %s" % (
                ln[:6], serial, ln[12:21], newch, rn, ln[27:].rstrip("\n"))
            out.append(nl + "\n")
            serial += 1
        out.append("TER\n")
    out.append("END\n")
    Path(out_pdb).write_text("".join(out))
    return len(resmap) if grp == "target" else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rfd-dirs", required=True, help="glob of RFd output dirs (quoted)")
    ap.add_argument("--out", default=str(REPO / "runs" / "rfd_af2in"))
    ap.add_argument("--work", default=str(REPO / "runs" / "rfd_mpnn"))
    ap.add_argument("--gpus", default="1-7")
    ap.add_argument("--nseq", type=int, default=2)
    ap.add_argument("--temperature", default="0.2")
    ap.add_argument("--variant", default="soluble", choices=["soluble", "vanilla"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--skip-mpnn", action="store_true")
    ap.add_argument("--mpnn-only", action="store_true",
                    help="run MPNN only; skip the backbone-only reorder (use build_fullatom_af2in.py next)")
    ap.add_argument("--binder-chain", default="B",
                    help="chain to redesign (the binder). de novo RFd output: B; partial diffusion: A")
    a = ap.parse_args()

    out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    work = Path(a.work).resolve(); work.mkdir(parents=True, exist_ok=True)
    stage = work / "stage"; stage.mkdir(exist_ok=True)
    gpus = parse_gpus(a.gpus)

    # 1. collect + stage uniquely-named symlinks
    pdbs = []
    for d in sorted(glob.glob(a.rfd_dirs)):
        tag = Path(d).name
        for p in sorted(glob.glob(os.path.join(d, "des_*.pdb"))):
            if "_traj" in os.path.basename(p):
                continue   # skip RFd trajectory frames (write_trajectory=True); keep only final des_<i>.pdb
            pdbs.append((os.path.abspath(p), f"{tag}__{Path(p).stem}"))
    if a.limit:
        pdbs = pdbs[:a.limit]
    if not pdbs:
        sys.exit(f"no des_*.pdb under {a.rfd_dirs}")
    staged = []
    for src, name in pdbs:
        link = stage / f"{name}.pdb"
        if not link.exists():
            try:
                link.symlink_to(src)
            except FileExistsError:
                pass
        staged.append(str(link))
    print(f"staged {len(staged)} RFd backbones", flush=True)

    ckpt = str(L / "model_params" / ("solublempnn_v_48_020.pt" if a.variant == "soluble" else "proteinmpnn_v_48_020.pt"))
    mtype = "soluble_mpnn" if a.variant == "soluble" else "protein_mpnn"
    ckflag = "--checkpoint_soluble_mpnn" if a.variant == "soluble" else "--checkpoint_protein_mpnn"
    py = str(L / ".venv" / "bin" / "python")

    if not a.skip_mpnn:
        # 2. shard across GPUs, one LigandMPNN process each
        k = len(gpus)
        per = (len(staged) + k - 1) // k
        shards = [staged[i * per:(i + 1) * per] for i in range(k)]
        procs = []
        for gpu, sh in zip(gpus, shards):
            if not sh:
                continue
            gwork = work / f"gpu{gpu}"; gwork.mkdir(exist_ok=True)
            js = gwork / "pdbs.json"
            json.dump(sh, open(js, "w"))
            env = {kk: vv for kk, vv in os.environ.items() if kk not in ("VIRTUAL_ENV", "CONDA_PREFIX")}
            env["CUDA_VISIBLE_DEVICES"] = str(gpu)
            cmd = [py, "run.py", "--model_type", mtype, ckflag, ckpt,
                   "--pdb_path_multi", str(js), "--chains_to_design", a.binder_chain,
                   "--out_folder", str(gwork), "--batch_size", str(a.nseq),
                   "--number_of_batches", "1", "--temperature", a.temperature]
            logf = open(work / f"mpnn_gpu{gpu}.log", "w")
            procs.append((subprocess.Popen(cmd, cwd=str(L), env=env, stdout=logf, stderr=subprocess.STDOUT), gpu, logf))
            time.sleep(1.0)
        print(f"launched {len(procs)} MPNN workers over {len(gpus)} gpus", flush=True)
        for pr, gpu, logf in procs:
            pr.wait(); logf.close()

    if a.mpnn_only:
        nbb = len(glob.glob(str(work / "gpu*" / "backbones" / "*.pdb")))
        print(f"MPNN-only: {nbb} designed backbones under {work}/gpu*/backbones\n"
              f"next: build_fullatom_af2in.py --mpnn '{work}/gpu*/backbones/*.pdb' --out <af2in_dir>", flush=True)
        return

    # 3. (legacy, backbone-only target -> af2ig mis-templates; prefer build_fullatom_af2in.py)
    n = 0
    for bb in glob.glob(str(work / "gpu*" / "backbones" / "*.pdb")):
        name = Path(bb).stem
        if reorder(bb, out / f"{name}.pdb") is not None or True:
            if os.path.getsize(out / f"{name}.pdb") > 0:
                n += 1
    print(f"reordered {n} designed complexes -> {out}\nnext: python scripts/run_af2ig.py --pdbdir {out} --out runs/rfd_af2ig --gpus 0-2 --mode fast", flush=True)


if __name__ == "__main__":
    main()
