#!/usr/bin/env python3
"""Objective-1 engine, NETWORK variant: install a *bidentate ionizable network* around His433.

The baseline obj1_switch.py installs a SINGLE binder carboxylate against His433. This explores the
Boyken-2019-style idea that a protonated His (+1 at tumor pH 6.5) can be co-coordinated by TWO binder
carboxylates simultaneously -> a cooperative, deeper proton-linked switch that is more robust than a
lone salt bridge (and harder to satisfy by a PROPKA desolvation artifact, since two real contacts must
land). For each binder backbone we scan positions near His433 and score, via the SAME PROPKA Wyman
linkage + conserved-sensor restriction + pKa clamps as obj1_switch.py:
  - single acids  D@p, E@p                        (baseline, for the network-vs-single comparison)
  - acid PAIRS    {D,E}@p1 + {D,E}@p2             (the network: two carboxylates co-facing His433)
We keep the best variant per design and flag whether the NETWORK beat the best single. A network win
must be RE-VALIDATED downstream (AF2-ig binding + MD occupancy on the engineered seq) before use --
installing two charges is more likely to perturb binding, so the gate is binding first.

Usage (driver shards across PyRosetta workers):
  envs/pyrosetta/bin/python scripts/obj1_switch_network.py --pdbs 'runs/rfd_passers_w2fa/*.pdb' \
      --outdir runs/obj1_net --out results/rfd_obj1_network.csv --workers 40
Run with envs/pyrosetta activated.
"""
from __future__ import annotations
import argparse, csv, glob, os, subprocess, sys, time
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import obj1_switch as O  # reuse the validated machinery

FIELDS = ["design", "best_tag", "best_is_network", "ddg_pH", "ddg_pH_single_best",
          "ddg_pH_total", "ddg_pH_his433", "dn_6.5", "dn_7.4", "his433_cx_pka", "his433_dpKa",
          "ph_drivers", "interface_ddg", "n_variants", "orig_ddg_pH",
          "binder_chain", "target_chain", "binder_len", "binder_seq", "best_pdb", "src_pdb", "error"]


def process_network(pyr, sf, pdb, outdir, acid_cut=10.0, max_acid=6, pair_top=5):
    name = Path(pdb).stem
    try:
        pose = O.load_and_relax(pyr, sf, pdb)
    except Exception as e:
        return {"design": name, "error": f"load:{repr(e)[:80]}"}
    sizes = O.chain_sizes(pose)
    if len(sizes) < 2:
        return {"design": name, "error": "fewer than 2 chains"}
    target_chain = max(sizes, key=sizes.get); binder_chain = min(sizes, key=sizes.get)
    his = O.find_res(pose, target_chain, O.HIS433_NUM)
    if his is None:
        return {"design": name, "error": "His433 not found"}
    hr = pose.residue(his)
    his_centers = [hr.xyz(a) for a in ("ND1", "NE2") if hr.has(a)]
    acid_pos = O.near_positions(pose, binder_chain, his_centers, acid_cut, max_acid)
    if len(acid_pos) < 2:
        return {"design": name, "error": "fewer than 2 acid positions near His433"}

    variants = [("orig", [], [his])]
    for (i, pnum, aa, d) in acid_pos:                       # singles (baseline)
        variants.append((f"D@{pnum}", [(i, "D")], [his]))
        variants.append((f"E@{pnum}", [(i, "E")], [his]))
    top = acid_pos[:pair_top]                               # pairs (the network)
    for (i1, p1, _, _), (i2, p2, _, _) in combinations(top, 2):
        for a1 in ("D", "E"):
            for a2 in ("D", "E"):
                variants.append((f"NET[{a1}@{p1}+{a2}@{p2}]",
                                 [(i1, a1), (i2, a2)], [his]))

    scored = []
    for tag, muts, foci in variants:
        try:
            p = O.make_variant(pyr, pose, sf, muts, foci)
        except Exception:
            continue
        tmp = os.path.join(outdir, f".{name}__{tag}.pdb".replace("@", "at").replace("[", "").replace("]", "").replace("+", "_"))
        p.dump_pdb(tmp)
        sc = O.ddg_ph_from_pdb(tmp, binder_chain, target_chain)
        os.remove(tmp)
        if sc is not None:
            scored.append((tag, muts, foci, sc))
    if not scored:
        return {"design": name, "error": "no scored variants"}

    best = max(scored, key=lambda s: s[3]["ddg_pH"])
    best_single = max((s for s in scored if not s[0].startswith("NET")),
                      key=lambda s: s[3]["ddg_pH"], default=None)
    tag, muts, foci, sc = best
    p = O.make_variant(pyr, pose, sf, muts, foci) if muts or tag == "orig" else pose
    outpdb = os.path.join(outdir, f"{name}__best.pdb")
    p.dump_pdb(outpdb)
    iddg = O.interface_ddg(pyr, p)
    seq = O.binder_seq(p, binder_chain)
    return {
        "design": name, "best_tag": tag, "best_is_network": tag.startswith("NET"),
        "ddg_pH": sc["ddg_pH"],
        "ddg_pH_single_best": (best_single[3]["ddg_pH"] if best_single else None),
        "ddg_pH_total": sc.get("ddg_pH_total"), "ddg_pH_his433": sc.get("ddg_pH_his433"),
        "dn_6.5": sc["dn_6.5"], "dn_7.4": sc["dn_7.4"],
        "his433_cx_pka": sc["his433_cx_pka"], "his433_dpKa": sc["his433_dpKa"],
        "ph_drivers": sc.get("ph_drivers", ""), "interface_ddg": iddg, "n_variants": len(scored),
        "orig_ddg_pH": next((s[3]["ddg_pH"] for s in scored if s[0] == "orig"), None),
        "binder_chain": binder_chain, "target_chain": target_chain,
        "binder_seq": seq, "binder_len": len(seq), "best_pdb": outpdb, "src_pdb": pdb, "error": "",
    }


def run_worker(shard_file, outdir, out_csv):
    pyr = O.init_pyrosetta(); sf = pyr.get_fa_scorefxn(); os.makedirs(outdir, exist_ok=True)
    pdbs = [l.strip() for l in open(shard_file) if l.strip()]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader()
        for pdb in pdbs:
            t0 = time.time()
            try:
                row = process_network(pyr, sf, pdb, outdir)
            except Exception as e:
                row = {"design": Path(pdb).stem, "error": f"proc:{repr(e)[:100]}"}
            w.writerow({k: row.get(k, "") for k in FIELDS}); fh.flush()
            print(f"[{Path(pdb).stem}] tag={row.get('best_tag')} net={row.get('best_is_network')} "
                  f"ddg_pH={row.get('ddg_pH')} single={row.get('ddg_pH_single_best')} "
                  f"err={row.get('error')} {time.time()-t0:.0f}s", flush=True)


def run_driver(pdbs, outdir, out_csv, workers):
    os.makedirs(outdir, exist_ok=True)
    shard_dir = os.path.join(outdir, "shards"); os.makedirs(shard_dir, exist_ok=True)
    workers = max(1, min(workers, len(pdbs)))
    shards = [[] for _ in range(workers)]
    for i, p in enumerate(pdbs): shards[i % workers].append(p)
    procs = []
    for k, sh in enumerate(shards):
        if not sh: continue
        sf_file = os.path.join(shard_dir, f"shard_{k:03d}.txt"); Path(sf_file).write_text("\n".join(sh) + "\n")
        pcsv = os.path.join(shard_dir, f"part_{k:03d}.csv")
        log = open(os.path.join(shard_dir, f"shard_{k:03d}.log"), "w")
        cmd = [sys.executable, __file__, "--worker", "--shard", sf_file, "--outdir", outdir, "--out", pcsv]
        procs.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), pcsv, log))
        time.sleep(0.3)
    print(f"launched {len(procs)} workers over {len(pdbs)} designs -> {outdir}", flush=True)
    for pr, _, _ in procs: pr.wait()
    rows = []
    for _, pcsv, log in procs:
        log.close()
        if os.path.exists(pcsv): rows += list(csv.DictReader(open(pcsv)))
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, "") for k in FIELDS})
    ok = [r for r in rows if not r.get("error")]
    nets = [r for r in ok if str(r.get("best_is_network")) == "True"]
    print(f"\nDONE {len(ok)}/{len(rows)} scored -> {out_csv} | network beat single in {len(nets)}/{len(ok)}", flush=True)
    for r in sorted(ok, key=lambda r: float(r["ddg_pH"] or -9), reverse=True)[:15]:
        print(f"  {r['design'][:34]:34s} tag={r['best_tag']:22s} ddg_pH={r['ddg_pH']} "
              f"(single {r['ddg_pH_single_best']}) dpKa={r['his433_dpKa']} len={r['binder_len']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true"); ap.add_argument("--shard")
    ap.add_argument("--pdbs", nargs="*"); ap.add_argument("--outdir", default=str(REPO / "runs" / "obj1_net"))
    ap.add_argument("--out", default=str(REPO / "results" / "rfd_obj1_network.csv"))
    ap.add_argument("--workers", type=int, default=40); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.worker:
        run_worker(a.shard, a.outdir, a.out); return
    pdbs = []
    for p in (a.pdbs or []):
        pdbs += sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p]
    pdbs = [p for p in pdbs if os.path.exists(p)]
    if a.limit: pdbs = pdbs[:a.limit]
    if not pdbs: sys.exit("no input PDBs")
    run_driver(pdbs, a.outdir, a.out, a.workers)


if __name__ == "__main__":
    main()
