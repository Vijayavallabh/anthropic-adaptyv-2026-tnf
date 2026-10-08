#!/usr/bin/env python3
"""Calibrated local novelty oracle for the Adaptyv/Proteinbase novelty score.

The portal folds each submitted SEQUENCE with a single-sequence model (ESMFold2) and
scores novelty = MMseqs2 identity vs SwissProt+PDB (sequence axis) + Foldseek TMalign-mode
max TM vs PDB+AFDB (structure axis). Pass = level >= 3/4.

The OLD proxy (novelty_proxy.py) reused the af2ig COMPLEX structure. That structure is too
PDB-canonical: its Foldseek TM reads ~0.05-0.1 high and CANNOT resolve the 2-vs-3 boundary
where every current design sits (ground truth: design 22 passed at af2ig-TM 0.89 while 24
FAILED at 0.892). This oracle fixes the one broken step: it folds the BARE SEQUENCE with a
single-sequence ESM folder (classic ESMFold facebook/esmfold_v1, cached locally), then runs
the identical Foldseek + MMseqs scoring. The TM threshold is calibrated against the 20 portal
ground-truth labels (--calibrate); classic ESMFold != ESMFold2 in absolute numerics, so the
empirical threshold on ground truth -- not the published 0.8 -- is what we trust, with a margin.

Usage:
  CUDA_VISIBLE_DEVICES=0 scripts/novelty_oracle.py \
      --seqs results/adaptyv_challenge1_resubmit_W19_FINAL.csv \
      --struct-dir runs/novelty_oracle/esmfold --out results/novelty_oracle_20.csv --calibrate
"""
import argparse, csv, os, sys, glob, time, warnings
warnings.filterwarnings("ignore")

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
sys.path.insert(0, os.path.join(REPO, "scripts"))
import novelty_proxy as npx  # reuse seq_identity(), struct_tm(), level(), classify()

# Portal ground truth from the 2026-10-05 novelty check of resubmit_W19_FINAL (user paste).
GT_PASS = {f"egfr_phsw_{i}" for i in ["02","04","06","08","09","12","16","19","22","23","27","29","32"]}
GT_FAIL = {f"egfr_phsw_{i}" for i in ["21","24","25","26","28","30","31"]}


def read_seqs(path):
    out = {}
    for r in csv.DictReader(open(path)):
        nm = r.get("name") or r.get("design")
        sq = (r.get("sequence") or r.get("binder_seq") or r.get("seq") or "").split(":")[0].strip().upper()
        if nm and sq:
            out[nm] = sq
    return out


def fold_esmfold(seqs, struct_dir, overwrite=False):
    """Fold each sequence with classic ESMFold (facebook/esmfold_v1) -> <struct_dir>/<name>.pdb."""
    os.makedirs(struct_dir, exist_ok=True)
    todo = {n: s for n, s in seqs.items()
            if overwrite or not os.path.exists(f"{struct_dir}/{n}.pdb")}
    if not todo:
        print(f"[fold] all {len(seqs)} structures present in {struct_dir}, skipping fold")
        return
    print(f"[fold] loading ESMFold facebook/esmfold_v1 ({len(todo)} to fold)...", flush=True)
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.cuda().eval()
    model.esm = model.esm.half()
    torch.backends.cuda.matmul.allow_tf32 = True
    try:
        model.trunk.set_chunk_size(64)
    except Exception:
        pass
    t0 = time.time()
    for i, (n, s) in enumerate(sorted(todo.items()), 1):
        try:
            with torch.no_grad():
                pdb = model.infer_pdb(s)
            with open(f"{struct_dir}/{n}.pdb", "w") as fh:
                fh.write(pdb)
            print(f"[fold] {i}/{len(todo)} {n} len={len(s)} {time.time()-t0:.0f}s", flush=True)
        except Exception as e:
            print(f"[fold] FAILED {n}: {repr(e)[:160]}", flush=True)
    del model
    torch.cuda.empty_cache()


def fold_esmfold2(seqs, struct_dir, overwrite=False):
    """Fold each sequence with the EXACT portal model, ESMFold2 (biohub/ESMFold2) -> <struct_dir>/<name>.pdb.
    Requires LD_LIBRARY_PATH to the env's nvidia wheel libs (libcublas.so.13 etc.) and HF_HOME set."""
    os.makedirs(struct_dir, exist_ok=True)
    todo = {n: s for n, s in seqs.items()
            if overwrite or not os.path.exists(f"{struct_dir}/{n}.pdb")}
    if not todo:
        print(f"[fold2] all {len(seqs)} structures present in {struct_dir}, skipping", flush=True); return
    print(f"[fold2] loading ESMFold2 biohub/ESMFold2 ({len(todo)} to fold)...", flush=True)
    import torch
    from esm.models.esmfold2 import EsmFold2Model
    model = EsmFold2Model.from_pretrained("biohub/ESMFold2").cuda().eval()
    t0 = time.time()
    for i, (n, s) in enumerate(sorted(todo.items()), 1):
        try:
            with torch.no_grad():
                pdb = model.infer_protein_as_pdb(s)
            with open(f"{struct_dir}/{n}.pdb", "w") as fh:
                fh.write(pdb)
            print(f"[fold2] {i}/{len(todo)} {n} len={len(s)} {time.time()-t0:.0f}s", flush=True)
        except Exception as e:
            print(f"[fold2] FAILED {n}: {repr(e)[:180]}", flush=True)
    del model; torch.cuda.empty_cache()


def calibrate(rows):
    labeled = [r for r in rows if r["name"] in GT_PASS or r["name"] in GT_FAIL]
    if not labeled:
        return
    print("\n===== CALIBRATION vs 20 portal labels =====")
    print(f"{'name':16s} {'seqid':>6s} {'TM':>6s} {'level':>5s} {'portal':>7s} {'proxyPass':>9s}")
    for r in sorted(labeled, key=lambda x: x["struct_tm"]):
        truth = "PASS" if r["name"] in GT_PASS else "FAIL"
        print(f"{r['name']:16s} {r['seqid']:6.1f} {r['struct_tm']:6.3f} {r['level']:5d} {truth:>7s} {str(r['pass']):>9s}")
    # best single-TM threshold (pass iff TM < thr), assuming seqid~0 for all
    tms = sorted(set(round(r["struct_tm"], 3) for r in labeled))
    cands = [(tms[i] + tms[i + 1]) / 2 for i in range(len(tms) - 1)] + [t + 0.001 for t in tms]
    best = None
    for thr in sorted(set(cands)):
        tp = sum(1 for r in labeled if r["name"] in GT_PASS and r["struct_tm"] < thr)
        tn = sum(1 for r in labeled if r["name"] in GT_FAIL and r["struct_tm"] >= thr)
        acc = tp + tn
        # also require no FAIL below thr (no false "novel") as tie-break -> conservative
        fp = sum(1 for r in labeled if r["name"] in GT_FAIL and r["struct_tm"] < thr)
        if best is None or (acc, -fp) > (best[1], -best[2]):
            best = (thr, acc, fp)
    thr, acc, fp = best
    npass, nfail = len(GT_PASS & {r["name"] for r in labeled}), len(GT_FAIL & {r["name"] for r in labeled})
    print(f"\nbest TM threshold = {thr:.3f} (pass iff TM < thr): {acc}/{len(labeled)} correct "
          f"({npass} pass / {nfail} fail labeled); false-novel (FAIL below thr) = {fp}")
    passing_tm = [r["struct_tm"] for r in labeled if r["name"] in GT_PASS]
    failing_tm = [r["struct_tm"] for r in labeled if r["name"] in GT_FAIL]
    if passing_tm and failing_tm:
        print(f"  PASS TM range [{min(passing_tm):.3f}, {max(passing_tm):.3f}]  "
              f"FAIL TM range [{min(failing_tm):.3f}, {max(failing_tm):.3f}]")
        print(f"  separation margin = {min(failing_tm) - max(passing_tm):+.3f} "
              f"(positive = cleanly separable by TM alone)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--struct-dir", default=f"{REPO}/runs/novelty_oracle/esmfold")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=f"{REPO}/runs/novelty_oracle/work")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--no-fold", action="store_true", help="structures already present")
    ap.add_argument("--folder", default="esmfold_v1", choices=["esmfold_v1", "esmfold2"],
                    help="esmfold_v1 = classic ESMFold (calibrated proxy); esmfold2 = the EXACT portal model")
    a = ap.parse_args()

    seqs = read_seqs(a.seqs)
    print(f"[oracle] {len(seqs)} sequences from {a.seqs} | folder={a.folder}")
    if not a.no_fold:
        (fold_esmfold2 if a.folder == "esmfold2" else fold_esmfold)(seqs, a.struct_dir, a.overwrite)

    seqid = npx.seq_identity(seqs, a.work, a.threads)
    tm = npx.struct_tm(a.struct_dir, a.work, a.threads)
    scored = [n for n in seqs if n in tm]
    if len(scored) < 0.5 * len(seqs):
        sys.exit(f"ORACLE_FAIL: foldseek scored only {len(scored)}/{len(seqs)} (tool failure?)")

    rows = []
    for n in seqs:
        s = seqid.get(n, 0.0); t = tm.get(n, 0.0); lv = npx.level(s, t)
        rows.append({"name": n, "seqid": round(s, 1), "struct_tm": round(t, 3),
                     "level": lv, "class": npx.classify(lv, t), "pass": lv >= 3})
    rows.sort(key=lambda r: (-r["pass"], r["struct_tm"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "seqid", "struct_tm", "level", "class", "pass"])
        w.writeheader(); w.writerows(rows)
    print(f"[oracle] wrote {a.out}: {len(rows)} designs | level>=3: {sum(r['pass'] for r in rows)}")
    if a.calibrate:
        calibrate(rows)


if __name__ == "__main__":
    main()
