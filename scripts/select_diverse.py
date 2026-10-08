#!/usr/bin/env python3
"""Diversity-aware final selection for the Adaptyv EGFR resubmit.

User directive (2026-10-03): "pick the 12 new ones to spread across distinct designs, not just the
highest-scoring cousins of the same few." Reviewers choose designs for the wet lab, so a set of
DISTINCT folds/lineages hedges risk better than near-clones of one winning seed. Partial diffusion
at T=20 stays near its seed, so same-seed neighbors are structurally redundant.

Input --pool CSV needs: name, and a score column (higher=better by default; --lower-better to flip).
'sequence' and any other columns are carried through. --pdbs <dir> of binder monomers enables
foldseek structural clustering (distinct folds); without it each seed lineage is its own cluster.

Greedy: sort by score, take the best candidate from each still-unrepresented (lineage AND cluster);
once coverage is exhausted, fill the rest by score while capping per lineage. Writes --out (the chosen
N, diversity-first) and <out>.pool.csv (full pool annotated with lineage/cluster/rank/selected) so you
can swap on the portal while keeping one-per-cluster.
"""
import argparse, csv, os, re, subprocess, collections

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FS = f"{REPO}/tools/novelty/foldseek/bin/foldseek"
LINEAGE_RE = re.compile(r"(egfr_phsw_\d+)")


def lineage(name):
    m = LINEAGE_RE.search(name)
    return m.group(1) if m else name


def struct_clusters(pdbs_dir, work, cluster_tm, threads):
    """foldseek easy-cluster -> {design_basename: cluster_rep}; {} on failure (caller falls back)."""
    if not pdbs_dir or not os.path.isdir(pdbs_dir):
        return {}
    os.makedirs(work, exist_ok=True)
    pref = f"{work}/clu"
    subprocess.run([FS, "easy-cluster", pdbs_dir, pref, f"{work}/tmp",
                    "--alignment-type", "1", "--tmscore-threshold", str(cluster_tm),
                    "-c", "0.7", "--threads", str(threads)],
                   cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = f"{pref}_cluster.tsv"
    clu = {}
    if os.path.exists(out):
        for line in open(out):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep = os.path.splitext(os.path.basename(parts[0]))[0]
            mem = os.path.splitext(os.path.basename(parts[1]))[0]
            clu[mem] = rep
    return clu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--pdbs", default=None)
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--max-per-lineage", type=int, default=2)
    ap.add_argument("--score-col", default="score")
    ap.add_argument("--lower-better", action="store_true")
    ap.add_argument("--cluster-tm", type=float, default=0.6)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--out", default=f"{REPO}/results/diverse_selected.csv")
    ap.add_argument("--work", default=f"{REPO}/runs/diverse_sel")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.pool)))
    sgn = 1.0 if a.lower_better else -1.0   # sort ascending by sgn*score -> best first
    kept = []
    for r in rows:
        try:
            r["_score"] = float(r.get(a.score_col, ""))
        except (TypeError, ValueError):
            continue                         # drop rows without a usable score
        r["_lin"] = lineage(r["name"])
        kept.append(r)

    clu = struct_clusters(a.pdbs, a.work, a.cluster_tm, a.threads)
    for r in kept:
        r["_clu"] = clu.get(r["name"], r["_lin"])   # fallback: lineage is its own cluster

    kept.sort(key=lambda r: sgn * r["_score"])

    chosen, picked = [], set()
    lin_ct = collections.Counter()
    seen_lin, seen_clu = set(), set()
    # round 1: diversity-first -- best from each still-unrepresented lineage OR cluster
    for r in kept:
        if len(chosen) >= a.n:
            break
        if r["_lin"] in seen_lin and r["_clu"] in seen_clu:
            continue
        if lin_ct[r["_lin"]] >= a.max_per_lineage:
            continue
        chosen.append(r); picked.add(id(r))
        lin_ct[r["_lin"]] += 1; seen_lin.add(r["_lin"]); seen_clu.add(r["_clu"])
    # round 2: fill remaining by score, still capping per lineage
    for r in kept:
        if len(chosen) >= a.n:
            break
        if id(r) in picked or lin_ct[r["_lin"]] >= a.max_per_lineage:
            continue
        chosen.append(r); picked.add(id(r)); lin_ct[r["_lin"]] += 1

    base = list(rows[0].keys()) if rows else ["name"]
    cols = base + [c for c in ("_lin", "_clu", "_score") if c not in base] + ["rank", "selected"]

    def write(path, items, mark):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for i, r in enumerate(items, 1):
                r2 = dict(r)
                r2["rank"] = i
                r2["selected"] = "yes" if (mark is None or id(r) in mark) else ""
                w.writerow(r2)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    write(a.out, chosen, None)
    write(a.out + ".pool.csv", kept, picked)

    nl = len({r["_lin"] for r in chosen})
    nc = len({r["_clu"] for r in chosen})
    print(f"selected {len(chosen)}/{a.n} | distinct lineages={nl} | distinct clusters={nc} "
          f"| clustering={'foldseek' if clu else 'lineage-fallback'}")
    for r in chosen:
        print(f"  {r['name']:42s} lin={r['_lin']:16s} clu={r['_clu']:30s} score={r['_score']:.3f}")


if __name__ == "__main__":
    main()
