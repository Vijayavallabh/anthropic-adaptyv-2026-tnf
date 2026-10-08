#!/usr/bin/env python3
"""All-vs-all structural diversity among candidate binder folds (local foldseek TMalign, no external DB).
Extracts each binder (chain A) from its complex, runs foldseek all-vs-all, greedily clusters at TM>=THR,
and reports clusters so selection can keep ONE design per distinct SHAPE (addresses 'N designs share one fold').
Usage: fold_diversity.py <name_TAB_pdb.tsv>... --out <dir> [--thr 0.6]
"""
import sys, os, csv, subprocess, argparse, collections
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FS = f"{R}/tools/foldseek_dl/foldseek/bin/foldseek"
ap = argparse.ArgumentParser(); ap.add_argument("tsvs", nargs="+"); ap.add_argument("--out", required=True)
ap.add_argument("--thr", type=float, default=0.6); ap.add_argument("--chain", default="A")
A = ap.parse_args()
bdir = os.path.join(A.out, "binders"); os.makedirs(bdir, exist_ok=True)
names = []
for tsv in A.tsvs:
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        n, pdb = r["name"], r["pdb"]
        out = os.path.join(bdir, n + ".pdb")
        with open(out, "w") as fh:
            for l in open(pdb):
                if l.startswith("ATOM") and l[21] == A.chain: fh.write(l)
            fh.write("TER\nEND\n")
        names.append(n)
m8 = os.path.join(A.out, "allvall.m8")
subprocess.run([FS, "easy-search", bdir, bdir, m8, os.path.join(A.out, "tmp"),
                "--alignment-type", "1", "--tmscore-threshold", "0.0", "-e", "inf",
                "--max-seqs", "2000", "--format-output", "query,target,qtmscore"],
               capture_output=True, text=True)
tm = collections.defaultdict(dict)
for l in open(m8):
    p = l.rstrip().split("\t")
    if len(p) < 3: continue
    q, t, s = p[0].replace(".pdb", ""), p[1].replace(".pdb", ""), float(p[2])
    tm[q][t] = s
# greedy clustering at THR (symmetric max)
def sim(a, b): return max(tm.get(a, {}).get(b, 0.0), tm.get(b, {}).get(a, 0.0))
clusters = []
for n in names:
    placed = False
    for c in clusters:
        if any(sim(n, m) >= A.thr for m in c): c.append(n); placed = True; break
    if not placed: clusters.append([n])
print(f"{len(names)} binders -> {len(clusters)} distinct-fold clusters at TM>={A.thr}\n")
for i, c in enumerate(sorted(clusters, key=len, reverse=True)):
    print(f"cluster {i+1} (n={len(c)}):")
    for n in c:
        short = n.split("__")[-1].replace("_af2pred", "")
        print(f"    {short}")
with open(os.path.join(A.out, "clusters.tsv"), "w") as fh:
    fh.write("cluster\tname\n")
    for i, c in enumerate(sorted(clusters, key=len, reverse=True)):
        for n in c: fh.write(f"{i+1}\t{n}\n")
print(f"\nwrote {A.out}/clusters.tsv  ({len(clusters)} distinct folds among {len(names)})")
