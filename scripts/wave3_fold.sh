#!/bin/bash
# Generic fold+ipSAE scorer. Usage: wave3_fold.sh <rundir> <tagpfx>  (default runs/wave3 wv3)
# Folds <rundir>/shard{0..3}.csv vs human+mouse TNF on 8 GPUs (tags <tag>h0-3 / <tag>m0-3),
# writes <rundir>/scored.json = [{variant,h,m,min,dHM}]. (dual/switch flags computed downstream.)
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="${1:-$R/runs/wave3}"; TAG="${2:-wv3}"
HT="runs/buried2/human_tnf.txt"; MT="runs/buried2/mouse_tnf_mature.txt"
echo "[$TAG $(date -u +%FT%TZ)] launching 8 Boltz jobs (dir $M)"
for k in 0 1 2 3; do
  nohup bash scripts/run_boltz_tnf.sh "$M/shard$k.csv" "$HT" ${TAG}h$k $k      > "$M/fold_${TAG}h$k.log" 2>&1 &
  nohup bash scripts/run_boltz_tnf.sh "$M/shard$k.csv" "$MT" ${TAG}m$k $((k+4)) > "$M/fold_${TAG}m$k.log" 2>&1 &
done
echo "[$TAG] waiting for 8 folds..."
until ( for t in ${TAG}h0 ${TAG}h1 ${TAG}h2 ${TAG}h3 ${TAG}m0 ${TAG}m1 ${TAG}m2 ${TAG}m3; do \
          grep -q "BOLTZ_${t}_DONE" "$M/fold_${t}.log" 2>/dev/null || exit 1; done ); do sleep 20; done
echo "[$TAG $(date -u +%FT%TZ)] folds done; scoring ipSAE"
"$R/third_party/LigandMPNN/.venv/bin/python" - "$M" "$TAG" <<'PY'
import os, subprocess, glob, csv, json, sys
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
M=sys.argv[1]; TAG=sys.argv[2]
def ipsae(did, tags):
    for tag in tags:
        b=f"{R}/runs/boltz_{tag}/out/by_seed/{did}/s0"; cif=f"{b}/{did}_model_0.cif"
        if not os.path.exists(cif): continue
        out=cif.replace(".cif","_10_10.txt")
        if not os.path.exists(out):
            npz=glob.glob(f"{b}/pae_{did}_model_0.npz")
            if npz: subprocess.run([PY,f"{R}/tools/ipsae.py",os.path.basename(npz[0]),f"{did}_model_0.cif","10","10"],capture_output=True,cwd=b)
        if not os.path.exists(out): continue
        a=[float(l.split()[5]) for l in open(out) if len(l.split())>6 and l.split()[4]=="asym"]
        if a: return min(a)
    return None
rows=list(csv.DictReader(open(f"{M}/variants.csv")))
HS=[f"{TAG}h{k}" for k in range(4)]; MS=[f"{TAG}m{k}" for k in range(4)]
out=[]
for r in rows:
    d=r["name"]; h=ipsae(d,HS); m=ipsae(d,MS)
    if h is None or m is None: continue
    out.append(dict(variant=d,h=h,m=m,min=min(h,m),dHM=abs(h-m)))
json.dump(out, open(f"{M}/scored.json","w"), indent=1)
print(f"[{TAG}] scored {len(out)}/{len(rows)}")
PY
echo "${TAG}_SCORE_DONE $(date -u +%FT%TZ)"
