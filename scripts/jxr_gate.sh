#!/bin/bash
# Wait for jxr human folds, run mouse folds, score ipSAE_min, gate vs des59_H48 baseline (min 0.577).
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
until grep -q "BOLTZ_jxrh0_DONE" runs/jxr/fold_h0.log 2>/dev/null && \
      grep -q "BOLTZ_jxrh1_DONE" runs/jxr/fold_h1.log 2>/dev/null && \
      grep -q "BOLTZ_jxrh2_DONE" runs/jxr/fold_h2.log 2>/dev/null; do sleep 30; done
echo "human folds done; launching mouse folds"
nohup bash scripts/run_boltz_tnf.sh runs/jxr/jxr_top_0.csv runs/buried2/mouse_tnf_mature.txt jxrm0 5 > runs/jxr/fold_m0.log 2>&1 &
nohup bash scripts/run_boltz_tnf.sh runs/jxr/jxr_top_1.csv runs/buried2/mouse_tnf_mature.txt jxrm1 6 > runs/jxr/fold_m1.log 2>&1 &
nohup bash scripts/run_boltz_tnf.sh runs/jxr/jxr_top_2.csv runs/buried2/mouse_tnf_mature.txt jxrm2 7 > runs/jxr/fold_m2.log 2>&1 &
wait
until grep -q "BOLTZ_jxrm0_DONE" runs/jxr/fold_m0.log 2>/dev/null && \
      grep -q "BOLTZ_jxrm1_DONE" runs/jxr/fold_m1.log 2>/dev/null && \
      grep -q "BOLTZ_jxrm2_DONE" runs/jxr/fold_m2.log 2>/dev/null; do sleep 30; done
echo "mouse folds done; scoring"
third_party/LigandMPNN/.venv/bin/python - <<'PY'
import os, subprocess, glob, csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
def ipsae(did, tags):
    for tag in tags:
        bs=glob.glob(f"{R}/runs/boltz_{tag}/out/by_seed/{did}/s0")
        if not bs: continue
        b=bs[0]; cif=f"{b}/{did}_model_0.cif"
        if not os.path.exists(cif): continue
        out=cif.replace(".cif","_10_10.txt")
        if not os.path.exists(out):
            npz=glob.glob(f"{b}/pae_{did}_model_0.npz")
            if npz: subprocess.run([PY,f"{R}/tools/ipsae.py",os.path.basename(npz[0]),f"{did}_model_0.cif","10","10"],capture_output=True,cwd=b)
        if not os.path.exists(out): continue
        a=[float(l.split()[5]) for l in open(out) if len(l.split())>6 and l.split()[4]=="asym"]
        if a: return min(a)
    return None
rows=list(csv.DictReader(open(f"{R}/runs/jxr/jxr_top.csv")))
res=[]
for r in rows:
    d=r["name"]; h=ipsae(d,["jxrh0","jxrh1","jxrh2"]); m=ipsae(d,["jxrm0","jxrm1","jxrm2"])
    if h is None or m is None: continue
    res.append((d,h,m,min(h,m)))
res.sort(key=lambda x:-x[3])
print(f"{'design':22}{'h':>7}{'m':>7}{'min':>7}")
for d,h,m,mn in res[:20]: print(f"{d:22}{h:>7.3f}{m:>7.3f}{mn:>7.3f}")
base=0.577
win=[r for r in res if r[3]>base]
print(f"\ndes59_H48 baseline min=0.577 ; v4 parent b1m9_A12H min=0.387")
print(f"BEAT 0.577: {len(win)} -> {[(r[0],round(r[3],3)) for r in win[:8]]}")
import json; json.dump([{'name':r[0],'h':r[1],'m':r[2],'min':r[3]} for r in res], open(f"{R}/runs/jxr/jxr_scored.json","w"), indent=1)
PY
echo "JXR_GATE_DONE $(date -u +%FT%TZ)"
