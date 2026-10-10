#!/bin/bash
# Mature bpair de-novo winners: fold all His-switch variants (runs/bpair_mature/shard{0..3}.csv) vs
# human + mouse mature TNF on all 8 GPUs (shard k human->gpu k, mouse->gpu k+4), score min-ipSAE, and
# flag variants that PRESERVE cross-react vs their own wt (switch candidates for MD/PROPKA/xtb confirm).
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/bpair_mature"
HT="runs/buried2/human_tnf.txt"; MT="runs/buried2/mouse_tnf_mature.txt"
echo "[bpm $(date -u +%FT%TZ)] launching 8 Boltz jobs (4 human gpu0-3, 4 mouse gpu4-7)"
for k in 0 1 2 3; do
  nohup bash scripts/run_boltz_tnf.sh "$M/shard$k.csv" "$HT" bpmh$k $k      > "$M/fold_h$k.log" 2>&1 &
  nohup bash scripts/run_boltz_tnf.sh "$M/shard$k.csv" "$MT" bpmm$k $((k+4)) > "$M/fold_m$k.log" 2>&1 &
done
echo "[bpm] waiting for 8 folds..."
until ( for t in bpmh0 bpmh1 bpmh2 bpmh3 bpmm0 bpmm1 bpmm2 bpmm3; do \
          grep -q "BOLTZ_${t}_DONE" "$M/fold_${t:3:1}${t:4}.log" 2>/dev/null || exit 1; done ); do sleep 20; done
echo "[bpm $(date -u +%FT%TZ)] all folds done; scoring min-ipSAE"
PYS="$R/third_party/LigandMPNN/.venv/bin/python"
"$PYS" - <<'PY'
import os, subprocess, glob, csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"; M=f"{R}/runs/bpair_mature"
def ipsae(did, tags):
    for tag in tags:
        b=f"{R}/runs/boltz_{tag}/out/by_seed/{did}/s0"
        cif=f"{b}/{did}_model_0.cif"
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
HT=["bpmh0","bpmh1","bpmh2","bpmh3"]; MT=["bpmm0","bpmm1","bpmm2","bpmm3"]
res=[]
for r in rows:
    d=r["name"]; h=ipsae(d,HT); m=ipsae(d,MT)
    if h is None or m is None: continue
    res.append((d,h,m,min(h,m)))
# per-design wt baseline
wt={d[:-3]:mn for (d,h,m,mn) in res if d.endswith("_wt")}
import json
out=[]
print(f"{'variant':20}{'h':>7}{'m':>7}{'min':>7}{'wt_min':>8}{'dmin':>7}  switch_cand")
for d,h,m,mn in sorted(res,key=lambda x:-x[3]):
    base=None
    for pfx,w in wt.items():
        if d.startswith(pfx+"_H") or d==pfx+"_wt": base=w; break
    dmin=(mn-base) if base is not None else None
    # switch candidate: a His variant that KEEPS cross-react (min within 0.08 of its wt, and >=0.45 absolute)
    is_his = ("_H" in d) and not d.endswith("_wt")
    cand = is_his and base is not None and mn>=0.45 and (mn-base)>=-0.08
    out.append({"variant":d,"h":h,"m":m,"min":mn,"wt_min":base,"dmin":dmin,"switch_cand":bool(cand)})
    bs=f"{base:.3f}" if base is not None else "  -  "
    ds=f"{dmin:+.3f}" if dmin is not None else "  -  "
    print(f"{d:20}{h:>7.3f}{m:>7.3f}{mn:>7.3f}{bs:>8}{ds:>7}  {'YES' if cand else ''}")
json.dump(out, open(f"{M}/scored.json","w"), indent=1)
cands=[o for o in out if o["switch_cand"]]
print(f"\nswitch candidates (keep cross-react, min>=0.45): {len(cands)}")
for o in sorted(cands,key=lambda x:-x['min'])[:12]: print(f"  {o['variant']:18} min={o['min']:.3f} (wt {o['wt_min']:.3f}, d{o['dmin']:+.3f})")
PY
echo "BPM_SCORE_DONE $(date -u +%FT%TZ)"
