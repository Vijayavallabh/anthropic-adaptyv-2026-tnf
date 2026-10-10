#!/bin/bash
# WAVE-2 maturation: fold all His-switch variants (runs/bpair_mature2/shard{0..3}.csv) vs human+mouse TNF
# on 8 GPUs, score min-ipSAE AND Δ_hm (|h-m|), and flag DUAL-STRONG candidates = keep cross-react vs wt
# AND balanced two-species binding (min>=0.52 AND Δ_hm<=0.13) per the stress-test gate. Then MD-confirm.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/bpair_mature2"; HT="runs/buried2/human_tnf.txt"; MT="runs/buried2/mouse_tnf_mature.txt"
echo "[bpm2 $(date -u +%FT%TZ)] launching 8 Boltz jobs"
for k in 0 1 2 3; do
  nohup bash scripts/run_boltz_tnf.sh "$M/shard$k.csv" "$HT" bpm2h$k $k      > "$M/fold_h$k.log" 2>&1 &
  nohup bash scripts/run_boltz_tnf.sh "$M/shard$k.csv" "$MT" bpm2m$k $((k+4)) > "$M/fold_m$k.log" 2>&1 &
done
echo "[bpm2] waiting for 8 folds..."
until ( for t in bpm2h0 bpm2h1 bpm2h2 bpm2h3 bpm2m0 bpm2m1 bpm2m2 bpm2m3; do \
          grep -q "BOLTZ_${t}_DONE" "$M/fold_${t:4:1}${t:5}.log" 2>/dev/null || exit 1; done ); do sleep 20; done
echo "[bpm2 $(date -u +%FT%TZ)] folds done; scoring"
PYS="$R/third_party/LigandMPNN/.venv/bin/python"
"$PYS" - <<'PY'
import os, subprocess, glob, csv, json
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"; M=f"{R}/runs/bpair_mature2"
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
HS=["bpm2h0","bpm2h1","bpm2h2","bpm2h3"]; MS=["bpm2m0","bpm2m1","bpm2m2","bpm2m3"]
res=[]
for r in rows:
    d=r["name"]; h=ipsae(d,HS); m=ipsae(d,MS)
    if h is None or m is None: continue
    res.append((d,h,m,min(h,m),abs(h-m)))
wt={d[:-3]:(mn,dl) for (d,h,m,mn,dl) in res if d.endswith("_wt")}
out=[]
print(f"{'variant':18}{'h':>7}{'m':>7}{'min':>7}{'dHM':>7}{'wt_min':>8}  DUAL-STRONG")
for d,h,m,mn,dl in sorted(res,key=lambda x:-x[3]):
    base=None
    for pfx,(wmn,wdl) in wt.items():
        if d.startswith(pfx+"_H") or d==pfx+"_wt": base=wmn; break
    is_his=("_H" in d) and not d.endswith("_wt")
    # DUAL-STRONG cross-react: balanced (dHM<=0.13) AND strong (min>=0.55) AND keeps binding vs wt
    dual = is_his and base is not None and mn>=0.55 and dl<=0.13 and (mn-base)>=-0.08
    out.append(dict(variant=d,h=h,m=m,min=mn,dHM=dl,wt_min=base,dual_strong=bool(dual)))
    bs=f"{base:.3f}" if base is not None else "  -  "
    print(f"{d:18}{h:>7.3f}{m:>7.3f}{mn:>7.3f}{dl:>7.3f}{bs:>8}  {'YES' if dual else ''}")
json.dump(out, open(f"{M}/scored.json","w"), indent=1)
cands=[o for o in out if o["dual_strong"]]
print(f"\nDUAL-STRONG candidates (balanced + strong, keep binding): {len(cands)}")
for o in sorted(cands,key=lambda x:-x['min'])[:20]: print(f"  {o['variant']:18} min={o['min']:.3f} dHM={o['dHM']:.3f}")
PY
echo "BPM2_SCORE_DONE $(date -u +%FT%TZ)"
