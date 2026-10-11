#!/bin/bash
# Post-MD harvest for v18: from a two-state-MD dir, take CONFIRMED switchers and fully re-validate each
# (cross-react Boltz human+mouse -> min-iptm; solubility; seq + structural novelty proxy) so they can be
# swapped into v16/v17 NONRESPOND slots (<=2 per lineage). Single-protomer bpair track (target 1TNF_B).
# Usage: harvest_v18.sh <md_verdict_dir> <variants.csv name,seq> <tag> "<gpus 'a b ..'>"
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
MDDIR=$1; VARCSV=$2; TAG=$3; IFS=' ' read -ra GA <<< "${4:-0 1 2 3 4 5 6 7}"
OUT="$REPO/runs/${TAG}"; mkdir -p "$OUT"
PY="$REPO/envs/uplifting-esmc/bin/python"; MDPY="$REPO/envs/md2/bin/python"
HUM="$REPO/runs/partial_denovo2/funnel/xreact/human_tnf.txt"; MOU="$REPO/runs/buried2/mouse_tnf_mature.txt"
# 1) CONFIRMED switchers from MD
"$MDPY" scripts/md_verdict.py "$MDDIR" 2>&1 | tee "$OUT/md_verdict.txt"
awk '$NF=="CONFIRMED"{print $1}' "$OUT/md_verdict.txt" > "$OUT/confirmed.txt"
NC=$(grep -c . "$OUT/confirmed.txt" 2>/dev/null || echo 0)
echo "[harvest] CONFIRMED switchers: $NC"
[ "$NC" -lt 1 ] && { echo "no confirmed switchers — nothing to harvest"; exit 0; }
# 2) build seqs CSV for confirmed
"$PY" - "$VARCSV" "$OUT/confirmed.txt" "$OUT/confirmed_seqs.csv" <<'PY'
import sys,csv
var,conf,out=sys.argv[1:4]
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open(var))}
names=[l.strip() for l in open(conf) if l.strip()]
w=csv.writer(open(out,"w")); w.writerow(["name","sequence"])
for n in names:
    if n in seqs: w.writerow([n,seqs[n]])
print("wrote",sum(1 for n in names if n in seqs),"confirmed seqs")
PY
# 3) cross-react Boltz human + mouse (sharded), min-iptm
"$PY" - "$OUT/confirmed_seqs.csv" "$OUT" "${GA[*]}" <<'PY'
import sys,csv
csvf,out,gpus=sys.argv[1],sys.argv[2],sys.argv[3].split()
rows=list(csv.DictReader(open(csvf))); w={}
for g in gpus:
    for sp in ("h","m"):
        c=csv.writer(open(f"{out}/xr_{sp}_{g}.csv","w")); c.writerow(["name","sequence"]); w[(g,sp)]=c
for i,r in enumerate(rows):
    g=gpus[i%len(gpus)]
    for sp in ("h","m"): w[(g,sp)].writerow([r["name"],r["sequence"]])
print("sharded",len(rows),"across",len(gpus),"gpus")
PY
pids=()
for g in "${GA[@]}"; do
  [ -s "$OUT/xr_h_$g.csv" ] || continue
  ( bash scripts/run_boltz_tnf.sh "$OUT/xr_h_$g.csv" "$HUM" "v18xr_h_$g" "$g" > "$OUT/boltz_h_$g.log" 2>&1
    bash scripts/run_boltz_tnf.sh "$OUT/xr_m_$g.csv" "$MOU" "v18xr_m_$g" "$g" > "$OUT/boltz_m_$g.log" 2>&1 ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"||true; done
# 4) solubility
uv run python scripts/solubility_predict.py "$OUT/confirmed_seqs.csv" "$OUT/solub.csv" 2>&1 | tail -2 || true
# 5) collect: dSwitch + xr(min iptm) + solub -> ranked table
"$PY" - "$OUT" "${GA[*]}" <<'PY'
import sys,glob,json,csv,os
out,gpus=sys.argv[1],sys.argv[2].split()
def iptm(tag,n):
    for g in gpus:
        c=glob.glob(f"/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/boltz_{tag}_{g}/out/**/confidence_{n}_model_0.json",recursive=True)
        if c: return round(json.load(open(c[0])).get("iptm",0),3)
    return None
# dSwitch from md_verdict.txt
dsw={}
for ln in open(f"{out}/md_verdict.txt"):
    f=ln.split()
    if len(f)>8 and f[-1] in ("CONFIRMED","NONRESPOND","MARGINAL"):
        try: dsw[f[0]]=float(f[5])
        except: pass
solub={}
if os.path.exists(f"{out}/solub.csv"):
    for r in csv.DictReader(open(f"{out}/solub.csv")):
        try: solub[r["name"]]=float(r.get("solubility") or r.get("scaled_solubility") or 0)
        except: pass
names=[r["name"] for r in csv.DictReader(open(f"{out}/confirmed_seqs.csv"))]
rows=[]
for n in names:
    h=iptm("v18xr_h",n); m=iptm("v18xr_m",n); xr=min(h,m) if (h is not None and m is not None) else None
    rows.append((n,dsw.get(n),h,m,xr,solub.get(n)))
rows.sort(key=lambda r:-(r[1] or -9))
print(f"{'confirmed switcher':40}{'dSwitch':>8}{'ipt_h':>7}{'ipt_m':>7}{'xr':>6}{'solub':>7}")
for n,ds,h,m,xr,s in rows:
    f=lambda v,p=3:(f'{v:.{p}f}' if isinstance(v,(int,float)) else 'NA')
    print(f"{n[:40]:40}{f(ds,2):>8}{f(h):>7}{f(m):>7}{f(xr):>6}{f(s,2):>7}")
csv.writer(open(f"{out}/v18_candidates.csv","w")).writerows(
    [["name","dSwitch","iptm_h","iptm_m","xr_min","solub"]]+[[n,ds,h,m,xr,s] for n,ds,h,m,xr,s in rows])
print(f"\n-> {out}/v18_candidates.csv  (then: swap xr>=0.70 + solub>=0.5 into NONRESPOND slots, <=2/lineage, novelty-check, validate)")
PY
echo "HARVEST_V18_DONE $(date -u +%FT%TZ)"
