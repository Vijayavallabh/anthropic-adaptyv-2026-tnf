#!/bin/bash
# v15 stage-3b: place SEPARATED multi-His (Baker SOTA) on cross-react survivors, Boltz-fold each variant
# vs HUMAN (linkage) + MOUSE (cross-react), score Wyman linkage (mh_screen). Survivors that yield a
# strong multi-His switch (linkage<=-0.7, nHis>=2) + mouse-binding advance to expression + MD.
# Usage: v15_multihis.sh [CAP=24] [gpus="1 3 4 5"]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
CAP="${1:-24}"; GPUS="${2:-1 3 4 5}"
BASE="${MHBASE:-$REPO/runs/partial_v15}"; case "$BASE" in /*) ;; *) BASE="$REPO/$BASE";; esac
XR="$BASE/xreact"; MHOUT="$BASE/multihis"; rm -rf "$MHOUT"; mkdir -p "$MHOUT"
POSEGLOB="${MHPOSEGLOB:-runs/partial_v15/af2ig/gpu*/pdbs}"
PY="$REPO/envs/uplifting-esmc/bin/python"; MPNNPY="$REPO/third_party/LigandMPNN/.venv/bin/python"
MOUSE_T="$REPO/runs/buried2/mouse_tnf_mature.txt"
# --- 1. subselect survivors (balanced per lineage, cap) + build MH_SEQS ---
"$PY" - "$XR/xreact_survivors.txt" "$XR/binders.csv" "$MHOUT/mh_seqs.csv" "$MHOUT/selected.txt" "$CAP" <<'PY'
import sys,csv,re
survf,binders,seqcsv,self,cap=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],int(sys.argv[5])
surv=[l.strip() for l in open(survf) if l.strip()]
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open(binders))}
bylin={}
for n in surv:
    lin=re.sub(r'_T(10|20)_pd.*','',n); bylin.setdefault(lin,[]).append(n)
# round-robin across lineages up to cap (balanced diversity)
sel=[]; i=0
while len(sel)<cap and any(bylin.values()):
    for lin in list(bylin):
        if bylin[lin]: sel.append(bylin[lin].pop(0))
        if len(sel)>=cap: break
w=csv.writer(open(seqcsv,'w')); w.writerow(["name","sequence"])
open(self,'w').write("")
for n in sel:
    if n in seqs: w.writerow([n,seqs[n]]); open(self,'a').write(n+"\n")
print(f"survivors={len(surv)} selected={len(sel)} (cap {cap}) across lineages {list(set(re.sub(r'_T.*','',s) for s in sel))}")
PY
NSEL=$(grep -c . "$MHOUT/selected.txt" 2>/dev/null || echo 0)
[ "$NSEL" -lt 1 ] && { echo "NO SURVIVORS to His-scan — stop"; exit 0; }
# --- 2. build separated multi-His variants (reads af2pred poses + seqs) ---
echo "[v15mh] building multi-His variants for $NSEL survivors $(date -u +%FT%TZ)"
MH_SEQS="$MHOUT/mh_seqs.csv" MH_POSEGLOB="$POSEGLOB" MHOUT="$MHOUT" \
  "$MPNNPY" scripts/bpair_multihis_build.py $(cat "$MHOUT/selected.txt") 2>&1 | tail -6
NV=$(($(wc -l < "$MHOUT/variants.csv" 2>/dev/null)-1)); echo "[v15mh] $NV variants to fold"
[ "$NV" -lt 1 ] && { echo "NO VARIANTS built — stop"; exit 0; }
# --- 3. Boltz fold every variant vs human + mouse (sharded) ---
"$PY" - "$MHOUT/variants.csv" "$MHOUT" "$GPUS" <<'PY'
import sys,csv
v,out,gpus=sys.argv[1],sys.argv[2],sys.argv[3].split()
rows=list(csv.DictReader(open(v))); w={}
for g in gpus:
    for sp in ("h","m"):
        c=csv.writer(open(f"{out}/csv_{sp}_{g}.csv","w")); c.writerow(["name","sequence"]); w[(g,sp)]=c
for i,r in enumerate(rows):
    g=gpus[i%len(gpus)]
    for sp in ("h","m"): w[(g,sp)].writerow([r["name"],r["sequence"]])
print("sharded",len(rows),"variants across",len(gpus),"gpus")
PY
pids=()
for g in $GPUS; do
  [ -s "$MHOUT/csv_h_$g.csv" ] || continue
  ( bash scripts/run_boltz_tnf.sh "$MHOUT/csv_h_$g.csv" "$XR/human_tnf.txt" "v15mh_h_$g" "$g" > "$MHOUT/boltz_h_$g.log" 2>&1
    bash scripts/run_boltz_tnf.sh "$MHOUT/csv_m_$g.csv" "$MOUSE_T"          "v15mh_m_$g" "$g" > "$MHOUT/boltz_m_$g.log" 2>&1 ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"||true; done
echo "[v15mh] boltz done $(date -u +%FT%TZ)"
# --- 4. Wyman linkage on human folds (mh_screen) ---
"$PY" scripts/mh_screen.py "$MHOUT/manifest.tsv" "runs/boltz_v15mh_h_*" "$MHOUT/linkage.tsv" 2>&1 | tail -20
# --- 5. join linkage + mouse cross-react iptm ---
"$PY" - "$MHOUT/linkage.tsv" "$MHOUT" "$GPUS" <<'PY'
import sys,glob,json,csv
link,out,gpus=sys.argv[1],sys.argv[2],sys.argv[3].split()
def miptm(name):
    for g in gpus:
        c=glob.glob(f"/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/boltz_v15mh_m_{g}/out/by_seed/{name}/s0/confidence_{name}_model_0.json")
        if c: return round(json.load(open(c[0])).get("iptm",0),3)
    return None
rows=[r for r in csv.DictReader(open(link),delimiter="\t")]
res=[]
for r in rows:
    mi=miptm(r["name"]); hi=float(r["iptm_h"] or 0); lk=float(r["linkage"] or 0); nh=int(r["nHis"] or 0)
    xr=min(hi,mi) if mi is not None else None
    res.append((r["name"],r["design"],nh,lk,hi,mi,xr))
res.sort(key=lambda x:x[3])  # most-negative linkage first
print(f"{'variant':34}{'nHis':>5}{'linkage':>9}{'hum_ipt':>8}{'mo_ipt':>8}{'xreact':>8}  STRONG")
strong=[]
for n,d,nh,lk,hi,mi,xr in res:
    s = "YES" if (lk<=-0.7 and nh>=2 and (xr or 0)>=0.70) else ""
    if s: strong.append((n,d,nh,lk,xr))
    print(f"{n[:34]:34}{nh:>5}{lk:>9.3f}{hi:>8.3f}{(mi if mi is not None else -1):>8.3f}{(xr if xr is not None else -1):>8.3f}  {s}")
open(f"{out}/switch_strong.txt","w").write("\n".join(f"{n}\t{d}\t{nh}\t{lk}\t{xr}" for n,d,nh,lk,xr in strong)+"\n")
print(f"\nSTRONG multi-His switch + cross-react: {len(strong)} / {len(res)}  -> {out}/switch_strong.txt")
PY
echo "V15_MULTIHIS_DONE $(date -u +%FT%TZ)"
