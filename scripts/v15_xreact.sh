#!/bin/bash
# v15 stage-3a: from stage-2 af2ig binders -> best seq per backbone -> top-N by binding ->
# Boltz-2 co-fold vs HUMAN + MOUSE TNF -> cross-react min(iptm) gate (>=0.70). The real TNF mouse
# target is runs/buried2/mouse_tnf_mature.txt (NOT the EGFR make_mouse_target.py).
# Usage: v15_xreact.sh [TOPN=120] [gpus space-list="1 3 4 5 6 7"]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
TOPN="${1:-120}"; GPUS="${2:-1 3 4 5 6 7}"
BASE="${XRBASE:-$REPO/runs/partial_v15}"; case "$BASE" in /*) ;; *) BASE="$REPO/$BASE";; esac
SC="$BASE/af2ig/out.sc"; AF2IN="$BASE/af2in"
OUT="$BASE/xreact"; mkdir -p "$OUT"
PY="$REPO/envs/uplifting-esmc/bin/python"
HUMAN_B="$REPO/data/targets/1TNF_B.pdb"; MOUSE_T="$REPO/runs/buried2/mouse_tnf_mature.txt"
# --- 1. select passing binders, best per backbone, top-N by pae_interaction; extract chain-A seq ---
"$PY" - "$SC" "$AF2IN" "$OUT/binders.csv" "$TOPN" <<'PY'
import sys,glob,os,re
sc,af2in,outcsv,topn=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4])
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def chainA_seq(pdb):
    seq={}
    for l in open(pdb):
        if l.startswith("ATOM") and l[21]=='A' and l[12:16].strip()=='CA':
            seq[int(l[22:26])]=AA3.get(l[17:20].strip(),'X')
    return ''.join(seq[k] for k in sorted(seq))
# parse out.sc: cols SCORE: rmsd pae_binder pae_int pae_target plddt_binder ... description(last)
rows=[]
for ln in open(sc):
    if not ln.startswith("SCORE:") or "description" in ln: continue
    f=ln.split()
    try: rmsd,pae_int,plddt_b,desc=float(f[1]),float(f[3]),float(f[5]),f[-1]
    except: continue
    if pae_int<10 and plddt_b>80 and rmsd<2: rows.append((pae_int,desc,rmsd,plddt_b))
# map description -> af2in stem (stem is the longest af2in basename that is a prefix of desc)
stems=sorted({os.path.basename(p)[:-4] for p in glob.glob(af2in+"/*.pdb")}, key=len, reverse=True)
def to_stem(desc):
    for s in stems:
        if desc.startswith(s): return s
    return None
best={}  # af2in stem -> (pae_int, rmsd, plddt_b)
for pae_int,desc,rmsd,plddt_b in rows:
    s=to_stem(desc)
    if s is None: continue
    if s not in best or pae_int<best[s][0]: best[s]=(pae_int,rmsd,plddt_b)
# backbone = stem minus trailing _<seqnum>; keep best stem per backbone
bb={}
for s,(pi,rm,pl) in best.items():
    b=re.sub(r'_\d+$','',s)
    if b not in bb or pi<bb[b][1]: bb[b]=(s,pi,rm,pl)
percap=int(os.environ.get("PERCAP","50"))
bylin={}
for v in sorted(bb.values(), key=lambda x:x[1]):
    lin=re.sub(r'_T(10|20)_pd.*','',v[0]); bylin.setdefault(lin,[]).append(v)
ranked=[]
for lin,items in bylin.items(): ranked+=items[:percap]
ranked=sorted(ranked,key=lambda x:x[1])[:topn]
print("per-lineage kept:",{k:len(v[:percap]) for k,v in bylin.items()})
import csv as _csv
w=_csv.writer(open(outcsv,'w')); w.writerow(["name","sequence","pae_int","rmsd","plddt_b"])
n=0
for s,pi,rm,pl in ranked:
    seq=chainA_seq(os.path.join(af2in,s+".pdb"))
    if len(seq)<60: continue
    w.writerow([s,seq,round(pi,2),round(rm,2),round(pl,1)]); n+=1
print(f"passing-stems={len(best)} backbones={len(bb)} -> wrote {n} (top {topn}) to {outcsv}")
PY
NB=$(($(wc -l < "$OUT/binders.csv")-1)); echo "[v15xr] $NB binder backbones to cross-fold $(date -u +%FT%TZ)"
[ "$NB" -lt 1 ] && { echo "NO BINDERS PASSED STAGE-2 — stop"; exit 0; }
# --- 2. shard binders.csv across GPUs, co-fold vs human + mouse ---
# build human target seq file from 1TNF_B.pdb chain B
"$PY" - "$HUMAN_B" "$OUT/human_tnf.txt" <<'PY'
import sys
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
pdb,out=sys.argv[1],sys.argv[2]; seq={}
for l in open(pdb):
    if l.startswith("ATOM") and l[21]=='B' and l[12:16].strip()=='CA': seq[int(l[22:26])]=AA3.get(l[17:20].strip(),'X')
open(out,'w').write(''.join(seq[k] for k in sorted(seq)))
print("human target len",len(seq))
PY
# shard binders.csv into per-GPU human+mouse CSVs (round-robin) — robust Python split
rm -f "$OUT"/csv_h_*.csv "$OUT"/csv_m_*.csv
"$PY" - "$OUT/binders.csv" "$OUT" "$GPUS" <<'PY'
import sys,csv
binders,out,gpus=sys.argv[1],sys.argv[2],sys.argv[3].split()
rows=list(csv.DictReader(open(binders)))
w={}
for g in gpus:
    for sp in ("h","m"):
        fh=open(f"{out}/csv_{sp}_{g}.csv","w"); c=csv.writer(fh); c.writerow(["name","sequence"]); w[(g,sp)]=c
for i,r in enumerate(rows):
    g=gpus[i%len(gpus)]
    for sp in ("h","m"): w[(g,sp)].writerow([r["name"],r["sequence"]])
print("sharded",len(rows),"binders across",len(gpus),"gpus")
PY
pids=()
for g in $GPUS; do
  [ -s "$OUT/csv_h_$g.csv" ] || continue
  ( bash scripts/run_boltz_tnf.sh "$OUT/csv_h_$g.csv" "$OUT/human_tnf.txt" "v15xr_h_$g" "$g" > "$OUT/boltz_h_$g.log" 2>&1
    bash scripts/run_boltz_tnf.sh "$OUT/csv_m_$g.csv" "$MOUSE_T"          "v15xr_m_$g" "$g" > "$OUT/boltz_m_$g.log" 2>&1 ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"||true; done
echo "[v15xr] boltz done $(date -u +%FT%TZ)"
# --- 3. collect iptm, compute cross-react ---
"$PY" - "$OUT" "$GPUS" <<'PY'
import sys,glob,json,os
OUT=sys.argv[1]; gpus=sys.argv[2].split()
def iptm(tag,name):
    g=glob.glob(f"/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/boltz_{tag}/out/by_seed/{name}/s0/confidence_{name}_model_0.json")
    if not g: g=glob.glob(f"/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/boltz_{tag}/out/**/confidence_{name}_model_0.json",recursive=True)
    return round(json.load(open(g[0])).get("iptm",0),3) if g else None
import csv
names=[r["name"] for r in csv.DictReader(open(f"{OUT}/binders.csv"))]
rows=[]
for n in names:
    h=next((iptm(f"v15xr_h_{g}",n) for g in gpus if iptm(f"v15xr_h_{g}",n) is not None), None)
    m=next((iptm(f"v15xr_m_{g}",n) for g in gpus if iptm(f"v15xr_m_{g}",n) is not None), None)
    xr=min(h,m) if (h is not None and m is not None) else None
    rows.append((n,h,m,xr))
rows.sort(key=lambda r:-(r[3] or -1))
print(f"{'design':44}{'hum_iptm':>9}{'mouse_iptm':>11}{'xreact':>8}  gate")
surv=[]
for n,h,m,xr in rows:
    g="XREACT" if (xr or 0)>=0.70 else ""
    if g: surv.append(n)
    print(f"{n[:44]:44}{(h if h is not None else -1):>9.3f}{(m if m is not None else -1):>11.3f}{(xr if xr is not None else -1):>8.3f}  {g}")
open(f"{OUT}/xreact_survivors.txt","w").write("\n".join(surv)+"\n")
print(f"\nCROSS-REACTIVE (min iptm>=0.70): {len(surv)} / {len(rows)}")
PY
echo "V15_XREACT_DONE $(date -u +%FT%TZ)"
