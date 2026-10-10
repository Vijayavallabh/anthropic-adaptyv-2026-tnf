#!/bin/bash
# Generalized His-position scan on a cross-react winner: install His at each binder interface position
# within 9A of a target cation, fold human+mouse (Boltz fast), score ipSAE_min. GPU-gated: waits for 2
# free GPUs before folding. Usage: xhiscan.sh <design_name_in_top_seqs> <gpu_h> <gpu_m>
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
D="$1"; GH="${2:-}"; GM="${3:-}"
OUT="runs/xhiscan/$D"; mkdir -p "$OUT"
# build His-variant CSV from the design's af2pred complex
.venv/bin/python -I - "$D" "$OUT" <<'PY'
import sys,glob,numpy as np,csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; D=sys.argv[1]; OUT=sys.argv[2]
seqs={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/partial4/top_seqs.csv"))}
s=seqs[D]
# af2pred complex matching this exact design name
cands=glob.glob(f"{R}/runs/partial4/af2ig/gpu*/pdbs/{D}_af2pred.pdb")
f=cands[0]
bCB={}; tcat=[]
for l in open(f):
    if l[:4]!="ATOM": continue
    ch=l[21]; res=l[17:20].strip(); at=l[12:16].strip(); rn=int(l[22:26])
    c=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    if ch=="A" and at=="CB": bCB[rn]=(c,res)
    if ch=="B" and ((res=="LYS" and at=="NZ")or(res=="ARG" and at=="CZ")): tcat.append(c)
pos=sorted(rn for rn,(c,res) in bCB.items() if res not in("GLY","PRO","HIS") and any(np.linalg.norm(c-tc)<9 for tc in tcat))
tag=D.split("_des_")[-1]
rows=[(f"{tag}_wt",s)]+[(f"{tag}_H{p}", s[:p-1]+"H"+s[p:]) for p in pos if p<=len(s)]
with open(f"{OUT}/scan.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in rows]
print(f"{D}: af2pred={f.split('/')[-1]} positions={pos} -> {len(rows)} seqs")
PY
# GPU-gate: wait for 2 idle GPUs if not specified
pick_free(){ nvidia-smi --query-gpu=index,memory.used --format=csv,noheader | awk -F, '$2+0<500{print $1}'; }
if [ -z "$GH" ]; then
  echo "waiting for 2 free GPUs..."
  until [ "$(pick_free | wc -l)" -ge 2 ]; do sleep 60; done
  mapfile -t FG < <(pick_free); GH=${FG[0]}; GM=${FG[1]}
fi
echo "folding $D human=gpu$GH mouse=gpu$GM"
TAG=$(basename "$D" | sed 's/.*_des_/xh/')
nohup bash scripts/run_boltz_tnf.sh "$OUT/scan.csv" runs/buried2/human_tnf.txt ${TAG}h "$GH" > "$OUT/fold_h.log" 2>&1 &
nohup bash scripts/run_boltz_tnf.sh "$OUT/scan.csv" runs/buried2/mouse_tnf_mature.txt ${TAG}m "$GM" > "$OUT/fold_m.log" 2>&1 &
wait
echo "XHISCAN_FOLD_DONE $D $(date -u +%FT%TZ)"
