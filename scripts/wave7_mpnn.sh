#!/bin/bash
# WAVE-7: His-biased SOLUBLE MPNN on c19-partial-diffused backbones (GPUs 3,4,5). 0-Cys. Bias His near
# conserved cations (runs/wave7/bias_AA_per_residue_multi.json). Extract binder seqs -> variants.csv.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/wave7"; L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"
CKPT="$L/model_params/solublempnn_v_48_020.pt"; BIAS="$M/bias_AA_per_residue_multi.json"
GPUS=(3 4 5); NSEQ=${1:-3}
"$PY" -I - "$M/pdbs.json" "$M" 3 <<'PY'
import json,sys
pdbs=json.load(open(sys.argv[1])); work=sys.argv[2]; ng=int(sys.argv[3])
import os
for i in range(ng):
    os.makedirs(f"{work}/gpu{i}",exist_ok=True)
    json.dump(pdbs[i::ng], open(f"{work}/gpu{i}/multi.json","w"))
print("sharded",len(pdbs))
PY
pids=()
for i in 0 1 2; do g=${GPUS[$i]}
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
      --pdb_path_multi "$M/gpu$i/multi.json" --chains_to_design A \
      --bias_AA_per_residue_multi "$BIAS" --bias_AA "C:-100.0" \
      --out_folder "$M/gpu$i" --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 --seed 19 ) \
      > "$M/mpnn_gpu$i.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[wv7] MPNN done; extracting"
"$PY" -I - "$M" <<'PY'
import glob,csv,os,sys
M=sys.argv[1]; rows=[]
for fa in sorted(glob.glob(f"{M}/gpu*/seqs/*.fa")):
    scaf=os.path.basename(fa)[:-3]
    lines=[l for l in open(fa).read().splitlines() if l.strip()]
    i=0
    for k in range(0,len(lines)-1,2):
        h=lines[k]
        if not h.startswith(">") or "id=" not in h: continue
        seq=lines[k+1].split(":")[0].strip()
        if "C" in seq: continue
        i+=1; rows.append((f"{scaf}_s{i}",seq))
seen=set(); ur=[]
for n,s in rows:
    if n in seen or not s: continue
    seen.add(n); ur.append((n,s))
with open(f"{M}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in ur]
NSH=3; sh=[[] for _ in range(NSH)]
for idx,r in enumerate(ur): sh[idx%NSH].append(r)
for kk,s in enumerate(sh):
    with open(f"{M}/shard{kk}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in s]
print(f"[wv7] {len(ur)} variants -> variants.csv + {NSH} shards")
PY
echo "WV7_MPNN_DONE"
