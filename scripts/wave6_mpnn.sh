#!/bin/bash
# WAVE-6 (Campaign H): LigandMPNN redesigns the shell around FIXED 3-His separated-cation networks
# (runs/wave6/jobs.tsv) -> extract binder seqs -> variants.csv + shards. GPUs 3,4,5 ONLY (0,1,2,6,7 foreign).
# 0-Cys enforced (--omit_AA C). Uses the proven vanilla proteinmpnn ckpt (gave robust c44/b14).
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/wave6"; LM="$R/third_party/LigandMPNN"; CKPT="$LM/model_params/proteinmpnn_v_48_020.pt"
GPUS=(3 4 5)
mkdir -p "$M/mpnn"
echo "[wv6 $(date -u +%FT%TZ)] MPNN shell-adapt around fixed 3-His networks"
pids=(); gi=0
while IFS=$'\t' read -r scaffold pdb redes; do
  [ -z "$scaffold" ] && continue
  g=${GPUS[$gi]}
  CUDA_VISIBLE_DEVICES=$g nohup "$LM/.venv/bin/python" "$LM/run.py" --model_type protein_mpnn \
    --checkpoint_protein_mpnn "$CKPT" --pdb_path "$pdb" --redesigned_residues "$redes" \
    --omit_AA C --temperature 0.2 --number_of_batches 2 --batch_size 4 \
    --out_folder "$M/mpnn/$scaffold" --seed 17 \
    > "$M/mpnn_$scaffold.log" 2>&1 &
  pids+=($!); gi=$(( (gi+1)%${#GPUS[@]} ))
done < "$M/jobs.tsv"
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[wv6] MPNN done; extracting binder seqs"
"$LM/.venv/bin/python" -I - <<'PY'
import glob, csv, os
M="/home/prachh/v/anthropic-adaptyv-2026-tnf/runs/wave6"
rows=[]
for fa in sorted(glob.glob(f"{M}/mpnn/*/seqs/*.fa")):
    scaf=os.path.basename(fa)[:-3]
    lines=[l for l in open(fa).read().splitlines() if l.strip()]
    i=0
    for k in range(0,len(lines)-1,2):
        hdr=lines[k]
        if not hdr.startswith(">") or "id=" not in hdr: continue
        seq=lines[k+1].split(":")[0].strip()
        if "C" in seq: continue  # hard 0-Cys guard
        i+=1; rows.append((f"{scaf}_w6_{i}", seq))
seen=set(); ur=[]
for n,s in rows:
    if n in seen or not s: continue
    seen.add(n); ur.append((n,s))
with open(f"{M}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in ur]
NSH=3
sh=[[] for _ in range(NSH)]
for idx,r in enumerate(ur): sh[idx%NSH].append(r)
for kk,s in enumerate(sh):
    with open(f"{M}/shard{kk}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in s]
print(f"[wv6] {len(ur)} variants -> variants.csv + {NSH} shards")
PY
echo "WV6_MPNN_DONE $(date -u +%FT%TZ)"
