#!/bin/bash
# WAVE-3B orchestrator: LigandMPNN (vanilla proteinmpnn ckpt) redesigns the interface shell AROUND fixed
# multi-His on each scaffold (runs/wave3b/jobs.tsv) -> extract binder seqs -> fold human+mouse -> analyze.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/wave3b"; LM="$R/third_party/LigandMPNN"; CKPT="$LM/model_params/proteinmpnn_v_48_020.pt"
mkdir -p "$M/mpnn"
echo "[wv3b $(date -u +%FT%TZ)] MPNN interface-adapt on scaffolds"
pids=(); g=0
while IFS=$'\t' read -r scaffold pdb redes; do
  [ -z "$scaffold" ] && continue
  CUDA_VISIBLE_DEVICES=$g nohup "$LM/.venv/bin/python" "$LM/run.py" --model_type protein_mpnn \
    --checkpoint_protein_mpnn "$CKPT" --pdb_path "$pdb" --redesigned_residues "$redes" \
    --temperature 0.2 --number_of_batches 4 --batch_size 4 --out_folder "$M/mpnn/$scaffold" --seed 17 \
    > "$M/mpnn_$scaffold.log" 2>&1 &
  pids+=($!); g=$(( (g+1)%8 ))
done < "$M/jobs.tsv"
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[wv3b] MPNN done; extracting binder seqs"
"$R/.venv/bin/python" -I - <<'PY'
import glob, csv, os
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; M=f"{R}/runs/wave3b"
wt={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/bpair/xreact_full2.csv"))}
rows=[]
for fa in sorted(glob.glob(f"{M}/mpnn/*/seqs/*.fa")):
    scaf=os.path.basename(fa)[:-3]
    if scaf in wt: rows.append((f"{scaf}_wt", wt[scaf]))
    lines=[l for l in open(fa).read().splitlines() if l.strip()]
    i=0
    for k in range(0,len(lines)-1,2):
        hdr=lines[k]
        if not hdr.startswith(">") or "id=" not in hdr: continue
        seq=lines[k+1].split(":")[0].strip()
        i+=1; rows.append((f"{scaf}_b{i}", seq))
seen=set(); ur=[]
for n,s in rows:
    if n in seen or not s: continue
    seen.add(n); ur.append((n,s))
with open(f"{M}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in ur]
sh=[[] for _ in range(4)]
for idx,r in enumerate(ur): sh[idx%4].append(r)
for kk,s in enumerate(sh):
    with open(f"{M}/shard{kk}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in s]
print(f"[wv3b] {len(ur)} variants ({sum(1 for n,_ in ur if '_b' in n)} designs) -> variants.csv + shards")
PY
echo "[wv3b] folding + analyzing"
bash scripts/wave3_fold.sh "$M" wv3b
"$R/.venv/bin/python" -I scripts/wave3b_analyze.py
echo "WV3B_ALL_DONE $(date -u +%FT%TZ)"
