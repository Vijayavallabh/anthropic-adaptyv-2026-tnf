#!/bin/bash
# WAVE-4B: His-biased DEEPER interface-adapt MPNN on the 5 CONFIRMED-lineage scaffolds (reuse wave3b prepped
# PDBs+shells). Goal: stronger switches (more His-cation contacts = bigger fold-change, Baker) + novelty-diverse
# ALTERNATES per lineage (hedge the des_4098_1 novelty risk). 48 samples/scaffold, bias_AA H:1.5.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/wave4b"; LM="$R/third_party/LigandMPNN"; CKPT="$LM/model_params/proteinmpnn_v_48_020.pt"
mkdir -p "$M/mpnn"
grep -P "^(des_4098_1|des_3050_6|des_180_6|des_818_4|des_875_5)\t" "$R/runs/wave3b/jobs.tsv" > "$M/jobs.tsv"
[ -s "$M/jobs.tsv" ] || { echo "ERROR: jobs.tsv empty (runs/wave3b/jobs.tsv missing or scaffold names changed)"; exit 1; }
# GPU idle guard (CLAUDE.md: never clobber a foreign job) — abort if any GPU is busy
busy=$(nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader,nounits | awk -F', ' '$2>10||$3>2000{print $1}' | tr '\n' ' ')
[ -n "$busy" ] && { echo "ERROR: GPU(s) busy [$busy]; refusing to clobber a foreign job"; exit 1; }
echo "[wv4b $(date -u +%FT%TZ)] His-biased MPNN on $(wc -l <"$M/jobs.tsv") confirmed-lineage scaffolds"
pids=(); g=0
while IFS=$'\t' read -r scaffold pdb redes; do
  [ -z "$scaffold" ] && continue
  CUDA_VISIBLE_DEVICES=$g nohup "$LM/.venv/bin/python" "$LM/run.py" --model_type protein_mpnn \
    --checkpoint_protein_mpnn "$CKPT" --pdb_path "$pdb" --redesigned_residues "$redes" \
    --bias_AA "H:1.5" --temperature 0.25 --number_of_batches 12 --batch_size 4 --out_folder "$M/mpnn/$scaffold" --seed 29 \
    > "$M/mpnn_$scaffold.log" 2>&1 &
  pids+=($!); g=$(( (g+1)%8 ))
done < "$M/jobs.tsv"
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[wv4b] MPNN done; extracting binder seqs"
"$R/.venv/bin/python" -I - <<'PY'
import glob, csv, os
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; M=f"{R}/runs/wave4b"
wt={r["name"]:r["sequence"] for r in csv.DictReader(open(f"{R}/runs/bpair/xreact_full2.csv"))}
rows=[]
for fa in sorted(glob.glob(f"{M}/mpnn/*/seqs/*.fa")):
    scaf=os.path.basename(fa)[:-3]
    if scaf in wt: rows.append((f"{scaf}_wt", wt[scaf]))
    lines=[l for l in open(fa).read().splitlines() if l.strip()]; i=0
    for k in range(0,len(lines)-1,2):
        if not lines[k].startswith(">") or "id=" not in lines[k]: continue
        i+=1; rows.append((f"{scaf}_c{i}", lines[k+1].split(":")[0].strip()))
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
print(f"[wv4b] {len(ur)} variants -> variants.csv + shards")
PY
echo "[wv4b] folding + analyzing + triple-confirming"
bash scripts/wave3_fold.sh "$M" wv4b
"$R/.venv/bin/python" -I scripts/wave4b_analyze.py
PYTHONUNBUFFERED=1 "$R/.venv/bin/python" -I scripts/wave4b_triple.py
echo "WV4B_ALL_DONE $(date -u +%FT%TZ)"
