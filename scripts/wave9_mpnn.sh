set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
M="$R/runs/wave9"; L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"
CKPT="$L/model_params/solublempnn_v_48_020.pt"; BIAS="$M/bias.json"
GPUS=(3 4 5); gi=0; pids=()
while IFS=$'\t' read -r name pdb redes; do
  [ -z "$name" ] && continue
  g=${GPUS[$gi]}
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn \
      --checkpoint_protein_mpnn "$CKPT" --pdb_path "$pdb" --redesigned_residues "$redes" \
      --bias_AA_per_residue_multi "$BIAS" --bias_AA "C:-100.0" \
      --out_folder "$M/mpnn/$name" --batch_size 4 --number_of_batches 24 --temperature 0.25 --seed 29 ) \
      > "$M/mpnn_$name.log" 2>&1 &
  pids+=($!); gi=$(( (gi+1)%3 ))
done < "$M/jobs.tsv"
for p in "${pids[@]}"; do wait "$p" || true; done
echo "[wv9] MPNN done; extracting"
"$PY" -I - "$M" <<'PY'
import glob,csv,os,sys
M=sys.argv[1]; rows=[]
for fa in sorted(glob.glob(f"{M}/mpnn/*/seqs/*.fa")):
    scaf=os.path.basename(fa)[:-3]
    lines=[l for l in open(fa).read().splitlines() if l.strip()]
    i=0
    for k in range(0,len(lines)-1,2):
        h=lines[k]
        if not h.startswith(">") or "id=" not in h: continue
        s=lines[k+1].split(":")[0].strip()
        if "C" in s: continue
        i+=1; rows.append((f"{scaf}_w9_{i}",s))
seen=set(); ur=[]
for n,s in rows:
    if s in seen or not s: continue
    seen.add(s); ur.append((n,s))
with open(f"{M}/variants.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in ur]
sh=[[] for _ in range(3)]
for idx,r in enumerate(ur): sh[idx%3].append(r)
for kk,s in enumerate(sh):
    with open(f"{M}/shard{kk}.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in s]
print(f"[wv9] {len(ur)} unique variants -> variants.csv + 3 shards")
PY
echo WV9_MPNN_DONE
