#!/bin/bash
# WAVE-10 confirm: for each winner variant (args = names), rigorously verify the switch is FOLD-ROBUST and
# the binder is CROSS-REACTIVE, apples-to-apples with the portfolio:
#   - Chai-1 fold vs human  (independent 2nd folder) + Boltz already has human
#   - PROPKA+xTB triple-confirm on BOTH the Boltz AND Chai human complexes -> fold-robust His = intersection
#   - Boltz fold vs mouse -> cross-react iptm (binds both species?)
# Only a design with >=2 fold-robust His (both folds) + iptm>=0.5 human AND mouse + novel fold is a 2nd-flagship
# candidate. Usage: wave10_confirm.sh <name1> <name2> ...   (GPUs: Chai 3-5, Boltz mouse 0-2)
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
PY="$R/envs/uplifting-esmc/bin/python"
W="$R/runs/wave10"; mkdir -p "$W/confirm"
names=("$@"); [ ${#names[@]} -eq 0 ] && { echo "no winners passed"; exit 0; }
# build a winners CSV from variants.csv
"$PY" - "${names[@]}" <<'PY'
import sys,csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
want=set(sys.argv[1:])
rows=[r for r in csv.DictReader(open(f"{R}/runs/wave10/variants.csv")) if r["name"] in want]
with open(f"{R}/runs/wave10/winners.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow([r["name"],r["sequence"]]) for r in rows]
print("winners.csv:",len(rows))
PY
echo "[w10c] Chai human (gpu3-5) + Boltz mouse (gpu0) ..."
nohup bash scripts/run_chai.sh "$W/winners.csv" w10_ch data/targets/human_tnf.seq "3-5" > "$W/chai.log" 2>&1 &
nohup bash scripts/run_boltz_tnf.sh "$W/winners.csv" runs/buried2/mouse_tnf_mature.txt w10m 0 > "$W/mouse.log" 2>&1 &
wait
echo "[w10c] folds done; triple-confirm both folds"
for n in "${names[@]}"; do
  s=$(awk -F, -v n="$n" '$1==n{print $2}' "$W/winners.csv")
  BZ=$(ls runs/boltz_w10h*/out/by_seed/$n/s0/${n}_model_0.cif 2>/dev/null | head -1)
  CH=$(ls runs/chai_w10_ch/out/fast/$n/seed_0/pred.model_idx_0.cif 2>/dev/null | head -1)
  [ -n "$BZ" ] && OMP_NUM_THREADS=8 "$PY" scripts/triple_confirm_mut.py "${n}_bz" "$BZ" "$s" > "$W/confirm/${n}_bz.log" 2>&1
  [ -n "$CH" ] && OMP_NUM_THREADS=8 "$PY" scripts/triple_confirm_mut.py "${n}_ch" "$CH" "$s" > "$W/confirm/${n}_ch.log" 2>&1
  echo "  $n: boltz=$(ls runs/hhead/${n}_bz_triple.json 2>/dev/null && echo ok) chai=$(ls runs/hhead/${n}_ch_triple.json 2>/dev/null && echo ok)"
done
echo "WAVE10_CONFIRM_DONE $(date -u +%FT%TZ)"
