#!/bin/bash
# Joint multi-state (human<->mouse) tied MPNN on the 3 top switch designs: redesign interface only,
# tie binder chain A<->C, keep switch-His + core fixed. One sequence optimized for BOTH orthologs.
# Vanilla ProteinMPNN (keeps interface hydrophobics). Extract chain-A binder -> jxr.csv.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
OUT="$R/runs/jxr"; W="$OUT/mpnn"; rm -rf "$W"; mkdir -p "$W"
GPUS=(5 6 7); i=0; pids=()
for d in des59_H48 b1m9_A12H b1m7_A66H b4m2_A57H; do
  g=${GPUS[$((i%3))]}; i=$((i+1))
  PDB=$("$PY" -c "import json;print(json.load(open('$OUT/meta.json'))['$d']['pdb'])")
  RED=$("$PY" -c "import json;print(json.load(open('$OUT/meta.json'))['$d']['redesigned'])")
  SYM=$("$PY" -c "import json;print(json.load(open('$OUT/meta.json'))['$d']['symmetry'])")
  WTS=$("$PY" -c "import json;print(json.load(open('$OUT/meta.json'))['$d']['weights'])")
  mkdir -p "$W/$d"
  # NOTE: --chains_to_design splits on COMMA not space; redesigned_residues already names chains -> omit it.
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn \
      --checkpoint_protein_mpnn "$CKPT" --pdb_path "$PDB" \
      --redesigned_residues "$RED" --symmetry_residues "$SYM" --symmetry_weights "$WTS" \
      --omit_AA C --out_folder "$W/$d" --batch_size 8 --number_of_batches 8 --temperature 0.3 \
  ) > "$W/$d.log" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
# extract chain-A binder sequence (field 0 before first ':'), dedupe
"$PY" - <<'PY'
import glob, csv, os
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; W=f"{R}/runs/jxr/mpnn"
rows=[]; seen=set()
for d in ("des59_H48","b1m9_A12H","b1m7_A66H","b4m2_A57H"):
    for fa in glob.glob(f"{W}/{d}/seqs/*.fa"):
        lines=[l.strip() for l in open(fa)]
        for h,s in zip(lines[::2],lines[1::2]):
            if h.startswith(">") and ", id=" in h:
                binder=s.split(":")[0]
                if binder and binder not in seen:
                    seen.add(binder); idx=len([r for r in rows if r[0].startswith("jxr_"+d)])
                    rows.append((f"jxr_{d}_{idx}", binder))
with open(f"{R}/runs/jxr/jxr.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in rows]
print("jxr seqs:", len(rows))
PY
echo "JXR_MPNN_DONE $(date -u +%FT%TZ)"
