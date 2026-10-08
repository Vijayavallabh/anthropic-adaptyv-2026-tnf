#!/bin/bash
# For every candidate complex in runs/pool/cands: emit His-switch variants (scan_switch.py),
# af2ig all variants (binding retained?), propka-score release linkage, write a ranked table.
set -uo pipefail; cd /home/prachh/v/anthropic-adaptyv-2026-tnf
SCAN=runs/pool/switch_scan; AIG=runs/pool/switch_af2ig; rm -rf "$SCAN" "$AIG"; mkdir -p "$SCAN"
for c in runs/pool/cands/*.pdb; do
  uv run python scripts/scan_switch.py "$c" "$SCAN" 6 2>&1 | grep -v VIRTUAL_ENV | head -1
done
echo "[bsw] His variants: $(ls "$SCAN"/*.pdb 2>/dev/null|wc -l)"
uv run python scripts/run_af2ig.py --pdbdir "$SCAN" --out "$AIG" --gpus 0-7 --mode fast > runs/pool/switch_af2ig_run.log 2>&1
echo "[bsw] scored $(( $(wc -l < "$AIG/out.sc" 2>/dev/null||echo 1) -1 )) variants; propka-scoring release linkage..."
# propka release score for each folded variant with binding retained (pae<14)
: > runs/pool/switch_table.tsv
echo -e "pae_switch\tplddt\trmsd\tddg_pH\tvariant" >> runs/pool/switch_table.tsv
while read -r line; do
  case "$line" in SCORE:*) ;; *) continue;; esac
  set -- $line; rmsd=$2; pae=$4; plb=$6; desc=${!#}
  case "$pae$plb$rmsd" in *nan*) continue;; esac   # any nan field -> skip (don't let it mis-gate)
  awk "BEGIN{exit !($pae<14 && $plb>80 && $rmsd<2.5)}" || continue
  v=$(find "$AIG" -name "${desc}.pdb" | head -1); [ -n "$v" ] || continue
  dg=$(uv run python scripts/ph_release_score.py "$v" 2>/dev/null | grep 'ddg_pH (sum' | grep -oE '[-+][0-9.]+')
  echo -e "${pae}\t${plb}\t${rmsd}\t${dg}\t${desc}" >> runs/pool/switch_table.tsv
done < "$AIG/out.sc"
echo "=== SWITCH TABLE (binding retained pae<14, sorted by release ddg_pH) ==="
(head -1 runs/pool/switch_table.tsv; tail -n +2 runs/pool/switch_table.tsv | sort -t$'\t' -k4 -g) | column -t
echo BATCH_SWITCH_DONE
