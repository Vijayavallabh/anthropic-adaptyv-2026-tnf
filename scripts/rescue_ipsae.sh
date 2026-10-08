#!/bin/bash
# Score ipSAE on rescued Boltz folds (human + mouse); emit ipsae_h, ipsae_m, ipsae_wc per variant.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf
PY="uv run --no-sync python"
cd "$R"
declare -A H M
for sp in h m; do
  for d in runs/boltz_rescue_$sp/out/by_seed/*/s0; do
    [ -d "$d" ] || continue
    name=$(basename "$(dirname "$d")")
    cif="$d/${name}_model_0.cif"; pae="$d/pae_${name}_model_0.npz"
    [ -f "$cif" ] && [ -f "$pae" ] || { echo "MISS $name $sp"; continue; }
    $PY tools/ipsae.py "$pae" "$cif" 10 10 >/dev/null 2>&1
    txt="$d/${name}_model_0_10_10.txt"
    val=$(awk '$5=="max"{print $6; exit}' "$txt" 2>/dev/null)
    [ "$sp" = h ] && H[$name]=$val || M[$name]=$val
  done
done
printf "%-12s %8s %8s %8s\n" variant ipsae_h ipsae_m ipsae_wc
for k in "${!H[@]}"; do
  h=${H[$k]:-NA}; m=${M[$k]:-NA}
  wc=$(python3 -c "print(f'{min($h,$m):.3f}')" 2>/dev/null || echo NA)
  printf "%-12s %8s %8s %8s\n" "$k" "$h" "$m" "$wc"
done | sort
