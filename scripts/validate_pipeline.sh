#!/bin/bash
# Wait for the current RFd gen sweep to finish, then run MPNN+af2ig on the batch and
# report raw-contact distribution + af2ig pass rate (pae_interaction<10). One-shot validation.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
echo "[val] waiting for gen sweep to finish $(date -u +%T)"
for i in $(seq 1 180); do
  grep -q 'sweep done' runs/gen_test_noise0.log 2>/dev/null && { echo "[val] gen done"; break; }
  sleep 20
done
N=$(ls runs/rfd_gen/g*/des_*.pdb 2>/dev/null | wc -l)
echo "[val] $N RFd designs on disk"
echo "=== RAW CONTACTS (noise=0, aromatic hotspots) ==="
uv run python scripts/rawcontacts2.py "runs/rfd_gen/g*/des_*.pdb" 2>&1 | grep -v VIRTUAL_ENV
echo "=== MPNN + af2ig (gpus 0-3, 4 seq/bb) $(date -u +%T) ==="
bash scripts/seq_af2ig_tnf.sh 0-3 4 'g*' 2>&1 | grep -viE 'VIRTUAL_ENV|oneDNN|LEVER|deprecat|map_fn|fn_output'
echo "=== AF2IG PASS RATE ==="
SC=runs/seq/af2ig/out.sc
if [ -s "$SC" ]; then
  awk 'NR>1 && $4!="nan"{n++; if($4+0<10 && $6+0>80 && $2+0<2)p++; if($4+0<10)pl++} END{printf "scored=%d  pass(pae<10 & plddt_b>80 & rmsd<2)=%d  pae<10(any)=%d\n", n, p+0, pl+0}' "$SC"
  echo "--- top 12 by pae_interaction ---"
  (head -1 "$SC"; awk 'NR>1 && $4!="nan"' "$SC" | sort -k4 -n | head -12)
else echo "NO af2ig scores"; fi
echo "VALIDATE_DONE $(date -u +%FT%TZ)"
