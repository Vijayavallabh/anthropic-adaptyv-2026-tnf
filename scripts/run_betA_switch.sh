#!/bin/bash
# After AF2-ig fold-filter: collect RFd binder passers -> canonicalize (His433 -> B/100)
# -> install + score the His433 pH switch (PROPKA Wyman linkage). CPU/PyRosetta (runs alongside GPU work).
# Usage: run_betA_switch.sh <tag> [max_pae] [workers]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
export TMPDIR="$REPO/runs/tmp"; mkdir -p "$TMPDIR"   # /tmp is read-only/stale here; obj1_switch writes temp PDBs
export PATH="$REPO/envs/pyrosetta/bin:$PATH"          # propka3 lives here; obj1_switch execs it by name
TAG=${1:-w1}; MAXPAE=${2:-12}; WORKERS=${3:-48}
MD="$REPO/envs/md/bin/python"; PR="$REPO/envs/pyrosetta/bin/python"
echo "=== [$TAG] collect AF2-ig passers (pae<$MAXPAE) $(date -u +%FT%TZ) ==="
"$MD" scripts/collect_passers.py --runs "runs/rfd_af2ig_$TAG" --outdir "runs/rfd_passers_$TAG" \
    --manifest "results/rfd_af2_passers_$TAG.csv" --max-pae "$MAXPAE"
npass=$(ls runs/rfd_passers_$TAG/*.pdb 2>/dev/null | wc -l)
echo "passers: $npass"
[ "$npass" -eq 0 ] && { echo "no passers; stopping"; exit 0; }
echo "=== [$TAG] canonicalize His433->B/100 ==="
"$MD" scripts/canon_renumber.py --indir "runs/rfd_passers_$TAG" --outdir "runs/rfd_passers_${TAG}_canon"
echo "=== [$TAG] obj1_switch: install + score His433 pH switch ($npass designs, $WORKERS workers) ==="
"$PR" scripts/obj1_switch.py --pdbs runs/rfd_passers_${TAG}_canon/*.pdb \
    --outdir "runs/rfd_obj1_$TAG" --out "results/rfd_obj1_switch_$TAG.csv" --workers "$WORKERS"
echo "BETA_SWITCH_${TAG}_DONE $(date -u +%FT%TZ)"
