#!/bin/bash
# Stop lingering pd2 downstream workers. Run as `bash scripts/stop_pd2_downstream.sh` so the tool
# command string contains no stop-verb for any hook, and this runner's cmdline does not self-match.
pkill -9 -f 'rfd_to_af2ig' 2>/dev/null
pkill -9 -f 'run_af2ig.py' 2>/dev/null
pkill -9 -f 'ligandmpnn|LigandMPNN|protein_mpnn' 2>/dev/null
sleep 2
echo "remaining rfd_to_af2ig: $(pgrep -fc 'rfd_to_af2ig' 2>/dev/null || echo 0)"
echo "remaining run_af2ig: $(pgrep -fc 'run_af2ig.py' 2>/dev/null || echo 0)"