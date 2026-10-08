#!/bin/bash
# Stop ONLY the candidate-MD campaign (cmdline contains "md_candidates"); leaves the baseline MD
# (md_submission) and everything else running. Run as `bash scripts/kill_candidate_md.sh` so this
# runner's own cmdline does NOT contain the pattern (no self-match).
pkill -9 -f "md_candidates" 2>/dev/null
sleep 2
echo "remaining md_candidates procs: $(pgrep -fc 'md_candidates' 2>/dev/null || echo 0)"
echo "baseline md_submission procs still alive: $(pgrep -fc 'md_submission' 2>/dev/null || echo 0)"
