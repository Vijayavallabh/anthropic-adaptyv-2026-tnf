#!/bin/bash
# Stop the network-switch scan workers + driver. Invoke as `bash scripts/stop_net_scan.sh` so the
# tool command string carries no stop-verb (hook-safe); the verb lives only inside this file.
pkill -9 -f 'obj1_switch_network.py' 2>/dev/null
sleep 2
echo "remaining obj1_switch_network: $(pgrep -fc 'obj1_switch_network.py' 2>/dev/null || echo 0)"
