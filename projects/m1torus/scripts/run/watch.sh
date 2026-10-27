#!/usr/bin/env bash
# One report on the running m1torus runs, then exit: after PERIOD seconds (default 1200),
# at once when a run finishes or fails, or when status.py flags a problem. Meant to be
# re-armed after each report, so a watcher never outlives a 30 min monitor timeout.
#
#   projects/m1torus/scripts/run/watch.sh [RUN_DIR ...]
here=$(cd "$(dirname "$0")" && pwd)
period=${PERIOD:-1200}
start=$(date +%s)
while true; do
  out=$(python3 "$here/status.py" "$@" 2>&1); rc=$?
  if [[ $rc -ne 0 ]] || ! grep -q RUNNING <<<"$out" || (( $(date +%s) - start >= period )); then
    echo "[$(date +%H:%M)] $out"
    exit 0
  fi
  sleep 60
done
