#!/usr/bin/env bash
# Continue a finished or stopped m1torus run in place, from its last full dump.
#
#   projects/m1torus/scripts/run/continue.sh RUN_DIR N_ORBITS [--threads N] [--fg]
#
# Every full phantom dump is a complete restart state (double precision positions and
# velocities, sinks, header), so no separate restart file is needed. N_ORBITS is the new
# end time; phantom starts a new <prefix>NN.ev and the live analysis appends to
# m1torus_m1.ev; the readers keep the last value at each time.
set -euo pipefail
run=$(cd "$1" && pwd); orbits=$2; shift 2
threads=8; fg=no
while [[ $# -gt 0 ]]; do
  case $1 in
    --threads) threads=$2; shift 2 ;;
    --fg) fg=yes; shift ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
cd "$run"
[[ $(cat status) == RUNNING ]] && { echo "$run is still running" >&2; exit 1; }
last=$(ls m1torus_[0-9][0-9][0-9][0-9][0-9] | tail -1)
tmax=$(python3 -c "import math; print(f'{$orbits*2*math.pi:.10g}')")
sed -i -E "s|^( *dumpfile *= *)[^!]*(!.*)?$|\1$last   \2|; s|^( *tmax *= *)[^!]*(!.*)?$|\1$tmax   \2|" m1torus.in
echo "$(date -Iseconds) continued from $last to $orbits orbits" >> continued.log
launch='export OMP_NUM_THREADS='$threads'; echo RUNNING > status
  if ./phantom m1torus.in >> run.log 2>&1; then echo COMPLETED > status; else echo FAILED > status; fi'
if [[ $fg == yes ]]; then
  bash -c "$launch"; echo "$(cat status): $run"
else
  setsid nohup bash -c "$launch" > /dev/null 2>&1 &
  echo "continued: $run from $last to $orbits orbits"
fi
