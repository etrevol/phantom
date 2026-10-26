#!/usr/bin/env bash
# Standard figures of one run: history, modes against theory, mid-plane at the last
# frame (plain and deviation), meridional cuts. Movies are separate (vis2d.py --animate).
#
#   scripts/ppi/vis/gallery.sh RUN_DIR
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
run=${1:-.}
python3 "$here/vishst.py" --run "$run"
python3 "$here/vismodes.py" --run "$run" --m 1,2,3,4 | tee "$run/modes.txt"
python3 "$here/vis2d.py" --run "$run" --mode all
