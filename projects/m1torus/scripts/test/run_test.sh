#!/bin/bash
# Build one unit test (default test_cloudcollisions) against phantom's objects and run it.
#   run_test.sh [test_cloudcollisions|test_cloudfeed]
# Needs a previous `make SETUP=m1torus` (any IND_TIMESTEPS) in this tree; test_cloudfeed
# needs it with SRCINJECT=inject_cloudfeed.f90.
set -euo pipefail
name=${1:-test_cloudcollisions}
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
build=$root/build
out=${TMPDIR:-/tmp}/$name
mkdir -p "$out"
# every object except the ones holding a main program (phantom, setup, moddump, ...)
objs=$(for o in "$build"/*.o; do nm "$o" 2>/dev/null | grep -qE ' T main$' || echo "$o"; done)
gfortran -O2 -fdefault-real-8 -fdefault-double-8 -fopenmp -I"$build" -J"$out" \
   "$root/projects/m1torus/scripts/test/$name.f90" $objs -o "$out/$name"
"$out/$name"
