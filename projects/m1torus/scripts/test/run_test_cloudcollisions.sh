#!/bin/bash
# Build the cloud-collision test against phantom's objects and run it.
# Needs a previous `make SETUP=m1torus` (any IND_TIMESTEPS) in this tree.
set -euo pipefail
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
build=$root/build
out=${TMPDIR:-/tmp}/test_cloudcollisions
mkdir -p "$out"
# every object except the main programs
objs=$(ls "$build"/*.o | grep -v -E '/(phantom|phantomsetup|phantommoddump|phantomanalysis)\.o$')
gfortran -O2 -fdefault-real-8 -fdefault-double-8 -fopenmp -I"$build" -J"$out" \
   "$root/projects/m1torus/scripts/test/test_cloudcollisions.f90" $objs -o "$out/test_cloudcollisions"
"$out/test_cloudcollisions"
