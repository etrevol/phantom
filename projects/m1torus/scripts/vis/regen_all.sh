#!/usr/bin/env bash
# Regenerate every figure and movie of the finished m1torus runs, after a change of the
# naming and style rules (projects/m1torus/FIGURES.md). Slow: the movies take most of it.
#
#   projects/m1torus/scripts/vis/regen_all.sh [--no-movies]
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../../.." && pwd)
A=$root/results/m1torus/A/runs
B=$root/results/m1torus/B/runs
movies=yes
[[ ${1:-} == --no-movies ]] && movies=no
py() { nice -n 10 python3 "$here/$1" "${@:2}" > /dev/null 2>&1 || echo "FAILED: $*"; }

# collisionless runs and their saturated windows [orbits]
declare -A WIN=(
  [20261027-095659-canon-16k]=110,150      [20261027-141227-canon-16k-eps0]=110,150
  [20261027-215121-canon-64k]=110,150      [20261028-121511-thin30-16k]=110,150
  [20261028-150605-canon-16k-fixeps001]=110,150
  [20261029-100846-mtor03-16k]=80,120      [20261029-122207-mtor006-16k]=200,250
  [20261028-144810-weak45-16k]=350,400
)
for r in "${!WIN[@]}"; do
  [[ -f $A/$r/status ]] && grep -q COMPLETED "$A/$r/status" || continue
  py vism1.py --run "$A/$r" --window "${WIN[$r]}"
  # orbit figures need the full dumps, which some early runs no longer keep
  compgen -G "$A/$r/m1torus_0*" > /dev/null && py vis_orbits.py --run "$A/$r" --window "${WIN[$r]}"
  echo "A $r"
done

# cross-run comparisons of the collisionless model
c=$A/20261027-141227-canon-16k-eps0
py compare_runs.py --runs "$c" "$A/20261027-215121-canon-64k" --window 110,150 --out "$A/compare_16k_64k.png"
py compare_runs.py --runs "$A/20261027-095659-canon-16k" "$c" --window 110,150 --out "$A/compare_eps_16k.png"
py compare_runs.py --runs "$c" "$A/20261028-150605-canon-16k-fixeps001" --window 110,150 --out "$A/compare_softening.png"
py compare_runs.py --runs "$c" "$A/20261028-121511-thin30-16k" --window 110,150 --out "$A/compare_thin30.png"
py compare_runs.py --runs "$A/20261029-122207-mtor006-16k" "$c" "$A/20261029-100846-mtor03-16k" \
   --window 100,120 --out "$A/compare_mass.png"
rm -f "$A/compare_mtor.png" "$A/compare_mtor006.png"     # superseded by compare_mass.png
# the collisionless model on NGC 613, as the reference of B
py compare_obs.py "$c/m1torus_00100" --out "$root/results/m1torus/B/compare_A3"

# dissipative runs: time series, mode panels, synthetic maps of the last dump
for d in "$B"/2*; do
  [[ -f $d/status ]] && grep -qE "COMPLETED|FAILED" "$d/status" || continue
  last=$(ls "$d"/m1torus_[0-9][0-9][0-9][0-9][0-9] | tail -1)
  t1=$(python3 -c "import sys; sys.path.insert(0,'$root/projects/common'); import phantom_io as p; print(round(p.read_dump('$last')['header']['time']/6.2831853))")
  py vis_b.py --run "$d"
  py vism1.py --run "$d" --window "$((t1 - 3)),$t1"
  py compare_obs.py "$last"
  py compare_circinus.py "$last"
  echo "B $(basename "$d")"
done

[[ $movies == yes ]] || exit 0
for r in "${!WIN[@]}"; do
  [[ -f $A/$r/status ]] && grep -q COMPLETED "$A/$r/status" || continue
  compgen -G "$A/$r/m1torus_0*" > /dev/null || continue
  ev=2; [[ $r == *64k* ]] && ev=3
  py anim_particles.py --run "$A/$r" --every $ev
  py anim_particles.py --run "$A/$r" --every $ev --corotate
  echo "movie A $r"
done
for d in "$B"/2*; do
  [[ -f $d/status ]] && grep -qE "COMPLETED|FAILED" "$d/status" || continue
  py anim_particles.py --run "$d"
  echo "movie B $(basename "$d")"
done
