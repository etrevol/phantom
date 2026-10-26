#!/usr/bin/env bash
# Build, set up and launch one PPI torus run.
#
#   scripts/ppi/run/run.sh --name unseeded -- --np 1e5 --n-orbits 40
#   scripts/ppi/run/run.sh --name m2 --kind tests --fg -- --pert-kind single --pert-amp 1e-3
#
# Everything after -- goes to scripts/ppi/theory/torus_model.py (see its --help).
# The run lands in results/<kind>/<YYYYMMDD-HHMMSS>-<name>/ with its own binaries, the
# .setup and .in, provenance.json, run.log and a status file (RELAXING, RUNNING,
# COMPLETED, FAILED). With relax = T (the default) stage 1 runs in relax/, and
# phantommoddump spins the last relaxed dump up into the run directory at t = 0.
set -euo pipefail

name=run
kind=runs
threads=8          # 8 physical cores; 16 hyperthreads are no faster
ind=yes
fg=no
maxp=
dry=no
while [[ $# -gt 0 ]]; do
  case $1 in
    --name) name=$2; shift 2 ;;
    --kind) kind=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    --ind) ind=$2; shift 2 ;;
    --maxp) maxp=$2; shift 2 ;;
    --fg) fg=yes; shift ;;
    --dry-run) dry=yes; shift ;;
    --) shift; break ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done

root=$(cd "$(dirname "$0")" && while [[ ! -f build/Makefile_setups && $PWD != / ]]; do cd ..; done; pwd)
[[ -f $root/build/Makefile_setups ]] || { echo "phantom root not found" >&2; exit 1; }

dir=$root/results/$kind/$(date +%Y%m%d-%H%M%S)-$name
setup_args=("$@")

if [[ $dry == yes ]]; then
  python3 "$root/scripts/ppi/theory/torus_model.py" -o /dev/null "${setup_args[@]}"
  echo "would run in $dir (ind=$ind, threads=$threads)"
  exit 0
fi

# Phantom's Makefile cleans by itself when SETUP or a preprocessor flag changes.
export SYSTEM=${SYSTEM:-gfortran}
( cd "$root" && make SETUP=ppitorus IND_TIMESTEPS=$ind >/dev/null \
  && make SETUP=ppitorus IND_TIMESTEPS=$ind setup >/dev/null \
  && make SETUP=ppitorus IND_TIMESTEPS=$ind moddump >/dev/null ) \
  || { echo "build failed: cd $root && make SETUP=ppitorus IND_TIMESTEPS=$ind" >&2; exit 1; }

mkdir -p "$dir"
cp "$root/bin/phantom" "$root/bin/phantomsetup" "$root/bin/phantommoddump" "$dir/"
cd "$dir"
python3 "$root/scripts/ppi/theory/torus_model.py" -o ppitorus.setup "${setup_args[@]}" > model.txt
np_target=$(awk '$1=="np"{print $3}' ppitorus.setup)
relax=$(awk '$1=="relax"{print $3}' ppitorus.setup)
maxp=${maxp:-$(( np_target * 12 / 10 + 1000 ))}
if [[ $relax == T ]]; then
  mkdir -p relax
  cp ppitorus.setup relax/
  ( cd relax && ../phantomsetup ppitorus --maxp="$maxp" ) > setup.log 2>&1 \
    || { echo FAILED > status; tail -20 setup.log; exit 1; }
else
  ./phantomsetup ppitorus --maxp="$maxp" > setup.log 2>&1 \
    || { echo FAILED > status; tail -20 setup.log; exit 1; }
fi

sha=$(git -C "$root" rev-parse --short HEAD)
dirty=$(git -C "$root" status --porcelain -- src build scripts/ppi | wc -l)
cat > provenance.json <<EOF
{
  "created": "$(date -Iseconds)",
  "git": "$sha",
  "dirty_files": $dirty,
  "make": "make SETUP=ppitorus IND_TIMESTEPS=$ind SYSTEM=$SYSTEM",
  "threads": $threads,
  "maxp": $maxp,
  "setup_args": "$(printf '%s ' "${setup_args[@]}")"
}
EOF

launch='export OMP_NUM_THREADS='$threads'
  if [[ -d relax ]]; then
    echo RELAXING > status
    ( cd relax && ../phantom ppitorus.in > run.log 2>&1 ) || { echo FAILED > status; exit 1; }
    last=$(ls relax/ppitorus_[0-9][0-9][0-9][0-9][0-9] | tail -1)
    ./phantommoddump "$last" ppitorus_00000 0.0 --maxp='$maxp' > moddump.log 2>&1 \
      || { echo FAILED > status; exit 1; }
  fi
  echo RUNNING > status
  if ./phantom ppitorus.in > run.log 2>&1; then echo COMPLETED > status; else echo FAILED > status; fi
  OPENBLAS_NUM_THREADS=1 '"$root"'/scripts/ppi/vis/gallery.sh . > figs.log 2>&1'
if [[ $fg == yes ]]; then
  bash -c "$launch"
  echo "$(cat status): $dir"
else
  setsid nohup bash -c "$launch" > /dev/null 2>&1 &
  echo "launched: $dir"
  echo "watch:    tail -f $dir/run.log"
fi
