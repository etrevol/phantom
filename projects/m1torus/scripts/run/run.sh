#!/usr/bin/env bash
# Build, set up and launch one m1torus run.
#
#   projects/m1torus/scripts/run/run.sh --name canon --campaign A -- np=131072 n_orbits=300
#   projects/m1torus/scripts/run/run.sh --name n16k --kind tests --fg -- np=16384 n_orbits=5
#
# key=value pairs after -- replace entries of the .setup that phantomsetup writes with its
# defaults (the paper's canonical run 0.5-60-0.3); phantom .in options go after --in, e.g.
# --in tree_accuracy=0.3. The run lands in results/m1torus/<campaign>/<kind>/<stamp>-<name>/
# with its own binaries, inputs, provenance.json, run.log and a status file.
# --restart DUMP continues from a full dump with that run's .setup and .in, for
# numerical tests on an evolved state; give the new end time with --in tmax=...
# --make VAR=VALUE passes a build flag (e.g. IND_TIMESTEPS=no; it forces a rebuild).
# --mod key=value (with --restart) passes the dump through moddump_m1torus first, e.g.
# --mod mode=fixed (sink -> fixed point mass) or --mod mode=soften --mod eps_c=0.05.
# --maxp N caps the particle arrays of a run that adds particles (cloud inflow): with
# injection phantom allocates 5.2M otherwise; kept in ./maxp for continue.sh.
set -euo pipefail

name=run
campaign=A
kind=runs
threads=8          # 8 physical cores; 16 hyperthreads are no faster
fg=no
restart=
maxp=
make_opts=()
mod_opts=()
in_opts=()
while [[ $# -gt 0 ]]; do
  case $1 in
    --name) name=$2; shift 2 ;;
    --campaign) campaign=$2; shift 2 ;;
    --kind) kind=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    --fg) fg=yes; shift ;;
    --in) in_opts+=("$2"); shift 2 ;;
    --make) make_opts+=("$2"); shift 2 ;;
    --mod) mod_opts+=("$2"); shift 2 ;;
    --maxp) maxp=$2; shift 2 ;;
    --restart) restart=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); shift 2 ;;
    --) shift; break ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
setup_opts=("$@")

root=$(cd "$(dirname "$0")" && while [[ ! -f build/Makefile_setups && $PWD != / ]]; do cd ..; done; pwd)
[[ -f $root/build/Makefile_setups ]] || { echo "phantom root not found" >&2; exit 1; }
dir=$root/results/m1torus/$campaign/$kind/$(date +%Y%m%d-%H%M%S)-$name

export SYSTEM=${SYSTEM:-gfortran}
( cd "$root" && make SETUP=m1torus "${make_opts[@]}" >/dev/null \
  && make SETUP=m1torus "${make_opts[@]}" setup >/dev/null \
  && make SETUP=m1torus "${make_opts[@]}" moddump >/dev/null ) \
  || { echo "build failed: cd $root && make SETUP=m1torus ${make_opts[*]:-}" >&2; exit 1; }

mkdir -p "$dir"
cp "$root/bin/phantom" "$root/bin/phantomsetup" "$root/bin/phantommoddump" "$dir/"
cd "$dir"
[[ -n $maxp ]] && echo "$maxp" > maxp

# replace key = value in an input file; a key that is not there is an error, since
# phantom would ignore it silently
set_opt() {
  local file=$1 key=${2%%=*} val=${2#*=}
  grep -qE "^ *$key *=" "$file" || { echo "no option '$key' in $file" >&2; exit 1; }
  sed -i -E "s|^( *$key *= *)[^!]*(!.*)?$|\1$val   \2|" "$file"
}

if [[ -n $restart ]]; then
  src=$(dirname "$restart")
  cp "$src/m1torus.setup" "$src/m1torus.in" .
  if [[ ${#mod_opts[@]} -gt 0 ]]; then
    num=${restart##*_}
    cp "$restart" "src_$num"
    cp m1torus.in src.in
    printf '%s\n' "${mod_opts[@]}" | sed 's/=/ = /' > moddump_m1torus.opts
    rm m1torus.in
    # an .in from a build without the new options (e.g. inflow) is rewritten with
    # their defaults on the first call; the second call then goes through
    md() { ./phantommoddump "src_$num" "m1torus_$num" ${maxp:+--maxp=$maxp} > setup.log 2>&1; }
    md || { grep -q "try again using revised input file" setup.log && md; } \
      || { echo FAILED > status; tail -20 setup.log; exit 1; }
    rm "src_$num" src.in
  else
    cp "$restart" .
    echo "restart from $restart" > setup.log
  fi
  set_opt m1torus.in "dumpfile=$(basename "$restart")"
else
  ./phantomsetup m1torus > /dev/null 2>&1 || true        # writes the default .setup
  for kv in "${setup_opts[@]}"; do set_opt m1torus.setup "$kv"; done
  np=$(awk '$1=="np"{print $3}' m1torus.setup)
  ./phantomsetup m1torus --maxp=$(( np + 1000 )) > setup.log 2>&1 \
    || { echo FAILED > status; tail -20 setup.log; exit 1; }
fi
for kv in "${in_opts[@]}"; do set_opt m1torus.in "$kv"; done

cat > provenance.json <<EOF
{
  "created": "$(date -Iseconds)",
  "git": "$(git -C "$root" rev-parse --short HEAD)",
  "dirty_files": $(git -C "$root" status --porcelain -- src build projects/m1torus | wc -l),
  "make": "make SETUP=m1torus ${make_opts[*]:-} SYSTEM=$SYSTEM",
  "threads": $threads,
  "setup_opts": "${setup_opts[*]:-}",
  "restart": "${restart}",
  "mod_opts": "${mod_opts[*]:-}",
  "maxp": "${maxp}",
  "in_opts": "${in_opts[*]:-}"
}
EOF

launch='export OMP_NUM_THREADS='$threads'; echo RUNNING > status
  if ./phantom m1torus.in $([[ -f maxp ]] && echo --maxp=$(cat maxp)) > run.log 2>&1; then echo COMPLETED > status; else echo FAILED > status; fi'
if [[ $fg == yes ]]; then
  bash -c "$launch"
  echo "$(cat status): $dir"
else
  setsid nohup bash -c "$launch" > /dev/null 2>&1 &
  echo "launched: $dir"
fi
