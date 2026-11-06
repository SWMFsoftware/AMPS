#!/usr/bin/env bash
# Run the Stage-9 transported-daughter/general-NTC compatibility gate at all
# accepted decompositions.  The caller supplies one already configured,
# immutable executable.  This runner never edits input, compiles AMPS, removes
# build products, or overwrites evidence; failed case directories remain
# available for diagnosis.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: run_stage9_integrated_tests.sh \
  --executable /absolute/path/to/stage9-amps \
  --output-root /new/path/to/stage9-collision-validation \
  [--mpiexec mpiexec] [--chemistry-validation-events 200000]

Runs --validate-transported-daughter-collision at exactly 1, 4, and 12 MPI
ranks. The executable must be the DIM1 Stage-9 build with individual particle
corrections, Stage-5 windows, Stage-6 packets, Stage-7 channel importance,
Stage-8 total-reaction forcing, and the production general NTC collision model.

The output root must not exist. Each case receives its own command.txt,
run.log, error.log, mesh, and focused runtime artifacts. A successful process
exit is not sufficient: the strict parser independently checks transported
pointer identity, collision feedback, corrections, completed trajectories,
probabilities, represented rate, and momentum closure.
EOF
}

executable=
output_root=
mpi_launcher=mpiexec
chemistry_validation_events=200000

while (( $# > 0 )); do
  case "$1" in
    --executable)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      executable=$2
      shift 2
      ;;
    --output-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      output_root=$2
      shift 2
      ;;
    --mpiexec)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      mpi_launcher=$2
      shift 2
      ;;
    --chemistry-validation-events)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      chemistry_validation_events=$2
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      printf 'unknown Stage-9 integrated-runner option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -n "$executable" && -n "$output_root" ]] || {
  printf '%s\n' 'both --executable and --output-root are required' >&2
  usage >&2
  exit 2
}
[[ "$chemistry_validation_events" =~ ^[1-9][0-9]*$ ]] || {
  printf 'invalid --chemistry-validation-events value: %s\n' \
    "$chemistry_validation_events" >&2
  exit 2
}
(( chemistry_validation_events >= 1000 )) || {
  printf 'Stage-9 chemistry-validation-events must be at least 1000: %s\n' \
    "$chemistry_validation_events" >&2
  exit 2
}
[[ -f "$executable" && -x "$executable" ]] || {
  printf 'Stage-9 executable is not executable: %s\n' "$executable" >&2
  exit 2
}
resolved_launcher=$(command -v -- "$mpi_launcher" || true)
[[ -n "$resolved_launcher" ]] || {
  printf 'MPI launcher is not on PATH: %s\n' "$mpi_launcher" >&2
  exit 2
}
[[ ! -e "$output_root" ]] || {
  printf 'Stage-9 output root already exists: %s\n' "$output_root" >&2
  exit 2
}

executable=$(readlink -f -- "$executable")
resolved_launcher=$(readlink -f -- "$resolved_launcher")
mkdir -p -- "$output_root"
output_root=$(cd -- "$output_root" && pwd -P)

run_case() {
  local ranks=$1
  local case_dir="$output_root/mpi${ranks}/transported-daughter-collision"
  mkdir -p -- "$case_dir"
  {
    printf 'OMP_NUM_THREADS=1'
    printf ' %q' "$resolved_launcher" -n "$ranks" "$executable" \
      --validate-transported-daughter-collision \
      --chemistry-validation-events "$chemistry_validation_events"
    printf '\n'
  } >"$case_dir/command.txt"

  printf 'START Stage-9 transported-daughter-collision ranks=%s directory=%s\n' \
    "$ranks" "$case_dir"
  if ! (
    cd -- "$case_dir"
    OMP_NUM_THREADS=1 "$resolved_launcher" -n "$ranks" "$executable" \
      --validate-transported-daughter-collision \
      --chemistry-validation-events "$chemistry_validation_events" \
      >run.log 2>error.log
  ); then
    printf 'FAIL Stage-9 transported-daughter-collision ranks=%s; inspect %s\n' \
      "$ranks" "$case_dir" >&2
    return 1
  fi
}

for ranks in 1 4 12; do
  run_case "$ranks"
done

parser_arguments=()
for ranks in 1 4 12; do
  parser_arguments+=(
    --collision-log "$ranks"
    "$output_root/mpi${ranks}/transported-daughter-collision/run.log"
  )
done
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage9_transported_collision_logs.py" \
  "${parser_arguments[@]}"

printf 'PASS Stage-9 transported-daughter collision runtime archive=%s\n' \
  "$output_root"
