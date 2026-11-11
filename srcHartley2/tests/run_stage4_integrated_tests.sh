#!/usr/bin/env bash
# Launch the already-built correction-on comet executable at the three Stage-4
# MPI decompositions, then validate its compact records.  This script never
# configures or compiles AMPS and never edits an input file or build tree.
# Keeping build selection outside the runner makes executable provenance and
# the clean-build procedure explicit responsibilities of the caller.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: run_stage4_integrated_tests.sh \
  --executable /absolute/path/to/amps \
  --output-root /new/path/to/stage4-validation \
  [--mpiexec mpiexec] [--validation-particles 200000] \
  [--chemistry-validation-events 200000]

Runs the local-timestep, weight-adaptation/core-compatibility, and
chemistry-ensemble validators at exactly 1, 4, and 12 MPI ranks. The
executable must already have individual particle-weight correction,
species-local timesteps, species-local weights, and local-timestep face
resampling compiled in. The output root must not exist: archived evidence is
never overwritten or silently reused.

This runner does not invoke make, ampsConfig.pl, or modify build/, input/, or
the repository-root generated input files.
EOF
}

executable=
output_root=
mpi_launcher=mpiexec
validation_particles=200000
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
    --validation-particles)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      validation_particles=$2
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
      printf 'unknown Stage-4 runner option: %s\n' "$1" >&2
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
[[ "$validation_particles" =~ ^[1-9][0-9]*$ ]] || {
  printf 'invalid --validation-particles value: %s\n' \
    "$validation_particles" >&2
  exit 2
}
[[ "$chemistry_validation_events" =~ ^[1-9][0-9]*$ ]] || {
  printf 'invalid --chemistry-validation-events value: %s\n' \
    "$chemistry_validation_events" >&2
  exit 2
}
[[ -f "$executable" && -x "$executable" ]] || {
  printf 'Stage-4 executable is not an executable regular file: %s\n' \
    "$executable" >&2
  exit 2
}

resolved_launcher=$(command -v -- "$mpi_launcher" || true)
[[ -n "$resolved_launcher" ]] || {
  printf 'MPI launcher is not on PATH: %s\n' "$mpi_launcher" >&2
  exit 2
}
executable=$(readlink -f -- "$executable")
resolved_launcher=$(readlink -f -- "$resolved_launcher")

# Refuse an existing path even when it is empty. A second attempt belongs in a
# new archive so a failed case cannot be mistaken for evidence from a retry.
[[ ! -e "$output_root" ]] || {
  printf 'Stage-4 output root already exists: %s\n' "$output_root" >&2
  exit 2
}
mkdir -p -- "$output_root"
output_root=$(cd -- "$output_root" && pwd -P)

run_case() {
  local ranks=$1
  local name=$2
  shift 2
  local case_dir="$output_root/mpi${ranks}/$name"
  mkdir -p -- "$case_dir"

  # command.txt is provenance only. Use shell-escaped arguments so spaces in
  # an executable/output path cannot make the record ambiguous.
  {
    printf 'OMP_NUM_THREADS=1'
    printf ' %q' "$resolved_launcher" -n "$ranks" "$executable" "$@"
    printf '\n'
  } >"$case_dir/command.txt"

  printf 'START Stage-4 %s ranks=%s directory=%s\n' \
    "$name" "$ranks" "$case_dir"
  if ! (
    cd -- "$case_dir"
    OMP_NUM_THREADS=1 "$resolved_launcher" -n "$ranks" "$executable" "$@" \
      >run.log 2>error.log
  ); then
    printf 'FAIL Stage-4 %s ranks=%s; inspect %s\n' \
      "$name" "$ranks" "$case_dir" >&2
    return 1
  fi
}

run_expected_failure() {
  local ranks=$1
  local name=$2
  local expected_message=$3
  shift 3
  local case_dir="$output_root/mpi${ranks}/$name"
  mkdir -p -- "$case_dir"
  {
    printf 'OMP_NUM_THREADS=1'
    printf ' %q' "$resolved_launcher" -n "$ranks" "$executable" "$@"
    printf '\n'
  } >"$case_dir/command.txt"

  printf 'START Stage-4 expected-failure %s ranks=%s directory=%s\n' \
    "$name" "$ranks" "$case_dir"
  if (
    cd -- "$case_dir"
    OMP_NUM_THREADS=1 "$resolved_launcher" -n "$ranks" "$executable" "$@" \
      >run.log 2>error.log
  ); then
    printf 'FAIL Stage-4 %s ranks=%s unexpectedly returned success\n' \
      "$name" "$ranks" >&2
    return 1
  fi
  if ! grep -Fq -- "$expected_message" "$case_dir/run.log" \
      "$case_dir/error.log"; then
    printf 'FAIL Stage-4 %s ranks=%s did not emit expected setup error\n' \
      "$name" "$ranks" >&2
    return 1
  fi
}

for ranks in 1 4 12; do
  run_case "$ranks" local-timestep \
    --validate-local-timestep \
    --validation-particles "$validation_particles" \
    --validation-mover both
  run_case "$ranks" local-timestep-unit \
    --validate-local-timestep \
    --validation-particles "$validation_particles" \
    --validation-mover both \
    --validation-unit-corrections
  run_case "$ranks" weight-adaptation --validate-weight-adaptation
  run_case "$ranks" chemistry \
    --validate-chemistry \
    --chemistry-validation-events "$chemistry_validation_events"
done

run_expected_failure 1 unsupported-simplified-ntc \
  "simplified NTC does not support individual particle-weight corrections" \
  --validate-unsupported-simplified-ntc

python_args=()
for ranks in 1 4 12; do
  python_args+=(
    --local-log "$ranks" "$output_root/mpi${ranks}/local-timestep/run.log"
    --unit-local-log "$ranks" \
      "$output_root/mpi${ranks}/local-timestep-unit/run.log"
    --adaptation-log "$ranks" "$output_root/mpi${ranks}/weight-adaptation/run.log"
    --chemistry-log "$ranks" "$output_root/mpi${ranks}/chemistry/run.log"
  )
done

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage4_integrated_logs.py" "${python_args[@]}"

printf 'PASS Stage-4 correction-on runtime archive=%s\n' "$output_root"
