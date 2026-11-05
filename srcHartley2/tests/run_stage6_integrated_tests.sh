#!/usr/bin/env bash
# Run the initialized-mesh Stage-6 packet validator and its Stage-4/5
# regression dependencies at every accepted MPI decomposition.  The caller
# supplies one frozen correction-on executable; this runner never configures,
# compiles, edits an input, or overwrites an existing evidence archive.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: run_stage6_integrated_tests.sh \
  --executable /absolute/path/to/stage6-amps \
  --output-root /new/path/to/stage6-validation \
  [--mpiexec mpiexec] [--chemistry-validation-events 200000]

Runs correction-window, exact chemistry, and daughter-packet validators at
exactly 1, 4, and 12 MPI ranks. The executable must be a DIM1 radial build
with individual corrections, Bernoulli chemistry, Stage-5 windows, and
Stage-6 packet mode enabled with a maximum of at least eight. The output root
must not exist; failed or partial evidence is never overwritten.
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
      printf 'unknown Stage-6 runner option: %s\n' "$1" >&2
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
  printf 'Stage-6 chemistry-validation-events must be at least 1000: %s\n' \
    "$chemistry_validation_events" >&2
  exit 2
}
[[ -f "$executable" && -x "$executable" ]] || {
  printf 'Stage-6 executable is not executable: %s\n' "$executable" >&2
  exit 2
}
resolved_launcher=$(command -v -- "$mpi_launcher" || true)
[[ -n "$resolved_launcher" ]] || {
  printf 'MPI launcher is not on PATH: %s\n' "$mpi_launcher" >&2
  exit 2
}
[[ ! -e "$output_root" ]] || {
  printf 'Stage-6 output root already exists: %s\n' "$output_root" >&2
  exit 2
}

executable=$(readlink -f -- "$executable")
resolved_launcher=$(readlink -f -- "$resolved_launcher")
mkdir -p -- "$output_root"
output_root=$(cd -- "$output_root" && pwd -P)

run_case() {
  local ranks=$1
  local name=$2
  shift 2
  local case_dir="$output_root/mpi${ranks}/$name"
  mkdir -p -- "$case_dir"
  {
    printf 'OMP_NUM_THREADS=1'
    printf ' %q' "$resolved_launcher" -n "$ranks" "$executable" "$@"
    printf '\n'
  } >"$case_dir/command.txt"

  printf 'START Stage-6 %s ranks=%s directory=%s\n' \
    "$name" "$ranks" "$case_dir"
  if ! (
    cd -- "$case_dir"
    OMP_NUM_THREADS=1 "$resolved_launcher" -n "$ranks" "$executable" "$@" \
      >run.log 2>error.log
  ); then
    printf 'FAIL Stage-6 %s ranks=%s; inspect %s\n' \
      "$name" "$ranks" "$case_dir" >&2
    return 1
  fi
}

for ranks in 1 4 12; do
  # Stage 6 relies on both exact corrected chemistry (Stage 4) and the normal
  # post-step correction controller (Stage 5), so repeat their initialized
  # production-path validators with the very same frozen executable.
  run_case "$ranks" correction-windows --validate-correction-windows
  run_case "$ranks" chemistry --validate-chemistry \
    --chemistry-validation-events "$chemistry_validation_events"
  run_case "$ranks" daughter-packets --validate-daughter-packets \
    --chemistry-validation-events "$chemistry_validation_events"
done

# Reuse the Stage-5 strict parser rather than weakening those already accepted
# gates in a second implementation.
stage5_args=()
packet_args=(--expected-events "$chemistry_validation_events")
for ranks in 1 4 12; do
  stage5_args+=(
    --window-log "$ranks" "$output_root/mpi${ranks}/correction-windows/run.log"
    --chemistry-log "$ranks" "$output_root/mpi${ranks}/chemistry/run.log"
  )
  packet_args+=(
    --packet-log "$ranks" "$output_root/mpi${ranks}/daughter-packets/run.log"
  )
done
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage5_integrated_logs.py" "${stage5_args[@]}"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage6_daughter_packet_logs.py" "${packet_args[@]}"

printf 'PASS Stage-6 daughter-packet runtime archive=%s\n' "$output_root"
