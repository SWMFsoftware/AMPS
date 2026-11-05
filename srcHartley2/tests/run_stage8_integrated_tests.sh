#!/usr/bin/env bash
# Run the mover-integrated Stage-8 deterministic forced-segment validator and
# every Stage-4/5/6/7 dependency at the accepted MPI decompositions.  The
# caller supplies one frozen Stage-8-enabled executable.  This script never
# edits an input, configures AMPS, compiles, removes a build, or reuses an
# evidence directory, so a validation archive remains attributable to one
# immutable executable and generated policy set.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: run_stage8_integrated_tests.sh \
  --executable /absolute/path/to/stage8-amps \
  --output-root /new/path/to/stage8-validation \
  [--mpiexec mpiexec] [--validation-particles 200000] \
  [--chemistry-validation-events 200000]

Runs the Stage-4 local-timestep mover, unit-correction mover, particle-weight
adaptation, correction-window, exact-chemistry, daughter-packet,
channel-importance, and mover-integrated total-reaction-forcing validators at
exactly 1, 4, and 12 MPI ranks. The executable must be a DIM1 radial build with
individual corrections, Bernoulli chemistry, Stage-5 windows, Stage-6 packets,
Stage-7 channel importance, and Stage-8 deterministic two-branch forcing
enabled.

The output root must not exist. Failed/partial evidence is retained in place
and is never overwritten or mixed with a retry. Each case receives separate
command, stdout, stderr, and working directories.
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
    --chemistry-validation-events)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      chemistry_validation_events=$2
      shift 2
      ;;
    --validation-particles)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      validation_particles=$2
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      printf 'unknown Stage-8 runner option: %s\n' "$1" >&2
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
[[ "$validation_particles" =~ ^[1-9][0-9]*$ ]] || {
  printf 'invalid --validation-particles value: %s\n' \
    "$validation_particles" >&2
  exit 2
}
(( chemistry_validation_events >= 1000 )) || {
  printf 'Stage-8 chemistry-validation-events must be at least 1000: %s\n' \
    "$chemistry_validation_events" >&2
  exit 2
}
[[ -f "$executable" && -x "$executable" ]] || {
  printf 'Stage-8 executable is not executable: %s\n' "$executable" >&2
  exit 2
}
resolved_launcher=$(command -v -- "$mpi_launcher" || true)
[[ -n "$resolved_launcher" ]] || {
  printf 'MPI launcher is not on PATH: %s\n' "$mpi_launcher" >&2
  exit 2
}
[[ ! -e "$output_root" ]] || {
  printf 'Stage-8 output root already exists: %s\n' "$output_root" >&2
  exit 2
}

executable=$(readlink -f -- "$executable")
resolved_launcher=$(readlink -f -- "$resolved_launcher")
mkdir -p -- "$output_root"
output_root=$(cd -- "$output_root" && pwd -P)

# Run the independent mathematical oracle before the executable.  A failure
# here means the host test environment cannot reproduce the formulas used by
# the strict log checker, so MPI evidence must not be interpreted.
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage8_forced_segment_reference.py" \
  >"$output_root/analytic-reference.log" 2>"$output_root/analytic-reference.err"

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

  printf 'START Stage-8 %s ranks=%s directory=%s\n' \
    "$name" "$ranks" "$case_dir"
  if ! (
    cd -- "$case_dir"
    OMP_NUM_THREADS=1 "$resolved_launcher" -n "$ranks" "$executable" "$@" \
      >run.log 2>error.log
  ); then
    printf 'FAIL Stage-8 %s ranks=%s; inspect %s\n' \
      "$name" "$ranks" "$case_dir" >&2
    return 1
  fi
}

for ranks in 1 4 12; do
  # Stage 8 must compose with every accepted correction layer.  Re-running
  # the Stage-4 mover/adaptation matrix with the same frozen binary is
  # essential: installing Stage-8 trajectory callbacks changes the production
  # first-order mover globally, so lower-stage face, MPI-transfer, continuation,
  # correction, collision, and adaptation contracts must still pass.  These
  # dependency replays do not replace the Stage-8 marker-specific fixture; they
  # prove that enabling it did not regress an existing transport path.
  run_case "$ranks" local-timestep --validate-local-timestep \
    --validation-particles "$validation_particles" \
    --validation-mover both
  run_case "$ranks" local-timestep-unit --validate-local-timestep \
    --validation-particles "$validation_particles" \
    --validation-mover both \
    --validation-unit-corrections
  run_case "$ranks" weight-adaptation --validate-weight-adaptation
  run_case "$ranks" correction-windows --validate-correction-windows
  run_case "$ranks" chemistry --validate-chemistry \
    --chemistry-validation-events "$chemistry_validation_events"
  run_case "$ranks" daughter-packets --validate-daughter-packets \
    --chemistry-validation-events "$chemistry_validation_events"
  run_case "$ranks" channel-importance --validate-channel-importance \
    --chemistry-validation-events "$chemistry_validation_events"
  run_case "$ranks" total-reaction-forcing \
    --validate-total-reaction-forcing \
    --chemistry-validation-events "$chemistry_validation_events"
done

# Preserve the accepted lower-stage parsers as the only authorities for their
# schemas.  The new Stage-8 checker owns homogeneous/piecewise segment
# placement, exponential survival, occurrence/channel/packet composition,
# capacity rollback, remaining-time motion, no-double-chemistry, and rank
# identity; duplicating lower-stage checks here would risk weaker criteria.
stage4_args=()
stage5_args=()
packet_args=(--expected-events "$chemistry_validation_events")
importance_args=(--expected-events "$chemistry_validation_events")
forcing_args=(--expected-events "$chemistry_validation_events")
for ranks in 1 4 12; do
  stage4_args+=(
    --local-log "$ranks" "$output_root/mpi${ranks}/local-timestep/run.log"
    --unit-local-log "$ranks" \
      "$output_root/mpi${ranks}/local-timestep-unit/run.log"
    --adaptation-log "$ranks" \
      "$output_root/mpi${ranks}/weight-adaptation/run.log"
    --chemistry-log "$ranks" "$output_root/mpi${ranks}/chemistry/run.log"
  )
  stage5_args+=(
    --window-log "$ranks" "$output_root/mpi${ranks}/correction-windows/run.log"
    --chemistry-log "$ranks" "$output_root/mpi${ranks}/chemistry/run.log"
  )
  packet_args+=(
    --packet-log "$ranks" "$output_root/mpi${ranks}/daughter-packets/run.log"
  )
  importance_args+=(
    --importance-log "$ranks" \
      "$output_root/mpi${ranks}/channel-importance/run.log"
    --diagnostics-root "$ranks" \
      "$output_root/mpi${ranks}/channel-importance"
  )
  forcing_args+=(
    --forcing-log "$ranks" \
      "$output_root/mpi${ranks}/total-reaction-forcing/run.log"
  )
done

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage4_integrated_logs.py" "${stage4_args[@]}"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage5_integrated_logs.py" "${stage5_args[@]}"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage6_daughter_packet_logs.py" "${packet_args[@]}"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage7_channel_importance_logs.py" \
  "${importance_args[@]}"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage8_total_reaction_forcing_logs.py" \
  "${forcing_args[@]}"

printf 'PASS Stage-8 total-reaction forcing runtime archive=%s\n' "$output_root"
