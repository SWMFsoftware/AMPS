#!/usr/bin/env bash
# Execute the initialized-mesh Stage-5 controller and exact-chemistry
# validators at all accepted MPI decompositions.  The executable is supplied
# by the caller: this runner never reconfigures AMPS, removes build/, or edits
# a science input, so archived evidence remains tied to one frozen binary.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: run_stage5_integrated_tests.sh \
  --executable /absolute/path/to/correction-on-amps \
  --output-root /new/path/to/stage5-validation \
  [--mpiexec mpiexec] [--chemistry-validation-events 200000]

Runs correction-window and exact chemistry validators at exactly 1, 4, and
12 MPI ranks. The executable must be a DIM1 radial build with individual
corrections and ParticleWeightCorrectionWindowMode=ON. The output root must
not already exist; failed or partial archives are never overwritten.
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
      printf 'unknown Stage-5 runner option: %s\n' "$1" >&2
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
[[ -f "$executable" && -x "$executable" ]] || {
  printf 'Stage-5 executable is not executable: %s\n' "$executable" >&2
  exit 2
}
resolved_launcher=$(command -v -- "$mpi_launcher" || true)
[[ -n "$resolved_launcher" ]] || {
  printf 'MPI launcher is not on PATH: %s\n' "$mpi_launcher" >&2
  exit 2
}
[[ ! -e "$output_root" ]] || {
  printf 'Stage-5 output root already exists: %s\n' "$output_root" >&2
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

  printf 'START Stage-5 %s ranks=%s directory=%s\n' \
    "$name" "$ranks" "$case_dir"
  if ! (
    cd -- "$case_dir"
    OMP_NUM_THREADS=1 "$resolved_launcher" -n "$ranks" "$executable" "$@" \
      >run.log 2>error.log
  ); then
    printf 'FAIL Stage-5 %s ranks=%s; inspect %s\n' \
      "$name" "$ranks" "$case_dir" >&2
    return 1
  fi
}

for ranks in 1 4 12; do
  # Individual corrections and the current owner-block systematic residual
  # state are intentionally unsupported together.  Override only the residual
  # allocation policy; the reaction network and physical RNG seed are frozen.
  run_case "$ranks" correction-windows \
    --chemistry-residual-sampling bernoulli \
    --validate-correction-windows
  run_case "$ranks" chemistry --validate-chemistry \
    --chemistry-residual-sampling bernoulli \
    --chemistry-validation-events "$chemistry_validation_events"
done

python_args=()
for ranks in 1 4 12; do
  python_args+=(
    --window-log "$ranks" "$output_root/mpi${ranks}/correction-windows/run.log"
    --chemistry-log "$ranks" "$output_root/mpi${ranks}/chemistry/run.log"
  )
done
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage5_integrated_logs.py" "${python_args[@]}"

printf 'PASS Stage-5 correction-window runtime archive=%s\n' "$output_root"
