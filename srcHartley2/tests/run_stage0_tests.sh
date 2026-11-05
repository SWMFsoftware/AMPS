#!/usr/bin/env bash
# Run every focused test required by srcHartley2/CHEMISTRY.md Stage 0.  The long Halley
# archive remains external to source control; when it is present at the
# documented validation_runs path (or HARTLEY2_STAGE0_BASELINE_ROOT is set),
# the Python suite also verifies every preserved checksum and C2 acceptance
# value.  Otherwise that one archive test is reported as skipped.  When
# HARTLEY2_STAGE0_PRODUCTION_ROOT names an instrumented production run, the
# suite also validates its exact cell measures, CSV/JSON accounting
# invariants, and an observed represented C2 carrier.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
stage0_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage0.XXXXXX")
trap 'rm -rf -- "$stage0_tmp"' EXIT

# The all-stage runner may select a non-default MPI launcher.  Keep the
# historical command name when no override is supplied so direct invocations
# retain their original behavior.  Resolve it before compiling anything;
# otherwise a missing launcher would leave a misleading partial test run.
stage0_mpi_launcher=${HARTLEY2_MPIEXEC:-mpiexec}
stage0_mpi_launcher=$(command -v -- "$stage0_mpi_launcher" || true)
if [[ -z "$stage0_mpi_launcher" ]]; then
  echo "Stage-0 MPI launcher is not on PATH: ${HARTLEY2_MPIEXEC:-mpiexec}" >&2
  exit 2
fi

mpicxx -std=c++11 -O2 -Wall -Wextra -DOMPI_SKIP_MPICXX=1 \
  -I"$source_root/src" \
  "$script_dir/test_stage0_chemistry_diagnostics.cpp" \
  -o "$stage0_tmp/test_stage0_chemistry_diagnostics"

reference_signature=""
for ranks in 1 4 12; do
  output=$("$stage0_mpi_launcher" -n "$ranks" --oversubscribe \
    "$stage0_tmp/test_stage0_chemistry_diagnostics")
  printf '%s\n' "$output"
  signature=$(printf '%s\n' "$output" | sed -n 's/^stage0-mpi-signature /stage0-mpi-signature /p')
  if [[ -z "$signature" ]]; then
    echo "missing MPI diagnostic signature for $ranks ranks" >&2
    exit 1
  fi
  if [[ -z "$reference_signature" ]]; then
    reference_signature=$signature
  elif [[ "$signature" != "$reference_signature" ]]; then
    echo "MPI diagnostic totals changed with rank count" >&2
    echo "one rank:  $reference_signature" >&2
    echo "$ranks ranks: $signature" >&2
    exit 1
  fi
done

PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_stage0_minor_species.py"

if [[ -n "${HARTLEY2_STAGE0_PRODUCTION_ROOT:-}" ]]; then
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage0_production_report.py" \
    "$HARTLEY2_STAGE0_PRODUCTION_ROOT" --require-observed-c2
else
  echo "SKIP Stage-0 production report: set HARTLEY2_STAGE0_PRODUCTION_ROOT"
fi
