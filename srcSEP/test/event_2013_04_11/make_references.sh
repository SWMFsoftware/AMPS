#!/bin/sh
# Make the DiffNum references of one 2013-04-11 real-event Table test from a
# finished run of that test.
#
# Usage: make_references.sh RUN_DIR TEST_NAME REFERENCE_DIR
#   RUN_DIR        run_test_<TEST_NAME> after the test's run stage
#   TEST_NAME      e.g. SEP--event-2013-04-11--Parker
#   REFERENCE_DIR  the checkout's output/test_SEP3D (the SWMF_data tree)
#
# Writes REFERENCE_DIR/<TEST_NAME>--<label>.ref_np<N>.gz for the labels
# density, flux, energy and pitch_angle, where N is the number of MPI ranks of
# the run (its PT/thread* files) and each file is the product assembled by
# regression_products.sh, so the Table Check compares like with like.
# gzip -n leaves the name and time stamp out of the archive, so a reference
# made twice from the same run is byte-identical.
#
# Regenerate references only for an intended change of the results, on the
# reference machine and compiler, and commit them to SWMF_data with the AMPS
# commit and the run's command line in the message.
set -eu
run=$1
name=$2
refdir=$3
here=$(cd "$(dirname "$0")" && pwd)
np=$(ls -d "$run"/PT/thread* 2>/dev/null | wc -l | tr -d ' ')
if [ "$np" -eq 0 ]; then
  echo "make_references: $run has no PT/thread* files (did the run stage finish?)" >&2
  exit 1
fi
products=$(mktemp -d "${TMPDIR:-/tmp}/sep-event-products.XXXXXX")
trap 'rm -rf "$products"' EXIT
sh "$here/regression_products.sh" "$run" "$products"
mkdir -p "$refdir"
for label in density flux energy pitch_angle; do
  gzip -9 -n -c "$products/$label.dat" > "$refdir/$name--$label.ref_np$np.gz"
  echo "wrote $refdir/$name--$label.ref_np$np.gz"
done
