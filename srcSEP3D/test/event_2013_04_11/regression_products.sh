#!/bin/sh
# Assemble the regression products of one srcSEP3D 2013-04-11 run.
#
# Usage: regression_products.sh RUN_DIR OUTPUT_DIR
#
# The run publishes under RUN_DIR/PT/plots/sep3d (output.directory of the
# deck).  This script writes OUTPUT_DIR/{observers,radial_energy,ledger}.dat,
# each a whitespace-separated copy (commas become spaces, so DiffNum.pl
# compares every number) of srcSEP3D output tables, in a fixed order:
#   observers      every observer row of every publication's spacecraft.csv,
#                  prefixed with the publication sequence number (the
#                  sequences in increasing order; rows as published)
#   radial_energy  <prefix>.radial-energy.csv: represented particles and
#                  macroparticles per heliocentric-radius and energy bin at
#                  every publication (observer-flux-tube runs)
#   ledger         the last publication's shocks.csv: injected, escaped,
#                  absorbed (inner boundary and returns to the shock) and
#                  failed macroparticles per step
# The Table Check compares each with its reference
# output/test_SEP3D/SEP3D--event-2013-04-11--<mover>[-long]--<label>.ref_np<N>.gz
# using share/Scripts/DiffNum.pl; make_references.sh builds the references
# from the same files.  A missing product fails the script (non-zero exit),
# which the Table turns into a non-empty .diff.
set -eu
run=$1
out=$2
plots="$run/PT/plots/sep3d"
mkdir -p "$out"
# Byte-order sorting and C number formatting, independent of the locale of
# the machine that runs the check.
LC_ALL=C
export LC_ALL

# Publication directories are <prefix>.<zero-padded sequence number>, e.g.
# event3d-nightly.00000001; the zero padding makes the byte order the time
# order.  The radial-energy histogram (<prefix>.radial-energy.csv) and the
# other files at this level do not match the pattern.  `|| true` lets an
# empty listing reach the explicit error below instead of stopping the
# script under set -e.
publications=$(ls -1d "$plots"/*.[0-9]* 2>/dev/null | sort || true)
if [ -z "$publications" ]; then
  echo "regression_products: no publication under $plots" >&2
  exit 1
fi

# observers: drop each spacecraft.csv header line and prefix every row with
# the publication sequence number, so a row is identified by its time and
# its position in the published order.
: > "$out/observers.dat"
for dir in $publications; do
  seq=${dir##*.}
  tail -n +2 "$dir/spacecraft.csv" | sed "s/^/$seq,/" | tr ',' ' ' >> "$out/observers.dat"
done

# radial_energy: the run's single histogram file, which holds the rows of
# every publication (its first column is the simulation time).
histogram=$(ls -1 "$plots"/*.radial-energy.csv 2>/dev/null | head -1 || true)
if [ -z "$histogram" ]; then
  echo "regression_products: no radial-energy histogram under $plots" >&2
  exit 1
fi
tr ',' ' ' < "$histogram" > "$out/radial_energy.dat"

# ledger: shocks.csv holds one row per completed step and species
# (step, species, injected_count, escaped_count, absorbed_count,
# failed_count, shock_crossings_count) for every step so far, so the last
# publication's copy covers the whole run.
last=$(echo "$publications" | tail -n 1)
if [ ! -f "$last/shocks.csv" ]; then
  echo "regression_products: no ledger in $last" >&2
  exit 1
fi
tr ',' ' ' < "$last/shocks.csv" > "$out/ledger.dat"
