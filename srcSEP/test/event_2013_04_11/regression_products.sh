#!/bin/sh
# Assemble the regression products of one 2013-04-11 real-event run.
#
# Usage: regression_products.sh RUN_DIR OUTPUT_DIR
#
# Writes OUTPUT_DIR/{density,flux,energy,pitch_angle}.dat, each the plain
# concatenation (no reformatting) of native AMPS output tables written under
# RUN_DIR/PT/plots, in a fixed lexicographic order:
#   density      every observer's sample.density history
#   flux         every observer's sample.flux history
#   energy       the final EnergyRSample table (energy versus radius)
#   pitch_angle  the pitch-angle distributions of the broadband observers:
#                the inner regression observers nightly-r*au (event.nightly.in
#                only), then the SOHO observer soho-spectrum
# The same four labels are used by the neighbouring SEP--Parker_spiral--*
# Table tests.  The Table Check compares each file with its reference
# output/test_SEP3D/SEP--event-2013-04-11--<mover>[--long]--<label>.ref_np<N>.gz
# using share/Scripts/DiffNum.pl, and the reference package is made from the
# same files, so both sides always use one ordering.  A missing product makes
# this script fail (non-zero exit), which the Table turns into a non-empty
# .diff.
set -eu
run=$1
out=$2
plots="$run/PT/plots"
mkdir -p "$out"
LC_ALL=C
export LC_ALL

collect() {
  label=$1
  shift
  : > "$out/$label.dat"
  found=0
  for file in "$@"; do
    if [ -f "$file" ]; then
      cat "$file" >> "$out/$label.dat"
      found=$((found + 1))
    fi
  done
  if [ "$found" -eq 0 ]; then
    echo "regression_products: no $label product under $plots" >&2
    exit 1
  fi
}

collect density $(ls -1 "$plots"/sample.density/observer=*.dat | sort)
collect flux $(ls -1 "$plots"/sample.flux/observer=*.dat | sort)
last=$(ls -1 "$plots"/EnergyRSample/cnt=*.dat | sed 's/.*cnt=\([0-9]*\)\.dat/\1 &/' | sort -n | tail -1 | cut -d' ' -f2)
collect energy "$last"
# Unmatched patterns stay literal and are skipped by collect's file test; the
# shell expands each pattern in sorted order (LC_ALL=C above).
collect pitch_angle "$plots"/sample.pitch_angle_distribution/observer=nightly-r*.dat \
  "$plots"/sample.pitch_angle_distribution/observer=soho-spectrum.*.dat
