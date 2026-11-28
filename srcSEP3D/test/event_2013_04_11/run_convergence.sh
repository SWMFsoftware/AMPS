#!/bin/sh
# Collection-radius convergence study of the srcSEP3D 2013-04-11 test
# (README.md, section 6): one long-tier run per collection radius, then
# convergence_report.py.
#
# Usage: run_convergence.sh AMPS_BINARY MOVER WORK_DIR [RADIUS_TAG ...]
#   AMPS_BINARY  an amps built with ./Config.pl -application=sep3d
#   MOVER        parker, fte-dmumu or fte-mfp
#   WORK_DIR     gets one run directory convergence-<RADIUS_TAG>-<MOVER> per radius
#   RADIUS_TAG   r0p005au r0p01au r0p02au r0p04au by default
# Environment: MPIRUN (default "mpirun -np 16", the long-tier rank count). A run directory that already
# holds a finished run (its run.status says exit=0) is reused, so the study
# can be resumed or extended.
#
# Each run directory receives the same inputs as a Table run directory: all
# event3d.*.in decks, the body ephemerides, and the shared 1-D event file
# with its WSO asset.  amps is a symbolic link to AMPS_BINARY.  The run's
# stdout/stderr go to run.log/run.err, its exit code to run.status.  After the
# last run the script runs convergence_report.py on all radii; its exit
# status (0 selected, 1 none, 2 bad input) is the script's exit status.  A
# failed run stops the script with exit 1 before the report.
#
# Example: MPIRUN="mpirun -np 16" sh run_convergence.sh ../../../amps fte-mfp /scratch/conv
set -eu
# here: this directory (decks, ephemerides, report); e1: the 1-D event
# directory, which holds event.event (#included by the decks) and its WSO
# asset.
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
e1=$(CDPATH= cd -- "$here/../../../srcSEP/test/event_2013_04_11" && pwd)
# The binary must be an executable file; resolve it to an absolute path so
# the link in each run directory works from there.
[ -x "$1" ] && [ -f "$1" ] || { echo "no executable $1" >&2; exit 2; }
amps=$(CDPATH= cd -- "$(dirname -- "$1")" && pwd)/$(basename -- "$1")
mover=$2
work=$3
shift 3
# Default radius tags: the four radii of event3d_settings.json "convergence".
[ $# -gt 0 ] || set -- r0p005au r0p01au r0p02au r0p04au
case $mover in parker|fte-dmumu|fte-mfp) ;; *) echo "unknown mover $mover" >&2; exit 2;; esac
[ -x "$amps" ] || { echo "no executable $amps" >&2; exit 2; }
mpirun=${MPIRUN:-mpirun -np 16}
mkdir -p "$work"
# pairs: the RADIUS_AU=RUN_DIR arguments of convergence_report.py.
pairs=""
for tag in "$@"; do
  deck="event3d.convergence-$tag.$mover.in"
  [ -f "$here/$deck" ] || { echo "no deck $here/$deck" >&2; exit 2; }
  dir="$work/convergence-$tag-$mover"
  # Reuse a finished run; otherwise start this radius from a clean directory
  # (an unfinished or failed earlier attempt is discarded).
  if ! grep -qs '^exit=0$' "$dir/run.status"; then
    rm -rf "$dir"
    mkdir -p "$dir/assets" "$dir/ephemeris"
    cp "$here"/event3d.*.in "$dir"/
    cp "$here"/ephemeris/*.csv "$dir"/ephemeris/
    cp "$e1/event.event" "$dir"/
    cp "$e1/assets/wso_cr2135_model_frame.harmonics.csv" "$dir"/assets/
    ln -s "$amps" "$dir/amps"
    echo "convergence: $tag $mover started $(date)"
    # MPIRUN is expanded unquoted so that "mpirun -np 16" splits into words.
    status=0
    (cd "$dir" && $mpirun ./amps -input "$deck" > run.log 2> run.err) || status=$?
    echo "exit=$status" > "$dir/run.status"
    echo "convergence: $tag $mover exit=$status $(date)"
    [ "$status" -eq 0 ] || exit 1
  fi
  # Tag -> radius in au: r0p005au -> 0.005.
  radius=$(echo "$tag" | sed 's/^r//; s/au$//; s/p/./')
  pairs="$pairs $radius=$dir"
done
# shellcheck disable=SC2086
python3 "$here/convergence_report.py" $pairs
