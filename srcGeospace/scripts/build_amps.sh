#!/bin/bash
#======================================================================================
# build_amps.sh -- build an `amps` executable of the srcGeospace application in a
# private, isolated AMPS root ("stage").
#
# The repository's own build products (AMPS/build, AMPS/amps, AMPS/share, Makefile.local,
# .last_input_info, ...) are never written, and nothing is written into srcGeospace.
# The stage directory contains:
#
#   symlinks (read only)   : src srcInterface utility species.input exosphere.pl ampsConfig.pl
#                            ampsConfigLib.pm Makefile Makefile.conf Makefile.test
#   private copies         : share input Makefile.def Makefile.local .amps.conf
#                            .general.conf (from the repository)
#   srcGeospace input      : srcGeospace/input/geospace.input -> <stage>/geospace.input,
#                            srcGeospace/input/srcGeospace.input -> <stage>/input/
#
# srcGeospace/input/srcGeospace.input selects the application with
# ProjectSourceDirectory=srcGeospace (correct for a build in the repository root).  In
# the stage copy that entry is replaced by the absolute path of --app, because
# ampsConfig.pl copies the directory with `cp -r` and a relative path would have to be a
# symbolic link, which `cp -r` copies as a link.
#
# Each invocation first REMOVES the whole previous stage (build tree, executable,
# configuration, logs, GPU library), i.e. every recompilation is a clean build.
#
# Usage:
#   srcGeospace/scripts/build_amps.sh [--app DIR] [--stage DIR] [-j N]
#                                     [--extra-cxx "FLAGS"] [--float]
#                                     [--gpu] [--gpu-target gpu|multicore]
#                                     [--gpu-opts "ARGS"] [--make-var VAR=VALUE]...
#
#   --app DIR       application sources (default: the srcGeospace directory of this script)
#   --stage DIR     private AMPS root (default: $GEOSPACE_RUNS/builds/<cpu|gpu|multicore>_<double|float>)
#   -j N            parallel make jobs (default: 16)
#   --extra-cxx F   appended to EXTRACOMPILEROPTIONS (e.g. -DGEOSPACE_PRECISION_FLOAT)
#   --float         Real=float engine build: same as --extra-cxx
#                   -DGEOSPACE_PRECISION_FLOAT (also selects the float device library
#                   with --gpu)
#   --gpu           build the device trajectory pool gpu/EngineDevice.cpp with nvc++
#                   (scripts/build_gpu_engine.sh) into <stage>/gpu_engine/
#                   libGeospaceEngine.so and link amps against it: make is called with
#                   GEOSPACE_DEVICE_OBJ= (the g++ EngineDevice.o is not built) and
#                   EXTRALINKEROPTIONS="<lib> -Wl,-rpath,<lib dir> -Wl,-rpath,<SDK lib>".
#                   No repository file is modified.  Nothing is executed on a GPU.
#   --gpu-target T  gpu (default, -acc=gpu) or multicore (-acc=multicore)
#   --gpu-opts A    extra arguments for build_gpu_engine.sh (e.g. "--arch cc80 -O 3")
#   --make-var V=X  additional make command-line variable (repeatable)
#
# Environment:
#   GEOSPACE_RUNS   root of builds and test results, outside srcGeospace
#                   (default: <repository>/../AMPS_runs/srcGeospace)
#
# Without --float/--gpu/--make-var the make invocation is exactly
#   make -j N amps EXTRACOMPILEROPTIONS="<--extra-cxx>"
#
# Result: <stage>/amps
#======================================================================================
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO=$(cd "$SRC_DIR/.." && pwd)
RUNS="${GEOSPACE_RUNS:-$REPO/../AMPS_runs/srcGeospace}"

APP="$SRC_DIR"
STAGE=""
JOBS=16
EXTRA_CXX=""
FLOAT=0
GPU=0
GPU_TARGET=gpu
GPU_OPTS=""
MAKE_VARS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --app) APP=$(cd "$2" && pwd); shift 2;;
    --stage) STAGE="$2"; shift 2;;
    -j) JOBS="$2"; shift 2;;
    --extra-cxx) EXTRA_CXX="$2"; shift 2;;
    --float) FLOAT=1; shift;;
    --gpu) GPU=1; shift;;
    --gpu-target) GPU_TARGET="$2"; shift 2;;
    --gpu-opts) GPU_OPTS="$2"; shift 2;;
    --make-var) MAKE_VARS+=("$2"); shift 2;;
    -h|--help) sed -n 2,58p "$0"; exit 0;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done
if [ -z "$STAGE" ]; then
  KIND=cpu; [ "$GPU" = 1 ] && KIND="$GPU_TARGET"
  PRECISION=double
  case " $EXTRA_CXX " in *GEOSPACE_PRECISION_FLOAT*) PRECISION=float;; esac
  [ "$FLOAT" = 1 ] && PRECISION=float
  STAGE="$RUNS/builds/${KIND}_${PRECISION}"
fi
mkdir -p "$STAGE"
STAGE=$(cd "$STAGE" && pwd)

case "$STAGE" in
  "$REPO"|"$REPO/build"|"$REPO/build/"*) echo "refusing to stage into the repository build tree" >&2; exit 2;;
  /|"$HOME"|"$SRC_DIR"|"$SRC_DIR/"*|"$APP"|"$APP/"*) echo "refusing to use $STAGE as a stage" >&2; exit 2;;
esac
case "$REPO/" in
  "$STAGE/"*) echo "refusing to use $STAGE as a stage (it contains the repository)" >&2; exit 2;;
esac

# Every build starts from an empty stage: the previous AMPS build tree, executable,
# configuration (.last_input_info, Makefile.def, deck), GPU library and logs are all
# removed, so nothing from an earlier build can be linked into the new executable.
# Only directories created by this script (marker file) or empty ones are wiped.
STAGE_MARKER=".srcGeospace_stage"
if [ -n "$(ls -A "$STAGE" 2>/dev/null)" ] && [ ! -e "$STAGE/$STAGE_MARKER" ]; then
  echo "refusing to remove $STAGE: not created by build_amps.sh (no $STAGE_MARKER)" >&2
  exit 2
fi
rm -rf "$STAGE"
mkdir -p "$STAGE"
touch "$STAGE/$STAGE_MARKER"

# Precision: --float, or the define given directly through --extra-cxx.
case " $EXTRA_CXX " in
  *GEOSPACE_PRECISION_FLOAT*) FLOAT=1;;
  *) [ "$FLOAT" = 1 ] && EXTRA_CXX="${EXTRA_CXX:+$EXTRA_CXX }-DGEOSPACE_PRECISION_FLOAT";;
esac

echo "[build_amps] repo  : $REPO"
echo "[build_amps] app   : $APP"
echo "[build_amps] stage : $STAGE"

# --- optional nvc++ device library (built first: fail before the long AMPS build) ---
if [ "$GPU" = 1 ]; then
  PREC=double
  [ "$FLOAT" = 1 ] && PREC=float
  # shellcheck disable=SC2086
  "$SCRIPT_DIR/build_gpu_engine.sh" --app "$APP" --out "$STAGE/gpu_engine" \
      --precision "$PREC" --target "$GPU_TARGET" $GPU_OPTS
  MAKE_VARS+=("GEOSPACE_DEVICE_OBJ=" "EXTRALINKEROPTIONS=$(cat "$STAGE/gpu_engine/link_flags.txt")")
  echo "[build_amps] gpu   : $GPU_TARGET/$PREC device pool $STAGE/gpu_engine/libGeospaceEngine.so"
fi

# --- read-only links ---------------------------------------------------------------
for f in src srcInterface utility species.input exosphere.pl ampsConfig.pl ampsConfigLib.pm \
         Makefile Makefile.conf Makefile.test; do
  ln -sfn "$REPO/$f" "$STAGE/$f"
done

# --- private copies ----------------------------------------------------------------
rm -rf "$STAGE/share" "$STAGE/input"
# share/ and input/ must be real directories in the stage: the build regenerates files in
# share/ and this script writes the deck into input/.  "cp -a" copies a symbolic link as
# a link, so if $REPO/share or $REPO/input is itself a link (e.g. a root assembled from
# links) the stage would write through it into the original directory.  Copy the
# directory contents (trailing "/.") and verify the result.
mkdir -p "$STAGE/share" "$STAGE/input"
cp -a "$REPO/share/." "$STAGE/share/"
cp -a "$REPO/input/." "$STAGE/input/"
for d in share input; do
  if [ -L "$STAGE/$d" ] || [ ! -d "$STAGE/$d" ]; then
    echo "[build_amps] ERROR: $STAGE/$d is not a private directory" >&2; exit 1
  fi
done
cp -f "$APP/input/geospace.input" "$STAGE/geospace.input"
cp -f "$APP/input/srcGeospace.input" "$STAGE/input/srcGeospace.input"
cp -f "$REPO/Makefile.local" "$STAGE/Makefile.local"
for f in .amps.conf .general.conf; do
  [ -e "$REPO/$f" ] && cp -f "$REPO/$f" "$STAGE/$f"
done
sed -e "s#^PTDIR *=.*#PTDIR    = $STAGE#" -e "s#^MYDIR *=.*#MYDIR    = $STAGE#" \
    "$REPO/Makefile.def" > "$STAGE/Makefile.def"

# The repository Makefile.local/.amps.conf describe whichever application was configured
# last in the repository root (Config.pl -application=...).  The stage always builds
# srcGeospace, so its application entries and the build options it needs are set here:
#   InputFileAMPS/APPLICATION  the stage deck geospace.input
#   INTERFACE=on               Fortran Geopack/IGRF/Tsyganenko interface (src/interface)
#   OPENMP=on                  -fopenmp (Mode3D density loops; MPI compilation mode kept)
# Site settings (SPICE, SPICEKERNELS, compilers) still come from the repository files.
set_conf() {   # set_conf FILE KEY VALUE: replace every KEY= line, or append one
  if grep -q "^$2=" "$1" 2>/dev/null; then sed -i "s#^$2=.*#$2=$3#" "$1"; else echo "$2=$3" >> "$1"; fi
}
set_conf "$STAGE/Makefile.local" InputFileAMPS geospace.input
set_conf "$STAGE/Makefile.local" APPLICATION geospace
set_conf "$STAGE/Makefile.local" INTERFACE on
set_conf "$STAGE/Makefile.local" OPENMP on
touch "$STAGE/.amps.conf"
set_conf "$STAGE/.amps.conf" InputFileAMPS geospace.input
set_conf "$STAGE/.amps.conf" APPLICATION geospace
if grep -q '^SPICE=nospice' "$STAGE/Makefile.local" || ! grep -q '^SPICE=' "$STAGE/Makefile.local"; then
  echo "[build_amps] WARNING: SPICE is not configured (repository Makefile.local); runs that use" >&2
  echo "[build_amps]          driver files or SPICE frames will fail (see srcGeospace/BUILD.md)" >&2
fi

# Point the deck at the application sources (absolute path, see the header).  The deck
# name is taken from the first #include of geospace.input.
DECK=$(sed -n 's/^#include[ \t]\+\(input\/[^ \t]*\).*/\1/p' "$STAGE/geospace.input" | head -1)
[ -n "$DECK" ] || { echo "cannot find the #include deck in geospace.input" >&2; exit 2; }
grep -q '^ProjectSourceDirectory=' "$STAGE/$DECK" || { echo "no ProjectSourceDirectory in $DECK" >&2; exit 2; }
sed -i "s#^ProjectSourceDirectory=.*#ProjectSourceDirectory=$APP#" "$STAGE/$DECK"
echo "[build_amps] deck  : $DECK -> ProjectSourceDirectory=$APP"

# (The stage was emptied above, so the AMPS source tree build/ is created fresh from
# --app by ampsConfig.pl.)

cd "$STAGE"
if [ ${#MAKE_VARS[@]} -gt 0 ]; then
  printf '[build_amps] make  : %s\n' "${MAKE_VARS[@]}"
fi
make -j "$JOBS" amps EXTRACOMPILEROPTIONS="$EXTRA_CXX" ${MAKE_VARS[@]+"${MAKE_VARS[@]}"} \
    > "$STAGE/build.log" 2>&1 || {
  echo "[build_amps] FAILED; see $STAGE/build.log" >&2; tail -40 "$STAGE/build.log" >&2; exit 1; }
ls -la "$STAGE/amps"
if [ "$GPU" = 1 ]; then
  # The device pool must come from the nvc++ library, not from a g++ object.
  if nm --defined-only "$STAGE/amps" 2>/dev/null | grep ' T eg_dev_create$' >/dev/null; then
    echo "[build_amps] ERROR: amps defines eg_dev_create itself (g++ EngineDevice.o linked)" >&2
    exit 1
  fi
  readelf -d "$STAGE/amps" | grep -E 'libGeospaceEngine|RUNPATH|RPATH' | sed 's/^/[build_amps]   /' || true
fi
echo "[build_amps] done: $STAGE/amps"
