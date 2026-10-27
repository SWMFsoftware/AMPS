#!/bin/bash
#======================================================================================
# build_gpu_engine.sh -- build the srcGeospace device trajectory pool
# (gpu/EngineDevice.cpp) with the NVIDIA HPC SDK compiler as a shared library that
# the g++-built AMPS executable links against.
#
# Only the device pool is compiled by nvc++.  Everything else (AMPS, the application
# drivers, BatchTracer, fibers, the host pool) stays g++/OpenMP.  The two sides talk
# through the extern "C" API of gpu/EngineDeviceAPI.h, which carries only POD types,
# so no C++ library object crosses the compiler boundary.
#
# Usage:
#   srcGeospace/scripts/build_gpu_engine.sh [--app DIR] [--out DIR]
#        [--precision double|float] [--target gpu|multicore] [--arch ccNN]
#        [--cuda X.Y] [--rdc] [--fma] [-O N] [--nvhpc-env FILE]
#        [--extra-nvcxx "FLAGS"]
#
#   --app DIR        application sources (default: the srcGeospace directory of this script)
#   --out DIR        output directory (default:
#                    $GEOSPACE_RUNS/builds/gpu_engine/<target>-<precision>)
#   --precision P    double (default; Real=double, SI units) or float
#                    (-DGEOSPACE_PRECISION_FLOAT: Real=float, normalized geospace
#                    units).  The AMPS side must be built with the same setting
#                    (build_amps.sh --gpu --float does this).
#   --target T       gpu (default: -acc=gpu) or multicore (-acc=multicore)
#   --arch ccNN      GPU compute capability (default: cc86 = RTX A6000)
#   --cuda X.Y       CUDA version bundled with the SDK (default: 12.4)
#   --rdc            relocatable device code (default: -gpu=nordc; the whole device
#                    pool is one translation unit, so no device link is needed)
#   --fma            allow nvc++ to contract a*b+c into fused multiply-adds.  Default is
#                    -gpu=nofma (gpu) / -Mnofma (multicore): nvc++ 24.5 contracts at -O2
#                    on both targets, while the g++ AMPS build (x86-64 baseline, no
#                    -march) never does, so nofma keeps the double build comparable with
#                    the legacy tracer (AMPS_EARTH_ENGINE_VERIFY=1).  See doc/BUILD.md.
#   -O N             optimization level (default: 2)
#   --nvhpc-env F    script that puts nvc++ on PATH
#                    (default: $NVHPC_ENV or /data/vtenishe/openacc/env.sh)
#   --extra-nvcxx F  additional nvc++ flags
#
# Result:
#   <out>/libGeospaceEngine.so      SONAME libGeospaceEngine.so, RUNPATH = SDK libraries
#   <out>/EngineDevice.o            object compiled with -fPIC
#   <out>/compile.log, link.log     nvc++ output (-Minfo=accel)
#   <out>/link_flags.txt            linker options for the AMPS executable, e.g.
#                                   make amps GEOSPACE_DEVICE_OBJ= \
#                                        EXTRALINKEROPTIONS="$(cat <out>/link_flags.txt)"
#
# This script only compiles and links; it never runs code on a GPU.
#======================================================================================
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO=$(cd "$SRC_DIR/.." && pwd)
RUNS="${GEOSPACE_RUNS:-$REPO/../AMPS_runs/srcGeospace}"

APP="$SRC_DIR"
OUT=""
PRECISION=double
TARGET=gpu
ARCH=cc86
CUDA=12.4
RDC=0
FMA=0
OPT=2
NVHPC_ENV_FILE="${NVHPC_ENV:-/data/vtenishe/openacc/env.sh}"
EXTRA_NVCXX=""
while [ $# -gt 0 ]; do
  case "$1" in
    --app) APP=$(cd "$2" && pwd); shift 2;;
    --out) OUT="$2"; shift 2;;
    --precision) PRECISION="$2"; shift 2;;
    --target) TARGET="$2"; shift 2;;
    --arch) ARCH="$2"; shift 2;;
    --cuda) CUDA="$2"; shift 2;;
    --rdc) RDC=1; shift;;
    --fma) FMA=1; shift;;
    -O) OPT="$2"; shift 2;;
    --nvhpc-env) NVHPC_ENV_FILE="$2"; shift 2;;
    --extra-nvcxx) EXTRA_NVCXX="$2"; shift 2;;
    -h|--help) sed -n 2,51p "$0"; exit 0;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done

case "$PRECISION" in
  double) PREC_FLAGS="";;
  float)  PREC_FLAGS="-DGEOSPACE_PRECISION_FLOAT";;
  *) echo "--precision must be double or float" >&2; exit 2;;
esac
case "$TARGET" in
  gpu)
    GPU_OPTS="-gpu=${ARCH},cuda${CUDA}"
    [ "$RDC" = 1 ] || GPU_OPTS="${GPU_OPTS},nordc"
    [ "$FMA" = 1 ] || GPU_OPTS="${GPU_OPTS},nofma"
    ACC_FLAGS="-acc=gpu $GPU_OPTS";;
  multicore)
    ACC_FLAGS="-acc=multicore"
    [ "$FMA" = 1 ] || ACC_FLAGS="$ACC_FLAGS -Mnofma";;
  *) echo "--target must be gpu or multicore" >&2; exit 2;;
esac

SRC="$APP/gpu/EngineDevice.cpp"
[ -f "$SRC" ] || { echo "missing $SRC" >&2; exit 2; }
[ -n "$OUT" ] || OUT="$RUNS/builds/gpu_engine/${TARGET}-${PRECISION}"
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)

if ! command -v nvc++ >/dev/null 2>&1; then
  [ -f "$NVHPC_ENV_FILE" ] || { echo "nvc++ not on PATH and $NVHPC_ENV_FILE missing" >&2; exit 2; }
  # shellcheck disable=SC1090
  source "$NVHPC_ENV_FILE"
fi
NVCXX=$(command -v nvc++)
NVHPC_BIN=$(dirname "$NVCXX")
NVHPC_LIB=$(cd "$NVHPC_BIN/../lib" && pwd)

CXXFLAGS="-std=c++17 -O${OPT} -fPIC ${ACC_FLAGS} ${PREC_FLAGS} -Minfo=accel ${EXTRA_NVCXX}"
echo "[build_gpu_engine] nvc++   : $NVCXX ($("$NVCXX" --version | sed -n 2p))"
echo "[build_gpu_engine] source  : $SRC"
echo "[build_gpu_engine] flags   : $CXXFLAGS"
echo "[build_gpu_engine] output  : $OUT"

rm -f "$OUT/EngineDevice.o" "$OUT/libGeospaceEngine.so" "$OUT/link_flags.txt"
# shellcheck disable=SC2086
"$NVCXX" $CXXFLAGS -c "$SRC" -o "$OUT/EngineDevice.o" > "$OUT/compile.log" 2>&1 || {
  echo "[build_gpu_engine] compile FAILED; see $OUT/compile.log" >&2
  grep -E 'fatal|NVC\+\+-[SF]-|nvlink error|[Uu]ndefined reference|Unresolved|error:' \
      "$OUT/compile.log" | head -20 >&2
  exit 1; }

# The SDK driver adds a RUNPATH to its own lib directories, so the library finds the
# OpenACC/CUDA runtime (libacchost, libaccdevice, libnvomp, libnvc, ...) by itself.
# shellcheck disable=SC2086
"$NVCXX" -shared -fPIC ${ACC_FLAGS} -Wl,-soname,libGeospaceEngine.so \
    -o "$OUT/libGeospaceEngine.so" "$OUT/EngineDevice.o" > "$OUT/link.log" 2>&1 || {
  echo "[build_gpu_engine] link FAILED; see $OUT/link.log" >&2
  head -20 "$OUT/link.log" >&2
  exit 1; }

NEXPORT=$(nm -D --defined-only "$OUT/libGeospaceEngine.so" | grep -c ' T eg_dev_' || true)
[ "$NEXPORT" -ge 14 ] || { echo "[build_gpu_engine] expected 14 eg_dev_* exports, found $NEXPORT" >&2; exit 1; }

# Linker options for the AMPS executable: the library by absolute path plus a RUNPATH
# entry for its directory and for the SDK libraries.
echo "$OUT/libGeospaceEngine.so -Wl,-rpath,$OUT -Wl,-rpath,$NVHPC_LIB" > "$OUT/link_flags.txt"
{
  echo "precision=$PRECISION"
  echo "target=$TARGET"
  echo "flags=$CXXFLAGS"
  echo "nvcxx=$("$NVCXX" --version | sed -n 2p)"
  echo "source=$SRC"
  echo "date=$(date -Is)"
} > "$OUT/build_info.txt"

ls -la "$OUT/libGeospaceEngine.so"
readelf -d "$OUT/libGeospaceEngine.so" | grep -E 'NEEDED|RUNPATH|SONAME' | sed 's/^/[build_gpu_engine]   /'
echo "[build_gpu_engine] exports : $NEXPORT eg_dev_* functions"
echo "[build_gpu_engine] link    : $(cat "$OUT/link_flags.txt")"
echo "[build_gpu_engine] done"
