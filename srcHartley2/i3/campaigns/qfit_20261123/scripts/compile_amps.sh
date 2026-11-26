#!/bin/bash
# Build the AMPS executable of the I3 campaign qfit_20261123 in a separate tree (not in the checkout that
# holds this script), from a given commit and the campaign's build input ../i3.conf.
#
#   compile_amps.sh --dest DIR [--commit REV] [--repo URL|PATH] [--jobs N] [--compiler LIST]
#                   [--share-repo URL|PATH] [--util-repo URL|PATH] [--share-commit REV] [--test-mode]
#
#   --dest DIR        new directory for the build tree (must not exist)
#   --commit REV      AMPS commit to build (default 848e9c5: the exact code of the campaign executable,
#                     c42b555 + the DIM1 domain-boundary fix; see ../README.md)
#   --repo URL|PATH   AMPS repository to clone (default: the repository that contains this script)
#   --jobs N          parallel make jobs (default 8)
#   --compiler LIST   Config.pl -compiler value (default gfortran,gcc_mpicc, as in the campaign)
#   --share-repo, --util-repo
#                     SWMF share and util repositories (default git@github.com:SWMFsoftware/share.git and
#                     .../util.git; a local clone path also works)
#   --share-commit REV  share commit (default c5f8beef, the campaign build's; '' = the repository's default
#                     branch)
#   --test-mode       configure with -amps-test=on (nightly-test build; NOT the production build of the
#                     campaign, which is -amps-test=off)
#
# Steps: clone, check out REV, clone share and util, copy ../i3.conf to input/i3.conf,
# Config.pl -install and -application=i3, rm -rf -- build, make EXTRACOMPILEROPTIONS=-O3, then check the
# compiled block layout (_BLOCK_CELLS_X_ 12, _MAX_REFINMENT_LEVEL_ 30) and write DIR/build_provenance.txt
# and DIR/amps.build.sha256. Nothing is written outside DIR.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CAMPAIGN=$(dirname "$HERE")
DEST="" COMMIT=848e9c5 REPO="" JOBS=8 COMPILER=gfortran,gcc_mpicc TESTMODE=off
SHARE_REPO=git@github.com:SWMFsoftware/share.git UTIL_REPO=git@github.com:SWMFsoftware/util.git
SHARE_COMMIT=c5f8beef
while [[ $# -gt 0 ]]; do
  case $1 in
    --dest) DEST=$2; shift 2 ;;
    --commit) COMMIT=$2; shift 2 ;;
    --repo) REPO=$2; shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;;
    --compiler) COMPILER=$2; shift 2 ;;
    --share-repo) SHARE_REPO=$2; shift 2 ;;
    --util-repo) UTIL_REPO=$2; shift 2 ;;
    --share-commit) SHARE_COMMIT=$2; shift 2 ;;
    --test-mode) TESTMODE=on; shift ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option $1 (see --help)"; exit 2 ;;
  esac
done
[[ -n $DEST ]] || { echo "--dest DIR is required"; exit 2; }
[[ ! -e $DEST ]] || { echo "refusing to use existing $DEST"; exit 2; }
[[ -n $REPO ]] || REPO=$(git -C "$HERE" rev-parse --show-toplevel)
( cd "$CAMPAIGN" && sha256sum -c --quiet SHA256SUMS ) || { echo "campaign files differ from SHA256SUMS"; exit 3; }

echo "## clone $REPO at $COMMIT into $DEST"
git clone --quiet "$REPO" "$DEST"
DEST=$(cd "$DEST" && pwd)
git -C "$DEST" -c advice.detachedHead=false checkout --quiet --detach "$COMMIT"
git clone --quiet "$SHARE_REPO" "$DEST/share"
[[ -z $SHARE_COMMIT ]] || git -C "$DEST/share" -c advice.detachedHead=false checkout --quiet --detach "$SHARE_COMMIT"
git clone --quiet "$UTIL_REPO" "$DEST/util"

echo "## configure (build input: $CAMPAIGN/i3.conf, -amps-test=$TESTMODE)"
cp "$CAMPAIGN/i3.conf" "$DEST/input/i3.conf"
cd "$DEST"
./Config.pl -install -compiler="$COMPILER" > install.log 2>&1
./Config.pl -application=i3 -spice-path=nospice -spice-kernels=nospice -amps-test=$TESTMODE > configure.log 2>&1
cmp input/i3.conf "$CAMPAIGN/i3.conf"
[[ ! -e i3.conf ]] || cmp i3.conf "$CAMPAIGN/i3.conf"
grep -q "^TESTMODE=$TESTMODE" .amps.conf || { echo ".amps.conf does not record TESTMODE=$TESTMODE"; exit 4; }

echo "## compile (make -j$JOBS EXTRACOMPILEROPTIONS=-O3, log: $DEST/compile.log)"
rm -rf -- build
make -j"$JOBS" EXTRACOMPILEROPTIONS=-O3 > compile.log 2>&1 || { tail -30 compile.log; exit 5; }
[[ -x amps ]] || { echo "no amps executable"; exit 5; }
grep -q '^#define _BLOCK_CELLS_X_ 12$' build/meshAMR/meshAMRdef.h || { echo "_BLOCK_CELLS_X_ is not 12"; exit 6; }
grep -q '^#define _MAX_REFINMENT_LEVEL_ 30$' build/meshAMR/meshAMRdef.h || { echo "_MAX_REFINMENT_LEVEL_ is not 30"; exit 6; }

sha256sum amps | cut -d' ' -f1 > amps.build.sha256
{
  echo "date=$(date -Is) host=$(hostname)"
  echo "amps_commit=$(git rev-parse HEAD) requested=$COMMIT repo=$REPO"
  echo "share_commit=$(git -C share rev-parse HEAD) util_commit=$(git -C util rev-parse HEAD)"
  echo "config: Config.pl -install -compiler=$COMPILER; Config.pl -application=i3 -spice-path=nospice" \
       "-spice-kernels=nospice -amps-test=$TESTMODE; make -j$JOBS EXTRACOMPILEROPTIONS=-O3"
  echo "build_input=$CAMPAIGN/i3.conf sha256=$(sha256sum "$CAMPAIGN/i3.conf" | cut -d' ' -f1)"
  echo "compilers: $( (mpicxx --version 2>/dev/null || true) | head -1) | $( (gfortran --version 2>/dev/null || true) | head -1)" \
       "| $( (mpiexec --version 2>/dev/null || true) | head -1)"
  echo "git status (tracked files):"
  git status --short --untracked-files=no
  sha256sum amps build/pic/picGlobal.dfn build/pic/pic_mover.cpp
} > build_provenance.txt
cat build_provenance.txt
# Fingerprints of the campaign build (validation_runs/i3_fixed_campaign_20261122, job 137954): the generated
# compile-time configuration (without its two build-date comment lines) and the mover. A different amps
# SHA-256 alone is expected (paths, compiler).
want_dfn=d82c0e8a6046026bdb1e2d5bb0ad2628f8db6ceebe0a08f276fe688072bd5063
want_mover=c9e0b74548d67248f4b49bddb08ef111ec32f2dc9a5cab1050e5d0e64620a151
have_dfn=$(grep -v 'are added while reading input file\. \[' build/pic/picGlobal.dfn | sha256sum | cut -d' ' -f1)
have_mover=$(sha256sum build/pic/pic_mover.cpp | cut -d' ' -f1)
for item in "build/pic/picGlobal.dfn (without date lines):$have_dfn:$want_dfn" "build/pic/pic_mover.cpp:$have_mover:$want_mover"; do
  IFS=: read -r name have want <<< "$item"
  if [[ $have == "$want" ]]; then echo "same as the campaign build: $name"
  else echo "DIFFERENT from the campaign build: $name (expected for --test-mode or another commit)"; fi
done
echo "BUILD DONE: $DEST/amps"
