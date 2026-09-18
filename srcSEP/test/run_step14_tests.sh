#!/bin/sh
set -eu

# Step 14 is primarily a source and documentation boundary.  These checks use
# exact file lists and stable symbol patterns so historical names remain legal
# in migration documents while being impossible to compile or dispatch.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
sep_common_dir=${SEP_COMMON_DIR:-"$src_root/../src/models/sep_common"}
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step14.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# DOC01: every documented dependency-light command must exist as a make target,
# and structured-report examples must use options implemented by the parser.
for target in test-cli-unit test-state-unit test-geometry-source-unit \
  test-mover-api-unit test-field-line-scope-unit test-transport-common-unit \
  test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit \
  test-coefficients-unit test-turbulence-core-unit \
  test-reproducibility-unit test-wp01-wp10-unit \
  test-wp11-wp20-unit test-wp21-wp30-unit test-wp31-wp41-unit \
  test-wp42-wp64-experimental test-wp59-wp64-native \
  test-package-hygiene-unit test-sep-common-ownership-unit \
  test-swcme-relocation-unit test-python-runner-unit \
  test-acceptance-unit test-documentation-unit \
  test-controlled-analytical test-scientific-validation \
  test-native-amps-validation test-swmf-validation \
  test-observational-validation test-stochastic-repeat; do
  grep -q "^$target:" "$src_root/makefile" || {
    echo "FAIL DOC01: README target '$target' is missing" >&2
    exit 1
  }
done
grep -q -- '--test-json' "$src_root/util/sep_cli.cpp"
grep -q -- '--test-junit' "$src_root/util/sep_cli.cpp"
echo "PASS DOC01: documented commands map to implemented targets and options"

# DOC02: the public symbol inventory and runtime registry expose exactly the
# three canonical movers.  Retired names may appear only in Markdown migration
# history and in this negative regression pattern.
for file in "$src_root/sep.h" "$src_root/production_mover_runtime.cpp" \
            "$src_root/makefile"; do
  if grep -Eq 'ParticleMoverPtr|ParticleMover_Droge_2009_AJ|ParticleMover_Tenishev_2005_FL|ParticleMover_He_2011_AJ|ParticleMover_MeanFreePathScattering|ParticleMover_Parker_MeanFreePath|ParticleMover_Parker_Dxx|ParticleMover_FTE\>|ParticleMover_ParkerEquation' "$file"; then
    echo "FAIL DOC02: retired public mover surface remains in $file" >&2
    exit 1
  fi
done
test ! -e "$src_root/mover.cpp"
test ! -e "$src_root/fte_mover.cpp"
test ! -e "$src_root/fte_mover_dmumu.cpp"
# The compiled Step 4 test asserts Registry().size()==3 and checks discovery,
# parsing, capability metadata, and rejection of every retired alias.
"$src_root/test/run_step4_tests.sh" >/dev/null
echo "PASS DOC02: public source contains exactly the canonical mover surface"

# DOC03: delivery archives must contain sources, not products from an earlier
# build or another delivery nested inside this one.
if find "$src_root" -type f \( -name '*.o' -o -name '*.a' -o \
     -name '*.gcda' -o -name '*.gcno' -o -name '*.profraw' -o \
     -name '*.pyc' -o \
     -name 'controlled-*-results.json' -o \
     -name 'controlled-*-results.xml' -o \
     -name 'step15-*-results.json' -o -name 'step15-*-results.xml' -o \
     -name '*.tar' -o -name '*.tar.gz' -o -name '*.tgz' -o \
     -name '*.zip' \) -print -quit | grep -q .; then
  echo "FAIL DOC03: generated binary, coverage, or archive product is present" >&2
  exit 1
fi
if find "$src_root" -type d -name '__pycache__' -print -quit | grep -q .; then
  echo "FAIL DOC03: Python bytecode cache is present" >&2
  exit 1
fi
grep -q '^\*\.tar\.gz$' "$src_root/.gitignore"
echo "PASS DOC03: source tree is free of generated and nested delivery products"

# WARN-FULL-SOURCE: these two host-dependent translation units cannot be
# compiled by the source-only archive because their PIC/field-line generated
# headers come from the enclosing AMPS build.  Guard the exact regressions that
# previously escaped dependency-light compilation: a jump across initialized
# C++ objects and an int/size_type variadic format mismatch.  The authoritative
# completion gate remains the native AMPS warning-clean build.
if grep -Eq 'goto[[:space:]]+end[[:space:]]*;' "$src_root/sampling.cpp"; then
  echo "FAIL WARN-FULL-SOURCE: sampling output jumps across initialized objects" >&2
  exit 1
fi
if grep -Eq 'I=%i.*field_line->size\(\)' "$src_root/mesh.cpp"; then
  echo "FAIL WARN-FULL-SOURCE: field-line size is printed with an int format" >&2
  exit 1
fi
grep -Eq 'I=%zu.*field_line->size\(\)' "$src_root/mesh.cpp" || {
  echo "FAIL WARN-FULL-SOURCE: field-line size has no size_t-safe output format" >&2
  exit 1
}
echo "PASS WARN-FULL-SOURCE: sampling control flow and mesh size format are safe"

# WARN01: compile every dependency-light production translation unit touched by
# Step 13/14 under the project's strictest portable C++11 warning policy.
for source in sep_test_registry sep_transport_common sep_coefficient_physics \
  sep_species_source sep_coefficient_registry sep_background_snapshot; do
  ${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
    -I"$sep_common_dir" -fsyntax-only "$sep_common_dir/$source.cpp"
done
# Application-owned units are checked separately.  This split is deliberate:
# a Step 14 run must compile the canonical shared source, never resurrect a
# private srcSEP/util implementation merely to satisfy a documentation gate.
for source in sep_focused_transport_core sep_focused_transport_mfp_core \
  sep_acceptance_cases sep_production_mover sep_cli; do
  definitions=
  case "$source" in
    sep_cli) definitions=-DSEP_CLI_PARSE_ONLY ;;
  esac
  ${CXX:-c++} -std=c++11 -Wall -Wextra -Wpedantic -Werror $definitions \
    -I"$sep_common_dir" -I"$src_root/util" \
    -fsyntax-only "$src_root/util/$source.cpp"
done
# The controlled component catalogs are now linked into the production test
# registry. Compile both under the same strict warning policy; fsyntax-only is
# sufficient here because their focused sanitizer targets perform full links.
${CXX:-c++} -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -I"$sep_common_dir" -I"$src_root" -I"$src_root/util" -fsyntax-only \
  "$src_root/util/sep_mover_validation.cpp" \
  "$src_root/util/sep_turbulence_validation.cpp"
echo "PASS WARN01: acceptance and registered validation sources compile warning-free"

# A compiler analyzer is used when available; absence is reported explicitly
# without pretending that the optional external cppcheck tool ran.
if ${CXX:-c++} --help=common 2>/dev/null | grep -q -- '-fanalyzer'; then
  ${CXX:-c++} -std=c++11 -Wall -Wextra -fanalyzer -DSEP_CLI_PARSE_ONLY \
    -I"$sep_common_dir" -I"$src_root/util" -fsyntax-only \
    "$sep_common_dir/sep_test_registry.cpp" \
    "$src_root/util/sep_acceptance_cases.cpp" \
    "$src_root/util/sep_cli.cpp"
  echo "PASS STATIC01: compiler static analyzer accepted new infrastructure"
else
  echo "SKIP STATIC01: compiler does not advertise -fanalyzer"
fi

# SAN01 source-only subset: Step 13 compiles the exact new callbacks with ASan
# and UBSan. Earlier sanitizer targets are composed by `make test-sanitizer` so
# this script stays non-recursive when invoked from that aggregate target.
"$src_root/test/run_step13_tests.sh"
echo "PASS SAN01: Step 13/14 acceptance infrastructure passed ASan and UBSan"
