#!/usr/bin/env bash
# Exercise Stage-3 input records in an isolated generated-source tree.  This
# runner intentionally reuses the current clean build as a source template;
# it never edits build/ and performs no compilation.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
template="$source_root/build"
assembled="$source_root/halley1d.input.Assembled.Block"
species="$source_root/halley1d.input.Assembled.Species"
if [[ ! -d "$template/main" || ! -f "$assembled" || ! -f "$species" ]]; then
  echo "FAIL Stage-3 parser tests require a configured clean build" >&2
  exit 1
fi

test_root=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage3-parser.XXXXXX")
trap 'rm -rf -- "$test_root"' EXIT

prepare_case() {
  local name=$1
  local case_root="$test_root/$name"
  mkdir -p "$case_root/work/main" "$case_root/work/pic" \
    "$case_root/work/meshAMR" "$case_root/work/models/dust" \
    "$case_root/work/models/exosphere"
  cp "$template/main"/*.cpp "$template/main"/*.h \
    "$template/main"/*.dfn "$template/main/makefile" "$case_root/work/main/"
  cp "$template/pic/pic.h" "$template/pic/picGlobal.dfn" \
    "$case_root/work/pic/"
  cp "$template/meshAMR/meshAMR_UserDefinitions.h" \
    "$template/meshAMR/meshAMRdef.h" "$case_root/work/meshAMR/"
  cp "$template/models/dust/Dust.cpp" "$template/models/dust/Dust.dfn" \
    "$case_root/work/models/dust/"
  cp "$template/models/exosphere/Exosphere_Chemistry.cpp" \
    "$case_root/work/models/exosphere/"

  # The configured build is useful as a complete parser fixture, but its
  # application definitions already contain the values selected by the real
  # Halley input.  Restore the two Stage-3 definition files from the source
  # defaults before every case.  This makes the omitted-record test meaningful
  # and prevents one configured build from silently supplying a seed or a
  # sampling mode that the test case did not request.
  cp "$source_root/srcHartley2/CometModel.cpp" \
    "$source_root/srcHartley2/CometChemistry.cpp" \
    "$case_root/work/main/"
  cp "$assembled" "$case_root/case.Block"
  cp "$species" "$case_root/case.Species"
  printf '%s\n' "$case_root"
}

run_success() {
  local case_root=$1
  perl "$source_root/input/hartley2.pl" "$case_root/case.Block" \
    "$case_root/work" "$source_root" >"$case_root/parser.log" 2>&1
}

run_failure() {
  local case_root=$1
  local expected=$2
  if run_success "$case_root"; then
    echo "FAIL parser accepted invalid case $(basename "$case_root")" >&2
    exit 1
  fi
  if ! rg -q "$expected" "$case_root/parser.log"; then
    echo "FAIL parser emitted the wrong error for $(basename "$case_root")" >&2
    cat "$case_root/parser.log" >&2
    exit 1
  fi
}

case_root=$(prepare_case valid)
run_success "$case_root"
rg -q 'std::uint64_t Comet::Model::RandomSeed=104729;' \
  "$case_root/work/main/CometModel.cpp"
rg -q 'bool Comet::Model::RandomSeedConfigured=true;' \
  "$case_root/work/main/CometModel.cpp"
rg -q 'ResidualSampling::RandomPhaseSystematic' \
  "$case_root/work/main/CometChemistry.cpp"

case_root=$(prepare_case omitted_defaults)
sed -i '/^167 halley1d.input$/,+1d;/^172 halley1d.input$/,+1d' \
  "$case_root/case.Block"
run_success "$case_root"
rg -q 'bool Comet::Model::RandomSeedConfigured=false;' \
  "$case_root/work/main/CometModel.cpp"
rg -q 'ResidualSampling::IndependentBernoulli;' \
  "$case_root/work/main/CometChemistry.cpp"

case_root=$(prepare_case invalid_seed)
sed -i 's/CometRandomSeed=104729/CometRandomSeed=0/' "$case_root/case.Block"
run_failure "$case_root" 'CometRandomSeed must be an integer'

case_root=$(prepare_case duplicate_seed)
sed -i '/CometRandomSeed=104729/a 999 stage3-parser-test\nCometRandomSeed=5' \
  "$case_root/case.Block"
run_failure "$case_root" 'duplicate CometRandomSeed'

case_root=$(prepare_case invalid_mode)
sed -i 's/ChemistryResidualSampling=systematic/ChemistryResidualSampling=invalid/' \
  "$case_root/case.Block"
run_failure "$case_root" 'must be BERNOULLI or SYSTEMATIC'

case_root=$(prepare_case duplicate_mode)
sed -i '/ChemistryResidualSampling=systematic/a 999 stage3-parser-test\nChemistryResidualSampling=bernoulli' \
  "$case_root/case.Block"
run_failure "$case_root" 'duplicate ChemistryResidualSampling'

case_root=$(prepare_case missing_seed)
sed -i '/^167 halley1d.input$/,+1d' "$case_root/case.Block"
run_failure "$case_root" 'SYSTEMATIC requires an explicit CometRandomSeed'

case_root=$(prepare_case chemistry_off)
sed -i 's/ChemistryMode=on/ChemistryMode=off/' "$case_root/case.Block"
run_failure "$case_root" 'ChemistryResidualSampling requires ChemistryMode=ON'

echo "PASS Stage-3 parser valid=1 omitted-defaults=1 invalid=6"
