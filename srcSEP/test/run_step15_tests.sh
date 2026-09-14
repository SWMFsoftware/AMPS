#!/bin/sh
set -eu

# Run all source-distribution Step 15 cases and assemble their evidence in one
# disposable directory.  Real SWMF and observational inputs are not bundled;
# the expected local outcome is therefore INCOMPLETE, not a false release pass.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step15-campaign.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
evidence_dir="$build_dir/evidence"
campaign_dir="$build_dir/campaign"
mkdir -p "$evidence_dir" "$campaign_dir"

STEP15_REPORT_DIR="$evidence_dir" "$src_root/test/run_step15_numerical_tests.sh"
STEP15_REPORT_DIR="$evidence_dir" "$src_root/test/run_step15_swcme_tests.sh"
PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/validation/run_campaign.py" \
  --numerical-json "$evidence_dir/step15-numerical-results.json" \
  --swcme-json "$evidence_dir/step15-swcme-results.json" \
  --source-root "$src_root" --output-dir "$campaign_dir"

# The release switch is a negative test here.  Without genuine SWMF and
# spacecraft manifests it must refuse the campaign even though VAL01-VAL04 pass.
if PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/validation/run_campaign.py" \
  --numerical-json "$evidence_dir/step15-numerical-results.json" \
  --swcme-json "$evidence_dir/step15-swcme-results.json" \
  --source-root "$src_root" --output-dir "$campaign_dir" --release; then
  echo "FAIL VAL-RELEASE-GATE: incomplete campaign was accepted" >&2
  exit 1
fi

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$src_root/test/step15/test_campaign_governance.py" \
  "$campaign_dir/step15-validation-campaign.json" "$src_root"
for report in "$campaign_dir/step15-validation-campaign.json" \
              "$campaign_dir/step15-validation-campaign.md"; do
  test -s "$report" || { echo "FAIL VAL-CAMPAIGN: missing $report" >&2; exit 1; }
done

# Operators may request copies in an explicit location.  The default continues
# to remove all generated campaign artifacts when the focused test exits.
if test -n "${STEP15_REPORT_DIR:-}"; then
  mkdir -p "$STEP15_REPORT_DIR"
  cp "$evidence_dir"/* "$campaign_dir"/* "$STEP15_REPORT_DIR/"
fi
echo "PASS STEP15: source-only evidence passes; external release gates remain explicit"
