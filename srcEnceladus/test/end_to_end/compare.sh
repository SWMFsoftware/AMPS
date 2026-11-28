#!/bin/sh
# Model-observation comparison of the Enceladus end-to-end Table tests
# (srcEnceladus/test/end_to_end/README.md). Never gating: the Table runs it
# after check.sh with its output in the test log.
#
# usage: sh srcEnceladus/test/end_to_end/compare.sh RUN_DIR REQUEST ZONES
# Run from the AMPS root. Writes RUN_DIR/comparison/:
#   observation.REQUEST.comparison.png   last-window model vs observation
#                                        (srcEnceladus/test/plot_observation_comparisons.py)
#   observation_comparison_summary.json  its metrics
#   flyby_density.E3.png                 E3 only: model density along the dense E3 track
# and prints the metrics: INMS chi-square and model/observation ratios; UVIS
# column-to-optical-depth shape fit and model peak column / published peak.
# These are diagnostics of the fitted source transferred into AMPS (gravity,
# photolysis, finite AMR box), not a refit, convergence or validation result.
#
# Python: the first of $ENCELADUS_PYTHON, python3, /usr/bin/python3 that has
# matplotlib; without one the plots are skipped and the reason is printed.

D=$1; REQUEST=$2; ZONES=$3
if [ -z "$D" ] || [ -z "$REQUEST" ] || [ -z "$ZONES" ]; then
    echo "usage: compare.sh RUN_DIR REQUEST ZONES"; exit 0
fi
PY=""
for P in "$ENCELADUS_PYTHON" python3 /usr/bin/python3; do
    [ -n "$P" ] || continue
    if "$P" -c 'import matplotlib' >/dev/null 2>&1; then PY=$P; break; fi
done
if [ -z "$PY" ]; then
    echo "comparison skipped: no python3 with matplotlib (set ENCELADUS_PYTHON)"; exit 0
fi
OBS=$D/PT/plots/observation.$REQUEST.dat
if [ ! -f "$OBS" ]; then echo "comparison skipped: $OBS is missing"; exit 0; fi
OUT=$D/comparison
mkdir -p "$OUT"
# Name the run's source distribution in the plot titles: the first .in asset of
# run_manifest.json, e.g. tenishev2014_hybrid_8__independent_each_observation__E3__H2O__absolute.in
# -> "tenishev2014_hybrid_8 fit to E3 alone" (another strategy: "<basis> <strategy> fit for <target>").
SOURCE=$(grep -o '"path": *"[^"]*\.in"' "$D/run_manifest.json" 2>/dev/null | head -1 | sed 's/.*\///; s/\.in"$//')
LABEL=$(echo "$SOURCE" | awk -F'__' 'NF == 5 && $2 == "independent_each_observation" { print $1 " fit to " $3 " alone"; next }
                                      NF == 5 { print $1 " " $2 " fit for " $3; next } { print }')
[ -n "$LABEL" ] || LABEL="source distribution of this run"
"$PY" -B srcEnceladus/test/plot_observation_comparisons.py "$OBS" \
    --output-dir "$OUT" --expected-zones "$ZONES" --source-label "$LABEL" || { echo "comparison plot failed"; exit 0; }
echo "COMPARISON source $LABEL"
FLYBY=$D/PT/plots/flyby_density.E3.dat
if [ -f "$FLYBY" ]; then
    "$PY" -B srcEnceladus/test/plot_flyby_density.py "$FLYBY" --output "$OUT/flyby_density.E3.png" \
        || echo "flyby density plot failed"
fi
"$PY" -B - "$OUT/observation_comparison_summary.json" <<'EOF'
import json, sys
summary = json.load(open(sys.argv[1]))
for item in summary["comparisons"]:
    head = "COMPARISON {} final window {} of {}".format(item["instrument"], item["final_output"] + 1, item["zones"])
    if item["instrument"] == "UVIS":
        fit, scale = item["shape_fit"], item["absolute_column"]
        print("{}: shape fit chi2={:.4g} over {} samples, column-to-optical-depth={:.4g} m^2, baseline={:.4g}; "
              "model peak column {:.4g} m^-2 = {:.3g} x published".format(
                  head, fit["chi2"], item["samples"], fit["column_to_optical_depth"], fit["signed_baseline"],
                  scale["model_peak_m-2"], scale["model_peak_over_published"]))
    else:
        print("{}: {}".format(head, json.dumps(item["metrics"], sort_keys=True)))
    print("COMPARISON plot {}".format(item["png"]))
EOF
exit 0
