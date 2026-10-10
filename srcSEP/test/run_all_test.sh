#!/bin/bash
# run_all_test.sh -- build srcSEP natively, run all of its test procedures, and
# report how many tests ran/passed/failed/were skipped, with the list of failed
# tests (name, runner, message, test-log location).
#
# This wrapper keeps the historical command line; the implementation is the
# shared orchestrator tools/sep_test_orchestrator.py, whose srcSEP step table
# lists every step (build, test runners, model suites) and the result reader
# that interprets each runner's report.  All options are passed through:
#   --skip-build --quick --only STEP --list --dry-run --output-dir DIR
#   --stream --log FILE --junit FILE --ranks N --ranks-large N   (see --help)
#
# Policy: every new srcSEP test, test target, or runner must be added to that
# step table, with a result reader, when it is developed (see
# srcSEP/test/README.md).
#
# exec replaces this shell, so the orchestrator's exit status (0 = all checks
# passed and no test failed/errored, 1 = failures, 2 = usage/preflight error)
# is the script's exit status.  Python reads the whole orchestrator before
# running it, so editing it during a run cannot affect that run.
AMPS_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
exec python3 "$AMPS_ROOT/tools/sep_test_orchestrator.py" --app srcsep "$@"
