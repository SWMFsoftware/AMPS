#!/bin/sh
# AMPS nightly testing, script 1 of 3: run the test procedure on this machine.
# (Replaces ~/bin/run_test_amps.sh of the old routine.)
#
# Clones/updates AMPS master (and share, util, BATL) from GitHub into local mirrors, takes one commit for
# the whole night, installs it in a fresh workspace for every runner and every configuration of this
# machine's file (machines/<machine>.json), and runs the runners of runners.json one after another:
# src/test core suites, srcEnceladus, srcSEP, srcSEP3D, srcMoon, srcEarth, Makefile.test, ...
# Results and logs stay on this machine in <workdir>/nights/<YYYYMMDD-HHMM>/runs/<column>/ until
# nightly_pack.sh sends them to the server.
#
#   nightly_run.sh                         the whole procedure (about 7 h on passion)
#   nightly_run.sh --only moon             one runner (see runners.json for the ids)
#   nightly_run.sh --dry-run               print the commands without running anything
#   nightly_run.sh --help                  every option
#
# crontab of a test machine (passion is also the server, so it adds script 3):
#   00 20 * * *  D=<AMPS>/utility/nightly; $D/nightly_run.sh; $D/nightly_pack.sh; $D/nightly_update.sh
here=$(cd "$(dirname "$0")" && pwd)
exec "${AMPS_NIGHTLY_PYTHON:-python3}" "$here/amps-nightly" run "$@"
