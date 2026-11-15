#!/bin/sh
# AMPS nightly testing, script 3 of 3 (on the server, passion): take the tarballs and update the results
# database, the log archive and the web pages.
# (Replaces ~/bin/process_amps_tests.sh of the old routine.)
#
# Unpacks every <run_id>.tar.gz in ~/AMPSTEST/nightly/incoming/ (or the tarballs given as arguments),
# checks it, adds its results to ~/AMPSTEST/nightly/results.db, moves its logs to
# www/archive/YYYY/MM/DD/<machine>-<config>/ (with a night summary index.html), rebuilds www/index.html and
# the family pages, and gzips archived logs after 14 nights. A tarball that fails the checks goes to
# rejected/ with REJECTED.txt.
#
#   nightly_update.sh                                   everything that arrived in incoming/
#   nightly_update.sh /path/to/<run_id>.tar.gz          a tarball brought by hand
#
# crontab of the server, for test machines other than passion (passion runs it right after its own pack):
#   30 10 * * *  <AMPS>/utility/nightly/nightly_update.sh
here=$(cd "$(dirname "$0")" && pwd)
exec "${AMPS_NIGHTLY_PYTHON:-python3}" "$here/amps-nightly" update "$@"
