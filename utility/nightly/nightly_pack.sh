#!/bin/sh
# AMPS nightly testing, script 2 of 3: collect the results of the night, tar them and send them to the server.
# (Replaces ~/bin/scp_test_amps.sh of the old routine.)
#
# For every run of nightly_run.sh not sent yet (one per configuration = one column of the report), writes
# <workdir>/outbox/<run_id>.tar.gz with everything the web page shows: run.json (machine, compiler, commit,
# times), results.jsonl (every test: status, stage, reason), runners.json (what each runner did) and logs/
# (each runner's console log, one log per test, the files the runners write, the orchestrator log).
# Then copies the tarball to the "upload" destination of the machine file: a local directory on the server
# (passion) or [user@]host:directory through rsync. A tarball that could not be sent stays in the outbox
# and is sent the next time. A run whose nightly_run.sh died (reboot, killed) is completed from the
# progress it saved after each runner and sent marked "incomplete".
#
#   nightly_pack.sh                        pack and send
#   nightly_pack.sh --no-send              only write the tarballs (in <workdir>/outbox/)
here=$(cd "$(dirname "$0")" && pwd)
exec "${AMPS_NIGHTLY_PYTHON:-python3}" "$here/amps-nightly" pack "$@"
