#!/bin/sh
# AMPS nightly testing: the script cron starts on a test machine.
# Its installed copy lives in ~/bin; this file, utility/nightly/bin/amps_nightly_test.sh in AMPS, is the
# version under version control. Install it, or update the copy in ~/bin, with
#     <AMPS>/utility/nightly/amps-nightly install
# which also writes ~/bin/amps_nightly_test.conf and the crontab line (it shows them and asks first).
#
# Each night it
#   1. takes a lock, so that only one night runs at a time on this machine;
#   2. downloads AMPS: git clone --depth 1 of $BRANCH into $TESTDIR/AMPS (the previous night's copy
#      is kept as $TESTDIR/AMPS.prev, and used if the download fails);
#   3. pulls SWMF_data: the nightly's own copy, at the "swmf_data" var of this machine's file
#      (cloned the first time, then updated to the newest commit of its branch);
#   4. runs the orchestrator of the downloaded copy on the downloaded commit:
#      $TESTDIR/AMPS/utility/nightly/nightly_run.sh --commit <sha>
#   5. packs the results and sends them to the server: nightly_pack.sh
#   6. on the results server only (UPDATE_SITE=yes): nightly_update.sh (results.db, archive, web pages).
#
# Where AMPS is installed: the download in $TESTDIR/AMPS only provides the orchestrator and its files;
# Config.pl is never run in it. The orchestrator installs AMPS once for every runner (runners.json) and
# every compiler of the machine file, as the old routine's $TESTDIR/GNU/AMPS:
#     $TESTDIR/<COMPILER>/<runner>/AMPS     e.g. ~/Tmp_AMPS_test/GNU/legacy/AMPS, ~/Tmp_AMPS_test/GNU/core/AMPS
# a copy of the night's commit with share, util and BATL and a link to SWMF_data, in which it runs the
# "config_pl" calls (./Config.pl -install -compiler=..., ...) and "after_install" lines (make share_lib) of
# machines/<machine>.json before the runner starts. Every night's log lists these commands (step 1), and
# each runner's log (logs/<runner>/runner.log, "### configure") shows them running.
# The orchestrator, runners.json and machines/ come from the download, so a change to them pushed to $BRANCH
# takes effect the next night. This script does not replace itself: when the downloaded copy has a
# different version of it, the night's log and the night's summary page say so, and
# 'amps-nightly install' updates ~/bin.
#
# Settings: amps_nightly_test.conf next to this script (written by 'amps-nightly install'):
#   AMPSTEST=/data/vtenishe/AMPSTEST   the nightly's logs/ and SWMF_data/ (on the server also its results)
#   TESTDIR=~/Tmp_AMPS_test            where the tests run: the "workdir" of the machine file
#   BRANCH=master                      the AMPS branch that is downloaded and tested
#   UPDATE_SITE=no                     yes on the results server (passion)
#   PYTHON=python3                     the Python (3.8 or newer) that runs amps-nightly
#   rarely needed: REPOSITORY, SWMF_DATA_REPOSITORY, SWMF_DATA_BRANCH, KEEP_LOG_DAYS (defaults below)
# Arguments are passed on to nightly_run.sh, e.g. 'amps_nightly_test.sh --only moon --no-upload' for a trial.
# Log of each night: $AMPSTEST/logs/nightly-YYYYMMDD-HHMM.log, kept KEEP_LOG_DAYS days.

here=$(cd "$(dirname "$0")" && pwd)
self="$here/$(basename "$0")"
conf=${AMPS_NIGHTLY_CONF:-$here/amps_nightly_test.conf}
if [ ! -r "$conf" ]; then
    echo "amps_nightly_test.sh: settings file $conf not found; run 'amps-nightly install'" >&2
    exit 2
fi
. "$conf"
if [ -z "${AMPSTEST:-}" ]; then
    echo "amps_nightly_test.sh: AMPSTEST is not set in $conf" >&2
    exit 2
fi
TESTDIR=${TESTDIR:-$AMPSTEST/tool}
BRANCH=${BRANCH:-master}
UPDATE_SITE=${UPDATE_SITE:-no}
PYTHON=${PYTHON:-python3}
REPOSITORY=${REPOSITORY:-git@github.com:SWMFsoftware/AMPS.git}
SWMF_DATA_REPOSITORY=${SWMF_DATA_REPOSITORY:-git@github.com:SWMFsoftware/SWMF_data.git}
SWMF_DATA_BRANCH=${SWMF_DATA_BRANCH:-master}
KEEP_LOG_DAYS=${KEEP_LOG_DAYS:-60}
PATH=/usr/local/bin:/usr/bin:/bin:$PATH          # cron's PATH is only /usr/bin:/bin
AMPS_NIGHTLY_PYTHON=$PYTHON                     # used by nightly_run.sh, nightly_pack.sh, nightly_update.sh
export PATH AMPS_NIGHTLY_PYTHON

mkdir -p "$AMPSTEST/logs" "$TESTDIR" || exit 2
log="$AMPSTEST/logs/nightly-$(date +%Y%m%d-%H%M).log"
exec >>"$log" 2>&1

say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }
notes=""
warn() { say "WARNING: $*"; notes="${notes:+$notes; }WARNING: $*"; }

# 1. the lock: flock on descriptor 9, held until this script and the orchestrator have finished
exec 9>>"$AMPSTEST/nightly.lock"
"$PYTHON" -c 'import fcntl, sys
try:
    fcntl.flock(9, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    sys.exit(75)'
case $? in
    0) ;;
    75) say "another night is still running on this machine (lock $AMPSTEST/nightly.lock); nothing done"
        exit 0 ;;
    *)  say "ERROR: cannot take the lock $AMPSTEST/nightly.lock with $PYTHON"
        exit 2 ;;
esac
say "AMPS nightly test on $(hostname); settings $conf; script $self"

# 2. download AMPS
tool=$TESTDIR
rm -rf "$tool/AMPS.new"
say "downloading AMPS ($BRANCH) from $REPOSITORY"
if git clone --quiet --depth 1 --branch "$BRANCH" "$REPOSITORY" "$tool/AMPS.new"; then
    rm -rf "$tool/AMPS.prev"
    if [ -d "$tool/AMPS" ]; then
        mv "$tool/AMPS" "$tool/AMPS.prev"
    fi
    mv "$tool/AMPS.new" "$tool/AMPS"
elif [ -d "$tool/AMPS" ]; then
    rm -rf "$tool/AMPS.new"
    warn "downloading AMPS failed; this night uses the copy downloaded the night before"
else
    say "ERROR: downloading AMPS failed and there is no earlier copy; nothing can run"
    exit 1
fi
nightly=$tool/AMPS/utility/nightly
if [ ! -f "$nightly/amps-nightly" ]; then
    say "ERROR: $nightly/amps-nightly is missing: utility/nightly is not on branch $BRANCH"
    exit 1
fi
commit=$(git -C "$tool/AMPS" rev-parse HEAD)
say "AMPS $BRANCH is at $commit"
if ! cmp -s "$self" "$nightly/bin/amps_nightly_test.sh"; then
    warn "$self differs from utility/nightly/bin/amps_nightly_test.sh on $BRANCH; update it with: $nightly/amps-nightly install"
fi

# 3. pull SWMF_data
data=$("$PYTHON" "$nightly/amps-nightly" machines --var swmf_data 2>/dev/null) || data=""
if [ -z "$data" ]; then
    say "this machine's file sets no swmf_data: SWMF_data is not updated"
elif [ -d "$data/.git" ]; then
    before=$(git -C "$data" rev-parse HEAD 2>/dev/null)
    if git -C "$data" fetch --quiet --depth 1 origin "$SWMF_DATA_BRANCH" &&
       git -C "$data" reset --quiet --hard FETCH_HEAD; then
        git -C "$data" clean -q -f -d
        after=$(git -C "$data" rev-parse HEAD)
        if [ "$before" != "$after" ]; then
            git -C "$data" reflog expire --expire=now --all
            git -C "$data" gc --quiet --prune=now
            say "SWMF_data ($data) updated from $before to $after"
        else
            say "SWMF_data ($data) is up to date at $after"
        fi
    else
        warn "updating SWMF_data in $data failed; the tests use it as it was"
    fi
else
    say "cloning SWMF_data into $data"
    mkdir -p "$(dirname "$data")"
    if ! git clone --quiet --depth 1 --branch "$SWMF_DATA_BRANCH" "$SWMF_DATA_REPOSITORY" "$data"; then
        warn "cloning SWMF_data into $data failed"
    fi
fi

# 4.-6. run, pack, update
say "the night's plan: where each runner is installed, the install commands and each runner's command"
"$PYTHON" "$nightly/amps-nightly" run --dry-run "$@" | sed 's/^/    /'
say "step 1/3: nightly_run.sh --commit $commit $*"
if [ -n "$notes" ]; then
    "$nightly/nightly_run.sh" --commit "$commit" --note "$notes" "$@"
else
    "$nightly/nightly_run.sh" --commit "$commit" "$@"
fi
run_status=$?
say "step 2/3: nightly_pack.sh"
"$nightly/nightly_pack.sh"
pack_status=$?
update_status=0
if [ "$UPDATE_SITE" = yes ]; then
    say "step 3/3: nightly_update.sh"
    "$nightly/nightly_update.sh"
    update_status=$?
else
    say "step 3/3 (nightly_update.sh) runs on the results server, not here (UPDATE_SITE=$UPDATE_SITE)"
fi
find "$AMPSTEST/logs" -name 'nightly-*.log' -mtime +"$KEEP_LOG_DAYS" -exec rm -f {} + 2>/dev/null
say "finished: nightly_run.sh $run_status, nightly_pack.sh $pack_status, nightly_update.sh $update_status"
[ "$run_status" -eq 0 ] && [ "$pack_status" -eq 0 ] && [ "$update_status" -eq 0 ]
