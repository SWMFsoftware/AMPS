#!/bin/sh
# Run an AMPS executable under valgrind memcheck with MPI and summarize the errors.
#
# Usage: run_memcheck.sh [-n NP] [-o LOGDIR] [-t SECONDS] [-c CPULIST] [-x "VALGRIND OPTIONS"] [--] ./amps [ARGS]
#
#   -n NP        number of MPI ranks (default 2); NP=0 runs the executable without mpirun
#   -o LOGDIR    directory for the per-rank logs memcheck.<rank>.log and summary.txt
#                (default ./memcheck; it is emptied first)
#   -t SECONDS   stop the run after SECONDS (SIGTERM); the errors found until then stay in the logs
#   -c CPULIST   run on these CPUs only (taskset -c CPULIST); the ranks are never bound (--bind-to none)
#   -x OPTIONS   extra valgrind options, e.g. "--gen-suppressions=all"
#   -p           do not preload the checks of MPI send buffers and fwrite() data (memcheck_preload.c)
#
# The run is made in the current directory (an AMPS run directory with its input files).
# Build the executable with debug information for file:line in the reports; -g does not change
# the code that gcc generates, e.g.
#   make test_Moon_compile EXTRACOMPILEROPTIONS=-g
# Memcheck does not follow data from one rank to another, and it reports undefined file data only where
# stdio flushes; memcheck_preload.c, which the script builds into LOGDIR and preloads, checks the send
# buffers of the MPI calls and the data of every fwrite() at the call.
# Suppressions: amps-memcheck.supp next to this script, the Open MPI file when it is installed,
# and every file in $MEMCHECK_SUPP (space-separated).
#
# Exit status: 0 no memcheck error, 1 errors were found (see LOGDIR/summary.txt), 2 usage error.
# The exit status of the program itself is in LOGDIR/exit_status.

here=$(cd "$(dirname "$0")" && pwd)
np=2
logdir=./memcheck
limit=
cpus=
extra=
preload=yes

while [ $# -gt 0 ]; do
  case "$1" in
    -n) np=$2; shift 2 ;;
    -o) logdir=$2; shift 2 ;;
    -t) limit=$2; shift 2 ;;
    -c) cpus=$2; shift 2 ;;
    -x) extra=$2; shift 2 ;;
    -p) preload=no; shift ;;
    --) shift; break ;;
    -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
    -*) echo "run_memcheck.sh: unknown option $1" >&2; exit 2 ;;
    *) break ;;
  esac
done

[ $# -ge 1 ] || { echo "run_memcheck.sh: no executable given (try -h)" >&2; exit 2; }
command -v valgrind >/dev/null || { echo "run_memcheck.sh: valgrind is not installed" >&2; exit 2; }

mkdir -p "$logdir" || exit 2
rm -f "$logdir"/memcheck.*.log "$logdir"/summary.txt "$logdir"/exit_status
logdir=$(cd "$logdir" && pwd)

supp="--suppressions=$here/amps-memcheck.supp"
for f in /usr/share/openmpi/openmpi-valgrind.supp /usr/share/openmpi/openmpi-valgrind.supp.txt $MEMCHECK_SUPP; do
  [ -f "$f" ] && supp="$supp --suppressions=$f"
done

# --leak-check=no: leaks are not the subject here and the MPI libraries report many of them.
# --expensive-definedness-checks=yes: gcc -O3 turns strcmp(s,"LITERAL")==0 into 8-byte compares that read past
#   the end of s; without the option every such comparison is reported ("Conditional jump ... depends on
#   uninitialised value" in the function that calls strcmp), although its result does not depend on those bytes.
# %q{OMPI_COMM_WORLD_RANK} gives one log per rank (Open MPI); with -n 0 the log is memcheck.0.log.
vgopts="--tool=memcheck --track-origins=yes --expensive-definedness-checks=yes --error-limit=no --num-callers=40 --leak-check=no \
 --read-inline-info=yes $supp $extra"

lib=
if [ $preload = yes ]; then
  if mpicc -shared -fPIC -O1 -g "$here/memcheck_preload.c" -o "$logdir/libmemcheck_preload.so" -ldl; then
    lib=$logdir/libmemcheck_preload.so
  else
    echo "run_memcheck.sh: cannot build memcheck_preload.c, MPI send buffers and fwrite() data are not checked" >&2
  fi
fi

pre=
[ -n "$cpus" ] && pre="taskset -c $cpus"
to=
[ -n "$limit" ] && to="timeout -s TERM $limit"

if [ "$np" = 0 ]; then
  OMPI_COMM_WORLD_RANK=0 LD_PRELOAD=$lib $to $pre valgrind $vgopts --log-file="$logdir/memcheck.%q{OMPI_COMM_WORLD_RANK}.log" "$@"
else
  # ranks are not bound: valgrind runs a rank's threads one at a time anyway, and bound ranks of two jobs would
  # share cores (Open MPI binds every job from core 0)
  bind="--bind-to none"
  xpre=
  [ -n "$lib" ] && xpre="-x LD_PRELOAD=$lib"
  $to $pre mpirun -np "$np" $bind $xpre valgrind $vgopts --log-file="$logdir/memcheck.%q{OMPI_COMM_WORLD_RANK}.log" "$@"
fi
echo $? > "$logdir/exit_status"

/usr/bin/python3 "$here/summarize_memcheck.py" "$logdir"/memcheck.*.log > "$logdir/summary.txt"
status=$?
cat "$logdir/summary.txt"
echo "logs: $logdir (program exit status $(cat "$logdir/exit_status"))"
exit $status
