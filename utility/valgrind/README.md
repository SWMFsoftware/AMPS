# Uninitialized-memory sweep with valgrind memcheck

Tools for finding reads of uninitialized memory in AMPS runs: values that can
change results, make two runs of the same executable differ, or leak into
output files.

| file | purpose |
|------|---------|
| `run_memcheck.sh` | runs an executable under memcheck, with MPI, one log per rank, then summarizes |
| `summarize_memcheck.py` | groups the errors of memcheck logs by kind, first AMPS frame and origin |
| `memcheck_preload.c` | library preloaded into every rank: reports undefined bytes in MPI send buffers and in the data of every `fwrite()` |
| `amps-memcheck.supp` | suppressions for MPI/system noise only (no AMPS frames) |

## Quick start

```sh
# 1. build the test with debug information (gcc -g does not change the generated code); -check_macro=off keeps the
#    build/ copies of the sources line-for-line equal to src/ (Config.pl -check_macro=on afterwards restores it)
./Config.pl -check_macro=off
make test_Moon_compile FLAGC_EXTRA=-g
make test_Moon_rundir

# 2. run it in a scratch copy of the run directory, 2 ranks, at most 2 hours, on CPUs 32-35
mkdir -p /tmp/vg-Moon && cp run_test_Moon/* /tmp/vg-Moon/ && cd /tmp/vg-Moon
$AMPS/utility/valgrind/run_memcheck.sh -n 2 -t 7200 -c 32-35 -o vg ./amps

# 3. read vg/summary.txt; vg/memcheck.<rank>.log has the full stacks
$AMPS/utility/valgrind/summarize_memcheck.py --all-frames vg/memcheck.*.log
```

Use `FLAGC_EXTRA=-g` rather than `EXTRACOMPILEROPTIONS=-g`: on the command line the latter replaces the values that
`Config.pl -compiler-option=...` stores (e.g. the gtest include path of the unit tests). With CheckMacro on (the
default) the file:line of a report refers to the `build/` copy, which has the `AMPS_CHECKMACRO` blocks inserted.

Arguments after the executable are passed to it (e.g. `./amps -mode gridless -i file.in`).
`-n 0` runs without `mpirun` (a single process).  The exit status is 1 when memcheck reported an error.

## What the options do

* `--track-origins=yes`: every report names the allocation (heap or stack frame) of the uninitialized bytes.
* `--expensive-definedness-checks=yes`: gcc -O3 compiles `strcmp(s,"LITERAL")==0` into 8-byte compares that
  read past the end of `s`.  Without this option memcheck reports each such comparison
  (in the caller of `strcmp`, e.g. the input parsers), although the result does not depend on those bytes.
* `--leak-check=no`: leaks are a different subject, and the MPI libraries report many.
* Memcheck tracks definedness inside one process only; bytes sent to another rank arrive there as
  "defined".  `memcheck_preload.c` (built by the script into the log directory and preloaded with
  `mpirun -x LD_PRELOAD=...`) checks the send buffer of `MPI_Send`/`Isend`/`Sendrecv`/`Bcast`/`Reduce`/
  `Allreduce`/`Gather(v)`/`Allgather(v)`/`Scatter(v)`/`Alltoall(v)`/`Scan`; an undefined byte is reported as
  `Uninitialised byte(s) found during client check request` with the AMPS caller and the origin.
  Only contiguous data are checked (datatypes whose size equals their extent).
* Memcheck reports undefined bytes written to a file only when stdio flushes its buffer, with the stack of the
  flush, and only once per flush site.  The same library wraps `fwrite()` and checks the data at the call, so
  each `fwrite()` that writes undefined bytes is reported at its own line.  `-p` disables both checks.
* With `-c CPULIST` the run is confined to those CPUs (`taskset`) and the ranks are not bound
  (`--bind-to none`), so that a shared machine is not loaded beyond them.

## Reading the reports

| report | meaning |
|--------|---------|
| `Conditional jump or move depends on uninitialised value(s)` | a branch on garbage: can change results |
| `Use of uninitialised value of size N` | garbage used as an address/index (e.g. a table lookup) |
| `Syscall param write(buf) points to uninitialised byte(s)` | garbage written to a file; the stack is where stdio flushes, the origin says whose bytes they are |
| `Uninitialised byte(s) found during client check request` | garbage in an MPI send buffer or in `fwrite()` data (from `memcheck_preload.c`; the caller is the frame after `MPI_*`/`fwrite`) |

Classify each finding as (a) affects results, (b) affects output files only (padding or unset fields
written), or (c) harmless/noise, and fix (a) and (b) at the allocation or in the constructor.  Add a
suppression only for MPI/system noise; `-x --gen-suppressions=all` prints the suppression text.

## Making runs short

valgrind is 20-50 times slower than a native run.  Cut-cell volumes are the slowest part of the start of
many tests: run the executable once natively and copy its `amr.sig=*.CenterNodeMeasure.v*.bin` into the
memcheck run directory, which then reads the volumes instead of computing them.  Shorten runs in the scratch copy of the run directory
(fewer particles, points, energies, iterations in the input files), never in `MakefileTest/Table`, or stop
them with `-t SECONDS`: the errors found until then stay in the logs.

## In the nightly tests

`MakefileTest/Table` has two memcheck variants of ordinary tests; the ordinary tests are unchanged:

| Table name | shown as | what it does | memcheck time |
|------------|----------|--------------|---------------|
| `Mars-ions--valgrind` | `Mars-ions(valgrind)` | `-g` build, a native 2-rank run that writes the cell-measure cache, then `run_memcheck.sh -n 2` | ~10 min |
| `unit_test_field-lines--valgrind` | `unit_test_field-lines(valgrind)` | `-g` build of the GoogleTest unit tests, `run_memcheck.sh -n 1` | ~5 min |

A variant PASSes when memcheck reported no error and `amps` ended normally (its result file then is empty);
otherwise the result file holds the error summary, and the run directory `run_test_<name>--valgrind/memcheck/` the
per-rank logs. The builds use `FLAGC_EXTRA=-g` and `Config.pl -check_macro=off` (set back to on after the build) and
remove `build/` first.

Make targets cannot contain parentheses (make reads `x(y)` as an archive member), so the Table name carries the
variant as the suffix `--valgrind`; the nightly tool (`utility/nightly`, `amps_nightly/table.py: VARIANTS`) shows
it as `<test>(valgrind)` in the results, on the pages and in the log name (`logs/legacy/tests/<test>(valgrind).log`).
A variant runs only on machines whose file `utility/nightly/machines/<machine>.json` declares the capability
`valgrind` (passion does); elsewhere it is reported as SKIPPED: "<machine> does not declare capability valgrind".
Run by hand (`make test_Mars-ions--valgrind`) on a machine without valgrind, its result file says so.

## Mesh files

`cMeshAMRgeneric::saveMeshFile()` output is byte-identical between two runs of the same executable
(on the same day: the default mesh name contains the date).  To check:

```sh
cmp run1/mesh.msh run2/mesh.msh && echo identical
```
