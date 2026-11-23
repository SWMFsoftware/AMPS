# UTestRunner: test-runner regression suite

Run from the AMPS repository root:

```bash
./srcEarth/test/UTestRunner/run_test.sh
```

The suite uses only the Python standard library. It imports the production
`srcEarth/test/test_runner.py` directly and does not need AMPS, MPI, SPICE,
Geopack, a compiler, or observational data.

## Permanent gates

### UTR-F01 — exclusive-directive parser and provenance

This test verifies that `! runner: exclusive`:

- applies to a scalar entry without changing its command, expected P/F state,
  or `last pass:` value;
- applies to every expansion when placed before a `for` declaration;
- does not change ordinary comments or ordinary entries;
- rejects misspelled, detached, and dangling directives instead of silently
  scheduling a different test; and
- survives a `last pass:` update byte-for-byte except for the intended commit
  value.

### UTR-F02 — exclusive pool ownership

This test runs the real asynchronous dispatcher with timed in-memory jobs. Two
ordinary jobs must overlap, proving that the test is not merely exercising a
serial runner. The following exclusive job must wait for both to finish, and
the next ordinary job must wait for the exclusive job. Consequently, a failure
detects either side of the required invariant:

1. an exclusive test starts only when the runner-owned pool is empty; and
2. no other test launches while an exclusive test is running.

### UTR-F03 — child-environment exclusivity handshake

The production dispatcher is exercised through its real subprocess path. An
exclusive child must receive exactly
`AMPS_TEST_RUNNER_EXCLUSIVE=1`; an ordinary child must see the variable unset.
The test deliberately seeds a bogus value in the parent environment and also
verifies that the parent mapping remains unchanged. This prevents a stale shell
export or shared-dictionary mutation from making an unprotected launch appear
exclusive.

### UTR-F04 — active C19 resource contract

The active C19 command in `test/list` must carry both the exclusive directive
and `--require-runner-exclusive`. Its field-initialization flag must occur once,
with `-np 4 -nt 8`. Mode3D uses the caller in addition to the requested eight
temporary workers, so the active configuration has `4 * (8 + 1) = 36`
participants. This matches the 36-CPU allocation seen in the production log and
eliminates the former `4 * (33 + 1) = 136`-participant oversubscription. The
change affects execution parallelism only; C19 mesh resolution, trajectory
sampling, reference observations, thresholds, and pass/fail gates are unchanged.

### UTR-F05 — deferred last-pass persistence

This end-to-end CLI test uses an isolated list containing an exclusive scalar
entry and a two-variant loop. It verifies three persistence contracts:

1. a normal run leaves the list unchanged and atomically saves actual results
   in `.list.last-pass-results.json`; a later `--commit-last-pass` updates only
   the actual passes, preserves the failed loop variant's old provenance,
   removes the cache, and does not execute any test command again;
2. editing the list between the run and commit makes the SHA-256 guard reject
   stale results while retaining both the list and pending cache; and
3. `--update-last-pass` still performs the immediate update and leaves no
   deferred cache.

The cached result also has to retain the scalar entry's `exclusive` provenance.
This makes deferred persistence exercise the current runner schema without
changing command text, expected P/F state, or validation behavior.

These CLI tests intentionally use the runner's normal unlimited per-test
runtime. They therefore also exercise the child-completion heartbeat: on hosts
where an asyncio child-watcher notification is queued without waking the
selector, an exited child must still produce the normal completion record
promptly. The heartbeat never prints or streams child output; stdout/stderr
remain in the per-test log, and the terminal continues to show only launch and
completion progress.

### UTR-F06 — Earth executable identity and transactional status

The C9/C10 shared runtime guard is exercised with temporary executable scripts.
An executable whose `-h` output carries the Earth application banner and command
contract must be accepted and SHA-256 fingerprinted; an executable for another
AMPS application must fail before MPI or physics starts.  The same test verifies
that status JSON is replaced atomically and that a preceding top-level result is
archived recoverably before a new physical invocation.  This is a provenance
gate only: it never substitutes for, edits, or relaxes a scientific comparison.

### UTR-F07 — per-entry time limit

`! runner: timeout SECONDS` gives one entry its own time limit. The parser must
attach it to the following entry only (in either order with `! runner:
exclusive`, and to every expansion of a loop) and reject a missing, non-numeric,
negative, zero, duplicate, detached, or dangling limit. Through the real
subprocess path, the entry's limit must replace a shorter `--timeout`, apply
when the runner has no `--timeout`, and leave other entries on `--timeout`; each
log records which limit applied. The active SA QUICK entries of `test/list` must
be exclusive and have a limit above the nightly's `--timeout 5400`.

### UTR-F08 — an interrupted runner still reports every entry

The AMPS nightly stops a runner at its time limit with SIGINT. The production
runner is started as a subprocess on a three-entry list with `-j 1` and is
signalled (SIGINT, and separately SIGTERM) while the second entry runs. It must
exit with status 130 and still write its reports: the first entry as it
finished, the second as actual F with `interrupted` (its log says so, and its
command, which runs in its own session, must be gone), and the third as actual
F with `not_started`. Before this, the runner printed "Interrupted by user"
and wrote no report, and the nightly lost every Earth result of the night.

Two more cases (2026-11-23) deliver the signal from inside the running entry's
coroutine. A thin launcher loads the unchanged production runner and wraps
`wait_for_process_completion()`; once the second entry's command has started,
the wrapper signals its own process. Before 2026-11-23 the runner's SIGTERM
handler (and Python's default SIGINT handler) raised `KeyboardInterrupt` in
whatever Python frame was executing. When that frame was the running entry's,
the exception skipped its cancellation handler: the log had only its header (no
"# Interrupted" record, no footer), and the command was not stopped and outlived
the runner. That was the nightly UTestRunner failure of 2026-11-22. The signals
are now delivered by the event loop (`run_all_tests_interruptibly`), which
cancels every running entry and awaits its cleanup. Against the old runner both
new cases fail under Python 3.8 (SIGTERM also under 3.14, where `asyncio.run`
already handles SIGINT); against the current runner all four pass.

### UTR-F09 — a failed entry says why it failed

A report line such as `exit=139` alone does not identify a failure. For every
actual F the runner now records two more fields (JSON, CSV columns
`exit_status` and `failure_cause`, the text reports, and the terminal):

* `exit_status`, the exact decoded status. `-N` means the command was killed
  by signal N. `128 + N` means, by the shell/MPI-launcher convention, that a
  process it ran was killed by signal N. For example, Open MPI `mpirun` returns
  139 when a rank receives SIGSEGV. The status is reported as a convention,
  not a certainty, because a program can exit with such a value deliberately.
* `failure_cause`, one line from the test log, chosen in this order:
  1. the wrapper's own `FAILURE CAUSE:` line;
  2. a launch error;
  3. the first crash message of the most specific kind: an MPI launcher line
     naming the rank and signal, then a UCX/Open MPI/C++-runtime/abort message,
     then a generic one;
  4. the exception line of a Python traceback;
  5. the last FAIL/ERROR line;
  6. the last output line.

  Interrupted, not-started and timed-out entries say so.

The C1-C4 wrappers call `earth_test_runtime.report_amps_failure` when AMPS
fails. It prints and appends to the AMPS log:

* the decoded status;
* the command, working directory and log path;
* host, `ulimit -s`, `ulimit -c` and the `OMP_*` variables;
* every crash message, with its line number, and the log tail.

It ends with the `FAILURE CAUSE:` line. The wrapper still returns the AMPS
status unchanged.

The to-address list no longer calls an interrupted or not-started expected-F
entry an "unexpected pass". Expected-F entries that ended by a signal are listed
separately, because a crash is not necessarily their documented failure. Their
P/F comparison is unchanged. The tests use a recorded Open MPI 4.1.7/UCX
SIGSEGV transcript and a stand-in `mpirun`, so they need neither MPI nor AMPS.
They also check that the runner and the wrapper helper decode statuses and
classify crash messages identically.

### Also run by this suite: C9/C10 storm-extremum timing criteria

`run_test.sh` also runs two offline test files of the shared evaluator
`earth_test_runtime.evaluate_extremum_time()`. They are listed here without a
UTR identifier:

- `srcEarth/test/C9/test_minimum_time.py`: the C9 minimum-cutoff timing
  criterion.
- `srcEarth/test/C10/tests/test_expansion_time.py`: the C10
  maximum-expansion timing criterion.

The tests check that near-equal modeled extrema are treated as ties, so the
-O0 and -O3 builds give the same status. A clear model extremum at the wrong
time still fails. An observed extremum that is not significant (below 2
combined sigma) is `NOT_EVALUATED`, never PASS. Missing inputs are handled
fail-closed. The fixtures are the 2006 December storm values and need no AMPS,
MPI, aacgmv2 or matplotlib. They are run from here so that the main list
executes them without a separate list entry.

The suite does not weaken or replace any scientific validation. Exclusivity is
runner-only scheduling metadata; AMPS commands, numerical inputs, observations,
reference solutions, expected P/F states, and acceptance gates remain unchanged.
