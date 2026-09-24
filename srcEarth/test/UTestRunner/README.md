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

The suite does not weaken or replace any scientific validation. Exclusivity is
runner-only scheduling metadata; AMPS commands, numerical inputs, observations,
reference solutions, expected P/F states, and acceptance gates remain unchanged.
