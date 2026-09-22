# UTestRunner: deferred last-pass regression tests

Run this dependency-free Python suite from any directory:

```bash
./srcEarth/test/UTestRunner/run_test.sh
```

The suite invokes the public `srcEarth/test/test_runner.py` command-line
interface against isolated temporary test lists. It checks three contracts:

1. A normal run leaves the list unchanged and atomically saves scalar and
   loop-expanded results in `.list.last-pass-results.json`. A later
   `--commit-last-pass` updates only actual passes, preserves the failed loop
   variant's previous provenance, deletes the cache, and does not execute any
   test command again.
2. Editing the list between the run and commit causes the SHA-256 guard to
   reject the stale results. The list and pending cache remain available and no
   command is re-executed.
3. `--update-last-pass` still performs the immediate update and leaves no
   deferred cache.

The fixture includes one scalar pass plus a two-variant loop whose outcomes are
an expected pass and an expected failure. This verifies that provenance follows
the saved **actual** result independently for each source entry and loop
variant.
