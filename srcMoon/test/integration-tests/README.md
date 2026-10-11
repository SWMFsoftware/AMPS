# M0 linked integration tests

Each `I01`–`I10` directory is a frozen test ID with an individual runner and
machine-readable contract. The runners call the actual configured Moon
executable with one MPI rank and one OpenMP thread. The executable records raw
SI quantities after full `amps_init()`; `integration_common.py` evaluates
independent analytical, invariant, or high-resolution numerical references.

Run a correctly configured executable with:

```sh
python3 srcMoon/test/integration-tests/I02/test.py \
  --amps ./amps --output-dir test_output/srcMoon/M0
```

Every run preserves `stdout.log`, `stderr.log`, `receipt.txt`, `result.json`,
and `run_manifest.json` below `<output-dir>/IXX/`. `PASS`, `FAIL`, `ERROR`, and
`SKIPPED` have the meanings defined in the top-level test README. A kernel or
preflight PASS never replaces a linked-runtime PASS.

I05 requires the no-SPICE neutral-Na build. I08 requires the dedicated
no-SPICE `Na,NA_PLUS` fixture. The other dynamic tests require the DE421 SPICE
build. I04 consumes the canonical `make -j test_Moon` zero-byte diff. The
campaign runner handles these configuration boundaries and rebuilds generated
trees rather than mixing object files from different fixtures.

