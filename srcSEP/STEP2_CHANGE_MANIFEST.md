# Step 2 archive change manifest

This archive was rebuilt after completing Step 2, **Create an immutable
background snapshot and one simulation clock**. It contains the prior Step 1
component-test registry plus all Step 2 changes.

## New Step 2 files

- `BACKGROUND_STATE.md`
- `STEP2_CHANGE_MANIFEST.md`
- `test/run_step2_tests.sh`
- `test/step2/test_background_snapshot.cpp`
- `util/sep_background_runtime.cpp`
- `util/sep_background_runtime.h`
- `util/sep_background_snapshot.cpp`
- `util/sep_background_snapshot.h`

## Existing files modified by Step 2

- `README.md`
- `main.cpp`
- `main_lib.cpp`
- `makefile`
- `sep.cpp`
- `sep.h`
- `shock_analytical_model.cpp`
- `sw1d.cpp`
- `test/README.md`

## Quick verification after extraction

From the extracted `srcSEP` directory:

```sh
test -f util/sep_background_snapshot.cpp
test -f util/sep_background_runtime.cpp
test -f test/step2/test_background_snapshot.cpp
make test-state-unit
```

The focused test should print:

```text
Step 2 background-snapshot tests: PASS
Step 2 single-clock source contract: PASS
```

The linked `make -j test` gate additionally requires the enclosing AMPS source
tree, `../../Makefile.conf`, MPI/PIC dependencies, and the production executable.
