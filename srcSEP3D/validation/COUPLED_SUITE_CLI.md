# Run all coupled SEP + corona tests

After overlaying the updated source package and rebuilding `amps` with
`srcSEP3D`, run this command from the AMPS root:

```bash
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/coupled-sep-corona/native.json --artifact-directory test_output/coupled-sep-corona/artifacts
```

This selects every descriptor registered in the native `sep-corona` suite.
Today it selects SCCM3D01–07; the command automatically grows with the registry.
It initializes the real distributed AMPS application once, captures state at a
collective boundary, evaluates every selected case, writes the existing native
JSON format and state artifact, and prints per-case results and totals. It
preserves the production Parker corridor, halo, and solar exclusion behavior.

`--test-steps 0` is initialization-only validation. To check advancement, set
an appropriate positive horizon within `run.maximum_time_steps`; use a deck
and horizon suitable for the tests' prerequisites. Suite selection does not
turn a skipped check into a passing check or generate missing physical inputs.
For full-domain initialization, use `sep3d_analytic_parker.in` instead.

`./amps --list-tests` lists the complete native registry and each descriptor's
suite. `--all-tests` runs that entire registry, including general AMPS and
MPI/restart cases. Existing repeated `--test ID` selection is unchanged.
These three execution selectors and `--list-tests` are mutually exclusive.
Unknown/missing suite names, duplicate suite selectors, missing input decks,
and incompatible dry-run/initialization-only switches fail before AMPS setup.

Exit codes: 0 means no FAIL/ERROR, 1 means at least one FAIL, and 2 means ERROR
or invalid invocation. Inspect SKIP counts separately when qualifying a run.
All MPI ranks receive the same test exit code. Particle stepping and full
AMPS/MPI compilation remain to be verified on the configured installation.

## Register future tests

Add a descriptor to `CoronalCmeNativeTests()` in
`validation/coronal_cme_application_test.cpp`, with the fourth field set to
`"sep-corona"`, and implement its production-state evaluator:

```cpp
{"NEW_TEST_ID", "Test name", "Acceptance contract", "sep-corona"},
```

Rebuild the executable. The suite selector traverses the entire native registry;
it has no copied ID list, fixed SCCM prefix/range, or fixed case count. Future
tests requiring new live state must also extend collective state capture and
its evaluator. The separate scientific matrix/global Python runner registries
continue to have their own explicit evidence/prerequisite policies; this
change targets native executable suite selection.

## Build and portable checks

From the AMPS root, use the same compiler/MPI configuration as your existing
build:

```bash
make -C src/models/sep_coronal_cme -j8 lib
make clean
make -j8
python3 srcSEP3D/test/run_native_boundary_regression.py
```

The CLI parser now lives in `runtime/standalone_command_line.cpp`, which is
included in the production/test runtime object lists. This permits testing the
actual parser without the omitted SWCME configuration resolver. Existing
configuration loading remains in `configuration_io.cpp`.

Portable validation passed 45 C++ checks (21 suite/CLI and the prior 24 native
boundary/geometry checks), native JSON evidence verification, and the 19
existing SEP3D Python runner tests. These are component checks, not a substitute
for executing the linked command above.
