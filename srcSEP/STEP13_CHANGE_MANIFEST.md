# Step 13 change manifest: acceptance CLI and evidence

Step 13 makes the test CLI a fail-closed acceptance boundary rather than a
console-only diagnostic launcher.

## Implemented behavior

- `--test-json <path>` writes deterministic `srcsep-component-tests-v1` JSON.
- `--test-junit <path>` writes CI-compatible JUnit XML.
- Both options require an executing `--test`, `--test-group`, or `--all-tests`
  selector and report-write failures return process status 2.
- MPI ranks execute each selected callback; the most severe status and longest
  duration are reduced before root writes evidence. Root also gathers every
  rank's pre-reduction status, message, seed, configuration, metrics, and
  artifact paths, so a non-root failure remains diagnosable.
- A callback returning `PASS` with a positive or non-finite
  `assertion_failures` metric is forcibly converted to `FAIL`.
- Complete seed, configuration, metric, artifact, message, and duration data
  are retained in both report formats.
- `BG01` supplies standalone analytic and SWCME snapshots; `BG02` supplies a
  mock, read-only SWMF import and explicit local-evolution handoff.
- `CROSS01` compares `fte-dmumu` and `fte-mfp` in the exact matched ballistic
  limit. `CROSS02` round-trips the isotropic closure shared by Parker `kappa`,
  focused-transport `Dmumu`, and event-driven `lambda`.

## Files

- `util/sep_test_registry.*`: result audit and structured report writers.
- `util/sep_acceptance_cases.*`: reusable background and cross-mover fixtures.
- `util/sep_cli.*`, `main.cpp`, `tests.h`, `component_tests.cpp`: selector,
  report-path, MPI, and production-registry integration.
- `test/step13/test_acceptance_cases.cpp`, `test/run_step13_tests.sh`: strict
  source-only ASan/UBSan validation of the exact linked callbacks.

Run `make test-acceptance-unit` and `make test-stochastic-repeat`. The full
linked/MPI gate is available through `make test-mpi` (configure `MPIEXEC`,
`MPI_NP`, and `SEP_EXECUTABLE` as needed); expensive linked diagnostics are
explicitly selected by `make test-stress`.
