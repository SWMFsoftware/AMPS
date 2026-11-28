# Native Enceladus startup adapter

This exclusive suite performs a clean native `enceladus` build through the
mandatory repository preflight and imports six unchanged gates from
`srcEnceladus/test/run_native_startup.py`. The gates cover CLI/default-deck
validation, rejected inward stripe geometry, rejected inaccessible plasma
cache, and generated pointer-mover/particle-processing hooks. They stop before
MPI initialization. Six additional groups from
`srcEnceladus/test/run_case_native_tests.py` validate the actual case CLI,
one/four-rank input agreement and failure paths, all four source bases with
H2O/CO2 transport, both UVIS operators, noncubic runtime AMR and negative
domain/capability gates. Each synthetic transport run advances four steps
and closes an exact global 40-primary-particle ledger; these are integration
checks, not a converged E3/Cassini comparison.

```bash
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001
```

The suite retains the configured executable and its hash, build transcript,
local receipts, per-command MPI/run logs, source snapshots and startup/case
transcripts. Its fixture publication is retained beside the binary so the
compiled model-data root remains inspectable. The fixture generator uses
only Python's standard library; no real fitting/observational data is needed.
Machine `tools.mpiexec` and `execution.mpi_rank_flag` configure launchers.
Missing MPI marks dependent groups SKIP; it never upgrades them to PASS.

The owner confirms the root and idle native state, removes only generated
`build`, runs `Config.pl -application=enceladus -model-data-path=FIXTURE_ROOT`,
generates hooks with `ampsConfig.pl -input enceladus.input -no-compile`, then
compiles through top-level `env MAKEFLAGS=-j16 make amps`. It holds and restores
existing `srcEnceladus/test_output` outside the copied source tree during
configuration. It changes the checkout-wide generated tree and therefore
always runs alone. After qualification, an operator can reconfigure the data
root for real case bundles; see `srcEnceladus/RUNTIME_INPUT.md`.

The exact 12 stable IDs and their purpose are in
`srcEnceladus/test/README.md` and `src/test/TEST_CATALOG.md`. Agents adding a
test must update the registry, traceability and generated catalog in the same
change. The historical suite ID retains `NATIVE-STARTUP` for compatibility.
