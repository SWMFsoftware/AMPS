# Run all shared-model and native SEP + corona tests

The complete entry point runs both the shared release registry and the linked
executable's native `sep-corona` registry, without individual IDs:

```bash
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0
```

Run this Python command once from the AMPS root; it launches MPI itself. It
builds the shared test binary/adapters, runs every shared gate (including
Stage-12 observation preprocessing and Stage-13/14 release/research Python), then launches the requested live suite. Currently
there are 222 shared and 7 native cases, totaling 229. The actual shared
`--list` and executable `--list-tests` outputs own selection; this driver
maintains no copied ID catalogue or fixed count.

Use `--launcher 'srun -n {ranks}'` on a suitable allocation, `--jobs 16` to
change shared compilation parallelism, and `--timeout SECONDS` for a longer
phase allowance. `--list` discovers both scopes without executing cases.
`--model-only` explicitly omits live AMPS and labels the report shared-only.
`--require-no-skips` makes any SKIP produce a nonzero exit.

The supplied active-tube deck selects schema 4, analytic Parker background and
SWCME shock. Its native passes inspect host initialization/corridor/species/
storage/MPI and shared readiness contracts. They do not qualify an installed
PFSS/SCS coronal provider or Stage-11 mover; current host authority enums lack
those selectors. Every combined row retains `shared-model` or `native-amps`
scope. A full coronal runtime adapter still needs its own live qualification.

## Live progress

The runner prints START/END records for discovery, shared compilation, shared
tests and native MPI tests. Build/test output streams as the child produces it,
instead of appearing only after the phase exits. Python test children inherit
unbuffered output. Each recognized result advances a phase-local counter, with
the discovered total, percentage and elapsed seconds:

```text
[progress] shared-model-tests completed=200/222 (90.1%) elapsed=5.7s result=SHEATH3D01 PASS
```

Every 15 seconds, including during quiet compilation or AMPS initialization,
a RUNNING record reports elapsed time and seconds since the last child output.
Append `--progress-interval 5` for five-second updates. The interval must be
positive and finite. Test percentages describe completed cases, not mesh
construction or physical time advancement; native initialization can remain
at zero completed cases while its output and heartbeat continue.

Console results are provisional until the fresh JSON report is checked. The
runner announces verified case counts and aggregate report publication, then
prints the existing final summary. Each subprocess log retains its original
bytes without the runner's progress messages or duplicate output. A timed-out
phase retains its diagnostics and is reported as ERROR.

This progress change is Python-only: replacing the runner needs no AMPS
rebuild. The earlier provider-report changes still require an updated linked
executable when upgrading from a source version predating those changes.

## Aggregate reports and failure handling

`test_output/coupled-sep-corona/summary.json` and `junit.xml` contain individual
results/totals. Fresh `runs/TIMESTAMP-ID/` directories retain original reports,
discovery/build/execution logs and native artifacts. An old passing JSON can
never satisfy a crashed invocation. Commands, durations, executable/deck
checksums, configured provider names, MPI ranks and completed steps are stored.

The final `failed_test_summary` prints every FAIL/ERROR with its scope, ID,
short diagnostic reason and absolute file paths. Infrastructure errors are
listed even if registry discovery failed before any tests were known. Strict
SKIPs are listed separately when `--require-no-skips` makes them fail the run.
An all-passing run still prints the log directory and phase-log paths.

| File beneath the selected output directory | Contents |
| --- | --- |
| `failures.txt` | Latest human-readable failure/log summary |
| `runs/TIMESTAMP-ID/failures.txt` | Summary for that invocation, preserved across later runs |
| `runs/TIMESTAMP-ID/failures/shared-model/ID.log` | Complete captured output of a failed shared case |
| `runs/TIMESTAMP-ID/failures/native-amps/ID.log` | Failed native case message, metrics and artifact references |
| `runs/TIMESTAMP-ID/shared.log` | Complete shared-test console output |
| `runs/TIMESTAMP-ID/native.log` | Complete native MPI console output, including initialization |
| `runs/TIMESTAMP-ID/shared-build.log` | Shared compilation output |
| `runs/TIMESTAMP-ID/shared-discovery.log`, `native-discovery.log` | Registry discovery output |
| `runs/TIMESTAMP-ID/shared/results.json`, `native/native.json` | Original case reports |

The shared case diagnostic is not clipped when its console reason is shortened.
Native per-case diagnostics reference the MPI suite log; they do not fabricate
separate per-case MPI stdout. Failed aggregate rows include `diagnostic_log`,
and aggregate JSON includes the run-specific `failure_summary` path. Log paths
are printed only for files that were generated. Discovery/build failure may
prevent later execution logs or source reports from being created.

For older reports created before this summary feature, identify failed cases
and print their stored diagnostics without rerunning the suite:

```bash
python3 -c 'import json; d=json.load(open("test_output/coupled-sep-corona-all/summary.json")); print("logs:",d["run_directory"]); [print(r["scope"],r["id"],r["status"],r.get("message") or r.get("output") or "",sep="\n") for r in d["results"] if r["status"] in ("FAIL","ERROR")]'
```

Adjust that JSON path to the selected `--output-dir`. If all seven native cases
pass but the aggregate reports one FAIL, that failure belongs to the shared
suite; its ID and captured assertion/traceback are in `summary.json`, with the
original output in the run's `shared.log` and `shared/results.json`.

If `ARCHSCCM01` reports a UTF-8 decode error in the neutral `sep_*` scan,
upgrade the architecture audit and its regression entry point. The old glob
included adjacent compiled objects such as `sep_*.o`; the corrected audit
decodes only regular C/C++ sources/headers. Genuine invalid UTF-8 sources
remain failures with named paths. Forbidden dependencies, archive symbols and
the external public-header consumer are still checked. No AMPS rebuild is
needed for this Python audit fix. To rerun just the affected shared gate from
the AMPS root:

```bash
python3 src/models/sep_coronal_cme/test/run_tests.py --test ARCHSCCM01
```

Every discovered case must appear exactly once. Missing/duplicate/extra cases,
invalid schemas/statuses, wrong MPI ranks, absent provider identities and
report/process exit disagreement are ERROR. A native suite with both PASS and
FAIL rows retains both statuses; its exit is reconciled with all rows together.
A failed phase does not suppress the independent phase. Discovery errors are
also exposed in JUnit, preventing surviving shared PASS rows from hiding a
missing native scope. Build errors create explicit unexecuted ERROR records.

Exit 0 means no FAIL/ERROR (and no SKIP with `--require-no-skips`), 1 means
failed cases or a strict SKIP policy, and 2 means execution/evidence errors.
Increasing `--test-steps` changes only the host horizon, not suite membership.
The aggregate total is shared verification plus native host evidence; it is
not a claim that every shared test executes through the live AMPS mover.

## Native-only executable command

After overlaying the updated source package and rebuilding `amps` with
`srcSEP3D`, run this command from the AMPS root:

```bash
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/coupled-sep-corona/native.json --artifact-directory test_output/coupled-sep-corona/artifacts
```

This selects only descriptors registered in the native `sep-corona` suite.
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
boundary/geometry checks), native JSON provider/active-region evidence, and 29
SEP3D Python runner tests (including ten aggregate-runner regressions).
The real shared aggregate passed 206/206. These are component checks, not a substitute
for executing the linked command above.

Stage-13 release machinery and Stage-14 research/protocol cases are discovered
from the same shared registry. A model-only run now verifies 222 canonical
records. These are portable software checks; the production release profile
separately requires coronal runtime, clean application/core, MPI/convergence,
D1--D10 and registered campaign evidence. Synthetic EVT/XMD/SLM checks do not
mean that an observational campaign passed. See the shared Stage-13/14 guides.

Reserved campaign records EVT3D01, XMD3D01 and SLM3D01 now retain passing
synthetic contract verification but report SKIP for absent actual campaign
evidence. Shared full selection is 222: 219 PASS, 3 SKIP, 0 FAIL. The
aggregate preserves shared SKIPs and enforces `--require-no-skips` across both
scopes. The Stage-13 baseline remains 209 PASS.
