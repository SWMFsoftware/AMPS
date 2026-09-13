# srcSEP standalone component tests

Step 1 provides one catalog and result contract for standalone component tests.
The catalog lives in `component_tests.cpp`; generic deterministic selection,
metadata validation, timing, status handling, and output formatting live in
`util/sep_test_registry.*`.  The production CLI remains `util/sep_cli.*`.

## Running tests

Use the linked executable directly or the Make wrappers:

```sh
# Dependency-light parser/registry contract test; does not need AMPS or MPI.
make test-cli-unit

# Dependency-light immutable-background and single-clock contract tests.
make test-state-unit

# Discover the catalog without initializing the model.
make test-list

# Run one stable ID, a group, or all bounded routine tests.
make test-case CASE=DXX01
make test-group GROUP=parker
make test-turbulence
make -j test
```

The makefile expects the linked executable at `../amps`.  If the enclosing AMPS
application writes it elsewhere, add
`SEP_EXECUTABLE=/absolute/path/to/executable`.  `CASE` and `GROUP` are required
for their respective targets, and any child failure is returned by Make.

Equivalent CLI examples are:

```sh
../amps --list-tests
../amps --test DXX01
../amps --test=DXX01 --test=TURB01
../amps --test-group parker
../amps --all-tests
```

`make -j test` intentionally sequences shared-state component execution even
when Make is given `-j`, then runs the embedded SWCME suite using its own native
parallel target.  Expensive extended tests are not part of `--all-tests`.

## Focused Step 2 state and clock tests

`test/run_step2_tests.sh` builds the exact production
`util/sep_background_snapshot.cpp` implementation in a disposable directory
with C++11, `-Wall -Wextra -Werror`, and pthread support. It verifies:

- `BKG01`: required metadata, validity domains, read-only SWMF ownership, and
  compile-time snapshot non-assignability;
- `BKG02`: one snapshot per particle read phase, mover acquisition, stale-time
  rejection, and no publication while particles are active;
- `BKG03`: cross-provider overwrite rejection and an explicit, provenance-bearing
  SWMF-to-local handoff with a new field-line generation;
- `BKG04`: concurrent scheduler threads acquire the identical const generation;
- `BKG05`: deterministic configuration fingerprints change when background
  configuration changes;
- `CLK01`: production code contains no standalone elapsed-time/launch clock and
  reads `PIC::SimulationTime::Get()` only through the runtime clock adapter.

These checks need neither AMPS nor MPI. The linked native executable is still
required to prove the snapshot boundary around the real `PIC::TimeStep()` in a
standalone and SWMF-coupled run. See
[../BACKGROUND_STATE.md](../BACKGROUND_STATE.md) for the full runtime contract.

## Registered tests

| ID | Group | Class | Initialization | What is asserted | State/artifacts |
|---|---|---|---|---|---|
| `DXX01` | `diffusion` | routine | field-line model | `GetDxx` agrees with the constant-coefficient analytical result and the existing million-panel independent quadrature at relative tolerance `1e-5`. | Temporarily replaces and restores the pitch-angle diffusion function pointer; no artifact. |
| `FTE01` | `transport` | routine | field-line model | The focused-transport mover follows the expected field-line displacement while preserving velocity in a static-plasma fixture, using the existing `1e-2`/`1e-5` checks. | Uses fixed registry seed 1002, restores vertex data and diffusion pointer, deletes its particle, and clears test lists. |
| `PARKER01` | `parker` | routine | field-line model | Parker convection keeps the line coordinate stationary in the fixture and matches the analytical density-driven momentum update within `1e-5`. | Uses fixed registry seed 1001, restores vertex data and diffusion pointer, deletes its particle, and clears lists. |
| `PARKER02` | `parker` | extended | field-line model | Four million legacy stochastic trials produce at least one in-range displacement sample and rank 0 successfully writes the histogram.  This is an execution/output assertion, not yet a Gaussian-shape validation. | Uses fixed registry seed 1003; rank 0 writes `dxParker.dat`; particle/lists are cleaned. |
| `SCAT01` | `scattering` | extended | field-line model | Currently reports `SKIP`: the legacy return-probability diagnostic has no approved reference/tolerance and contains a singular zero-energy case. | Registry path performs no mutation. Historical TestManager may write `rmax-E=...` and `time-E=...` files. |
| `TURB01` | `turbulence` | routine | none | The production 1-AU helper equals the independently evaluated magnetic-pressure closure `delta_B^2/(2 mu_0)` within 32 machine epsilons. | Pure deterministic calculation; no RNG, model state, or artifact. |

The list printed by `--list-tests` is authoritative and also includes supported
build modes, seed policy, and state/isolation notes.  Entries are sorted by ID
regardless of their construction order.

## Results, exit codes, and MPI

Each test returns one of `PASS`, `FAIL`, `SKIP`, or `ERROR`, plus a message,
duration, optional seed, metrics/tolerances, and artifact paths.  Exit status is:

- `0`: every requested result is `PASS` or `SKIP`;
- `1`: one or more scientific/test assertions are `FAIL`;
- `2`: one or more tests encounter `ERROR`.

In MPI execution every required rank calls the adapter.  Status severity and
duration are reduced deterministically; the most severe rank status and maximum
duration are reported.  Only rank 0 prints the registry summary.  Shared output
from `PARKER02` is also root-only.  Legacy test bodies may still print detailed
per-rank diagnostics; their removal requires a later test refactor.

Selectors are case-insensitive and de-duplicated.  Unknown IDs/groups, missing
values, `--list-tests` combined with execution, or `--all-tests` combined with
explicit selectors fail before model initialization.  `--help` and
`--list-tests` are pre-initialization success paths.  A new execution selector
always exits before the production loop.

## Focused Step 1 contract tests

`test/run_step1_tests.sh` builds in a temporary directory using `-Wall -Wextra
-Werror`.  It compiles the exact production CLI source with
`SEP_CLI_PARSE_ONLY`, which excludes only the function that writes AMPS globals,
and exercises:

- `CLI01`: help/list mode selection;
- `CLI02`: both `--test` forms and exactly-once execution;
- `CLI03`: deterministic ordering, de-duplication, groups, routine/extended;
- `CLI04`: malformed/conflicting/unknown selection rejection;
- `CLI05`: result exit propagation and frozen legacy TestManager parsing;
- registry completeness: required metadata, callbacks, and unique IDs.

The linked-host targets are still required for final evidence that early exits
precede real AMPS initialization and that test-only execution never enters the
production loop.  A source-only archive cannot substitute static scans for that
runtime evidence.

## Adding a component test

1. Write an adapter returning `SEP::Testing::Result`; never call `exit()` or
   infer PASS from printed output.
2. Add a descriptor in `ComponentTestRegistry()` with a unique stable ID,
   nonempty name/group/description, exact initialization level, supported build
   modes, runtime class, seed policy, and state-isolation contract.
3. Return metrics with comparison operators, tolerances, and units.  List every
   created artifact and make its shared writer MPI-root-only.
4. Restore production globals, function pointers, field-line data, particle
   lists, output settings, and deterministic random state with RAII where
   available.  Otherwise isolate the test and document the limitation.
5. Add focused selection/status tests, run `make test-cli-unit`, then run the
   linked individual/group target and the complete bounded `make -j test` suite.

Self-consistent Alfvén turbulence remains a production provider/evolution
subsystem, not a particle mover.  Its registry group is only a component-test
entry point and does not replace the later dedicated turbulence verification and
validation campaign.
