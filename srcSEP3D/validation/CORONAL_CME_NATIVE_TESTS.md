# Native AMPS tests for `sep_coronal_cme`

## Purpose and scope

`srcSEP3D` is a production SEP application with native checks of shared
`src/models/sep_coronal_cme` initialization/API contracts. They observe the
selected production providers and do not install a coronal background. The
supplied decks select analytic Parker/SWCME; current host selectors do not
activate PFSS/SCS or Stage-11 providers. The native test mode does not substitute a
small mock mesh or a second physics driver. It parses one complete srcSEP3D
input deck, constructs the ordinary immutable configuration, executes
`amps_init_mesh()` and `amps_init()`, optionally executes the requested number
of ordinary `amps_time_step()` calls, and then observes the joined application
state without mutating it.

This distinction is important. The standalone model tests prove equations,
units, parsers, and deterministic algorithms without AMPS. The tests in this
document inspect generic host contracts on the real generated species table,
AMR allocation, MPI decomposition, center-node storage, spherical boundary,
source provider, particle numerics, halo exchange, and Tecplot writer. Neither
class of evidence replaces the other. Native initialization PASS is not full
coronal-provider or Stage-11 mover integration evidence.

Normal production commands are unchanged:

```bash
mpiexec -n 8 ./amps --input srcSEP3D/examples/sep3d_analytic_parker.in
```

## Executable CLI

The linked executable provides allocation-free discovery:

```bash
./amps --list-tests
```

Run every currently registered native SEP + corona host check without IDs:

```bash
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/coupled-sep-corona/native.json --artifact-directory test_output/coupled-sep-corona/artifacts
```

The selector reads native descriptor suite membership directly. It currently
selects SCCM3D01–07 and includes future `suite="sep-corona"` descriptors after
rebuilding. Discovery displays suite membership. Results print one status per
case plus total PASS/FAIL/SKIP/ERROR counts. See [COUPLED_SUITE_CLI.md](COUPLED_SUITE_CLI.md).

To include all shared-model Stage-0--12 tests in the same invocation, use:

```bash
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0
```

This runs both authorities and combines their reports (currently 206 shared
plus 7 native cases). Rebuild after updating the native report boundary:
`providers` and `completed_steps` must be captured from production state. An
obsolete report fails the aggregate runner rather than receiving an inferred
provider identity.

Run one test or the complete native registry with an explicit immutable input:

```bash
mpiexec -n 4 ./amps \
  --test SCCM3D01 \
  --test-input srcSEP3D/examples/sep3d_analytic_parker.in \
  --test-steps 0 \
  --expect-mpi-ranks 4 \
  --test-json test_output/SCCM3D01/native.json \
  --artifact-directory test_output/SCCM3D01/artifacts

mpiexec -n 4 ./amps \
  --all-tests \
  --test-input srcSEP3D/examples/sep3d_analytic_parker.in \
  --test-steps 1 \
  --expect-mpi-ranks 4 \
  --test-json test_output/native-all.json \
  --artifact-directory test_output/native-all
```

`--test` is repeatable. `--test`, `--test-suite`, `--all-tests` and
`--list-tests` are mutually exclusive. `--test-suite` accepts `sep-corona`
(case insensitive), exactly once; missing/unknown names are usage errors. `--test-input` is mandatory for execution; `--input` may be
given as an alias only when both paths are textually identical. Test execution
cannot be combined with `--dry-run` or `--initialization-only`, because the
native callback must observe a completed real AMPS initialization. A zero
`--test-steps` value tests initialization only; a positive value calls the
production timestep function exactly that many times unless the application's
normal end-of-simulation condition occurs first. The value cannot exceed
`run.maximum_time_steps` in the immutable deck.

`--expect-mpi-ranks` is optional but recommended for MPI evidence. It makes a
launcher error a test failure instead of allowing a one-rank invocation to be
mistaken for the requested matrix row. The output paths are created only by
rank zero after all ranks have joined.

The preferred public runner remains:

```bash
python3 test/run_tests.py --test SCCM3D01 \
  --amps /path/to/amps \
  --validation-input /path/to/reviewed-sep3d.in \
  --validation-launch-prefix "mpiexec -n 4"
```

The runner first requires `--list-tests` to advertise the selected ID, removes
stale JSON, passes the exact deck to the executable, verifies the report
schema and exit status, and hashes the executable, deck, JSON, and artifacts.

## Coronal-CME native cases

| ID | Real AMPS state observed | Shared-model contract |
|---|---|---|
| `SCCM3D01` | mesh, solar sphere, background, turbulence, shock/no-shock state, halo, compiled species, block numerics, observers, output dictionary | `ValidateInitializationForOutput` accepts only the complete ten-bit initialization ledger and a joined background generation |
| `SCCM3D02` | every generated AMPS species, global base weight/time step, and every owner-local block value | `ValidateAllSpeciesWeights` plus positive finite per-species timesteps |
| `SCCM3D03` | exact compiled slot, AMPS chemical symbol, mass, charge, and source-enabled state | `ValidateSourceSpecies` covers every compiled slot; no proton-at-slot-zero assumption is permitted |
| `SCCM3D04` | registered `R_sun` internal sphere, replicated active-use plan, and allocated active leaves | installed full-domain/corridor plan, verified owner-block counts, nonempty physical region, consistent pruning evidence, and solar-interior exclusion |
| `SCCM3D05` | immutable background generation and shock state observed after publication at one joined boundary | background/shock generation coherence; explicit no-shock transport is a complete state, not missing data |
| `SCCM3D06` | mesh, Parker-line, and all per-species AMPS initialization Tecplot products | every product exists, is nonempty, and contains no numeric NaN or infinity token |
| `SCCM3D07` | configuration identity, collective mesh counts, generation, initialization mask, and all species numerics | identical deterministic fingerprint on every rank; optional exact rank-count assertion |

The pre-existing linked `NAT3D01–03/09–12` and `MPI3D01–02` IDs are exposed
by the same executable registry. This removes the former disconnect in which
the Phase-V runner specified a native CLI that the actual srcSEP3D driver
rejected.

## State-capture boundary

`CaptureNativeApplicationState()` is implemented in `main_lib.cpp`, the only
translation unit allowed to see both private srcSEP3D application state and
AMPS. Every rank calls it after initialization or the requested test steps.
It performs only reads and MPI reductions:

1. owner-local background samples are checked with the production completeness
   validator; all magnetic/velocity derivatives and focusing quantities must
   be finite;
2. the installed turbulence provider is evaluated at the identical immutable
   snapshot positions; a physical sample must have finite positive spectral
   support, while an explicitly selected ballistic result must contain finite
   zero-valued sentinels rather than NaN;
3. generated AMPS species identities and global/block-local weights and time
   steps are enumerated over `PIC::nTotalSpecies`;
4. blocks and physical cells are summed collectively, while readiness flags
   use collective logical AND;
5. rank zero checks the already-closed initialization files, then broadcasts
   the result;
6. each rank hashes only immutable or globally reduced quantities, and MPI
   minimum/maximum equality proves cross-rank identity.

The capture never calls `Prepare()`, rebuilds a provider, changes a cell,
injects or resamples a particle, clears sampling, advances a random stream, or
rewrites a Tecplot product. Thus a passing native test describes the state
consumed by production, not a separately prepared test state.

## Species and electron handling

AMPS `SpeciesList` is a compile-time authority. Native tests enumerate all
compiled slots through the same read-only table used by production and never
use `_H_SPEC_`, `_ELECTRON_SPEC_`, or another optional macro as a count or
slot-zero assumption. Mass and charge must be finite physical AMPS values.

The shared source identity record permits `nucleonCount = 0` to mean “not
applicable or unavailable.” This is the physically correct representation for
an electron. A downstream diagnostic that reports energy per nucleon must
separately require a positive nucleon count; source identity validation must
not invent `A=1` for leptons.

## Build and evidence limits

The srcSEP3D makefile builds `sep_coronal_cme` through its independent model
makefile and flattens its objects into `mainlib.a`, because the enclosing AMPS
link consumes application archives and cannot extract objects recursively from
a nested static archive. The production audit requires the shared
`ValidateInitializationForOutput` symbol and exactly one copy of every listed
member.

A source-only checkout can run the portable model and srcSEP3D tests, including
`BLDL3D11`, but it cannot claim native AMPS evidence. `BLDL3D01` reports `SKIP`
when no configured `Makefile.conf`, generated PIC headers, and linked AMPS
executable are present. Final qualification requires running `SCCM3D01–07` on
the target AMPS/MPI toolchain with a reviewed complete deck.
