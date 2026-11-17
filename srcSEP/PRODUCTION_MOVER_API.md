# Production field-line mover API

Step 4 establishes one public registry for the three supported production
movers. The registry describes physics and state capabilities independently of
implementation addresses; the PIC callback enters through one validating
adapter.

See [STATUS.md](STATUS.md) for current qualification. The registry and retained
coefficient interfaces exist, but the adopted four-spectrum joint scattering
law is still being connected to native particle/wave transactions. Advertised
streaming capability does not establish physical qualification; dependent
production paths reject explicitly while controlled diagnostic tests remain
available.

## Canonical movers

| CLI name | Public enum | Production implementation | Coefficient contract |
|---|---|---|---|
| `parker` | `ProductionMover::Parker` | `ParticleMover_Parker` | `SpatialDiffusionProvider`: `kappa_parallel` [m²/s], `d(kappa_parallel)/ds` [m/s], and provenance |
| `fte-dmumu` | `ProductionMover::FocusedTransportDiffusion` | `ParticleMover_FocusedTransport_Dmumu` | `PitchAngleDiffusionProvider`: `Dmumu` [s^-1], `dDmumu/dmu` [s^-1], provenance, and turbulence-state identity |
| `fte-mfp` | `ProductionMover::FocusedTransportMeanFreePath` | `ParticleMover_FocusedTransport_EventDriven` | Required adopted law: same polarized spectrum-derived joint operator as `fte-dmumu`; derive `lambda_parallel` [m] where defined. The retained `MeanFreePathProvider`/isotropic `nu=v/lambda` path is not its qualified replacement. |

Self-consistent turbulence is deliberately absent from this table. Step 11
selects its source and integrated/spectral authority through
`Turbulence::Configuration`; it remains a coefficient/state subsystem consumed
by these movers, never a fourth dispatch entry. Step 12 changes only the
deterministic reduction of deferred particle-to-wave contributions.

The selected event process carries residual optical depth and draws a branch
mark only when its combined physical hazard expires. A scalar mean free path
does not determine that angular law. If the summed pitch diffusion vanishes on
a positive-width interior interval, the usual fixed-momentum diffusion-limit
integral diverges: local scattering can still occur elsewhere, but a finite
global diffusion mean free path and its Parker reduction are unavailable.
Preserving that gap forbids a numerical scattering floor or isotropic reset.

Native particle state retains all three operator streams, cached normal state,
event history and split serial. Schema 4 is being integrated with an explicit
source-driver window and persisted first-movement duration; only fresh source
or initial-condition APIs may create that history. Phase space, RNG and pending
birth consumption commit together. Completed checkpoints reject unfinished
source windows/births. Build 6b passes BIRTHN01–06 on one and two ranks after
correcting the fixture's explicit static epoch; snapshot ownership assertions
remain intact. Fresh-process and complete source/coupling qualification remain
separate, as recorded in STATUS.md.

These are one-dimensional transport models along field lines embedded in 3-D
space. The field geometry remains three-dimensional. Fully 3-D particle
trajectories, drift/Boris variants, and their Cartesian particle state are not
compiled into srcSEP after Step 5. Direct-wave experimental movers and
historical ambiguous field-line scattering functions were deleted in Step 14.

## Selection and discovery

```sh
../amps --list-movers
../amps --particle-mover parker
../amps --particle-mover=fte-dmumu
../amps --mover fte-mfp
```

A runtime input file (`-input FILE`, or the synonym `--input FILE`) in the
shared SEP content schema 6 (`../src/INPUT_FILE.md` section 5) sets the mover
with the required key `application.sep.transport`, whose spelling is shared
with srcSEP3D:

```text
#begin application: sep
  schema_version = 6
  transport = parker      ! parker | focused-diffusion (fte-dmumu) | focused-scattering (fte-mfp)
  #begin turbulence
    ...                   ! the 28 srcSEP turbulence keys (README.md)
  #end turbulence
#end application
```

The retired srcSEP key `application.sep.particle_mover` (schema 4) is an
unknown key, and the `fte-*` names are not accepted in the file.
Precedence is `--particle-mover` > `application.sep.transport` > default
`fte-dmumu` (the default applies only when no input file is read); with an input file and no `--particle-mover`, the
mover-dependent CLI checks run after the file value is applied
(`SEP::Util::CLI::ValidateMoverDependentOptions`). At `amps_init()` srcSEP
registers `SEP::Mover::DispatchProductionMover` as generic PIC's
user-defined mover (the deck sets
`_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_ = _PIC_MODE_OFF_`); it dispatches to the
selected mover.

`--help` and `--list-movers` advertise the complete accepted set. Step 14 ended
the transition interval, so former short/legacy aliases now fail before model
initialization. Names such as `coupled-fte`, `focused-transport-wave-scattering`,
`parker-mean-free-path`, `mean-free-path-scattering`, `tenishev-2005-fl`, 3-D,
drift, and Boris names are rejected. They differ physically or are ambiguous
and therefore cannot be silently mapped. See
[MIGRATION_MANIFEST.md](MIGRATION_MANIFEST.md) for explicit replacements.

## Capabilities and callback contract

Each descriptor reports whether the mover requires field-line attachment, uses
pitch-angle state, accumulates streaming for the turbulence manager, evolves
wave state directly, and which coefficient contract it consumes. Main-loop
turbulence dispatch queries `CurrentCapabilities()`; it no longer compares raw
function addresses.

All PIC mover calls route through `DispatchProductionMover()`. Before invoking
the selected implementation, the adapter verifies:

- a finite nonnegative time interval and a valid particle record;
- enabled field-line attachment;
- a valid field-line ID and attached segment;
- a finite field-line coordinate; and
- finite parallel and normal particle velocities.

The adapter translation unit additionally requires `_PIC_FIELD_LINE_MODE_` and
`_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_` at compile time. Supported movers
commit directly to a field-line segment temporary list; there is no AMR-node
particle-list branch or Cartesian post-mover scattering pass.

The three shells share `transport_common.*`. It converts physical displacement
to the host field-line coordinate exactly once through `FieldLine::move`, reads
focusing and flow derivatives from one immutable background snapshot, and owns
the only validation/commit/segment-attachment sequence. The coefficient-driven
movers additionally use explicit keyed random streams whose keys exclude
MPI rank and worker identity.

For native integration evidence, the adapter increments a process-local atomic
counter only after the selected concrete mover returns. The value is not a
particle number, flux, or source normalization: a particle that is advanced to
and deleted at a valid boundary still completed one dispatch. Mover selection
resets the count. After stepping stops, the standalone driver reduces the count
across MPI ranks and prints `completed_particle_dispatches`; D03 requires a
positive value for every mover/decomposition/restart case. Fatal preflight or
physics errors do not return and therefore cannot satisfy that requirement.

Startup output records the canonical mover, representation requirements,
turbulence streaming behavior, coefficient contract, and active provider.
Self-consistent Alfvén turbulence remains a separate production subsystem; all
three movers provide manager-consumed streaming terms and none directly owns
wave-state evolution. Parker uses the established pitch-angle-averaged
streaming closure. `fte-dmumu` associates every kick with the provider's
turbulence identity and defers coupling records until the particle survives and
is reattached. Post-step sorting precedes the legacy G+/G- accumulation.

Step 10 routes all three shells through `coefficient_providers.*`. The runtime
registry rejects recursive conversions and source/background ownership
mismatches. The canonical `fte-mfp` adapter replaces the former `fte_mover.cpp`
production object; Step 14 deleted that file and the other legacy monoliths.

## Focused tests

Run `make test-mover-api-unit` for:

- `MOVCLI01`: canonical parsing, help, and exactly three discovery entries;
- `MOVCLI02`: rejection of every retired mover alias;
- `MOVCLI03`: capability reporting and absence of main-loop pointer comparison;
- `MOVCLI04`: rejection of unsupported or ambiguous historical movers.

Run `make test-transport-common-unit`, `make test-parker-unit`,
`make test-fte-dmumu-unit`, `make test-fte-mfp-unit`, and
`make test-coefficients-unit` for the Step 6–10 numerical contracts. Their stable
IDs and tolerances are documented in [test/README.md](test/README.md).

The complete native AMPS regression gate exercises at least one completed
dispatch through each adapter mapping, confirms startup/final metadata, and
compares serial/MPI/restart artifacts from the linked executable. Its campaign
inputs remain site-owned and that dependency is not present in a source-only
`srcSEP` archive. See [NATIVE_INTEGRATION.md](NATIVE_INTEGRATION.md).

## WP34 and WP35 native evidence boundary

`sep_configuration_matrix.*` forms the 90-row mover/coefficient-source/
turbulence-source/coupling product from the production mover registry and
classifies every row before initialization. `make print-configuration-matrix`
prints that exact registry output. Unsupported rows return stable `CFG-*`
diagnostics with an alternative and allocate no particle state.

`sep_evidence.*` defines the observation record a linked harness must return:
registered production mover entered, configured coefficient adapter entered,
common turbulence driver entered, exactly one queue-flush owner, executable and
configuration identity, source generation, thread count, and rank count. The
validator rejects any source-only attempt to claim `native-amps`. Completion of
the actual PIC fixture remains a linked AMPS gate, not a dependency-light test.
