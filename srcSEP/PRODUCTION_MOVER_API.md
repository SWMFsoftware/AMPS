# Production field-line mover API

Step 4 establishes one public registry for the three supported production
movers. The registry describes physics and state capabilities independently of
implementation addresses; the PIC callback enters through one validating
adapter.

## Canonical movers

| CLI name | Public enum | Production implementation | Coefficient contract |
|---|---|---|---|
| `parker` | `ProductionMover::Parker` | `ParticleMover_Parker` | `SpatialDiffusionProvider`: `kappa_parallel` [m²/s], `d(kappa_parallel)/ds` [m/s], and provenance |
| `fte-dmumu` | `ProductionMover::FocusedTransportDiffusion` | `ParticleMover_FocusedTransport_Dmumu` | `PitchAngleDiffusionProvider`: `Dmumu` [s^-1], `dDmumu/dmu` [s^-1], provenance, and turbulence-state identity |
| `fte-mfp` | `ProductionMover::FocusedTransportMeanFreePath` | `ParticleMover_FocusedTransport_EventDriven` | `MeanFreePathProvider`: `lambda_parallel` [m], validity, provenance, turbulence identity, and exact `nu=v/lambda` events |

Self-consistent turbulence is deliberately absent from this table. Step 11
selects its source and integrated/spectral authority through
`Turbulence::Configuration`; it remains a coefficient/state subsystem consumed
by these movers, never a fourth dispatch entry. Step 12 changes only the
deterministic reduction of deferred particle-to-wave contributions.

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

The complete native AMPS regression gate must additionally exercise one
particle through each adapter mapping and confirm startup metadata on the linked
executable. That dependency is not present in a source-only `srcSEP` archive.
