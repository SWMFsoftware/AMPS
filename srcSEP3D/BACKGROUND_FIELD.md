# Phase B Background Providers and Snapshots

Phase B separates ambient-field acquisition from AMPS storage. Providers and
snapshots are AMPS/MPI-free; `main_lib.cpp` is the only layer that copies a
validated sample into an AMPS center-node buffer.

## Complete background sample

Every published sample is SI and contains:

- magnetic field, magnitude, unit direction, and optional analytic gradient;
- `div(b_hat)`, focusing length, and curvature;
- plasma velocity, gradient, divergence, and field-aligned strain;
- number density, temperature, pressure, and Alfvén speed;
- snapshot generation and configuration digest.

Completeness validation rejects invalid flags, zero/non-finite required
primitives, non-unit directions, non-finite advertised derivatives, and a zero
generation. Required values are never replaced by guessed defaults.

## Installation into AMPS native storage

The production executable has two independent center-node layouts. The frozen
srcSEP3D layout is the versioned mover/restart authority, while
`PIC::CPLR::DATAFILE` owns the fields printed by AMPS' native Tecplot callback
and read by native background accessors. Publishing a `BackgroundSnapshot`
does not populate that second layout automatically. `main_lib.cpp` therefore
performs one explicit bridge operation for every validated physical cell:

| Canonical `BackgroundSample` value | Native AMPS value | Definition |
|---|---|---|
| `B` | `MagneticField` | identical Cartesian SI vector [T] |
| `U` | `PlasmaBulkVelocity` | identical Cartesian SI vector [m/s] |
| `numberDensityM3` | `PlasmaNumberDensity` | SWCME electron density [m^-3] |
| `temperatureK` | `PlasmaTemperature` | SWCME proton temperature [K] |
| `pressurePa` | `PlasmaIonPressure` | canonical total thermal pressure [Pa] |
| `divU` | `PlasmaDivU`, when allocated | identical analytic divergence [s^-1] |
| `gradB` | `MagneticFieldGradient`, when allocated | identical row-major `dB_i/dx_j` [T/m] |
| `U`, `B` | `ElectricField` | ideal-MHD `E=-U x B` [V/m] |
| `gradB` | `Current`, when allocated | stationary Ampere law `J=curl(B)/mu0` [A/m^2] |

If AMPS allocates a separate electron-pressure slot, proton-only closure writes
zero because that closure explicitly excludes electron pressure; multi-species
closure writes the exact `ne*kB*Te`. The source configuration currently defines
one canonical solar-wind state, so a DATAFILE layout with `nIonFluids != 1`
fails before output instead of assigning undocumented fluid identities.

Before physical cells are copied, every owner-local native record is zeroed.
Cells inside the inner boundary or outside the spherical outer boundary remain
finite placeholders and are marked `background_valid=0`. After the background
and directional turbulence values are complete, AMPS exchanges block halos;
only then is the native background marked ready. The same bridge and halo
boundary run on subsequent background generations, so AMPS-native accessors do
not lag the immutable srcSEP3D snapshot.

`sep3d-initialization-data.dat` is written only after this bridge, runtime
snapshot publication, turbulence preparation, all-species weight/time-step
initialization, mover-context installation, and optional restart restoration.
The earlier `sep3d-initialization-mesh.dat` remains deliberately geometry-only.

## SWCME-backed analytic Parker provider

The standalone provider is analytic in space, but its ambient physics is not a
second srcSEP3D solar-wind approximation. `configuration_io.cpp` first resolves
the complete canonical SWCME3D layer. `[background.parker]` is a fail-closed
human-review cross-check; after agreement, the provider receives the exact
resolved wind, Parker, density, composition, temperature, and thermodynamic
values from SWCME. `Prepare(epoch)` builds one immutable
`swcme::solarwind::PreparedState` before publishing generation metadata.

SWCME declares total `B1AU` at a reference `sin(theta)`, while the application
cross-check declares positive radial `Br` at an arbitrary reference radius and
keeps polarity separate. The provider exactly inverts the one-AU pitch
normalization, then calls the shared SWCME Cartesian evaluator. The resulting
field can be written

\[
\mathbf B=C r^{-3}\left[\mathbf x-k(r-r_0)
(\hat{\mathbf a}\times\mathbf x)\right],
\]

with the implemented sign expressed as
`x - winding*(r-r0)*(axis × x)`. This avoids singular spherical basis vectors
on the rotation axis. The provider differentiates the exact same expression
analytically for `grad(B)` and derives focusing and curvature from it. The wind
is radial with `grad(U)=V/r*(I-r_hat r_hat)` and `div(U)=2V/r`.

### Density and thermodynamics

`numberDensityM3` is SWCME electron density. It uses the normalized Leblanc,
Dulk & Bougeret law

\[
n_e(r)=C_2/r^2+C_4/r^4+C_6/r^6,
\]

where SWCME scales all coefficients together so the configured density at one
AU is exact. The density reference is deliberately distinct from the arbitrary
magnetic `Br` reference. The positive fourth- and sixth-power terms make the
near-Sun state physically different from the retired pure `r^-2` profile.

The selected canonical SWCME closure determines the remaining plasma state:

- `PROTON_ONLY`: `np=ne`, `rho=mp*np`, and `p=np*kB*Tp`;
- `MULTI_SPECIES`: for `f=na/np`, charge neutrality gives
  `np=ne/(1+2f)` and `na=f*np`; mass density is
  `mp*np+m_alpha*na`, and pressure is
  `kB*(np*Tp+ne*Te+na*Talpha)`.

In both cases `temperatureK` is proton temperature and
`alfvenSpeedMpS=|B|/sqrt(mu0*rho)`. A nonzero alpha abundance therefore changes
the Alfvén speed through mass density rather than being ignored.

The curve and tangent are not duplicated here. `core/parker_geometry.cpp` owns
the same source geometry, wind, rotation rate, and rotation axis used by
Phase-M tube refinement. Polarity reverses `B` and the focused-transport pitch
orientation while leaving the mesh centerline and tube distance unchanged.
`CFG3D04` verifies that identity. `BGP3D07` proves exact one-AU normalization,
the near-Sun Leblanc correction, and multi-species pressure/Alfvén speed.

## Reserved Python interpolation provider

`background.provider = python-interpolator` is recognized by the schema and by
the typed `BackgroundAuthority`/`ProviderKind` vocabulary, but returns
`ReservedFeature` before mesh allocation. There is not yet a reviewed contract
for process lifetime, source-data identity, interpolation order, units, missing
points, timeout, or Python failure propagation. The selection never falls
through to `analytic-parker` or `swmf`.

The future bridge must implement `BackgroundProvider`, not call a script from a
particle mover. At each background cadence it will:

1. freeze epoch, coordinate frame, data identity, and configuration;
2. send owner-local cell positions in bounded SI batches;
3. receive a complete record and per-point status for every coordinate;
4. validate units, finite values, ordering, derivatives, and coverage in
   temporary storage; and
5. enter the existing all-rank staged publication only if every point/rank
   succeeds.

No partial result, stale generation, pointwise subprocess fallback, or guessed
value will be publishable.

## SWMF/AWSoM ambient import

`SwmfAwsomProvider` is an import adapter, not an SWMF reader. The host supplies
a complete candidate with frame, epoch, generation, ownership, and records.
Two input unit systems are supported:

| Quantity | SI input | AWSoM coupling input |
|---|---|---|
| position | m | solar radii |
| magnetic field | T | nT |
| velocity | m/s | km/s |
| number density | m⁻³ | cm⁻³ |
| pressure | Pa | nPa |
| magnetic gradient | T/m | nT/R_sun |
| velocity gradient | s⁻¹ | (km/s)/R_sun |
| curvature, `div(b_hat)` | m⁻¹ | 1/R_sun |

Imported ownership must be read-only. Every record must have the declared
epoch, and the coordinate frame must match the expected
`HCI-like-inertial` frame unless the host has already applied an explicit
transformation. Conversion occurs in temporary records. A failed late field,
mixed epoch, bad unit value, or frame mismatch leaves the prior generation
unchanged.

## Immutable construction and batch semantics

`BackgroundSnapshotBuilder` requires a prepared provider, finite positions,
valid metadata, successful per-point statuses, complete samples, and matching
sample/provider generations. It assigns the output shared pointer only after
all validation succeeds.

`EvaluateBatchDetailed` has transactional per-point output: successful index
`i` replaces `out[i]`; a failure writes only `status[i]` and leaves `out[i]`
unchanged. The aggregate return reports the first failure so a builder can
reject the candidate while preserving the exact failing positions.

## Current/next time buffer

`SnapshotBuffer` accepts a next generation only when provider, ownership,
frame, configuration identity, spatial grid, and sample count match the
current generation and both generation and epoch increase. `SnapshotAt(t)`:

- returns an endpoint snapshot at its epoch;
- linearly interpolates complete fields inside the two-epoch bracket;
- rejects requests outside the bracket;
- rejects interpolation between different provider configuration digests.

There is no temporal extrapolation.

## Production and coupling path

The standalone driver configures SWCME-backed analytic authority. During
`amps_init()`, the
application gathers owner-local cell centers outside the inner sphere, builds
one immutable Parker snapshot, writes every complete field at the Phase-M
offset, and publishes it with `StandaloneAdapter`.

A coupled host configures SWMF authority and calls
`InstallBackgroundSnapshot()` with a complete read-only snapshot in the same
deterministic owner-cell order before `amps_init()`. The application validates
the authority, frame, grid, and fields, writes the same layout, and publishes
through `SwmfAdapter`. Both adapters invoke the same `Runtime` transition from
`MeshReady` to `WaitingForSnapshot` to `SnapshotReady`.

## R03 joined-boundary update transaction

Subsequent coupled or analytic generations use a separate update state:

`Idle -> Requested -> Filling -> Staged -> Idle`.

`RequestSnapshotUpdate` freezes the requested epoch/generation;
`BeginSnapshotFill` authorizes provider work; `StageSnapshot` validates the
complete descriptor without replacing the active descriptor; and
`PublishStagedSnapshot(true)` performs the commit. A provider or rank failure
enters `Failed`. `AcknowledgeSnapshotFailure` returns to `Idle`, leaving the
previous active generation and its validity interval unchanged.

`main_lib.cpp` builds the candidate background and evaluates the matching
turbulence provider at every owner-local physical cell. An `MPI_Allreduce` of
the readiness flag precedes publication. Only after all ranks agree does the
application swap both immutable shared pointers. The AMPS associated-data
bytes are a cache/diagnostic copy; mover resolution reads the immutable active
snapshot, so partially filled next-generation cache bytes are never physics
authority. Snapshot provenance records provider, frame, configuration digest,
epoch, validity interval, and monotonically increasing generation and is
written into every R07 checkpoint.

## Evidence

| IDs | Contract |
|---|---|
| `BGP3D01–03` | divergence-free field, component laws, field-line tangency |
| `BGP3D04–07` | focusing derivative, radial-wind derivatives, finite poles, SWCME ambient closure |
| `SNAP3D01–02` | required fields and finite-value policy |
| `SNAP3D03–05` | unit equivalence, epoch consistency, atomic failed update |
| `SNAP3D06` | bracketed time interpolation and extrapolation rejection |
| `SNAP3D07–08` | per-sample batch status and coordinate-frame rejection |
| `R3D03` | requested/filling/staged/failed transitions and collective atomic publication |

Run `test/run_tests.py --suite phase-b --rebuild`.
