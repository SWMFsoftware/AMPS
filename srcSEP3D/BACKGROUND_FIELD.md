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

## Analytic Parker provider

The magnetic field is evaluated in Cartesian form:

\[
\mathbf B=C r^{-3}\left[\mathbf x-k(r-r_0)
(\hat{\mathbf a}\times\mathbf x)\right],
\]

with the implemented sign expressed directly in
`bg_parker.cpp` as `x - spiralFactor*(axis × x)`. This avoids singular
spherical basis vectors on the rotation axis. The provider differentiates the
same expression analytically for `grad(B)` and derives focusing and curvature
from it. The solar wind is radial with analytic
`grad(U)=V/r*(I-r_hat r_hat)` and `div(U)=2V/r`.

Density scales as `r^-2`; pressure is `n k_B T`; Alfvén speed is
`|B|/sqrt(mu0 m_p n)`. `Prepare(epoch)` freezes only validity/generation
metadata; the analytic configuration is immutable.

The curve and tangent are not duplicated here. `core/parker_geometry.cpp`
owns the same source longitude, source colatitude, source radius, wind speed,
rotation rate, and rotation axis used by Phase-M tube refinement. The provider
computes the unsigned tangent first and applies magnetic polarity only to the
field vector. Thus `polarity=-1` reverses `B` and the focused-transport pitch
orientation while leaving the mesh centerline and tube distance unchanged.
`CFG3D04` verifies both alignments at the same physical point.

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

The standalone driver configures analytic authority. During `amps_init()`, the
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
| `BGP3D04–06` | focusing derivative, radial-wind derivatives, finite poles |
| `SNAP3D01–02` | required fields and finite-value policy |
| `SNAP3D03–05` | unit equivalence, epoch consistency, atomic failed update |
| `SNAP3D06` | bracketed time interpolation and extrapolation rejection |
| `SNAP3D07–08` | per-sample batch status and coordinate-frame rejection |
| `R3D03` | requested/filling/staged/failed transitions and collective atomic publication |

Run `test/run_tests.py --suite phase-b --rebuild`.
