# SWMF-coupled field snapshots, trajectories, and directional access

## Roadmap Step 9: coherent live state, freeze, export, and replay

Before every requested backward-product callback,
`ReadyForBackwardProductCalculation()` collectively verifies that the mesh, first
coupling receive, magnetic-field offset, and bulk-velocity offset are available on all
ranks.  The callback then also checks that ranks agree on the PT simulation time,
coupler call counter, domain, offsets, and electric-field mode.  Time must be finite,
non-negative, and not older than the previously accepted state; equality is allowed for
an exact retry/restart. The absolute UTC
epoch is the configured reference epoch plus that authoritative PT time; no product is
allowed to apply the offset a second time. A coupled build without SPICE may process
only the zero-offset reference state; a nonzero PT time fails rather than fabricating
calendar/leap-second arithmetic.

`PrepareGlobalSWMFCoupledMagneticFieldForCutoff()` assigns canonical AMR block IDs and
packs exactly one value for every used interior owner cell.  It reads B and plasma bulk
velocity directly from the live owner buffer, converts neither quantity because the
coupler contract is already tesla and metre/second, and rejects missing, duplicate, or
non-finite contributions.  After MPI assembly it rereads each local owner cell and
compares it component-by-component with the compact value.  This direct parity gate
detects a wrong data offset, accidental ghost-cell use, or a source buffer that changed
during assembly even when all ranks would otherwise reduce to a plausible array.

Publication records source `PIC::CPLR:SWMF`, GSM coordinates, SI units, domain,
authoritative simulation time, absolute epoch, AMR `mesh_revision`, complete
`content_fingerprint`, and immutable field-generation ID.  Every rank constructs the
portable representation independently and must agree on the fingerprint.  Cutoff,
directional access, flux, and spectrum then run inside a frozen-field lease; modifying
the global compact arrays while the lease is active is a fatal lifecycle error.  A
receive that arrives while products are running is left for the next scheduler
callback, so one output batch can contain only one MHD epoch.

The released default is magnetic-only:

```text
SWMF_DERIVED_ELECTRIC_FIELD OFF
```

Bulk velocity remains in the compact state and exported file, but E is zero and marked
unavailable to the particle solver.  `EXPERIMENTAL` is the only accepted opt-in; it
derives and verifies the explicit ideal-MHD convention `E=-u x B`.  This label is part
of the content identity and replay contract.  It does not enable a time-dependent
electromagnetic characteristic: the existing mover gates still reject electric-field,
time-dependent, and physical-backward-time requests.

For reproducible offline diagnosis, configure:

```text
SWMF_SNAPSHOT_EXPORT        T
SWMF_SNAPSHOT_EXPORT_PREFIX swmf_field_snapshot
```

Rank zero writes the exact frozen generation before any product is evaluated.  The
strict v1 CSV contains all metadata plus canonical owner cells `(block,i,j,k,x,B,u)`;
diagnostic E is reconstructed from the declared convention rather than stored twice.
Rank-zero export/status I/O success is broadcast before any rank continues, preventing
a local filesystem error from leaving peers blocked in a later solver collective.
Standalone Mode3D replays that file with `FIELD_MODEL SWMF_SNAPSHOT`; it rejects a
schema, unit, frame, epoch, mode, domain, topology, cell-centre, mesh-revision,
fingerprint, or ID mismatch instead of interpolating or silently using the last valid
state.  A fail-closed status JSON is written for every attempted coupled batch.

The phenomenological Step-7 models remain compile-time isolated from this path.  The
shared Step-4 through Step-6 trajectory/access/product kernels are unchanged, as are
all existing validation thresholds.  See `../examples/standalone_step9_swmf_replay.in.template`
for replay controls and `../test/USWMFSnapshot/README.md` for the Step-9 tests.

Roadmap Step 6 requires no SWMF-specific spectrum implementation. After the current
coupler state is frozen and published, `Mode3DForwardSWMF.cpp` calls
`Earth::Mode3D::RunDensityAndFlux()`. That path uses the same
`BoundaryProducts.h` kernel as standalone Mode3D and gridless runs, including explicit
MeV/particle or MeV/nucleon coordinates, log-intensity time-table selection, boundary
uncertainty, normalized/raw PAD and spatial modes, differential spectra, density,
omnidirectional and planar flux, channels, detector rates, and unresolved bounds.

This is an instantaneous quasi-static fold at each published SWMF snapshot. It does
not yet transport a distribution across changing coupled states or enable electric
acceleration. Production Step 6 deliberately selects the static-magnetic
`J_local=A*J_boundary` mapping; the general `j/p^2` mapping remains a tested kernel for
the later electromagnetic mover. Time-history products remain later-roadmap work.

The coupled output schemas and append-only compatibility rules are those documented in
`../3d/README.md`. Validate the backend-independent product physics with
`./srcEarth/test/UBoundaryProducts/run_test.sh`; retain the existing coupled and C/F
tests for field assembly, trajectory access, and observation comparisons.

Assembly rejects unavailable buffers, incomplete/duplicate ownership, and non-finite
owner or reduced values. Those errors are not mapped to forbidden particle access.
Step 3 does not change the coupled field values, AMR interpolation, trajectory mover,
termination policy, or test thresholds.

Roadmap Step 4 uses that published generation through the same backend-neutral
`TrajectoryRequest`/`TrajectoryResult` contract as standalone Mode3D and gridless.
Every request fingerprint is checked against the active coupled snapshot before a
trajectory starts, and every result records the verified fingerprint, mover, backward
convention, exact termination, retry/extension provenance, and optional complete outer-
boundary phase-space state. A coupled snapshot replacement therefore cannot be mixed
silently into an in-progress cutoff/flux batch.

The released Step-4 backtracer remains frozen magnetic-only. Electric or explicitly
time-dependent characteristics and `PhysicalBackwardTime` requests fail fast; even an
explicitly experimental-E snapshot is not permission to use the static antiparticle
shortcut.

Roadmap Step 5 is coupled through the same Mode3D `DIRECT_ACCESS` driver; there is no
SWMF-specific reconstruction path. For each current quasi-static SWMF snapshot, the
writer saves exact solid-angle weights, the three-state `A(E,Omega)` classification,
complete allowed exit states, and adaptive error/work-limit metadata. The snapshot
fingerprint in each trajectory result is checked against the generation published for
that SWMF callback, preventing access rows from two coupling times from being mixed.

The Step-5 product remains instantaneous/quasi-static. It maps an incident boundary
distribution through one magnetospheric state; it does not yet model acceleration,
loss, or trapping evolution across a time sequence of SWMF snapshots. Adaptive
refinement controls and the output schema are identical to standalone Mode3D and are
documented in `../README.md` and `../3d/README.md`.

The dependency-free shared contract tests are:

```bash
./srcEarth/test/UTrajectoryCore/run_test.sh
./srcEarth/test/UDirectionalAccess/run_test.sh
./srcEarth/test/USWMFSnapshot/run_test.sh
```
