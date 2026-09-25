# SWMF-coupled field snapshots, trajectories, and directional access

Before each requested backward-product callback,
`PrepareGlobalSWMFCoupledMagneticFieldForCutoff()` gathers authoritative owner-cell
magnetic field and plasma velocity into the compact Mode3D arrays. Electric field is
derived cell-by-cell as

```text
E = -v × B.
```

The arrays are published with one `FieldProvider.h` metadata record containing source
`PIC::CPLR:SWMF`, the configured reference epoch plus authoritative PT simulation-time
offset, GSM/SI units, domain, derived-E interpolation mode, and a deterministic
generation ID. Cutoff and density/flux retain that ID and fail if a replacement is
published between products.

Roadmap Step 6 requires no SWMF-specific spectrum implementation. After the current
coupler state is frozen and published, `Mode3DForwardSWMF.cpp` calls
`Earth::Mode3D::RunDensityAndFlux()`. That path uses the same
`BoundaryProducts.h` kernel as standalone Mode3D and gridless runs, including explicit
MeV/particle or MeV/nucleon coordinates, log-intensity time-table selection, boundary
uncertainty, normalized/raw PAD and spatial modes, differential spectra, density,
omnidirectional and planar flux, channels, detector rates, and unresolved bounds.

This is an instantaneous quasi-static fold at each published SWMF snapshot. It does
not yet transport a distribution across changing coupled states or enable electric
acceleration. Even though the compact snapshot includes derived `E=-v×B`, production
Step 6 deliberately selects the static-magnetic `J_local=A*J_boundary` mapping; the
general `j/p^2` mapping remains a tested kernel for the later electromagnetic mover.
Full cadence synchronization and time-history products remain Roadmap Step 11 work.

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

Although the coupled snapshot contains the derived electric field, the released Step 4
backtracer remains frozen magnetic-only. Electric or explicitly time-dependent
characteristics and `PhysicalBackwardTime` requests fail fast; the code does not treat
the presence of an E array as permission to use the static antiparticle shortcut.

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
```
