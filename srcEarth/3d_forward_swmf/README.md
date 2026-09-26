# SWMF-coupled field snapshots, trajectories, and directional access

## Roadmap Step 11: coupled flux/spectrum callback

The backward-product callback now completes flux and spectrum products from the same
immutable snapshot used by coupled cutoff/access. Once Step 9 has assembled and frozen
the current SWMF B/u arrays, the callback makes a value-owned parameter block at the
snapshot's absolute UTC and calls `Earth::Mode3D::RunDensityAndFlux()` directly. This
is the shared Step-6 backward-characteristic product integrator. The forward
`cDensity3D` and `cSphereFlux3D` samplers remain confined to historical forward mode
and are not used to validate or fill a Step-11 result.

### Time, ephemeris, and control synchronization

For POINTS and SHELLS, the configured geometry is evaluated at the frozen field UTC.
For TRAJECTORY, `BuildSynchronizedCoupledProductParam_()` retains only spacecraft
samples matching that UTC within 1 ms, rebuilds the flattened GSM point list, and
remaps LOCATION-qualified directional apertures. No match is fatal; the code never
forms an unintended field-epoch × trajectory-time Cartesian product.

Before any trajectory collective, every rank independently fingerprints and compares:

- species charge and mass;
- energy range, spacing, access scan, and direction budget;
- the raw boundary-spectrum definition, basis, units, mass number, and uncertainty;
- ordered integral-flux channels and detector-response support/geometric factors; and
- coordinate frame plus selected point/trajectory/shell state and ephemeris UTC.

The boundary spectrum is then selected once at the final snapshot UTC. Every numeric
file records that requested UTC, the active/interpolated table epoch, temporal status,
gap flag/fraction, and the complete product/response/observation fingerprints.

### Fail-closed completion

Coupled mode requires the termination summary even if an input would normally suppress
that optional diagnostic; this does not change trajectory outcomes or thresholds. The
Mode3D writer records files only after flush and close. Rank zero then requires a local
spectrum, density plus flux (or a combined shell density/flux file), termination
summary, exact termination-count closure, nonempty artifacts, matching snapshot and
spectrum epochs, and the configured unresolved tolerance. It writes
`swmf_flux_spectrum_manifest<SUFFIX>.json` only after those checks. The frozen field is
released and `swmf_product_status<SUFFIX>.json` changes to PASS only after the manifest
itself closes successfully.

Run the portable tests with:

```bash
./srcEarth/test/USWMFCoupledProducts/run_test.sh
```

Use `../examples/swmf_step11_flux_spectrum.in.template` for live runs and the Step-9
snapshot replay template for the offline reference. Exact multi-artifact comparison,
1x1/2x8/8x16 and scheduler matrices, response closure, and O3 requirements are in
`../test/USWMFCoupledProducts/README.md`. This path is Phase-1 quasi-static and
magnetic-only; enabling experimental E does not release a time-dependent product map.

## Roadmap Step 10: production cutoff/access callback

The cutoff/access hook is now a collective, content-identified transaction rather than
a rank-local callback counter. After `ReadyForBackwardProductCalculation()` confirms a
complete receive, `ShouldRunBackwardProductCalculation()` compares the authoritative PT
time, configured `FIELD_UPDATE_DT`, and last-completed scheduler state on every MPI
rank. All ranks either skip or run. An in-process clock rollback or divergent cadence
state is fatal; a skipped callback never enters a solver collective and never consumes
the cadence slot. The current epoch is committed only after every requested product
returns successfully. A newly started/restarted process evaluates its first complete
epoch again, which provides the required restart-parity result.

Final cutoff/access filenames are constructed only after Step-9 B/u assembly creates
the content-derived snapshot ID:

```text
.swmf_t0000003600.000000000s_sidfield-v1-...
```

Callback number, MPI decomposition, and thread count are deliberately absent from this
scientific identity. The diagnostic mesh dump uses exactly the same suffix. Every
cutoff/access Tecplot file embeds `AUXDATA` for snapshot ID, absolute epoch, mesh
revision, content fingerprint, and boundary policy. Rank zero records files only after
`fclose` succeeds, verifies that each is nonempty, and writes
`swmf_cutoff_access_manifest*.json`. The manifest has `RESULT: PASS` only after the
complete artifact set exists. The older fail-closed `swmf_product_status*.json` remains
`FAILED` until all requested products finish. Before field assembly can determine a
snapshot ID, a clearly labelled `swmf_attempt_n...` status protects early fatal exits;
after assembly the final snapshot-named FAILED status is durably written before that
provisional marker is removed.

Step 10 activates the parsed outer-boundary policy:

- `BOUNDARY_TYPE BOX` preserves the historical six-face physical escape box.
- `BOUNDARY_TYPE SHUE` uses the Shue magnetopause and `DOMAIN_X_MIN` nightside tail
  cap. `SHUE_R0 AUTO` and `SHUE_ALPHA AUTO` use configured PDYN [nPa] and IMF Bz [nT].
  Numeric r0 is in Re and alpha is dimensionless.

For SHUE, a particle that reaches another computational face while still inside the
magnetopause terminates as `INVALID_FIELD`; it is not counted as allowed access. An
in-domain missing AMR leaf or incomplete compact interpolation row receives the same
unresolved termination. BOX follows its pre-Step-10 event code unchanged.

`POINTS`, `TRAJECTORY`, and `SHELLS` continue to call the common Mode3D movers,
adaptive/fixed access searches, termination accounting, and MPI schedulers. Progress
and filesystem writes remain root-only. Phase-1 output is instantaneous/quasi-static:
the frozen B snapshot does not represent a time-dependent characteristic.

Use `../examples/swmf_step10_cutoff_access.in.template` for the live configuration and
`../examples/standalone_step9_swmf_replay.in.template` for exact offline replay. The
numerical tests, live/replay comparator, linked 1x1/2x8/8x16 acceptance matrix, and
failure cases are specified in `../test/USWMFCoupledAccess/README.md`.

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
