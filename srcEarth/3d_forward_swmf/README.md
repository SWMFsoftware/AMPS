# SWMF-coupled field snapshots for cutoff, flux, and spectra

The coupled backward-product path implements roadmap Step 3 through the same
`util/FieldProvider.h` contract used by standalone fields. It does not call or link the
standalone DIPOLE/IGRF/Tsyganenko providers for production SWMF calculations.

## Source and freeze boundary

`PrepareGlobalSWMFCoupledMagneticFieldForCutoff()` is called only after every MPI rank
reports a mesh, valid magnetic-field and bulk-velocity offsets, and at least one
completed SWMF receive. It then:

1. visits authoritative owner-rank interior cells;
2. copies `PIC::CPLR::SWMF::MagneticFieldOffset` values;
3. reads `BulkVelocityOffset` and derives `E=-v×B`;
4. rejects missing, duplicate, or non-finite owner values;
5. MPI-replicates only compact B, E, and presence arrays; and
6. publishes those arrays with one immutable metadata record and snapshot ID.

No trajectory samples a live owner buffer after this boundary. Cutoff/access and
density/flux/spectrum run consecutively on the compact arrays and each checks that the
snapshot ID is unchanged. The output suffix and diagnostic mesh dump use the same SWMF
callback time.

## Metadata contract

| Item | Coupled value |
|---|---|
| Source | `PIC::CPLR:SWMF` |
| Model | `SWMF` |
| Epoch | configured reference epoch plus authoritative PT simulation-time offset |
| Frame | GSM |
| Units | position m, B T, E V/m |
| Interpolation | decomposition-independent cell-centred linear row stencil |
| Domain | configured Mode3D box plus the used AMR-tree validity check |
| Time semantics | immutable/quasi-static for one trajectory batch |

This Phase-1 implementation supports snapshot access calculations. It does not claim a
fully time-dependent Lorentz characteristic, trapped-population evolution, or local
acceleration/loss physics.

## Steps 4 and 5 behavior in coupled runs

After the compact SWMF snapshot is frozen, coupled Mode3D invokes the same common
trajectory and `DIRECT_ACCESS` code as standalone Mode3D. Thus termination reasons,
retry/extension limits, outer-boundary exit state, angular weights, adaptive rigidity
criteria, and saved `A(E,Omega)` schema are identical. The snapshot fingerprint carried
by the trajectory result provides the cutoff/flux synchronization provenance.

Although the compact coupled snapshot contains derived `E=-v×B`, the released Step-4
characteristic intentionally remains magnetic-only and quasi-static. Turning on E or
explicit time dependence requires `PhysicalBackwardTime` and an implemented
electromagnetic mover; validation fails before integration until that later roadmap
stage is delivered. This prevents a static-field reversal from being silently applied
to time-dependent SWMF physics.

## Verification

The dependency-free contract test is:

```bash
./srcEarth/test/UFieldProvider/run_test.sh
./srcEarth/test/UTrajectoryCore/run_test.sh
./srcEarth/test/UDirectionalAccess/run_test.sh
```

A configured SWMF/AMPS validation build must additionally compare selected owner-cell
coupler values with compact cell access and interpolated samples, verify explicit
pre-receive/out-of-domain failures, and run a combined cutoff plus density/flux case to
confirm one snapshot ID is retained through both products.
