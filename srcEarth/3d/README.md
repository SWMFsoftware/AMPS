# Mode3D cutoff, flux, and spectrum backend

`DensityMode3D.cpp` computes backward-access transmission, local differential spectra,
number density, total omnidirectional flux, and configured energy-channel fluxes using
the compact AMR magnetic/electric-field snapshot.

The same code path serves two configurations:

- standalone Mode3D, after `Mode3DPrepareMagneticFieldSnapshot()`; and
- SWMF coupling, after `PrepareGlobalSWMFCoupledMagneticFieldForCutoff()` assembles the
  current coupler snapshot.

The field evaluator is the only intended numerical difference from the gridless solver.
Energy/rigidity conversion, energy coordinates, angular sampling, channel clipping,
quadrature, transmission diagnostics, and unresolved bounds come from
`../util/FluxNumerics.h` in both backends.

## Step 4 common trajectory request/result

`TraceTrajectoryMesh(prm, request)` consumes the same `TrajectoryRequest` and returns
the same `TrajectoryResult` as the direct gridless backend. The request-specific mover
is passed through the full integration, numerical retry, and unresolved-extension
stack; it is not re-read from `gDefaultMover` inside a trajectory. Time, step, and path
budgets are applied to an isolated `AmpsParam` copy, leaving the run configuration and
compact field immutable.

Before the first step, Mode3D hashes the currently published
`GlobalMagneticField::CurrentSnapshotMetadata().snapshotId` and compares it with every
nonzero request fingerprint. A mismatch or missing required identity is fatal rather
than a forbidden trajectory. This applies equally to standalone compact fields and to
SWMF-published generations.

An allowed captured trajectory returns position and SI momentum at the same exact
outer-box chord fraction, velocity direction, pitch cosine, event time, exit rigidity,
and `valid=true`. The common `TrajectoryContract.h` interpolation helpers are shared
with gridless, while only the field evaluation differs. All non-allowed terminations
leave the exit state invalid.

The released solver remains frozen-B only. Compact snapshots may contain electric
field data, but setting `electricFieldEnabled`, `fieldTimeDependent`, or
`PhysicalBackwardTime` in a request fails before integration because no physical-
backward electromagnetic mover has been released. Adaptive access refinement is a
separate Step 5 concern.

## Step 3 immutable compact-field generation

`GlobalMagneticField` publishes compact B/E arrays and `FieldProvider.h` metadata as a
single generation. Assembly proceeds in this order:

1. invalidate the previous generation;
2. gather exactly one owner value for every used interior cell;
3. reject missing, duplicate, or non-finite B/E values;
4. validate the GSM/SI/source/epoch/interpolation metadata;
5. publish arrays, metadata, snapshot ID, and generation together;
6. run requested products against read-only arrays; and
7. verify that the snapshot ID is unchanged after each product.

An `IFieldSnapshot` view captures the publication generation. If a new field is
assembled or the arrays are cleared, an old view returns `STALE_EPOCH`; it never reads
new arrays under an old identity. Standalone metadata includes the existing Step 2
field and electric drivers plus mesh/domain controls. Step 3 does not add new empirical
field models or modify interpolation, movers, trace limits, or validation tolerances.

## Unresolved trajectories

`TraceTrajectoryMesh()` returns a structured termination code.  Physical inner losses
and validated trapping are resolved forbidden states.  Time, step, distance, invalid
field, invalid time step, and numerical failures remain unresolved.  For each energy,
Mode3D now retains sampled, resolved, allowed, retried, and per-termination counts.

Outputs use these values as follows:

- `T` uses resolved trajectories only and is `NaN` if none resolves;
- `T_lower` assumes all unresolved directions are forbidden;
- `T_upper` assumes all unresolved directions are maximally allowed;
- density, local spectrum, total flux, and channel flux repeat the same nominal/lower/
  upper convention.

`mode3d_termination_summary*.dat` contains the accounting when
`DS_SAVE_TERMINATION_SUMMARY=T`.  Set `DS_FAIL_ON_UNRESOLVED=T` in validation cases to
abort when any location/energy exceeds `DS_UNRESOLVED_TOL`; leave it false for survey
runs that intentionally continue with explicit uncertainty bounds.

## Linear-grid correction

The former Mode3D LINEAR branch evaluated `Emin + a*(Emax-Emin)*a`.  The shared builder
uses `Emin + a*(Emax-Emin)`, pins both endpoints exactly, and is covered by U-F04 in
`../test/UFluxNumerics`.

## Verification

Run the dependency-free shared numerical suite from `srcEarth`:

```bash
./test/UFluxNumerics/run_test.sh
./test/UFieldProvider/run_test.sh
./test/UTrajectoryCore/run_test.sh
```

A full Mode3D or SWMF build still requires the parent AMPS build tree and its configured
MPI, SPICE, and model dependencies.
