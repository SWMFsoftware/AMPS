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

## Step 7 standalone product orchestration

Standalone Mode3D now accepts cutoff-only, density/flux/spectrum-only, and combined
targets. For every requested epoch it materializes one compact field generation,
validates a shared `StandaloneProducts::RunPlan`, writes
`standalone_mode3d_manifest[_snapshot].json`, executes the selected products, and
checks that the published snapshot ID is unchanged after each one. The manifest and all
product suffixes therefore refer to the same epoch and field generation.

The standalone field list is DIPOLE, IGRF, T96, T01, T05/TS05, TA15N, TA15B, and
TA16. Step 7 adds explicit Mode3D initialization and evaluation for T01 and both TA15
variants: T01 receives PDYN, DST, BY, BZ and G1..G3; TA15 receives PDYN, BY, BZ and
XIND after selecting its N or B coefficient set. Driver columns, units, monotonic time,
and inclusive epoch coverage are checked before interpolation and mesh construction.

All empirical wrappers remain in the existing non-SWMF compile-time branch. An SWMF
build obtains its field from the coupler and cannot accidentally link or dispatch a
standalone T01/TA15 model. Step 7 makes no change to coupled cadence or coupled physics.
Run `./test/UStandaloneProducts/run_test.sh` for the common startup and isolation
contract, then retain the existing Mode3D C/F tests for linked numerical validation.

## Step 6 common boundary and product kernel

`DensityMode3D.cpp` now sends the complete access curve to
`../util/BoundaryProducts.h`, the same production integrator used by gridless mode.
The result contains number density; total, configured-channel, and
isotropic-equivalent one-way planar flux; optional detector-response rates; and the
differential boundary/local/omnidirectional/planar spectra with nominal/lower/upper
bounds. Boundary uncertainty and unresolved access are combined monotonically without
changing the underlying trajectory classification.

Energy is interpreted in the spectrum's declared coordinate. For a per-nucleon ion
spectrum, the grid and channel bounds remain MeV/nucleon, whereas rigidity and speed
receive total particle kinetic energy. Rigidity-scan grids are constructed uniformly
in physical log rigidity and converted back to the declared output coordinate. The
legacy per-particle behavior is unchanged by default.

Point spectrum and flux outputs append Step-6 fields after their historical prefixes.
Shell density/flux output appends planar and detector products, and each shell also
writes `mode3d_shell_*km_spectrum.dat`, with one structured zone per energy. These rows
are sufficient to reproduce every reported trapezoidal integral. All files include
`AUXDATA` for mapping, energy basis, mass number, differential-intensity unit,
uncertainty, temporal selection, and planar convention.

Mode3D `DIRECT_ACCESS` likewise retains the exact 45-column legacy/Step-5 prefix and
appends nine boundary-factor/intensity fields. Static magnetic characteristics use
`J_local=A*J_boundary`. The shared phase-space mapping is tested but is not selected by
production until the electric/time-dependent mover is validated.

Input keys and detector blocks are documented in `../README.md` and
`../boundary/README.md`. Run `./test/UBoundaryProducts/run_test.sh` for the focused
closed-form references; existing F/C gates are unchanged.

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
separate opt-in product and does not change scalar cutoff behavior.

## Step 5 directional-access product

With `CUTOFF_SEARCH_ALGORITHM DIRECT_ACCESS`, Mode3D writes
`cutoff_3d_dir_access_loc_######.dat`. The product now contains the complete
three-state `A(E,Omega)` contract rather than only a binary access flag:

- exact regular-cell `direction_weight_sr` and `weighted_access_sr`;
- stable termination, retry, extension, and trap diagnostics;
- the full outer-boundary position, SI momentum, unit velocity, pitch cosine, event
  time, and exit rigidity for every allowed characteristic; and
- adaptive refinement/error/support metrics plus explicit target/sample-limit status.

The angular grid continues to use `DIRMAP_LON_RES`/`DIRMAP_LAT_RES` and is independent
of the rigidity refinement. Adaptive mode evaluates every `CUTOFF_RIGIDITY_LIST_GV`
seed, performs the configured guard probes, and refines all visible state changes
without assuming monotonic access. Its absolute/relative tolerances and sample cap are
read from the common Step-5 controls documented in `../README.md`.

Standalone Mode3D and coupled SWMF calls execute this same code after publishing their
respective immutable Step-3 snapshot. The coupled path therefore changes only the
field snapshot; it does not have a separate access schema or convergence algorithm.
Dense mode remains available for convergence/reference runs. Existing scalar cutoff
searches, trace policies, C/F thresholds, and reference tables are unchanged.

Run the dependency-free shared Step-5 tests with:

```bash
./test/UDirectionalAccess/run_test.sh
```

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
./test/UDirectionalAccess/run_test.sh
```

A full Mode3D or SWMF build still requires the parent AMPS build tree and its configured
MPI, SPICE, and model dependencies.
