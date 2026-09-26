# Shared numerical utilities

## `SWMFCoupledProductsContract.h` — Roadmap Step 11

This dependency-free header defines the scientific identity and completion contract
for live SWMF density, local-spectrum, integral-flux/channel, detector-rate, and
termination products. Production and tests use the same functions; the suite does not
mock the rules in a separate implementation.

`ProductControl` contains every setting that changes the Step-6 fold: output mode,
species, energy/access sampling, boundary mode, raw spectrum key/value definition,
energy basis and units, spectrum uncertainty, ordered energy channels, detector
support/geometric factors, coordinate frame, ephemeris/point state, and shell geometry.
Canonical fingerprints sort only the map-shaped spectrum keys; channel, response, and
observation order is retained because it changes output schema and row meaning.
Dedicated channel, response, and observation fingerprints make a mismatch diagnosable
without weakening the complete control identity.

`ProductRunSummary` carries the actual spectrum evaluation and active-table epochs,
temporal interpolation status, location/energy/direction sizes, sampled/retried/
resolved/allowed counts, the complete termination vector, maximum unresolved fraction,
the configured unresolved tolerance, and the close-verified artifact inventory.
`ValidateRunSummary()` requires counts to close and, for a PASS transaction, requires
spectrum, density, flux, and termination roles plus the existing unresolved gate.
`BuildProductsManifestJson()` additionally requires complete Step-9 provenance, an
exact time/snapshot suffix, valid BOX/SHUE policy, and equality between field UTC and
boundary-spectrum evaluation UTC. It emits the explicitly limited
`INSTANTANEOUS_QUASI_STATIC`/`STATIC_MAGNETIC` schema.

Run the analytic, fixed-contract, negative, and exact-replay tests with:

```bash
./test/USWMFCoupledProducts/run_test.sh
```

The companion `compare_flux_products.py` validates a complete manifest-enrolled set,
including units and response/observation identity, before numeric comparison. Its
default tolerance is exact; any nonzero campaign tolerance must be passed explicitly.
The linked F6/F7/F13/F17, I-F03/I-F06/I-F07/I-F10, and O3 procedures remain separate
release evidence and are documented in the test README.

## `SWMFCoupledAccessContract.h` — Roadmap Step 10

This dependency-free header defines the decisions that must be identical in the live
SWMF component, standalone replay, and numerical tests:

- `CadenceGate` schedules distinct authoritative PT epochs, skips duplicates, rejects
  an in-process time rollback, and advances only after a successful product commit.
- `BuildProductSuffix()` binds nine-decimal simulation time and the complete Step-9
  snapshot ID into a restart/layout-independent filename.
- `BoundaryPolicy` separates a physical BOX or Shue+tail-cap escape from unavailable
  computational data. The Shue AUTO coefficients use explicit PDYN and IMF Bz inputs;
  missing AMR coverage inside the physical surface is never allowed access.
- `BuildAccessManifestJson()` emits the self-identifying, quasi-static Step-10 artifact
  transaction with an explicit `RESULT`. A PASS requires complete field provenance,
  finite time, a matching time/snapshot suffix, a valid domain/policy, and at least one
  distinct nonempty closed artifact.

The production Mode3D tracer calls these functions directly; the test harness does not
copy their formulas. Run:

```bash
./test/USWMFCoupledAccess/run_test.sh
```

The companion comparator requires identical snapshot/content/mesh provenance before
it compares numeric cutoff/access rows. Its default numeric tolerance is exact; a
campaign tolerance must be supplied explicitly and does not alter any C/F gate.

## `SWMFSnapshotContract.h` — Roadmap Step 9

`SWMFSnapshotContract.h` is the dependency-free interchange and lifecycle contract
shared by live SWMF publication, coupled export, standalone Mode3D replay, and the
Step-9 unit suite.  It defines the strict
`sep-in-geospace/swmf-field-snapshot/v1` CSV schema in GSM/SI:
position in metres, B in tesla, plasma bulk velocity in metre/second, and optional
experimental ideal-MHD E in volt/metre.  The default `MAGNETIC_ONLY` mode retains u for
provenance but makes E unavailable; `EXPERIMENTAL_IDEAL_MHD` must be selected
explicitly and requires every stored E vector to equal `-u x B` within the declared
numerical contract.

Identity is intentionally content-derived.  Canonical sorting by `(block,i,j,k)` makes
the mesh revision and state fingerprint independent of MPI decomposition or record
order.  The mesh revision covers topology and exact SI cell centres.  The content
fingerprint additionally covers epoch, authoritative simulation time, domain, units,
mode, B, u, and the derived-E convention; the final snapshot ID binds both. Readers
recompute all three and
reject metadata-only edits, stale geometry, missing/duplicate cells, non-finite values,
wrong units/frame/mode, negative time, and truncated files.

`Compare()` supplies the live-versus-replay numerical comparison with separate
position, B, u, and E tolerances and reports maxima instead of a boolean smoke result.
`PublicationQueue` models the fail-closed lifecycle: unavailable/stale/corrupt receives
cannot reuse the previous state, a valid receive behind a frozen batch is queued, and
release promotes only that validated next state.  `BuildProductStatusJson()` records a
machine-readable PASS/FAILED result for every attempted product batch.

The production parser includes `pic.h` directly before evaluating coupler-mode macros;
this is required so standalone and SWMF translation units cannot obtain different
preprocessor meanings through accidental include order.  Its Step-9 controls are
`SWMF_SNAPSHOT_FILE`, `SWMF_SNAPSHOT_EXPORT`, `SWMF_SNAPSHOT_EXPORT_PREFIX`, and the
strict `SWMF_DERIVED_ELECTRIC_FIELD OFF|EXPERIMENTAL` selector.

Run the fixed numerical identities, analytic `-u x B`, serialization, corruption,
tolerance, lifecycle, and production-wiring gates with:

```bash
./test/USWMFSnapshot/run_test.sh
```

## `StandaloneProductContract.h` — Roadmap Step 7

`StandaloneProductContract.h` is the dependency-free contract used by both standalone
dispatchers before either launches trajectories. It canonicalizes field aliases,
strictly parses single or combined product targets, validates POINTS/TRAJECTORY/SHELLS,
defines the exact external-driver columns and native units for each released field
model, binds all inputs and outputs to one epoch, and builds the machine-readable run
manifest. Unknown model, product, or domain tokens are rejected rather than accepted by
substring matching.

`RunPlan::Validate()` is intentionally separate from the numerical and scientific
acceptance gates. It proves startup provenance: immutable snapshot identity, driver
schema and unit validation, inclusive driver coverage of the requested epoch, required
Geopack initialization, and field validity. `RequireSameSnapshot()` from
`FieldProvider.h` is then called after every requested product. Neither routine alters
cutoff tolerances, unresolved classification, mover budgets, quadrature, or spectra.

`TsDriverTable` carries the provenance consumed by this plan. JSON metadata propagates
`UNIT`/`UNITS` to scalar and vector elements. Simple AMPS-wizard headers may write
bracketed units such as `By[nT]` and `Pdyn[nPa]`; old bracket-free files are accepted
only as the fixed documented wizard schema. Required row values must be complete finite
numbers, rows must have strictly increasing epochs, and production callers test
inclusive table coverage before interpolation. A wrong declared unit, missing or
malformed model driver, duplicate timestamp, or uncovered epoch fails.

Run the strict positive/negative contract and source-wiring references with:

```bash
./test/UStandaloneProducts/run_test.sh
```

## `BoundaryProducts.h` — Roadmap Step 6

`BoundaryProducts.h` is a header-only C++11 production kernel between trajectory
access and science products. It has no AMPS, MPI, SPICE, field-model, or mesh types, so
gridless, standalone Mode3D, and SWMF-coupled Mode3D cannot acquire different spectrum
normalizations or quadrature rules.

The unit contract is explicit:

| Quantity | Unit/convention |
|---|---|
| public energy coordinate | MeV/particle or MeV/nucleon, declared by `SpectrumUnits` |
| spectrum callback energy | coordinate energy converted to joule-equivalent units |
| callback intensity | m^-2 s^-1 sr^-1 per coordinate joule |
| saved differential intensity | m^-2 s^-1 sr^-1 per declared coordinate MeV |
| omnidirectional/planar integral flux | m^-2 s^-1 |
| number density | m^-3 |
| detector geometric factor | m2 sr |
| detector folded rate | s^-1 |

`BuildEnergyCoordinateGridMeV()` performs rigidity scans in total particle energy and
returns nodes in the declared coordinate. `IntegrateDensityWithUnits()` integrates over
that coordinate but evaluates relativistic speed from total particle energy. This is
the required Jacobian for MeV/nucleon spectra; treating the coordinate as particle MeV
would shift ion access and density by factors involving mass number.

`SelectTemporalSpectrum()` requires a rectangular, strictly increasing time table with
positive finite intensities and interpolates their logarithms. It reports exact,
interpolated, gap-interpolated, gap-held, clamped, or zero-outside status. Gap and
out-of-range choices are explicit; the caller can choose a deterministic earlier row
on an exact HOLD_NEAREST tie or fail.

`PadWeight()` and `SpatialWeight()` support `Raw` and `UnitMean`; their analytic
normalizers are also used when constructing conservative unresolved upper bounds.
`DirectionalAccessBounds()` maps allowed/forbidden/unresolved states to exact or bounded
access without inventing a nominal value for unresolved trajectories.

`EvaluateIsotropicProducts()` emits one `ProductSet`: density, omnidirectional and
isotropic-equivalent one-way planar flux, clipped channels, detector rates, and every
differential sample needed to reproduce the integrals. `FoldDirectionalDifferential()`
performs an explicit solid-angle and projected-area fold when directional samples are
available. Static production uses `J_local=A*J_boundary`; the separately tested general
phase-space branch multiplies by `(p_local/p_boundary)^2` for a future validated
electromagnetic mover.

Run its strict references with:

```bash
./test/UBoundaryProducts/run_test.sh
```

## `TrajectoryContract.h` — Roadmap Step 4

`TrajectoryContract.h` owns the common gridless/Mode3D/SWMF characteristic interface.
It is header-only C++11 and intentionally independent of AMPS, MPI, field-model
libraries, and mesh classes.

### Request and result

| Record | Required information |
|---|---|
| `Trajectory::Request` | SI launch position and unit direction, rigidity, physical charge/mass, per-request mover, explicit backward convention, time/step/path budgets, exit capture, field-use flags, reduced-orbit validity, snapshot fingerprint |
| `Trajectory::Result` | terminal reason, total time/path/steps, numerical retry and unresolved-extension provenance, mirror/bounce/drift evidence, mover and backward convention, verified snapshot fingerprint |
| `Trajectory::ExitState` | outer-event position, SI momentum, velocity direction, pitch cosine, event time, event rigidity, and a validity bit |

The request validator rejects non-finite state, non-unit direction, invalid species or
rigidity, missing budgets/identity, unknown enum values, undeclared reduced-orbit use,
and incompatible electromagnetic/backward-time combinations. Production adapters add
two backend checks: requested species must match the parsed run species, and every
nonzero snapshot fingerprint must equal the active immutable field generation.

`PopulateExitKinematics()` interpolates boundary position and momentum at the same
first-event fraction. `CompleteExitPitchAngle()` then uses a finite nonzero boundary
field to finish the record. Both direct and mesh backends call these helpers, so flux
and spectrum code cannot receive backend-dependent asymptotic state definitions.

### Retry and unresolved extension

`RetryPolicy` keeps two mechanisms distinct:

1. only an invalid timestep, invalid field, or numerical failure can receive the
   bounded smaller-step recovery retry;
2. only `TIME_LIMIT` or `STEP_LIMIT` can receive an explicitly configured larger-time
   convergence pass, restarted from the original seed.

`DISTANCE_LIMIT` is never expanded by the time-convergence policy. Primary and final
terminations, budgets, and pass counts remain in the result. Gridless and Mode3D call
the same `ShouldRetryNumerical()`, `ShouldExtendUnresolved()`,
`ExtensionTimeBudget()`, and `ScaledStepBudget()` functions.

Run the strict reference/invariant suite with:

```bash
./test/UTrajectoryCore/run_test.sh
```

See `test/UTrajectoryCore/README.md` for the U-F10 through U-F13 gates.

## `AdaptiveDirectAccess.h` and `DirectionalAccess.h` — Roadmap Step 5

`AdaptiveDirectAccess.h` constructs one deterministic geometric-rigidity candidate
tree per direction. Every requested seed is evaluated. Configured guard levels probe
midpoints even when the endpoints agree, which can expose a hidden non-monotone access
or forbidden pocket. Visible state-change and resolved/unresolved brackets are then
refined until

```text
Delta R <= max(absoluteTolerance_GV, relativeTolerance * R)
```

or until a hard maximum depth/sample count is reached. The returned report distinguishes
successful convergence from depth or work-budget exhaustion and carries the maximum
ambiguous width, summed bracket-width support, and response-weighted unresolved
support. It never relabels an unresolved trajectory to make the target pass. The legacy
depth-only API remains as a zero-tolerance compatibility wrapper.

`DirectionalAccess.h` defines one backend-neutral saved `A(E,Omega)` sample and the
field-independent reconstruction of lower, effective, and upper cutoff plus penumbra.
An allowed sample is valid only when it contains a finite, internally consistent Step-4
exit position, momentum, velocity direction, pitch cosine, event time, and rigidity.
Physical-forbidden and unresolved samples must not carry a valid exit state. The exact
regular lon/lat cell-area helper closes to `4*pi`, including polar half cells.

The production gridless and Mode3D/SWMF writers expose identical columns and repeat the
per-direction convergence report on every realized sparse row. Consequently a consumer
can reconstruct the cutoff diagnostics and angular integral from the saved file alone.

Run the strict analytic, invariant, negative, deterministic, and C19-reader tests with:

```bash
./test/UDirectionalAccess/run_test.sh
```

See `test/UDirectionalAccess/README.md` for the U-F14 through U-F17 references and
acceptance conditions.

## `FieldProvider.h` — Roadmap Step 3

`FieldProvider.h` is a dependency-free C++11 contract that separates a trajectory
solver from the origin of its background field:

```text
IFieldProvider.CreateSnapshot(request) -> immutable IFieldSnapshot
IFieldSnapshot.Sample(query)           -> B, E, status, snapshot ID
```

Every valid snapshot states its source/model, frozen epoch, coordinate frame, SI units,
validity domain, interpolation method, B/E availability, and deterministic physical-
state identity. `MakeSnapshotId()` uses stable FNV-1a over canonical source, epoch, and
field-defining state. Request labels and output filenames are deliberately excluded.

`ValidateQuery()` keeps failure modes separate: `OUTSIDE_DOMAIN`, `STALE_EPOCH`,
`SOURCE_UNAVAILABLE`, `INTERPOLATION_FAILURE`, and `NONFINITE_VALUE` cannot silently
become magnetic shielding. `RequireSameSnapshot()` prevents cutoff and flux/spectrum
products from being combined after a field generation changes.

Step 3 does not change `IGridlessFieldEvaluator`, particle movers, trajectory limits,
or access classification. The direct evaluator implements the provider contract by
multiple inheritance internally, while movers continue using the original `GetB_T()`
call. Mode3D publishes metadata only after complete finite compact arrays are assembled.

Run the strict contract and closed-form reference suite with:

```bash
./test/UFieldProvider/run_test.sh
```

## `FluxNumerics.h`

`FluxNumerics.h` is a header-only C++11 module shared by the standalone gridless and
Mode3D/SWMF-backed flux solvers.  Header-only placement is deliberate: its unit suite
can compile without the AMPS link graph, while production translation units execute the
same inline implementation.

### Units and public operations

| Operation | Input | Output |
|---|---|---|
| `RigidityFromEnergyGV` | kinetic energy J, `abs(q)` C, mass kg | GV |
| `EnergyFromRigidityMeV` | GV, `abs(q)` C, mass kg | MeV |
| `RelativisticSpeed` | kinetic energy J, mass kg | m/s |
| `BuildEnergyGridMeV` | MeV bounds and scan controls | MeV nodes |
| `BuildEqualSolidAngleDirections` | `N_mu`, `N_phi` | unit vectors and sr weights |
| `IntegrateFlux` | MeV nodes, transmission, spectrum per J | m^-2 s^-1 |
| `IntegrateDensity` | MeV nodes, transmission, mass, spectrum per J | m^-3 |

`IntegrateFlux` clips a requested energy channel to the solver grid, inserts channel
endpoints, linearly interpolates transmission there, evaluates the spectrum at the
actual endpoints, and applies trapezoidal quadrature.  Density uses the same energy
coordinate and relativistic speed, so the only physics difference is the expected
`1/v(E)` factor.

### Structured access accounting

`AccessAccumulator` records every sampled direction, including its termination class
and retry status.  `ResolveAccess()` returns:

```text
T       = weighted allowed / resolved
T_lower = weighted allowed / sampled
T_upper = (weighted allowed + unresolved * maximum_allowed_weight) / sampled
```

The nominal value is `NaN` when `resolved==0`.  This is intentional: an exhausted
numerical budget is missing information, not evidence for zero access.  The lower and
upper curves remain finite and can be integrated to form conservative product bounds.

`TrajectoryTermination.h` defines which outcomes are physically resolved.  Do not add
time, step, distance, invalid-field, or mover failures to `IsResolvedTermination()` to
make a validation pass; adjust the trace controls or explicitly accept the resulting
bounds instead.

### Adding a numerical kernel

Keep this module independent of MPI, AMPS mesh classes, field models, and the global
spectrum object.  Pass spectra as callables and backend results as plain values.  Add a
case to `test/UFluxNumerics/test_flux_numerics.cpp`, compile with the provided runner,
and then route both production backends through the new function.
