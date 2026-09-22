# Shared numerical utilities

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
