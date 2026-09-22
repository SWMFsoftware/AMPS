# Shared numerical utilities

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
