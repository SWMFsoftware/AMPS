# Native spherical dust-injection qualification

This application is selected by `Config.pl -application=test/model-dust-injection`.
It is a one-rank T7 fixture for the actual
`ElectricallyChargedDust::DustInjection__Sphere` path, not a standalone dust
science model.

The fixture creates a real 3-D AMR mesh and non-origin spherical boundary,
allocates native DUST records with individual statistical weights, invokes the
production injection and boundary mover, and verifies:

- the final particle belongs to a locally allocated AMR owner block;
- radius, per-grain mass, compiled dust species, outward velocity, and
  individual weight correction survive insertion;
- one represented grain closes two half-grain physical source intervals via
  the signed residual mass reservoir; and
- the sampled represented mass rate equals the inserted mass divided by the
  production time step; and
- the general gas NTC kernel leaves dust records unchanged even with a
  `1e40` contrast of individual corrections. The test supplies two valid
  cloned grains, invokes the production kernel and compares their complete
  record bytes. Dust must bypass gas prefix sums and gas collisions; its
  steep size/weight law is not an invalid gas sampling distribution.

The deterministic callback fixes radius and direction so the oracle tests
native ownership/accounting rather than relying on a stochastic tolerance.
The dependency-light T7 cases separately qualify the steep power-law radius
sampler and isotropic spherical direction law. Multi-rank source arbitration,
migration, restart, and MPI reduction remain T10 responsibilities.

The registered runner performs the mandatory clean native preflight. Direct
manual builds must use the same sequence from the repository root:

```sh
rm -rf -- build
./Config.pl -application=test/model-dust-injection -spice-path=nospice
env MAKEFLAGS=-j16 make amps
mpiexec -n 1 ./amps
```

The manufactured sphere has no external ephemeris dependency. The adapter
explicitly resets global SPICE selection, verifies effective compiler options,
and checks actual clean-build commands before MPI. Its
`dust-injection-native-configuration.json`/`.txt` receipts distinguish a
configuration ERROR from the native physical ledger result. This prevents a
previous application from silently adding a library/kernel requirement.
