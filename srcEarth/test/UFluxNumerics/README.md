# Shared flux-numerics unit tests

This directory tests the backend-independent kernels in `util/FluxNumerics.h` without
requiring MPI, AMPS, SPICE, Geopack, or an SWMF installation.

Run from any directory:

```bash
./test/UFluxNumerics/run_test.sh
```

From the AMPS repository root it is also executed automatically as the first entry in
the normal validation list:

```bash
python3 srcEarth/test/test_runner.py -j 6 srcEarth/test/list
```

Its independent `last pass:` line is intentionally left empty until a successful run
is committed with the test runner's last-pass update workflow.

The runner compiles with C++11 plus `-Wall -Wextra -Werror -pedantic`, then exercises:

- **U-F03** — relativistic energy/rigidity round trips and the speed bound;
- **U-F04** — LINEAR, LOG, and log-rigidity grids, including the Mode3D linear-grid
  regression (the second and middle nodes must be linear, not proportional to `a*a`);
- **U-F05** — trapezoidal total-flux and clipped-channel quadrature, and a density
  unit/finite-value check;
- **U-F06** — unit direction vectors, angular weights summing to `4*pi`, and repeatable
  deterministic subsampling;
- **U-F07** — the local phase-space-density to differential-intensity identity and
  calibrated, nonuniform angular integration used by observational product adapters;
- **U-F12** — structured termination accounting, conservative unresolved bounds,
  retry counts, tolerance decisions, and the all-unresolved `NaN` nominal result.

U-F07 deliberately does not integrate over energy: an energy-centroid array is not
an energy-bin-edge definition. Instrument adapters must supply documented masses and
calibrated angular weights, and invalid values fail instead of being replaced by an
isotropic approximation.

These tests validate numerical infrastructure, not magnetic-field physics. Full
gridless/Mode3D comparison tests and observational validation still require the normal
AMPS runtime and data dependencies.
