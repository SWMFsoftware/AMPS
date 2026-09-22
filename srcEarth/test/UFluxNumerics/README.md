# Shared flux-numerics unit tests

This directory tests the backend-independent kernels in `util/FluxNumerics.h` without
requiring MPI, AMPS, SPICE, Geopack, or an SWMF installation.

Run from any directory:

```bash
./test/UFluxNumerics/run_test.sh
```

The runner compiles with C++11 plus `-Wall -Wextra -Werror -pedantic`, then exercises:

- **U-F03** — relativistic energy/rigidity round trips and the speed bound;
- **U-F04** — LINEAR, LOG, and log-rigidity grids, including the Mode3D linear-grid
  regression (the second and middle nodes must be linear, not proportional to `a*a`);
- **U-F05** — trapezoidal total-flux and clipped-channel quadrature, and a density
  unit/finite-value check;
- **U-F06** — unit direction vectors, angular weights summing to `4*pi`, and repeatable
  deterministic subsampling;
- **U-F12** — structured termination accounting, conservative unresolved bounds,
  retry counts, tolerance decisions, and the all-unresolved `NaN` nominal result.

These tests validate numerical infrastructure, not magnetic-field physics. Full
gridless/Mode3D comparison tests and observational validation still require the normal
AMPS runtime and data dependencies.
