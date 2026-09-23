# UTrajectoryCore — strict Roadmap Step 4 tests

This dependency-free suite validates the shared trajectory request/result contract and
the production full-orbit mover dispatch without MPI, AMPS, SPICE, Geopack, or SWMF.
Run it from the repository root or from `srcEarth`:

```bash
./srcEarth/test/UTrajectoryCore/run_test.sh
# or, from srcEarth:
./test/UTrajectoryCore/run_test.sh
```

The runner first compiles the gridless and Mode3D public headers together using the
production C++17 standard and strict warnings. It then compiles the real
`gridless/GridlessParticleMovers.cpp` plus the dependency-free contract test with
C++11 and `-Wall -Wextra -Werror -pedantic`. It does not substitute a permissive test
mover.

## Gates and reference solutions

| ID | Gate | Independent reference or invariant |
|---|---|---|
| U-F10 | BORIS and RK4 mover convergence | Closed-form relativistic helix in uniform `B`; three step sizes, monotone error, order-ratio floors, absolute finest-step error ceilings, and Boris momentum conservation |
| U-F11 | Request/backward-time contract | Explicit negative-input matrix; exact `dp/dt=qE` forward/backward cancellation; Liouville `p²` intensity relation |
| U-F12 | Complete outer-boundary state | Analytic linear event-fraction reference for position, momentum, direction, time, rigidity, and pitch angle; invalid event/field paths must fail closed |
| U-F13 | Termination and bounded retry policy | Exact retry/extension counts, no distance-cap relaxation, exact scaled step budget, invalid-policy rejection, and closure over every termination enum |

The fixed snapshot-fingerprint assertion is a literal FNV-1a reference value. This
prevents a platform-specific `std::hash` or an accidental identity-algorithm change
from weakening provenance checks.

No tolerance, trajectory limit, mover, or expected value from an existing C/F test is
changed by this suite. A failure should be fixed in the implementation or explained as
an intentional physics-contract change; do not relax these gates to accommodate a
regression.
