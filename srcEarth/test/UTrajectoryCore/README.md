# Shared trajectory-core unit tests

This standalone suite validates the Step 4 trajectory contract without MPI, AMPS,
SPICE, Geopack, or SWMF. Run it from any directory:

```bash
./srcEarth/test/UTrajectoryCore/run_test.sh
```

It compiles the production `GridlessParticleMovers.cpp` with C++11 and strict warnings.
The checks are:

- **U-F10 — analytic uniform-field orbit:** BORIS and RK4 are compared with the closed-
  form relativistic helix in a uniform magnetic field over three timesteps. The test
  requires decreasing phase-space error, the expected convergence trend, and BORIS
  momentum-magnitude conservation. This is a reference-solution comparison.
- **U-F11 — backward-time/E-field semantics:** request validation rejects static-B
  conventions when E or explicit time dependence is enabled, rejects unreleased
  electromagnetic and reduced-orbit combinations, checks an exact uniform-E momentum
  forward/backward update, and checks the Liouville momentum factor.
- **U-F12 — deterministic termination/retry policy:** numerical retry and unresolved
  extension budgets are bounded and distinct; every termination reason is accumulated
  once and the category counts must close to the sample total.

The generated AMPS build normally supplies `constants.h`. This directory contains a
test-local file defining the exact SI speed of light so this unit can remain standalone.
