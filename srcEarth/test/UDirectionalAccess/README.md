# Adaptive directional-access unit tests

Run from any directory:

```bash
./srcEarth/test/UDirectionalAccess/run_test.sh
```

The standalone **U-F13** suite validates Step 5 without AMPS, MPI, or a field model:

- an adaptive access curve is compared with an analytic 4 GV step and must meet its
  requested bracket-width tolerance;
- a narrow non-monotone forbidden island with allowed seed endpoints must be found by
  the guard probe and refined on both sides;
- a five-transition analytic multi-band/oscillatory penumbra must retain every island;
- response-weighted unresolved support must be bounded, and exhausting `MAX_SAMPLES`
  must produce an explicit non-converged report;
- lower, effective, and upper cutoff plus penumbra width are reconstructed using only
  saved `A(E,Omega)` samples, including their required allowed boundary exit states.

These are reference-function and invariant tests. Full dipole/SWMF trajectory and
observation comparisons remain in the C/F validation suites and require a configured
AMPS executable.
