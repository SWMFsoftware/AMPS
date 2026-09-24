# UDirectionalAccess — strict Roadmap Step 5 tests

Run from the AMPS repository root:

```bash
./srcEarth/test/UDirectionalAccess/run_test.sh
```

The suite compiles the production `AdaptiveDirectAccess.h` and
`DirectionalAccess.h` headers without AMPS, MPI, SPICE, or a field-model library. It
contains two groups of tests:

- **U-F14 — adaptive directional access:** compares the sampler with exact analytic
  Heaviside, hidden-pocket, and five-transition access functions; checks absolute and
  relative error targets; verifies exact `4*pi` angular closure; exercises response-
  weighted unresolved support; and requires depth/sample exhaustion to remain an
  explicit non-converged result.
- **U-F15 — saved-product reconstruction:** reconstructs lower, effective, and upper
  cutoff plus penumbra width from saved `A(E,Omega)` rows, compares unresolved bounds
  with a hand-derived reference, and rejects missing/inconsistent exit states,
  termination mismatches, and non-increasing rigidity coordinates.
- **U-F16 — observation-reader schema:** feeds valid and deliberately corrupted Step-5
  Tecplot rows through the production C19 reader. It rejects inconsistent angular
  weights, incomplete or nonphysical exit states, partial schemas, and convergence
  metadata that changes within one direction; it also requires both production writers
  to expose the same schema and shared adaptive-control markers. Archived pre-Step-5
  seven-column files remain explicitly readable.
- **U-F17 — CLI contract:** compiles the production `cutoff_cli.cpp` and verifies every
  new convergence option, no-override default, help entry, and negative/non-finite
  rejection. Only the AMPS fatal-exit hook is replaced by an exception for the test.
- **U-F18 — C8 schema and complete-cube guard:** runs the production C8 self-test. It
  proves that the legacy 24-column, corrected append-only 45-column, and historical
  named 45-column layouts recover identical core samples and reductions; malformed
  widths/names/blanks fail explicitly; full polar coverage and C8-G11 are enforced;
  and a partial cube cannot reach the FOV/comparison reductions that consume numeric
  fractions. The physical and sign-mirrored East–West references remain discriminating.

These tests include independent reference solutions; they are not format-only or
self-consistency-only checks. Full mover/field comparisons remain in UTrajectoryCore,
F1–F6, and C1–C19. Observation comparisons remain in C8/C9/C10/C19 and are not replaced
by this dependency-free suite.
