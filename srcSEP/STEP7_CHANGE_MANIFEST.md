# Step 7 change manifest: canonical Parker solver

> Step 14 supersession: `mover.cpp` and the historical symbols listed below
> have now been deleted. See `MIGRATION_MANIFEST.md` for final replacements.

## Implemented

- Added `SpatialDiffusionProvider` and its SI/provenance result contract in
  `util/sep_transport_coefficients.h`.
- Added the dependency-light canonical Itô kernel in
  `util/sep_parker_core.*`, including convection, the positive
  `d(kappa)/ds` drift, Gaussian diffusion, and exact adiabatic momentum.
- Added `parker_mover.cpp`, whose production provider delegates to configured
  `GetDxx`, applies named segment-scale timestep limits, advances through the
  common metric adapter, samples mean free path consistently, and queues only
  the declared pitch-angle-averaged turbulence feedback.
- Updated the production mover registry, default callback, makefile, and legacy
  component-test call sites to use `ParticleMover_Parker` as the sole public
  Parker implementation.
- Added `PARK01`–`PARK07`, `PARK-SOURCE`, and `make test-parker-unit`.
  `PARK07` compares absorbing-boundary exit probability and mean first-passage
  time with the exact drift-free Brownian solution.

## Public boundary

The CLI remains frozen at `parker`; historical source functions may remain for
migration reference but are absent from the public registry and production
object selection. The production SDE consumes one coefficient-provider result
per substep and contains no hidden random-number source.

| Former symbol/location | Production replacement | Compatibility |
|---|---|---|
| `ParticleMover_Parker_Dxx`, `mover.cpp` | `ParticleMover_Parker`, `parker_mover.cpp` | Historical implementation deleted by Step 14. |
| `ParticleMover_ParkerEquation`, `mover.cpp` | `ParticleMover_Parker`, `parker_mover.cpp` | Historical implementation deleted by Step 14. |
| `ParticleMover_Parker_MeanFreePath`, `mover.cpp` | `parker` for pitch-angle-averaged transport, or `fte-mfp` for an event-driven MFP closure | No ambiguous CLI compatibility alias; rejected names require an explicit canonical choice. |

## Remaining native gate

The focused suite verifies manufactured solutions, ensemble moments, automatic
refinement order, and first-passage statistics in the exact pure production
core. A full AMPS checkout must still verify a real PIC
particle, segment-boundary flux, restart/campaign seed configuration,
serial/OpenMP/MPI distribution equivalence, and coupled background updates.
