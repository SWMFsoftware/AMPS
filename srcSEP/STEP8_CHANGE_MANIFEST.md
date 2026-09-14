# Step 8 change manifest: coefficient-driven focused transport

> Step 14 supersession: the historical `mover.cpp` and
> `fte_mover_dmumu.cpp` implementations described below have now been deleted.

## Implemented

- Added `PitchAngleDiffusionProvider` with `Dmumu`, `dDmumu/dmu`, SI units,
  status, coefficient provenance, and turbulence-state identity.
- Added `util/sep_focused_transport_core.*` with symmetric deterministic /
  stochastic / deterministic splitting, midpoint deterministic pitch-angle
  evolution, exact cooling, Itô derivative drift and noise, midpoint streaming,
  and multi-crossing reflective pitch-angle boundaries.
- Added `focused_transport_dmumu.cpp` as the canonical `fte-dmumu` adapter. It
  uses the configured coefficient function through the Step 10 shared provider,
  provides a documented numerical
  derivative mode, applies four independent timestep constraints, and advances
  exclusively through the common field-line adapter.
- Corrected the Kolmogorov QLT normalization and full SI `Dmumu` expression in
  `QLT.cpp` and `diffusion_dmumu.cpp`.
- Added thread-local, post-step deterministic wave-coupling queues. A mover
  retains records until successful particle reattachment, preventing absorbed
  particles from leaving dangling deferred handles.
- Updated both standalone and coupled timestep drivers to flush contributions
  after particle workers finish, and removed the experimental direct-wave
  Dmumu mover from production objects. Its historical source file
  `fte_mover_dmumu.cpp` remains available for migration review but is not linked.
- Added `FTED01`–`FTED07`, `FTED-SOURCE`, and
  `make test-fte-dmumu-unit`.
- The controlled-validation enhancement moves these cases into the common
  `--test` registry and adds `FTED08`, which checks analytical Legendre `P2`
  decay for `Dmumu=D0(1-mu^2)`. The focused Make runner and linked CLI now call
  the same callbacks.

## Coupling contract

The coefficient and the numerical deposition record carry the same immutable
turbulence identity. The particle mover does not evolve wave arrays directly;
it uses the existing reduced feedback interface after deterministic queue
ordering. Boundary-exit wave deposition is intentionally omitted while the
legacy callback requires a live particle record.

| Former symbol/location | Production replacement | Compatibility |
|---|---|---|
| `ParticleMover_FTE`, `mover.cpp` | `ParticleMover_FocusedTransport_Dmumu`, `focused_transport_dmumu.cpp` | Historical implementation deleted by Step 14. |
| `ParticleMover_FTE_DmuMu`, `fte_mover_dmumu.cpp` | `ParticleMover_FocusedTransport_Dmumu`, `focused_transport_dmumu.cpp` | Historical source deleted by Step 14. |
| hard-coded `QLT::calculateDmuMu(v,mu,r)` mover path | configured `PitchAngleDiffusionProvider` | No fallback alias; provider choice is explicit in runtime configuration. |

## Remaining native and scientific gates

The source-only suite establishes kernel behavior, units, QLT reconstruction,
refinement, and source dispatch. It does not establish native PIC compilation,
MPI/restart identity, long-run particle-plus-wave conservation, distributional
equivalence across decompositions, or agreement with observational data.
