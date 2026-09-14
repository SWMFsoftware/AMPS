# Step 9 change manifest: event-driven focused transport

> Step 14 supersession: the unlinked `fte_mover.cpp` migration copy has now
> been deleted; `focused_transport_mfp.cpp` is the only event-driven adapter.

## Implemented

- Added `MeanFreePathProvider` with SI lambda, momentum validity, provenance,
  and turbulence-state identity.
- Added `util/sep_focused_transport_mfp_core.*`: exact exponential event
  sampling (`nu=v/lambda`), infinite-lambda ballistic transport, focusing and
  exact cooling between events, counter-propagating Alfvén branch selection,
  isotropic wave-frame redistribution, and exact Lorentz transforms.
- Added the thin production adapter `focused_transport_mfp.cpp`. It composes
  named limits for the global end, segment scale, focus, cooling, shock
  crossing, and immutable snapshot validity.
- Removed `fte_mover.o` from `MAINLIBOBJ`; Step 14 later deleted
  `fte_mover.cpp`. The public `fte-mfp` registry entry continues to resolve
  to `ParticleMover_FocusedTransport_EventDriven`, now defined by the thin
  adapter.
- Replaced legacy speed and lambda floors/clamps with structured errors or the
  explicit `ballistic` policy. Invalid samples and ballistic substitutions are
  counted.
- Added `FTEM01`–`FTEM08`, `FTEM-SOURCE`, and `make test-fte-mfp-unit`.
  `FTEM08` compares the finite-time persistent-random-flight MSD and its
  diffusion limit with the independent velocity-correlation solution.

## Fixed-seed validation

The dependency-light suite uses deterministic keyed streams and checks the
exponential and Poisson distributions, infinite-lambda motion, wave-frame
invariance, manufactured focusing/cooling solutions, and refinement under
combined operators. Refinement cases report automatically calculated observed
order rather than a hard-coded ratio. Native PIC field-line attachment, MPI
decomposition, and long campaign conservation remain host-integration gates.
