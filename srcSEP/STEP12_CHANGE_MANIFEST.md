# Step 12 change manifest: deterministic parallel reductions

## Implemented

- Added `util/sep_reproducible_reduction.*` with worker-local contribution
  buffers and physical keys containing field line, segment, branch, particle,
  timestep, and purpose.
- Added a complete-key canonical sort and fixed long-double accumulation order
  for wave energy, streaming, and resonant-event counts.
- Added a gather-then-canonical-reduce MPI policy helper whose results are
  independent of rank partitioning.
- Kept MPI rank and worker identity out of stochastic keys and exposed
  campaign/particle/step/purpose-separated random streams.
- Added atomic integer warning/limiter counters with explicit reset semantics.
- Added portable evidence hashes that avoid host padding and endian dependence.
- Added PAR01–PAR05, source gates, and `make test-reproducibility-unit`.

Canonical particle-to-wave reductions and keyed draws are bitwise identical
for a fixed physical contribution set and seed across the tested worker and
synthetic MPI layouts. Native collectives that do not promise a fixed
floating-point tree require a documented numerical-tolerance policy and
process-count evidence.
