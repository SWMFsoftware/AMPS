# Exosphere finite surface inventories

`Exosphere_FiniteInventory.cpp` implements the November 12 user-approved
continuous finite-reservoir release representation. An application opts in
through `UseFluxBalancedSurfaceInventory`; the default null callback and
prescribed reservoirs retain the existing source path.

For a synchronized species step with initial surface atoms N and constant
first-order coefficients k_j [1/s], the production kernel uses
`released_j=N*(-expm1(-sum(k_j)*dt))*k_j/sum(k_j)` and leaves the remainder
on the surface. PSD and thermal coefficients come from their existing
production rate/lifetime callbacks; no physical coefficient is fitted or
introduced. Temperature/illumination are those supplied by the production
operator split, held constant during this release interval. Incoming adsorbed
material enters the next synchronized reservoir update. Resolving coefficient
changes and within-step deposition still requires time-step convergence.

Each process/element's released mass is represented by full base-weight
packets plus a final fractional-weight packet. This is a continuous-inventory
Monte Carlo quadrature: atom totals are deterministic conditional on the
initial state; positions, velocities and departure times are sampled. It
does not reproduce microscopic atom-count shot noise. Individual particle
weights must be compiled ON; the dedicated Moon M2 builder enforces this.
The existing PSD/thermal generators retain their energy and angular laws.
Departure times follow the conditional exponential CDF, rather than uniform
times; stiff Ar release therefore occurs near the beginning of the interval.

Logical surface populations are replicated across MPI ranks. Rank zero draws
each candidate using the normal generator with an explicit source element,
then broadcasts it. The normal terrain postprocessor resolves final geometry;
exactly one AMR owner must allocate it. Actual accepted packet weights credit
gas; after all packets, every rank stores the analytic continuous remainder
directly. Releases are not also queued in `DesorptionFluxUP`, which would debit
twice. The independent ledger bounds the difference between the analytic
remainder and the sum of represented packet weights. Incoming adsorption stays
pending for the ordinary surface exchange. This avoids subtractive depletion
roundoff without inventing replacement atoms or clamping negative stocks.

The scoped `PIC::Mover::cMovingInjectionEpochScope` carries the sampled departure
offset [s] into the normal terrain postprocessor and moving-boundary tracer.
Moon evaluates the exact SPICE body-to-LSO transform at that departure epoch;
ownership is determined after final terrain positioning. The same epoch then
starts the shortened particle flight. This prevents combining an end-step
surface position with a different partial-flight surface epoch. The current
collective finite source path explicitly requires `OMP_NUM_THREADS=1`.

Duplicate/missing owners, unsupported weight/time
modes, invalid rates/populations and unrepresentable packet counts terminate
with explicit diagnostics. Packet counts are not silently truncated; process
memory/wall/disk guards remain the runner's operational responsibility.

Population-independent sputtering cannot be converted into an exponential
hazard by dividing by N. This unqualified coupling is explicitly rejected in
the finite path; other unsupported active depleting sources also fail. The
legacy prescribed-reservoir configuration remains available. Tests use the
same production kernel through U12 and real source/mover/surface callbacks in
the supporting I19 ledger. Formal I19 and observational qualification retain
all roadmap prerequisites and acceptance criteria.
