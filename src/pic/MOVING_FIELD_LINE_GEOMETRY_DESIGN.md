# Callback-driven material field-line geometry: implementation design

Status: the generic container transaction, field tracer, ownership-map refresh
and native segment-exposure sampling are implemented. The application-specific
analytical state callback, moving-geometry checkpoint continuation and real MPI/
coupled qualification remain incomplete. See README.md for the actual API,
registered test commands and precise limits; remaining requirements below are
not claims that those capabilities have already been qualified.
This extends the existing `PIC::FieldLine::cFieldLine` representation rather
than introducing a second field-line or particle model.

## Agreed physical contract

* Existing vertices, including the first and last, move according to a
  user-supplied motion function. The field-line container owns connectivity,
  not a solar-wind model or a velocity integrator.
* Particle motion relative to the line excludes bulk solar-wind displacement.
  Motion of the vertices supplies that displacement once. Physical U and its
  gradients remain available to the existing momentum/scattering operators.
* If the moved first vertex is beyond the minimum heliocentric radius, extend
  toward the beginning along B by the current first segment's length. Repeat
  until the first vertex is at or inside that radius. U is not the extension
  direction: U and B need not be parallel, particularly on a Parker spiral.
* Remove trailing segments only when they are completely outside the maximum
  radius. Delete their attached particles once and record physical losses.
* New segments start empty. Existing configured source rules alone create
  particles; geometric insertion never copies or manufactures a population.
* Sampling length belongs to each segment. The AMPS physical timestep and the
  global output schedule remain unchanged.

## Interfaces and responsibilities

The public entry points are on `cFieldLine`. The implemented API also accepts
an explicit edge-state initializer when stored segment physics is allocated:

```cpp
Geometry::Status cFieldLine::PrepareGeometryStep(
    const Geometry::Context& context,
    const Geometry::Motion& moveVertex,
    const Geometry::InnerTrace& traceInner,
    const VertexStateInitializer& initializeVertex,
    GeometryPlan* plan,
    const SegmentStateInitializer& initializeSegment);

Geometry::Status cFieldLine::CommitGeometryStep(GeometryPlan* plan);
```

`Geometry::Context` supplies departure and arrival epochs [s], expected geometry
generation, heliocentric center [m], minimum and maximum radii [m], a spatial
tolerance [m] and an insertion limit. Stable identity is registered separately
on the line. The owner binds an immutable background to its callbacks; the
generic context does not itself store a background snapshot. Epoch, radius,
generation and owner checks reject before topology mutation.

`moveVertex` receives the old vertex index, a const position and both epochs,
and returns its new Cartesian position [m] or an explicit failure. The owner
can derive the material label from `FirstMaterialLabel + index`.
Returning a position permits an exact flow map, a controlled numerical
integrator, or an externally supplied position. The container must not assume
Euler motion `X += U*dt`, select an integrator, query a different provider, or
infer a timestep from a mutable PIC counter. The callback cannot mutate live
vertices, particles, sampling data or shared provider state.

`traceInner` evaluates the owner's B field at the target geometry epoch and
traces toward the beginning by a specified arc length [m]. Choose the sign of
B from continuity with the oriented first segment; do not assume +B is the
increasing-index direction. The tracer owns curvature resolution and must
finish the requested length or return an error. A straight displacement along the endpoint
tangent is not a generally valid curved-field trace. Tracing subdivisions are
spatial numerical work; they never shorten an AMPS particle timestep.

`initializeVertex` supplies the actual current/previous B, U, density, thermal,
wave and other required datum state, with provider identity/epochs/validity.
It must distinguish an existing advected material record from a newly entering
record. New geometry cannot obtain physical state by copying the nearest
vertex or by zero-filling required primitives. Current and previous positions
may have to be supplied for time-dependent derivatives; a new record cannot
pretend that it existed in an earlier material population. A model lacking
required history or primitive coverage returns an explicit unsupported status.

All four concerns remain separate: container topology, vertex trajectories,
magnetic tracing, and physical state ownership. The same container can then be
used by analytical and imported data owners without deciding their physics.
SWMF-imported positions must not also receive a local wind displacement.

## Prepare, validate, then commit

1. Enter a host boundary after workers, particle migration, sampling and read
   guards have joined. Reject closed loops for this open-line boundary policy.
2. Capture the entire old line and material segment identities. Require valid
   particle S/segment ownership and empty temporary migration lists.
3. Evaluate `moveVertex` once for every old vertex, including both endpoints,
   into candidate storage. Validate finite positions, nondegenerate segments,
   ordering and provider epoch consistency. No old vertex is changed yet.
4. On candidate geometry, repeatedly trace an inner node along B using the
   current first segment's length as the requested arc length. New nodes and
   segments receive unique, deterministic material identities and empty
   particle/sample buffers. Validate the full candidate physical state.
5. Plan trailing removals using the full segment geometry, not its endpoint
   radii. Stage particle-loss and retired-sampling records before deletion.
6. Stage all surviving particle coordinates/Cartesian views and rebuild length,
   tangent, pointer tables and ownership plans. Validate every candidate before
   committing any part of the line. Allocate required capacity and detached
   records before publication; provide cleanup for failed preparation.
7. In MPI execution, all participating ranks agree on readiness and the same
   topology/generation before commit. A rank-local callback failure cannot
   leave another rank holding a different segment count.
8. Commit the prepared positions, topology, particle-coordinate mapping and
   associated state together; publish the new geometry generation. Deletion
   follows the staged loss ledger exactly once. A plan is single-use and must
   reject a changed expected generation.

The current implementation's per-vertex mutation, conditional call inside
`PIC::TimeStep()`, and last-vertex early exit must not remain an additional
update path for a line managed by this owner. Use an explicit owner/registration
mechanism, not function-pointer identity, to prevent double advection. Other
applications keep their existing path until explicitly migrated.

## Boundary geometry and failure cases

Let r(X)=|X-center|. For a straight trailing segment A--B, calculate

```
D = B-A
q = clamp(dot(center-A,D)/dot(D,D), 0, 1)
r_nearest = |A + q*D - center|
```

The segment is wholly outside the maximum sphere only if its nearest radius
is outside under the declared boundary tolerance. Both endpoints can lie
outside while the segment crosses the sphere. Retain a segment that intersects
the supported sphere; reject a zero-length segment before this calculation.

These are coverage bounds, not exact particle-absorption surfaces: the first
vertex may lie below the minimum radius, and a retained last segment may
straddle the maximum radius. If particles must be absorbed exactly at either
sphere, that is a separate mover/source boundary contract, not an incidental
effect of deleting whole segments.

Full-length inner extension can overshoot below the minimum radius. It must
remain inside the field provider's actual valid domain. Trace coverage, null B,
an unresolved orientation, an unreachable sphere, an extension that makes no
inward progress, or a step that passes through the sphere and ends outside
again must return a diagnostic instead of looping indefinitely. A numerical
resource limit is explicit and produces failure; it never silently clips an
extension or changes the requested length. Reaching an exact physical inner
surface instead would require a separately declared clipping policy.

If trimming removes the complete line, return a structured inactive/empty-line
result with the staged losses. The application must handle its sources,
observers and shock queries explicitly; it must not recreate a line or index
a missing segment. Repeated geometry updates must not silently grow storage
without limit: boundary placement, segment-resolution growth and resource
limits must be diagnosed. General splitting/coalescing is a separate remeshing
task and cannot inherit sampling histories by copying buffers.

## Material particles, identities and physical ledgers

After k prefix insertions, a surviving particle changes from `S` to `S+k`.
Its segment identity and fraction remain unchanged, and its Cartesian position
is recomputed from that segment's moved endpoints. Preserve species, momentum,
pitch angle, weight, identity, birth history and RNG event counters. Renumbering
is not a physical cooling, acceleration or scattering operation. Any physical
momentum response remains in the documented particle operators and must not be
applied a second time during geometry commit.

Delete particles only on removed segments. Record line/material-segment IDs,
species, macroparticle counts, represented particle counts and relevant
energy/momentum losses in their declared frames. Retire sampling products
before recycling records. A vacuum segment has zero particle loss but may
still have finite sampling exposure and stored wave energy.

Surviving segment-local data remain attached to their material identity.
Index-based MPI ownership, edge/vertex exchange sizes, turbulence arrays,
shock segment indices, source handles and observer connectivity must all be
rebuilt or invalidated. Extensive wave energy and intensive B/U/n are different
quantities; geometry deformation must use the provider's conservative flux,
mass, area and energy contracts. The container cannot invent magnetic or
thermodynamic evolution by rescaling everything with length.

There is also a concrete source-scheduler prerequisite in
`srcSEP/field_line.cpp`: `FixedNativeLineLayout()` currently hashes the complete
non-SWMF coordinate sequence and segment count as a restricted fixed-geometry
source namespace. Motion/insertion would correctly trip that guard. Replace
that restricted namespace with explicitly registered stable physical line IDs
and immutable source configuration; do not remove the guard or reset event
counters after a geometry change. A topology generation identifies changing
geometry and is distinct from physical source identity. The analytical shock's
launch radius/epoch also remain their original trajectory data: adding a new
first vertex must not relaunch the shock or reset its source history.

## Per-segment sampling exposure

Each segment owns collecting and completed sample counts (checked uint64), collecting and
completed exposure [s], and the corresponding window identity. Counters swap,
flush, checkpoint and migrate with their sample buffers. They are initialized
for every allocated/recycled segment, including zero-sized datum layouts.

Increment exposure once for every eligible sampling call in which the segment
exists, including calls with zero particles. Never increment per particle,
species loop, mover substep or MPI contribution. For fixed dt and full-step
sampling, exposure is N*dt; segment age is not a reason to change dt. A segment
created between steps starts with zero exposure and first participates in the
next scheduled sample. Zero exposure means no measurement, distinguishable in
output validity metadata from a measured zero density.

For a timed average, divide a completed segment numerator by that segment's
completed count/exposure according to the existing numerator's units. Fluence
is an integral and must not acquire that division. Particle-weighted moments
retain their particle-weight denominator. Different diagnostic classes cannot
share one normalization formula just because their buffers are adjacent.

Preserve the existing estimator for changing volume: if each sample deposits
W/V(t), output its time average. Dividing accumulated W by the final volume,
or replacing it with sum(W)/sum(V), estimates a different quantity. Actual
segment volume/area and physical metric must correspond to each sample epoch.

The present field-line sampler deposits many quantities directly to vertices,
and output uses the global `SampleCycleCounter` even when reading segments.
Therefore a new segment counter alone is insufficient. Segment accumulators
must preserve separate histories before output projects normalized segment
quantities to vertices. The projection must state its spatial weights and
conservation identity; histories with unequal exposure cannot first be summed
and then divided by one neighbor's age. Tests must preserve the equal-exposure
static limit and verify total represented particle density under the chosen
projection. Fixed-radius observers retain their own exposure; material-segment
averages must not be relabeled as observations at a fixed spatial point.

MPI reduces sampled numerators under the existing ownership rules. Exposure
is replicated geometry metadata or supplied once by its owner: summing an
identical count over all ranks would incorrectly multiply sampling time.
Validate rank agreement before publication. Completed output also retains the
geometry/window identity it describes; a later insertion cannot relabel it.

## Host order and restart

The host chooses one full dt and a single geometry owner. At departure, line
geometry, physical primitives, shock trajectory/intersection, source window
and frozen background metadata all describe t_n. Execute particle operators
and sampling with this immutable realization. Once they have joined, advance
geometry to t_(n+1), perform boundary topology changes and publish coherent
target-epoch state. Then recompute shock intersections and distances, source
and observer mappings for t_(n+1), and publish the next frozen realization.
No geometry writer runs inside the particle read phase. Any chosen departure-
frozen temporal splitting and its accuracy limit must remain documented and
convergence-tested; endpoint agreement alone is not temporal accuracy.

Finish old-topology particle-wave reductions and any old-generation operator
before insertion/removal, or explicitly remap its pending contributions by
stable segment identity. Existing host post-particle turbulence code cannot
receive newly indexed arrays containing contributions from old indices. This
ordering must preserve the approved physical operator splitting and keep
outgoing wave/particle losses in the same ledger.

Checkpoint schema must include geometry epoch/generation, stable material
labels, all collecting/completed counters and windows, and next-identity state.
Restore exactly; do not fabricate old history from a global output counter or
from which segments happen to be nonempty. Pre-change restarts need an explicit
compatibility policy and cannot be silently treated as having valid exposure.

## Acceptance gates and implementation order

1. Add production geometry/whole-segment-domain kernels and tests: uniform
   translation, every endpoint moved, curved field tracing, multiple prefix
   insertions, antiparallel B orientation, unreachable/degenerate boundaries,
   straddling segments and both-endpoints-outside chords.
2. Add native prepare/commit tests: callback failure leaves old state intact;
   material coordinates/identities survive insertion; removal deletes only
   attached particles; retired samples/loss ledgers balance; repeated commits
   reject; closed loops and complete removal have explicit outcomes.
3. Implement collecting/completed segment exposure, sampler/output migration
   and tests. Test different birth times, empty segments, multiple species,
   all buffer swaps/resets, zero exposure and repeated output. An independent
   constant-density oracle must give identical density for old/new segments.
4. Connect the analytical host owner; verify shock-surface residual and the
   material-relative radial crossing speed V_sh-U against independent oracles.
   Test actual update ordering, not only the pure intersection kernel.
5. Qualify MPI ownership changes, OpenMP behavior, checkpoints/restarts and
   imported-geometry ownership. Exposures and physical losses cannot depend on
   decomposition; keyed particle randomness cannot change because of indices.
6. Run strict warnings/sanitizers and the bounded transport, sampling, native
   SEP and BATS-R-US regression gates. Preserve references and tolerances;
   numerical differences are evidence to inspect, not permission to bless a
   new reference automatically.

New code must carry detailed comments at the callback boundary, staged commit,
inner trace, outer nearest-point test, S remap, particle deletion, sampling
increment/swap/reduction and restart decoding. Explain ownership, SI units,
frame/epoch, conservation and failure behavior there; do not merely narrate
loops. Keep this design and the relevant README sections synchronized with the
implemented and actually tested scope.
