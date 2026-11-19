# M3 production trajectory exposure contract

This package implements item 3 of the user-approved M3 completion work. The
approved numerical protocol is recorded separately; scientific source/state
admission and formal I21–I24 predecessors remain required.

| Field | Contract |
| --- | --- |
| ID/name | Supporting extension of existing U14/U15/U16/U17 and I21–I24 coverage; no new U/I ID |
| Purpose | Integrate source hazards over the actual accepted force-tracer arc, including point-Sun terrain and opaque-sphere visibility transitions |
| Production entries | `TrajectoryExposure::BuildPhotonHistory`; explicit `MoonReactionController` trajectory provider; existing `TerrainRayTracer::IntersectRay` |
| Configuration | Frozen external geometry; explicit epoch/frame, reference distance, rate table, numerical depth/time/work bounds; legacy constant-history mode unchanged |
| Inputs | Actual SI position, velocity and starting acceleration; complete channel rates; actual mesh face geometry and production ray callback; explicit Sun/Earth/sphere positions |
| Oracle | Closed-form constant/linear-flight inverse-square integral; analytic shadow intervals; independent conservative interpolation-error identity; native same-callback accepted-flight checks |
| Units | m, m/s, m/s², s; rates s^-1; exposure dimensionless |
| Frame/time | All geometry on one explicitly named frozen Cartesian frame; transport epoch remains the real particle epoch; body rotation supplied explicitly by the caller |
| Controls | No RNG during geometry, integration, preview or retries; one channel draw only after accepted exposure crosses the preserved optical depth |
| Procedure | Enumerate all polynomial silhouette crossings before integrating; classify open intervals through the production ray service; adapt only the smooth inverse-square factor; use identical preview history on accepted prefix |
| Metrics | Integrated hazard, analytic affine-interpolation bound, shadow-root time envelope, instantaneous channel rates, accepted optical-clock continuation and draw counts |
| Acceptance | Existing U/I gates unchanged; requested numerical error/work bounds must be met or evaluation refuses; analytical errors must lie in declared numerical bounds |
| Expected PASS | Only evaluated numerical/software criteria pass; no physical profile or observational claim is inferred |
| Failure modes | Missed short shadow, physical rate frozen along flight, interpolation across visibility jump, incomplete channel family, preview RNG, accepted/debited arc mismatch, exhausted resource bound |
| Required data | None for analytical fixtures; actual production terrain/source profile for physical runs |
| Data hashes | Recorded in the caller's manifest; unavailable data has no invented hash |
| Predecessors | Local software controls can run early; formal I21–I24 retain original same-campaign gates |
| Artifacts | Source/build identities, exact commands, stdout/stderr, result/manifest, per-case numerical comparisons; native package resource/restore receipts |
| Semantics | PASS criterion evaluated/met; FAIL evaluated mismatch; ERROR build/runtime/corruption; SKIPPED unavailable dependencies/data/provider |

The force tracer already represents each accepted position arc as
`x(t)=x0+v0*t+a0*t²/2`. The exposure service uses that same representation;
it neither adds a force law nor replaces the first/second-order velocity
quadrature. Photon visibility is point-Sun as approved P01, not finite-disc
illumination. Smooth inverse-square exposure uses an affine approximation
whose interpolation error is bounded from the analytic second derivative.
Shadow root time uncertainty is reported separately, never hidden in the
interpolation tolerance. No empirical physics correction is introduced.

The controller may retain disposable numerical preview scratch. Scratch is
not particle/history state, contains no random draw and is never migrated.
Only accepted flight updates the optical clock/interval/ledger. An accepted
prefix must consume the exact history produced by its final preview; adaptive
retries may replace scratch but cannot debit exposure. An inconsistent or
missing preview is an infrastructure error, not permission to recompute a
different hazard after acceptance.

### Review corrections and numerical admission

The final root implementation uses `TrajectoryIntervals.h`: every arithmetic
operation rounds outward, all geometry coefficients retain enclosures, and a
complement interval is excluded only by a strict interval sign. Unresolved
brackets, including near-tangent/coefficient-rounding cases, are unioned and
their whole possible photon exposure is charged to the bound. The initial
root-count heuristic is rejected by preserved independent counterexamples.
The maximum geometry-interval visit count is an explicit required control.
The affine interpolation bound also includes endpoint projection rounding.
Actual BVH and mesh accuracy remain separately qualified by terrain gates.

The controller's provider returns `ExposureHistory` with an explicit error
bound (omission is NaN and refuses). Accepted prefixes conservatively debit
the full preview bound once. Exact compensated clock updates expose their
storage error bound; both cumulative numerical bounds persist through MPI and
outer-interval reopening, are passed to the physical event context and are
stored in the accepted `CommittedReaction` record. No preview publishes them.
Scalar `ExactPendingDepth` remains a legacy explicit projection check; actual
variable location consumes both components through an exact arithmetic receipt.
The native weak-hazard fixture is a regression for the earlier realistic
second-flight failure.

### Represented midpoint projection admission — independent review correction

The exact quadratic arc and the binary64 `Vector` passed to a visibility
callback are different numerical objects. A root-free cell may reuse its one
callback classification only when outward evaluations prove the same strict
sign for every terrain/sphere silhouette constraint at the exact midpoint and
at the represented callback point. If a sign is unresolved or differs, the
ENTIRE cell duration is added to the shadow-error union; no guessed position
epsilon or root-count estimate is used. Already charged root brackets are not
charged twice. The unchanged total absolute/relative budget must admit both
this exposure uncertainty and smooth quadrature error, otherwise evaluation
refuses transactionally. A translated, sub-ULP transverse-flight screen fixture
checks both a permissive-budget bound and a restrictive-budget refusal. This
accounts for projection into the declared callback coordinates; unknown BVH,
ray-intersection or caller frame-transform error remains a separate production
geometry gate and is not claimed covered by the temporal quadrature bound.
