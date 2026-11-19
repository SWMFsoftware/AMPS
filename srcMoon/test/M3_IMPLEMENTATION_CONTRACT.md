# Approved M3 implementation and local tracing contract

2026-11-17; starting revision `74974c2b2fc1b4cf09a73c320064a22370c0ca94`.
The user approved implementing `M3_COMPLETION_PROPOSAL.md`, then explicitly
required **first- and second-order particle-tracing movers with local time**,
cut-cell collisions, lunar/Earth gravity, radiation pressure, Lorentz and
Coriolis forces. These additions supersede any common-step-only restriction
on the movers in the proposal. Missing authoritative product kinetics remain
missing; approval does not fabricate those inputs.

## Authorities and implementation scope

Source roadmap SHA-256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
External roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source copy adds approved development-order/resource amendments; native
data labels/external AGENTS control data interpretation, and unchanged formal
qualification predecessors still control milestone PASS. This work develops
production capability and local/supporting verification before M1/M2 closure.

The source inventory, entry points, configurations, data/oracles and per-ID
contracts are in the approved proposal. Before individual production changes,
the owning component also records its detailed contract. New rate and event
services use the same production implementations exercised by tests. Existing
Na behavior, numerical gates, raw observations and reviewed references remain.

Parallel ownership: force-cut-cell mover/configuration/native fixtures; frozen
reaction rate registry and CEX kernels; conservative reaction transactions;
parent integration/configuration, gas geometry, global runners and evidence.
No concurrent configuration/build of the active checkout is allowed. Build
and full-core verification use isolated directories under `/tmp`.

## Local tracing and physical time

Each mover accepts a caller-supplied physical interval in seconds. Elapsed,
remaining and absolute epoch information survives face crossings, reflections,
reactions and MPI continuation without being rescaled as a stationary flux
weight. A daughter receives the remainder of its parent's actual interval;
its next reaction clock is fresh, not a copy of consumed optical depth.

Normal surface-source boundary injection must also reach the selected Moon
tracer/controller. The Moon boundary wrapper preserves the exact old function,
flag and optional face for legacy receipts. New tracers use actual entering-
root classification, including a curved return during the first segment,
rather than skip every internal boundary for that whole segment. A new source
uses its scoped moving-boundary departure epoch, or step-start plus
`local_species_dt - supplied_remaining_seconds` for the stationary generic
source. No second birth offset is drawn. Already initialized histories retain
their epoch. Native tests must invoke the actual configured boundary macro.

Steady/frozen-environment local tracing and a common-time transient campaign
are distinct configurations. No evolving environment may silently evaluate
different particle epochs at one unrelated global SPICE epoch. Expose explicit
time context to force/event providers; refuse unsupported temporal modes
before stepping. Force tests isolate terms analytically, then exercise the
actual configured Moon sum and selected mover. Neither order duplicates the
Moon force laws. Numerical controls, curved-path intersections, ownership and
local-time continuation have separate independent tests and documented errors.

The initial Moon adapter also requires `m2_environment=disabled`. Current
frozen ephemeris updates have zero global thermal epoch duration, while local
particle histories continue by their own physical intervals. No shared dynamic
thermal/reservoir clock mapping exists for those asynchronous histories.
Parser and direct installation must therefore reject the combination rather
than infer a prescribed-stock subset or silently freeze a transient state.
This does not remove the legacy M2 path or modify its numerical/physics gates.

## U09 gas photon visibility extension: pre-implementation specification

- **ID/name:** U09, volume-point finite-disc photon visibility sub-contract.
- **Purpose/capability:** correct gas absorption geometry using the production
  terrain/Earth masks, distinct from facet irradiance or projected disc area.
- **Entry points:** `SurfaceEnvironment::EvaluateSolarDisc`, its production
  `TerrainRayTracer`, and the proposed explicit volume-origin request mode.
- **Configuration/inputs:** deterministic uniform-radiance apparent Sun;
  explicit radii, distances and quadrature; volume point starts exactly at the
  particle coordinate, no facet offset/horizon or ignored BVH face. Old surface
  mode unchanged. A downward disc through a genuine triangle must be fully
  blocked; supplying a finite self-face is refused transactionally.
- **Independent reference:** solid-angle cap `2*pi*(1-cos(alpha))`; full/half
  visible angular symmetry; concentric occultation
  `(cos(beta)-cos(alpha))/(1-cos(alpha))`. No observation is manufactured.
- **Units/frame/time:** SI metres, steradians and dimensionless fractions;
  one fixed Cartesian frame and explicitly fixed geometry, no SPICE dependency.
- **Deterministic/stochastic:** no random draws; predeclared analytical cases.
- **Procedure:** evaluate full disc, half-disc surface horizon, wide-disc
  concentric Earth mask and volume origin, including reversed dummy normal;
  verify invalid mode/offset combinations refuse transactionally.
- **Metrics/acceptance:** full/half fractions and solid-angle formula within
  `2e-13`; concentric scalar fraction absolute error `<=2e-6` at explicit
  512-by-128 quadrature and smaller error than a 32-ring grid. Original U09
  gates are unchanged. Positive difference from projected-area fraction must
  be resolved for the deliberately wide numerical disc.
- **Expected PASS:** gas fraction follows integral of radiance over solid
  angle; default surface results and self-intersection policy remain intact.
- **Failure modes:** raw area substituted for gas flux, missing Jacobian,
  independently multiplied masks, facet cosine/albedo, displaced volume point,
  bad normalization/domain or state mutation after failure.
- **Required data/hashes:** no external observations; analytical fixture and
  production source/build hashes preserved in the normal U09 artifacts.
- **Predecessors:** original U09 dependencies retained; this local extension
  does not qualify D02, D04, I15 or I21 without their separate gates.
- **Artifacts:** actual production probe metrics, stdout/stderr, result and
  source/build manifest; comparison/reference values retained by the U owner.
- **Status:** PASS only evaluated local geometry; FAIL evaluated disagreement;
  ERROR build/runtime/infrastructure; SKIPPED missing configured dependency.

## Execution and evidence

One-hour development packages include setup/builds/scoring/plots/restoration;
default linked fixtures use eight MPI ranks/OMP=1. Reviewed nightly remains
four ranks. Remove `build/` before required normal Config.pl / `make -j`.
Require successful `make -j test_Moon` with zero-byte diff. Core changes must
preserve previous PASS identities under the full core owner. Longer formal
ensembles require the separately declared cluster tier; resource shortage
causes a detailed ERROR, never shortened physics or weakened confidence.

Keep STATUS.md and relevant READMEs current. No commit, push, staging, raw-data
modification or automatic reference replacement is authorized in this package.
