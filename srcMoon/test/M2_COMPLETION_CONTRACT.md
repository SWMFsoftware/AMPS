# M2 remaining implementation contract — 2026-11-16

The user requested completion and testing of unfinished M2 code. This authorizes
implementation and supporting verification under the existing development-order
exception. It does not waive M1/I13, invent D02 generation settings, approve an
I17 threshold, expose December residuals, or qualify a missing observation.
The initial source/data roadmap hashes and source revision are retained in the
new `test_output/srcMoon/m2-completion-*/initial_state.json`. Source amendments
control execution; the common scientific gates and native data labels remain
authoritative. Raw and canonical processed validation data are read-only here.

## Inventory before changes

* `TerrainRayTracer::IsSunVisible` is a point-direction BVH query. Production
  thermal forcing additionally computes opaque-Earth circular-disc overlap,
  but terrain remains a binary central ray. The existing I15 wrapper therefore
  stops before runtime. D02 exact temporal/height/ephemeris/Earth conventions
  and an independent epoch benchmark are unavailable.
* `RunI17`, `prepare_i17_input`, and `evaluate_m2` implement independently spun-up
  June-only evolution, three-step numerical qualification, diagnostics and a
  review-only threshold candidate. They intentionally cannot score December.
  No reviewed frozen observational threshold exists.
* `run_m2_cluster.py` already executes the approved fixed 288-case numerical
  design and scores I19 with an independent evaluator. It preflights same-run
  predecessors, thermal initialization, statistical evidence and cost. Normal
  I19 remains a capability SKIP; I20 has no independent polar comparison path.
* Current local/linked tests use compiled production objects. U09 covers point
  rays, U11 conduction/eclipse limits, U12 kinetics, U13 retention. Existing
  protocol fixtures cover campaign/scorer/cost behavior, not observations.

## Test specification

| Field | Contract |
|---|---|
| IDs/name | Extend U09 and existing I15/I17/I19/I20; no new permanent IDs. |
| Purpose | Finish production finite-disc terrain illumination and fail-closed formal orchestration/scoring. |
| Production entries | TerrainRayTracer; production thermal illumination; RunI15/RunI17; normal production surface/mover/source/cohort paths; existing long numerical executor/scorer. |
| Configuration | Explicit optional `.in` controls; absent controls preserve M0/M1 and existing demonstrations. Eight MPI ranks/OpenMP=1 for selected linked development cases. |
| Inputs | Existing native LOLA/SPICE, immutable D02/D03 metadata, reviewed external contracts where actually available; independently declared analytical geometry and parser fixtures for local verification. |
| Independent oracle | Closed-form unobstructed/half-disc/horizon/eclipse geometry, conservation identities, existing independent numerical scorer; native labelled observations only for formal scoring. Infrastructure fixtures are never observational evidence. |
| Units | Metres, radians, dimensionless visibility and incidence, TDB seconds/explicit UTC, K brightness temperature, physical atoms and kg. |
| Frame/time | MOON_ME_DE421 terrain; explicit SPICE aberration and averaging epochs in reference contract. No inference from label acquisition dates. |
| Deterministic/stochastic controls | Deterministic disc quadrature with explicit resolution and refinement diagnostics; preserve approved numerical seeds, weights, confidence and physical horizons. |
| Procedure | Implement/test local geometry; validate/stage external contracts without scoring; wire linked exports and independent evaluators; exercise missing/corrupt/wrong-selection/predecessor rejection; run bounded production development and necessary baseline regressions in an isolated checkout. |
| Metrics | Disc area and incidence separately from any-disc occurrence; averaging counts; numerical refinement, exact site/time identity; existing thermal/numerical/inventory metrics; explicit missing evidence. |
| Acceptance | Existing scientific thresholds unchanged. Local analytical error tolerances explicitly numerical. A finite ray grid is a resolved numerical approximation, not exact any-disc visibility; formal use needs its own resolution qualification. |
| Expected PASS | Only the declared evaluated analytical/production or infrastructure contract. Formal I15/I17/I19/I20 additionally require every original gate. |
| Failure modes | Area/occurrence confusion, correlated terrain/Earth masks multiplied independently, false hold-out exposure, stale/fabricated freeze, mismatched sites/times/builds, repeated-state mutation, numerical-grid/cost/metadata gaps. |
| Required data/hashes | D01, D02, D03 and SPICE as applicable; retain exact files/hashes in each manifest. Missing reference contracts are SKIPPED; altered or internally contradictory evidence is ERROR. |
| Predecessors | Supporting tests use local/M0 contracts; formal results require same-campaign I13 and the individual I-series predecessor chain. |
| Artifacts | Input/source/build/executable/data identities, streams/receipts, per-stage numerical CSV, exact scored rows and PNG/EPS, resource logs, all attempts, pending questions and completion report. |
| Status semantics | PASS evaluated/met; FAIL evaluated/missed; ERROR parser/build/runtime/corruption; SKIPPED unavailable scientific input/predecessor. Never promote support/preflight to formal PASS. |

## External evidence and remaining investigator decisions

1. Product-specific D02 timing/geometry and independent epoch truth cannot be
   manufactured. Implement explicit interfaces and leave those runs SKIPPED
   until authoritative evidence is supplied.
2. The user approved independent observation-free December periodic spin-up.
   Its separate explicit UTC interval is still required before execution;
   no new interval, state or tolerance is silently chosen.
3. I17 observational acceptance requires later investigator review of June.
   Implement validation of that frozen artifact; do not create an approval or
   populate its scientific limits automatically.

The development campaign must include necessary builds/plots/reporting within
its selected one-hour limit. Long qualification stays explicit cluster work;
do not launch it without measured feasibility and all scientific inputs.

The final native campaign completed on November 16 in 37.4 minutes, with
58 supporting PASS / two SKIPPED / no FAIL or ERROR. All ten M0 gates PASS;
the reviewed nightly diff remains zero bytes. The final Python-only export
extension was exercised in a separate 22.64-second eight-rank source-only
run and 27 passing protocol tests. Formal M2 remains unqualified. STATUS.md
section 5.14 and the completion report retain separate counts and exact scope.

## Production PSD wiring extension identified during review

The core PSD kernel currently multiplies the actual surface population,
species cross section and 1-AU photon flux by a logical radial cosine and
Moon-centred inverse square. The source CDF, finite-inventory hazard and
parallel surface exchange share that kernel, but the total-rate selector has
a cache. Consequently a correct thermal shadow or closed ledger alone does
not establish correct illumination of the PSD source (R2.2).

The minimal correction is an optional NULL-default runtime direct-photon
factor hook in the SAME production PSD kernel, returning projected incident
photon flux / its 1-AU normalization. The active M2 provider uses the shared
terrain/opaque-Earth finite-disc geometry and actual surface-point distance,
without albedo or a duplicate PSD law. Runtime totals must refresh through
the same per-element kernel so cache/CDF/finite hazards agree. Preserve the
complete legacy NULL arithmetic and all physical coefficients. Coefficients
are evaluated at the explicit source-BEGIN epoch for the approved operator
split; per-particle departure geometry does not silently change the source
CDF/hazard during that interval. Register only after real terrain/logical
state and start epoch exist, and test initial/cache/frozen/evolving behavior.

Existing U09/U22 and linked I18/I19 support can legitimately cover this
illumination/source wiring without repurposing an ID. Independent limits are
N*sigma*Phi with exact full/half/zero projected illumination, NULL legacy
identity, changed-population total refresh, and invalid-factor refusal.
M0 must retain its zero-byte nightly diff. Required shared-core regressions
must preserve previous PASS identities; source-hook support cannot qualify a
missing D02 comparison. No external-data or acceptance change is authorized.

Lifecycle invariant: disabling/reinitializing the production M2 state must
remove its own PSD provider before discarding that provider's cache. The
independent oracle is NULL dispatch and the original cached-total behavior.
U09 exercises this public initialization path without SPICE or observations;
an unrelated application's provider must not be removed by a Moon reset.

I20 surface-kinetics export extension: a checkpoint-disabled native run must
report each actual sampled surface element's Na/Ar residence time and its
classification against the explicit `cold_trap_retention_timescale_s`. The
operator-support runner independently checks the published Na Arrhenius and
Ar modified-Arrhenius limits at the exported kinetic temperature, using the
preserved input coefficients and SI constants. Relative 1e-11 bounds cover
floating evaluation only, not physical calibration uncertainty; no science
acceptance limit is changed. Positive infinity and zero retain their exact
cold/instantaneous kernel meanings. Subsurface rows must have empty kinetics.
Malformed fields are ERROR; finite computed-law or classification misses
are FAIL; these analytical checks do not qualify polar observations.

Approved follow-up: Diviner maps explicitly select `isothermal-grey`, and
production converts their bolometric brightness to kinetic K using the
configured emissivity. No unique physical inversion of an anisothermal
footprint is claimed. All prescribed boundary modes refresh the logical-cell
kinetic temperatures used by desorption/diagnostics at their declared state
epoch; impacts evaluate the SAME prescribed boundary function at their point.
Preserve the existing debug cosine law and cylindrical Earth-shadow kernel,
factoring its position arguments so an explicit epoch does not borrow stale
global midpoint states. The legacy wrapper retains its exact arithmetic.
Independent checks are brightness/kinetic roundtrip, selected-map benchmarks,
nonthermal cache/callback/source consistency, and unchanged nightly baseline.

## Actual initial thermal state and formal input preservation

The final code audit found that I19 accepted generic qualification evidence
while its native columns always started from the declared uniform scalar.
Formal staging also reused the I18 helper that replaces AMR/source controls.
Neither is proof that the reviewed physical problem was actually executed.

Production will expose named `declared-uniform` and `periodic-spinup`
initialization modes. The latter applies the EXISTING observation-free
periodic solver before source/weight queries and requires the explicitly
declared cycle end to equal simulation t0. No cycle, tolerance, physical
material or qualifying initial condition is selected automatically. A repeated
I17/I20 spin-up request may reuse only the same unadvanced completed state.
Native I19 exports the actual ordered full column state before cohort seeding
when requested, and records its applied mode/epoch/metrics. The reviewed
initial-state document must bind those actual bytes, executable, physics and
qualification evidence; unattached or unapplied evidence cannot PASS.

Formal case preparation preserves the reviewed AMR/source/thermal/clock
configuration, rebases paths without changing meaning, and varies only the
approved numerical controls, seed and diagnostic outputs. Missing initial
qualification stays SKIPPED; changed or contradictory state/physics evidence
is ERROR. Source and conservation/numerical criteria remain unchanged.
Initialization compute/storage cost is part of measured cluster feasibility.

The bounded native I20 fixture selects periodic initialization on its explicit
1-ms replay, sets simulation t0 to that replay end, and checks the actual full
column export and initialization metadata before/after its two queries. This
is production dispatch/provenance verification, not a lunar spin-up or a
qualified initial state. Model-state CSV files are model inputs/diagnostics,
never independent observations or physical validation truth. All scientific
predecessors and actual reviewed-state requirements still apply.
