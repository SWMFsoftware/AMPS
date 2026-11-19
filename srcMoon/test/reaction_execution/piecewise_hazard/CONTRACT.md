# Explicit piecewise-affine hazard history: pre-implementation contract

Written 2026-11-18 against source revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`, before this service or its probe.
This is supporting numerical coverage of the existing chemistry meanings,
not a new U/I identifier, a physical provider, or formal M3 qualification.

## Inventory and authority

`species/ReactionExecution` already implements compensated optical clocks,
constant competing hazards, explicit channel selection and conservative event
plans. `MoonReactionController` currently accepts constant histories only;
neither module provides a trajectory-dependent exposure law. The new pure
`species/PiecewiseHazard` service integrates an **explicitly supplied** affine
representation without changing either existing module, PIC, an observation,
a source, a clock or a random generator. A future caller must independently
qualify the physical rate/exposure provider and its coupling to trajectories.

Source roadmap SHA-256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
External roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source adds approved bounded-development and development-order amendments;
external labels/AGENTS govern real data. Original formal predecessors and
physics thresholds remain unchanged. No validation package is consumed here.

## Representation and exact endpoint convention

Each supplied segment has finite `0 <= beginS < endS` and an explicit,
nonempty list of unique reaction IDs with nonnegative finite endpoint rates.
Every segment retains the same ID order. Segments must be ordered without
overlap; gaps are permitted in the representation but never filled. Between
the TWO endpoints of ONE segment the caller explicitly declares an affine
hazard. This is a mathematical provider representation, not permission to
interpolate measurements, bridge missing records, or infer rates.

Rates may jump at a shared boundary. At an event whose first-hitting time is
exactly a segment's right endpoint, its **LEFT-limit** rates own the event.
The parent approved this deterministic representation convention; it is not
a claim that discontinuous physical laws have a uniquely defined point value.
An exact endpoint with zero left-limit total still returns the event and all
zero rates. A downstream channel sampler must explicitly refuse that case;
this module neither substitutes right rates nor normalizes an empty CDF.
A strictly larger optical depth continues into the next declared segment.

The entire requested window must be covered before a result is published,
even if a possible event precedes a later gap. Missing start/end support or a
gap returns `IncompleteSupport`, not no-event. Zero-length windows require a
covered point and have zero exposure. No extension beyond the final segment
is inferred. Inputs/results are caller-owned and never modified on invalid,
incomplete or numerical-error returns; complete no-event is a distinct result.

## Test specification

| Required field | Contract |
|---|---|
| ID/name | Supporting piecewise-affine exposure/event locator under existing chemistry verification; no permanent ID. |
| Purpose/capability | Integrate declared nonnegative competing hazards and locate the first optical crossing, independently of PIC, clocks and RNG. |
| Production entry points | `PiecewiseHazard::{ValidateHistory,Integrate,LocateEvent}` in the maintained production module. |
| Configuration | Pure C++17; explicit caller `NumericControls.maximumRootRefinements`, with no default or event-acceptance tolerance. No Config.pl or native application rebuild. |
| Inputs | Explicit segment boundaries, channel IDs/endpoint rates, covered relative-time query window, strictly positive finite long-double optical depth. |
| Independent oracle | Closed-form constant/rising/falling quadratics and independent Python Decimal analytical integrals/root evaluation. Explicit analytical fixtures, never observations. |
| Units | Physical relative seconds; rates s^-1; integrated exposure/optical depth dimensionless. |
| Frame/time | Caller-declared relative physical-time coordinate. No UTC/TDB conversion, trajectory, global PIC epoch or spatial frame is inferred. |
| Deterministic controls | No stochastic operation. All fixtures and mathematical values predeclared; no channel selected and no optical clock debited. |
| Procedure | Compile actual production source with a thin I/O adapter; compare exposure, event time, segment ownership and instantaneous rates; inspect negative-control output preservation; retain source/compiler/executable identities. |
| Metrics | Exact status/segment/ID decisions; dimensionless exposure, event seconds, instantaneous per-channel/total rates; partition equality; unchanged outputs/history on refusals. |
| Acceptance | Exact decision/endpoint/gap/refusal identities, without tolerance. Nonzero numerical values within `128 * native long-double epsilon` relative to the independent analytical reference; exact zeros remain exact. This is a declared arithmetic bound, not a physical uncertainty. The locator certifies the first representable long-double timestamp that reaches the analytic exposure, with no time/depth epsilon. |
| Expected PASS | Every declared mathematical/transactional control evaluated and met. Photon/electron/CEX physical providers and formal I21-I24 remain unqualified. |
| Failure modes | Negative/nonfinite rates, duplicate/changed channel IDs, overlap/reversed windows, missing support, invalid depth/controls, arithmetic range loss, exhausted root certification. No partial output publication. |
| Required data/hashes | No observations or external package. Production/test source hashes, roadmap hashes and compiler/executable identity retained. |
| Predecessor gates | Existing development-order exception permits this local work. Native controller/mover/provider and formal predecessor gates are separate. |
| Artifacts | `result.json`, `run_manifest.json`, compile/stdout/stderr logs, fixture inputs, per-case comparison CSV and mathematical-reference PNG/EPS. |
| Status semantics | PASS evaluated local criterion; FAIL evaluated disagreement; ERROR compiler/runtime/parser/corrupt evidence; SKIPPED missing compiler/dependency. Optional missing Matplotlib only skips visualization. |

## Numerical method and limitations

Exposure is the analytical affine trapezoid, not time quadrature. Positive
factor products are exponent-scaled to avoid avoidable intermediate overflow
or underflow; returned areas use compensated summation. All supplied finite
binary values are exactly rational. A private standard-library-only positive
integer/rational layer therefore compares cumulative exposure and optical
depth EXACTLY. This is necessary at a falling-to-zero endpoint: a rounded
`1-epsilon^2` must not become an early event at depth1. No external arithmetic
library, tolerance, negative-rate clipping or empirical correction is used.

An analytical quadratic seed is certified with those exact comparisons;
termination requires an exact rational root or adjacent representable
timestamps. Exact comparisons also preserve a positive segment exposure
below long-double range while locating an event in a later segment. Only
the seed/display accumulator may round; it never decides event membership.
An unrepresentable complete returned integral or instantaneous rate is
refused. The explicit refinement limit bounds work; reaching it is a numerical
error, never permission to move an event earlier. An overflowing whole-window
integral can still have a representable finite event, which remains locatable.

The local campaign, including implementation, compilation, checks and plots,
must finish within one hour. No physical law, observation, measured-data
interpolation rule, callback activation, PIC scheduling or conservation
criterion is introduced by this module.
