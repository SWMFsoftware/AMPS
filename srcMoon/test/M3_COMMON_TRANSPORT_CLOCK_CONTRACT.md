# Explicit common-physical-time transport clock contract

Pre-implementation inventory and test contract, 2026-11-18. This is approved
software infrastructure, not a daughter-kinematics, source, sampling or
observational qualification. It introduces no permanent U/I ID.

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the external roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Their dated development/resource/status amendments differ; the approved
source M3 development exception controls this supporting implementation.
Common unchanged physics, acceptance and predecessor requirements, and the
external data rules, retain control of formal qualification.

## Inventory and missing capability

Canonical `input/moon.input` currently selects Na only and spatially global
per-species dt/weight. `PIC::TimeStep` samples once per outer iteration before
injection/motion. Density/LOS averages are iteration averages, with no
per-packet residence integration. The optional physical tracer advances
packets by their actual caller interval and retains epoch/history through
AMR/MPI, while refusing evolving external geometry. `amps_time_step` freezes
`OrbitalMotion::et` in the selected frozen environment. The existing
`ProductionInjectionEpochS` uses that external step epoch plus a departure
offset, so it does not supply an advancing common clock for repeated births.

Production entry points: `Moon::MoveParticle`, both actual source-boundary
macros/wrapper, `SurfaceEnvironment::BeginProductionParticleStep`,
`Production{Source,Injection}EpochS`, `amps_time_step`, actual PIC block/species
dt getters, first/second force tracers, normal ownership/continuation lists.

## Explicit mode and physical scope

`tracing_clock_mode = independent-history | common-time` is optional.
Absent means `independent-history`, preserving every existing receipt and
legacy behavior. `common-time` requires one of the two new Moon tracers,
`simulation_environment_mode=frozen`, `tracing_field_frame=lso`, M2 disabled
and the existing unsupported force/kinematic guards. It does not install
reaction providers, change fields, advance thermal state, move the lunar
pose or reinterpret frozen external coefficients as evolving ephemerides.

The common transport origin is the explicitly configured UTC interpreted by
normal CSPICE initialization. It is not J2000 zero or an inferred observation
epoch. The transport clock thereafter advances relative elapsed SI seconds;
external ephemeris/source-CDF/thermal/geometry coefficients remain at the
frozen external epoch. Physical birth time and the coefficient-evaluation
epoch are distinct named values.

Every owned allocated block/species dt must be EXACTLY the same positive
initialized production reference dt, on every rank, before each normal step.
This admits SingleGlobal or an actually identical per-species/local table;
a compile-time label alone is not proof. Mixed local/species steps are
refused before any step. The initial runtime scope also keeps that actual
common duration constant over outer steps; this avoids silently presenting
unequal-duration iteration averages as a physical time average. The pure
clock itself supports explicit varying intervals for numerical testing.
No minimum gyro/flight step or physical accuracy is selected by this check:
those remain separately declared convergence criteria.

## Clock, source and transaction semantics

The PIC-independent clock is initialized once, begins one immutable step,
commits it once after completed physical work, and can abort only an
unaccepted begin. A token binds sequence, relative begin/end, duration,
explicit origin and projected absolute ET. Relative elapsed time is retained
with compensation; repeatedly adding dt to a large absolute ET is forbidden.
Births are calculated relative to the origin BEFORE one projection to ET.
Nonfinite/negative/oversize source offsets, a nonpositive/unrepresentable
outer step, duplicate/mismatched commit, nested begin and chronology reversal
are errors with no clock-state mutation.

Normal `amps_time_step` begins before source callbacks and commits only after
`PIC::TimeStep` returns. Normal source-boundary wrappers use actual scoped
Poisson departure offsets or the supplied remaining interval (no second RNG
draw). Existing finite birth epochs must agree with that declared interval.
No unspecified/ambiguous initial population time is assigned automatically.

Common-only opaque Moon packet metadata is requested before buffer allocation
and copied by normal MPI/cloning operations; no shared PIC record ABI changes.
The metadata binds the current common target sequence/relative endpoint.
Conservative daughters inherit that target with the parent record, while the
backend/core retain their actual event epoch and remaining interval.
Commit requires no pending continuations, a valid common target on every
owned particle, and a completed actual core interval ending at the common
endpoint. Only the derived binary64 epoch-projection bound can align the
stored absolute endpoint; that bound never consumes optical depth or makes a
reaction due. Invalid completion terminates the run; it is not rolled back
or silently advanced. Checkpoint/restart of this new state is NOT provided.

Callbacks perform no MPI collective. Begin/Commit perform explicit normal
outer-step collectives, including ranks with empty particle/block sets. Any
configuration, source, resource or clock failure produces a detailed
`exit(__LINE__,__FILE__,message)` with the violated invariant and corrective
configuration/instructions. No automatic retry hides partially accepted work.

## Test specification

| Field | Contract |
|---|---|
| Name/ID | Supporting M3 common transport clock, no permanent new ID. |
| Purpose | Separate advancing transport/birth chronology from frozen external coefficients; prove explicit common actual dt and complete endpoint semantics. |
| Production capability/entries | New PIC-independent clock and normal Moon Begin/Commit/source/mover hooks listed above; same production object linked into normal application and thin probes. |
| Configuration | Explicit common-time plus new first/second tracer, frozen LSO prescribed fields, M2 disabled. Native fixture uses actual tiny PIC mesh/allocator/normal mover/owned lists on1/8MPI. Existing independent-history fixtures retain their mode. |
| Inputs | Explicit analytical start/step/source offsets, actual compiled species/masses/charges and actual per-block dt tables. Invalid controls are parser/temporal fixtures, never observations. |
| Independent oracle | Closed-form relative chronology n*dt and t_birth=t_begin+offset, exact endpoint complements, immutable external epoch; integer counters and independent long-double arithmetic for accumulated drift. No saved model output is physical truth. |
| Units/frame/time | Relative SI seconds plus projected SPICE ET [TDB s past J2000] from explicit start UTC; actual tracer LSO state; frozen external frame/pose unchanged. |
| Deterministic/stochastic controls | Fixed source offsets and declared cohorts; zero RNG needed by the clock. No physical channel/kinematic activation. |
| Procedure | Pure class and native parser/installer controls first; isolated normal clean SPICE build (remove build, Config, make-j) and thin native1/8MPI/order1/2 source/mover chronology checks if prerequisites permit. Original aggregate3600s includes setup/build/retries/plots/report/cleanup. |
| Metrics/criteria | Exact sequence/token/default-mode and dt/refusal semantics; compensated relative chronology; projected source/endpoint values within stated arithmetic projection bounds; every normal-step survivor/child target completes, no pending physical interval; external epoch/fields remain unchanged. No sampling/physical threshold is invented. |
| Expected PASS | Declared clock/configuration/native wiring controls evaluated and met; existing independent/default behavior preserved. |
| Failure modes | Reused frozen birth epoch, stale/mixed common dt, incomplete/MPI-stranded packet, ambiguous source epoch, stale/forged token, duplicate commit, numerical epoch drift, hidden rollback, half-step sampled as complete, incompatible restart. |
| Required data/hashes | No observations or external rate tables. Normal SPICE build dependencies; native explicit ephemeris snapshots where declared. Retain all actual source/compiler/config/input/generated/library/executable identities. |
| Predecessors | Existing local-time core tracing and packet/backend controls. Formal M1/M2/I21-I24 order unchanged. |
| Artifacts/status | Per-case stdout/stderr, result/run manifests, source/config snapshots, CSV/PNG/EPS analytical chronology comparisons, budget/resource receipts. PASS/FAIL/ERROR/SKIPPED retain normal meanings; no source/preflight-only result is linked PASS. |

## Remaining scope

This common-time path addresses an explicit equal-physical-step population;
it does not qualify arbitrary asynchronous sampling. Physical source rates,
external drivers, recoil/bath laws, observation operators and numerical /
statistical convergence still require their existing gates. No observational
claim, hidden time origin, checkpoint support or physical accuracy decision
is added by selecting common-time.


## Analytical controller/daughter composition extension — 2026-11-18

Extend the existing native clock fixture, without a new ID or production
change: keep its first three source/chronology steps inert via an explicitly
resolved numerical Na state. A separate fourth step starts at relative3s,
creates tag7 at offset0.25s and gives it one immutable constant1s^-1 channel
with optical depth0.4. Expected event ET1003.65, unused0.35s and common end
ET1004 are analytical control values from declared origin1000.

Use the actual `Controller`, `MakeEventPlan`, `ParticleTransaction` and queue.
The declared nominal-heavy-mass Na+He+→Na++He zero-Q mechanical control has
explicit co-moving incoming bath velocity and signed mass/momentum/energy
moments BEFORE products are supplied. It installs no physical CEX/recoil law.
The source has represented N=1; exactly one Na+ and one He carrier must each
represent N=1, inherit common target4 through the actual byte clone, complete
unused time and pass common Commit. Exact nuclei/charge and accepted event
count remain integer criteria. Existing backend floating reconstruction bounds
apply only to rounded W*(N/W); the configured1e-12s interval consistency and
derived absolute epoch projection envelope are unchanged, not fitted to output.

This adds one named global criterion per native job:20×(1/8MPI×bothorders)=80
control evaluations. Per-rank and aggregate daughter CSVs must contain two
actual global rows per successful job with finite matching state/time/number
fields. Original76-PASS evidence remains unchanged and is not relabelled.
All production defaults/providers and normal32-control startup coverage stay
unchanged. No new build or physical dataset is needed for this thin extension.
