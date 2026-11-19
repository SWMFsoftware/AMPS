# M3 independent software and model/protocol preparation contract

2026-11-18, before implementation. Source basis
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus the existing working changes.
The user approves continuing independent software and preparing concrete
model/protocol specifications. Existing approvals are retained; no physical
product distribution, mass/Q convention or observation is created by this
software authorization.

## Controlling contracts and approved scope

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the external roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source adds the approved early-development/resource amendments. Both
retain the same formal chemistry and predecessor requirements. Native data
labels and external AGENTS control data interpretation. The canonical data
repository is read-only. No roadmap milestone, permanent ID, acceptance
relaxation or automatic physical-provider approval is introduced.

Current production inventory: `ReactionExecution` owns integer material,
compensated optical clocks, constant-rate decisions and packet plans;
`MoonReactionAdapter` owns actual PIC atomic publication;
`MoonReactionController` accepts constant-history hazards only. The core
first/second tracers expose physical elapsed/remaining time. Frozen external
step epochs currently repeat source birth epochs. Normal sampling uses loop
counts; general asynchronous science remains unqualified. Existing M0,
controller and 843-case core preservation evidence must retain its exact basis.

This package develops (1) declared piecewise-affine hazard mathematics,
(2) a bridge into the same production optical-clock/channel selector,
(3) an opt-in advancing common transport clock independent of frozen external
ET, and (4) concrete physical-model and qualification-preparation documents.
It does not install a physical variable-rate provider, infer spectra/recoil,
enable missing channels, or claim formal I21-I24 execution.

## Shared test specification

| Field | Contract |
|---|---|
| ID/name | Supporting extensions of existing U14-U17 and M3 event/time meanings; no new permanent U/I ID. |
| Purpose | Verify independent software for time-dependent exposure, conditional process selection and common chronology; make unresolved physical inputs and qualification choices machine-readable. |
| Production capability | New standalone exposure/clock services linked into normal Moon library; existing optical-clock and CDF decision reused. Common-time source/step hooks opt-in. Missing physical providers remain disabled. |
| Entry points | `PiecewiseHazard::{ValidateHistory,Integrate,LocateEvent}`, `ReactionExecution` integrated-clock/selector APIs, `MoonTransportClock` and normal Moon source/step wrappers; preparation validator reads maintained model/protocol JSON. |
| Configuration | Named numerical controls; explicit common-time versus independent-history selection; frozen external LSO fields/pose, M2 disabled and one fixed actual species/owned-block step for initial common-time runtime. Independent mode retains its old layout/behavior. |
| Inputs | Declared affine segments with explicit boundaries/channel identities; long-double nonnegative rates [s^-1]; pending dimensionless optical depth; explicit clock origin/duration; versioned documents with source hashes and unresolved fields. |
| Independent oracle | Closed-form constant/linear integrals and event roots; high-precision evaluator; branch probabilities from instantaneous rates at the event; exact binary interval identities, immutable snapshot and chronology tests. No science output is used as truth. |
| Units | Relative physical seconds, external/transport ET seconds, rate s^-1 and dimensionless exposure; species/bath physical inputs retain their documented SI/state/frame conventions. |
| Frame/time | Relative arithmetic before projection onto absolute ET. Frozen external epoch is distinct from advancing common transport epoch. No external ephemeris or material motion inferred from advancing particle time. |
| Stochastic controls | Supplied open uniforms; pure queries draw nothing. Selection occurs once after located crossing. Protocol plans declare original-parent trials/streams, never MPI ownership or sibling records as independent trials. |
| Procedure | Contract first, standalone compiled software/control tests, real native common-clock source/step tests, registered global support, isolated normal clean configuration/build and unchanged M0/nightly checks. |
| Metrics | Integral/root/channel error; exact event/no-event/strict-future classification; clock freshness/consumption; unchanged outputs after refusal; all-rank chronology; correct schema/blocking/integrity status. |
| Acceptance | Existing physics gates unchanged. Analytical software comparisons use predeclared representation-level bounds; event crossing uses no optical-depth tolerance. Numerical exhaustion/gaps/domain violations refuse rather than fabricate coverage. |
| Expected PASS | Declared local/native mechanical criteria actually evaluated/met. Actual physical-model/preparation readiness remains SKIPPED where inputs/providers are absent. |
| Failure modes | Cumulative instead of instantaneous branch weights; erased tiny exposure; repeated clock draws; gap interpolation; incorrect exact-boundary ownership; frozen birth clock reused as transport time; noncommon runtime dt; fabricated readiness. |
| Required data | Standalone mechanics require no observations. Model specification audits actual selected D04/D05 source files and prior CEX references read-only; missing joint kinetics remain null/BLOCKED. |
| Data hashes | Baseline/evidence manifests plus exact model-spec source identities; no fabricated hash for unavailable arrays. |
| Predecessors | Current local clock/packet mechanics; formal M1/M2 and I21-I24 gates unchanged. Common-time verification does not waive thermal/surface or asynchronous observation-operator requirements. |
| Artifacts | stdout/stderr, result/manifest, source/compiler/input hashes, numerical reference/candidate CSV/PNG/EPS, resource histories, actual protocol/model documents and explicit missing-input reasons. |
| Status semantics | PASS evaluated software criterion; FAIL evaluated mismatch; ERROR build/runtime/corruption/contradiction; SKIPPED absent support/data/provider/qualification. Structural preflight success never becomes physical validation. |

## Exposure and clock bridge decisions

An affine history is an explicitly declared mathematical representation, not
permission to interpolate observations or missing driver samples. Segment
jumps are retained. At an exact first-hitting segment endpoint the declared
representation uses the left-limit rates; a zero total instantaneous rate
cannot be normalized into an invented channel. Unsupported intervals/gaps and
unrepresentable roots preserve outputs and report their actual cause.

The bridge reads the existing compensated optical depth in extended precision,
locates the event with the production history service and uses one production
instantaneous-rate CDF selector. The original constant-hazard service delegates
to that same selector without changing its existing arithmetic/order. No
integrated per-channel probability is substituted for lambda_i(t_event)/sum.
Clock mutation remains an explicit accepted-result operation, not preview.

## Protocol and physical-model documents

Physical models name their actual source states, energy/mass/bath conventions,
validity and missing numerical distributions. BLOCKED/null is required when
authoritative information is absent. A rate table or `qualified=true` alone is
not readiness. The preparation validator checks schemas, hashes and frozen
existing criteria but cannot manufacture an installed physical executor.

A concrete statistical/numerical case design can be proposed before outputs;
new precision/power/error-allocation choices not previously specified remain
marked PROPOSED pending review. Existing 2%/0.1%/3% and 95%/99% gates are not
changed. Common physical time is the approved initial qualification direction;
full H2O/OH and unresolved U24/forcing dependencies remain explicit.

## Builds, preservation and evidence

Use separate <=3600-second end-to-end native packages including setup, clean
builds, comparisons and restoration, eight MPI/OMP1 with explicit one-rank
controls. Remove `build/` before Config.pl/make. Preserve exact source and
ABI identity, old failed attempts and zero-byte reviewed `test_Moon.diff`.
No active-root configuration/build, staged file, commit, push, raw-data edit
or reviewed-reference replacement is authorized here. Shared PIC code remains
unchanged; new Moon-only source must retain current shared-core proof scope.
