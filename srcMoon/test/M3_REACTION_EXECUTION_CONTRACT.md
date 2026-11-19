# M3 conservative reaction primitives — implementation contract

Written before production changes, 2026-11-17, source `74974c2b2fc1b4cf09a73c320064a22370c0ca94`.
This work implements reusable production primitives, not complete I21–I24 or
normal-mover chemistry. The approved M3 proposal controls design; source/data
roadmaps retain SHA-256 `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.

| Specification field | Contract |
|---|---|
| ID / name | Supporting conservative transaction subcontracts of existing U15/U17/U24 and I21–I24; no new permanent test ID. |
| Purpose | Verify integer material stoichiometry/state identity, explicit bath bookkeeping, represented-number transfer, local-time continuation and capacity-safe transaction ownership. |
| Production capability / entry | New PIC-independent `Moon::Species::ReactionExecution` services in `srcMoon/species/ReactionExecution.h/.cpp`; parent PIC adapter will consume the same event plans. |
| Inventory | Upstream `ChemicalSolver::cReactionNetwork` validates structural term/species bounds but allows noninteger coefficients, indexes only unit-unimolecular reactions and rejects signed excess energy. `cChemistryDiagnostics` explicitly documents number-rate moments; it must not silently receive transient-number moments. |
| Reuse | Use the upstream network as a structural stoichiometric view; enforce integer chemistry and tracked-parent/bath distinctions separately. Keep dedicated transient-number ledger; an eventual diagnostics adapter must explicitly label any number-to-rate conversion. No generic continuum rate formula or product distribution is duplicated here. |
| Configuration | Registry and caller-supplied per-event kinematics/bath conservation contract. No built-in lunar recoil/angular distribution, solar spectrum, bath density or cross-section table. |
| Inputs | Explicit species/state identities, integer nuclei and elementary charge, tracked/bath reaction terms, per-parent physical interval, already evaluated nonnegative hazards, open-interval draws, positive reference/individual weights, qualified product kinematics and signed bath exchange. |
| Oracle | Integer stoichiometric identities; exponential optical-depth/competing-hazard equations; physical-number identity W*c; analytic closed-system and explicitly declared bath fixtures; transaction state snapshots. Fixtures are controls, never observations or kinetic evidence. |
| Units | Physical epoch/duration/elapsed in seconds; hazards s^-1; optical depth dimensionless; represented number in physical particles; mass kg, velocity m/s, momentum kg m/s, energy J. |
| Frame/time | Every transport kinematics contract names one event frame. Each immutable trace interval retains its own physical start, duration and elapsed time. No global/common dt assumption; multi-species asynchronous scheduling policy belongs to the adapter. Products inherit the event position/physical interval, never their nominal species dt as a new birth interval. |
| Stochastic controls | No RNG in the service. Caller supplies strict 0<u<1 draws and fresh daughter clocks. Same physical continuation preserves residual optical depth; real products get fresh clocks. No draw occurs on capacity retries. |
| Procedure | Validate catalog/channel with overflow guards; advance constant hazards analytically or consume a caller-integrated no-event exposure; create a complete immutable plan; reserve all storage; stage/commit through an atomic backend contract; rollback reservation/staging on refusal/error. |
| Metrics | Exact nuclei/charge/state identities, unchanged outputs on invalid input, stable probability/time/selection, local interval preservation, packet-weight and moment residuals, no plan changes/redraw on retry, no double commit. |
| Acceptance | Integer composition/charge closure exact; event quantities match independent analytic calculations within explicitly declared FP comparison bounds. All transport product occurrences present and correlated kinematics caller-qualified. Closure bounds are explicit inputs, not tuned defaults. Failure cannot partially consume a parent or diagnostics. |
| Expected PASS | Only evaluated primitive/API/transaction controls; no real lunar kinetics, PIC integration, physical event-chain or formal M3 PASS. |
| Failure modes | Missing/duplicate state, invalid/overflowing stoichiometry, missing kinematic authority/frame/bath contract, zero/nonfinite weights, wrong product multiplicity/reference weight, consumed daughter clock, stale local interval, invalid uniform/hazard, missing reserve capacity, backend staging/commit refusal. |
| Data / hashes | No authoritative kinetic/observational table is consumed. Record exact production/header/upstream-oracle/test/compiler hashes in supporting evidence. External raw data are untouched. |
| Predecessors | Approved early M3 implementation/local verification. Formal M1/M2 and I21–I24 predecessors remain unchanged and unqualified. |
| Artifacts | Standalone production-service compiler log, probe stdout/stderr, named control results and source/compiler identity; README and retained supporting result.json. |
| Status semantics | PASS evaluated primitive control; FAIL evaluated mismatch; ERROR compile/runtime/parser/infrastructure; SKIPPED unavailable compiler. Missing kinetic policy explicitly refuses transport, while rate-only plans cannot commit particles. |

## Atomic backend contract

The PIC-independent layer cannot undo arbitrary undocumented mutations by an
external backend. The required interface therefore separates reserve, staging,
atomic commit and rollback. Reserve/staging may hold private slots but cannot
publish particles or change parent/physical ledgers. A failed commit must be
rollback-safe. The service calls rollback for any deferred/failed attempt,
retains the immutable plan, and marks committed only after backend success.
The parent PIC adapter must verify this contract through actual PIC allocation,
lists and diagnostics before any production reaction activation. A successful
in-memory control backend is not PIC integration qualification.

## Scope boundaries

This layer does not infer kinematics from total rates/mean excess energies.
Caller-qualified classical product velocities and explicit signed bath moment/
energy/model-mass changes are required for transport plans. Those approval
identities are supplied by an independently validated production policy, not
made authoritative by a boolean alone. Rate-only planning is explicitly unable
to allocate/remove model particles. Classical moment bookkeeping itself is not
a new recoil distribution; a future relativistic treatment requires an explicit
extension rather than silent use outside the caller-qualified domain.
