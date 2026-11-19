# Opt-in constant-history reaction controller contract

Pre-implementation inventory/contract,2026-11-17. This module connects the
approved production primitives, actual PIC transaction backend and existing
trajectory hooks. **It adds no physical recoil distribution, rate table,
`.in` transport selector, observation, permanent test ID or default activation.**

## Inventory and implementation boundary

`ReactionExecution` already owns exact topology, compensated optical clocks,
combined constant hazards, represented-number plans and explicit bath ledgers.
`MoonReactionAdapter::ParticleTransaction` owns actual private PIC allocation,
atomic publication and detached retained-parent ownership. Core
`QueryTrajectoryEvent` previews without side effects; `AdvanceTrajectoryEvent`
debits only accepted physical flight; `ProcessTrajectoryEvent` commits an
event and returns Delete/Continue. The force tracer owns parent deletion,
remaining physical time, MPI exchange and continuation queues. Normal Moon
photochemistry is currently a separate whole-step list callback, which must
be bypassed/reinserted when the new controller is explicitly installed.

The controller is Moon-only. Existing source/runtime rate registries do not
provide complete daughter kinematics and are not attached automatically.
The first supported hazard scope is immutable **constant over each physical
particle/species-state history**, including across AMR/MPI/outer-step
continuations. A trajectory/illumination/time-dependent hazard request is
rejected, not frozen at its initial position. Its geometry/exposure-integral
provider remains a separately unimplemented physical contract.

## Test specification

| Field | Contract |
|---|---|
| ID/name | Supporting M3 combined constant-history native controller, extending existing chemistry/mover meanings without a new ID. |
| Purpose | Verify one combined optical clock and channel decision through actual normal Moon mover/PIC hooks and accepted conservative product publication. |
| Production entries | `MoonReactionController::{Install,PrepareParticle,PrepareProduct,Query,Advance,Process}`; `ReactionExecution::{ConstantHazardEventDelay,AdvanceConstantHazardClock,AdvanceConstantHazards,MakeEventPlan}`; `ParticleTransaction`; normal `Moon::MoveParticle`, `Moon::Chemistry::ProcessParticle`, core first/second tracer and detached queue. |
| Configuration | Default DISABLED; installation before buffer allocation, explicit immutable catalog/reactions/constant rates/evidence IDs/source resolver/open-uniform provider/kinematic-policy provider/ledger. Only selected first/second physical-time Moon tracers, MPI-only/OMP1, individual weightsON, generic PIC electron executorOFF, complete empty trajectory-hook set. No existing hook overwritten. |
| Inputs | Catalog states remain explicit even when several map one PIC species. Initial source resolver must select the state; fresh source clock is drawn once. Annotated daughters inherit adapter state/clock/interval and never use the source resolver or a newborn-dt factor. |
| Independent oracle | Closed-form constant competing hazards with known open uniforms/event times/channel intervals; exact integer nuclei/charge/represented-number/bath budgets; independent mechanical fixture kinematics; actual list/slot/state/epoch observations. No fabricated observations. |
| Units/frame/time | Rates s^-1, optical depth dimensionless, physical local time seconds, x/v/SI units and explicit event-frame identity. Actual force-tracer physical context controls event epoch; no global epoch substituted. |
| Deterministic controls | Explicit fixed analytical source/channel/product clock draws and counters. Query/refinement has no draw or mutation. Accepted crossing selects once and caches consumed clock/channel in opaque trivially-copyable particle bytes before MPI exchange. |
| Procedure | Separate bounded≤3600s package, clean isolated normal eight-species Moon build, actual normal PIC walker/tracer at1/8MPI. Two competing channels, three explicit source states, source/channel counts, partial-time product continuation, no-double-deletion/whole-step photolytic bypass, local-dt/weight differences, MPI crossing before an event, exact interval endpoint versus strictly future deadline, and expected resource-fatal subprocess. Exact geometry/endpoint coincidence requires its own explicit control; an earlier MPI crossing does not establish that coverage. |
| Metrics/acceptance | Exact topology/channel/source/accepted-event counts; recorded clocks/draw counts preserved. Represented N and mechanical ledgers within explicit rounding-level bounds; local physical event/remaining times within declared tracer time tolerance plus binary64 clock representation. No science acceptance threshold is relaxed. |
| Required data/hashes | No observations or external rate/kinematic data. Retain exact source/compiler/MPI/input/generated/library/executable identities and original manifests. Analytical topology is storage/controller-only, not real Na/CEX qualification. |
| Predecessors | Local primitive controls and real backend controls. Formal I21–I24/predecessor order remain unchanged. |
| Artifacts | stdout/stderr, result/run manifest, exact input/source snapshots, build/probe identities, per-rank states/draw/channel/event-budget/remaining-time CSVs, expected-fatal logs, comparison diagnostics. No earlier meaningful result overwritten. |
| Status | PASS declared local/native controller controls; evaluated invariant missFAIL; parser/runtime/resource/infrastructure failureERROR; unavailable capability/dependencySKIPPED. Expected resource-refusal control passes only on its specific detailed failure. Missing physical kinetic evidence remains NOT QUALIFIED. |

## Provider and transaction decisions

The controller owns the catalog/topology and calls `MakeEventPlan` itself.
An event provider supplies a caller-qualified `KinematicPolicy` at the actual
event, including every fresh product clock and explicit bath moments. It
cannot substitute an arbitrary event plan or silently change the selected
registered reaction. No kinetic policy means refusal, not zero recoil.

`PrepareParticle` runs only after Moon has assigned a finite particle epoch.
New sources require an explicit state resolver and one fresh clock draw.
Survivors retain their optical clock when a completed outer interval is
replaced by the caller's next physical interval. Queued daughters initialize
controller history from `ReadReactionState` before entering the core queue,
because that continuation may bypass the outer Moon wrapper.

Pure `ConstantHazardEventDelay` reads the compensated pending depth and the
existing combined-rate sum, with no channel draw. The event is represented
at the first binary64 residual second not earlier than exact H/Σλ. Zero total rate
returns infinity. Unrepresentable positive event time is an error; no
zero-time fudge is applied. The controller retains an immutable interval-start
clock and the mover's actual initial remaining seconds. A compensated
three-term sum locates the residual against the mover's subtractive remaining
coordinate; an exact endpoint maps to that remaining value byte-for-byte.
Accepted flight reconstructs TOTAL exposure from this origin using
`AdvanceConstantHazardClock`, also used by the original interval primitive.
It does not independently sum rounded adaptive dt values. Selected channel
and consumed clock are stored before ownership changes. Pure repeated
previews cannot consume optical depth or RNG.

At Process, the actual parent x/v/weights/state/time are validated, a complete
plan is made and the real backend commits detached children while retaining
the mover-owned parent. Every child is initialized/queued for the actual
remaining physical interval, then Delete makes the core delete the parent
exactly once. No callback performs an MPI collective. Insufficient complete-
packet capacity terminates with `exit(__LINE__,__FILE__,message)` explaining
capacity, why products cannot be discarded and how to increase capacity or
reduce model-particle resolution. It never redraws/rescales the event.

No optical-depth/event tolerance is permitted. An exactly located caller
endpoint consumes its original full physical duration; a deadline at the
next binary64 value beyond that duration remains pending. A native follow-on
interval then consumes the SAME tiny residual without another source draw.
Fresh daughters anchor to their actual queued unused seconds, rather than a
rounded subtraction from the old parent's duration. An unrepresentable
positive event residual is ERROR, never a forced zero-time reaction. Physical
epoch consistency checks retain the declared representation/time bound;
that bound is not an event-acceptance rule.

The first native campaign preserved a real FAIL: second-order adaptive
partitions left about1e-16s of optical lifetime after the mover reached its
declared endpoint. First order passed. This motivated the origin-based
reconstruction above; neither the exact endpoint criterion nor the future
deadline control was weakened. The original53PASS/8FAIL result is retained.

An additional analytical straight-line control places an exact1s reaction at
an internal AMR face AND the caller endpoint. It uses the prescribed zero
acceleration vector only for this geometry control,1m/s over1m, with no fit to
a lunar trajectory. Geometry retains deterministic priority, but the pending
event must publish exactly once on a legitimate owner at that endpoint, with
two children and zero remaining seconds. A child's exact-plane owner is not
prescribed. Physical surfaces and external-domain exits retain their declared
priority; this internal-face control does not qualify those distinct ties.

## Source-birth lifecycle audit and boundary controls

The real upstream allocator's `ResetParticleRecordForAllocation` clears every
`ParticleDataLength` byte after saving the free-list successor; both MPI and
hybrid `GetNewParticle` call it. Generic Exosphere injection also zeroes its
entire temporary record before cloning. The native test deletes/reuses the
same slot and verifies one genuinely new resolved state/clock. Commented
legacy main_lib surface-clone examples are not active capabilities.

The actual normal and FACE boundary-injection macros are separately invoked
with a new unprepared controller source and an explicitly supplied finite
birth epoch. They must route through the selected Moon wrapper, initialize
the clock once and complete daughters over only the supplied0.5s interval.
Each independent boundary case has its own accepted-budget check before its
fixture state is cleared. Its candidate diagnostics carry a `case` column;
the subsequent normal-walker cohort is separate. No caller silently assumes
that a linked legacy force-free boundary mover invokes Moon preparation.
Active generic particle splitting is refused until its stochastic-history
clone policy is qualified.
