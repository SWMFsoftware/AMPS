# Force-aware local-time cut-cell tracers: implementation contract

## Internal-AMR endpoint event amendment (2026-11-18)

The new explicitly bounded work package begins when this amendment is written;
its original source/results remain preserved. A due event at the caller's
exact endpoint must not disappear merely because that point is also an
administrative AMR face. After committing flight and projecting the exact
shared plane, a zero-horizon query may promote a due event ONLY when the
outgoing endpoint has an internal neighbour. It executes on the current
owner before the zero-time hand-off. Continue then performs the ordinary
handoff using the actual post-event species/velocity; Delete remains owned
by the mover, and detached products keep zero remaining seconds without an
extra full step. No event is moved earlier and no hazard/time epsilon is used.

Physical cut/sphere/moving-surface and external-domain priorities remain
unchanged. Their simultaneous endpoint chemistry conventions are NOT
qualified by this amendment. Tests extend the existing core fixture with
both orders/1 and8 ranks, exact binary steps, Continue and Delete with a
zero-duration detached child. The real Moon controller additionally verifies
complete conservative zero-duration products. No public API/record ABI changes.

Before any midpoint force evaluation, the proposed interval is clipped to
the earliest cut face, analytic sphere, moving surface, AMR face or reaction
event. Adaptive two-half-step estimates also preview their own second arc;
an earlier collision there causes refinement, never a force sample inside
the obstacle. The native fixture's force callback explicitly refuses states
inside either solid, so a collision PASS exercises this guard in production.
Integrated impact velocity must also enter the selected surface in its rest
frame. If midpoint velocity and the quadratic position extension disagree at
an accelerated grazing contact, adaptive mode refines BEFORE any physical
clock/state/ledger mutation; fixed mode terminates with a reduce-substep
diagnostic. The explicit time-dependent fixture is
`z(t)=-0.79-0.1*t+(50/6)*t^3` [m], interval 0.11 s, and stays above the
z=-0.8m floor; its coarse frozen-a arc predicts a spurious interior impact
at 0.1 s. Fixed-mode refusal and adaptive no-impact are separate controls.

The zero-force service uses an explicit physical-velocity overload. It
rejects entry inferred only from rounded displacement D=v*h when actual
stored velocity is leaving/tangent; wider arithmetic and the ordinary
callback dot product must both classify entry. No angular tolerance is
fitted. The original five-argument service uses NULL filtering and keeps its
original behavior. The archived adversarial vector has exact-double Decimal
v·N=+2.8699027721572765e-18 although ordinary double gives zero; this is local
numerical evidence, never an observational datum.

At AMR crossings, bounded root-rounding residuals in outward crossed
coordinates are projected onto the exact shared plane before owner
selection. Lookup probes use machine epsilon times domain/block scale, so
x=0 does not reduce to ineffective subnormal nextafter increments. Physical
time and arbitrary interior coordinates are unchanged.

The final opaque tracing record is 64 bytes, alignment 8 (earlier development
record 56 bytes). `lastFaceIndexPlusOne=0` means none. A stable global triangle
index survives zero-time local/MPI hand-offs; positive flight or fresh-history
reset clears it. A rank-local face pointer is never serialized. Producer,
generated layout/header, library and consumer must be built together; old/new
record layouts cannot be mixed. Geometry topology replacement during a live
trace is not qualified by this contract.

Written before implementation against HEAD
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`. The user explicitly requested
TWO particle-tracing movers, first and second order, based on the committed
`PIC::Mover::ForceFreeCutSurface`, with local caller-selected physical time,
all production forces and chemistry trajectory events. This supersedes the
proposal's optional common-step restriction for the mover itself. Chemistry
and evolving environmental clocks still require their own explicit policy.

## Inventory and production path

`pic_mover_force_free_cut_surface.cpp` traces a straight segment through every
crossed AMR leaf using conservative triangle lists and watertight tests. It
handles the earliest entering cut face, surface callbacks, domain boundaries
and temporary-list insertion. It does not integrate acceleration or invoke
trajectory-event callbacks, and a final unallocated leaf is a numerical loss.
`CutSurface::FindFirstCrossing` exposes its segment geometry to other movers.

The normal selected mover uses `_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_`;
Moon's configured callback already contains lunar gravity, Sun/Earth
differential gravity, rotating-frame terms, Na radiation pressure and Lorentz
force. The new core must call that SAME callback, including its velocity
argument; it must contain no Moon force formula or assumed E/B/plasma value.
An optional epoch-aware runtime acceleration callback is provided for explicit
local-time forcing and controlled analytical verification fields.

Existing Query/Advance/Process trajectory-event hooks are implemented only by
the first-order historical tracer. Existing local-time continuation fields
encode a W/dt flux-resampling invariant. New physical-seconds tracing must
neither reset a remaining interval to a destination dt nor silently resample
particle weights. New default-OFF storage/continuation is therefore required.

## Public API and switches

Two explicit compile-input selectors:

```text
MoverIntegratorMode=ForceCutSurfaceFirstOrder
MoverIntegratorMode=ForceCutSurfaceSecondOrder
```

They select corresponding `PIC::Mover` functions with `(ptr,dtSeconds,node)`.
All historical selection/defaults remain unchanged. A dedicated default-OFF
macro automatically reserves tracer state. Runtime `.in` applications may
instead call idempotent `RegisterStorage()` before particle-buffer initialization;
the physical-time continuation manager checks actual registration at runtime.
No registration means no extra record bytes or changed default ordering.
The new mode initially requires 3-D MPI particle-list tracing and
refuses CUDA, field-line mode, hybrid OpenMP and simultaneous legacy flux
resampling rather than claiming those unverified configurations.

`PIC::Mover::ForceCutSurface` exposes `cTimeContext` with `InitialEpochS`,
`ElapsedTimeS`, `RequestedTimeS`, `RemainingTimeS`; optional acceleration
callback `(a,species,ptr,x,v,node,const cTimeContext&)`; `SetParticleEpochS`,
`GetParticleEpochS`, `GetTimeContext`, and `ResetParticleState(byte*)`.
Acceleration NULL calls the existing configured force macro unchanged.
Absolute epoch can be unavailable (NaN) only for time-independent fixtures or
explicit frozen forcing. Evolving Moon adapters must supply actual per-particle
TDB epoch and use the context, not infer a global step time.

Fresh allocation already zeros the complete particle record. Reset means NEW
physical history: invoke after cloning a newly born physical particle, then
set its birth epoch. MPI/continuation clones preserve all bytes and MUST NOT
reset state. An inactive completed particle starts a new caller interval at
its accumulated current epoch; an active migrated particle resumes its exact
remaining seconds and original elapsed time. Requested dt is not truncated.

Named numerical settings include relative state tolerance, absolute position
[m] and velocity [m/s] tolerances, geometry/time tolerance, min/max substep [s],
maximum refinements/substeps/interactions and same-point event guard. Defaults
are numerical, not physical constants; all settings must be recordable and
adjustable by an application. Exhaustion is a detailed fatal ERROR with
settings/remaining time and remedies, never deletion or a false complete step.

## Integration and curved-path geometry

Both modes trace piecewise constant-acceleration position arcs
`x(t)=x0+v0*t+a0*t*t/2` over a numerically controlled segment. First order
freezes acceleration for the velocity update. Second order evaluates the
production acceleration at the predicted midpoint (including midpoint
velocity) for its velocity correction. For smooth position/velocity dependent
forces these have first/second global order respectively; uniform constant
acceleration is exact in both. Step-doubling supplies a state error estimate;
fixed-segment mode is explicitly available for analytical order tests.

Limit segments at the earliest physical crossing of the current leaf box.
Thus every tested arc lies inside a leaf whose conservative cut-face list is
known; no whole-step chord can bypass an intervening block/triangle. Solve
the quadratic arc/face-plane equation for ALL admissible roots, including
turn-back trajectories. Test the closed triangle with robust projected
orientation predicates; retain the earliest entering root. Edge/vertex ties
use a declared head-on/face-index order inherited from ForceFreeCutSurface.
Bound and report any rounding correction, never conceal a geometric leak.

The current moving-surface callback receives this same position arc and the
physical elapsed time. A Moon adapter must combine its explicit trace origin
with elapsed time: the old callback's shared `ProductionInjectionEpochS()`
alone is not an evolving local particle clock. Recompute acceleration after
reflection, species change, event or leaf crossing.

An event competes with cut surface, moving surface, domain and block-face
times. Geometry wins exact ties, matching existing event ownership. Advance
the event clock only by accepted physical time and before any MPI clone.
Process a due event once, sync phase space/species to the record, and continue
with actual remaining seconds. Callback deletion/list ownership stays in the
mover. New physical-time mode does not accept a hidden W/dt rescaling of its
remaining interval. Sibling/cascade chemistry is separately owned by the Moon
transaction service; the mover must not fabricate products.

At an ownership boundary, save remaining/elapsed seconds and stage in the
allocated destination ghost leaf. The normal exchange transfers ownership;
the selected continuation manager advances ONLY unfinished records and
exchanges after each pass, including just-completed remote records. Final
lists contain each survivor once. MPI migration is not physical deletion.

## Focused native verification specification

These are existing core-test extensions and supporting Moon mover coverage,
not new permanent U/I IDs and not observational validation.

| Field | Contract |
| --- | --- |
| Purpose | Verify both production tracers, force callback, geometry, physical-time ownership and event continuation. |
| Entry points | Actual new mover functions, settings/callback/state APIs, actual cut-face/moving-surface/domain callbacks and production list/exchange manager. |
| Inputs | Declared analytical constant-force, harmonic/linear velocity, Lorentz and Coriolis fields; triangulated plane/closed polyhedra; explicit local dt, seeds and individual weights. No observations. |
| Independent oracles | Exact constant-acceleration parabolic impact, exact uniform-B/Coriolis rotation, harmonic solution, independent high-precision root/ODE reference; list/number/weight identities. |
| Units / frame | m, s, m/s, m/s2; local Cartesian test frame; explicit optional absolute epoch and physical elapsed time. Moon force tests use its production LSO callback separately. |
| Procedure | Zero-force correspondence, curved-hit missed by a chord, turn-back/two roots, edge/vertex/grazing and start-on-face cases; reflection/absorption/domain behavior; dt refinement and error-control refusal; event before/after face and at exact tie; source birth/epoch; multi-leaf and 1/8-rank continuation with unequal local dt and unchanged physical weights. |
| Metrics | State/impact-time error and measured order, callback evaluation context, root/plane/barycentric residuals, event counts/position/species/remainder, list ownership, final physical clock and weighted inventory, numerical error counters. |
| Acceptance | Uniform-force trajectory/root within configured roundoff/geometry tolerances; deterministic fixtures conserve declared quantities; fixed-step convergence approaches factors2/4 for first/second order without relaxing gates after results; temporal endpoint invariant within declared floating accumulation bound; no missing/duplicated particles/events or numerical deletions. |
| Expected failures | Invalid settings/dt/state/epoch; nonfinite acceleration; incompatible execution mode; unavailable leaf/continuation storage; inward/unresolved callback; capacity/event loop exhaustion; incomplete event callback set. All are ERROR, never science PASS. |
| Data / hashes | No external observations. Record source/compiler/config/executable/fixture hashes. Production Moon-linked checks additionally hash actual furnished SPICE. |
| Controls | Deterministic force/intersection cases; OMP=1; explicit1/8 MPI; no bitwise cross-rank demand for stochastic surface laws. |
| Predecessors | Existing cut-surface/core baseline; local force kernels; normal M0; full same-source core preservation after shared edits. Formal Moon M1/M2/M3 predecessors remain unchanged. |
| Artifacts | Native stdout/stderr, per-case result/metrics/reference CSV, source/config/build identity, raw state/event/list audits and resource log; visualization when useful. |
| Status | Supporting local/native PASS only for evaluated criteria. Missing dependency SKIPPED; computed disagreement FAIL; parser/build/runtime/corruption ERROR. No qualification inferred from source inspection. |

Do not configure/build the active checkout while other work is using it.
Parent coordinates isolated normal configuration, targeted native tests, M0
and full core preservation. Remove its generated build before recompilation.
No commit, push, staging or reference replacement is authorized.
