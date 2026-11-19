# Combined in-flight constant-history controller

Read the [pre-implementation contract](../M3_REACTION_CONTROLLER_CONTRACT.md).
`MoonReactionController` connects the same conservative primitive/backend
services to the actual core Query/Advance/Process trajectory hooks. It is
DISABLED by default; there is no `.in` physics activation flag or automatic
attachment of D04/D05/CEX rate registries. Explicit programmatic providers
must declare immutable constant rates, source state, open uniforms and
qualified actual-event kinematics/bath moments. Missing data do not create
zero recoil or an isotropic angular distribution.

Trajectory-dependent/illumination-dependent rates are rejected. Freezing a
photorate at injection would not integrate exposure along a moving particle.
The controller's first scope is constant over each species/state physical
history, including outer-step, AMR/MPI and numerical-refinement continuation.
Query is side-effect-free; accepted Advance reconstructs compensated optical
depth from an immutable origin and the actual mover remaining-time coordinate,
and stores a single selected channel before migration. This prevents adaptive
substep summation from omitting a mathematically exact endpoint event. A
strictly future nextafter deadline remains pending, and continues into the
next interval without a source-clock redraw. Process constructs the
registered conservative plan at the actual event, queues every detached
child for remaining physical time and returns Delete for core-owned parent
deletion. Whole-step Moon photochemistry only reinserts while the controller
is active, avoiding a second lifetime/removal decision.

Initial states are selected explicitly even when they share a compiled PIC
species. Allocator reuse clears all opaque bytes; ordinary migration clones
retain history, while genuine new source clones must call
`ResetForNewParticle`. Generic Exosphere injection starts from a zeroed full
temporary record. Active particle splitting is rejected until a separately
qualified history/RNG cloning policy exists. Actual controller products read
their fresh clocks/state from the native backend; source clocks are not
redrawn for them.

Build the all-species [native backend](NATIVE_BACKEND.md) in an isolated
checkout, then before restoring its generated ABI:

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_reaction_controller.py \
  --producer-manifest /absolute/path/to/PASS-backend/run_manifest.json \
  --output-dir /absolute/path/to/new/controller-evidence \
  --wall-seconds 600
```

API: `run(output, producer_manifest, wall_seconds=600)`. The parent/global
coordinator invokes it after the actual Moon-force tracer runner, consuming
the SAME producer/source/header/library identities and remaining aggregate
budget. This thin runner performs no application configuration/build. A
separate explicit development test package must keep clone/setup/clean
builds/native work/cleanup within3600s; do not reset an earlier package's clock.

The numerical controls use three explicit Na source states, two competing
constant rates0.25/0.75s^-1 per state, predetermined open clock/channel draws
and a nominal-heavy-mass zero-Q Na+He+→Na++He topology. Products use explicit
analytical mechanical bath/velocity policies, not physical CEX data. Real
normal PIC walking and both Moon tracers run at1/8MPI. Checks cover pure
repeated previews, source-clock initialization/recycled slots, channel/source
identity, accepted event epochs/remaining intervals, pending-clock migration,
fresh product clocks/continued motion, unequal species dt/reference weights,
whole-step photolytic bypass and exact accepted material/number budgets.
Both actual normal/FACE source-boundary macros are also invoked without
direct controller preparation; their finite supplied birth epoch and0.5s
remaining interval must survive the normal wrapper and child continuation.
Separate `case` labels distinguish these independent cohorts in diagnostics.

The final fixture has18 named controls per MPI/order job (four jobs), plus
one expected capacity refusal:73 evaluations. An exact1m/1m-per-second
analytical zero-acceleration cohort verifies the internal-face/event/endpoint
tie and zero-time children without prescribing a particular exact-plane owner.
The other cohorts exercise actual Moon forces. A previous unchanged-core
run retained69PASS/4FAIL, isolating this missing endpoint dispatch; its results
remain separate from the corrected producer. Physical-surface/external-domain
coincidence follows its existing priority and is not qualified by this control.

A dedicated1MPI negative subprocess fills the actual buffer so the complete
packet cannot fit. It passes only on the controller's specific detailed
resource termination with the probe's documented exit1, including instructions to increase capacity or reduce
model-particle resolution. A generic launcher/runtime failure is ERROR. No
event is redrawn or partial packet accepted to make it fit. Signals and
unrelated launcher codes remain ERROR even if buffered stdout contains the
expected diagnostic. Protocol-only checks are executable with:

```bash
/usr/bin/python3 -B -m unittest discover -s srcMoon/test \\
  -p test_m3_reaction_controller.py -v
```

`result.json` and `run_manifest.json` retain producer/compiler/generated ABI,
source/test executable identities, actual commands and all logs. Each layout
has named check CSV and per-rank candidate-event/product-state CSV. Candidate
rows are work diagnostics until paired with accepted ledger entries by the
native budget checks; they are not independently called accepted events.
`event_comparison.csv`, `.png`, `.eps` compare actual and closed-form event
delays with the analytical-only label. Missing Matplotlib is SKIPPED before
execution; no figure is claimed. Exit codes0/1/2/77 mean PASS/FAIL/ERROR/SKIPPED.
Formal chemistry/predecessor gates and physical kinetic qualification remain
separate; no permanent test ID is introduced.
