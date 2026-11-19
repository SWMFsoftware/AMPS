# Conservative reaction primitive controls

This supporting runner compiles **production** `species/ReactionExecution.cpp`
with the upstream `ChemicalSolver::cReactionNetwork` structural validator. It
then checks the same event planner/clock/weight/transaction APIs intended for
the PIC adapter. It does not configure/build AMPS, enable a reaction, implement
a daughter distribution, run a lunar trajectory or close a formal U/I gate.
No permanent test ID is added; these are supporting subcontracts of the
existing U15/U17/U24 and I21–I24 meanings.

Read [the pre-implementation contract](../M3_REACTION_EXECUTION_CONTRACT.md).
From the repository root:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_reaction_execution.py \
  --output-dir test_output/srcMoon/m3-reaction-execution-NEW
```

Use a new/empty output directory. `--cxx` selects a C++17 compiler executable.
The generated application `build/` and active machine configuration are not
modified. Exit codes are 0 PASS, 1 evaluated FAIL, 2 ERROR, 77 unavailable
compiler/SKIPPED. The parent/global M3 coordinator must invoke this actual
runner and retain its result; a zero launch exit or file presence is not a
replacement for the structured status.

The mechanical labels and convenient kg masses are **analytical controls,
not physical lunar species or observations**. Actual controls cover exact
integer nuclei/charge/state identities and overflow refusal; local physical
interval preservation and memcpy-safe interval/clock storage; open uniforms,
competing hazards, stable small probabilities, tiny accumulated exposure and
end-of-interval events; explicit bath electron multiplicities; independent
reference weights; repeated daughter carriers and fresh clocks; missing
kinematic/bath policy refusal; capacity deferral, staging/commit refusal,
backend exceptions, rollback, identical retry and double-commit refusal.

The fixed-seed 100000-draw channel control uses a predeclared conservative
99% simultaneous Hoeffding bound for probabilities 1/4 and 3/4. It is not a
cross-section model, an observation, or the complete U17 probabilistic-channel
qualification. Kinematics are supplied as an explicitly named analytical
mechanical policy, never inferred from a total rate or a scalar excess energy.

Artifacts: compiler/source/executable identities in `run_manifest.json`,
`compile.log`, `stdout.log`, `stderr.log`, named `result.json` outcomes and the
independent-count comparison `branch_control.csv`. No observational residuals
or model calibration is produced. A negative fixture succeeds only when the
production service refuses it without publication.

The backend fixture stages private state and publishes through allocation-free
swaps. A PIC backend must separately prove the same contract with actual
particle slots, lists, diagnostics, MPI ownership and mover continuations.
Rate-only plans cannot commit. Heap allocation exceptions before publication
are resource failures for the caller to handle, not permission to drop
products or redraw physical decisions. Catalog/runtime policy must remain
frozen while events are pending; plans retain resolved parent/product state
snapshots to prevent accidental identity reinterpretation.

The separate [native PIC backend runner](NATIVE_BACKEND.md) exercises the
actual allocator/node/cell-list adapter at1/8MPI. Its analytical nominal-mass
topology does not supply missing physical daughter distributions or promote
the local controls to formal chemistry qualification.

The [actual Moon force/tracer runner](MOON_FORCE_TRACERS.md) subsequently
uses that same matched generated library/ABI to check both normal physical-
time tracer orders against independently integrated production acceleration.
Its explicit analytical snapshot and shared-force RK4 reference verify
integration/wiring only; existing independent force oracles remain required.

The [combined constant-history controller](CONTROLLER.md) then connects the
same planner/backend to actual in-flight hooks. It requires explicit source/
hazard/kinematic providers, defaults disabled, and refuses trajectory-dependent
rates until their exposure service is qualified. Its analytical native runner
uses the same successful all-species build and never enables a physical lunar
channel by reading a rate table.
