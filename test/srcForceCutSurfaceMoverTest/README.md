# Force-aware physical-time cut-surface tracers

This target verifies the production first- and second-order tracers in
`src/pic/pic_mover_force_cut_surface.cpp`. It is analytical numerical/core
verification, not lunar observational validation or M3 qualification. The
contract is [M3_FORCE_CUT_CELL_CONTRACT.md](../../srcMoon/test/M3_FORCE_CUT_CELL_CONTRACT.md).
Existing AMPS defaults remain selected unless the application opts in.

The native target constructs a closed 3 m × 3 m × 0.1 m triangulated slab
whose top is z = −0.8 m, plus a genuine analytic internal sphere centred at
(1, 0, 0.5) m with radius 0.1 m. It uses the normal NASTRAN reader, AMR mesh,
default execution track, mover, boundary callbacks, temporary lists, MPI
owner exchange and continuation manager. Both orders use the same fixture.

Local base particle weights and time steps differ between neighbouring
leaves. The oracle preserves represented **N = W × correction**, and the
original caller interval in physical seconds across every crossing. This
is separate from steady-flux W/dt resampling. A source resets a fresh record
and assigns its physical birth epoch; MPI clones preserve continuation bytes.

Controls cover flight across owners, constant acceleration, a parabolic
floor entry whose interval endpoints are above the floor, gravity and
reflection, analytic-sphere and domain reflections, velocity-dependent
cyclotron/Coriolis form, an exact-time species event, and a detached sibling.
An accelerated sphere exercises the genuine quartic geometry. A moving
plane exercises the actual moving-surface finder/interaction callbacks and
relative-velocity reflection. An x=0 roundoff case guards owner oscillation.
The time-dependent grazing trajectory remains above the floor although a
coarse position arc predicts impact; adaptive execution must converge without
impact, and fixed execution must refuse an inconsistent leaving velocity.
The sibling is an additional independent queue control, not an atom-budget
reaction. Its epoch and remaining interval are checked independently.

`numerics_test.cpp` directly exercises the same production stepper/root
solver without AMPS. Independent oracles are closed-form acceleration,
quadratic intersections and exact unit-frequency rotation. The 64/128-step
rotation sequence checks first-/second-order convergence ratios; no present
production output is called physical truth.
There are 27 primitive checks, 124 native fixed-mode checks and 134 adaptive
checks. The native counts include BOTH integration orders. The arbitrary
fixture epoch1000s is a verification origin, not a UTC/SPICE observation.

An exact endpoint at an internal AMR face also exercises a due event before
the administrative owner handoff. Both orders check Continue with a species
change and Delete with a detached zero-duration child: the event commits once,
the physical epoch remains 1001 s, and no child receives another flight. An
internal mesh face has no material interaction. The zero-duration child may
remain in the current closed leaf until its next actual flight; it is neither
lost nor moved twice. This rule does not qualify an unresolved chemistry tie
with a material surface or external domain boundary; their existing geometric
priority is preserved. The test uses an exactly representable 0.125 s substep
and an exactly due event, without a time or optical-depth acceptance window.

Fast numerical controls, without application configuration:

```sh
g++ -std=c++17 -O2 -I src/pic test/srcForceCutSurfaceMoverTest/numerics_test.cpp -o /tmp/amps-force-cut-numerics-test
/tmp/amps-force-cut-numerics-test
```

Native controls, in an isolated checkout to preserve a configured Moon build:

```sh
rm -rf build
./Config.pl -application=test/force-cut-surface-mover -amps-test=off
make -j
mkdir -p PT/plots
mpiexec -n 1 ./amps
mpiexec -n 8 ./amps
mpiexec -n 8 ./amps adaptive
```

The committed thin runner runs the complete matrix, the expected-fatal
fixed-mode guard, the mathematical oracle and diagnostic comparison figures:

```sh
python3 -B test/srcForceCutSurfaceMoverTest/run.py --amps ./amps --output-root ./test_output/force-cut --mpi-ranks 8
```

It creates a new timestamped directory with `result.json`, `run_manifest.json`,
per-case stdout/report/CSV/resource receipts and the rotation-error comparison
in PNG/CSV/EPS. This host uses `/usr/bin/python3` with Matplotlib, plus g++ and
MPI. Missing optional Matplotlib is recorded as visualization SKIPPED while
the numerical kernel status remains separate. Compiler, runtime or missing/
truncated result failures are ERROR; evaluated numerical disagreements are
FAIL. A native exit1 is not sufficient evidence of a completed FAIL report.

The normal nightly registration is `unit_test_force_cut_surface_mover` in
`MakefileTest/Table`. In an isolated, installed AMPS checkout regenerate the
normal test makefile, then execute:

```sh
./Config.pl -install
make test_unit_test_force_cut_surface_mover
```

No Moon U/I identifier is reassigned. The runner's nonempty incomplete marker
prevents a killed job or missing report from leaving a stale empty test diff.

Results are `force-cut-surface.np<N>.report`; generated slab/mesh files are in
the run directory. Preserve stdout/stderr, executable/source hashes, input
and compiler/MPI identity with the report. An evaluated failed criterion
exits 1; build/fatal runtime failures are ERROR. Source implementation alone
does not claim these linked gates have run.

An application selects `MoverIntegratorMode=ForceCutSurfaceFirstOrder` or
`ForceCutSurfaceSecondOrder`, or calls `ForceCutSurface::RegisterStorage()`
before particle-buffer allocation and dispatches either function at runtime.
Configuration validation refuses DIM other than 3, hybrid/OpenMP, CUDA,
field-line lists, periodic wrapping, backward integration, any active
self-consistent PIC field solver/current-charge deposition and concurrent
local-time flux resampling. Prescribed Lorentz acceleration is supported;
self-consistent charge/current deposition is a different production path.

The source also compiles when MPI uses C++ worker threads. Its temporary-list
insertion selects the existing atomic-exchange protocol in that configuration
and retains the plain list head in serial MPI builds. Every worker must own a
distinct particle, and traversal/commit must wait for all insertion workers:
the head is published before its next link is completed. Initialize `prev`
before publishing so a later inserter's predecessor update is preserved.
The two threaded fast-wave regressions exercise build compatibility with this
source while selecting the ECSIM/Lapenta mover. They retain their physical
workload and numerical gates; they do not qualify threaded force/surface
callbacks or extend the native tracer matrix above.

Reproduce a native unsupported-mode refusal in an isolated checkout (replace
`periodic` with `backward` or `field-solver` for the other controls):

```sh
python3 -B test/srcForceCutSurfaceMoverTest/make_guard_input.py --mode periodic --output input/test/force-cut-refuse-periodic.input
rm -rf build
./Config.pl -application=test/force-cut-refuse-periodic -amps-test=off
make -j
mkdir -p refusal-run/PT/plots
tracer_amps_path="$PWD/amps"
(cd refusal-run; mpiexec -n 1 "$tracer_amps_path") > refusal-run/native.log 2>&1
tracer_refusal_code=$?
python3 -B test/srcForceCutSurfaceMoverTest/make_guard_input.py --mode periodic --check-log refusal-run/native.log --exit-code "$tracer_refusal_code"
```

These builds must exit through the declared production diagnostic before
geometry/particle execution. A PASS means the unsupported configuration was
refused. It never qualifies periodic tracing, backward clocks or a field
solver. The input generator refuses to overwrite an existing input file.

Public `ForceCutSurface::Settings` names the relative tolerance, position
[m]/velocity [m/s] absolute tolerances, geometry [m]/time [s] tolerances,
minimum/maximum substeps [s], adaptive mode and iteration/interaction limits.
Exhaustion terminates with an explanatory diagnostic, without discarding
unused physical time. The fixed 0.01 s native fixtures test accuracy
separately from adaptive acceptance. All constants are numerical fixture
definitions, not undocumented lunar physical parameters.

The force callback receives phase space, interpolation leaf and `cTimeContext`
with initial epoch and elapsed/remaining seconds. A NULL callback calls the
application's production total-acceleration macro; applications requiring an
evolving environment must install a force/geometry adapter that uses the
actual per-particle epoch rather than a shared step-start ephemeris.

The optional tracer's final opaque record occupies64bytes, aligned8. It
stores a stable global cut-face index through zero-time AMR/MPI continuations;
fresh reset and positive flight clear that history. Library, generated
particle layout and consumer must be rebuilt together. Default unregistered
movers retain their old record layout and behavior.

An event callback never frees/list-inserts the mover-owned current particle.
It returns Delete for core deletion, or updates the state, consumes its due
event and returns Continue while preserving remaining physical seconds.
New siblings are allocated unlinked, reset and assigned the event epoch;
`QueueDetachedParticle` validates/stages their TEMP continuation. It does
not move them recursively or insert them into the outer old-list sweep.
