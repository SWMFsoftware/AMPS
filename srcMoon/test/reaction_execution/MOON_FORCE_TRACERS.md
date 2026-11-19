# Actual Moon forces through both physical-time tracers

Read the [frozen supporting contract](../M3_MOON_FORCE_TRACER_CONTRACT.md).
This thin native probe uses the actual generated eight-species Moon library,
normal `Moon::MoveParticle` selector and `Moon::TotalParticleAcceleration`.
Core `ForceCutSurface::Acceleration` stays NULL, selecting the normal
production force macro. No test force is installed into either tracer.

First run the [native backend producer](NATIVE_BACKEND.md) in an isolated
checkout. Before changing/restoring that checkout's generated ABI, execute:

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_moon_force_tracers.py \
  --producer-manifest /absolute/path/to/PASS-backend/run_manifest.json \
  --output-dir /absolute/path/to/new/moon-force-evidence \
  --wall-seconds 600
```

Python API: `run(output, producer_manifest, wall_seconds=600)`. Launch the
script from the same retained checkout as the producer. It verifies that
producer's successful actual backend result, library bytes, current/preserved
generated headers and actual production source identities before compiling
only the new thin probe. The original producer manifest/source inventory is
retained. New test-only source is hashed separately; adding this probe does
not pretend to have rebuilt unchanged production objects. The compile uses
the actual MPI preprocessing flag `-DMPI_ON` as well as MPI linking.

The runner executes **order1/order2 ×1/8MPI**, with one neutralNa and oneNa+
global tracer per case, fixed physical horizon4s and maximum substeps
1,.5,.25s. It installs explicit stationary LSO E/B and a declared analytical
Sun/Earth/orbit/angular-velocity snapshot. These are dimensioned controls,
not an actual ephemeris, ARTEMIS observations, a calibrated solar-wind state
or a lunar observing case. No SPICE kernels or observations are consumed.
The path remains in its original leaf, so this fixture does not claim MPI
face-migration coverage; that requires the separate core mover tests.

Initial-state independent vector bookkeeping checks lunar gravity,
differential Sun/Earth gravity, centrifugal/Coriolis terms, q/m(E+v×B) and
Na radiation geometry. The Na scalar radiation curve remains the existing
production table: its independent physical qualification belongs to U05.
All applicable contributions must be nonzero. First/second integration
orders are then checked against independently coded long-double RK4 which
calls the **same actual production total force**. This qualifies integration
and wiring; it cannot independently validate that force law. U03–U06 and
formal production/convergence/data gates remain distinct prerequisites.

The predeclared contract gives separate position/velocity order windows and
requires the two fine RK4 references to differ by at most0.1% of the finest
tracer error. An insufficient oracle is ERROR; no tolerance is relaxed.
Actual species, represented number, opaque incoming clock and physical
interval/epoch are checked. Species local dt differs from the requested4s.
MPI layout agreement uses a16-binary64-epsilon rounding bound per state
component; it is not an observational agreement threshold.

Results are `result.json`, `run_manifest.json`, source/producer snapshots,
actual test executable identity and `np1-order1/`, `np1-order2/`,
`np8-order1/`, `np8-order2/`. Each retains stdout/stderr and named-check,
force-term and phase-space/error CSVs.
The `moon-force.reference.csv` file retains each4096/8192-step RK4
position/velocity difference and its fraction of the finest tracer error,
including an insufficient-reference result before termination.
`moon-force-convergence.png` and `.eps`
show deterministic numerical differences with the supporting-only scope in
the figure title. Exit codes0/1/2/77 mean PASS/FAIL/ERROR/SKIPPED. No permanent
test ID or formal physical/observational qualification is assigned.
