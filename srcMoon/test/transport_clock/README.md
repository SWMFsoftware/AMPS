# Common transport time: supporting M3 controls

The [contract](../M3_COMMON_TRANSPORT_CLOCK_CONTRACT.md) was written before
implementation. `MoonTransportClock` is PIC-independent numerical timekeeping;
normal Moon wiring is in `MoonTracing.cpp`, `MoonSurfaceEnvironment.cpp` and
`main_lib.cpp`. No new U/I ID, physical rate, daughter recoil law or observation
is introduced. Formal chemistry and earlier milestone gates remain separate.

## Runtime selection and semantics

```text
tracing_clock_mode = independent-history
```

This is the default when absent and preserves existing legacy/tracing receipts.
Explicit `common-time` requires a selected first/second force tracer, frozen
external LSO coefficients and `m2_environment=disabled`. Normal SPICE startup
interprets `simulation_start_utc` and supplies the explicit transport origin.
The common clock advances compensated **relative seconds**; frozen external
ET, geometry, source-CDF/photon factors and fields remain at that origin.
A source birth uses advancing transport begin plus its actual departure offset,
or transport end minus the supplied unused interval. There is no second RNG
draw and no repeated reuse of frozen ET as the next step's birth time.

Every actual owned block/species timestep must be identical, positive and
constant across outer steps. Both a SingleGlobal build and an actually equal
per-species/local table are admitted; the native probe deliberately uses the
latter to test real getters. Canonical `moon.input` is Na-only SpeciesGlobal,
so it already has one spatially global step. The check does not select a
physically accurate dt: gyro/flight, chemistry and convergence criteria remain
separate. Mixed asynchronous dt and unequal-duration iteration averaging are
refused before a step. Arbitrary asynchronous density/brightness sampling is
not qualified by this work.

Normal `amps_time_step` begins before PIC sampling/injection and commits after
all normal motion/MPI/product continuations. Common-only opaque packet targets
are copied with the real MPI/product records. Commit rejects pending histories,
stale targets and incomplete endpoints on every rank before accepting the step.
It may align only the declared two-addition binary64 projection envelope; that
never changes optical exposure or makes an event due. Actual allocator reuse
clears all runtime bytes (`pic_pbuffer.cpp`, `ResetParticleRecordForAllocation`
and MPI `GetNewParticle`); generic Exosphere injection zeroes its full runtime
record before cloning. The native test deletes a completed source carrier,
reuses the same slot in the NEXT physical step and checks its new birth/target.
Migration/product clones retain their target and are never reset as new sources.

The clock has no checkpoint protocol. Particle/sampling restart read, autosave
and SWMF restart are refused in common mode. Invalid source epochs, offsets,
dt tables or partial completion terminate normal runs with instructions; an
accepted physical step cannot be retried by rolling back only its clock.

Examples:

- [first order common-time](../../examples/m3_common_time_first_order.in)
- [second order common-time](../../examples/m3_common_time_second_order.in)

These are complete coarse **wiring** receipts. They establish no physical or
statistical convergence and activate no unsupported reaction provider.

## Local arithmetic and native parser controls

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_transport_clock.py \
  --local-only --output-dir test_output/srcMoon/my-clock-local \
  --wall-seconds 600
```

API: `run_local(output, wall_seconds=600)`. This configures no application and
compiles the maintained production clock plus parser/installer. There are37
named arithmetic controls and14 parser/default/refusal cases. The runner
requires every name exactly once and hashes sources again after execution;
an empty/truncated successful process cannot report PASS. Inputs are numerical
protocol fixtures, never observations. The100000-step0.1s control compares
with independent `n * actual-binary-dt` arithmetic; repeated absolute ET
addition is only a diagnostic counterexample, not the production clock.

Outputs include commands/compiler logs, named outcomes, `result.json` and
`run_manifest.json`. Missing compiler is SKIPPED; bad compilation/parser or
corrupted output is ERROR; evaluated invariant failure is FAIL.

## Actual PIC/source/mover hooks

After a clean successful all-species backend build in an isolated checkout,
while retaining its generated ABI and library:

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_transport_clock.py \
  --producer-manifest /absolute/PASS-backend/run_manifest.json \
  --output-dir /absolute/new-clock-native-evidence --wall-seconds 600
```

API: `run(output, producer_manifest, wall_seconds=600)`. The runner performs no
configuration/application build. It hashes the actual successful producer,
source, library and generated headers, then compiles a thin adapter and runs
both orders at1/8MPI. Each of four jobs has20 named global checks (80 control
evaluations). Equal actual dt, rejected unequal/changing dt, source boundaries,
recycled later sources, MPI continuation, bad completion refusal and normal
Begin/Commit chronology are exercised on real PIC lists. The native code uses
explicit analytical ET1000 and actual normal Moon force snapshots. It calls the
same production hooks/source macros/mover; it does **not** call full
`amps_init` or full `amps_time_step` and does not qualify UTC startup/sampler
ordering. That separate coverage is below.

Retained birth CSVs, `birth_comparison.csv/.png/.eps`, per-rank logs and manifests
compare actual and analytical source times. These figures contain numerical
controls only, never observed lunar quantities.

The final added composition control installs an explicitly analytical
controller registry before the buffer: the first three chronology steps use
an inert numerical Na state. A separate fourth-step tag selects ONE constant
1s^-1 channel with H=0.4. The declared nominal-heavy-mass zero-Q
Na+He+→Na++He control supplies co-moving incoming bath moments BEFORE product
construction; it is not a physical lunar CEX model. Actual conservative
publication creates one Na+ and one He carrier, fresh daughter clocks and the
cloned common target. Both must finish unused0.35s at common endpoint1004,
with one accepted exact nuclei/charge ledger and represented N=1 per child.
`daughter_comparison.csv` and per-rank daughter CSVs retain actual states;
missing/nonfinite/wrong-state CSV evidence is ERROR even after native PASS.
Floating weight reconstruction uses the existing declared backend rounding
bounds; integer topology and event counts remain exact. The original76-PASS
chronology-only evidence is retained separately from this80-control extension.

## Full normal startup and two production steps

The second thin driver calls the actual production parser/installer,
`amps_init()` and exactly two actual `amps_time_step()` calls. It runs each
order at8MPI in a restored canonical Na-only SPICE build. It changes no source,
force, sampler or reaction policy. A passive diagnostic chains/restores any
existing production event observer, so normal ledgers remain active. Checks
cover real CSPICE UTC interpretation, sampling-before-source order, actual
source callback coverage, completed common endpoints and unchanged frozen ET.

API:

```python
run_startup(output, wall_seconds=600, build_evidence=actual_restore_receipt)
```

Or, with the parent's REAL successful normal restore receipt saved as JSON:

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_transport_clock.py \
  --normal-startup-only --build-evidence /absolute/normal-build-evidence.json \
  --output-dir /absolute/new-normal-clock-evidence --wall-seconds 600
```

The mapping requires `status=PASS`, actual `source_root`,
`production_source_hashes`, `library` and `executable` path/SHA descriptors,
`build_log` path/SHA, retained successful `build_commands`, and complete
`generated_snapshots` (header path→SHA). Missing proof is SKIPPED; mismatched
bytes/commands or incomplete output is ERROR. A source-only preflight does not
create a producer receipt. The global development coordinator supplies this
proof after its actual successful normal restore; it must account for this
extra compile/startup work within the ORIGINAL3600s package budget.

Each order has16 named controls;32 evaluations total. Logs are
`compile-startup.log`, `{first,second}-order/native.log`; complete inputs,
actual furnished SPICE files/hashes, compiler/source/archive/header identities,
`normal_startup_chronology.csv/.png/.eps` and result/manifests are retained.
Actual archive/kernel identities establish consumed dependencies, not physical
or observational validation.

## Scope and reproducibility

Use one end-to-end development deadline, including clone/setup, clean builds,
retries, plots, reports and cleanup. The global coordinator reuses the same
fresh all-species producer for linked numerical probes, then restores the
normal Moon build before the full-startup control. Never configure the active
checkout solely to reuse a test's generated headers. Keep failed attempts and
source identities separate; no threshold is changed following failure.

Only an explicitly common physical horizon is addressed here. Evolving
external coefficients, asynchronous residence-time estimators, predictive
thermal/M2 coupling, unsupported source clocks, physical product laws and
complete restart/UQ/release qualification remain outside this contract.
