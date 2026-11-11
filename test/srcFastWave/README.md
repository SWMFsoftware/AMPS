# Fast-wave native regression and performance

Run from an idle AMPS root:

```sh
make test_fast-wave -j
```

Keep other native campaigns and make invocations, including `make PDF`, out
of this interval: the native process guard conservatively rejects competing
make processes in the checkout before configuration/build work starts.

The maintained recipe is the `Name=test/fast-wave` block in
`MakefileTest/Table`; `Makefile.test` is its host-generated counterpart.
The target performs native preflight (root/idle check and removal of only
`build/`), selects the application through `Config.pl`, regenerates production
hooks, compiles through the top-level workflow with at most 16 build jobs,
creates `run_test_fast-wave`, runs the configured MPI launcher, and compares
the proton and electron outputs. The passion configuration uses four ranks.
The driver actually advances 60 iterations and produces three sampling files
per species. The timestep, weights, particles and sampling controls remain
defined by the existing driver and `input/test/fast-wave.input`.

## Implemented algorithm improvements

The subsequent P0–P4 implementation batches eligible periodic ghost records
per peer, builds complete uniform pointer stencils directly, reuses per-call
deposition geometry/live node addresses, accelerates the exact legacy CRC and
canonical particle copy, and scopes a compact scalar sparse view to each
linear solve. Arithmetic order, solver settings, repeatability checks,
physical workload and references remain unchanged. Custom/uncertain periodic
copies, incomplete/adaptive stencils and other backend paths retain their
documented general implementations.

The final exact native command passes in **79.387 s**, including **26.326 s**
simulation, compared with the already optimized 87.082 / 35.043 s baseline:
8.8% less complete-command time and 24.9% less simulation time. Every measured
stage reproduces all six full sampled files byte-for-byte. A second frozen
binary control takes 26.410 s. These are measured host-specific timings;
performance is not a numerical acceptance gate.

Detailed production contracts, per-stage evidence, the 18 permanent new test
cases, complete-core campaign status and limitations are in
[`STATUS.md`](STATUS.md). The implementation is documented beside the
[`periodic transport`](../../src/pic/parallel/README.md),
[`stencils`](../../src/pic/interpolation/README.md),
[`deposition`](../../src/pic/README.md),
[`CRC`](../../src/general/README.md) and
[`sparse solver`](../../srcInterface/README.md).

## Compiler fallback

The top-level Makefile defines `FAST_WAVE_CXXFLAGS`. For GCC host builds
with no explicit optimization level in the effective C/C++ flags, its default
is `-O1`. The recipe adds this through `FLAGC_EXTRA` only to its native AMPS
build. Debug symbols remain enabled when selected by the site. No site
configuration file, shared Fortran-library flags, other test recipe, input,
physics algorithm, comparison column or tolerance is changed.

Existing `-O0`, `-Og`, `-O1`, `-O2`, `-O3`, etc. choices are respected. There is
no automatic fallback for NVCC or another compiler family. To retain the
previous unoptimized policy, or use another explicitly selected extra flag:

```sh
make test_fast-wave -j FAST_WAVE_CXXFLAGS=
make test_fast-wave -j FAST_WAVE_CXXFLAGS=-O2
```

An explicit site optimization level can take precedence over extra flags;
these switches do not rewrite that site policy. The target prints its
additional flags and compiler invocations record the actual effective flags.
This change does not enable fast math, relaxed floating-point semantics,
additional MPI ranks, or reduced scientific work.

## Measured validation on passion, 2026-11-10

GNU C++ on an Intel Core i9-10980XE, four MPI ranks, clean native builds:

| Build policy | Preflight/configure/build/setup (s) | Simulation (s) | Comparison (s) | Total (s) |
|---|---:|---:|---:|---:|
| Original site flags, debug symbols and no optimization | 39.26 | 77.58 | 2.20 | 119.04 |
| Explicit `-O1` trial | 49.61 | 35.17 | 2.22 | 86.99 |
| Explicit `-O2` trial | 57.10 | 32.51 | 2.18 | 91.79 |
| Exact command with the new default | 49.84 | 35.04 | 2.20 | 87.08 |

The implemented default reduces complete command time by **26.8% (1.37x)**;
the simulation phase improves by **2.21x**. `-O2` runs slightly faster but
costs more compilation time for this mandatory clean-build workload.
The two `-O1` measurements agree within 0.1 s. These are host/load-dependent
observations, not a portable timing acceptance threshold. Phase boundaries
come from monotonic timestamps of stdout markers and are approximate.
Existing internal cumulative timers overlap; do not add them to obtain wall
time or interpret their inclusive labels as exclusive MPI costs.

All four completed native runs PASS: 60 iterations, four rank logs, empty
proton/electron/combined diff files under the unchanged `DiffNum.pl -a=2.0E-14`
contract. All **six** full sampled species files are byte-identical across
the baseline and candidates, including the column excluded by the historical
comparison. Twelve separate make probes PASS for fallback selection, explicit
site flags, compiler-family exclusions and user overrides. This validation
does not rerun the full core campaign or qualify other wave configurations.

Evidence directory:
`test_output/fast-wave-performance-20261110/`, including `summary.json`,
`optimization-policy-validation.json` and per-run `run.log`, `lines.jsonl`,
`time.txt`, `result.json`, sampled files, binary and diff files. An initial
sandbox attempt failed at MPI socket setup and is excluded from native
performance/physics results. The native runs then used normal local MPI.

The evidence driver copies references into the workspace, temporarily
redirects the legacy output symlink, and restores it after execution; original
reference hashes remain unchanged. The historical check recipe still edits
its reference files temporarily and can ignore child comparison failures.
For any future run, inspect both species diff files and the combined diff,
required output/rank logs and actual completion, alongside the make return
code. Missing files or a failed process are not agreement with the reference.

## Original individual-function profile

The subsequent four-rank function profile uses the same frozen validated
executable without an AMPS rebuild. A 100 Hz main-thread CPU sampler and PMPI
interposer collected 14,313 samples with no drops. The isolated control and
profile runs took 35.191 and 36.419 s; both passed the original comparisons
and reproduced all six sampled files byte-for-byte. Instrumentation overhead
was 3.49% in this measurement.

The leading measured targets are deposition, interpolation, sparse
matrix-vector multiplication, deterministic particle signatures, and periodic
ghost-block exchange. Timing-report collectives are low priority after actual
call-site measurement. That original profiling step changed no production
algorithm; the later implementation and validation are recorded above.

- [Measured findings and optimization strategy](PERFORMANCE_OPTIMIZATION_PLAN.md).
- [Profiler commands, machine requirements and interpretation](profiling/README.md).
- Raw measurements, JSON and plots:
  `test_output/fast-wave-function-profile-20261110/`.
