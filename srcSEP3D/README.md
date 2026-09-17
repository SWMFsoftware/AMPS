# srcSEP3D

`srcSEP3D` is the AMPS application for three-dimensional solar-energetic-
particle and energetic-electron transport in the heliosphere. Its target
physics is the Parker or focused transport equation in an analytic Parker or
coupled SWMF/AWSoM background, with Alfvén-turbulence scattering and SWCME
shock/source parameters.

## Implemented scope

The production tree implements the rebaseline and shared foundations (R0–R2),
**Phase M Mesh and Storage**, **Phase B Background Providers and Snapshots**,
**Phase T Turbulence and Scattering Inputs**, **Phase P Transport Cores**,
**Phase A AMPS Mover and Source Adapters**, and **Phase O Sampling, Output, and
Restart**.

`amps_time_step()` now enters the typed Runtime particle phase, calls the AMPS
step, and completes the Runtime cadence transition. Scientific production
still requires the configured AMPS mover macro and a host-installed local-state
resolver described in Phase A; no legacy mover, source, sampler, or fallback
physics is substituted when that coupling is absent.

### R0–R2 foundation

- The obsolete `SEP3D.cpp`, axisymmetric mover, legacy global sampler, wedge
  mesh, Maxwellian prepopulation, and placeholder output paths are absent.
- `src/models/sep_common/sep_common.a` is the single owner of general SEP
  transport/coefficient kernels used by both SEP applications.
- `src/models/swcme/swcme.a` is the single compiled SWCME implementation.
- srcSEP3D never inspects or requires the independent `srcSEP` application
  directory.
- `RunConfiguration3D` is immutable after validation and fingerprints every
  trajectory-relevant option. Output-only options remain outside that physics
  fingerprint.
- `Runtime` owns the lifecycle, active background generation, output cadence,
  checkpoint sequence, and restart counters.
- The standalone and SWMF adapters enter the same `Runtime` acquisition and
  publication path.
- AMPS mover integers are translated only in `amps/amps_mover_status.h`, where
  `static_assert` binds the adapter to the actual `pic.h` ABI.

### Phase M: mesh and storage

- Earth and Mars heliospheric domain presets are represented as Cartesian
  cubes containing a physical inner sphere and requested outer sphere.
- One resolution law provides radial refinement and an optional finite-width
  Parker-spiral tube with a continuous shoulder.
- The same AMPS-independent law drives the standalone octree verifier and the
  production `localResolution()` callback.
- The standalone octree enforces 2:1 face balance, deterministic global leaf
  and cell IDs, owner assignment, leaf histograms, and a pre-allocation memory
  estimate.
- The complete static/sampling byte layout is frozen in `RunConfiguration3D`
  before AMPS initializes its cell buffer. AMPS requests those exact bytes
  through its model allocation callbacks.
- Production cell population iterates only owner-local blocks. The physical
  cells inside the inner boundary remain allocated and zero-initialized; they
  are not presented as valid heliospheric background cells.
- Least-squares scalar/vector gradients support mixed coarse/fine neighbor
  distances and reject rank-deficient stencils.

See [MESH_STORAGE.md](MESH_STORAGE.md) for formulas, storage order, production
initialization order, and invariants.

### Phase B: backgrounds and immutable snapshots

- `BackgroundProvider` is active and independent of AMPS/MPI.
- `AnalyticParkerProvider` supplies Cartesian magnetic field, analytic
  magnetic/velocity gradients, focusing length, curvature, density,
  temperature, pressure, Alfvén speed, and validity metadata. Its Cartesian
  form has finite polar limits.
- `SwmfAwsomProvider` imports read-only SI or documented AWSoM coupling units,
  validates frame/epoch/ownership/completeness, and commits a generation only
  after the entire candidate succeeds.
- `BackgroundSnapshotBuilder` evaluates into temporary storage and publishes
  an immutable snapshot only after all points and fields validate.
- `SnapshotBuffer` supports current/next generations and linear interpolation
  only inside a compatible epoch bracket; extrapolation is rejected.
- Batch evaluation returns one status per point and never overwrites output for
  a failed point.
- The AMPS boundary stores every complete background field at its frozen
  offset, then publishes the same immutable object through either the
  standalone or SWMF `Runtime` adapter.

See [BACKGROUND_FIELD.md](BACKGROUND_FIELD.md) for units, field completeness,
coupling ownership, and atomicity rules.

### Phase T: turbulence and scattering inputs

- `TurbulenceProvider` is a separate authority from `BackgroundProvider`.
- The prescribed provider produces a normalized finite-band Kolmogorov
  spectrum with explicit amplitude, wave-number bounds, radial scaling, and
  correlation length.
- The AWSoM adapter accepts `w+` propagating along `+B` and `w-` propagating
  against `+B`, both in J/m³. It uses `δB²=μ₀w` and resolves outward/inward
  labels from the sign of `B·r`, so polarity reversals cannot silently swap
  coupling semantics.
- Out-of-band resonances use an explicit reject or power-law-extension policy;
  extension evaluates the physical resonant wave number and never clamps it to
  a band edge.
- Missing waves fail by default. Zero-rate ballistic transport exists only as
  an explicit configuration and is represented by a typed state.
- `CoefficientBridge` calls `sep_common` directly for Dμμ, mean-free-path,
  and parallel-diffusion conversions. srcSEP3D contains no copied coefficient
  formula.
- The provider API in `turbulence_models.h` is intentionally independent of
  `sep_common`. Coefficient users opt into `coefficient_bridge.h`; this keeps
  the AMPS-facing `main_lib.cpp` compilable under historic `Makefile.conf`
  rules that do not propagate application-local include flags.
- Self-consistent 3-D turbulence remains a hard reserved feature until a mesh
  wave-energy equation and conservation tests exist.

See [TURBULENCE_SCATTERING.md](TURBULENCE_SCATTERING.md) for conventions,
normalization, policies, and the shared-kernel boundary.

### Phase P: transport cores

- The Parker core advances the rank-one tensor
  `kappa_parallel * b * b` with the complete Itô drift, including the
  field-aligned coefficient gradient, field-line curvature, and `div(b)`.
- The focused core advances full gyrotropic focusing and flow coefficients
  with a symmetric split, reflecting pitch boundaries, and a declared
  Milstein or Euler–Maruyama stochastic scheme.
- Cell crossing, diffusion, focusing, cooling, background variation, shock
  crossing, and snapshot validity are separate named timestep limits.
- Counter-based random streams are keyed by campaign, particle, step,
  substep, and physical purpose, making histories independent of iteration
  order and worker ownership.
- Perpendicular diffusion and drifts must be exactly zero until their own
  physics and validation gates are implemented.

See [TRANSPORT_CORES.md](TRANSPORT_CORES.md) for the equations, splitting
algorithm, reproducibility contract, and Phase-P acceptance tests.

### Phase A: AMPS mover and SWCME source adapters

- One AMPS particle-buffer entry point validates and dispatches exactly the
  tensor Parker or split focused core selected by immutable configuration.
- A packed particle extension persists stable ID, stochastic step/substep,
  shock generation, momentum, pitch cosine, and gyrophase through migration
  and AMPS checkpointing.
- The adapter performs deterministic gyrotropic-to-Cartesian velocity
  reconstruction, exact destination-list insertion, and explicit terminal
  deletion/return-code mapping.
- Inner absorption, outer escape, invalid background, and failed transport are
  distinct semantic outcomes. Per-step/species integer ledgers require exact
  closure of active, injected, escaped, absorbed, and failed counts.
- Moving spherical shock crossings use the first analytic segment/surface root
  and are de-duplicated by shock generation.
- The SWCME adapter consumes the canonical common `SEPSourceState`, maps the
  DSA law to the shared `sep_common` injection sampler, and uses independent
  semantic random streams for momentum, pitch, gyrophase, and stable identity.

See [AMPS_ADAPTERS.md](AMPS_ADAPTERS.md) for buffer layout, dispatch, shock
geometry, DSA spectrum mapping, conservation, and host configuration.

### Phase O: sampling, output, and restart

- Read-only particle observations are sorted by stable ID and reduced with a
  specified compensated sum into cell moments, virtual-spacecraft spectra and
  anisotropy, field-line projections, and closed-ledger shock diagnostics.
- Output uses SI unit-bearing CSV schemas, per-artifact hashes, and a manifest
  containing configuration, code, and snapshot identities. A staging-directory
  rename publishes the complete sequence atomically.
- The independent parser verifies manifest keys, exact schemas, and hashes;
  corrupted or partial products never replace caller state.
- The restart codec writes a canonical versioned little-endian image rather
  than C++ object memory. It includes Runtime cadence/checkpoint counters,
  stochastic identity, all active particle state, snapshot/turbulence/source
  generations, sampling state, next stable ID, and closed ledger rows.
- Snapshot mismatch has an explicit reject or bounded-wait policy. Failed
  checkpoint writes roll Runtime back to `SnapshotReady` without incrementing
  the checkpoint sequence.

See [SAMPLING_OUTPUT_RESTART.md](SAMPLING_OUTPUT_RESTART.md) for algorithms,
file schemas, atomicity, restart contents, and lifecycle rules.

## Current limitations

The following are intentionally not enabled:

- perpendicular diffusion and gradient/curvature drifts;
- self-consistent 3-D turbulence evolution;
- external-script background providers;
- unconfigured direct access to mutable SWMF state from mover workers;
- full linked/MPI scientific validation and conservation campaigns.

For coupled operation, the host must configure SWMF authority and call
`InstallBackgroundSnapshot()` and, when selected,
`InstallTurbulenceProvider()` before `amps_init()`. A standalone analytic run
constructs its Parker snapshot and prescribed turbulence from the immutable
configuration. Before injecting particles, either host must install an
`AMPS::Movers::Context` whose resolver supplies coefficients from the pinned
snapshot. Output collection likewise must gather one globally stable-ID-ordered
observation set before invoking the Phase-O sampler.

## Source layout

```text
srcSEP3D/
├── core/                         semantic types and test registry
├── mesh/                         Phase-M resolution, octree, storage, gradients
├── background/                   Phase-B providers and immutable snapshots
├── turbulence/                   Phase-T providers, spectra, coefficient bridge
├── transport/                    Phase-P Parker/focused cores, timestep, RNG
├── adapters/                     Phase-A neutral dispatch, ledger, SWCME source
├── output/                       Phase-O sampling, publication, restart
├── runtime/                      immutable configuration and lifecycle
├── amps/                         AMPS-only ABI adapters
├── MESH_STORAGE.md
├── BACKGROUND_FIELD.md
├── TURBULENCE_SCATTERING.md
├── TRANSPORT_CORES.md
├── AMPS_ADAPTERS.md
├── SAMPLING_OUTPUT_RESTART.md
├── MIGRATION_MANIFEST.md
├── SEP3D.h                       production/coupling interface
├── main_lib.cpp                  AMPS mesh/storage/provider boundary
├── main.cpp                      standalone typed host and driver
├── makefile
└── test/
    ├── run_tests.py              srcSEP-style unified runner
    ├── stage1.cpp                AMPS/MPI-free C++ test executable
    ├── individual-test/          component acceptance callbacks
    └── frozen/                   reviewed byte-exact references
```

Generated objects, archives, binaries, reports, and `test_output/` are not
source deliverables. Overlaying this package on an older checkout cannot
delete stale files; remove an old `AMPS/srcSEP3D/SEP3D.cpp` explicitly if
`BLDL3D02` reports it.

## Layering contract

| Layer | Location | AMPS/MPI allowed? | Responsibility |
|---|---|---:|---|
| L0/L1 | `core/`, `mesh/` | No | types, resolution, standalone octree/storage, gradients |
| L1 | `background/` | No | analytic/imported ambient state and snapshots |
| L1 | `turbulence/` | No | scattering authority, spectra, AWSoM mapping, coefficient bridge |
| L1/L2 | `runtime/`, `transport/`, `adapters/`, `output/` | No | lifecycle, numerical transport, neutral coupling, diagnostics/restart |
| L2 | `amps/` | Yes | model-to-AMPS ABI translation |
| L3 | `SEP3D.h`, `main_lib.cpp`, `main.cpp` | Yes | AMPS allocation, owner-local filling, host entry points |

Every directory except `amps/` and L3 must not include AMPS/MPI headers or
refer to the AMPS namespace. `LAY01`, `LAY02`, and `BLD01` enforce this with a source
scan, negative control, AMPS-free link, and symbol-table inspection.

## Test runner

Run from `AMPS/srcSEP3D`:

```bash
# Discover all IDs, groups, and suites without building.
test/run_tests.py --list

# Fast development gate.
test/run_tests.py --routine --amps-source .. --rebuild

# One implemented phase.
test/run_tests.py --suite phase-m --rebuild
test/run_tests.py --suite phase-b --rebuild
test/run_tests.py --suite phase-t --rebuild
test/run_tests.py --suite phase-p --rebuild
test/run_tests.py --suite phase-a --rebuild
test/run_tests.py --suite phase-o --rebuild

# Complete source and configured-production evidence.
env MAKEFLAGS="-j16" test/run_tests.py --all \
  --amps-source .. --make-config ../Makefile.conf \
  --output-dir test_output/all --rebuild
```

`MAKEFLAGS` controls recursive GNU Make compilation, including the enclosing
AMPS build. Tests themselves remain intentionally sequential. Every run writes
JSON and JUnit summaries. A missing real `Makefile.conf` makes `BLDL3D01`
**SKIP**, never PASS.

For a detached application directory, provide the canonical shared dependency
explicitly:

```bash
test/run_tests.py --routine \
  --sep-common-dir /path/to/AMPS/src/models/sep_common \
  --sep-common-archive /path/to/AMPS/src/models/sep_common/sep_common.a \
  --amps-source /path/to/AMPS
```

See [test/README.md](test/README.md) for the complete evidence catalog,
selection rules, exit codes, and troubleshooting.

## Production build gate

Within a configured AMPS tree:

```bash
make strict-production
```

The application target delegates to the enclosing `make amps` workflow, then
audits `AMPS/build/main/mainlib.a` and `main.a`. This is required because AMPS
copies `srcSEP3D` to `build/main`; a direct source-directory compile lacks the
generated include/definition set and is not production evidence. All paths are
resolved from the active makefile, so source and copied locations find the
same `AMPS/Makefile.conf`, `src/models/sep_common`, and `src/models/swcme`.

## Implemented acceptance groups

| Group | Scope |
|---|---|
| `BLDL3D`, `ARCH3D`, `SWCME3D` | production routing, retired-symbol audit, canonical archives |
| `HARN`, `RUNNER`, `LAY`, `BLD`, `UTIL` | runner, layering, binary boundary, frozen common kernels |
| `LIFE3D01–04` | immutable configuration and complete lifecycle transition matrix |
| `MSH3D01–09` | resolution bounds/laws, tube geometry, balance, octrees, memory, ownership, presets, gradients |
| `BGP3D01–06` | analytic Parker identities, component laws, focusing, wind derivatives, polar limits |
| `SNAP3D01–08` | completeness, finite values, units, epochs, atomicity, interpolation, batch status, frame |
| `TUR3D01–04` | spectrum normalization, AWSoM mapping, resonance range, missing-data policy |
| `COEF3D01–02` | six-decade conversions and bitwise shared-kernel identity |
| `COEF3D03–05`, `PRK3D01–08` | tensor assembly/Itô drift and Parker transport behavior |
| `FTE3D01–07`, `RNG3D01–03` | focused transport, pitch boundaries, strong-scattering limit, keyed reproducibility |
| `ADP3D01`, `NAT3D04–05/08`, `SHK3D01–04` | mover dispatch, boundaries, ledger, moving shocks, common SWCME source |
| `NAT3D06–07`, `RST3D01–03` | sampling isolation, transactional output/schema, complete restart |

## Next release gate

The next work is configured-host integration and scientific validation: make
the generated AMPS mover macro see `AMPS::Movers::MoveParticle`, install the
analytic/SWMF local coefficient resolver, connect the coupled SWCME event
schedule and global observation gather, and run linked single-/multi-rank
conservation and restart campaigns. A configured `BLDL3D01` run on the target
AMPS checkout remains required after every production-boundary change.
