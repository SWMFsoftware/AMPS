# srcSEP3D

`srcSEP3D` is the AMPS application being developed for three-dimensional
solar-energetic-particle and energetic-electron transport in the heliosphere.
The target model will support Parker and focused transport in either an
analytic Parker spiral or a coupled SWMF/AWSoM background, with scattering
controlled by prescribed or coupled Alfvén-turbulence state and shock source
parameters supplied by SWCME.

## Current state: Phases R0–R2

Phase R0 (**Rebaseline the Production Tree**), Phase R1 (**Complete Shared
Kernel and SWCME Integration**), and Phase R2 (**Implement Runtime and
Lifecycle**) are implemented in this source package. This is now a tested
orchestration foundation, but it is not yet a runnable transport model: mesh,
background-field population, and mover phases still stop explicitly at their
unimplemented boundaries.

The retained R0 baseline is:

- the obsolete `SEP3D.cpp` translation unit was deleted;
- the retired axisymmetric mover and legacy global sampler are absent from all
  production sources and archives;
- the unused `amps/amps_sampling.h` declaration-only scaffold was deleted;
- the unintegrated `core/sep3d_energy_distribution.{h,cpp}` sampler replacement
  was also removed; sampling will return as one tested Phase O implementation;
- `main_lib.cpp` no longer creates the wedge mesh, prepopulates a Maxwellian,
  writes placeholder fields, or advertises a partially working calculation;
- `core/sep3d_types.h` contains a semantic `ParticleMotionOutcome`, not copied
  AMPS integer values;
- `amps/amps_mover_status.h` is the only AMPS mover-status translation boundary
  and verifies the real `pic.h` values with `static_assert`;
- the production manifest excludes every retired prototype object; the current
  manifest contains `main_lib.cpp`, `main.cpp`, and the R2 runtime sources;
- makefile-owned external paths are canonical absolute paths derived from the
  active makefile location, so the same file works in `AMPS/srcSEP3D` and after
  AMPS copies the application to `AMPS/build/main`;
- a configured enclosing-production target and archive-symbol audit were added;
- the test CLI was rebuilt to follow the `srcSEP` runner vocabulary and to
  distinguish standalone, source/ABI, and configured-AMPS evidence.

R1 adds one canonical build of each shared dependency:

- `AMPS/src/models/sep_common/sep_common.a` owns the seven general SEP kernels;
- `AMPS/src/models/swcme/swcme.a` owns the compiled 3-D SWCME implementation;
- both application archives consume those canonical objects rather than
  compiling private copies;
- srcSEP3D depends only on the canonical `src/models/swcme` interface and does
  not inspect or require the independent `srcSEP` application tree;
- `src/models/swcme/test/run_tests.py` provides the restored srcSEP-style CLI
  and executable common 1-D/3-D integration gates.

R2 adds an AMPS-independent lifecycle implementation:

- `RunConfiguration3D::Create` validates host options and publishes an
  immutable configuration with a physics fingerprint;
- the complete per-cell/sampling storage layout is frozen before mesh binding;
- `Runtime` implements typed, transactional transitions from `Created` through
  `Finalized` and owns output, step, checkpoint, and restart counters;
- standalone Parker and SWMF adapters call the identical acquisition and
  snapshot-publication API;
- a minimal analytic Parker snapshot reaches `SnapshotReady` without AMPS;
- coupled production entry points do not parse process arguments,
  `AMPS_PARAM.in`, streams, or environment variables.

The detailed file-by-file record is in [MIGRATION_MANIFEST.md](MIGRATION_MANIFEST.md).

## Current limitations

R2 defines snapshot ownership and lifecycle metadata; it does not yet populate
AMPS cells with Parker or SWMF fields. The Parker/focused movers, heliospheric
AMR mesh, detailed SWMF field import, turbulence closure/evolution, SWCME
particle injection, sampling, checkpoint serialization, and validation
campaigns remain later phases. Calling an AMPS hook that needs one of those
components aborts with the specific missing phase instead of substituting a
placeholder calculation. A successful R0–R2 build therefore proves the shared
ownership and lifecycle contracts, not a scientific simulation.

## Source layout

```text
srcSEP3D/
├── core/                         AMPS-independent numerical/data types
│   ├── sep3d_types.h
│   └── sep3d_test_registry.h
├── background/                   AMPS-independent provider contract
│   └── bg_provider.h
├── runtime/                      immutable R2 configuration and lifecycle
│   ├── run_configuration.{h,cpp}
│   ├── runtime.{h,cpp}
│   └── runtime_adapters.{h,cpp}
├── amps/                         AMPS-only adapter boundary
│   └── amps_mover_status.h
├── SEP3D.h                       production umbrella header
├── main_lib.cpp                  process-owned Runtime and AMPS phase gates
├── main.cpp                      standard AMPS driver archive source
├── makefile                      production and standalone build gates
└── test/
    ├── run_tests.py              unified srcSEP-style runner
    ├── stage1.cpp                dependency-light C++ registry executable
    ├── individual-test/          HARN, LAY, BLD, UTIL, and LIFE3D callbacks
    └── frozen/                   reviewed byte-exact kernel reference
```

`SEP3D.cpp` and the prebuilt `test/stage1` binary are intentionally absent.
Generated objects, archives, reports, and executables must not be committed to
the source package.

When this package is overlaid on an older checkout, extraction cannot delete
files that existed only in the old tree. Remove the obsolete
`AMPS/srcSEP3D/SEP3D.cpp` explicitly before running `BLDL3D02`; retaining it is
a failed R0 installation even though the current makefile no longer links it.

## Layering contract

| Layer | Location | AMPS/MPI allowed? | Responsibility through R2 |
|---|---|---:|---|
| L0 | `core/` | No | semantic types and dependency-free kernels |
| L1 | `background/` | No | provider interfaces expressed in L0 types |
| L1/L2 | `runtime/` | No | immutable configuration, layout, lifecycle, and provider-neutral adapters |
| L2 | `amps/` | Yes | translate model outcomes to the AMPS ABI |
| L3 | `SEP3D.h`, `main_lib.cpp`, `main.cpp` | Yes | AMPS application boundary and process-owned Runtime |

Lower-layer code must not include `pic.h`, `mpi.h`, or use `PIC::`. The rule is
checked by source scan, an AMPS-free link, and a binary symbol-table test.

## Test runner

Run commands from the `srcSEP3D` directory. The CLI mirrors the selectors used
by `srcSEP/test/run_tests.py`:

```bash
# Discover IDs, groups, and named suites without building.
python3 test/run_tests.py --list

# Bounded normal-development gate.
python3 test/run_tests.py --routine --amps-source /path/to/AMPS

# One test or one/more groups.
python3 test/run_tests.py --test BLDL3D02
python3 test/run_tests.py --group HARN --group BLDL3D \
  --amps-source /path/to/AMPS

# Every public R0-R2 test, isolated so later tests continue after a failure.
python3 test/run_tests.py --all --amps-source /path/to/AMPS \
  --output-dir test_output/all

# Configured enclosing AMPS build and archive audit.
python3 test/run_tests.py --suite production \
  --amps-source /path/to/AMPS \
  --make-config /path/to/AMPS/Makefile.conf \
  --output-dir test_output/production
```

The runner always writes `srcsep3d-tests.json` and `srcsep3d-tests.xml` beneath
the selected output directory. At completion it prints PASS, FAIL, SKIP, and
ERROR totals and lists failed/error IDs. A missing `Makefile.conf` causes
`BLDL3D01` to **SKIP**; it never becomes a false production-build PASS.

The standalone executable uses the shared `sep_common.a` test registry. In a
normal AMPS tree the runner discovers the archive and headers automatically.
For a detached source package, specify them explicitly:

```bash
python3 test/run_tests.py --routine \
  --sep-common-dir /path/to/AMPS/src/models/sep_common \
  --sep-common-archive /path/to/sep_common/sep_common.a \
  --amps-source /path/to/AMPS
```

See [test/README.md](test/README.md) for test definitions, evidence classes,
exit codes, and troubleshooting.

## Production build gate

Within a configured AMPS application tree:

```bash
make strict-production
```

The target delegates to the enclosing AMPS `make amps` workflow and then uses
`nm` to reject retired symbols in `AMPS/build/main/mainlib.a` and `main.a`.
This is important because a direct compile from `AMPS/srcSEP3D` does not inherit
the generated-header include list assembled by the top-level build and cannot
find `AMPS/build/pic/pic.h`. The production gate therefore tests the real
copied `build/main` application rather than an artificial source-directory
submake.

Paths may be overridden without editing the makefile:

```bash
make strict-production \
  AMPS_ROOT=/path/to/AMPS \
  AMPS_CONFIG=/path/to/AMPS/Makefile.conf \
  SEP_COMMON_DIR=/path/to/AMPS/src/models/sep_common \
  SWCME_DIR=/path/to/AMPS/src/models/swcme
```

By default, the makefile first resolves its own absolute directory. It then
locates the AMPS root one level above `srcSEP3D` or two levels above
`build/main`. Consequently, the included `Makefile.conf`, `sep_common` and
SWCME source directories, archives, and object paths do not depend on the
process working directory. `make print-layout-paths` prints the resolved values
for diagnosis.

## Phase R0–R2 acceptance tests

| ID | Evidence | Acceptance |
|---|---|---|
| `BLDL3D01` | configured enclosing AMPS build | top-level `make amps` succeeds and `build/main` archives pass retired-symbol plus exact shared-member audits |
| `BLDL3D02` | production source/manifest scan | no retired source file, symbol, wedge constant, prepopulation, or placeholder output remains |
| `BLDL3D03` | AMPS `pic.h` plus adapter inspection | AMPS values are 0/2/3 for deleted/left/finished and all conversions occur in the L2 adapter |
| `BLDL3D04` | compile-only AMPS macro regression | `sep3d_types.h` remains valid when legacy `constants.h` has defined `Pi` |
| `BLDL3D05` | source/build makefile relocation fixture | source `srcSEP3D` and copied `build/main` makefiles resolve identical absolute paths when invoked from another directory, and production orchestration delegates to enclosing `make amps` |
| `LAY01`–`LAY02` | source scan and negative control | core/background/runtime remain independent of AMPS/MPI and the guard detects a deliberate violation |
| `BLD01` | `nm -u test/stage1` | standalone binary has no AMPS/MPI symbols |
| `HARN01`–`HARN04` | registry self-tests | selection, status, exit-code, JSON, and JUnit contracts work |
| `RUN3D01` | Python runner unit test | srcSEP-style selection, de-duplication, usage errors, JSON, and JUnit work |
| `UTIL02` | byte-exact record | shared kernel behavior matches the reviewed frozen record |
| `ARCH3D02` | canonical archive build/member audit | one exact `sep_common.a` and `swcme.a`; srcSEP3D consumes their objects without inspecting or requiring srcSEP |
| `SWCME3D01` | relocated SWCME runner | common configuration, 1-D, 3-D, and cross-dimensional tests pass |
| `LIFE3D01` | canonical transition path | lifecycle reaches `Finalized`; cadence and checkpoint counters are Runtime-owned |
| `LIFE3D02` | exhaustive transition matrix | all 80 state/operation pairs accept or reject as specified without partial mutation |
| `LIFE3D03` | layout and fingerprint checks | storage offsets are frozen before binding; physics/output identity is separated; mismatches are atomic |
| `LIFE3D04` | adapter and source audit | analytic and SWMF hosts reach `SnapshotReady` through the same API and coupled hooks contain no internal parser |

## Next work package

The next implementation phase is the mesh/storage layer. It must allocate the
exact `StorageLayout` frozen by R2, bind it transactionally, and connect the
production AMPS callbacks without changing configuration ownership. The
external configured-build gate `BLDL3D01` must also pass in the complete AMPS
checkout before R1/R2 can be merged; this source subset has no `Makefile.conf`
or generated `pic.h`, so it cannot manufacture that evidence locally. No
transport mover should be enabled before both the configured link and mesh
layout gates pass.
