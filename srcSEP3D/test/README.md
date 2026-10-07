# srcSEP3D Testing Procedure

`test/run_tests.py` is the general application test interface. For combined
SEP/corona shared-model and native coverage, use `test/run_coupled_sep_corona.py`
as described at the end of this document. The general runner's selectors
match `srcSEP/test/run_tests.py` so the two applications can use the same
automation habits. The runner combines R0/R1/R2 foundation evidence with
Phase-M mesh, Phase-B background, Phase-T turbulence/coefficient, Phase-P
transport, Phase-A adapter, and Phase-O sampling/restart evidence without
requiring AMPS for dependency-free tests. Phase V adds controlled integration
and physics checks plus explicit linked/cross-model/observational evidence
gates. R01–R09 production-runtime improvements are native registry tests, not
an external checklist.

## Quick commands

```bash
python3 test/run_tests.py --list
python3 test/run_tests.py --routine --amps-source /path/to/AMPS
python3 test/run_tests.py --test LAY01
python3 test/run_tests.py --suite r1 --suite r2
python3 test/run_tests.py --suite improvements-c --rebuild
python3 test/run_tests.py --suite improvements-r --rebuild
python3 test/run_tests.py --suite improvements-v --rebuild
python3 test/run_tests.py --test V2D01 --sep1d-source /path/to/AMPS/srcSEP
python3 test/run_tests.py --suite phase-m --suite phase-b --suite phase-t
python3 test/run_tests.py --suite phase-p --suite phase-a --suite phase-o
python3 test/run_tests.py --suite phase-v --rebuild
python3 test/run_tests.py --suite phase-v --amps /path/to/amps \
  --validation-input /path/to/reviewed-sep3d.in \
  --validation-data /path/to/evidence \
  --validation-launch-prefix "mpiexec -n 8"
python3 test/run_tests.py --test SCCM3D01 --amps /path/to/amps \
  --validation-input /path/to/reviewed-sep3d.in \
  --validation-launch-prefix "mpiexec -n 4"
python3 test/run_tests.py --group HARN --group BLDL3D \
  --amps-source /path/to/AMPS
python3 test/run_tests.py --all --amps-source /path/to/AMPS \
  --output-dir test_output/all
```

The runner does not accept an implicit mode. Choose exactly one of `--list`,
`--test`/`--group`, `--routine`, `--all`, or `--suite`. This prevents an empty
or misspelled selection from exiting successfully.

### Reduced shock-front native gates

Run the entire selected reduced-profile campaign through one orchestrator:

```bash
cd /home/vtenishe/Mars2/AMPS
python3 srcSEP3D/test/run_reduced_shock_front.py
```

The default invocation uses the existing configured `./amps` and executes the
shared 54-assertion RSH harness, `ARCHCSWC01`, `RSHAPP01--03`, one- and
four-rank `RSH24--RSH27` smoke cases, an actual native uninterrupted/checkpoint/
resume matrix, and the separate four-rank `RSH24--RSH28` 1-AU campaign.  The
restart matrix compares four-rank uninterrupted with four-to-four resume, then
repeats the checkpoint at one rank and resumes it at four ranks.  It creates a
fresh directory below
`test_output/reduced-front/runner/` and writes `summary.txt`, `summary.json`,
one complete log per phase, native JSON receipts, products and artifacts.
Every FAIL/ERROR in the final summary carries its log path.  The expected
complete selected-profile aggregate is 93 PASS with no FAIL/SKIP/ERROR; IDs
repeated at different layers/rank counts remain distinct evidence records.
The 2026-10-05 native result is `93/0/0/0` under
`test_output/reduced-front/runner/20261005T-native-restart-qualification-02/`.
On the 2026-10-05 validation host, the non-restart native phases took
approximately 57 s (one-rank smoke), 36 s (four-rank smoke), and 1,750 s
(four-rank 1-AU).  The complete invocation therefore takes more than 30
minutes.  The runner emits a heartbeat
every 30 s by default and the growing `execution.log` records every committed
native tick; a quiet console during the long phase is not evidence of a hang.

Useful controls are:

```bash
python3 srcSEP3D/test/run_reduced_shock_front.py --list
python3 srcSEP3D/test/run_reduced_shock_front.py --dry-run
python3 srcSEP3D/test/run_reduced_shock_front.py --skip-native
python3 srcSEP3D/test/run_reduced_shock_front.py --rebuild-native
python3 srcSEP3D/test/run_reduced_shock_front.py \
  --mpi-launch-prefix 'srun -n {ranks}' \
  --output-dir test_output/reduced-front/runner/site-allocation-001
```

`--rebuild-native` is deliberately opt-in.  It requires execution from the
exact Mars2 AMPS root, rejects an active make/compiler/MPI/AMPS process and a
symlinked build target, removes only root `build` with `rm -rf -- build`, then
regenerates the srcSEP3D configuration/hooks and builds `amps` with `-j16`.
An existing output directory is rejected so older evidence cannot be
overwritten or mistaken for a fresh result.  `--skip-native` records every
omitted native ID as SKIP; it cannot produce the 93-PASS complete result.

#### Actual AMPS checkpoint/resume protocol

The restart gate is not the dependency-light shared-provider serialization
test.  The runner generates one evidence-local `restart-smoke.in` whose
`[output] checkpoint_cadence_steps = 5` and whose `[restart] output_path =
restart.chk`.  Every leg uses those same input bytes.  It then performs:

1. an uninterrupted four-rank run from tick 0 through tick 10;
2. a four-rank run through tick 5, which writes an actual AMPS checkpoint;
3. a four-rank `--restart <checkpoint>` run for five additional steps;
4. a one-rank checkpoint at tick 5; and
5. a four-rank resume of that one-rank checkpoint for five additional steps.

Each phase has its own working and output directory.  Relocation is necessary
because native publication refuses to overwrite an earlier process group's
products.  The checkpoint compatibility identity therefore freezes the
physics fingerprint, checksummed event/asset manifest, native storage layout,
output cadence and checkpoint cadence, while excluding only output,
initialization and checkpoint path names.  The complete provenance manifest
still records those path names.  A changed physical option, event asset,
storage layout or cadence is rejected transactionally.

At the common tick-10 boundary, `RSH29` and the Python comparator require exact
event identity, reduced-front generation/epoch/phase, every-record front-state
fingerprint, apex radius and normal speed, accepted/numerical area and shock
classification.  They also compare rank-independent XOR and modular-sum
fingerprints over actual owner-cell positions plus all ambient plasma/IMF
values after native readback.  Four-rank resumes must inspect nonzero received
physical ghosts, all owner/ghost/provider epoch checks must be true, and the
actual AMPS particle lists and source ledger must remain exactly zero.
Checkpoint source-rank count, input tick/generation, final checkpoint sequence
and final completed tick are recorded in each native JSON receipt.

The maintained loader selects deterministic stable-ID repartitioning.  The
one-to-four-rank case qualifies that mechanism for this background-only,
zero-particle profile.  It does not by itself qualify nonempty particle
repartitioning, which remains outside the reduced-model scope.

`RSHAPP01--03` are dependency-light application boundary tests.  Run their
binary from the `srcSEP3D` directory because other application fixtures use
paths relative to that documented working directory:

```bash
make -C srcSEP3D -j16 test/stage1
cd srcSEP3D
./test/stage1 --test-group RSHAPP
./test/stage1 --all-tests
```

These portable tests prove factory/adapter construction, immutable epochs,
rollback and deck/asset resolution from two working directories.  They cannot
prove AMPS owner storage or MPI synchronization.

The linked executable registers `RSH24--RSH29`.  `RSH24--RSH27` are members of
the `sep-corona` smoke suite; `RSH28` is long-only so the short handoff case
cannot accidentally pass or skip an arrival requirement, and `RSH29` is
explicit-only so a fresh smoke run cannot masquerade as restart evidence.
Native acceptance requires:

- `RSH24`: a rank-local injected candidate failure is rejected collectively
  without replacing the committed provider/front epoch;
- `RSH25`: actual owner fields and, on multiple ranks, received physical ghost
  fields and derivatives agree with the shared ambient epoch;
- `RSH26`: front/event identity, generation, apex state and absolute area
  ledgers agree on every rank and between the one-/four-rank smoke runs;
- `RSH27`: actual AMPS linked-list and source-ledger global counts remain zero;
- `RSH28`: a committed native horizon passes the exact continuous-time 1-AU
  observer root and reports geometric and accepted-shock arrival separately.
- `RSH29`: the process was actually launched with `--restart`, advanced beyond
  the serialized tick/generation, rebuilt owner and received-ghost state, and
  retained zero particles; the runner then compares it with the independent
  uninterrupted receipt at the identical final tick.

Use a unique `--output-dir` for every MPI process group.  JSON/artifact options
do not redirect observer products, and reusing the deck's default product path
correctly fails rather than overwriting evidence.  Exact smoke/long commands
and current result paths are in
`../../src/models/sep_corona_swcme/shock_front/README.md` and
`../../CODEX_REDUCED_SHOCK_PLAN.md`.

The long polar sensitivity fixture is expected to reach 1 AU geometrically
while returning `non-forward-inflow`; it must not be reported as accepted
shock arrival.  All reduced decks keep sources disabled and request zero
particle allocation.  None of these gates qualifies BG3D-4 downstream volume
physics.

### Parallel production compilation

The `BLDL3D01` gate delegates compilation to the enclosing AMPS GNU Make
build. Set `MAKEFLAGS` for this runner invocation to allow that build and its
recursive make operations to compile in parallel. The portable `env` form
works with `tcsh`, `csh`, `bash`, and `zsh`:

```bash
env MAKEFLAGS="-j16" test/run_tests.py --all --amps-source .. \
  --make-config ../Makefile.conf --output-dir test_output/all --rebuild
```

Adjust `16` to the CPU and memory available on the build host. Parallelism
applies to GNU Make compilation performed by the runner, including independent
standalone object compilation and the enclosing AMPS build. The tests
themselves remain sequential. The shorter `MAKEFLAGS="-j16" command`
assignment-prefix form is valid
in `bash` and `zsh`, but not in `csh` or `tcsh`; use the documented `env` form
for a shell-independent command.

## CLI correspondence with srcSEP

| srcSEP selector | srcSEP3D behavior |
|---|---|
| `--list` | prints stable IDs, groups, evidence kind, and suites |
| `--test ID` | runs one ID; repeatable and case-insensitive |
| `--group GROUP` | expands a registered group; repeatable |
| `--routine` | bounded normal-development set; excludes shell exit-code probes |
| `--all` | runs every public test separately and continues after failures |
| `--suite NAME` | runs one dependency/evidence class; repeatable |
| `--output-dir DIR` | writes the complete report bundle below `DIR` |
| `--amps PATH` | supplies the configured linked executable for Phase-V native cases |
| `--validation-input PATH` | supplies the immutable complete input deck used to construct linked native state |
| `--timeout SEC` | applies a per-command timeout |

Additional setup options are `--amps-source`, `--make-config`,
`--sep-common-dir`, `--sep-common-archive`, `--validation-data`,
`--validation-input`, and `--validation-launch-prefix`.

## Evidence classes

### Standalone C++ registry

`test/stage1` is compiled with no AMPS include path and no MPI library. It links
the Phase-M mesh model, Phase-B providers/snapshots, Phase-T turbulence bridge,
Phase-P transport cores, Phase-A neutral adapters, Phase-O output/restart,
Phase-V metrics/audits, R2 Runtime, test callbacks, and shared `sep_common.a`.
The runner adds
`-Wall -Wextra -Wpedantic -Werror`.

| Group | IDs | Purpose |
|---|---|---|
| `HARN` | `HARN01`–`HARN04` | registry selection, PASS/FAIL/SKIP/ERROR exits, JSON, and JUnit |
| `LAY` | `LAY01`, `LAY02` | AMPS-free core/background/runtime/mesh/turbulence/transport/adapters/output/validation rule and negative control |
| `BLD` | `BLD01` | `nm -u` confirms the standalone binary has no AMPS/MPI symbols |
| `UTIL` | `UTIL02` | byte-exact shared-kernel reference record |
| `LIFE3D` | `LIFE3D01`–`LIFE3D04` | immutable configuration, state machine, frozen layout, counters, adapter parity, and no-parser boundary |
| `R3D` | `R3D01`–`R3D09` | mover hook, requested-time loop, snapshot transaction, tick/events, source, observers, restart, canonical initialization source, finite empty-cell output |
| `CFG3D` | `CFG3D01`–`CFG3D13` | schema/CLI, typed contracts, domains, Parker geometry, mesh/memory preflight, finite-line/schema-3 initialization, complete compiled AMPS species binding, turbulence selection, CME/Parker linkage, schema-4 mover/coefficient and fixed/local source choices, active-corridor connectivity, corner geometry, and the shared application-section parser |
| `MSH3D` | `MSH3D01`–`MSH3D15` | resolution, exact Parker geometry, balance, octree budget/ownership, presets, gradients, finite line, initialization Tecplot output, conservative active-corridor classification, hole-free AMR topology, and fixed solar-boundary geometry |
| `BGP3D` | `BGP3D01`–`BGP3D06` | analytic Parker field/plasma identities and polar limits |
| `SNAP3D` | `SNAP3D01`–`SNAP3D08` | snapshot completeness, coupling conversion, atomicity, interpolation, batch/frame policy |
| `TUR3D` | `TUR3D01`–`TUR3D06` | spectrum, AWSoM convention, resonance, missing-data policy, selectable spectral/amplitude closures, and mandatory Tecplot energy |
| `COEF3D` | `COEF3D01`–`COEF3D07` | conversion/shared identity, tensor assembly, Itô drift, rejection, field-aligned gradient stencil, and selected-model/required-quantity isolation |
| `PRK3D` | `PRK3D01`–`PRK3D08` | Parker moments, characteristics, PDE/first passage, and named limits |
| `FTE3D` | `FTE3D01`–`FTE3D09` | focused streaming, focusing, pitch diffusion/boundaries, event-driven scattering, frame-energy invariants, momentum, and strong-scattering limit |
| `POP3D` | `POP3D01` | relativistic three-to-two weight, momentum, total-energy, and centroid conservation |
| `RNG3D` | `RNG3D01`–`RNG3D03` | worker/order independence and random-purpose isolation |
| `V1D` | `V1D01`–`V1D05` | V01 tensor, cross-field moments, drift, focused invariance, and timestep |
| `V2D` | `V2D01` | distinct compiled srcSEP/srcSEP3D production-core parity |
| `V5D` | `V5D01` | native profiles, deferred-R8 campaign block, and release governance |
| `ADP3D` | `ADP3D01` | exact three-core production registry and one validating dispatch |
| `NAT3D` | `NAT3D04`–`NAT3D08` | boundary outcomes, ledger closure, sampling isolation, output schema, shock crossing |
| `SHK3D` | `SHK3D01`–`SHK3D04` | common source identity, moving-sphere geometry, guards, and normalization |
| `RST3D` | `RST3D01`–`RST3D03` | full round trip, transactional rejection, and snapshot policy |
| `INT3D` | `INT3D01`–`INT3D03` | stable-ID rank merge, global conservation, and resource budgets |
| `VFY3D` | `VFY3D01`–`VFY3D06` | comparison metrics and analytical Parker/focused/SWCME/fixed-source validation |
| `RUNNER` | `RUN3D01` | Python selector, de-duplication, usage-error, JSON, and JUnit contract |

`HARN02-EXITCODE` and `HARN03-EXITCODE` are shell-level probes. They launch the
internal fail/skip beacons and verify the outer process status. The beacons are
not public `--all` tests because one is intentionally failing.

### Phase R0 source and ABI gates

#### BLDL3D01 — configured enclosing AMPS build

The runner locates or accepts a configured `Makefile.conf`, then executes:

```bash
make -f makefile strict-production \
  AMPS_ROOT=/path/to/AMPS \
  AMPS_CONFIG=/path/to/AMPS/Makefile.conf
```

That target invokes the top-level `make amps`, which owns generated headers,
include paths, compile definitions, libraries, and the final executable link.
It then audits `AMPS/build/main/mainlib.a` and `main.a` with `nm`/`ar`, including
exactly one member for each shared kernel, SWCME, and `sep_coronal_cme`
implementation object. A direct
`make lib` from `AMPS/srcSEP3D` is not production evidence because it does not
inherit the enclosing include configuration and cannot reliably locate
`build/pic/pic.h`. If the real configuration is absent, the result is SKIP. A
source-only compile or mock header is never reported as a production pass.

#### BLDL3D02 — retired-source exclusion

This check reads the live source tree and active makefile lines. It requires:

- no `SEP3D.cpp` file or makefile object;
- no retired mover, acceleration, sampler, or field-line injection symbol in
  production headers/sources;
- no former wedge bounds;
- no `PrepopulateDomain`, placeholder mesh output, or mesh-file write in
  `main_lib.cpp`;
- only the current R0–R2/M/B/T/P/A/O production manifest described in the root README.

The check deliberately scans production code, not documentation, because the
migration record must be allowed to name what was removed.

#### BLDL3D03 — AMPS mover status mapping

This check verifies both sides of the boundary:

1. `amps/amps_mover_status.h` contains the required `static_assert` and mapping
   branches.
2. The supplied AMPS `src/pic/pic.h` defines deleted-on-face `0`, left-domain
   `2`, and motion-finished `3`.

Use `--amps-source /path/to/AMPS` when the source is not in a discoverable
parent directory. Missing `pic.h` produces SKIP, not PASS.

#### BLDL3D04 — AMPS macro hygiene

AMPS `general/constants.h` defines `Pi` as a preprocessor macro. Namespaces do
not prevent macro substitution, so the former `Const::Pi` declaration broke
unrelated translation units after `pic.h` included the application header.
This compile-only regression defines the same macro before including
`core/sep3d_types.h` and requires the collision-safe `Const::kPi` API.

#### BLDL3D05 — copied-build makefile paths

AMPS copies the application into `AMPS/build/main` before compiling it. This
test constructs both `AMPS/srcSEP3D/makefile` and copied
`AMPS/build/main/makefile` fixtures, then invokes each with `make -f` from an
unrelated working directory. Both must resolve the same absolute `AMPS_ROOT`,
`AMPS_CONFIG`, `SEP_COMMON_DIR`, and `SWCME_DIR`; the active makefile directory
itself must match its source or copied location. This directly protects against the
`build/main: ../Makefile.conf: No such file or directory` failure.
The same fixture also supplies a synthetic enclosing `amps` target and verifies
that `strict-production` delegates to it and audits `build/main` archives rather
than attempting a bare compile in `srcSEP3D`. The synthetic build uses valid
ELF fixture objects, creates the complete application/shared member manifest,
and exports the normalized-domain and C05 mesh-preflight ABI names required by
the production audit. It therefore exercises the strict archive contract while
leaving real AMPS compilation coverage to BLDL3D01.

The fixture manifest is intentionally explicit and must be updated with the
production member list when a shared kernel is added. Stage 14 requires
`common_sep_coherent_transport.o` in this synthetic `mainlib.a`, matching the
coronal model archive that the application flattens into its production
archive. A missing-member failure under `BLDL3D05-layout/build/main` diagnoses
this fixture, not the real enclosing AMPS build. Retain the production audit:
it must still reject an archive with a missing or duplicated required member.

#### OUT3D01–02 and BLDL3D12 — excluded-volume output and header reentrancy

`OUT3D01` compiles the real `InterpolateInitializationCellData` body extracted
from `main_lib.cpp`, together with `output/sampling.cpp`. It exercises legal
zero-donor stencils, null unused arrays, optional gradient storage, unaligned
application offsets, canaries on neighbouring AMPS bytes, and malformed-input
rejection. `OUT3D02` additionally compiles the real print callback and checks
finite in-shell excluded rows, solar-interior rows, populated background with
empty/occupied particle windows, and owner-send/root-receive callback branches.
The host node/channel services are portable test doubles; these cases do not
claim native MPI transport or cut-cell geometry qualification. They are
selected by `--all`, `--group OUT3D`, and `--suite phase-o`.

`BLDL3D12` reads the permanent `src/pic/pic.h` and
`src/pic/ecsim/domain_bc.h` guards and compiles the actual PIC model-header
tail and `cDomainBC` declaration through distinct source/build copies. It makes
no repair-script call. Mesh metadata and recursive dependency stubs keep the
test portable; it is not a full PIC/MPI compilation. Negative controls remove
each guard correction and require class/default-argument compilation failures.
The package hygiene check requires both fixed headers. This case runs in
`--all`, `--group BLDL3D`, and `--suite production`.
The check lives inside `run_tests.py`; it does not require a separate script in
the AMPS `tools/` directory. `--amps-source` selects the checkout to inspect.
Generated probes and individual compiler logs remain in
`<output-dir>/pic-header-guards-probe/` for diagnosing failures.

```bash
python3 srcSEP3D/test/run_tests.py --group OUT3D --test BLDL3D12 --no-build --output-dir test_output/output-boundary
```

Run that command from the AMPS root. After applying the header repair and
rebuilding AMPS, repeat the multi-rank initialization output that originally
aborted. The seven native SCCM initialization contracts keep their existing IDs;
the new portable output cases belong to the application runner.

#### DOM3D01–04 and CFG3D12 — corner cube, solar neighbourhood and corridor

`DOM3D01` encloses densely sampled bent Parker curves and their complete
cross-sections under automatic/all eight explicit corners, translated origins,
rotated axes and strong winding. It checks endpoint scaling, plot-count
independence, complete solar-sphere containment and invalid controls.
`DOM3D02` checks off-corridor sphere/box intersection, face tangency, sphere
disablement and independent photospheric coarsening identities for all three
profiles. `DOM3D03` builds a balanced corner octree, proves dense line and
full sphere-surface coverage in physical-core leaves, checks a connected union,
retains topological halos and verifies pruning and photospheric preflight.

`DOM3D04` verifies the x-y corner mode with the Sun on the z midplane. It checks
automatic and all four x/y selections, whole tilted/polar corridors and their
cross-sections, full solar-sphere clearance, cubic extents, and rejection of a
conflicting z corner direction. For the equatorial case, it also builds the
real balanced mesh and active plan and checks reflected-z leaf bounds and
core/halo/inactive classes. It retains the complete sphere and source-to-endpoint
corridor while pruning unused leaves.

These four portable tests compile the actual `domain_geometry`, Parker and
mesh kernels without AMPS, MPI or replacement shared-model headers. Their
content-addressed executable lives under the runner output directory. They
are included by `--all`, `--group DOM3D`, and `--suite phase-m`.
`CFG3D12` additionally exercises the real input parser and immutable factory,
endpoint/arc-length normalization, geometry fingerprints, negative inputs and
dry-run reporting in the complete application Stage-1 build. It additionally
parses the x-y variant, checks its symmetric z bounds and distinct identity,
and rejects a nonzero z corner direction. The uploaded
overlay lacks canonical legacy sep_common/SWCME dependencies, so this latter
gate and a real linked MPI initialization require the enclosing AMPS checkout.

From `srcSEP3D`, the portable-only command is:

```sh
env MAKEFLAGS="-j16" python3 test/run_tests.py --group DOM3D --no-build --output-dir test_output/domain-geometry
```

See [`docs/DOMAIN_GEOMETRY.md`](../docs/DOMAIN_GEOMETRY.md) and the complete
`examples/sep3d_analytic_parker_corner_sphere.in` input for controls and usage.

#### BLDL3D06 — production transitive-header boundary

Every generic AMPS translation unit reaches `SEP3D.h` through generated
`pic.h`. Those compiler commands do not necessarily contain the canonical
`sep_common` or SWCME include paths. This regression rejects concrete
source/restart/particle-adapter includes from the umbrella header and rejects a
concrete `source_runtime.h` include from the AMPS adapter declaration. The
corresponding types must remain forward declarations until an implementation
file that is compiled with model include paths.

The same gate compiles `turbulence_models.h` with only the srcSEP3D include
root and requires it to be independent of `sep_common`. It then compiles the opt-in
`coefficient_bridge.h` with the canonical `SEP_COMMON_DIR` include path.

Finally, the gate creates a mock installed `Makefile.conf` whose generic
`main_lib.o` recipe intentionally ignores `CPPFLAGS`, `CXXFLAGS`, and
`INCLUDE`. The command has the same fixed shape seen in production AMPS logs.
It must nevertheless compile a translation unit that includes both
`sep_injection_spectrum.h` and `swcme_sep_source.hpp`, proving that the
makefile's target-scoped `CPLUS_INCLUDE_PATH` reaches `mpicxx`. The test thus
reproduces both former `sep_coefficient_physics.h` and
`sep_injection_spectrum.h: No such file or directory` failures without a
configured AMPS build.

#### BLDL3D07 — application-object ABI freshness

Reproducible source packages normalize file timestamps. If such a package is
overlaid on an existing configured tree, timestamp-only dependency checks can
mistake an older `build/main/mesh/mesh_model.o` for a current object. The
symptom is a final-link failure for the new normalized-domain or C05 preflight
symbols even though other, unchanged mesh symbols resolve.

The production makefile therefore forces every srcSEP3D-owned object to be
recompiled whenever its `lib` or `amps` target runs. It recreates indexed
archives, verifies exactly one copy of every application/shared member, and
uses `nm -C` to require the current `MakeDomain(RunConfiguration3DOptions)` and
`BuildRefinementPreflight` definitions plus the fixed-sphere
`MakeSolarBoundary(RunConfiguration3DOptions)` definition before AMPS reaches
`mpif90`. BLDL3D07 protects that makefile contract without requiring a
configured host.

#### BLDL3D08 — initialized native-background output ordering

This source gate protects the boundary between srcSEP3D's application-owned
background cache and AMPS' separate DATAFILE cache. It requires the production
driver to zero native padding, copy the validated density/velocity/temperature/
pressure/B/E/gradient fields, exchange block halos, and mark the installation
complete before the final `OutputDistributedDataTECPLOT` call. The call must
enable data printing for each compiled species and assemble the rank-local
cut-cell fragments into the requested initialization filename. This gate
prevents reverting to the whole-brick writer, which does not clip the solar
surface. It also requires the
physical-cell selection to enforce both spherical radii, unconditional wave
storage, immediate center-node turbulence write/readback, the registered AMPS
`InterpolateCenterNode` hook for application-owned bytes, and mandatory
total/directional SI turbulence-energy columns. The gate catches both observed
regressions: custom `B_x_T/B_y_T/B_z_T` initialized while native `Bx/By/Bz`
remained zero, and native plasma/IMF initialized while turbulence read as zero
from AMPS' temporary Tecplot vertex node. `TUR3D06` numerically checks that a
positive directional variance survives the same weighted interpolation. A configured
`BLDL3D01` run remains the compile/link authority for the AMPS API itself.

#### BLDL3D09 — active-region and population-control wiring

This routine source gate verifies the lifecycle facts that portable geometry
and conservation tests cannot establish alone. It requires the Parker-corridor
whole-mesh planner to consume AMPS' coarse/fine face/edge/corner links and feed
`SetTreeNodeActiveUseFlag` after `buildMesh()` and before load measurement,
distribution, and block allocation. It also requires an allocation audit after
`AllocateTreeBlocks()` so inactive leaves cannot retain storage. Inactive shock
patches must be excluded, the legacy automatic AMPS splitter must remain
disabled, and the SEP-aware relativistic controller must run after shock
injection but before observers and checkpoints. Finally, it checks the generic
controller's semantic ordinal is scoped to one cell/species population so MPI
block repartitioning cannot change post-resampling histories. It also checks
the generic AMPS split/merge entry points for the repaired non-positive-target,
empty-list, no-op, and singleton guards. `BLDL3D01` remains the configured
compile/link authority.

The guard audit is order-sensitive: the target check must precede linked-list
traversal/vector reservation, the population-size checks must precede the
first-record dereference, and an empty velocity-bin list must be handled before
requesting an iterator successor. A guard-looking token in a comment or after
the unsafe operation does not pass this gate.

The guarded core implementation is a required source-release member at
`src/pic/pic_particle_spliting.cpp`, not merely a prerequisite assumed to be
present in the destination checkout. Run the AMPS-level package audit before
publishing an overlay; it verifies that the exact file checked here is carried
with the SEP applications and prevents an older unguarded core implementation
from surviving installation.

The package policy therefore distinguishes **allowed** AMPS-level files from
**required** AMPS-level files. `TOP_LEVEL_ALLOWED` controls what may be copied;
`TOP_LEVEL_REQUIRED` additionally makes absence of the guarded splitter a
`missing-required` audit error. `BLDL3D09` parses that literal required-member
declaration rather than passing on a filename that appears only in a comment or
permissive allowlist. The hygiene self-test removes the fixture splitter and
requires the audit to fail, providing a negative control for this exact
packaging regression.

Because the policy is stored at the AMPS root, install a source archive from
the AMPS root and retain its `tools/` and `src/` members. Copying only
`srcSEP3D/` leaves both the old policy and potentially the old core splitter in
place, which is intentionally reported by this gate.

#### BLDL3D10 — solar internal-boundary wiring

This routine source gate covers the AMPS-only part of the photospheric
boundary that the dependency-light C++ registry cannot link. It requires
internal-boundary support and user-defined spherical callbacks at compile
time; exactly one all-rank
`Sphere::Init()`/`RegisterInternalSphere()` sequence after
`PIC::Init_BeforeParser()` and before mesh construction; the application
resolution callback; null injection hooks; and an absorbing callback that
returns `_PARTICLE_DELETED_ON_THE_FACE_` without double-deleting the particle.
It rejects a redundant direct `mesh->RegisterInternalBoundary()` call because
`RegisterInternalSphere()` already performs that registration.

The gate also requires one fixed `Core::Const::R_sun` geometry authority that
does not read `innerRadiusM`, rank agreement on center/radius, union of fully
solid leaves with the active-use mask, and post-`InitCellMeasure()` correction
of physical and ghost cells proven wholly inside the sphere. It checks that the
physics identity is `sep3d-physics-v8` and explicitly contains the
`amps-absorbing-sphere-v1` contract plus the fixed SI radius. `MSH3D15`
provides the numerical geometry/negative-control tests, while `BLDL3D01`
remains the configured AMPS compile/link authority.

#### BLDL3D11 — coronal-CME native-test wiring

This fast source gate requires the executable test CLI, immutable test-deck
loading, normal `amps_init_mesh()`/`amps_init()`/`amps_time_step()` calls,
collective read-only state capture, all seven `SCCM3D` evaluators, public
`sep_coronal_cme` validation calls, runner-compatible JSON, and flattened
model objects in `mainlib.a`. It also guards the physically correct
`nucleonCount = 0` source-identity representation used for electrons.

`BLDL3D11` cannot prove the generated AMPS ABI or MPI runtime. `BLDL3D01`
proves compilation/linkage, and `SCCM3D01–07` prove the initialized numerical
state on the configured executable. The detailed division of responsibility
is in `../validation/CORONAL_CME_NATIVE_TESTS.md`.

### Phase R1 shared-library gates

`ARCH3D02` runs both canonical archives' `verify` targets, compares exact `ar`
membership, and checks that the srcSEP3D makefile consumes the canonical object
lists. It deliberately does not inspect or require the independent `srcSEP`
application tree. `SWCME3D01` invokes the relocated SWCME runner and requires
`R1CFG01`, `R1D01`, `R3D01`, and `R13D01` to pass.

```bash
../src/models/swcme/test/run_tests.py --routine \
  --output-dir test_output/swcme-r1 --rebuild
```

### Phase R2 lifecycle gates

| ID | Acceptance contract |
|---|---|
| `LIFE3D01` | canonical legal path reaches `Finalized`; step/output/checkpoint counters remain inside `Runtime` |
| `LIFE3D02` | all 80 operation/state pairs match the transition table; rejected calls preserve state, counters, and snapshot generation |
| `LIFE3D03` | pre-mesh layout is deterministic; output-only settings do not change physics identity; mismatched binding is atomic |
| `LIFE3D04` | standalone Parker and SWMF adapters reach `SnapshotReady` through `Runtime`; coupled sources contain no process/file/environment parser |

```bash
python3 test/run_tests.py --suite r2 --rebuild \
  --output-dir test_output/r2
```

### R01–R09 production-runtime gates

| ID | Acceptance contract |
|---|---|
| `R3D01` | generated `picGlobal.dfn` hook names the exact mover signature and strict production depends on hook installation |
| `R3D02` | one AMPS request is consumed by multiple accepted substeps with fresh local resolution and exact final time |
| `R3D03` | failed fill preserves the active generation; a collectively accepted staged pair commits atomically |
| `R3D04` | integer clock/event schedule agrees with host time and restores without cadence drift |
| `R3D05` | physical source number survives rounding/cap, policies are counted, and successive ticks use distinct identities |
| `R3D06` | moving observer geometry, acceptance/uncertainty metadata, and commit-only accumulator reset are enforced |
| `R3D07` | schema-2 restart round-trips identities/layout, clocks/events, providers, shock, RNG tuple, ledgers, and sampling state |
| `R3D08` | canonical schema-3 provider preflights the first valid source surface and allocates the exact per-species, per-step particle count deterministically |
| `R3D09` | native Tecplot presentation distinguishes invalid background, no completed window, a valid empty particle cell, and an occupied cell without emitting `NaN` |

```bash
python3 test/run_tests.py --suite improvements-r --rebuild \
  --output-dir test_output/improvements-r
```

These tests are AMPS-independent and always participate in `--all`.
`BLDL3D01/03/05` remain the configured-host evidence for the generated mover
ABI and copied `build/main` production layout. On a configured tree the build
sequence is:

```bash
make -C srcSEP3D prepare-production
env MAKEFLAGS="-j16" srcSEP3D/test/run_tests.py --all \
  --amps-source . --make-config Makefile.conf \
  --output-dir srcSEP3D/test_output/all --rebuild
```

### Configuration, preflight, finite-line, and species-binding gates

| ID | Acceptance contract |
|---|---|
| `CFG3D01` | complete versioned input parses; file and typed construction have one fingerprint; documented CLI, early-error, and dry-run contracts hold |
| `CFG3D02` | output-only changes preserve physics identity; physical changes alter it; invalid shock/source intent fails; SWMF uses the same factory |
| `CFG3D03` | solar, one-AU, Mars, and explicit radii normalize exactly; invalid observers fail; the source shell may equal but not undercut the photosphere; transport-boundary status respects crossing direction |
| `CFG3D04` | mesh centerline/tangent and analytic field use one Parker geometry; polarity reverses `B` without moving the tube |
| `CFG3D05` | composite profiles are monotone, tube width scales from its reference, all memory categories/levels report, and an impossible level cap fails |
| `CFG3D06` | schema version 2 requires the complete finite Parker line and rejects a source-inconsistent initial point |
| `CFG3D07` | a complete mixed ion/electron table binds, while count mismatch, non-contiguous indices, duplicate symbols, invalid mass, neutral charge, out-of-range observers, and missing fingerprint state fail closed |
| `CFG3D08` | complete schema-3 SWCME input resolves while a missing canonical field, inconsistent weight, or skipped-step injection fails closed |
| `CFG3D09` | spectral and amplitude turbulence models parse only with consistent slopes and exactly one active amplitude normalization; Python background remains reserved |
| `CFG3D10` | `cme-launch-point` resolves the canonical SWCME launch apex and rejects radius/direction mismatches; explicit mode remains independent |
| `CFG3D11` | schema-4 active corridor/observer connectivity, population hysteresis, fixed/local source spectra, mover/coefficient compatibility, and dry-run output |
| `CFG3D12` | corner-domain modes, endpoint normalization, photospheric refinement identity, and conflicting controls fail closed |
| `CFG3D13` | shared sep3d sections, recursive includes, comments, continuations, CLI defaulting, immutable pre-mesh commit, and provenance-rich negative cases |

```bash
python3 test/run_tests.py --suite improvements-c --rebuild \
  --output-dir test_output/improvements-c
```

### Phase M mesh/storage gates

| ID | Acceptance contract |
|---|---|
| `MSH3D01` | one million deterministic points stay within the configured cell-size bounds |
| `MSH3D02` | linear radial surface, midpoint, and transition values match closed forms |
| `MSH3D03` | generated polarity-independent Parker centreline has negligible tube distance |
| `MSH3D04` | fast tube-distance approximation converges above second order |
| `MSH3D05` | composite tube mesh is 2:1 balanced and a deliberately illegal jump is detected |
| `MSH3D06` | co-rotation preserves resolution |
| `MSH3D07` | five octrees reproduce exact leaf/memory counts and reject non-owner writes |
| `MSH3D08` | Earth and Mars domains enclose exact declared outer spheres |
| `MSH3D09` | mixed coarse/fine gradients are linear-exact and rank-deficient stencils fail |
| `MSH3D10` | exact equal-arc Parker line preserves configured count/arc length and origin-relative refinement |
| `MSH3D11` | initialization Parker line is deterministic unit-labeled Tecplot data |
| `MSH3D12` | conservative finite Parker capsule intersects complete leaves and preserves full-domain mode |
| `MSH3D13` | analytic Parker derivative is parallel to the field tangent and arc-length inversion round-trips |
| `MSH3D14` | whole-octree mask covers a dense finite line, applies exact halo layers, prunes exterior leaves, fills cavities, and remains face-connected |
| `MSH3D15` | the internal sphere uses the fixed physical `R_sun`, remains distinct from the source shell, classifies solid boxes conservatively, and receives the clamped surface resolution |

```bash
python3 test/run_tests.py --suite phase-m --rebuild \
  --output-dir test_output/phase-m
```

### Phase B background/snapshot gates

| IDs | Acceptance contract |
|---|---|
| `BGP3D01–03` | analytic Parker divergence, components, and tangency |
| `BGP3D04–06` | focusing derivative, radial-wind derivatives, and finite polar limits |
| `SNAP3D01–02` | missing/non-finite required fields reject the candidate |
| `SNAP3D03–05` | coupling units, epoch consistency, and atomic failure |
| `SNAP3D06` | compatible snapshots interpolate only inside their time bracket |
| `SNAP3D07–08` | per-sample batch status and explicit coordinate-frame policy |

```bash
python3 test/run_tests.py --suite phase-b --rebuild \
  --output-dir test_output/phase-b
```

### Phase T turbulence/coefficient gates

| ID | Acceptance contract |
|---|---|
| `TUR3D01` | log-space quadrature recovers prescribed magnetic variance |
| `TUR3D02` | AWSoM energies use `deltaB²=mu0*w` and directions follow field polarity |
| `TUR3D03` | proton/electron resonances below, inside, and above the band apply the selected policy exactly |
| `TUR3D04` | incomplete waves fail unless ballistic mode is explicit and typed |
| `TUR3D05` | named spectral slopes, cross helicity, and directional SI wave-energy conversion agree |
| `TUR3D06` | direct wave-energy normalization follows its declared radial power and supplies mandatory total/directional Tecplot columns |
| `COEF3D01` | Dmumu/mean-free-path/kappa conversions round-trip below 1e-12 over six decades |
| `COEF3D02` | srcSEP3D bridge and direct `sep_common` Jokipii calls are bitwise identical |
| `COEF3D06` | centered and both one-sided stencils recover an exact nonzero linear `dKappa_parallel/ds`; no usable neighbor fails closed |
| `COEF3D07` | selector dispatch isolates unused physics and reproduces analytic constant-MFP, species-charge-aware radial-rigidity MFP, and constant-Dmumu limits |

```bash
python3 test/run_tests.py --suite phase-t --rebuild \
  --output-dir test_output/phase-t
```

`phase-t` deliberately contains `COEF3D01–02` and the host-neutral
`COEF3D06–07` selector/stencil tests; the later tensor/drift
coefficient tests belong to Phase P even though they share the `COEF3D` group.

### Phase P transport gates

| IDs | Acceptance contract |
|---|---|
| `COEF3D03–05` | rank-one tensor assembly, complete numerical/analytic Itô divergence, invalid/reserved coefficient rejection |
| `PRK3D01–08` | diffusion moments, advection/rotation, nonuniform equilibrium, cooling, first passage, radial PDE, all named limits |
| `FTE3D01–07` | ballistic characteristic, focusing, eigenmodes, reflecting boundaries, momentum, Parker limit, zero-perpendicular guard |
| `FTE3D08–09` | carried optical-depth event sequence and energy conservation in plasma/Alfvén scattering frames |
| `RNG3D01–03` | worker partition, list order, and future-purpose changes cannot alter keyed histories |
| `POP3D01` | relativistic 3-to-2 resampling conserves weight, momentum, energy, and position centroid |

```bash
python3 test/run_tests.py --suite phase-p --rebuild \
  --output-dir test_output/phase-p
```

### Phase A mover/source adapter gates

| ID | Acceptance contract |
|---|---|
| `ADP3D01` | registry contains Parker, focused diffusion, and focused scattering; all pass one validator/dispatch boundary |
| `NAT3D04` | inner absorption, outer escape, and invalid background remain distinct |
| `NAT3D05` | integer particle ledger closes exactly and mismatch is transactional |
| `NAT3D08` | first expanding-sphere root is recorded once per generation |
| `SHK3D01` | identical common SWCME records produce bitwise-identical dimensional source samples |
| `SHK3D02–04` | analytic moving shock, ownership/generation guards, and event-weight sum |

```bash
python3 test/run_tests.py --suite phase-a --rebuild \
  --output-dir test_output/phase-a
```

The standalone gate proves the AMPS-independent conversion and conservation
logic. `BLDL3D01` remains required to compile `amps_particle_adapter.cpp`
against the real particle-buffer/list ABI and the generated mover macro.

### Phase O sampling/output/restart gates

| ID | Acceptance contract |
|---|---|
| `NAT3D06` | stable-ID order and repeated sampling are bitwise identical; input remains untouched |
| `NAT3D07` | atomic output bundle has SI headers, complete identity manifest, and verified hashes |
| `RST3D01` | all state round-trips and future keyed normal draws are identical |
| `RST3D02` | fingerprint/checksum errors leave destination state unchanged |
| `RST3D03` | unavailable snapshot is rejected or awaited only through explicit bounded policy |

```bash
python3 test/run_tests.py --suite phase-o --rebuild \
  --output-dir test_output/phase-o
```

### Phase V integration/scientific-validation gates

The Phase-V suite combines immediately executable prerequisites with external
release evidence; their statuses must be interpreted separately.

| IDs | Acceptance contract |
|---|---|
| `INT3D01–03` | rank/order-independent stable-ID gather, exact global integer conservation, and explicit load/wall/memory budgets |
| `VFY3D01–02` | normalization, log-space metrics, coverage rejection, and malformed-evidence classification |
| `VFY3D03` | 3-D Parker projections reproduce the independent one-dimensional Gaussian Green function |
| `VFY3D04` | focused transport converges at second order to the exact focusing characteristic |
| `VFY3D05` | sampled SWCME/DSA momentum CDF and total represented event weight agree with their declared laws |
| `VFY3D06` | fixed phase-space q=5 overrides a q=4 shock patch and the sampled ensemble agrees with the independent `dN/dp proportional to p^-3` CDF |
| `NAT3D01–03/09–12` | configured AMPS mesh, storage, gradients, balance, budgets, coupled cadence, and output grammar |
| `MPI3D01–02` | multi-rank sampling and restart continuation are decomposition independent |
| `SCCM3D01–07` | shared-model initialization ledger, all compiled species' weights/time steps and source identities, solar/active mesh contract, provider generations, finite AMPS products, and collective identity |
| `XM3D01–06` | checksum-owned cross-model profiles, exact source identity, convergence, and independent PDE evidence |
| `OV3D01–04` | reviewed event comparisons, with release-gate and diagnostic roles retained in reports |
| `VALRUN3D01` | runner lists all classes, preserves SKIP, verifies checksums, evaluates convergence bundles, validates the OV3D01 known/unresolved-parameter blueprint, and proves that the native matrix rejects missing callbacks while hashing its exact input/executable |

```bash
# Runs controlled prerequisites now and records unavailable external work as SKIP.
python3 test/run_tests.py --suite phase-v --rebuild \
  --output-dir test_output/phase-v

# Adds independently exported cross-model/observational evidence.
python3 test/run_tests.py --suite phase-v \
  --validation-data /path/to/evidence \
  --output-dir test_output/phase-v-evidence

# Adds native tests from a configured linked application.
python3 test/run_tests.py --suite phase-v --amps /path/to/amps \
  --validation-input /path/to/reviewed-sep3d.in \
  --validation-launch-prefix "mpiexec -n 8" \
  --output-dir test_output/phase-v-linked
```

The launch prefix is parsed into process arguments and is never passed to a
shell. A linked binary must advertise the requested ID through `--list-tests`
and receives the deck only through `--test-input`; a missing deck, ordinary
production driver, or stale binary is `ERROR`. Scientific evidence must use the templates under
`validation/templates/`, remain inside `EVIDENCE_ROOT/CASE_ID`, and match every
declared SHA-256. Missing evidence is `SKIP`; checksum/schema/provenance failure
is `ERROR`; a valid metric outside tolerance is `FAIL`.

`XM3D03`, `OV3D03`, and `OV3D04` remain diagnostic. V01 now supplies controlled
perpendicular diffusion and guiding-centre drift, but these cases still need
reviewed evidence before cross-field transport can be called scientifically
validated. See `INTEGRATION_SCIENTIFIC_VALIDATION.md` for the full
equations, algorithms, case roles, and evidence schemas.

## Named suites

| Suite | Contents |
|---|---|
| `standalone` | C++ registry, shell exit-code probes, `RUN3D01`, and `VALRUN3D01` |
| `r0` | R0 source/ABI/production gates plus RUN3D01, LAY01, and BLD01 |
| `r1` | canonical shared-archive audit, relocated SWCME suite, and frozen common kernels |
| `r2` | LIFE3D01–LIFE3D04 immutable configuration and lifecycle gates |
| `improvements-c` | CFG3D01–CFG3D13 production configuration, preflight, finite-line/schema-3 initialization, species/turbulence/source selection, CME linkage, mover/coefficient selection, active-corridor, corner-domain, and shared-section parser gates |
| `improvements-r` | R3D01–R3D09 production runtime integration gates |
| `improvements-v` | V1D01–05 controlled physics, V2D01 true parity, and V5D01 governance |
| `phase-m` | MSH3D01–MSH3D15 mesh/storage, finite active-corridor, hole-free topology, fixed photosphere, and initialization-output gates |
| `phase-b` | BGP3D01–06 and SNAP3D01–08 background/snapshot gates |
| `phase-t` | TUR3D01–06, COEF3D01–02, and COEF3D06–07 turbulence/coefficient gates |
| `phase-p` | COEF3D03–05, PRK3D01–08, FTE3D01–09, RNG3D01–03, and POP3D01 |
| `phase-a` | ADP3D01, NAT3D04–05/08, SHK3D01–04 |
| `phase-o` | NAT3D06–07 and RST3D01–03 |
| `phase-v` | INT3D/VFY3D prerequisites, external NAT3D/MPI3D/XM3D/OV3D cases, and VALRUN3D01 |
| `production` | BLDL3D01–10 |

Suites can be repeated. Overlapping IDs are de-duplicated in stable order.

## Shared utility discovery

In the normal AMPS layout, the runner discovers `src/models/sep_common` from
the AMPS root and builds `sep_common.a` when needed. For a detached tree:

```bash
python3 test/run_tests.py --routine \
  --sep-common-dir /path/to/AMPS/src/models/sep_common \
  --sep-common-archive /path/to/sep_common/sep_common.a \
  --amps-source /path/to/AMPS
```

The archive must be the shared build used by the applications. The runner does
not compile private copies of shared kernel sources.

## Reports and exit codes

Every execution writes:

- `srcsep3d-tests.json` — commands, messages, elapsed times, totals, and final
  exit code;
- `srcsep3d-tests.xml` — JUnit representation for CI;
- one native `<ID>.json` for each C++ registry test.

Exit codes are:

| Code | Meaning |
|---:|---|
| 0 | all selected tests passed or were explicitly skipped |
| 1 | at least one test executed and failed an acceptance criterion |
| 2 | setup/usage error or at least one test could not be evaluated |

SKIP is an evidence limitation. It is acceptable in source-only development
but cannot close a release gate that requires the skipped test.

## Adding a test

For a dependency-free C++ test:

1. add a callback file under `test/individual-test/`;
2. declare its registration function in `core/sep3d_test_registry.h`;
3. register it in `test/stage1.cpp`;
4. add its source to the makefile and runner build command;
5. add one `TestDefinition` to `test/run_tests.py`;
6. document setup, independent reference, tolerance, and failure meaning here;
7. run `--list`, the individual ID, its group, `--routine`, and `--all`.

For an AMPS-linked test, do not add a mock replacement to the standalone
binary. Register it in the production executable and make the Python runner
invoke that exact linked callback, following the srcSEP pattern.

## Troubleshooting

| Symptom | Interpretation and action |
|---|---|
| `cannot locate sep_common directory` | pass `--sep-common-dir` |
| `cannot locate sep_common.a` | build the shared archive, then pass `--sep-common-archive` |
| `BLDL3D01 SKIP` | run in a configured AMPS checkout or pass `--make-config` |
| `BLDL3D01` reports `pic.h: No such file or directory` from `srcSEP3D` | update the makefile; the gate must delegate to top-level `make amps`, not compile `main_lib.cpp` directly in the source directory |
| `BLDL3D02` reports stale `SEP3D.cpp` | remove the obsolete file from the installed `srcSEP3D`; overlay extraction does not delete files left by an older version |
| `BLDL3D03 SKIP` | pass `--amps-source` pointing to a tree containing `src/pic/pic.h` |
| `build/main` cannot find `../Makefile.conf` | run `make print-layout-paths`; `AMPS_CONFIG` must resolve to the absolute `AMPS/Makefile.conf` path |
| `main_lib.cpp` reports `sep_injection_spectrum.h: No such file or directory` | install the updated srcSEP3D makefile in the source tree and refresh the copied `build/main`; `BLDL3D06` verifies that the fixed generic recipe receives the target-scoped canonical model search path |
| final link reports undefined `Mesh::MakeDomain(RunConfiguration3DOptions)`, `BuildRefinementPreflight`, or `MakeSolarBoundary` | stale application `mesh_model.o`; install the updated makefile, run `make clean`, and rebuild. BLDL3D07 prevents recurrence |
| standalone compile failure | rerun with `--rebuild --verbose` |
| Phase-V linked case `SKIP` | supply `--amps`; also supply `--validation-input`, and use `--validation-launch-prefix` when MPI launch arguments are required |
| Phase-V linked case reports that `--test-input` is required | pass the reviewed complete deck with `--validation-input`; hidden callback defaults are intentionally forbidden |
| `V5D01` cannot open `release/generate_release_evidence.py` | restore the complete `srcSEP3D/release/` source set (`generate_release_evidence.py`, `profiles.json`, `capabilities.json`, `README.md`, and `checklist.md`). These files are required test/governance inputs declared by `SOURCE_MANIFEST.json`, not generated build products. |
| XM3D/OV3D case `SKIP` | supply `--validation-data` containing `CASE_ID/manifest.json` and its declared artifacts |
| Phase-V checksum `ERROR` | regenerate the SHA-256 only after reviewing the changed evidence; never edit a hash merely to silence the gate |
| linked executable does not advertise a case | rebuild the configured AMPS application from this source tree and verify its production test registry |
| unknown test/group | use `--list`; unknown selectors are usage errors |
| report missing after a C++ test | treat as ERROR; inspect verbose subprocess output |

`CFG3D06`–`CFG3D10`, `TUR3D05`–`TUR3D06`, `MSH3D10`–`MSH3D15`, and `R3D08`–`R3D09` are routine C++ entries
in the runner manifest.
`CFG3D06` uses live negative controls for an omitted version-2 key and an
initial point inconsistent with the inner sphere. `MSH3D10` constructs the
configured number of vertices, sums every segment to the requested arc length,
and translates the origin/probe together to prove the AMR law is not tied to
coordinate zero. `CFG3D07` tests complete generated-table ownership with a
mixed positive-ion/negative-electron table and rejected count mismatch,
non-contiguous index, duplicate case-normalized symbol, zero mass, neutral
charge, out-of-range observer selection, and changed weight fingerprint.
`R3D05` additionally proves that identical total kinetic-energy bounds produce
different valid momentum intervals for proton and electron masses. These
extend the gates; none of the earlier CFG3D/MSH3D thresholds or negative
controls was relaxed.
`CFG3D08` exercises complete canonical input and three live schema-3 negative
controls. `MSH3D11` creates, verifies, and removes a real Tecplot product.
`MSH3D13` catches any future split between Parker points and tangents;
`MSH3D14` builds a mixed-level octree and rejects centerline gaps, approximate
halo depth, bounded inactive cavities, detached components, and continuation
beyond the configured finite line.
`MSH3D15` separates the fixed physical photosphere from the configurable
Parker/CME source shell and supplies live inside/surface/outside/malformed-box
controls for the predicate used by production leaf and cell masking.
`CFG3D09` and `TUR3D06` exercise both pre-existing turbulence-amplitude
prescriptions and the public total wave-energy output contract. `CFG3D10`
changes CME radius/direction independently to prove the optional launch-apex
link fails closed rather than moving the Parker tube implicitly.
`R3D08` constructs the canonical provider, checks delayed activation and the
full surface, and proves deterministic exact-count allocation plus downstream
no-cap behavior.
`R3D09` verifies finite serialization and the independent background,
sampling-window, and particle-occupancy flags for empty and occupied cells.
# Aggregate shared-model and live SEP + corona testing

From the AMPS root, after rebuilding the native executable:

```sh
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0
```

This runs both registries without naming individual tests (currently 222
shared-model plus 10 native cases), verifies complete fresh reports and writes
combined JSON/JUnit totals. `--list` discovers cases; `--model-only` explicitly
runs shared verification; `--require-no-skips` fails on unexercised prerequisites.
The supplied native deck uses Parker/SWCME, so generic host initialization PASS
does not establish the coronal runtime provider. See the
[coverage/CLI guide](../validation/COUPLED_SUITE_CLI.md).

The runner streams build/test output, prints each phase's start/end and counts
completed tests against the discovered registry, including percentages and
elapsed seconds. A heartbeat every 15 seconds keeps quiet phases visible;
append `--progress-interval 5` for more frequent updates. During MPI
initialization, this reports elapsed/idle time rather than estimating mesh
completion. Original subprocess logs and checked JSON remain authoritative.
Updating this Python progress reporting alone does not require rebuilding AMPS.

The final failure summary lists FAIL/ERROR cases by scope and ID, with a short
reason, a per-case diagnostic file and the execution log. It also prints all
generated phase-log paths and source-report paths. `failures.txt` beside the
latest `summary.json` is replaced each run; `runs/TIMESTAMP-ID/failures.txt`
and `failures/SCOPE/ID.log` preserve each invocation's summary and complete
failure details. Native case diagnostics contain messages/metrics/artifact
references, while full MPI output remains in `native.log`. Strict SKIPs and
infrastructure failures remain explicit. No AMPS rebuild is needed for this
Python reporting change. See the CLI guide for the complete log-path table
and a command to diagnose older reports without rerunning tests.

`test_coupled_sep_corona_runner.py` tests orchestration with deliberately
synthetic process fixtures: dynamic future membership, mixed suite exits,
missing/duplicate/extra/stale reports, provider/rank ownership, explicit scope
and SKIP policy. It also verifies that test progress appears before child exit,
quiet phases emit heartbeats, timeouts remain bounded and raw child logs remain
intact. Failure-summary coverage includes complete shared output, native
metrics/artifact references, generated-log paths, infrastructure errors,
strict SKIPs and replacing the latest summary while preserving older runs.
Those fixtures are never MPI/physical qualification evidence.

The current shared registry has 222 IDs through Stage 14; with seven generic
initialization and three SWCME mesh checks the aggregate has 232. The baseline Stage-13 subset has
209. Research campaign-protocol and portable kernel passes remain distinct
from actual observed-campaign/native-MPI qualification. No explicit test-ID
list is required to include future registered cases.

Reserved campaign records EVT3D01, XMD3D01 and SLM3D01 now retain passing
synthetic contract verification but report SKIP for absent actual campaign
evidence. Shared full selection is 222: 219 PASS, 3 SKIP, 0 FAIL. The
aggregate preserves shared SKIPs and enforces `--require-no-skips` across both
scopes. The Stage-13 baseline remains 209 PASS.


### Source-free shock propagation prerequisites

The standard runner automatically discovers `CME3D03` (schema-4 source-off
parser/factory, temporal-coverage conflicts, radius-stop boundaries and canonical
DBM versus an independent drag oracle) and `CME3D04` (contiguous telemetry,
clock/identity/MPI/source-off rejection and explicit close/overwrite behavior).
`CME3D01` includes launcher/input/checksum/raw-log protocol fixtures and PNG/EPS
rendering. They run with `--all` and `--suite phase-v`; no per-case names are
needed. They are portable software prerequisites, not evidence of native MPI
propagation. `CME3D02` remains the externally reviewed observational campaign.

Launch a native control through
`validation/run_swcme_coupling_validation.py --amps ../amps --ranks 10` from
srcSEP3D, or use the AMPS-root command in the validation README. The separate
`run_coupled_sep_corona.py` discovers shared coronal-model/native initialization
tests; it does not run this time-dependent SWCME observational campaign.

## SWCME mesh-background gates

The public application catalog includes eight `SWBG3D` portable cases: input and
frozen identity, exact RH vector/heated plasma endpoints and the inner handoff,
field evolution across a fixed cell, independent Cartesian derivative checks,
transactional invalid/empty-owner snapshots, future-model registration, polar
region snapshots and an independent analytic parallel RH limit.
They run with `--all`, `--group SWBG3D`, or `--suite phase-b`.

The portable `SWBG3D01–07` fixtures load
`examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in`. The former name
without `sphere` has been retired. Both the provider fixture and the parser's
text-mutation checks use one path constant in `test_swcme_background.cpp`;
`SWBG3D01` also verifies the canonical shape is `Sphere`. Install the renamed
example together with the test sources and rebuild `test/stage1`. These gates
test the spherical control, so substituting the finite-SSE deck changes their
intended geometry. `SWBG3D08` constructs its analytic RH control directly and
does not load either example.

The native sep-corona registry also contains `SWBGAMPS01–03`. Use
`examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in`, four ranks and
`--test-steps 2` with the coupled runner to check actual owner/native buffers,
cadence publication and received ghost blocks. These tests do not prepare or
write model fields during observation. The runner discovers them automatically;
an analytic-Parker input explicitly skips the SWCME-specific gates.
See [the update guide](../../SWCME_MESH_BACKGROUND_CHANGES.md) for full commands.

| Case | Oracle/evidence and prerequisite |
|---|---|
| `SWBG3D01` | Real schema-4 parser/factory; reject conflicting acceleration/source modes and a foreign canonical fingerprint |
| `SWBG3D02` | Exact canonical RH vector/pressure endpoint, declared heating partition and continuous ambient handoff; coupling verification, not an independent MHD solve |
| `SWBG3D03` | Fixed cell sampled before/after front passage; old snapshot stays immutable while U/B and generation advance |
| `SWBG3D04` | Smaller independent Cartesian stencil of canonical vectors and exact ambient `divU=2U/r` limit |
| `SWBG3D05` | Sentinel/pointer preservation on rejection, mixed valid/invalid points, empty owner snapshot and restored Parker counter |
| `SWBG3D06` | Registered extension constructed and published through the actual provider/builder/Runtime interfaces; duplicate/built-in protection |
| `SWBG3D07` | 1632 polar/oblique layer, sheath and ejecta samples at 0/60/120/4000 s through the actual snapshot builder, including Cartesian derivative stencils |
| `SWBG3D08` | Independent analytic switch-on compression/pressure/transverse-field magnitude; magnetic polarity and tangential boosts; conserved tangential fluxes |
| `SWBGAMPS01` | Live owner/application/native byte readback; requires SWCME mesh authority |
| `SWBGAMPS02` | Actual committed refreshes and installed epoch; additionally requires a crossed cadence |
| `SWBGAMPS03` | Live received-block representatives; additionally requires multiple ranks and received physical blocks |

Native checks cover mapped primitive/transport/E fields that are allocated and
both allocated DATAFILE time slots. Native current/electron-pressure slots are
filled but are not separately compared by these readback gates.
Remote coverage samples one physical center per received active block; it does
not claim exhaustive coverage of every ghost cell. The aggregate runner
discovers native descriptors after rebuilding the executable; portable SWBG3D
callbacks belong to the separate application `--all`/phase-b catalog.

```sh
python3 srcSEP3D/test/run_tests.py --suite phase-b --rebuild --output-dir test_output/background
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in --test-steps 2
```

See [background/README.md](../background/README.md) for the stencil and units,
[runtime/README.md](../runtime/README.md) for publication ordering, and the
[native guide](../validation/CORONAL_CME_NATIVE_TESTS.md) for evidence limits.

## Finite-SSE application acceptance

`python3 srcSEP3D/test/run_tests.py --group SSE3D --rebuild` builds and executes
the real application parser/providers/mover/checkpoint code without AMPS/MPI.
The nine cases are also included in the standalone, phase-b and phase-a suites.

| ID | Acceptance |
| --- | --- |
| SSE3D01 | Finite input, malformed axes/width, unsupported ellipsoid, tangent-flank inner handoff |
| SSE3D02 | 543 rotated/translated directions at three epochs against canonical radii/normals/speeds |
| SSE3D03 | Evolving canonical mesh primitives, ambient outside the cap and finite Cartesian derivatives |
| SSE3D04 | Rear/outside/duplicate rejection, moving center, oblique flank, tangent and SI-scale crossings |
| SSE3D05 | Actual canonical source preparation, physical patches confined to the cap, and a fixed probe inside the configured corridor |
| SSE3D06 | Schema-4 SSE geometry round trip and valid schema-3 spherical migration |
| SSE3D07 | Complete requested-time mover, expanding-cap subcycling and outside-cap exclusion |
| SSE3D08 | Exact mirrored failed coordinates and true weak layer through ten 60-s updates, including full derivative stencils |
| SSE3D09 | All canonical shapes: ambient support bypasses irrelevant RH; genuinely unresolved in-CME queries still reject without output writes |

`SSE3D05` converts the source-free propagation example into an injection
fixture. Injection requires an observer, but a fixed Earth position imported
from the analytic-Parker example need not lie inside this example's finite
active corridor. The fixture therefore computes a **test probe** at half the
configured Parker arc length, using the normalized source angles, inner
radius, rotation axis, wind speed and rotation rate. It adds the coordinate
origin and writes the fixed Cartesian position with 17-digit precision before
parsing the injection deck. This fixture works with both `parker-tube` and
`full-domain` allocation; it does not alter a real Earth's location or widen
the production corridor. With the default narrow corridor, a second probe in
the opposite direction must still be rejected by the configuration factory.

The propagation input remains source-free and needs no observers. For a real
injection campaign, choose the corridor's source angles and width to include
the actual fixed observer, or use `full-domain`; do not move an observational
position merely to satisfy an allocation check. The input's comments describe
how to switch the allocation mode independently of tube refinement.

To rerun just these fixtures from `AMPS/srcSEP3D`:

```sh
test/run_tests.py --group SSE3D --group SWBG3D --rebuild --output-dir test_output/swcme-fixtures
```

These are portable prerequisites. Re-run the native `sep-corona` suite on the
new SSE example to verify owner buffers, epoch scheduling and the actual MPI
receive halo. A portable PASS does not establish that native campaign, mesh
convergence or observational agreement.

The standalone halo-mask fixture can be run with
`python3 srcSEP3D/test/test_received_background_mask.py`. It compiles the actual
AMPS mask generator and native evidence-capture function, accepts current
received face layers in all six orientations, rejects stale received data, and
checks the null-mask full-block case. It does not start MPI. This specifically
regresses the earlier SWBGAMPS03 false failure caused by selecting an allocated
remote center whose packing-mask bit was zero.


SSE3D08 pins the failed `epoch_s=180` direction and a compression obtained
independently with 80-digit Decimal direct flux equations. Run the audit-only
reference verifier from the AMPS root:

```sh
python3 srcSEP3D/test/reference/verify_sse_weak_flank.py
```

This verifies a frozen reference and does not install or modify sources.
SSE3D09 deliberately retains a genuinely unresolved roundoff-scale jump;
ambient field sampling must succeed while the actual ICME layer and direct
shock diagnostic still fail. These tests complement the native ten-step MPI
run, which must be executed on the configured target build.
