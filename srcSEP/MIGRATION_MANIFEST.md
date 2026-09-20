# Step 14 mover migration manifest

## B01 source-tree rebaseline

The release manifest now matches the physical tree. Generated products are
excluded by `SOURCE_MANIFEST.json`, local ignore rules, and the AMPS-level
`tools/sep_package_hygiene.py` gate. The following unlinked historical files
were removed only after the makefiles and source references were audited:

| Removed source | Active owner or replacement |
| --- | --- |
| `mover.cpp` | `parker_mover.cpp`, `focused_transport_dmumu.cpp`, `focused_transport_mfp.cpp`, and `mover_state.cpp` |
| `fte_mover.cpp`, `fte_mover_dmumu.cpp` | the two canonical focused-transport adapters and tested pure cores |
| `drift.cpp` | outside the field-line-only `srcSEP` scope; reserved 3-D transport interfaces are owned by `srcSEP3D` |
| `output.cpp`, `sample3d*` | field-line `sampling.cpp`, `sampling_output.cpp`, `output_fl_background.cpp`, and the independent 3-D application's output layer |
| `demo1d.cpp`, `sw1d.cpp` | canonical SWCME examples under `src/models/swcme` and the private `adapters/swcme1d_adapter.cpp` consumer |

Object/dependency/archive files, native test executables, `test_output`, Python
caches, AMPS `amr.sig` data, and historical `PT` plot/log output were also
removed. They are generated evidence, not source inputs.

## B02 canonical SEP-common relocation

| Former application-local source | Canonical owner | Build/test behavior |
|---|---|---|
| `util/sep_transport_common.*` | `src/models/sep_common/sep_transport_common.*` | canonical public include; object inserted once into `mainlib.a` |
| `util/sep_coefficient_physics.*`, `util/sep_coefficient_registry.*` | matching files in `src/models/sep_common` | one coefficient implementation and registry across applications |
| `util/sep_background_snapshot.*` | matching files in `src/models/sep_common` | immutable snapshot ABI comes from the shared archive |
| `util/sep_test_registry.*` | matching files in `src/models/sep_common` | standalone and production reports use one result schema |
| `util/sep_injection_spectrum.*`, `util/sep_species_source.*` | matching files in `src/models/sep_common` | injection measure and species normalization cannot drift by application |

There are no forwarding headers or `.cpp` wrappers.  The makefile resolves the
canonical directory from the active source/copied-build location, builds its
archive, inserts its exact seven objects, and audits member counts plus strong
definitions. `test/run_sep_common_ownership_tests.sh` also links a public
consumer and constructs synthetic `srcSEP` and `build/main` layouts so a
working-directory-dependent path cannot pass accidentally.

`util/sep_common_header_path.h` is the sole include-path resolver, not a copy of
any canonical declaration. In an installed tree it prefers the canonical
AMPS-root path, preventing a stale short-name header from winning by include
order; detached component builds use `-I$(SEP_COMMON_DIR)`. When an AMPS sibling
library includes `build/main/sep.h` through `pic.h`, it has the global
`-I$(AMPS_ROOT)` but cannot inherit variables from the build/main child make;
the resolver therefore selects `src/models/sep_common/<header>` from that root.
Production translation units copied to `build/main` have the same constraint,
so `field_line.cpp`, the diffusion implementations, the private SWCME adapter,
and linked validation sources use the resolver as well. Focused test sources
may retain canonical short names because their compile commands explicitly add
`-I$(SEP_COMMON_DIR)`. The B02 gate compiles both contexts, rejects direct path
fallbacks anywhere outside this resolver, and rejects bare shared-header
includes from every non-test `.cpp` file.

## B04 external-suite result migration

The historical Python runner propagated a failed Make return code but omitted
that suite from merged JSON/JUnit when no child report existed. Every selected
source suite now owns a `SUITE-*` wrapper result with command, timestamps,
elapsed time, return code, and bounded diagnostic output. Wrapper and child
records are merged before any summary is written. Compile/assertion exits are
`FAIL`; launch failures, timeout, and signals are `ERROR`; explicit or blocked
skips remain `SKIP`. Direct launchers may return 77; Make-backed gates emit the
exact `SRCSEP_SUITE_RESULT=SKIP` line and return zero because Make rewrites a
recipe exit 77 as failure. Exit precedence is ERROR/status 2, then FAIL/status
1, then PASS-or-SKIP/status 0. Runner self-tests cover every state and a real
invalid-C++ compilation.

## B03 canonical SWCME relocation

| Former application-local surface | Canonical replacement | Compatibility |
|---|---|---|
| `srcSEP/swcme/` or transitional `srcSEP/swcme----moved-out/` | `AMPS/src/models/swcme/` | no application-local forwarding tree |
| public or application-local `#include "swcme/swcme1d.hpp"` | private `adapters/swcme1d_adapter.cpp` includes `swcme1d.hpp` with `-I$(SWCME_DIR)` | `sep.h` is provider-free; source and copied `build/main` layouts supported |
| `#include "../../swcme/swcme_sep_interface.hpp"` in VAL04 | `#include "swcme_sep_interface.hpp"` with the canonical test include root | no relative sibling dependency |
| `make -C swcme/test` | canonical `$(SWCME_DIR)/test/run_tests.py --routine` | test failures continue to propagate |
| root `srcSEP/demo1d.cpp` | `src/models/swcme/demo1d.cpp` | duplicate application copy removed |

The makefile now uses the same active-makefile path discovery as srcSEP3D,
builds `$(SWCME_DIR)/swcme.a`, and inserts canonical `swcme3d.o` exactly once
into `mainlib.a`. `test/run_swcme_relocation_tests.sh` enforces ownership,
the srcSEP3D-style private-adapter boundary, and source-versus-build/main path
equivalence.

## B05 WP42–WP64 disposition and promotion boundary

The previous overlay could compile prototypes while its prose implied that
they were active production algorithms. B05 replaces that ambiguity with
[`WP42_WP64_DISPOSITION.json`](WP42_WP64_DISPOSITION.json), an exact inventory
that must classify every WP42–WP64 item and state why it has that status, which
symbol or gate owns it, and what evidence is required for promotion.

| Classification | Work packages | Build and evidence behavior |
|---|---|---|
| `experimental-component` | WP42–WP58, WP61–WP62 | four extension `.cpp` files compile under strict warnings and sanitizers and execute controlled component tests, but are absent from `MAINLIBOBJ` and have no production selector/call path |
| `external-gate` | WP59, WP60, WP63, WP64 | native AMPS/SWMF, external-data, or scaling evidence is required; an unconfigured gate produces an explicit B04 `SKIP` record |

Two supporting contracts were accepted independently of the experimental
overlay. `ParkerMeasure` now names the distribution measure in the canonical
SEP-common transport header. Event-driven wave coupling now carries the
actual resonant branch and the pre/post-event momenta produced by the completed
MFP event; the active adapter no longer substitutes limiter-shell endpoints.
These changes are audited by source reachability plus the focused event tests.

`test-wp42-wp64-experimental` verifies the 23 component contracts and the
machine-readable disposition, including a deliberately corrupted negative
control. `test-wp59-wp64-native` invokes an explicitly supplied
`SRCSEP_NATIVE_GATE`, or emits the exact `SRCSEP_SUITE_RESULT=SKIP` marker when
no reviewed command exists. Neither result silently inserts an experimental
object into `mainlib.a`. Promotion requires a reviewed production selector,
restart behavior, native observation, and coordinated updates to the source
manifest, archive audit, disposition, and documentation.

Step 14 closes the temporary compatibility period and makes the source layout
match the public runtime contract. Only three particle movers are selectable:
`parker`, `fte-dmumu`, and `fte-mfp`. There is no fallback based on function
address, no mutable public mover pointer, and no retained legacy mover source.

## WP01--WP10 production-contract replacements

| Former production surface | Replacement | Compatibility |
|---|---|---|
| inline turbulence mutation sequence in `main.cpp` | `SEP::Turbulence::PICAdapter::Advance` | former block retained constant-false for review; remove after native parity gate |
| `QueueAveragedWaveContribution` / `QueueFocusedWaveContribution` with a PIC handle | typed `QueueWaveContribution(const CouplingRecord&)` | no live-pointer compatibility path in the production queue |
| `EvaluateLocalBackground(context, particle_dt, ...)` | `EvaluateLocalBackgroundAt(context, relative_arc_length_m, ...)` and zero-displacement wrapper | particle-timestep overload removed because its derivative semantics were invalid |
| sign(mu) / 50:50 event branch inference | provider `nuPlusPerS` / `nuMinusPerS` | lambda-only providers map explicitly to balanced rates |
| resampled waiting time at every limiter shell | carried `remainingOpticalDepth` and `nextEventIndex` | NaN state initializes a new/backward-compatible particle history |

Full details and validation boundaries are in
[WP01_WP10_IMPLEMENTATION.md](WP01_WP10_IMPLEMENTATION.md).

## WP11--WP20 physics and configuration replacements

| Former behavior | Required replacement |
|---|---|
| clipped or unconstrained stochastic pitch update | bounded reflecting Milstein scheme in `sep_focused_transport_core` |
| anonymous mover substep fractions | validated `NumericalTolerances` and named limiter diagnostics |
| coefficient inputs gathered from unrelated globals | immutable source/representation/generation/checksum `LocalInputView` |
| constant callback that could omit its configured value | pure `EvaluateConstantDmumu` with exact zero derivative |
| independently coded Jokipii value and derivative | one finite-band piecewise kernel with typed transition state |
| Florinskiy missing assignment/unbounded derivative/duplicated denominator | assigned output, branch-specific denominators, and bounded stencil |
| fixed Dmumu-to-kappa quadrature and numerical floors | adaptive error-controlled quadrature with typed gap policy |
| embedded proton constants and common species source | explicit `SpeciesProperties` and normalized `SpeciesSource` records |
| QLT1 sampling field and anonymous turbulence scales | selected local field plus named spectrum/correlation configuration |
| untyped `+infinity` reaching any mover | typed ballistic MFP and mover/provider preflight matrix |

See [WP11_WP20_IMPLEMENTATION.md](WP11_WP20_IMPLEMENTATION.md) for equations,
configuration names, focused evidence, and limitations.

## WP21–WP30 physics/output/run-control replacements

| Former behavior | Replacement |
|---|---|
| forward-Euler shock radius and static last-time | exact restartable knot trajectory in `sep_shock_source_core` |
| first-inside-only quadratic intersection | oriented full-polyline sphere intersection with typed policies |
| ambiguous log-momentum weights and Sokolov NaN normalizer | named normalized spectrum measure and inverse CDF |
| global source RNG / commutative XOR keys | tagged ordered semantic source keys |
| direct shock mutation of wave datum | queued `pendingShock` records consumed by the common ledger |
| one implicit tube reference area | per-line `MagneticFluxRecord` with generation/provenance |
| diagnostic superluminal cap / edge-bin clamping | strict exclusion and separate invalid/underflow/overflow counters |
| integral-then-peak normalized output | distinct physical products or explicitly unit-integral shape |
| shell mkdir, fixed `sprintf`, direct final writes | validated transactional output plus checksum completion manifest |
| scattered driver literals/static init flags | frozen, fingerprinted `RunConfiguration` and CLI controls |

See [WP21_WP30_IMPLEMENTATION.md](WP21_WP30_IMPLEMENTATION.md) and
[WP21_WP30_VALIDATION_REPORT.md](WP21_WP30_VALIDATION_REPORT.md).

## WP31–WP41 numerical and evidence replacements

| Former behavior or gap | Replacement |
| --- | --- |
| advection-only CFL with full-step local operators | `PlanAdvance` shared advection/source/reflection/cascade limits |
| silent wave floors, skips, and `mu` resonance epsilon | typed events, signed rejected source, explicit `mu=0` failure |
| fixed merge/split bin and threshold arguments | frozen `populationControl` policy passed to PIC calls |
| source-only checks described as production tests | fail-closed `NativeObservation` contract and evidence levels |
| selected compatibility examples | generated 90-row matrix and stable preflight codes |
| separate restart/remap/conservation assertions | checksummed global seam ledger checkpoint |
| one-seed stochastic claims | versioned domain-separated seed panels and confidence statistics |
| example-only edge tests | bounded IEEE/property cases and deterministic fault injection |
| normalized curves compared directly with instruments | SI response-folded count forward model |
| timing without algorithmic attribution | exact work counters plus environment-bound timing |
| README capability keywords as evidence | validated claim level, command, artifact, and status |

See [WP31_WP41_IMPLEMENTATION.md](WP31_WP41_IMPLEMENTATION.md) and
[WP31_WP41_VALIDATION_REPORT.md](WP31_WP41_VALIDATION_REPORT.md).

## Runtime-name replacements

| Retired input name | Required canonical replacement | Reason |
|---|---|---|
| `parker-dxx`, `parker-diffusion`, `dxx` | `parker` | The Parker mover's coefficient contract is spatial diffusion. |
| `fte`, `default`, `focused-transport`, `focused-transport-equation`, `legacy-fte` | `fte-dmumu` | The canonical name makes the pitch-angle-diffusion closure explicit. |
| `focused-transport-event-driven`, `event-driven`, `event-driven-fte`, `fte-event-driven` | `fte-mfp` | The canonical name makes the event-driven mean-free-path closure explicit. |
| `parker-mean-free-path` | Choose `parker` or `fte-mfp` explicitly | The old phrase mixed an averaged equation with an event-scattering closure and is intentionally ambiguous. |
| `mean-free-path-scattering`, `tenishev-2005-fl` | `fte-mfp` plus `--mean-free-path-provider <name>` | Algorithm selection and coefficient selection are independent. |
| `coupled-fte`, `focused-transport-wave-scattering` | `fte-dmumu` plus explicit turbulence source/coupling options | Self-consistent turbulence is a subsystem, not a fourth mover. |
| Parker3D, He-2019, Kartavykh-2016, Borovikov-2019, drift, or Boris spellings | Move the run to the separate Cartesian-transport application | Those algorithms require 3-D particle state deliberately absent from srcSEP. |

Unknown and retired names now return a CLI error. Migration is therefore
visible in batch logs instead of silently changing the requested physics.

## C++ symbol and source replacements

| Removed interface/source | Destination or replacement |
|---|---|
| `SEP::ParticleMoverPtr` and `SEP::fParticleMover` | `SEP::Mover::SelectProductionMover` and `SEP::Mover::DispatchProductionMover` |
| `ParticleMover_Parker_Dxx`, `ParticleMover_ParkerEquation` | `ParticleMover_Parker` in `parker_mover.cpp` |
| `ParticleMover_FTE`, `ParticleMover_Droge_2009_AJ` | `ParticleMover_FocusedTransport_Dmumu` in `focused_transport_dmumu.cpp` |
| `ParticleMover_MeanFreePathScattering`, `ParticleMover_Tenishev_2005_FL`, `ParticleMover_He_2011_AJ`, `ParticleMover_Parker_MeanFreePath` | `ParticleMover_FocusedTransport_EventDriven` in `focused_transport_mfp.cpp`, with an explicitly selected MFP provider where applicable |
| monolithic `mover.cpp` | the three canonical mover files plus shared `mover_state.cpp` |
| legacy `fte_mover.cpp` and `fte_mover_dmumu.cpp` | `focused_transport_mfp.cpp` and `focused_transport_dmumu.cpp` |
| `AccountTransportCoefficient`, `LimitScatteringUpcomingWave`, `NumericalScatteringEventMode`, `NumericalScatteringEventLimiter` | removed; these had no active canonical-mover behavior |

`MaxTurbulenceLevel`, `MaxTurbulenceEnforceLimit`, `LimitMeanFreePath`, and
`AccountAdiabaticCoolingFlag` remain in `mover_state.cpp` because canonical
coefficient/mover code still consumes them.

## Stage 1–2 runtime ownership migration

| Former behavior | Current owner and invariant |
|---|---|
| Standalone `main.cpp` advanced turbulence after `amps_time_step()` had already done so | `amps_time_step()` is the sole transaction owner for standalone and coupled execution; it advances exactly once after the immutable particle phase. |
| Shock injection history existed only in the standalone driver | `main_lib.cpp` keeps provider-neutral previous/current shock radii at the common transaction boundary. |
| `ModelInit::Init()` and a second wave initializer both wrote startup turbulence | Only self-consistent integrated/spectral sources receive the local initializer. Prescribed, SWMF-read-only, and SWMF-handoff inputs retain provider authority. |
| The radial density helper overwrote every provider | It is restricted to the analytic provider; SWCME and SWMF plasma state is preserved. |
| Native mover fixtures called production movers outside a snapshot read phase | `FTE01`, `PARKER01`, and `PARKER02` hold an immutable generation around mover calls and release it before restoring fixture data. |
| Detached CV01 compilation depended on incidental include state | `run_cv01_tests.sh` supplies the canonical repository root explicitly and uses a disposable object directory. |

The dependency-light enforcement command is
`make test-stage1-stage2-contracts`; the numerical companion gates are
`test-state-unit`, `test-turbulence-core-unit`, and `test-cv01-unit`.

## Verification

Run the dependency-light migration and documentation gate:

```sh
make test-documentation-unit
```

Then, in a configured AMPS checkout, build with strict warnings and run the
native acceptance suite:

```sh
make WARNINGS='-Wall -Wextra -Wpedantic -Werror' lib amps
make test SEP_EXECUTABLE=../amps
```

The source archive cannot truthfully claim the second gate when AMPS-generated
headers, MPI, or the linked executable are unavailable.
