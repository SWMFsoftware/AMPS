# Step 14 mover migration manifest

## Stage-7 material checkpoint fixture migration, 2026-11-22

| Former contract | Maintained replacement | Compatibility and reason |
| --- | --- | --- |
| CP05C01–05 on a frozen analytic line with nonzero U | Same stable IDs on an explicitly owned moving reduced-front material line | Real host steps, nonempty accelerated population, prescribed scattering and fixed dt 300 s; old bare startup is ERROR before allocation. |
| CP05F01–03 protocol 1 fixed-line transport | Same IDs, strict protocol 2 with explicit startup deck/event and write/resume roles | Separate processes compare complete final files. Nondefault flux, malformed state and all old continuation assertions remain. |
| Wave/codec checks mixed with fixed-line particle motion | CP05O01–05 / CP05OF01–03 in `test/native_checkpoint_operators.cpp` | Native isolated operator/serialization coverage, with no particle or host stepping. Operator role protocol 1 remains supported here. |
| Orchestrator bare `physics-checkpoint-native` group | `stage7-checkpoint-native` complete 44-job gate | Every old C/F assertion is included. Python `--suite physics-checkpoint-native` aliases the same maintained gate; bare native cases require their explicit fixture. |
| Only pre-decode corruption rejection | MCP12 rejection after valid nonempty detached preparation on MPI 1/8 | Rejection-only diagnostic callback exposes counts, is forbidden in production purpose, leaves state/pointers/counters unchanged and permits a normal retry. |
| Cartesian remap from the old fractional coordinate | Derive Cartesian X from the exact newly stored translated line coordinate | Integer-prefix insertion can round a binary64 fraction. Coordinate and X must remain one canonical physical position; no clipping, timestep change or added bulk particle advection. |
| Threaded field-line dt read through an Eulerian block | Shared configured-global timestep selection for segment workers | Field-line MPI owners need not allocate the containing AMR block. Same dt as the serial path; node-attached block-local stepping is retained. |
| Material samples accumulated in mover arrival order | Read-only canonical physical-input sampling index | Exact sampled continuation without changing particle/list state, formulas or exposure. Reversing identical records fails the old native sampler and passes the repaired path. |
| Material serialization in native arrival order | Complete normalized records ordered by immutable particle ID | Keyed transport has no physical arrival order; all record bytes remain exact. Fixed-line shared-stream operator ordering is preserved. |
| Main/debug object dependencies declared before `MAINOBJ` | Define `MAINOBJ` before expanding the complete application object list | Prevents a worker/PIC header change from leaving stale weak inline pool definitions. A real Make negative control fails with the old ordering. |

Native material checkpoint schema remains 5. The diagnostics API is not a
production configuration option and does not change the on-disk layout or the
normal load transaction. Intrinsic material primitive validation precedes the
mandatory exact core-source comparison, so corrupt data receive precise errors
without changing the admitted set of checkpoints. See the maintained
[qualification report](validation/STANDALONE_INTEGRATION_STATUS.md) for exact
source identities and execution results. Historical build8 schema-4 results
below remain historical and do not certify this material path.

## 2026-11-22 material-continuation interfaces, under qualification

| Previous state/interface | Replacement and compatibility |
| --- | --- |
| Corrected physics checkpoint schema 4 | Schema 5: fieldwise material topology/identity/exposure, source ledgers and output phase. Earlier schemas are rejected; no invented history migration. Fixed analytic capability remains bounded separately. |
| Function-local field-line `lastDataOutputFileNumber` | Explicit `PIC::FieldLine::LastDataOutputFileNumber`, persisted for the material continuation owner. |
| `StoredDataBytes()` derived from collecting offset | Stable physical prefix excluding both diagnostic buffers, independent of their current phase. Signature unchanged. |
| `cFieldLine::Output(FILE*,bool)` local formatting | Third optional `collective` argument; existing two-argument source calls preserve local behavior. Collective public output explicitly supplies it. Rebuild callers because the linked C++ signature changes. |
| Six-row source-blind conservative merger | Additional immutable birth-energy row when all candidates carry reduced-source history. Legacy candidates keep their previous six rows. Mixed or missing reduced histories reject. |

## 2026-11-16 corrective interfaces and stored-state revisions

| Former interface/state | Replacement and compatibility |
| --- | --- |
| Post-pass `PICAdapter::RekeyClonedParticles()` | Explicit core `PrepareExactSplit`/`PostExactSplit` lifecycle callbacks; only actual children receive new identities. Internal post-hoc rekeying is removed. Ordinary migration/restart clones remain byte-exact. |
| Implicit per-call random-stream initialization | Particle transport schema 4 retains three keyed streams, cached normals, interval/event state, split serial and explicit physical birth window/remaining first duration; earlier histories are rejected under the user’s corrected-checkpoints-only policy. |
| Keyed-stream algorithm schema 1 | Schema 2 has a strictly open uniform mapping with zero-mantissa rejection; old algorithm state is rejected instead of claiming identical continuation. |
| Core turbulence checkpoint schema 2 | Schema 3 additionally persists magnetic polarity relative to geometric line coordinate; earlier layouts are rejected. This pure-core codec is distinct from native field-line physics restart. |
| SEP use of the core velocity-average merger | Optional `FledLine::MergeController` dispatches conservative positive-weight resampling for SEP. Other applications retain their default core policy. Surviving phase points/history remain unchanged; preserved and unpreserved moments are explicit. |

Native object consumers must rebuild after header layout changes. The SEP
makefile now emits/includes transitive dependency files and bootstraps objects
that lack them; generated AMPS configuration is preserved. The adopted polarized
MHD collision core and application checkpoint are under qualification, with
their current evidence boundaries recorded in the requirement matrix.

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

## Stage 3 release and validation migration

Stage 3 closes the gap between this disposition record and the physical
release tree. `mover.cpp` is no longer merely absent from `MAINLIBOBJ`; it is
physically absent and listed as retired in `SOURCE_MANIFEST.json`. The
AMPS-level `tools/sep_package_hygiene.py` auditor reads the manifests for both
applications and both shared-model directories in one transaction. It rejects
retired, generated, unclassified, or unsafe-symlink paths and can create a
metadata-normalized archive whose exact members are verified after writing.

The same stage promotes the existing observational and ensemble case sources
to discoverable registry entries without overstating their evidence:

| IDs | Registry group | Evidence disposition |
|---|---|---|
| `OV01`–`OV02` | `observational-validation` | release gates using immutable publication inputs and one declared shared normalization |
| `OV03`–`OV05` | `observational-validation` | diagnostic only; documented 1-D/connectivity limitations remain |
| `EV01`–`EV02` | `campaign-evidence` | incomplete calibration/holdout pilots |

`test-stage3-contracts` audits these mappings across the JSON registry, Python
runner, native descriptors, Make targets, manifests, and the independent
`srcSEP3D` species boundary. `test-ov01-ov05-unit` and
`test-ev01-ev02-unit` execute the focused case-specific source checks. A linked
OV run still requires `SEP_EXECUTABLE`; its absence is reported as `SKIP` and
cannot close a release gate.

## Native source-time amendment: transport schema 4

This amendment adds explicit source timing without changing a distribution or
inventing an omitted physical birth epoch. User-selected compatibility is the
corrected current model only; old/incompatible particle layouts are rejected.

| Boundary | Current implementation | Compatibility and ownership |
|---|---|---|
| Particle birth history | `Offset::BirthRecord`, padding-free `Transport::Birth::Record` (64 bytes); transport schema `0x5352435345500004` | Copies include origin, frozen source window, remaining first duration, and pending/consumed flag. Zero/old history cannot be reconstructed from a coordinate or current clock. |
| Host particle step | `PICAdapter::FreezeNativeDriverWindow` / `CloseNativeDriverWindow` in `birth_window_native.cpp` | Collective, before `PIC::TimeStep` and after joined particle work/read guard respectively. Standalone single-global, or one-species species-global with counter ON; other compiled clock/coupling contracts reject before mutation. Failed preflight/close preserves the existing state. |
| Scheduled initialization | `InitializeParticleTransportStateWithIdentity(..., const Birth::Record&)` | Caller supplies a valid explicit source record and local-plasma-frame state; it draws no source distribution and performs no inferred velocity boost. |
| Two-argument fresh initializer | Existing `InitializeParticleTransportStateWithIdentity(particle,id)` | Explicit initial conditions/controlled fixtures only; rejected inside an active source window. Existing particle load has no fresh-state fallback. |
| Three native mover adapters | `PrepareParticleStep`, `ParticleSnapshotRemaining`, `SelectParticleSubstep`, `CommitAndAttach` | Resolve the remaining first duration before draws; commit consumption with the successful phase-space/RNG/hazard state. Source helper pre-drift is not added. |
| Core clone/migration | Generic PIC byte-copy payload semantics retained | Statistical splitting and migration preserve original physical birth provenance, including pending state. Neither creates an artificial age or consumes a first move. Split lineage/future RNG remain separate semantics. |
| Conservative population reduction | Existing weight-only phase-point survivors | Surviving birth records remain unchanged; no velocity/epoch/age average is invented. |
| Native physics checkpoint | Completed-step capability, schema/layout checks | Active windows/pending first movement reject before publication or restore. Live four-polarization authority remains unsupported until the complete four-state codec/transaction is implemented and qualified. |
| Flux-owner persistence | Inner physics-checkpoint schema 4; `FluxTubeState::State` schema 1 and `flux_tube_checkpoint_state.h` staged owner API | Stores actual seed area/Φ/generation/provenance and exact lazy slots; callback pointers reject. Explicit numeric vertex XYZ pins decoded fixed geometry and particle Cartesian views; nonpositive drawn hazards reject. No same-process warm cache is treated as reconstruction evidence. |
| Stale PB::X after field-line motion | `CommitAndAttach` derives/validates canonical Cartesian position before any writes | Coordinate is sole authority; PB::X commits with phase/RNG/birth/list state. Nonfinite geometry rejects atomically. Material-geometry-only motion needs its own refresh transaction. |
| Raw SWMF Wave1/Wave2 treated as magnetic or direction-unknown | `ReadStandardBatsrusWaveView`, `util/sep_standard_batsrus_waves.*` | User-selected standard BATSRUS totals convert privately by half to magnetic energies. Wave1 stays along B, Wave2 opposite B; raw shared PIC storage and other apps are unchanged. Contract/generation/raw identity bind metadata; L/R remain absent. Full imported provider/handoff qualification remains separate. |
| `WN02` synthetic magnetic copy labelled SWMF | Explicit `unidentified-already-magnetic-view` metadata | Same physical values, reference and tolerance; generic unknown-direction controls retained. New SBW/SBWN gates own the actual standard raw-total import conversion. |
| Coefficient bridge recomputed radius as `|globalX|`; optional library received untranslated global position | Provider-defined radius/origin/known/frame and solar-relative position in `LocalBackgroundView` | Analytic uses the immutable declared origin; known SWCME/reduced frames stay global; unknown imported origins reject. The library receives provider-translated position with B/U basis unchanged. One radius-frame contract enters background fingerprint/provenance, so identical wire layout does not hide the physical interpretation change. RADIUSN01/02 explicitly opt schema6 `--test-input` into the same immutable startup pipeline. Both old-producer callbacks fail as intended on two ranks in BUILD10a; post-fix native gates pending. |
| Thermal source APIs | Typed preflight in `source_capabilities.cpp` | Main Maxwellian beginning and optional solar-wind thermal helpers reject before counter/RNG/allocation. User selected drifting relativistic Maxwell–Jüttner flux. The separate full-particle material-plane kernel now exists; native proper-density/frame/guiding-center boundary/weight/count/phase/geometry bridge remains gated. |
| Relativistic material-plane flux kernel | `util/sep_maxwell_juttner_flux.*` | Explicit proper n/T/m, full U, exact boundary=U and inward normal. Exact gamma-two/three mixture and cosine hemisphere, full Lorentz momentum, six staged keyed draws, no tail cutoff/redraw. This is not a replacement registration for the legacy native helper; tilted-cap gyroaveraged matching is unresolved. |
| Analytic moving-shock source epoch | `MovingShockDeparturePreflight` before source event reservation | Rejects mutable trajectory-query epoch different from frozen StepStart departure; fixed-boundary sources retain their declared timing. Pure-at-epoch query and separate end-step ledger remain pending. |

`BIRTH01`–`BIRTH08` are dependency-light timing/value regressions. The new
`BIRTHN01`–`BIRTHN06` pass **6/6 on one/two ranks in build6b**. Initial build6
errors exposed invalid same-epoch fixture republication; original evidence is
retained and fixtures now declare a later physical epoch without weakening the
publication invariant. Actual mover-manager/clone/split/MPI paths are covered;
the compiled `PIC::TimeStep` driver-counter gate is separate. Build8 repeats
the expanded native birth/Cartesian checks on two ranks with **6/6 PASS**.
Generic Maxwellian `SOURCE02` geometry/owner passes do not establish source
physics. [STATUS.md](STATUS.md) records current gates.

`FLCP01`–`FLCP04` pass sanitizer value/codec tests. Hardened schema-4 native
`CP05C01`–`CP05C05` pass **5/5 on two ranks in build8**, including decoded
Cartesian and strict hazard rejection. All three separate-process
`CP05F01`–`CP05F03` writer/reader pairs pass with byte-identical final artifacts.
The paired harness uses existing `--test-input`/`--test-output-dir`; it adds no
production parser or resume CLI. It compares full uninterrupted/resumed final
checkpoint artifacts and deliberately uses nondefault conserved flux, so a
missing owner restore cannot pass through accidental startup reconstruction.

Expanded `MJ01`–`MJ06` pass **6/6** with strict warnings, ASan/UBSan and applied
campaign/repeat/reversed-order controls. Independent positive radial integrals
and 90-digit references verify SI normalization; moving-plane quadrature and
low-temperature flux moments remain distinct from native source evidence.
The old code failed independent binary64 cold-boost crossing-loss checks;
numeric conditioning/encoded-phase precision rejection now passes the expanded
pure gate and **6/6 native cases on two ranks in build8**. `SBW01`–`SBW05` and
actual vertex-accessor `SBWN01`–`SBWN02` likewise pass **7/7** in that build.
No native thermal source bridge or full imported SWMF provider is enabled.
See [MAXWELL_JUTTNER_FLUX_IMPLEMENTATION.md](validation/physics/MAXWELL_JUTTNER_FLUX_IMPLEMENTATION.md)
for numerical limits, physical scope and exact commands. No classical fallback
or general nonmaterial-plane reuse is introduced.

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
