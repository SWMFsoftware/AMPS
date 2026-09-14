# Step 14 mover migration manifest

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
