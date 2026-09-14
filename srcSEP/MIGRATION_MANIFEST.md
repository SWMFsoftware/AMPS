# Step 14 mover migration manifest

Step 14 closes the temporary compatibility period and makes the source layout
match the public runtime contract. Only three particle movers are selectable:
`parker`, `fte-dmumu`, and `fte-mfp`. There is no fallback based on function
address, no mutable public mover pointer, and no retained legacy mover source.

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
