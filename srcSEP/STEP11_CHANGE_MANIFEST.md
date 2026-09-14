# Step 11 change manifest: authoritative turbulence subsystem

## Implemented

- Added `util/sep_turbulence_core.*` with five source modes and an independently
  validated integrated-or-spectral authority.
- Added explicit restartable boundary objects, source/representation/coupling
  parsers, one installed active configuration, and local-evolution ownership
  gating in the standalone loop.
- Added derived density, amplitude, cross-helicity, resonance, `Dmumu`, mean
  free path, and spatial-diffusion views without secondary mutable ownership.
- Added the fixed source/coupling/advection/reflection/cascade/sync driver order,
  explicit CFL subcycling, positivity limiting, and structured status returns.
- Added a signed energy ledger, conservative arc-length remap, one-time SWMF
  handoff metadata, deterministic checkpoints, and stable state hashes.
- Replaced production reflection/cascade constants with installed CLI values
  and changed the invalid CFL fallback into an explicit fatal error.
- Added CLI controls for source, representation, coupling, operator flags,
  boundaries, coefficients, correlation length, spectral grid, cadence, CFL,
  and conservation tolerance.
- Added the original focused turbulence cases, a CLI validation test, source gates, and
  `make test-turbulence-core-unit`.
- The controlled-validation enhancement registers `TURB02`–`TURB20` under
  their stable IDs, retains source ownership as `TURBOWN01` to avoid collision
  with the established production `TURB01`, and adds `TURB21` translated-sine
  advection plus `TURB22` time-dependent analytical wave-energy growth. The
  source-only runner consumes these same descriptors and emits JSON/JUnit
  evidence. `TURB22` supplies a known source to isolate its application and
  ledger accounting; deriving a QLT source from particle distributions remains
  the responsibility of the native `growth_rate_validation_test.cpp` path.
- Added `TURB23`, which balances the weighted relativistic kinetic-energy
  change from controlled plus/minus-branch wave-frame scattering against the
  signed wave source and total-energy ledger.

No-argument behavior remains self-consistent integrated turbulence with the
existing operators enabled. The subsystem remains independent of the
three-mover registry.

The source archive lacks the enclosing AMPS `Makefile.conf`, PIC headers, and
MPI runtime. Native standalone/SWMF runs must still provide process-level
restart, remap, and long-duration conservation evidence.
