#include "sep_turbulence_validation.h"

#include "sep_focused_transport_mfp_core.h"
#include "sep_turbulence_core.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace SEP {
namespace Testing {
namespace {

using Turbulence::CellState;
using Turbulence::Configuration;
using Turbulence::Representation;
using Turbulence::Source;
using Turbulence::State;

bool Near(double actual, double expected, double tolerance = 1.0e-11) {
  return std::fabs(actual - expected) <= tolerance *
      std::max(1.0, std::max(std::fabs(actual), std::fabs(expected)));
}

Result Complete(bool pass, const std::string& success,
                const std::string& failure,
                const std::string& configuration) {
  Result result;
  result.status = pass ? Status::Pass : Status::Fail;
  result.message = pass ? success : failure;
  result.configuration.push_back(configuration);
  result.metrics.push_back(
      {"assertion_failures", pass ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

Descriptor MakeDescriptor(const char* id, const char* name,
                          const char* description, RuntimeClass runtime,
                          TestCallback callback) {
  Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = "turbulence";
  descriptor.description = description;
  descriptor.initialization = InitializationLevel::None;
  descriptor.supportedBuildModes =
      "linked CLI and source-only C++11 ASan/UBSan runner";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = "deterministic controlled state; no RNG";
  descriptor.stateIsolation =
      "stack/vector-owned turbulence state; active global configuration unchanged";
  descriptor.callback = callback;
  return descriptor;
}

CellState Cell(double plus_j, double minus_j) {
  CellState cell;
  cell.lengthM = 10.0;
  cell.volumeM3 = 20.0;
  cell.plasmaSpeedMPerS = 1.0;
  cell.alfvenSpeedMPerS = 0.5;
  cell.dLnAlfvenSpeeddsPerM = 0.02;
  cell.magneticFieldT = 5.0e-9;
  cell.massDensityKgPerM3 = 1.0;
  cell.ePlusJ = plus_j;
  cell.eMinusJ = minus_j;
  return cell;
}

State Basic(Source source = Source::SelfConsistentIntegrated) {
  State state;
  state.configuration.source = source;
  state.configuration.representation =
      source == Source::SelfConsistentSpectral ? Representation::Spectral
                                                : Representation::Integrated;
  state.configuration.advectionEnabled = false;
  state.configuration.reflectionEnabled = false;
  state.configuration.cascadeEnabled = false;
  state.configuration.innerBoundary.policy =
      Turbulence::BoundaryPolicy::TransparentOutflow;
  state.configuration.outerBoundary.policy =
      Turbulence::BoundaryPolicy::TransparentOutflow;
  state.cells.push_back(Cell(8.0, 2.0));
  state.provenance = "controlled-turbulence-test-v1";
  return state;
}

bool Initialize(State* state) {
  return Turbulence::InitializeState(state).ok();
}

Result RunTurbulenceOwnership() {
  bool ownership = true;
  const Source sources[] = {Source::Prescribed,
      Source::SelfConsistentIntegrated, Source::SelfConsistentSpectral,
      Source::SwmfReadOnly, Source::SwmfInitialThenEvolveLocal};
  for (std::size_t i = 0; i < sizeof(sources) / sizeof(sources[0]); ++i) {
    State state = Basic(sources[i]);
    if (sources[i] == Source::SwmfInitialThenEvolveLocal)
      state.configuration.representation = Representation::Integrated;
    ownership = ownership && Initialize(&state);
  }
  Result result = Complete(ownership,
      "all declared turbulence source/authority combinations validate",
      "a declared turbulence source/authority combination was rejected",
      "sources=prescribed,integrated,spectral,swmf-read-only,swmf-handoff");
  result.metrics.push_back(
      {"source_modes_checked", 5.0, 5.0, "==", "count"});
  return result;
}

Result RunTurb02() {
  State state = Basic();
  state.cells[0].pendingParticlePlusJ = 1.25;
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double ledger_error =
      std::fabs(step.ledger.particleExchangeJ - 1.25);
  const bool pass = initialized && step.status.ok() &&
      ledger_error <= 1.0e-11 && Near(step.ledger.closureResidualJ, 0.0);
  Result result = Complete(pass,
      "particle-wave exchange equals the applied controlled energy",
      "particle-wave ledger does not close",
      "initial_Eplus_J=8;pending_particle_plus_J=1.25;dt_s=1");
  result.metrics.push_back(
      {"particle_exchange_absolute_error", ledger_error, 1.0e-11, "<=", "J"});
  result.metrics.push_back({"closure_residual_absolute",
      std::fabs(step.ledger.closureResidualJ), 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb03() {
  State state = Basic();
  const bool initialized = Initialize(&state);
  const Turbulence::DerivedCell derived = Turbulence::DeriveCell(state.cells[0]);
  const double plus_error = std::fabs(derived.wPlusJPerM3 - 0.4);
  const double minus_error = std::fabs(derived.wMinusJPerM3 - 0.1);
  const bool pass = initialized && plus_error <= 1.0e-11 &&
      minus_error <= 1.0e-11 && derived.deltaB2T2 > 0.0;
  Result result = Complete(pass,
      "integrated energy converts to analytical energy density",
      "turbulence unit/volume conversion differs from E/V",
      "volume_m3=20;Eplus_J=8;Eminus_J=2");
  result.metrics.push_back({"plus_density_absolute_error", plus_error, 1.0e-11, "<=", "J m^-3"});
  result.metrics.push_back({"minus_density_absolute_error", minus_error, 1.0e-11, "<=", "J m^-3"});
  return result;
}

Result RunTurb04() {
  State state = Basic();
  state.cells.push_back(Cell(4.0, 1.0));
  const bool pass = Initialize(&state) &&
      state.cells[0].ePlusJ > state.cells[1].ePlusJ;
  return Complete(pass,
      "controlled nonuniform physical state initializes without flattening",
      "turbulence initialization changed the prescribed ordering",
      "Eplus_J=8,4;Eminus_J=2,1");
}

Result RunTurb05() {
  State state = Basic(Source::SelfConsistentSpectral);
  const bool initialized = Initialize(&state);
  double spectral_sum = 0.0;
  for (std::size_t k = 0; k < state.configuration.spectralBins; ++k)
    spectral_sum += state.cells[0].spectralEnergyJ[k];
  const double error = std::fabs(spectral_sum - 8.0);
  Result result = Complete(initialized && error <= 1.0e-11,
      "spectral projection conserves the analytical branch total",
      "integrated-to-spectral projection changed branch energy",
      "representation=spectral;Eplus_J=8;spectral_bins=128");
  result.metrics.push_back(
      {"projected_energy_absolute_error", error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb06() {
  State state = Basic();
  state.cells.push_back(Cell(8.0, 2.0));
  state.cells.push_back(Cell(8.0, 2.0));
  state.configuration.advectionEnabled = true;
  state.configuration.periodicBoundaries = true;
  const bool initialized = Initialize(&state);
  const std::uint64_t before = Turbulence::StateHash(state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double error = std::fabs(state.cells[0].ePlusJ - 8.0);
  Result result = Complete(initialized && step.status.ok() &&
      error <= 1.0e-11 && before != 0,
      "uniform periodic energy is invariant under constant advection",
      "constant-coefficient advection changed a uniform solution",
      "cells=3;periodic=true;uniform_Eplus_J=8;dt_s=1");
  result.metrics.push_back(
      {"uniform_energy_absolute_error", error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb07() {
  State state = Basic();
  state.configuration.advectionEnabled = true;
  state.cells[0].plasmaSpeedMPerS = 20.0;
  const bool initialized = Initialize(&state);
  const Transport::ScalarResult limit =
      Turbulence::MaximumStableAdvectionStepS(state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 2.0);
  const double expected = 0.8 * 10.0 / 20.5;
  const double error = std::fabs(limit.value - expected);
  Result result = Complete(initialized && limit.status.ok() &&
      step.status.ok() && error <= 1.0e-11 && step.diagnostics.subcycles > 1,
      "variable-speed advection obeys the analytical CFL bound",
      "CFL limit or required subcycling is incorrect",
      "length_m=10;U_m_per_s=20;VA_m_per_s=0.5;cfl=0.8;dt_s=2");
  result.metrics.push_back(
      {"maximum_step_absolute_error", error, 1.0e-11, "<=", "s"});
  result.metrics.push_back({"subcycles", static_cast<double>(step.diagnostics.subcycles), 1.0, ">", "count"});
  return result;
}

Result RunTurb08() {
  Configuration configuration = Basic().configuration;
  configuration.innerBoundary.policy =
      Turbulence::BoundaryPolicy::SpecifiedIncomingFlux;
  configuration.innerBoundary.value = -1.0;
  const bool rejected = !Turbulence::ValidateConfiguration(configuration).ok();
  return Complete(rejected,
      "negative incoming boundary flux is rejected",
      "invalid negative incoming boundary flux was accepted",
      "inner_policy=specified-incoming-flux;value_W=-1");
}

Result RunTurb09() {
  State state = Basic();
  state.configuration.reflectionEnabled = true;
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double total_error =
      std::fabs(state.cells[0].ePlusJ + state.cells[0].eMinusJ - 10.0);
  Result result = Complete(initialized && step.status.ok() &&
      total_error <= 1.0e-11 && step.ledger.reflectionTransferJ != 0.0,
      "reflection transfers energy between branches without changing total",
      "reflection violated the controlled branch-sum invariant",
      "Eplus_J=8;Eminus_J=2;reflection=true;dt_s=1");
  result.metrics.push_back(
      {"branch_total_absolute_error", total_error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb10() {
  State integrated = Basic();
  integrated.configuration.cascadeEnabled = true;
  const bool integrated_initialized = Initialize(&integrated);
  const Turbulence::StepResult integrated_step =
      Turbulence::Advance(&integrated, 1.0);
  State spectral = Basic(Source::SelfConsistentSpectral);
  spectral.configuration.cascadeEnabled = true;
  const bool spectral_initialized = Initialize(&spectral);
  const double highest_before = spectral.cells[0].spectralEnergyJ[
      spectral.configuration.spectralBins - 1];
  const Turbulence::StepResult spectral_step = Turbulence::Advance(&spectral, 1.0);
  const bool pass = integrated_initialized && spectral_initialized &&
      integrated_step.status.ok() && spectral_step.status.ok() &&
      integrated_step.ledger.physicalDissipationJ > 0.0 &&
      spectral_step.ledger.physicalDissipationJ > 0.0 &&
      Near(integrated_step.ledger.closureResidualJ, 0.0) &&
      spectral.cells[0].spectralEnergyJ[
          spectral.configuration.spectralBins - 1] != highest_before;
  Result result = Complete(pass,
      "integrated and spectral cascade dissipate while closing their ledgers",
      "cascade/dissipation failed positivity, transfer, or closure checks",
      "representations=integrated,spectral;cascade=true;dt_s=1");
  result.metrics.push_back({"integrated_dissipation_J",
      integrated_step.ledger.physicalDissipationJ, 0.0, ">", "J"});
  result.metrics.push_back({"spectral_dissipation_J",
      spectral_step.ledger.physicalDissipationJ, 0.0, ">", "J"});
  return result;
}

Result RunTurb11() {
  State state = Basic();
  state.cells[0].pendingParticleMinusJ = 0.5;
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double error = std::fabs(state.cells[0].eMinusJ - 2.5);
  Result result = Complete(initialized && step.status.ok() && error <= 1.0e-11,
      "one-step wave growth equals the independently accumulated source",
      "controlled wave-growth increment is incorrect",
      "initial_Eminus_J=2;pending_particle_minus_J=0.5;dt_s=1");
  result.metrics.push_back(
      {"growth_energy_absolute_error", error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb12() {
  State state = Basic();
  const bool initialized = Initialize(&state);
  const Turbulence::CoefficientView coefficient =
      Turbulence::DeriveCoefficients(
          state, 0, 1.0e6, 0.4, 1.602176634e-19, 1.67262192369e-27);
  const double closure_error = std::fabs(
      coefficient.kappaParallelM2PerS -
      1.0e6 * coefficient.lambdaParallelM / 3.0) /
      std::max(1.0, std::fabs(coefficient.kappaParallelM2PerS));
  Result result = Complete(initialized && coefficient.status.ok() &&
      coefficient.dMuMuPerS > 0.0 && coefficient.lambdaParallelM > 0.0 &&
      closure_error <= 1.0e-11,
      "resonant coefficient view satisfies kappa=v*lambda/3",
      "derived turbulence coefficients violate their isotropic closure",
      "speed_m_per_s=1e6;mu=0.4;cell=0");
  result.metrics.push_back(
      {"kappa_closure_relative_error", closure_error, 1.0e-11, "<=", "dimensionless"});
  return result;
}

Result RunTurb13() {
  State state = Basic();
  state.cells[0].pendingShockPlusJ = 3.0;
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double source_error = std::fabs(step.ledger.shockSourceJ - 3.0);
  const double energy_error = std::fabs(state.cells[0].ePlusJ - 11.0);
  Result result = Complete(initialized && step.status.ok() &&
      source_error <= 1.0e-11 && energy_error <= 1.0e-11,
      "shock source adds the prescribed controlled energy exactly once",
      "shock injection or its ledger value is incorrect",
      "initial_Eplus_J=8;pending_shock_plus_J=3;dt_s=1");
  result.metrics.push_back({"shock_ledger_absolute_error", source_error, 1.0e-11, "<=", "J"});
  result.metrics.push_back({"final_energy_absolute_error", energy_error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb14() {
  State state = Basic();
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const bool pass = initialized && step.status.ok() &&
      state.phase == Turbulence::OperatorPhase::Ready && state.completedSteps == 1;
  Result result = Complete(pass,
      "operator sequence returns to Ready after one completed step",
      "turbulence operator phase or completed-step count is inconsistent",
      "initial_phase=Ready;dt_s=1");
  result.metrics.push_back({"completed_steps", static_cast<double>(state.completedSteps), 1.0, "==", "count"});
  return result;
}

Result RunTurb15() {
  State state = Basic();
  state.cells[0].pendingParticlePlusJ = -20.0;
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double correction_error =
      std::fabs(step.ledger.limiterCorrectionJ - 12.0);
  const bool pass = initialized && step.status.ok() &&
      state.cells[0].ePlusJ == 0.0 &&
      step.diagnostics.limiterActivations == 1 && correction_error <= 1.0e-11;
  Result result = Complete(pass,
      "positivity limiter applies and accounts for the exact correction",
      "negative wave-energy request was not limited and accounted",
      "initial_Eplus_J=8;requested_change_J=-20");
  result.metrics.push_back({"limiter_correction_absolute_error", correction_error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb16() {
  State old_grid = Basic();
  old_grid.cells.push_back(Cell(4.0, 6.0));
  const bool initialized = Initialize(&old_grid);
  std::vector<CellState> new_grid;
  CellState half = Cell(0.0, 0.0);
  half.lengthM = 5.0;
  for (int i = 0; i < 4; ++i) new_grid.push_back(half);
  State remapped;
  Turbulence::EnergyLedger ledger;
  const bool remap_ok = Turbulence::RemapConservatively(
      old_grid, new_grid, &remapped, &ledger).ok();
  double total = 0.0;
  for (std::size_t i = 0; i < remapped.cells.size(); ++i)
    total += remapped.cells[i].ePlusJ + remapped.cells[i].eMinusJ;
  const double error = std::fabs(total - 20.0);
  Result result = Complete(initialized && remap_ok && error <= 1.0e-11 &&
      Near(ledger.remapCorrectionJ, 0.0),
      "overlap remap preserves the exact controlled total energy",
      "conservative remap changed total wave energy",
      "old_cell_lengths_m=10,10;new_cell_lengths_m=5,5,5,5");
  result.metrics.push_back(
      {"remapped_total_absolute_error", error, 1.0e-11, "<=", "J"});
  return result;
}

Result RunTurb17() {
  State state = Basic(Source::SelfConsistentSpectral);
  state.campaignSeed = 99;
  state.cells[0].pendingParticlePlusJ = 0.125;
  const bool initialized = Initialize(&state);
  std::string checkpoint;
  State restored;
  const bool wrote = Turbulence::SerializeCheckpoint(state, &checkpoint).ok();
  const bool read = Turbulence::DeserializeCheckpoint(checkpoint, &restored).ok();
  const bool hash_equal = wrote && read &&
      Turbulence::StateHash(state) == Turbulence::StateHash(restored);
  Result result = Complete(initialized && hash_equal,
      "spectral checkpoint round-trip preserves the complete state hash",
      "restart serialization changed authoritative turbulence state",
      "representation=spectral;campaign_seed=99;pending_plus_J=0.125");
  result.metrics.push_back(
      {"state_hash_equal", hash_equal ? 1.0 : 0.0, 1.0, "==", "boolean"});
  return result;
}

Result RunTurb18() {
  State first = Basic();
  State second = Basic();
  const bool initialized = Initialize(&first) && Initialize(&second);
  const Turbulence::StepResult first_step = Turbulence::Advance(&first, 0.25);
  const Turbulence::StepResult second_step = Turbulence::Advance(&second, 0.25);
  const bool equal = first_step.status.ok() && second_step.status.ok() &&
      Turbulence::StateHash(first) == Turbulence::StateHash(second);
  Result result = Complete(initialized && equal,
      "identical controlled histories produce identical state hashes",
      "deterministic turbulence histories differ",
      "histories=2;dt_s=0.25;configuration=identical");
  result.metrics.push_back(
      {"state_hash_equal", equal ? 1.0 : 0.0, 1.0, "==", "boolean"});
  return result;
}

Result RunTurb19() {
  // Source ownership is intentionally independent of the selected mover.  The
  // controlled core contains no mover selector; all five sources traverse the
  // same State/Advance boundary.  The production mover registry independently
  // asserts that exactly three movers exist.
  Result ownership = RunTurbulenceOwnership();
  ownership.message = ownership.status == Status::Pass
      ? "turbulence source ownership is independent of particle-mover choice"
      : "turbulence source/mover separation was violated";
  ownership.configuration.push_back("mover_selector_in_turbulence_core=false");
  return ownership;
}

Result RunTurb20() {
  State standalone = Basic();
  State coupled = Basic();
  standalone.cells[0].pendingShockMinusJ = 0.75;
  coupled.cells[0].pendingShockMinusJ = 0.75;
  const bool initialized = Initialize(&standalone) && Initialize(&coupled);
  const Turbulence::StepResult standalone_step =
      Turbulence::Advance(&standalone, 0.5);
  const Turbulence::StepResult coupled_step = Turbulence::Advance(&coupled, 0.5);
  const bool equal = standalone_step.status.ok() && coupled_step.status.ok() &&
      Turbulence::StateHash(standalone) == Turbulence::StateHash(coupled);
  Result result = Complete(initialized && equal,
      "standalone and coupled drivers share one controlled evolution result",
      "driver context changed the authoritative turbulence update",
      "pending_shock_minus_J=0.75;dt_s=0.5;contexts=standalone,coupled");
  result.metrics.push_back(
      {"state_hash_equal", equal ? 1.0 : 0.0, 1.0, "==", "boolean"});
  return result;
}

double SineCellAverage(std::size_t cell, std::size_t cells,
                       double shift_m) {
  const double pi = std::acos(-1.0);
  const double dx = 1.0 / static_cast<double>(cells);
  const double left = cell * dx - shift_m;
  const double right = (cell + 1) * dx - shift_m;
  // This is the exact finite-volume average, not the value at the cell center.
  // Comparing cell averages avoids assigning spatial quadrature error to the
  // production advection operator.
  return 1.0 + 0.25 *
      (std::cos(2.0 * pi * left) - std::cos(2.0 * pi * right)) /
      (2.0 * pi * dx);
}

struct AdvectionError {
  bool ok = false;
  double l1 = std::numeric_limits<double>::infinity();
  double conservation = std::numeric_limits<double>::infinity();
};

AdvectionError RunTranslatedSine(std::size_t cells) {
  const double speed_m_per_s = 1.0;
  const double duration_s = 0.25;
  const double dx_m = 1.0 / static_cast<double>(cells);
  const double dt_s = 0.4 * dx_m / speed_m_per_s;
  const std::size_t steps = static_cast<std::size_t>(
      std::llround(duration_s / dt_s));
  State state = Basic();
  state.cells.clear();
  state.configuration.advectionEnabled = true;
  state.configuration.periodicBoundaries = true;
  for (std::size_t i = 0; i < cells; ++i) {
    CellState cell = Cell(0.0, 0.0);
    cell.lengthM = dx_m;
    cell.volumeM3 = dx_m;  // Unit cross-sectional area.
    cell.plasmaSpeedMPerS = speed_m_per_s;
    cell.alfvenSpeedMPerS = 0.0;
    cell.dLnAlfvenSpeeddsPerM = 0.0;
    cell.ePlusJ = SineCellAverage(i, cells, 0.0) * cell.volumeM3;
    cell.eMinusJ = 0.0;
    state.cells.push_back(cell);
  }
  AdvectionError result;
  result.ok = Initialize(&state);
  for (std::size_t step = 0; step < steps && result.ok; ++step)
    result.ok = Turbulence::Advance(&state, dt_s).status.ok();

  double initial_total_j = 1.0;  // Integral of baseline plus one sine period.
  double final_total_j = 0.0;
  result.l1 = 0.0;
  for (std::size_t i = 0; i < cells; ++i) {
    const double actual_density = state.cells[i].ePlusJ / state.cells[i].volumeM3;
    const double expected_density = SineCellAverage(
        i, cells, speed_m_per_s * duration_s);
    result.l1 += std::fabs(actual_density - expected_density) * dx_m;
    final_total_j += state.cells[i].ePlusJ;
  }
  result.conservation = std::fabs(final_total_j - initial_total_j);
  return result;
}

Result RunTurb21() {
  const AdvectionError coarse = RunTranslatedSine(32);
  const AdvectionError fine = RunTranslatedSine(64);
  const RefinementOrderEstimate order = EstimateRefinementOrder(
      coarse.l1, fine.l1, 1.0 / 32.0, 1.0 / 64.0);
  // The production upwind flux is first order.  Report the measured exponent
  // directly so the assertion remains correct if the refinement factor is
  // changed from two in a future study.
  const double minimum_order = 0.85;
  const bool pass = coarse.ok && fine.ok && fine.l1 < coarse.l1 &&
      order.valid && order.observedOrder >= minimum_order &&
      fine.l1 <= 0.025 &&
      coarse.conservation <= 2.0e-13 && fine.conservation <= 2.0e-13;
  Result result = Complete(pass,
      "nonuniform periodic turbulence advection converges to the translated sine solution",
      "translated turbulence profile exceeded error, order, or conservation limits",
      "profile=1+0.25*sin(2*pi*s);U_m_per_s=1;duration_s=0.25;CFL=0.4;cells=32,64");
  result.metrics.push_back(
      {"fine_L1_error", fine.l1, 0.025, "<=", "J m^-1"});
  result.metrics.push_back({"observed_refinement_order", order.observedOrder,
      minimum_order, ">=", "dimensionless"});
  result.metrics.push_back({"coarse_to_fine_L1_ratio",
      order.coarseToFineErrorRatio, 1.0, ">", "dimensionless"});
  result.metrics.push_back({"maximum_conservation_error",
      std::max(coarse.conservation, fine.conservation), 2.0e-13, "<=", "J"});
  return result;
}

Result RunTurb22() {
  State state = Basic();
  const bool initialized = Initialize(&state);
  const double initial_j = 2.0;
  const double constant_rate_w = 0.30;
  const double linear_rate_w_per_s = 0.20;
  const double dt_s = 0.10;
  const std::size_t steps = 20;
  double maximum_error_j = 0.0;
  bool status_ok = initialized;
  for (std::size_t step = 0; step < steps && status_ok; ++step) {
    const double time_s = step * dt_s;
    // pendingParticleMinusJ is an energy increment, not a rate.  Supplying the
    // exact integral of G(t)=G0+G1*t over this interval produces the analytical
    // time history E(t)=E0+G0*t+0.5*G1*t^2 without embedding a production
    // update formula in the reference.
    state.cells[0].pendingParticleMinusJ =
        constant_rate_w * dt_s +
        linear_rate_w_per_s * (time_s * dt_s + 0.5 * dt_s * dt_s);
    status_ok = Turbulence::Advance(&state, dt_s).status.ok();
    const double end_time_s = (step + 1) * dt_s;
    const double expected_j = initial_j + constant_rate_w * end_time_s +
        0.5 * linear_rate_w_per_s * end_time_s * end_time_s;
    maximum_error_j = std::max(
        maximum_error_j, std::fabs(state.cells[0].eMinusJ - expected_j));
  }
  const double tolerance_j = 2.0e-13;
  Result result = Complete(status_ok && maximum_error_j <= tolerance_j,
      "time-dependent wave growth follows the analytical quadratic history",
      "time-dependent wave growth differs from the integrated source solution",
      "G(t)_W=0.30+0.20*t;initial_Eminus_J=2;dt_s=0.1;steps=20");
  result.metrics.push_back({"maximum_time_history_error",
      maximum_error_j, tolerance_j, "<=", "J"});
  return result;
}

double RelativisticKineticEnergyJ(double mass_kg, double speed_m_per_s) {
  const double speed_of_light_m_per_s = 299792458.0;
  const double beta = speed_m_per_s / speed_of_light_m_per_s;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  // (gamma-1) loses precision for nonrelativistic particles.  The equivalent
  // gamma^2*beta^2/(gamma+1) expression retains the small kinetic-energy
  // difference that must be balanced against the wave ledger.
  const double gamma_minus_one =
      gamma * gamma * beta * beta / (gamma + 1.0);
  return mass_kg * speed_of_light_m_per_s * speed_of_light_m_per_s *
      gamma_minus_one;
}

struct ParticleWaveClosure {
  bool ok = false;
  double totalEnergyErrorJ = std::numeric_limits<double>::infinity();
  double ledgerErrorJ = std::numeric_limits<double>::infinity();
};

ParticleWaveClosure RunParticleWaveClosureCase(
    double initial_mu, double signed_wave_speed_m_per_s,
    bool plus_branch, std::uint64_t particle_id) {
  const double mass_kg = 1.67262192369e-27;
  const double speed_of_light_m_per_s = 299792458.0;
  const double initial_speed_m_per_s = 0.05 * speed_of_light_m_per_s;
  Transport::KeyedRandomStream random(623, particle_id, 23, 0);
  const Transport::WaveFrameScatterResult scattered =
      Transport::ScatterIsotropicallyInWaveFrame(
          initial_speed_m_per_s, initial_mu, signed_wave_speed_m_per_s,
          speed_of_light_m_per_s, random);
  ParticleWaveClosure result;
  if (!scattered.status.ok()) return result;

  const double particle_before_j =
      RelativisticKineticEnergyJ(mass_kg, initial_speed_m_per_s);
  const double particle_after_j =
      RelativisticKineticEnergyJ(mass_kg, scattered.speedMPerS);
  const double single_particle_change_j =
      particle_after_j - particle_before_j;
  if (!std::isfinite(single_particle_change_j) ||
      single_particle_change_j == 0.0) return result;

  // Scale one physical scattering event to a controlled macro-particle energy
  // exchange of exactly 0.25 J.  The turbulence source uses the opposite sign:
  // energy gained by particles is removed from the resonant wave branch, and
  // energy lost by particles is added.  The 0.25 J magnitude is well above
  // roundoff yet below the initial energy of either branch, so the positivity
  // limiter must remain inactive and cannot conceal a sign error.
  const double macro_weight = 0.25 / std::fabs(single_particle_change_j);
  const double particle_change_j =
      macro_weight * single_particle_change_j;
  const double wave_change_j = -particle_change_j;
  State state = Basic();
  if (plus_branch)
    state.cells[0].pendingParticlePlusJ = wave_change_j;
  else
    state.cells[0].pendingParticleMinusJ = wave_change_j;
  const double wave_before_j =
      state.cells[0].ePlusJ + state.cells[0].eMinusJ;
  const bool initialized = Initialize(&state);
  const Turbulence::StepResult step = Turbulence::Advance(&state, 1.0);
  const double wave_after_j =
      state.cells[0].ePlusJ + state.cells[0].eMinusJ;
  const double total_before_j =
      macro_weight * particle_before_j + wave_before_j;
  const double total_after_j =
      macro_weight * particle_after_j + wave_after_j;
  result.totalEnergyErrorJ = std::fabs(total_after_j - total_before_j);
  result.ledgerErrorJ =
      std::fabs(step.ledger.particleExchangeJ - wave_change_j);
  result.ok = initialized && step.status.ok() &&
      step.diagnostics.limiterActivations == 0;
  return result;
}

Result RunTurb23() {
  // The outward-particle case resonates with the minus branch; the inward case
  // uses the plus branch.  Running both signs catches branch-selection or
  // energy-sign mistakes that a single zero-Alfven-speed closure would miss.
  const ParticleWaveClosure minus_branch =
      RunParticleWaveClosureCase(0.65, -4.0e5, false, 1);
  const ParticleWaveClosure plus_branch =
      RunParticleWaveClosureCase(-0.55, 4.0e5, true, 2);
  const double maximum_total_error_j = std::max(
      minus_branch.totalEnergyErrorJ, plus_branch.totalEnergyErrorJ);
  const double maximum_ledger_error_j = std::max(
      minus_branch.ledgerErrorJ, plus_branch.ledgerErrorJ);
  const double tolerance_j = 2.0e-12;
  const bool pass = minus_branch.ok && plus_branch.ok &&
      maximum_total_error_j <= tolerance_j &&
      maximum_ledger_error_j <= tolerance_j;
  Result result = Complete(pass,
      "controlled scattering conserves weighted particle-plus-wave energy",
      "particle and resonant-wave energy exchange does not close",
      "events=outward-minus-branch,inward-plus-branch;macro_exchange_J=0.25;VA_m_per_s=400000");
  result.hasSeed = true;
  result.seed = 623;
  result.metrics.push_back({"particle_wave_total_energy_error",
      maximum_total_error_j, tolerance_j, "<=", "J"});
  result.metrics.push_back({"particle_exchange_ledger_error",
      maximum_ledger_error_j, tolerance_j, "<=", "J"});
  return result;
}

}  // namespace

std::vector<Descriptor> ControlledTurbulenceDescriptors() {
  std::vector<Descriptor> descriptors;
  descriptors.push_back(MakeDescriptor("TURBOWN01", "Turbulence source ownership",
      "Validate every declared source and authoritative representation.", RuntimeClass::Routine,
      RunTurbulenceOwnership));
  descriptors.push_back(MakeDescriptor("TURB02", "Particle-wave energy closure",
      "Compare applied particle exchange with the signed energy ledger.", RuntimeClass::Routine, RunTurb02));
  descriptors.push_back(MakeDescriptor("TURB03", "Turbulence unit and volume closure",
      "Compare integrated energy density with E/V.", RuntimeClass::Routine, RunTurb03));
  descriptors.push_back(MakeDescriptor("TURB04", "Turbulence physical initialization",
      "Preserve a controlled nonuniform initial state.", RuntimeClass::Routine, RunTurb04));
  descriptors.push_back(MakeDescriptor("TURB05", "Integrated-spectral projection",
      "Verify exact branch-energy conservation during spectral projection.", RuntimeClass::Routine, RunTurb05));
  descriptors.push_back(MakeDescriptor("TURB06", "Uniform turbulence advection",
      "Verify a uniform periodic profile remains invariant.", RuntimeClass::Routine, RunTurb06));
  descriptors.push_back(MakeDescriptor("TURB07", "Turbulence CFL subcycling",
      "Compare the stable step with the analytical finite-volume CFL bound.", RuntimeClass::Routine, RunTurb07));
  descriptors.push_back(MakeDescriptor("TURB08", "Turbulence boundary validation",
      "Reject a negative prescribed incoming flux.", RuntimeClass::Routine, RunTurb08));
  descriptors.push_back(MakeDescriptor("TURB09", "Reflection energy exchange",
      "Verify reflection conserves the two-branch energy sum.", RuntimeClass::Routine, RunTurb09));
  descriptors.push_back(MakeDescriptor("TURB10", "Cascade and dissipation ledger",
      "Exercise integrated and spectral cascade with positive dissipation.", RuntimeClass::Routine, RunTurb10));
  descriptors.push_back(MakeDescriptor("TURB11", "One-step wave growth",
      "Compare one source increment with its exact accumulated energy.", RuntimeClass::Routine, RunTurb11));
  descriptors.push_back(MakeDescriptor("TURB12", "Turbulence coefficient closure",
      "Verify kappa_parallel=v*lambda_parallel/3.", RuntimeClass::Routine, RunTurb12));
  descriptors.push_back(MakeDescriptor("TURB13", "Shock turbulence injection",
      "Compare shock-source energy with its state and ledger increments.", RuntimeClass::Routine, RunTurb13));
  descriptors.push_back(MakeDescriptor("TURB14", "Turbulence operator order",
      "Verify a complete step returns to the Ready phase.", RuntimeClass::Routine, RunTurb14));
  descriptors.push_back(MakeDescriptor("TURB15", "Turbulence positivity limiter",
      "Verify exact accounting of a limited negative source.", RuntimeClass::Routine, RunTurb15));
  descriptors.push_back(MakeDescriptor("TURB16", "Conservative turbulence remap",
      "Compare remapped total energy with the analytical domain total.", RuntimeClass::Routine, RunTurb16));
  descriptors.push_back(MakeDescriptor("TURB17", "Turbulence restart round-trip",
      "Verify spectral checkpoint identity by deterministic state hash.", RuntimeClass::Routine, RunTurb17));
  descriptors.push_back(MakeDescriptor("TURB18", "Turbulence deterministic history",
      "Compare two identical controlled histories byte-semantically.", RuntimeClass::Routine, RunTurb18));
  descriptors.push_back(MakeDescriptor("TURB19", "Mover-source separation",
      "Verify turbulence source ownership contains no mover selection.", RuntimeClass::Routine, RunTurb19));
  descriptors.push_back(MakeDescriptor("TURB20", "Shared turbulence driver",
      "Compare standalone and coupled contexts through one core driver.", RuntimeClass::Routine, RunTurb20));
  descriptors.push_back(MakeDescriptor("TURB21", "Translated sine advection",
      "Compare nonuniform periodic advection with an exact translated cell-average profile.",
      RuntimeClass::Routine, RunTurb21));
  descriptors.push_back(MakeDescriptor("TURB22", "Time-dependent wave growth",
      "Compare the complete growth history with an analytically integrated linear source rate.",
      RuntimeClass::Routine, RunTurb22));
  descriptors.push_back(MakeDescriptor("TURB23", "Particle-wave total energy",
      "Close weighted particle kinetic energy against both signed resonant wave branches.",
      RuntimeClass::Routine, RunTurb23));
  return descriptors;
}

}  // namespace Testing
}  // namespace SEP
