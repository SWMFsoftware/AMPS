#include "sep_turbulence_core.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP {
namespace Turbulence {
namespace {

const double kMu0 = 1.25663706212e-6;
Configuration gActiveConfiguration;

Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                  message);
}

bool FiniteNonNegative(double value) {
  return std::isfinite(value) && value >= 0.0;
}

double CellEnergy(const CellState& cell) {
  return cell.ePlusJ + cell.eMinusJ;
}

double TotalEnergy(const State& state) {
  double result = 0.0;
  for (std::size_t i = 0; i < state.cells.size(); ++i)
    result += CellEnergy(state.cells[i]);
  return result;
}

// A stable, branch-major projection is used only when an integrated state is
// promoted to the spectral authority.  The k^(-2/3) energy-per-log-k shape is
// the Kolmogorov inertial-range convention; normalization makes the projection
// exactly conservative apart from floating-point summation roundoff.
void ProjectBranch(double total, const Configuration& config,
                   std::vector<double>* spectrum, std::size_t offset) {
  long double weightSum = 0.0L;
  std::vector<long double> weights(config.spectralBins, 0.0L);
  const long double logMin = std::log(config.spectralKMinPerM);
  const long double dlog =
      (std::log(config.spectralKMaxPerM) - logMin) /
      static_cast<long double>(config.spectralBins);
  for (std::size_t k = 0; k < config.spectralBins; ++k) {
    const long double center = std::exp(logMin + (k + 0.5L) * dlog);
    weights[k] = std::pow(center, -2.0L / 3.0L);
    weightSum += weights[k];
  }
  long double assigned = 0.0L;
  for (std::size_t k = 0; k + 1 < config.spectralBins; ++k) {
    const double value = static_cast<double>(total * weights[k] / weightSum);
    (*spectrum)[offset + k] = value;
    assigned += value;
  }
  // Assigning the remainder to the final bin preserves the input branch total
  // even when a long spectrum accumulates small per-bin rounding errors.
  (*spectrum)[offset + config.spectralBins - 1] =
      std::max(0.0, total - static_cast<double>(assigned));
}

// When a branch-integrated source changes a spectral state, preserve the
// existing shape by a conservative rescale. A previously empty branch receives
// the same normalized Kolmogorov initialization used at state construction.
void MatchSpectrumBranchToIntegrated(CellState* cell,
                                     const Configuration& config,
                                     bool plus) {
  const std::size_t offset = plus ? 0 : config.spectralBins;
  const double requested = plus ? cell->ePlusJ : cell->eMinusJ;
  long double oldTotal = 0.0L;
  for (std::size_t k = 0; k < config.spectralBins; ++k)
    oldTotal += cell->spectralEnergyJ[offset + k];
  if (oldTotal > 0.0L) {
    const double factor = requested / static_cast<double>(oldTotal);
    long double assigned = 0.0L;
    for (std::size_t k = 0; k + 1 < config.spectralBins; ++k) {
      cell->spectralEnergyJ[offset + k] *= factor;
      assigned += cell->spectralEnergyJ[offset + k];
    }
    cell->spectralEnergyJ[offset + config.spectralBins - 1] =
        std::max(0.0, requested - static_cast<double>(assigned));
  }
  else {
    ProjectBranch(requested, config, &cell->spectralEnergyJ, offset);
  }
}

void MatchSpectrumToIntegrated(State* state) {
  if (state->configuration.representation != Representation::Spectral) return;
  for (std::size_t i = 0; i < state->cells.size(); ++i) {
    MatchSpectrumBranchToIntegrated(&state->cells[i], state->configuration, true);
    MatchSpectrumBranchToIntegrated(&state->cells[i], state->configuration, false);
  }
}

void AddLedger(EnergyLedger* target, const EnergyLedger& value) {
  target->initialEnergyJ += value.initialEnergyJ;
  target->innerBoundaryJ += value.innerBoundaryJ;
  target->outerBoundaryJ += value.outerBoundaryJ;
  target->shockSourceJ += value.shockSourceJ;
  target->reflectionTransferJ += value.reflectionTransferJ;
  target->cascadeTransferJ += value.cascadeTransferJ;
  target->physicalDissipationJ += value.physicalDissipationJ;
  target->particleExchangeJ += value.particleExchangeJ;
  target->limiterCorrectionJ += value.limiterCorrectionJ;
  target->rejectedSourceJ += value.rejectedSourceJ;
  target->remapCorrectionJ += value.remapCorrectionJ;
  target->finalEnergyJ += value.finalEnergyJ;
  target->closureResidualJ += value.closureResidualJ;
}

std::vector<double> LedgerValues(const EnergyLedger& value) {
  std::vector<double> result;
  result.push_back(value.initialEnergyJ);
  result.push_back(value.innerBoundaryJ);
  result.push_back(value.outerBoundaryJ);
  result.push_back(value.shockSourceJ);
  result.push_back(value.reflectionTransferJ);
  result.push_back(value.cascadeTransferJ);
  result.push_back(value.physicalDissipationJ);
  result.push_back(value.particleExchangeJ);
  result.push_back(value.limiterCorrectionJ);
  result.push_back(value.rejectedSourceJ);
  result.push_back(value.remapCorrectionJ);
  result.push_back(value.finalEnergyJ);
  result.push_back(value.closureResidualJ);
  return result;
}

bool ReadLedger(std::istream& in, EnergyLedger* value) {
  return static_cast<bool>(
      in >> value->initialEnergyJ >> value->innerBoundaryJ >>
      value->outerBoundaryJ >> value->shockSourceJ >>
      value->reflectionTransferJ >> value->cascadeTransferJ >>
      value->physicalDissipationJ >> value->particleExchangeJ >>
      value->limiterCorrectionJ >> value->rejectedSourceJ >>
      value->remapCorrectionJ >>
      value->finalEnergyJ >> value->closureResidualJ);
}

bool IsImmutable(Source source) {
  return source == Source::Prescribed || source == Source::SwmfReadOnly;
}

// Apply an externally accumulated source while preventing negative wave
// energy.  The ledger records the *applied* exchange.  The difference from the
// request is separately exposed as a limiter correction and activation count.
double ApplyLimitedSource(const char* operatorName, double requested,
                          double* energy, EnergyLedger* ledger,
                          Diagnostics* diagnostics) {
  OperatorEvent event;
  event.operatorName = operatorName;
  event.requestedJ = requested;
  if (requested == 0.0) {
    // Exact zero is a physical no-op, not a missing or rejected contribution.
    // Recording it distinguishes a quiet source from an adapter that failed to
    // enqueue a record at all.
    event.disposition = UpdateDisposition::PhysicalZero;
    event.reason = "exact zero source";
    ++diagnostics->physicalZeroEvents;
    diagnostics->events.push_back(event);
    return 0.0;
  }
  const double applied = std::max(requested, -*energy);
  if (applied != requested) {
    ledger->limiterCorrectionJ += applied - requested;
    ledger->rejectedSourceJ += requested - applied;
    ++diagnostics->limiterActivations;
    event.disposition = UpdateDisposition::Corrected;
    event.reason = "positivity-preserving source limit";
  }
  event.appliedJ = applied;
  diagnostics->events.push_back(event);
  *energy += applied;
  return applied;
}

// Returns energy entering the domain.  For the outer boundary the caller
// supplies the inward characteristic sign, so this helper has no hidden
// orientation convention.
double IncomingBoundaryEnergy(const BoundaryCondition& boundary,
                              double reservoirEnergy, double speed,
                              double area, double dt) {
  switch (boundary.policy) {
    case BoundaryPolicy::SpecifiedIncomingEnergy:
      return boundary.value;
    case BoundaryPolicy::SpecifiedIncomingFlux:
      return boundary.value * dt;
    case BoundaryPolicy::FixedReservoir:
      return reservoirEnergy * std::fabs(speed) * area * dt;
    case BoundaryPolicy::TransparentOutflow:
      return 0.0;
  }
  return 0.0;
}

// First-order finite-volume advection is intentionally modest but auditable:
// every face flux is evaluated once, then applied with equal and opposite signs
// to its two neighboring cells.  Periodic tests therefore close to roundoff.
void AdvectBranch(State* state, bool plus, double dt, EnergyLedger* ledger,
                  Diagnostics* diagnostics) {
  const std::size_t n = state->cells.size();
  std::vector<double> delta(n, 0.0);
  const Configuration& config = state->configuration;

  for (std::size_t face = 0; face + 1 < n; ++face) {
    ++diagnostics->advectionFaceUpdates;
    const CellState& left = state->cells[face];
    const CellState& right = state->cells[face + 1];
    const double speedLeft = left.plasmaSpeedMPerS +
        (plus ? left.alfvenSpeedMPerS : -left.alfvenSpeedMPerS);
    const double speedRight = right.plasmaSpeedMPerS +
        (plus ? right.alfvenSpeedMPerS : -right.alfvenSpeedMPerS);
    const double speed = 0.5 * (speedLeft + speedRight);
    const CellState& upwind = speed >= 0.0 ? left : right;
    const double upwindEnergy = plus ? upwind.ePlusJ : upwind.eMinusJ;
    const double density = upwindEnergy / upwind.volumeM3;
    const double areaLeft = left.volumeM3 / left.lengthM;
    const double areaRight = right.volumeM3 / right.lengthM;
    const double area = 0.5 * (areaLeft + areaRight);
    const double transfer = speed * density * area * dt;
    delta[face] -= transfer;
    delta[face + 1] += transfer;
  }

  if (config.periodicBoundaries && n > 1) {
    ++diagnostics->advectionFaceUpdates;
    const CellState& left = state->cells[n - 1];
    const CellState& right = state->cells[0];
    const double speed = 0.5 *
        (left.plasmaSpeedMPerS + right.plasmaSpeedMPerS +
         (plus ? left.alfvenSpeedMPerS + right.alfvenSpeedMPerS
               : -left.alfvenSpeedMPerS - right.alfvenSpeedMPerS));
    const CellState& upwind = speed >= 0.0 ? left : right;
    const double energy = plus ? upwind.ePlusJ : upwind.eMinusJ;
    const double area = 0.5 * (left.volumeM3 / left.lengthM +
                               right.volumeM3 / right.lengthM);
    const double transfer = speed * energy / upwind.volumeM3 * area * dt;
    delta[n - 1] -= transfer;
    delta[0] += transfer;
  }
  else {
    diagnostics->advectionFaceUpdates += 2;
    const CellState& first = state->cells.front();
    const CellState& last = state->cells.back();
    const double innerSpeed = first.plasmaSpeedMPerS +
        (plus ? first.alfvenSpeedMPerS : -first.alfvenSpeedMPerS);
    const double outerSpeed = last.plasmaSpeedMPerS +
        (plus ? last.alfvenSpeedMPerS : -last.alfvenSpeedMPerS);
    const double firstEnergy = plus ? first.ePlusJ : first.eMinusJ;
    const double lastEnergy = plus ? last.ePlusJ : last.eMinusJ;

    if (innerSpeed > 0.0) {
      const double incoming = IncomingBoundaryEnergy(
          config.innerBoundary, firstEnergy / first.volumeM3, innerSpeed,
          first.volumeM3 / first.lengthM, dt);
      delta[0] += incoming;
      ledger->innerBoundaryJ += incoming;
    }
    else {
      const double outgoing = innerSpeed * firstEnergy / first.volumeM3 *
                              (first.volumeM3 / first.lengthM) * dt;
      delta[0] += outgoing;
      ledger->innerBoundaryJ += outgoing;
    }

    if (outerSpeed < 0.0) {
      const double incoming = IncomingBoundaryEnergy(
          config.outerBoundary, lastEnergy / last.volumeM3, outerSpeed,
          last.volumeM3 / last.lengthM, dt);
      delta[n - 1] += incoming;
      ledger->outerBoundaryJ += incoming;
    }
    else {
      const double outgoing = outerSpeed * lastEnergy / last.volumeM3 *
                              (last.volumeM3 / last.lengthM) * dt;
      delta[n - 1] -= outgoing;
      ledger->outerBoundaryJ -= outgoing;
    }
  }

  for (std::size_t i = 0; i < n; ++i) {
    double* energy = plus ? &state->cells[i].ePlusJ : &state->cells[i].eMinusJ;
    *energy += delta[i];
    // The CFL calculation guarantees positivity for the first-order scheme;
    // this branch catches only accumulated roundoff or malformed host input.
    if (*energy < 0.0) {
      OperatorEvent event;
      event.operatorName = plus ? "advection-plus" : "advection-minus";
      event.disposition = UpdateDisposition::Corrected;
      event.requestedJ = *energy;
      event.appliedJ = 0.0;
      event.reason = "roundoff-only conservative positivity correction";
      diagnostics->events.push_back(event);
      ++diagnostics->limiterActivations;
      ledger->limiterCorrectionJ -= *energy;
      *energy = 0.0;
    }
  }
}

void WriteBoundary(std::ostream& out, const BoundaryCondition& boundary) {
  out << static_cast<int>(boundary.policy) << ' '
      << std::setprecision(17) << boundary.value << '\n';
}

std::string HexEncode(const std::string& value) {
  static const char digits[] = "0123456789abcdef";
  // A leading marker keeps the encoded token non-empty, which allows the
  // whitespace-delimited checkpoint grammar to represent an empty string.
  std::string encoded("x");
  encoded.reserve(1 + 2 * value.size());
  for (std::size_t i = 0; i < value.size(); ++i) {
    const unsigned char byte = static_cast<unsigned char>(value[i]);
    encoded.push_back(digits[byte >> 4]);
    encoded.push_back(digits[byte & 0x0fU]);
  }
  return encoded;
}

bool HexDecode(const std::string& encoded, std::string* value) {
  if (!value || encoded.empty() || encoded[0] != 'x' ||
      encoded.size() % 2 != 1) return false;
  value->clear();
  value->reserve((encoded.size() - 1) / 2);
  for (std::size_t i = 1; i < encoded.size(); i += 2) {
    const std::size_t high = std::string("0123456789abcdef").find(encoded[i]);
    const std::size_t low = std::string("0123456789abcdef").find(encoded[i + 1]);
    if (high == std::string::npos || low == std::string::npos) return false;
    value->push_back(static_cast<char>((high << 4) | low));
  }
  return true;
}

bool ReadBoundary(std::istream& in, BoundaryCondition* boundary) {
  int policy = -1;
  if (!(in >> policy >> boundary->value) || policy < 0 || policy > 3)
    return false;
  boundary->policy = static_cast<BoundaryPolicy>(policy);
  return true;
}

}  // namespace

const char* SourceName(Source value) {
  switch (value) {
    case Source::Prescribed: return "prescribed";
    case Source::SelfConsistentIntegrated: return "self-consistent-integrated";
    case Source::SelfConsistentSpectral: return "self-consistent-spectral";
    case Source::SwmfReadOnly: return "swmf-read-only";
    case Source::SwmfInitialThenEvolveLocal: return "swmf-initial-then-local";
  }
  return "unknown";
}

const char* RepresentationName(Representation value) {
  return value == Representation::Integrated ? "integrated" : "spectral";
}

const char* BoundaryPolicyName(BoundaryPolicy value) {
  switch (value) {
    case BoundaryPolicy::SpecifiedIncomingEnergy: return "specified-incoming-energy";
    case BoundaryPolicy::SpecifiedIncomingFlux: return "specified-incoming-flux";
    case BoundaryPolicy::TransparentOutflow: return "transparent-outflow";
    case BoundaryPolicy::FixedReservoir: return "fixed-reservoir";
  }
  return "unknown";
}

const char* CouplingPolicyName(CouplingPolicy value) {
  return value == CouplingPolicy::Disabled ? "disabled" : "streaming-energy-exchange";
}

bool ParseSource(const std::string& text, Source* value) {
  if (!value) return false;
  if (text == "prescribed") *value = Source::Prescribed;
  else if (text == "self-consistent-integrated" || text == "integrated")
    *value = Source::SelfConsistentIntegrated;
  else if (text == "self-consistent-spectral" || text == "spectral")
    *value = Source::SelfConsistentSpectral;
  else if (text == "swmf-read-only") *value = Source::SwmfReadOnly;
  else if (text == "swmf-initial-then-local")
    *value = Source::SwmfInitialThenEvolveLocal;
  else return false;
  return true;
}

bool ParseBoundaryPolicy(const std::string& text, BoundaryPolicy* value) {
  if (!value) return false;
  if (text == "specified-incoming-energy")
    *value = BoundaryPolicy::SpecifiedIncomingEnergy;
  else if (text == "specified-incoming-flux")
    *value = BoundaryPolicy::SpecifiedIncomingFlux;
  else if (text == "transparent-outflow")
    *value = BoundaryPolicy::TransparentOutflow;
  else if (text == "fixed-reservoir") *value = BoundaryPolicy::FixedReservoir;
  else return false;
  return true;
}

bool ParseCouplingPolicy(const std::string& text, CouplingPolicy* value) {
  if (!value) return false;
  if (text == "disabled" || text == "off") *value = CouplingPolicy::Disabled;
  else if (text == "streaming-energy-exchange" || text == "streaming")
    *value = CouplingPolicy::StreamingEnergyExchange;
  else return false;
  return true;
}

Transport::Status ValidateConfiguration(const Configuration& c) {
  if (c.source == Source::SelfConsistentIntegrated &&
      c.representation != Representation::Integrated)
    return Error("self-consistent-integrated requires integrated authority");
  if (c.source == Source::SelfConsistentSpectral &&
      c.representation != Representation::Spectral)
    return Error("self-consistent-spectral requires spectral authority");
  if (!(c.cflSafety > 0.0 && c.cflSafety <= 1.0) ||
      !std::isfinite(c.cflSafety))
    return Error("turbulence CFL safety must be finite and in (0,1]");
  if (!(c.operatorAccuracySafety > 0.0 && c.operatorAccuracySafety <= 1.0) ||
      !(c.maximumSourceFraction > 0.0 && c.maximumSourceFraction <= 1.0) ||
      !(c.maximumCascadeFraction > 0.0 && c.maximumCascadeFraction <= 1.0) ||
      !(c.minimumSubstepS > 0.0) || !std::isfinite(c.minimumSubstepS) ||
      c.maximumSubsteps == 0)
    return Error("turbulence operator limits must be finite, positive, and bounded");
  if (!FiniteNonNegative(c.reflectionCoefficient) ||
      !FiniteNonNegative(c.cascadeCoefficient) ||
      !(c.perpendicularCorrelationLengthM > 0.0) ||
      !std::isfinite(c.perpendicularCorrelationLengthM))
    return Error("turbulence coefficients and correlation length are invalid");
  if (!(c.electronHeatingFraction >= 0.0 && c.electronHeatingFraction <= 1.0) ||
      !std::isfinite(c.electronHeatingFraction))
    return Error("electron heating fraction must be in [0,1]");
  if (!(c.spectralKMinPerM > 0.0) ||
      !(c.spectralKMaxPerM > c.spectralKMinPerM) || c.spectralBins == 0)
    return Error("spectral grid requires 0 < k_min < k_max and at least one bin");
  if (!FiniteNonNegative(c.conservationRelativeTolerance))
    return Error("conservation tolerance must be finite and non-negative");
  if ((c.innerBoundary.policy != BoundaryPolicy::TransparentOutflow &&
       !FiniteNonNegative(c.innerBoundary.value)) ||
      (c.outerBoundary.policy != BoundaryPolicy::TransparentOutflow &&
       !FiniteNonNegative(c.outerBoundary.value)))
    return Error("boundary values must be finite and non-negative");
  return Transport::Status::Ok();
}

Transport::Status SetActiveConfiguration(const Configuration& configuration) {
  const Transport::Status status = ValidateConfiguration(configuration);
  if (!status.ok()) return status;
  gActiveConfiguration = configuration;
  return Transport::Status::Ok();
}

const Configuration& ActiveConfiguration() {
  return gActiveConfiguration;
}

bool EvolvesLocally(Source source) {
  return source == Source::SelfConsistentIntegrated ||
         source == Source::SelfConsistentSpectral ||
         source == Source::SwmfInitialThenEvolveLocal;
}

Transport::Status InitializeState(State* state) {
  if (!state) return Error("turbulence state pointer is null");
  Transport::Status status = ValidateConfiguration(state->configuration);
  if (!status.ok()) return status;
  if (state->cells.empty()) return Error("turbulence state has no cells");
  const std::size_t bins = state->configuration.spectralBins;
  for (std::size_t i = 0; i < state->cells.size(); ++i) {
    CellState& cell = state->cells[i];
    if (!(cell.lengthM > 0.0) || !(cell.volumeM3 > 0.0) ||
        !std::isfinite(cell.lengthM) || !std::isfinite(cell.volumeM3) ||
        !FiniteNonNegative(cell.ePlusJ) || !FiniteNonNegative(cell.eMinusJ))
      return Error("cell geometry and branch energies must be finite and positive/non-negative");
    if (state->configuration.representation == Representation::Spectral) {
      if (cell.spectralEnergyJ.empty()) {
        cell.spectralEnergyJ.assign(2 * bins, 0.0);
        ProjectBranch(cell.ePlusJ, state->configuration,
                      &cell.spectralEnergyJ, 0);
        ProjectBranch(cell.eMinusJ, state->configuration,
                      &cell.spectralEnergyJ, bins);
      }
      if (cell.spectralEnergyJ.size() != 2 * bins)
        return Error("spectral cell storage does not match configured bin count");
      for (std::size_t k = 0; k < 2 * bins; ++k)
        if (!FiniteNonNegative(cell.spectralEnergyJ[k]))
          return Error("spectral energy must be finite and non-negative");
    }
    else if (!cell.spectralEnergyJ.empty()) {
      return Error("integrated authority cannot retain a competing spectral owner");
    }
  }
  status = SynchronizeDerivedRepresentation(state);
  if (!status.ok()) return status;
  state->phase = OperatorPhase::Ready;
  return Transport::Status::Ok();
}

Transport::Status HandoffImportedState(State* state, double epochS,
                                       const std::string& checksum) {
  if (!state || state->configuration.source != Source::SwmfInitialThenEvolveLocal)
    return Error("SWMF handoff requires swmf-initial-then-local source");
  if (state->handoffCompleted)
    return Error("SWMF handoff is a one-time ownership transfer");
  if (!std::isfinite(epochS) || checksum.empty())
    return Error("SWMF handoff requires a finite epoch and non-empty checksum");
  Transport::Status status = InitializeState(state);
  if (!status.ok()) return status;
  state->handoffCompleted = true;
  state->handoffEpochS = epochS;
  state->sourceChecksum = checksum;
  return Transport::Status::Ok();
}

Transport::Status SynchronizeDerivedRepresentation(State* state) {
  if (!state) return Error("turbulence state pointer is null");
  if (state->configuration.representation == Representation::Integrated)
    return Transport::Status::Ok();
  const std::size_t bins = state->configuration.spectralBins;
  for (std::size_t i = 0; i < state->cells.size(); ++i) {
    CellState& cell = state->cells[i];
    if (cell.spectralEnergyJ.size() != 2 * bins)
      return Error("spectral storage length changed outside the authority");
    long double plus = 0.0L, minus = 0.0L;
    for (std::size_t k = 0; k < bins; ++k) {
      plus += cell.spectralEnergyJ[k];
      minus += cell.spectralEnergyJ[bins + k];
    }
    cell.ePlusJ = static_cast<double>(plus);
    cell.eMinusJ = static_cast<double>(minus);
  }
  return Transport::Status::Ok();
}

DerivedCell DeriveCell(const CellState& cell) {
  DerivedCell result;
  if (!(cell.volumeM3 > 0.0)) return result;
  result.wPlusJPerM3 = cell.ePlusJ / cell.volumeM3;
  result.wMinusJPerM3 = cell.eMinusJ / cell.volumeM3;
  result.deltaB2T2 = 2.0 * kMu0 *
      (result.wPlusJPerM3 + result.wMinusJPerM3);
  const double total = cell.ePlusJ + cell.eMinusJ;
  result.crossHelicity = total > 0.0 ?
      (cell.ePlusJ - cell.eMinusJ) / total : 0.0;
  return result;
}

Transport::ScalarResult MaximumStableAdvectionStepS(const State& state) {
  Transport::ScalarResult result;
  result.status = ValidateConfiguration(state.configuration);
  if (!result.status.ok()) return result;
  result.value = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const CellState& cell = state.cells[i];
    const double maxSpeed = std::max(
        std::fabs(cell.plasmaSpeedMPerS + cell.alfvenSpeedMPerS),
        std::fabs(cell.plasmaSpeedMPerS - cell.alfvenSpeedMPerS));
    if (maxSpeed > 0.0)
      result.value = std::min(result.value,
          state.configuration.cflSafety * cell.lengthM / maxSpeed);
  }
  return result;
}

AdvancePlan PlanAdvance(const State& state, double dtS) {
  AdvancePlan plan;
  plan.requestedStepS = dtS;
  plan.status = ValidateConfiguration(state.configuration);
  if (!plan.status.ok()) return plan;
  if (!(dtS > 0.0) || !std::isfinite(dtS) || state.cells.empty()) {
    plan.status = Error("turbulence plan requires a finite positive step and cells");
    return plan;
  }

  const double infinity = std::numeric_limits<double>::infinity();
  plan.advectionLimitS = infinity;
  plan.sourceLimitS = infinity;
  plan.reflectionLimitS = infinity;
  plan.cascadeLimitS = infinity;

  if (state.configuration.advectionEnabled) {
    const Transport::ScalarResult advection = MaximumStableAdvectionStepS(state);
    if (!advection.status.ok() || !(advection.value > 0.0)) {
      plan.status = advection.status.ok()
          ? Transport::Status::Error(Transport::StatusCode::StepUnderflow,
                                     "advection has no positive stable step")
          : advection.status;
      return plan;
    }
    plan.advectionLimitS = advection.value;
  }

  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const CellState& cell = state.cells[i];
    const double energy = CellEnergy(cell);
    const double sourceMagnitude =
        std::fabs(cell.pendingParticlePlusJ) +
        std::fabs(cell.pendingParticleMinusJ) +
        std::fabs(cell.pendingShockPlusJ) +
        std::fabs(cell.pendingShockMinusJ);
    if (!std::isfinite(sourceMagnitude)) {
      plan.status = Error("nonfinite pending turbulence source");
      return plan;
    }
    if (sourceMagnitude > 0.0 && energy > 0.0) {
      // Pending values are amounts over dtS.  Limiting their fractional change
      // defines an equivalent source timescale without inventing a fallback
      // rate when the host has supplied an empty cell.
      plan.sourceLimitS = std::min(plan.sourceLimitS,
          dtS * state.configuration.maximumSourceFraction * energy /
          sourceMagnitude);
    }
    if (state.configuration.reflectionEnabled) {
      const double rate = 0.5 * state.configuration.reflectionCoefficient *
          std::fabs(cell.alfvenSpeedMPerS * cell.dLnAlfvenSpeeddsPerM);
      if (!std::isfinite(rate)) {
        plan.status = Error("nonfinite turbulence reflection rate");
        return plan;
      }
      if (rate > 0.0)
        plan.reflectionLimitS = std::min(plan.reflectionLimitS,
            state.configuration.operatorAccuracySafety / rate);
    }
    if (state.configuration.cascadeEnabled && energy > 0.0) {
      if (!(cell.volumeM3 > 0.0) || !(cell.massDensityKgPerM3 > 0.0) ||
          !std::isfinite(cell.volumeM3) ||
          !std::isfinite(cell.massDensityKgPerM3)) {
        plan.status = Error("cascade timescale requires positive finite volume and density");
        return plan;
      }
      const double rate = state.configuration.cascadeCoefficient *
          std::sqrt(energy / (cell.volumeM3 * cell.massDensityKgPerM3)) /
          state.configuration.perpendicularCorrelationLengthM;
      if (!std::isfinite(rate)) {
        plan.status = Error("nonfinite turbulence cascade rate");
        return plan;
      }
      if (rate > 0.0)
        plan.cascadeLimitS = std::min(plan.cascadeLimitS,
            state.configuration.maximumCascadeFraction / rate);
    }
  }

  const double limits[] = {plan.advectionLimitS, plan.sourceLimitS,
                           plan.reflectionLimitS, plan.cascadeLimitS};
  const char* names[] = {"advection", "source", "reflection", "cascade"};
  double selected = dtS;
  plan.limitingOperator = "requested-step";
  for (std::size_t i = 0; i < 4; ++i) {
    if (limits[i] < selected) {
      selected = limits[i];
      plan.limitingOperator = names[i];
    }
  }
  if (!(selected > 0.0) || !std::isfinite(selected)) selected = dtS;
  if (selected < state.configuration.minimumSubstepS) {
    plan.status = Transport::Status::Error(
        Transport::StatusCode::StepUnderflow,
        "combined turbulence operator limit is below minimumSubstepS");
    return plan;
  }
  const long double count = std::ceil(
      static_cast<long double>(dtS) / static_cast<long double>(selected));
  if (!(count >= 1.0L) ||
      count > static_cast<long double>(state.configuration.maximumSubsteps)) {
    plan.status = Transport::Status::Error(
        Transport::StatusCode::StepUnderflow,
        "combined turbulence operator plan exceeds maximumSubsteps");
    return plan;
  }
  plan.substeps = static_cast<std::uint64_t>(count);
  plan.selectedSubstepS = dtS / static_cast<double>(plan.substeps);
  plan.status = Transport::Status::Ok();
  return plan;
}

CoefficientView DeriveCoefficients(const State& state, std::size_t cellIndex,
                                   double speed, double mu, double charge,
                                   double mass) {
  CoefficientView result;
  if (cellIndex >= state.cells.size() || !(speed > 0.0) ||
      !(mass > 0.0) || charge == 0.0 || !std::isfinite(mu)) {
    result.status = Error("coefficient query contains an invalid particle or cell");
    return result;
  }
  const CellState& cell = state.cells[cellIndex];
  const DerivedCell derived = DeriveCell(cell);
  if (!(derived.deltaB2T2 > 0.0) || !(cell.magneticFieldT > 0.0)) {
    result.status = Transport::Status::Error(
        Transport::StatusCode::InvalidCoefficient,
        "wave amplitude and magnetic field must be positive");
    return result;
  }
  const double absMu = std::fabs(mu);
  if (absMu == 0.0) {
    result.status = Transport::Status::Error(
        Transport::StatusCode::OutOfDomain,
        "mu=0 has no finite resonant wave number; no hidden pitch floor is applied");
    return result;
  }
  result.resonantWaveNumberPerM =
      std::fabs(charge) * cell.magneticFieldT / (mass * speed * absMu);
  result.resonantBranch = mu >= 0.0 ? 1 : 0;
  if (state.configuration.representation == Representation::Spectral) {
    const double fraction = std::log(result.resonantWaveNumberPerM /
                                     state.configuration.spectralKMinPerM) /
        std::log(state.configuration.spectralKMaxPerM /
                 state.configuration.spectralKMinPerM);
    if (!(fraction >= 0.0 && fraction < 1.0)) {
      result.status = Transport::Status::Error(
          Transport::StatusCode::OutOfDomain,
          "resonant wave number lies outside the spectral grid");
      return result;
    }
    result.resonantBin = std::min(state.configuration.spectralBins - 1,
        static_cast<std::size_t>(fraction * state.configuration.spectralBins));
  }
  result.lambdaParallelM = state.configuration.perpendicularCorrelationLengthM *
      cell.magneticFieldT * cell.magneticFieldT / derived.deltaB2T2;
  result.dMuMuPerS = 0.5 * speed / result.lambdaParallelM *
                     std::max(0.0, 1.0 - mu * mu);
  result.kappaParallelM2PerS = speed * result.lambdaParallelM / 3.0;
  result.status = Transport::Status::Ok();
  return result;
}

StepResult Advance(State* state, double dtS) {
  StepResult result;
  if (!state || !(dtS > 0.0) || !std::isfinite(dtS)) {
    result.status = Error("turbulence advance requires a finite positive timestep");
    return result;
  }
  result.status = ValidateConfiguration(state->configuration);
  if (!result.status.ok()) return result;
  result.status = SynchronizeDerivedRepresentation(state);
  if (!result.status.ok()) return result;
  if (state->configuration.source == Source::SwmfInitialThenEvolveLocal &&
      !state->handoffCompleted) {
    result.status = Error("local evolution cannot start before the SWMF handoff");
    return result;
  }

  const AdvancePlan plan = PlanAdvance(*state, dtS);
  if (!plan.status.ok()) {
    ++result.diagnostics.rejectedUpdates;
    result.status = plan.status;
    return result;
  }
  const std::uint64_t operatorSubsteps =
      std::max<std::uint64_t>(UINT64_C(1), plan.substeps);
  result.diagnostics.subcycles = operatorSubsteps;

  result.ledger.initialEnergyJ = TotalEnergy(*state);
  if (IsImmutable(state->configuration.source)) {
    // Prescribed and read-only SWMF states are immutable by contract.  Pending
    // inputs are left untouched so a caller cannot mistake a discarded source
    // for an applied conservative exchange.
    result.ledger.finalEnergyJ = result.ledger.initialEnergyJ;
    result.status = Transport::Status::Ok();
    return result;
  }

  for (std::size_t i = 0; i < state->cells.size(); ++i) {
    CellState& cell = state->cells[i];
    // The additive source map is exact over the host interval, so it does not
    // need an Euler stability substep.  Its diagnosed timescale still controls
    // the following noncommuting operators, reducing Lie-splitting error when a
    // strong source changes the wave state seen by reflection or cascade.
    result.ledger.shockSourceJ += ApplyLimitedSource(
        "shock-source-plus", cell.pendingShockPlusJ, &cell.ePlusJ,
        &result.ledger, &result.diagnostics);
    result.ledger.shockSourceJ += ApplyLimitedSource(
        "shock-source-minus", cell.pendingShockMinusJ, &cell.eMinusJ,
        &result.ledger, &result.diagnostics);
    cell.pendingShockPlusJ = cell.pendingShockMinusJ = 0.0;
  }
  result.diagnostics.sourceSubcycles = 1;
  MatchSpectrumToIntegrated(state);
  state->phase = OperatorPhase::SourcesApplied;

  if (state->configuration.coupling == CouplingPolicy::StreamingEnergyExchange) {
    for (std::size_t i = 0; i < state->cells.size(); ++i) {
      CellState& cell = state->cells[i];
      result.ledger.particleExchangeJ += ApplyLimitedSource(
          "particle-coupling-plus", cell.pendingParticlePlusJ, &cell.ePlusJ,
          &result.ledger, &result.diagnostics);
      result.ledger.particleExchangeJ += ApplyLimitedSource(
          "particle-coupling-minus", cell.pendingParticleMinusJ, &cell.eMinusJ,
          &result.ledger, &result.diagnostics);
      cell.pendingParticlePlusJ = cell.pendingParticleMinusJ = 0.0;
    }
  }
  MatchSpectrumToIntegrated(state);
  state->phase = OperatorPhase::CouplingApplied;

  if (state->configuration.advectionEnabled) {
    Transport::ScalarResult stable = MaximumStableAdvectionStepS(*state);
    if (!stable.status.ok() || !(stable.value > 0.0)) {
      ++result.diagnostics.rejectedUpdates;
      result.status = Transport::Status::Error(
          Transport::StatusCode::StepUnderflow,
          "no finite positive turbulence advection stability limit");
      return result;
    }
    const std::uint64_t advectionCycles = std::isfinite(stable.value) ?
        static_cast<std::uint64_t>(std::ceil(dtS / stable.value)) : 1;
    const std::uint64_t cycles = std::max(operatorSubsteps, advectionCycles);
    const double substep = dtS / static_cast<double>(std::max<std::uint64_t>(1, cycles));
    const BoundaryCondition advanceInner = state->configuration.innerBoundary;
    const BoundaryCondition advanceOuter = state->configuration.outerBoundary;
    for (std::uint64_t cycle = 0; cycle < std::max<std::uint64_t>(1, cycles); ++cycle) {
      // Specified energy is an amount for the complete driver update, whereas
      // specified flux is a rate multiplied by each substep in the face helper.
      if (advanceInner.policy == BoundaryPolicy::SpecifiedIncomingEnergy)
        state->configuration.innerBoundary.value =
            advanceInner.value * substep / dtS;
      if (advanceOuter.policy == BoundaryPolicy::SpecifiedIncomingEnergy)
        state->configuration.outerBoundary.value =
            advanceOuter.value * substep / dtS;
      if (state->configuration.representation == Representation::Integrated) {
        AdvectBranch(state, true, substep, &result.ledger,
                     &result.diagnostics);
        AdvectBranch(state, false, substep, &result.ledger,
                     &result.diagnostics);
      }
      else {
        // Every logarithmic bin is transported independently. The compact
        // branch totals are temporary work registers, never a second owner.
        const std::size_t bins = state->configuration.spectralBins;
        const BoundaryCondition savedInner = state->configuration.innerBoundary;
        const BoundaryCondition savedOuter = state->configuration.outerBoundary;
        std::vector<double> innerBins(2 * bins, 0.0), outerBins(2 * bins, 0.0);
        if (savedInner.policy == BoundaryPolicy::SpecifiedIncomingEnergy ||
            savedInner.policy == BoundaryPolicy::SpecifiedIncomingFlux)
          ProjectBranch(savedInner.value, state->configuration, &innerBins, 0);
        if (savedOuter.policy == BoundaryPolicy::SpecifiedIncomingEnergy ||
            savedOuter.policy == BoundaryPolicy::SpecifiedIncomingFlux)
          ProjectBranch(savedOuter.value, state->configuration, &outerBins, 0);
        for (std::size_t k = 0; k < bins; ++k) {
          // A configured total incoming energy/flux is distributed once across
          // log-k. Fixed-reservoir input already obtains the current bin's
          // density from the temporary work register and needs no rescaling.
          if (savedInner.policy == BoundaryPolicy::SpecifiedIncomingEnergy ||
              savedInner.policy == BoundaryPolicy::SpecifiedIncomingFlux)
            state->configuration.innerBoundary.value = innerBins[k];
          if (savedOuter.policy == BoundaryPolicy::SpecifiedIncomingEnergy ||
              savedOuter.policy == BoundaryPolicy::SpecifiedIncomingFlux)
            state->configuration.outerBoundary.value = outerBins[k];
          for (std::size_t i = 0; i < state->cells.size(); ++i) {
            state->cells[i].ePlusJ = state->cells[i].spectralEnergyJ[k];
            state->cells[i].eMinusJ = state->cells[i].spectralEnergyJ[bins + k];
          }
          AdvectBranch(state, true, substep, &result.ledger,
                       &result.diagnostics);
          AdvectBranch(state, false, substep, &result.ledger,
                       &result.diagnostics);
          result.diagnostics.spectralBinUpdates +=
              UINT64_C(2) * static_cast<std::uint64_t>(state->cells.size());
          for (std::size_t i = 0; i < state->cells.size(); ++i) {
            state->cells[i].spectralEnergyJ[k] = state->cells[i].ePlusJ;
            state->cells[i].spectralEnergyJ[bins + k] = state->cells[i].eMinusJ;
          }
        }
        state->configuration.innerBoundary = savedInner;
        state->configuration.outerBoundary = savedOuter;
        result.status = SynchronizeDerivedRepresentation(state);
        if (!result.status.ok()) return result;
      }
    }
    state->configuration.innerBoundary = advanceInner;
    state->configuration.outerBoundary = advanceOuter;
    result.diagnostics.advectionSubcycles = std::max<std::uint64_t>(1, cycles);
  }
  state->phase = OperatorPhase::AdvectionApplied;

  if (state->configuration.reflectionEnabled) {
    const double reflectionDtS = dtS / static_cast<double>(operatorSubsteps);
    for (std::uint64_t stage = 0; stage < operatorSubsteps; ++stage) {
      for (std::size_t i = 0; i < state->cells.size(); ++i) {
        ++result.diagnostics.reflectionCellUpdates;
        CellState& cell = state->cells[i];
        const double rate = 0.5 * state->configuration.reflectionCoefficient *
            std::fabs(cell.alfvenSpeedMPerS * cell.dLnAlfvenSpeeddsPerM);
        const std::size_t count = state->configuration.representation ==
            Representation::Spectral ? state->configuration.spectralBins : 1;
        for (std::size_t k = 0; k < count; ++k) {
          if (state->configuration.representation == Representation::Spectral)
            result.diagnostics.spectralBinUpdates += 2;
          double* plus = state->configuration.representation == Representation::Spectral
              ? &cell.spectralEnergyJ[k] : &cell.ePlusJ;
          double* minus = state->configuration.representation == Representation::Spectral
              ? &cell.spectralEnergyJ[state->configuration.spectralBins + k]
              : &cell.eMinusJ;
          const double sum = *plus + *minus;
          const double difference = (*plus - *minus) *
                                    std::exp(-2.0 * rate * reflectionDtS);
          const double oldPlus = *plus;
          *plus = 0.5 * (sum + difference);
          *minus = 0.5 * (sum - difference);
          result.ledger.reflectionTransferJ += *plus - oldPlus;
        }
      }
    }
    result.diagnostics.reflectionSubcycles = operatorSubsteps;
    result.status = SynchronizeDerivedRepresentation(state);
    if (!result.status.ok()) return result;
  }
  state->phase = OperatorPhase::ReflectionApplied;

  if (state->configuration.cascadeEnabled) {
    const double cascadeDtS = dtS / static_cast<double>(operatorSubsteps);
    for (std::uint64_t stage = 0; stage < operatorSubsteps; ++stage) {
      for (std::size_t i = 0; i < state->cells.size(); ++i) {
        ++result.diagnostics.cascadeCellUpdates;
        CellState& cell = state->cells[i];
        const double total = CellEnergy(cell);
        if (total == 0.0) {
          // Empty wave state is a valid physical zero.  Count it once per
          // selected stage so diagnostics distinguish it from a skipped cell.
          ++result.diagnostics.physicalZeroEvents;
          continue;
        }
      // The nonlinear turnover rate is delta-v/lambda_perp.  Using the cell's
      // mass density keeps the units explicit: sqrt((J/m^3)/(kg/m^3)) is m/s.
      if (!(cell.massDensityKgPerM3 > 0.0)) {
        result.status = Error("cascade requires positive cell mass density");
        return result;
      }
      const double rate = state->configuration.cascadeCoefficient *
          std::sqrt(total / (cell.volumeM3 * cell.massDensityKgPerM3)) /
          state->configuration.perpendicularCorrelationLengthM;
      if (state->configuration.representation == Representation::Integrated) {
        const double remaining = total / (1.0 + rate * cascadeDtS);
        const double factor = remaining / total;
        cell.ePlusJ *= factor;
        cell.eMinusJ *= factor;
        result.ledger.physicalDissipationJ += total - remaining;
      }
      else {
        // The backward bin sweep prevents a packet from traversing several
        // bins in one update. Transfers between represented bins cancel; only
        // the final-bin flux is classified as physical dissipation.
        const double fraction = rate * cascadeDtS /
                                (1.0 + rate * cascadeDtS);
        const std::size_t bins = state->configuration.spectralBins;
        for (int branch = 0; branch < 2; ++branch) {
          const std::size_t offset = branch == 0 ? 0 : bins;
          for (std::size_t reverse = bins; reverse > 0; --reverse) {
            ++result.diagnostics.spectralBinUpdates;
            const std::size_t k = reverse - 1;
            const double transfer = cell.spectralEnergyJ[offset + k] * fraction;
            cell.spectralEnergyJ[offset + k] -= transfer;
            if (k + 1 < bins)
              cell.spectralEnergyJ[offset + k + 1] += transfer;
            else
              result.ledger.physicalDissipationJ += transfer;
          }
        }
      }
      // Cascade transfer is internal to the represented spectrum. It is zero
      // in branch-integrated mode and cancels across bins in spectral mode.
        result.ledger.cascadeTransferJ += 0.0;
      }
    }
    result.diagnostics.cascadeSubcycles = operatorSubsteps;
    result.status = SynchronizeDerivedRepresentation(state);
    if (!result.status.ok()) return result;
  }
  state->phase = OperatorPhase::CascadeApplied;

  if (state->configuration.representation == Representation::Spectral) {
    result.status = SynchronizeDerivedRepresentation(state);
    if (!result.status.ok()) return result;
  }
  state->phase = OperatorPhase::Synchronized;
  result.ledger.finalEnergyJ = TotalEnergy(*state);
  result.ledger.closureResidualJ = result.ledger.finalEnergyJ -
      (result.ledger.initialEnergyJ + result.ledger.innerBoundaryJ +
       result.ledger.outerBoundaryJ + result.ledger.shockSourceJ +
       result.ledger.particleExchangeJ - result.ledger.physicalDissipationJ);
  AddLedger(&state->accumulatedLedger, result.ledger);
  state->epochS += dtS;
  ++state->completedSteps;
  state->phase = OperatorPhase::Ready;
  result.status = Transport::Status::Ok();
  return result;
}

Transport::Status RemapConservatively(const State& oldState,
                                      const std::vector<CellState>& newGeometry,
                                      State* remapped, EnergyLedger* ledger) {
  if (!remapped || !ledger || newGeometry.empty())
    return Error("remap requires output, ledger, and non-empty new geometry");
  *remapped = oldState;
  remapped->cells = newGeometry;
  const std::size_t bins = oldState.configuration.spectralBins;
  for (std::size_t j = 0; j < remapped->cells.size(); ++j) {
    remapped->cells[j].ePlusJ = remapped->cells[j].eMinusJ = 0.0;
    remapped->cells[j].pendingParticlePlusJ = 0.0;
    remapped->cells[j].pendingParticleMinusJ = 0.0;
    remapped->cells[j].pendingShockPlusJ = 0.0;
    remapped->cells[j].pendingShockMinusJ = 0.0;
    if (oldState.configuration.representation == Representation::Spectral)
      remapped->cells[j].spectralEnergyJ.assign(2 * bins, 0.0);
    else remapped->cells[j].spectralEnergyJ.clear();
  }
  std::vector<double> oldStart(oldState.cells.size() + 1, 0.0);
  std::vector<double> newStart(newGeometry.size() + 1, 0.0);
  for (std::size_t i = 0; i < oldState.cells.size(); ++i)
    oldStart[i + 1] = oldStart[i] + oldState.cells[i].lengthM;
  for (std::size_t j = 0; j < newGeometry.size(); ++j)
    newStart[j + 1] = newStart[j] + newGeometry[j].lengthM;
  const double domain = std::min(oldStart.back(), newStart.back());
  if (!(domain > 0.0)) return Error("remap domains must have positive length");
  for (std::size_t i = 0; i < oldState.cells.size(); ++i) {
    for (std::size_t j = 0; j < newGeometry.size(); ++j) {
      const double overlap = std::max(0.0, std::min(oldStart[i + 1], newStart[j + 1]) -
                                           std::max(oldStart[i], newStart[j]));
      if (overlap == 0.0) continue;
      const double fraction = overlap / oldState.cells[i].lengthM;
      remapped->cells[j].ePlusJ += oldState.cells[i].ePlusJ * fraction;
      remapped->cells[j].eMinusJ += oldState.cells[i].eMinusJ * fraction;
      if (oldState.configuration.representation == Representation::Spectral)
        for (std::size_t k = 0; k < 2 * bins; ++k)
          remapped->cells[j].spectralEnergyJ[k] +=
              oldState.cells[i].spectralEnergyJ[k] * fraction;
    }
  }
  ledger->remapCorrectionJ = TotalEnergy(*remapped) - TotalEnergy(oldState);
  ++remapped->fieldLineGeneration;
  return SynchronizeDerivedRepresentation(remapped);
}

Transport::Status SerializeCheckpoint(const State& state, std::string* text) {
  if (!text) return Error("checkpoint output pointer is null");
  std::ostringstream out;
  out << "SEP_TURBULENCE_CHECKPOINT 2\n" << std::setprecision(17)
      << static_cast<int>(state.configuration.source) << ' '
      << static_cast<int>(state.configuration.representation) << ' '
      << static_cast<int>(state.configuration.coupling) << '\n';
  WriteBoundary(out, state.configuration.innerBoundary);
  WriteBoundary(out, state.configuration.outerBoundary);
  out << state.configuration.advectionEnabled << ' '
      << state.configuration.reflectionEnabled << ' '
      << state.configuration.cascadeEnabled << ' '
      << state.configuration.shockInjectionEnabled << ' '
      << state.configuration.periodicBoundaries << '\n'
      << state.configuration.cflSafety << ' '
      << static_cast<int>(state.configuration.operatorSplitting) << ' '
      << state.configuration.operatorAccuracySafety << ' '
      << state.configuration.maximumSourceFraction << ' '
      << state.configuration.maximumCascadeFraction << ' '
      << state.configuration.minimumSubstepS << ' '
      << state.configuration.maximumSubsteps << ' '
      << state.configuration.reflectionCoefficient << ' '
      << state.configuration.cascadeCoefficient << ' '
      << state.configuration.perpendicularCorrelationLengthM << ' '
      << state.configuration.electronHeatingFraction << '\n'
      << state.configuration.spectralKMinPerM << ' '
      << state.configuration.spectralKMaxPerM << ' '
      << state.configuration.spectralBins << ' '
      << state.configuration.diagnosticCadence << ' '
      << state.configuration.conservationRelativeTolerance << '\n'
      << state.epochS << ' ' << state.fieldLineGeneration << ' '
      << state.campaignSeed << ' ' << state.completedSteps << ' '
      << static_cast<int>(state.phase) << ' ' << state.handoffCompleted << ' '
      << state.handoffEpochS << '\n'
      << HexEncode(state.sourceChecksum) << '\n'
      << HexEncode(state.provenance) << '\n';
  const std::vector<double> ledger = LedgerValues(state.accumulatedLedger);
  for (std::size_t i = 0; i < ledger.size(); ++i)
    out << ledger[i] << (i + 1 == ledger.size() ? '\n' : ' ');
  out << state.cells.size() << '\n';
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const CellState& c = state.cells[i];
    out << c.lengthM << ' ' << c.volumeM3 << ' ' << c.plasmaSpeedMPerS << ' '
        << c.alfvenSpeedMPerS << ' ' << c.dLnAlfvenSpeeddsPerM << ' '
        << c.magneticFieldT << ' ' << c.massDensityKgPerM3 << ' '
        << c.ePlusJ << ' ' << c.eMinusJ << ' '
        << c.pendingParticlePlusJ << ' ' << c.pendingParticleMinusJ << ' '
        << c.pendingShockPlusJ << ' ' << c.pendingShockMinusJ << ' '
        << c.spectralEnergyJ.size();
    for (std::size_t k = 0; k < c.spectralEnergyJ.size(); ++k)
      out << ' ' << c.spectralEnergyJ[k];
    out << '\n';
  }
  *text = out.str();
  return Transport::Status::Ok();
}

Transport::Status DeserializeCheckpoint(const std::string& text, State* state) {
  if (!state) return Error("checkpoint state pointer is null");
  std::istringstream in(text);
  std::string magic;
  int version = 0, source = 0, representation = 0, coupling = 0, phase = 0;
  int splitting = 0;
  State parsed;
  if (!(in >> magic >> version) || magic != "SEP_TURBULENCE_CHECKPOINT" || version != 2)
    return Error("unsupported turbulence checkpoint schema");
  if (!(in >> source >> representation >> coupling) || source < 0 || source > 4 ||
      representation < 0 || representation > 1 || coupling < 0 || coupling > 1)
    return Error("checkpoint contains invalid source/representation metadata");
  parsed.configuration.source = static_cast<Source>(source);
  parsed.configuration.representation = static_cast<Representation>(representation);
  parsed.configuration.coupling = static_cast<CouplingPolicy>(coupling);
  if (!ReadBoundary(in, &parsed.configuration.innerBoundary) ||
      !ReadBoundary(in, &parsed.configuration.outerBoundary))
    return Error("checkpoint boundary metadata is malformed");
  if (!(in >> parsed.configuration.advectionEnabled >>
        parsed.configuration.reflectionEnabled >> parsed.configuration.cascadeEnabled >>
        parsed.configuration.shockInjectionEnabled >> parsed.configuration.periodicBoundaries >>
        parsed.configuration.cflSafety >> splitting >>
        parsed.configuration.operatorAccuracySafety >>
        parsed.configuration.maximumSourceFraction >>
        parsed.configuration.maximumCascadeFraction >>
        parsed.configuration.minimumSubstepS >>
        parsed.configuration.maximumSubsteps >>
        parsed.configuration.reflectionCoefficient >>
        parsed.configuration.cascadeCoefficient >>
        parsed.configuration.perpendicularCorrelationLengthM >>
        parsed.configuration.electronHeatingFraction >> parsed.configuration.spectralKMinPerM >>
        parsed.configuration.spectralKMaxPerM >> parsed.configuration.spectralBins >>
        parsed.configuration.diagnosticCadence >>
        parsed.configuration.conservationRelativeTolerance >> parsed.epochS >>
        parsed.fieldLineGeneration >> parsed.campaignSeed >> parsed.completedSteps >> phase >>
        parsed.handoffCompleted >> parsed.handoffEpochS))
    return Error("checkpoint header is malformed");
  if (splitting != static_cast<int>(OperatorSplitting::LieFirstOrder))
    return Error("checkpoint operator-splitting scheme is unsupported");
  parsed.configuration.operatorSplitting =
      static_cast<OperatorSplitting>(splitting);
  std::string checksumHex, provenanceHex;
  if (!(in >> checksumHex >> provenanceHex) ||
      !HexDecode(checksumHex, &parsed.sourceChecksum) ||
      !HexDecode(provenanceHex, &parsed.provenance))
    return Error("checkpoint provenance fields are malformed");
  if (phase < 0 || phase > 6) return Error("checkpoint operator phase is invalid");
  parsed.phase = static_cast<OperatorPhase>(phase);
  if (!ReadLedger(in, &parsed.accumulatedLedger))
    return Error("checkpoint ledger is malformed");
  std::size_t count = 0;
  if (!(in >> count) || count == 0) return Error("checkpoint cell count is invalid");
  parsed.cells.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    CellState& c = parsed.cells[i];
    std::size_t spectralSize = 0;
    if (!(in >> c.lengthM >> c.volumeM3 >> c.plasmaSpeedMPerS >>
          c.alfvenSpeedMPerS >> c.dLnAlfvenSpeeddsPerM >> c.magneticFieldT >>
          c.massDensityKgPerM3 >> c.ePlusJ >> c.eMinusJ >>
          c.pendingParticlePlusJ >> c.pendingParticleMinusJ >>
          c.pendingShockPlusJ >> c.pendingShockMinusJ >> spectralSize))
      return Error("checkpoint cell record is malformed");
    c.spectralEnergyJ.resize(spectralSize);
    for (std::size_t k = 0; k < spectralSize; ++k)
      if (!(in >> c.spectralEnergyJ[k])) return Error("checkpoint spectrum is malformed");
  }
  Transport::Status status = ValidateConfiguration(parsed.configuration);
  if (!status.ok()) return status;
  // Validate cell ownership and dimensions on a temporary copy. The published
  // object preserves the exact saved phase; InitializeState resets that phase
  // and is therefore never applied directly to the restored object.
  State validationCopy = parsed;
  status = InitializeState(&validationCopy);
  if (!status.ok())
    return Error("checkpoint state failed validation: " + status.message);
  std::string trailing;
  if (in >> trailing) return Error("checkpoint contains trailing records");
  *state = parsed;
  return Transport::Status::Ok();
}

std::uint64_t StateHash(const State& state) {
  std::string checkpoint;
  if (!SerializeCheckpoint(state, &checkpoint).ok()) return 0;
  std::uint64_t hash = 1469598103934665603ULL;
  for (std::size_t i = 0; i < checkpoint.size(); ++i) {
    hash ^= static_cast<unsigned char>(checkpoint[i]);
    hash *= 1099511628211ULL;
  }
  return hash;
}

}  // namespace Turbulence
}  // namespace SEP
