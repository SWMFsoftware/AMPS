#include "sep_coronal_cme/mhd_jump_solver.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

double SafeScale(double a, double b = 0.0) {
  return std::max({std::abs(a), std::abs(b), 1.0e-300});
}

Vec3 Tangential(Vec3 value, Vec3 normal) {
  return value - Dot(value, normal) * normal;
}

double EnergyFlux(const MhdPrimitiveState& state, Vec3 normal,
                  double shockSpeedMPerS, double gammaAdiabatic) {
  const Vec3 shockVelocity = state.velocityMPerS -
      shockSpeedMPerS * normal;
  const double normalVelocity = Dot(shockVelocity, normal);
  const double normalField = Dot(state.magneticFieldT, normal);
  const double magneticSquared = Dot(state.magneticFieldT,
                                     state.magneticFieldT);
  const double totalSpecificFlux =
      0.5 * state.massDensityKgM3 * Dot(shockVelocity, shockVelocity) +
      gammaAdiabatic / (gammaAdiabatic - 1.0) * state.pressurePa +
      magneticSquared / Constants::kVacuumPermeabilityHPerM;
  return normalVelocity * totalSpecificFlux - normalField *
      Dot(shockVelocity, state.magneticFieldT) /
      Constants::kVacuumPermeabilityHPerM;
}

struct Candidate {
  bool valid = false;
  MhdPrimitiveState downstream;
  double energyResidual = 0.0;
};

Candidate DownstreamForCompression(
    const MhdPrimitiveState& upstream, Vec3 normal, double shockSpeedMPerS,
    double gammaAdiabatic, double compression) {
  Candidate result;
  if (!(compression > 1.0)) return result;
  const double mu0 = Constants::kVacuumPermeabilityHPerM;
  const double labNormal1 = Dot(upstream.velocityMPerS, normal);
  const double normalVelocity1 = labNormal1 - shockSpeedMPerS;
  const Vec3 tangentialVelocity1 = Tangential(upstream.velocityMPerS, normal);
  const double normalField = Dot(upstream.magneticFieldT, normal);
  const Vec3 tangentialField1 = Tangential(upstream.magneticFieldT, normal);
  const double massFlux = upstream.massDensityKgM3 * normalVelocity1;
  if (massFlux == 0.0) return result;
  const double normalVelocity2 = normalVelocity1 / compression;
  const double magneticTerm = normalField * normalField / (mu0 * massFlux);
  const double denominator = normalVelocity2 - magneticTerm;
  const double numerator = normalVelocity1 - magneticTerm;
  if (std::abs(denominator) <= 1.0e-14 *
      std::max(std::abs(normalVelocity1), 1.0)) return result;
  const Vec3 tangentialField2 =
      (numerator / denominator) * tangentialField1;
  const Vec3 tangentialVelocity2 = tangentialVelocity1 +
      normalField / (mu0 * massFlux) *
      (tangentialField2 - tangentialField1);
  const double pressure2 = upstream.pressurePa +
      upstream.massDensityKgM3 * normalVelocity1 * normalVelocity1 *
          (1.0 - 1.0 / compression) +
      (Dot(tangentialField1, tangentialField1) -
       Dot(tangentialField2, tangentialField2)) / (2.0 * mu0);
  if (!(Finite(pressure2) && pressure2 > 0.0)) return result;
  result.downstream.massDensityKgM3 =
      compression * upstream.massDensityKgM3;
  result.downstream.pressurePa = pressure2;
  result.downstream.velocityMPerS = tangentialVelocity2 +
      (normalVelocity2 + shockSpeedMPerS) * normal;
  result.downstream.magneticFieldT = tangentialField2 +
      normalField * normal;
  result.energyResidual = EnergyFlux(result.downstream, normal,
      shockSpeedMPerS, gammaAdiabatic) -
      EnergyFlux(upstream, normal, shockSpeedMPerS, gammaAdiabatic);
  result.valid = Finite(result.energyResidual);
  return result;
}

RankineHugoniotResiduals ComputeResiduals(
    const MhdPrimitiveState& upstream,
    const MhdPrimitiveState& downstream, Vec3 normal,
    double shockSpeedMPerS, double gammaAdiabatic) {
  const double mu0 = Constants::kVacuumPermeabilityHPerM;
  const Vec3 velocity1 = upstream.velocityMPerS - shockSpeedMPerS * normal;
  const Vec3 velocity2 = downstream.velocityMPerS - shockSpeedMPerS * normal;
  const double vn1 = Dot(velocity1, normal), vn2 = Dot(velocity2, normal);
  const double bn1 = Dot(upstream.magneticFieldT, normal);
  const double bn2 = Dot(downstream.magneticFieldT, normal);
  const Vec3 vt1 = Tangential(velocity1, normal);
  const Vec3 vt2 = Tangential(velocity2, normal);
  const Vec3 bt1 = Tangential(upstream.magneticFieldT, normal);
  const Vec3 bt2 = Tangential(downstream.magneticFieldT, normal);
  RankineHugoniotResiduals result;
  const double mass1 = upstream.massDensityKgM3 * vn1;
  const double mass2 = downstream.massDensityKgM3 * vn2;
  result.mass = std::abs(mass2 - mass1) / SafeScale(mass1, mass2);
  result.normalMagnetic = std::abs(bn2 - bn1) /
      SafeScale(Norm(upstream.magneticFieldT), Norm(downstream.magneticFieldT));
  const Vec3 electric1 = vn1 * bt1 - bn1 * vt1;
  const Vec3 electric2 = vn2 * bt2 - bn2 * vt2;
  result.tangentialElectric = Norm(electric2 - electric1) /
      SafeScale(Norm(electric1), Norm(electric2));
  const double normalMomentum1 = upstream.massDensityKgM3 * vn1 * vn1 +
      upstream.pressurePa + Dot(bt1, bt1) / (2.0 * mu0);
  const double normalMomentum2 = downstream.massDensityKgM3 * vn2 * vn2 +
      downstream.pressurePa + Dot(bt2, bt2) / (2.0 * mu0);
  result.normalMomentum = std::abs(normalMomentum2 - normalMomentum1) /
      SafeScale(normalMomentum1, normalMomentum2);
  const Vec3 tangentialMomentum1 = mass1 * vt1 -
      bn1 / mu0 * bt1;
  const Vec3 tangentialMomentum2 = mass2 * vt2 -
      bn2 / mu0 * bt2;
  result.tangentialMomentum = Norm(tangentialMomentum2 -
      tangentialMomentum1) /
      SafeScale(Norm(tangentialMomentum1), Norm(tangentialMomentum2));
  const double energy1 = EnergyFlux(upstream, normal,
                                    shockSpeedMPerS, gammaAdiabatic);
  const double energy2 = EnergyFlux(downstream, normal,
                                    shockSpeedMPerS, gammaAdiabatic);
  result.totalEnergy = std::abs(energy2 - energy1) /
      SafeScale(energy1, energy2);
  result.maximum = std::max({result.mass, result.normalMagnetic,
      result.tangentialElectric, result.normalMomentum,
      result.tangentialMomentum, result.totalEnergy});
  return result;
}

bool StrictlyIncreasing(const std::vector<double>& values) {
  for (std::size_t index = 1; index < values.size(); ++index) {
    if (!(Finite(values[index]) && values[index] > values[index - 1]))
      return false;
  }
  return !values.empty() && Finite(values.front());
}

std::size_t LowerInterval(const std::vector<double>& grid, double value) {
  const auto upper = std::upper_bound(grid.begin(), grid.end(), value);
  if (upper == grid.begin()) return 0;
  if (upper == grid.end()) return grid.size() - 2;
  return static_cast<std::size_t>(upper - grid.begin() - 1);
}

}  // namespace

Core::Result<MhdCharacteristicSpeeds> EvaluateMhdCharacteristics(
    const MhdPrimitiveState& state, Vec3 normal, double gammaAdiabatic) {
  normal = Unit(normal);
  if (!(state.massDensityKgM3 > 0.0 && state.pressurePa > 0.0 &&
        gammaAdiabatic > 1.0 && Norm(normal) > 0.0 &&
        Finite(state.velocityMPerS.x) && Finite(state.velocityMPerS.y) &&
        Finite(state.velocityMPerS.z) && Finite(state.magneticFieldT.x) &&
        Finite(state.magneticFieldT.y) && Finite(state.magneticFieldT.z))) {
    return Core::Result<MhdCharacteristicSpeeds>::Failure(
        Core::StatusCode::InvalidState,
        "MHD characteristics require positive rho/p/gamma and finite vectors");
  }
  const double fieldSquared = Dot(state.magneticFieldT, state.magneticFieldT);
  const double normalField = Dot(state.magneticFieldT, normal);
  const double soundSquared = gammaAdiabatic * state.pressurePa /
      state.massDensityKgM3;
  const double alfvenSquared = fieldSquared /
      (Constants::kVacuumPermeabilityHPerM * state.massDensityKgM3);
  const double normalAlfvenSquared = normalField * normalField /
      (Constants::kVacuumPermeabilityHPerM * state.massDensityKgM3);
  const double sum = soundSquared + alfvenSquared;
  const double discriminant = std::max(0.0,
      sum * sum - 4.0 * soundSquared * normalAlfvenSquared);
  MhdCharacteristicSpeeds result;
  result.soundSpeedMPerS = std::sqrt(soundSquared);
  result.alfvenSpeedMPerS = std::sqrt(alfvenSquared);
  result.normalAlfvenSpeedMPerS = std::sqrt(normalAlfvenSquared);
  result.fastSpeedMPerS = std::sqrt(0.5 * (sum + std::sqrt(discriminant)));
  result.slowSpeedMPerS = std::sqrt(std::max(0.0,
      0.5 * (sum - std::sqrt(discriminant))));
  if (fieldSquared > 0.0) {
    result.obliquityRad = std::acos(std::max(0.0, std::min(1.0,
        std::abs(normalField) / std::sqrt(fieldSquared))));
    result.plasmaBeta = 2.0 * Constants::kVacuumPermeabilityHPerM *
        state.pressurePa / fieldSquared;
  } else {
    result.obliquityRad = 0.0;
    result.plasmaBeta = std::numeric_limits<double>::infinity();
  }
  return Core::Result<MhdCharacteristicSpeeds>::Success(result);
}

Core::Result<MhdShockSolution> SolveObliqueFastShock(
    const MhdPrimitiveState& upstream, Vec3 outwardNormal,
    double shockNormalSpeedMPerS, double gammaAdiabatic,
    double residualTolerance) {
  const Vec3 normal = Unit(outwardNormal);
  if (!(Norm(normal) > 0.0 && Finite(shockNormalSpeedMPerS) &&
        gammaAdiabatic > 1.0 && residualTolerance > 0.0)) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "shock solve requires a normal, gamma>1, and positive tolerance");
  }
  auto upstreamCharacteristics = EvaluateMhdCharacteristics(
      upstream, normal, gammaAdiabatic);
  if (!upstreamCharacteristics.ok()) {
    return Core::Result<MhdShockSolution>::Failure(
        upstreamCharacteristics.status.code,
        upstreamCharacteristics.status.message);
  }
  const double inflow = shockNormalSpeedMPerS -
      Dot(upstream.velocityMPerS, normal);
  if (!(inflow > 0.0)) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::InvalidState,
        "candidate front is an expansion/outflow in the shock frame");
  }
  const double fastMach = inflow /
      upstreamCharacteristics.value.fastSpeedMPerS;
  if (!(fastMach > 1.0)) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::OutOfDomain,
        "candidate front is sub-fast and has no fast-shock downstream state");
  }
  const double maximumCompression = (gammaAdiabatic + 1.0) /
      (gammaAdiabatic - 1.0);
  const double lowerCompression = 1.0 + 1.0e-6;
  const double upperCompression = maximumCompression * (1.0 - 1.0e-10);
  const double energyScale = SafeScale(EnergyFlux(
      upstream, normal, shockNormalSpeedMPerS, gammaAdiabatic));

  bool bracketed = false;
  double left = 0.0, right = 0.0, fLeft = 0.0, fRight = 0.0;
  Candidate previous;
  double previousCompression = 0.0;
  constexpr int kScanIntervals = 4096;
  for (int index = 0; index <= kScanIntervals; ++index) {
    const double compression = lowerCompression +
        (upperCompression - lowerCompression) * index / kScanIntervals;
    const Candidate candidate = DownstreamForCompression(
        upstream, normal, shockNormalSpeedMPerS, gammaAdiabatic, compression);
    if (!candidate.valid) continue;
    const double normalized = candidate.energyResidual / energyScale;
    if (previous.valid && normalized * fLeft <= 0.0) {
      left = previousCompression;
      right = compression;
      fRight = normalized;
      bracketed = true;
      break;
    }
    previous = candidate;
    previousCompression = compression;
    fLeft = normalized;
  }
  if (!bracketed) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::NumericalFailure,
        "compressive fast-shock energy root is not bracketed");
  }
  Candidate root;
  for (int iteration = 0; iteration < 160; ++iteration) {
    const double middle = 0.5 * (left + right);
    root = DownstreamForCompression(upstream, normal, shockNormalSpeedMPerS,
                                    gammaAdiabatic, middle);
    if (!root.valid) {
      return Core::Result<MhdShockSolution>::Failure(
          Core::StatusCode::NumericalFailure,
          "fast-shock bracket crossed a singular/nonphysical branch");
    }
    const double fMiddle = root.energyResidual / energyScale;
    if (std::abs(fMiddle) <= residualTolerance * 0.1 ||
        right - left <= 1.0e-12 * middle) {
      left = right = middle;
      break;
    }
    if (fLeft * fMiddle <= 0.0) {
      right = middle;
      fRight = fMiddle;
    } else {
      left = middle;
      fLeft = fMiddle;
    }
  }
  (void)fRight;
  const double compression = 0.5 * (left + right);
  root = DownstreamForCompression(upstream, normal, shockNormalSpeedMPerS,
                                  gammaAdiabatic, compression);
  if (!root.valid) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::NumericalFailure,
        "bracketed shock root did not produce a downstream state");
  }
  auto downstreamCharacteristics = EvaluateMhdCharacteristics(
      root.downstream, normal, gammaAdiabatic);
  if (!downstreamCharacteristics.ok()) {
    return Core::Result<MhdShockSolution>::Failure(
        downstreamCharacteristics.status.code,
        downstreamCharacteristics.status.message);
  }
  const double downstreamInflow = shockNormalSpeedMPerS -
      Dot(root.downstream.velocityMPerS, normal);
  const double entropy = std::log(root.downstream.pressurePa /
      upstream.pressurePa) - gammaAdiabatic * std::log(compression);
  if (!(compression > 1.0 && compression < maximumCompression &&
        entropy > 0.0 && downstreamInflow > 0.0 &&
        downstreamInflow < downstreamCharacteristics.value.fastSpeedMPerS)) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::InvalidState,
        "root fails compression, entropy, or fast-characteristic ordering");
  }
  const RankineHugoniotResiduals residuals = ComputeResiduals(
      upstream, root.downstream, normal, shockNormalSpeedMPerS,
      gammaAdiabatic);
  if (!(Finite(residuals.maximum) && residuals.maximum <= residualTolerance)) {
    return Core::Result<MhdShockSolution>::Failure(
        Core::StatusCode::NumericalFailure,
        "accepted branch exceeds normalized Rankine-Hugoniot tolerance");
  }
  MhdShockSolution result;
  result.downstream = root.downstream;
  result.upstreamCharacteristics = upstreamCharacteristics.value;
  result.downstreamCharacteristics = downstreamCharacteristics.value;
  result.compressionRatio = compression;
  result.upstreamFastMach = fastMach;
  result.upstreamAlfvenMach = upstreamCharacteristics.value.alfvenSpeedMPerS > 0.0
      ? inflow / upstreamCharacteristics.value.alfvenSpeedMPerS
      : std::numeric_limits<double>::infinity();
  result.entropyLogIncrement = entropy;
  result.residuals = residuals;
  result.branch = Norm(Tangential(upstream.magneticFieldT, normal)) < 1.0e-14 *
      std::max(Norm(upstream.magneticFieldT), 1.0e-30)
      ? "regular-parallel-fast" : "regular-oblique-fast";
  return Core::Result<MhdShockSolution>::Success(std::move(result));
}

Core::Result<CriticalMachQuery> QueryCriticalMach(
    const CriticalMachTable& table, double beta, double obliquityRad,
    double gammaAdiabatic, MachConvention requestedConvention,
    bool exactNormalFieldIsZero) {
  if (table.version.empty() || table.checksum.empty() ||
      table.betaGrid.size() < 2 || table.obliquityGridRad.size() < 2 ||
      !StrictlyIncreasing(table.betaGrid) ||
      !StrictlyIncreasing(table.obliquityGridRad) ||
      table.values.size() != table.betaGrid.size() *
          table.obliquityGridRad.size() ||
      std::abs(table.obliquityGridRad.front()) > 1.0e-14 ||
      std::abs(table.obliquityGridRad.back() - 0.5 * Constants::kPi) >
          1.0e-14) {
    return Core::Result<CriticalMachQuery>::Failure(
        Core::StatusCode::DataIntegrityFailure,
        "critical-Mach table is malformed or lacks [0,pi/2] coverage");
  }
  for (double value : table.values) {
    if (!(Finite(value) && value > 0.0)) {
      return Core::Result<CriticalMachQuery>::Failure(
          Core::StatusCode::DataIntegrityFailure,
          "critical-Mach table contains a nonpositive/nonfinite value");
    }
  }
  if (std::abs(table.gammaAdiabatic - gammaAdiabatic) > 1.0e-12 ||
      table.convention != requestedConvention) {
    return Core::Result<CriticalMachQuery>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "critical-Mach EOS or Mach convention does not match the run");
  }
  if (!(Finite(beta) && beta >= 0.0 && Finite(obliquityRad) &&
        obliquityRad >= 0.0 && obliquityRad <= 0.5 * Constants::kPi)) {
    return Core::Result<CriticalMachQuery>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "critical-Mach query has invalid beta or obliquity");
  }
  if (requestedConvention == MachConvention::NormalAlfven &&
      exactNormalFieldIsZero) {
    return Core::Result<CriticalMachQuery>::Success(
        {CriticalQueryValidity::NormalAlfvenInapplicable, 0.0,
         "exact B_n=0 has no finite normal-Alfven Mach query"});
  }
  if (beta < table.betaGrid.front() || beta > table.betaGrid.back()) {
    return Core::Result<CriticalMachQuery>::Success(
        {CriticalQueryValidity::BetaOutsideCoverage, 0.0,
         "upstream beta is outside the versioned table"});
  }
  const std::size_t i = LowerInterval(table.betaGrid, beta);
  const std::size_t j = LowerInterval(table.obliquityGridRad, obliquityRad);
  const double betaWeight = (beta - table.betaGrid[i]) /
      (table.betaGrid[i + 1] - table.betaGrid[i]);
  const double thetaWeight =
      (obliquityRad - table.obliquityGridRad[j]) /
      (table.obliquityGridRad[j + 1] - table.obliquityGridRad[j]);
  const std::size_t columns = table.obliquityGridRad.size();
  const auto at = [&](std::size_t row, std::size_t column) {
    return table.values[row * columns + column];
  };
  const double lower = (1.0 - thetaWeight) * at(i, j) +
      thetaWeight * at(i, j + 1);
  const double upper = (1.0 - thetaWeight) * at(i + 1, j) +
      thetaWeight * at(i + 1, j + 1);
  return Core::Result<CriticalMachQuery>::Success(
      {CriticalQueryValidity::Valid,
       (1.0 - betaWeight) * lower + betaWeight * upper, ""});
}

Core::Result<CriticalityDecision> ApplyCriticalityPolicy(
    double matchingMachNumber, const CriticalMachQuery& query,
    bool sourceRequiresSupercritical, CriticalityMissPolicy policy) {
  if (!(Finite(matchingMachNumber) && matchingMachNumber > 0.0)) {
    return Core::Result<CriticalityDecision>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "criticality classification requires a positive matching Mach number");
  }
  CriticalityDecision result;
  if (query.validity == CriticalQueryValidity::Valid) {
    result.criticalityAvailable = true;
    result.supercritical = matchingMachNumber > query.value;
    result.sourceEligible = !sourceRequiresSupercritical || result.supercritical;
    return Core::Result<CriticalityDecision>::Success(result);
  }
  if (!sourceRequiresSupercritical &&
      policy == CriticalityMissPolicy::DiagnosticOnlyFast) {
    result.sourceEligible = true;
    return Core::Result<CriticalityDecision>::Success(result);
  }
  if (policy == CriticalityMissPolicy::ExcludeSourceBudgeted) {
    result.excludedWithoutRenormalization = true;
    return Core::Result<CriticalityDecision>::Success(result);
  }
  return Core::Result<CriticalityDecision>::Failure(
      Core::StatusCode::OutOfDomain,
      "critical-Mach coverage miss is fatal under the selected policy");
}

Core::Result<CriticalCoverageLedger> PreflightCriticalCoverage(
    const CriticalMachTable& table,
    const std::vector<CriticalCoverageSample>& samples,
    double gammaAdiabatic, MachConvention convention,
    double maximumAreaFraction, double maximumNumberFraction,
    double maximumEnergyFraction) {
  if (!(maximumAreaFraction >= 0.0 && maximumAreaFraction <= 1.0 &&
        maximumNumberFraction >= 0.0 && maximumNumberFraction <= 1.0 &&
        maximumEnergyFraction >= 0.0 && maximumEnergyFraction <= 1.0)) {
    return Core::Result<CriticalCoverageLedger>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "critical-coverage budgets must lie in [0,1]");
  }
  CriticalCoverageLedger result;
  for (const auto& sample : samples) {
    if (!(sample.areaM2 >= 0.0 && sample.incidentNumberRatePerS >= 0.0 &&
          sample.incidentKineticEnergyRateW >= 0.0)) {
      return Core::Result<CriticalCoverageLedger>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "critical-coverage physical measures must be nonnegative");
    }
    result.candidateAreaM2 += sample.areaM2;
    result.candidateNumberRatePerS += sample.incidentNumberRatePerS;
    result.candidateEnergyRateW += sample.incidentKineticEnergyRateW;
    auto query = QueryCriticalMach(table, sample.beta, sample.obliquityRad,
        gammaAdiabatic, convention, sample.exactNormalFieldIsZero);
    if (!query.ok()) {
      return Core::Result<CriticalCoverageLedger>::Failure(
          query.status.code, query.status.message);
    }
    if (query.value.validity != CriticalQueryValidity::Valid) {
      result.unavailableAreaM2 += sample.areaM2;
      result.unavailableNumberRatePerS += sample.incidentNumberRatePerS;
      result.unavailableEnergyRateW += sample.incidentKineticEnergyRateW;
    }
  }
  result.noCandidateSupport = result.candidateAreaM2 == 0.0 &&
      result.candidateNumberRatePerS == 0.0 &&
      result.candidateEnergyRateW == 0.0;
  result.areaFraction = result.candidateAreaM2 > 0.0
      ? result.unavailableAreaM2 / result.candidateAreaM2 : 0.0;
  result.numberFraction = result.candidateNumberRatePerS > 0.0
      ? result.unavailableNumberRatePerS / result.candidateNumberRatePerS : 0.0;
  result.energyFraction = result.candidateEnergyRateW > 0.0
      ? result.unavailableEnergyRateW / result.candidateEnergyRateW : 0.0;
  if (result.areaFraction > maximumAreaFraction ||
      result.numberFraction > maximumNumberFraction ||
      result.energyFraction > maximumEnergyFraction) {
    return Core::Result<CriticalCoverageLedger>::Failure(
        Core::StatusCode::OutOfDomain,
        "complete-history critical-Mach exclusion budget was exceeded");
  }
  return Core::Result<CriticalCoverageLedger>::Success(result);
}

} }  // namespace SEP::CoronalCME
