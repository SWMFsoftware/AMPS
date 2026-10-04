#include "sep_coronal_cme/shock_provider.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

bool QualifiedProductionInterface(const PatchInterfaceEvidence& evidence,
                                  std::uint64_t generation) {
  if (evidence.backgroundGeneration != generation ||
      !evidence.provenanceComplete || !evidence.topologyPassed ||
      !evidence.kinematicPassed || !evidence.balancePassed ||
      evidence.interfaceIdentity.empty()) return false;
  return evidence.policy == InterfacePolicy::BoundedApproximation ||
      evidence.policy == InterfacePolicy::StationaryTangentialDiscontinuity;
}

Core::Result<BudgetRatio> Ratio(double excluded, double candidate,
                                double maximum) {
  return EvaluateBudgetRatio(excluded, candidate, maximum);
}

}  // namespace

Core::Result<std::vector<ShockPatchInput>> SplitShockPatch(
    const ShockPatchInput& parent, const std::vector<double>& areaFractions,
    const std::vector<PatchInterfaceEvidence>& childEvidence) {
  if (areaFractions.size() < 2 ||
      areaFractions.size() != childEvidence.size() ||
      !(parent.areaM2 > 0.0) || !Finite(parent.areaM2)) {
    return Core::Result<std::vector<ShockPatchInput>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "patch split requires positive parent area and matching child arrays");
  }
  double sum = 0.0;
  for (double fraction : areaFractions) {
    if (!(fraction > 0.0) || !Finite(fraction)) {
      return Core::Result<std::vector<ShockPatchInput>>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "patch split fractions must be finite and positive");
    }
    sum += fraction;
  }
  if (std::abs(sum - 1.0) > 1.0e-12) {
    return Core::Result<std::vector<ShockPatchInput>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "patch split fractions do not close to one");
  }

  std::vector<ShockPatchInput> children;
  children.reserve(areaFractions.size());
  for (std::size_t index = 0; index < areaFractions.size(); ++index) {
    ShockPatchInput child = parent;
    // IDs are a deterministic lineage encoding and do not depend on rank,
    // vector order outside this parent, or mesh ownership.
    child.parentId = parent.stableId;
    child.stableId = parent.stableId * 1024U + index + 1U;
    child.areaM2 *= areaFractions[index];
    child.incidentNumberRatePerS *= areaFractions[index];
    child.incidentKineticEnergyRateW *= areaFractions[index];
    child.interface = childEvidence[index];
    children.push_back(std::move(child));
  }
  return Core::Result<std::vector<ShockPatchInput>>::Success(
      std::move(children));
}

Core::Result<double> LocateFirstFastCrossing(
    double beginTimeS, double endTimeS,
    const std::function<double(double)>& fastMach, double timeToleranceS,
    int maximumIterations) {
  if (!(Finite(beginTimeS) && Finite(endTimeS) && endTimeS > beginTimeS &&
        timeToleranceS > 0.0 && Finite(timeToleranceS) &&
        maximumIterations > 0 && fastMach)) {
    return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "fast crossing requires a finite ordered bracket and tolerance");
  }
  double left = beginTimeS, right = endTimeS;
  double fLeft = fastMach(left) - 1.0;
  double fRight = fastMach(right) - 1.0;
  if (!(Finite(fLeft) && Finite(fRight)) || fLeft > 0.0 || fRight < 0.0) {
    return Core::Result<double>::Failure(
        Core::StatusCode::OutOfDomain,
        "fast crossing bracket must start sub-fast and end fast");
  }
  if (fLeft == 0.0) return Core::Result<double>::Success(left);
  for (int iteration = 0; iteration < maximumIterations; ++iteration) {
    const double middle = 0.5 * (left + right);
    const double fMiddle = fastMach(middle) - 1.0;
    if (!Finite(fMiddle)) {
      return Core::Result<double>::Failure(
          Core::StatusCode::NumericalFailure,
          "fast Mach function became nonfinite inside event bracket");
    }
    if (fMiddle >= 0.0) right = middle;
    else left = middle;
    if (right - left <= timeToleranceS)
      return Core::Result<double>::Success(0.5 * (left + right));
  }
  return Core::Result<double>::Failure(
      Core::StatusCode::NumericalFailure,
      "fast crossing did not converge within the iteration limit");
}

Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>
TransactionalShockProvider::Prepare(
    double timeS, std::uint64_t backgroundGeneration,
    const std::vector<ShockPatchInput>& inputs,
    const ShockPreparationOptions& options) {
  if (!Finite(timeS) || backgroundGeneration == 0 || inputs.empty() ||
      !Finite(options.gammaAdiabatic) || options.gammaAdiabatic <= 1.0 ||
      options.minimumFastAreaFraction < 0.0 ||
      options.minimumFastAreaFraction > 1.0) {
    return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "shock preparation requires time, background generation, and patches");
  }

  ShockSurfaceSnapshot candidate;
  candidate.generation = CurrentGeneration() + 1;
  candidate.backgroundGeneration = backgroundGeneration;
  candidate.timeS = timeS;
  candidate.patches.reserve(inputs.size());

  for (const ShockPatchInput& input : inputs) {
    if (input.stableId == 0 || !(input.areaM2 > 0.0) ||
        !Finite(input.areaM2) || input.incidentNumberRatePerS < 0.0 ||
        input.incidentKineticEnergyRateW < 0.0 ||
        !Finite(input.incidentNumberRatePerS) ||
        !Finite(input.incidentKineticEnergyRateW)) {
      return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
          Core::StatusCode::InvalidState,
          "shock patch has invalid identity, area, or physical source measure");
    }
    if (options.productionIntent &&
        !QualifiedProductionInterface(input.interface,
                                      backgroundGeneration)) {
      return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
          Core::StatusCode::DataIntegrityFailure,
          "production patch lacks qualified one-sided interface evidence");
    }

    auto characteristics = EvaluateMhdCharacteristics(
        input.upstream, input.outwardNormal, options.gammaAdiabatic);
    if (!characteristics.ok()) {
      return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
          characteristics.status.code, characteristics.status.message);
    }
    const Vec3 normal = Unit(input.outwardNormal);
    const double inflow = input.shockNormalSpeedMPerS -
        Dot(input.upstream.velocityMPerS, normal);
    const double mach = inflow / characteristics.value.fastSpeedMPerS;
    if (!Finite(mach)) {
      return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
          Core::StatusCode::NumericalFailure,
          "patch fast Mach number is nonfinite");
    }

    ShockPatchSnapshot patch;
    patch.stableId = input.stableId;
    patch.parentId = input.parentId;
    patch.areaM2 = input.areaM2;
    patch.fastMach = mach;
    patch.fast = mach > 1.0;
    patch.supercritical = patch.fast &&
        (!options.requireSupercritical ||
         mach >= options.criticalFastMach);
    patch.sourceTerminated = input.sourceTerminated;
    patch.interface = input.interface;
    candidate.measures.geometricAreaM2 += input.areaM2;

    if (patch.fast) {
      candidate.measures.fastAreaM2 += input.areaM2;
      auto jump = SolveObliqueFastShock(input.upstream, normal,
          input.shockNormalSpeedMPerS, options.gammaAdiabatic);
      if (!jump.ok()) {
        return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
            jump.status.code, jump.status.message);
      }
      patch.jump = jump.value;
    } else {
      patch.rejectionReason = "sub-fast";
    }
    if (patch.supercritical) candidate.measures.supercriticalAreaM2 += input.areaM2;

    patch.sourceEligibleBeforeClearance = input.sourceEnabled && patch.fast &&
        patch.supercritical && !input.sourceTerminated;
    if (patch.sourceEligibleBeforeClearance) {
      candidate.measures.counterfactualAreaM2 += input.areaM2;
      candidate.measures.counterfactualNumberRatePerS +=
          input.incidentNumberRatePerS;
      candidate.measures.counterfactualKineticEnergyRateW +=
          input.incidentKineticEnergyRateW;
    }
    if (patch.sourceEligibleBeforeClearance &&
        input.intersectsTransitionClearance) {
      patch.excludedByTransitionClearance = true;
      patch.rejectionReason = "transition-clearance";
      candidate.measures.excludedAreaM2 += input.areaM2;
      candidate.measures.excludedNumberRatePerS +=
          input.incidentNumberRatePerS;
      candidate.measures.excludedKineticEnergyRateW +=
          input.incidentKineticEnergyRateW;
    }
    patch.sourceActive = patch.sourceEligibleBeforeClearance &&
        !patch.excludedByTransitionClearance;
    if (patch.sourceActive) {
      patch.physicalNumberRatePerS = input.incidentNumberRatePerS;
      patch.physicalKineticEnergyRateW = input.incidentKineticEnergyRateW;
      candidate.measures.sourceActiveAreaM2 += input.areaM2;
    }
    if (input.sourceTerminated)
      candidate.measures.sourceTerminatedAreaM2 += input.areaM2;
    candidate.patches.push_back(std::move(patch));
  }

  const double fastFraction = candidate.measures.fastAreaM2 /
      candidate.measures.geometricAreaM2;
  if ((options.initialGate == InitialShockGate::AnyFastPatch &&
       candidate.measures.fastAreaM2 == 0.0) ||
      (options.initialGate == InitialShockGate::MinimumFastAreaFraction &&
       fastFraction < options.minimumFastAreaFraction)) {
    return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
        Core::StatusCode::InvalidState,
        "prepared surface does not satisfy the configured initial fast gate");
  }

  auto area = Ratio(candidate.measures.excludedAreaM2,
      candidate.measures.counterfactualAreaM2,
      options.maximumTransitionAreaFraction);
  auto number = Ratio(candidate.measures.excludedNumberRatePerS,
      candidate.measures.counterfactualNumberRatePerS,
      options.maximumTransitionNumberFraction);
  auto energy = Ratio(candidate.measures.excludedKineticEnergyRateW,
      candidate.measures.counterfactualKineticEnergyRateW,
      options.maximumTransitionEnergyFraction);
  const auto exceeds = [](const Core::Result<BudgetRatio>& ratio) {
    return ratio.ok() && ratio.value.validity == MeasureValidity::Valid &&
        !ratio.value.withinBound;
  };
  if (!(area.ok() && number.ok() && energy.ok()) || exceeds(area) ||
      exceeds(number) || exceeds(energy)) {
    return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Failure(
        Core::StatusCode::InvalidState,
        "transition-clearance exclusion exceeds a registered source budget");
  }
  candidate.measures.areaFractionValidity = area.value.validity;
  candidate.measures.numberFractionValidity = number.value.validity;
  candidate.measures.energyFractionValidity = energy.value.validity;

  std::shared_ptr<const ShockSurfaceSnapshot> publication =
      std::make_shared<const ShockSurfaceSnapshot>(std::move(candidate));
  current_ = publication;
  return Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>>::Success(
      std::move(publication));
}

} }  // namespace SEP::CoronalCME
