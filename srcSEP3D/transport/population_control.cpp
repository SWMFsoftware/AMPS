#include "population_control.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Transport {
namespace {

Core::Status Invalid(const char* message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double RelativeResidual(double observed, double reference) {
  return std::fabs(observed - reference) /
      std::max(std::fabs(reference), std::numeric_limits<double>::min());
}

}  // namespace

double RelativisticTotalEnergyJ(const Core::Vec3& momentumKgMPerS,
                                double massKg) {
  if (!Finite(momentumKgMPerS) || !std::isfinite(massKg) || massKg <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double pc = momentumKgMPerS.Norm() * Core::Const::c;
  const double rest = massKg * Core::Const::c * Core::Const::c;
  return std::hypot(pc, rest);
}

RelativisticMergeResult MergeRelativisticThreeToTwo(
    const std::array<WeightedPhasePoint, 3>& input,
    double massKg, const Core::Vec3& direction) {
  RelativisticMergeResult result;
  if (!std::isfinite(massKg) || massKg <= 0.0 || !Finite(direction) ||
      std::fabs(direction.Norm() - 1.0) > 1.0e-12) {
    result.status = Invalid("relativistic merge mass or direction is invalid");
    return result;
  }

  double totalWeight = 0.0;
  double totalEnergy = 0.0;
  double momentumScale = 0.0;
  Core::Vec3 totalMomentum;
  Core::Vec3 weightedPosition;
  double largestMomentum = 0.0;
  for (const WeightedPhasePoint& point : input) {
    if (!std::isfinite(point.weight) || point.weight <= 0.0 ||
        !Finite(point.positionM) || !Finite(point.momentumKgMPerS)) {
      result.status = Invalid("relativistic merge phase point is invalid");
      return result;
    }
    const double energy = RelativisticTotalEnergyJ(
        point.momentumKgMPerS, massKg);
    if (!std::isfinite(energy)) {
      result.status = Invalid("relativistic merge input energy is invalid");
      return result;
    }
    totalWeight += point.weight;
    totalMomentum += point.weight * point.momentumKgMPerS;
    weightedPosition += point.weight * point.positionM;
    totalEnergy += point.weight * energy;
    momentumScale += point.weight * point.momentumKgMPerS.Norm();
    largestMomentum = std::max(largestMomentum,
                               point.momentumKgMPerS.Norm());
  }
  if (!std::isfinite(totalWeight) || totalWeight <= 0.0 ||
      !std::isfinite(totalEnergy) || totalEnergy <= 0.0) {
    result.status = Invalid("relativistic merge moments overflowed");
    return result;
  }

  const Core::Vec3 meanMomentum = totalMomentum / totalWeight;
  result.outputPositionM = weightedPosition / totalWeight;
  result.outputWeight = 0.5 * totalWeight;
  const double targetPairEnergy = 2.0 * totalEnergy / totalWeight;
  auto pairEnergy = [&](double q) {
    return RelativisticTotalEnergyJ(meanMomentum + q * direction, massKg) +
           RelativisticTotalEnergyJ(meanMomentum - q * direction, massKg);
  };

  double qLower = 0.0;
  double qUpper = std::max(massKg * Core::Const::c, largestMomentum);
  if (!(qUpper > 0.0)) qUpper = std::numeric_limits<double>::min();
  const double minimumPairEnergy = pairEnergy(0.0);
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
      targetPairEnergy;
  if (!std::isfinite(minimumPairEnergy) ||
      targetPairEnergy + tolerance < minimumPairEnergy) {
    result.status = Invalid(
        "relativistic merge violates the convex energy lower bound");
    return result;
  }
  if (targetPairEnergy > minimumPairEnergy + tolerance) {
    unsigned expansion = 0;
    while (pairEnergy(qUpper) < targetPairEnergy && expansion++ < 256)
      qUpper *= 2.0;
    if (!std::isfinite(qUpper) || pairEnergy(qUpper) < targetPairEnergy) {
      result.status = Invalid(
          "relativistic merge could not bracket its energy root");
      return result;
    }
    // A fixed iteration count makes the result reproducible across tolerance
    // library implementations while driving the bracket below double ulps.
    for (unsigned iteration = 0; iteration < 160; ++iteration) {
      const double middle = 0.5 * (qLower + qUpper);
      if (pairEnergy(middle) < targetPairEnergy)
        qLower = middle;
      else
        qUpper = middle;
    }
  } else {
    qUpper = 0.0;
  }

  const double q = 0.5 * (qLower + qUpper);
  result.firstMomentumKgMPerS = meanMomentum + q * direction;
  result.secondMomentumKgMPerS = meanMomentum - q * direction;
  const Core::Vec3 outputMomentum = result.outputWeight *
      (result.firstMomentumKgMPerS + result.secondMomentumKgMPerS);
  const double outputEnergy = result.outputWeight *
      (RelativisticTotalEnergyJ(result.firstMomentumKgMPerS, massKg) +
       RelativisticTotalEnergyJ(result.secondMomentumKgMPerS, massKg));
  result.relativeWeightResidual = RelativeResidual(
      2.0 * result.outputWeight, totalWeight);
  result.relativeMomentumResidual =
      (outputMomentum - totalMomentum).Norm() /
      std::max(momentumScale, std::numeric_limits<double>::min());
  result.relativeEnergyResidual = RelativeResidual(outputEnergy, totalEnergy);
  if (!Finite(result.firstMomentumKgMPerS) ||
      !Finite(result.secondMomentumKgMPerS) ||
      !std::isfinite(result.relativeWeightResidual) ||
      !std::isfinite(result.relativeMomentumResidual) ||
      !std::isfinite(result.relativeEnergyResidual)) {
    result.status = Invalid("relativistic merge produced invalid output");
    return result;
  }
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Transport
}  // namespace SEP3D
