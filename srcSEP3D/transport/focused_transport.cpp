#include "focused_transport.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Transport {
namespace {

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double AdvanceDeterministicMu(double mu, double speedMPerS,
                              const FocusedLocalState& local,
                              double intervalS) {
  // Explicit midpoint is second order for one immutable local background.
  const double first = FocusedPitchDriftPerS(mu, speedMPerS, local);
  const double midpoint = ReflectPitchAngle(mu + 0.5 * intervalS * first);
  return ReflectPitchAngle(
      mu + intervalS * FocusedPitchDriftPerS(midpoint, speedMPerS, local));
}

}  // namespace

double RelativisticSpeed(double momentumKgMPerS, double massKg) {
  if (!std::isfinite(momentumKgMPerS) || momentumKgMPerS < 0.0 ||
      !std::isfinite(massKg) || massKg <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double mc = massKg * Core::Const::c;
  return Core::Const::c * momentumKgMPerS /
      std::hypot(momentumKgMPerS, mc);
}

double FocusedPitchDriftPerS(double mu, double speedMPerS,
                             const FocusedLocalState& local) {
  const double oneMinusMu2 = std::max(0.0, 1.0 - mu * mu);
  return 0.5 * oneMinusMu2 *
      (speedMPerS * local.divBhatPerM +
       mu * (local.divUPerS - 3.0 * local.fieldAlignedStrainPerS));
}

double FocusedLogMomentumRatePerS(double mu,
                                  const FocusedLocalState& local) {
  return -0.5 * ((1.0 - mu * mu) * local.divUPerS +
      (3.0 * mu * mu - 1.0) * local.fieldAlignedStrainPerS);
}

double ReflectPitchAngle(double mu, unsigned* reflections) {
  if (reflections != nullptr) *reflections = 0;
  if (!std::isfinite(mu) || (mu >= -1.0 && mu <= 1.0)) return mu;

  // The period-four triangular wave is equivalent to repeated specular
  // reflection but remains O(1) for an arbitrarily large stochastic kick.
  const double cell = std::floor((mu + 1.0) / 2.0);
  if (reflections != nullptr)
    *reflections = static_cast<unsigned>(std::fabs(cell));
  double folded = std::fmod(mu + 1.0, 4.0);
  if (folded < 0.0) folded += 4.0;
  return folded <= 2.0 ? folded - 1.0 : 3.0 - folded;
}

FocusedStepResult AdvanceFocused(const FocusedParticleState& initial,
                                  const FocusedLocalState& local,
                                  double massKg,
                                  double dtS,
                                  KeyedRandomStream* random) {
  FocusedStepResult result;
  result.state = initial;
  const double bNorm = local.bHat.Norm();
  if (!Finite(initial.positionM) ||
      !std::isfinite(initial.momentumKgMPerS) ||
      initial.momentumKgMPerS < 0.0 || !std::isfinite(initial.mu) ||
      initial.mu < -1.0 || initial.mu > 1.0 ||
      !Finite(local.bulkVelocityMPerS) || !Finite(local.bHat) ||
      std::fabs(bNorm - 1.0) > 1.0e-12 ||
      !std::isfinite(local.divBhatPerM) || !std::isfinite(local.divUPerS) ||
      !std::isfinite(local.fieldAlignedStrainPerS) ||
      !std::isfinite(local.dMuMuPerS) || local.dMuMuPerS < 0.0 ||
      !std::isfinite(local.dDmuMuDmuPerS) ||
      !std::isfinite(local.kappaPerpendicularM2PerS) ||
      local.kappaPerpendicularM2PerS != 0.0 ||
      !Finite(local.driftVelocityMPerS) ||
      local.driftVelocityMPerS.NormSq() != 0.0 ||
      !std::isfinite(massKg) || massKg <= 0.0 ||
      !std::isfinite(dtS) || dtS < 0.0) {
    result.status = Core::Status(
        Core::StatusCode::InvalidInput,
        "invalid focused state or nonzero reserved perpendicular/drift input");
    return result;
  }
  if (local.dMuMuPerS > 0.0 && random == nullptr) {
    result.status = Core::Status(
        Core::StatusCode::InvalidInput,
        "scattering step requires a keyed pitch-angle random stream");
    return result;
  }

  const double initialSpeed = RelativisticSpeed(
      initial.momentumKgMPerS, massKg);
  if (!std::isfinite(initialSpeed)) {
    result.status = Core::Status(Core::StatusCode::InvalidInput,
                                 "focused speed is invalid");
    return result;
  }

  unsigned reflections = 0;
  double muHalf = AdvanceDeterministicMu(
      initial.mu, initialSpeed, local, 0.5 * dtS);
  muHalf = ReflectPitchAngle(muHalf, &reflections);
  result.pitchReflections += reflections;

  const double firstLogP =
      0.5 * dtS * FocusedLogMomentumRatePerS(initial.mu, local);
  double momentumHalf = initial.momentumKgMPerS * std::exp(firstLogP);
  if (!std::isfinite(momentumHalf)) {
    result.status = Core::Status(Core::StatusCode::InvalidInput,
                                 "focused half-step momentum overflowed");
    return result;
  }

  const double normal = local.dMuMuPerS > 0.0 ? random->Normal01() : 0.0;
  const double dW = std::sqrt(dtS) * normal;
  result.deterministicMuIncrement = local.dDmuMuDmuPerS * dtS;
  result.stochasticMuIncrement = local.dMuMuPerS > 0.0
      ? std::sqrt(2.0 * local.dMuMuPerS) * dW : 0.0;
  if (local.scheme == PitchAngleScheme::ReflectingMilstein &&
      local.dMuMuPerS > 0.0) {
    // b=sqrt(2D) gives 0.5*b*b' = 0.5*dD/dmu.
    result.stochasticMuIncrement += 0.5 * local.dDmuMuDmuPerS *
        (dW * dW - dtS);
  }
  double muKicked = muHalf + result.deterministicMuIncrement +
                    result.stochasticMuIncrement;
  muKicked = ReflectPitchAngle(muKicked, &reflections);
  result.pitchReflections += reflections;

  const double secondLogP =
      0.5 * dtS * FocusedLogMomentumRatePerS(muKicked, local);
  const double finalMomentum = momentumHalf * std::exp(secondLogP);
  const double finalSpeed = RelativisticSpeed(finalMomentum, massKg);
  if (!std::isfinite(finalMomentum) || !std::isfinite(finalSpeed)) {
    result.status = Core::Status(Core::StatusCode::InvalidInput,
                                 "focused final momentum is invalid");
    return result;
  }
  result.logMomentumIncrement = firstLogP + secondLogP;

  double muFinal = AdvanceDeterministicMu(
      muKicked, finalSpeed, local, 0.5 * dtS);
  muFinal = ReflectPitchAngle(muFinal, &reflections);
  result.pitchReflections += reflections;

  const double midpointMu = 0.5 * (initial.mu + muFinal);
  const double midpointSpeed = 0.5 * (initialSpeed + finalSpeed);
  result.displacementM =
      (local.bulkVelocityMPerS + local.bHat * (midpointMu * midpointSpeed)) *
      dtS;
  result.state.positionM += result.displacementM;
  result.state.momentumKgMPerS = finalMomentum;
  result.state.mu = muFinal;
  if (!Finite(result.state.positionM) || !std::isfinite(result.state.mu) ||
      result.state.mu < -1.0 || result.state.mu > 1.0) {
    result.status = Core::Status(Core::StatusCode::InvalidInput,
                                 "focused step produced an invalid state");
    return result;
  }
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Transport
}  // namespace SEP3D
