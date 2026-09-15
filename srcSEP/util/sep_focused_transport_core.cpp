#include "sep_focused_transport_core.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace Transport {
namespace {

double DeterministicMuRate(double mu, double speedMPerS,
                           const FocusedTransportBackground& background) {
  return FocusedPitchDriftPerS(mu, speedMPerS, background);
}

ScalarResult ApplyFocusedMomentum(double momentumKgMPerS, double mu,
                                  const FocusedTransportBackground& background,
                                  double dtS) {
  ScalarResult result;
  const double rate = FocusedLogMomentumRatePerS(mu, background);
  const double exponent = rate * dtS;
  if (!std::isfinite(momentumKgMPerS) || momentumKgMPerS < 0.0 ||
      !std::isfinite(exponent) ||
      exponent > std::log(std::numeric_limits<double>::max())) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
        "focused momentum update overflowed or received invalid state");
    return result;
  }
  result.value = momentumKgMPerS * std::exp(exponent);
  result.status = std::isfinite(result.value)
      ? Status::Ok()
      : Status::Error(StatusCode::InvalidParticleState,
                      "focused momentum result is not finite");
  return result;
}

double AdvanceDeterministicMu(double mu, double speedMPerS,
                              const FocusedTransportBackground& background,
                              double intervalS) {
  // Explicit midpoint is second order for a frozen background snapshot.  It is
  // used for each half of the symmetric split so focusing does not collapse to
  // two first-order half-Euler updates when D_mumu is zero.
  const double first = DeterministicMuRate(mu, speedMPerS, background);
  const double midpoint = ReflectPitchAngle(mu + 0.5 * intervalS * first, NULL);
  return ReflectPitchAngle(
      mu + intervalS * DeterministicMuRate(midpoint, speedMPerS, background),
      NULL);
}

}  // namespace

double FocusedPitchDriftPerS(
    double mu, double speedMPerS,
    const FocusedTransportBackground& background) {
  const double oneMinusMu2 = std::max(0.0, 1.0 - mu * mu);
  const double focusing = -0.5 * oneMinusMu2 * speedMPerS *
                          background.dLnAbsBdsPerM;
  if (background.equationMode ==
      FocusedEquationMode::ReducedFieldAligned1D) {
    // Compatibility form used by pre-WP07 campaigns.  It is available only
    // behind the explicit mode and is never silently selected from geometry.
    return focusing - 0.5 * mu * oneMinusMu2 *
        background.parallelVelocityGradientPerS;
  }

  // Full gyrotropic coefficient in the local plasma frame:
  // dmu/dt=(1-mu^2)/2[-v dln|B|/ds + mu(divU-3 bb:gradU)].
  return focusing + 0.5 * mu * oneMinusMu2 *
      (background.velocityDivergencePerS -
       3.0 * background.fieldAlignedStrainPerS);
}

double FocusedLogMomentumRatePerS(
    double mu, const FocusedTransportBackground& background) {
  if (background.equationMode ==
      FocusedEquationMode::ReducedFieldAligned1D) {
    return -background.velocityDivergencePerS / 3.0;
  }
  // dp/dt=-p/2[(1-mu^2) divU +(3mu^2-1) bb:gradU].  This
  // pitch-dependent form reduces to Parker cooling after isotropic averaging.
  return -0.5 * ((1.0 - mu * mu) * background.velocityDivergencePerS +
      (3.0 * mu * mu - 1.0) * background.fieldAlignedStrainPerS);
}

void ThreadLocalWaveAccumulator::Deposit(
    const WaveContribution& contribution) {
  if (!contribution.turbulenceStateIdentity.empty() &&
      std::isfinite(contribution.signedStreaming)) {
    contributions_.push_back(contribution);
  }
}

std::vector<WaveContribution> ThreadLocalWaveAccumulator::DeterministicReduce(
    const std::vector<const ThreadLocalWaveAccumulator*>& locals) {
  std::vector<WaveContribution> reduced;
  for (const ThreadLocalWaveAccumulator* local : locals) {
    if (!local) continue;
    reduced.insert(reduced.end(), local->contributions_.begin(),
                   local->contributions_.end());
  }
  std::sort(reduced.begin(), reduced.end(),
            [](const WaveContribution& left, const WaveContribution& right) {
              if (left.turbulenceStateIdentity != right.turbulenceStateIdentity)
                return left.turbulenceStateIdentity < right.turbulenceStateIdentity;
              return left.signedStreaming < right.signedStreaming;
            });
  return reduced;
}

double ReflectPitchAngle(double mu, unsigned* reflections) {
  if (reflections) *reflections = 0;
  if (!std::isfinite(mu)) return mu;
  if (mu >= -1.0 && mu <= 1.0) return mu;

  const double cell = std::floor((mu + 1.0) / 2.0);
  if (reflections) *reflections = static_cast<unsigned>(std::fabs(cell));
  double folded = std::fmod(mu + 1.0, 4.0);
  if (folded < 0.0) folded += 4.0;
  return folded <= 2.0 ? folded - 1.0 : 3.0 - folded;
}

FocusedTransportIncrement AdvanceFocusedTransportDmumu(
    const FocusedTransportState& initial,
    const FocusedTransportBackground& background,
    double massKg,
    double speedOfLightMPerS,
    double dtS,
    const PitchAngleDiffusionProvider& provider,
    RandomStream& random,
    ThreadLocalWaveAccumulator* waveAccumulator) {
  FocusedTransportIncrement result;
  result.state = initial;
  if (!std::isfinite(initial.arcLengthM) ||
      !std::isfinite(initial.momentumKgMPerS) ||
      initial.momentumKgMPerS < 0.0 || !std::isfinite(initial.mu) ||
      initial.mu < -1.0 || initial.mu > 1.0 ||
      !std::isfinite(massKg) || massKg <= 0.0 ||
      !std::isfinite(speedOfLightMPerS) || speedOfLightMPerS <= 0.0 ||
      !std::isfinite(dtS) || dtS < 0.0 ||
      !std::isfinite(background.dLnAbsBdsPerM) ||
      !std::isfinite(background.plasmaAdvectionMPerS) ||
      !std::isfinite(background.parallelVelocityGradientPerS) ||
      !std::isfinite(background.velocityDivergencePerS) ||
      !std::isfinite(background.fieldAlignedStrainPerS)) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "invalid focused-transport input state");
    return result;
  }

  const ScalarResult initialSpeed = SpeedFromMomentum(
      initial.momentumKgMPerS, massKg, speedOfLightMPerS);
  if (!initialSpeed.status.ok()) {
    result.status = initialSpeed.status;
    return result;
  }

  // First half of the deterministic operator.
  double muHalf = AdvanceDeterministicMu(
      initial.mu, initialSpeed.value, background, 0.5 * dtS);
  unsigned reflectionCount = 0;
  muHalf = ReflectPitchAngle(muHalf, &reflectionCount);
  result.pitchAngleReflections += reflectionCount;

  const ScalarResult halfMomentum = ApplyFocusedMomentum(
      initial.momentumKgMPerS, initial.mu, background, 0.5 * dtS);
  if (!halfMomentum.status.ok()) {
    result.status = halfMomentum.status;
    return result;
  }

  // Coefficients are sampled at a predicted physical midpoint, not at a
  // hard-coded adapter coordinate.  This is second-order accurate for smooth
  // profiles and ensures crossing a segment changes the sampled turbulence.
  const double predictedMidpointM = initial.arcLengthM + 0.5 * dtS *
      (background.plasmaAdvectionMPerS + initialSpeed.value * muHalf);
  const PitchAngleDiffusionSample coefficient = provider.Evaluate(
      predictedMidpointM, halfMomentum.value, muHalf);
  if (!coefficient.status.ok() || !std::isfinite(coefficient.dMuMuPerS) ||
      !std::isfinite(coefficient.dDmuMuDmuPerS) ||
      coefficient.dMuMuPerS < 0.0 || coefficient.provenance.empty() ||
      coefficient.turbulenceStateIdentity.empty()) {
    result.status = Status::Error(StatusCode::InvalidCoefficient,
        coefficient.status.message.empty()
            ? "invalid pitch-angle-diffusion provider result"
            : coefficient.status.message);
    return result;
  }
  result.coefficientProvenance = coefficient.provenance;

  // Ito drift includes dD_mumu/dmu.  The square-root amplitude has the units
  // required for a dimensionless pitch-angle increment because D_mumu is s^-1.
  const double normal = coefficient.dMuMuPerS == 0.0
      ? 0.0 : random.Normal01();
  const double dW = std::sqrt(dtS) * normal;
  result.deterministicMuIncrement = coefficient.dDmuMuDmuPerS * dtS;
  result.stochasticMuIncrement = coefficient.dMuMuPerS == 0.0
      ? 0.0 : std::sqrt(2.0 * coefficient.dMuMuPerS) * dW;

  if (background.pitchAngleScheme ==
      BoundedPitchAngleScheme::ReflectingMilstein &&
      coefficient.dMuMuPerS > 0.0) {
    // For b(mu)=sqrt(2D), the scalar Milstein correction is
    // 0.5*b*b'*(dW^2-dt)=0.5*D'*(dW^2-dt).  Adding it to the Ito drift D'dt
    // yields the requested higher weak-order treatment without evaluating D
    // at a stochastic predictor that might lie outside [-1,1].  The subsequent
    // mirror map is the discrete realization of the zero-normal-flux boundary.
    result.stochasticMuIncrement +=
        0.5 * coefficient.dDmuMuDmuPerS * (dW * dW - dtS);
  }
  double muKicked = muHalf + result.deterministicMuIncrement +
                    result.stochasticMuIncrement;
  muKicked = ReflectPitchAngle(muKicked, &reflectionCount);
  result.pitchAngleReflections += reflectionCount;

  const ScalarResult finalMomentum = ApplyFocusedMomentum(
      halfMomentum.value, muKicked, background, 0.5 * dtS);
  if (!finalMomentum.status.ok()) {
    result.status = finalMomentum.status;
    return result;
  }
  const ScalarResult finalSpeed = SpeedFromMomentum(
      finalMomentum.value, massKg, speedOfLightMPerS);
  if (!finalSpeed.status.ok()) {
    result.status = finalSpeed.status;
    return result;
  }

  double muFinal = AdvanceDeterministicMu(
      muKicked, finalSpeed.value, background, 0.5 * dtS);
  muFinal = ReflectPitchAngle(muFinal, &reflectionCount);
  result.pitchAngleReflections += reflectionCount;

  const double midpointMu = 0.5 * (initial.mu + muFinal);
  result.displacementM =
      (background.plasmaAdvectionMPerS + finalSpeed.value * midpointMu) * dtS;
  result.state.arcLengthM += result.displacementM;
  result.state.momentumKgMPerS = finalMomentum.value;
  result.state.mu = muFinal;
  if (!std::isfinite(result.state.arcLengthM) ||
      !std::isfinite(result.state.mu)) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "focused-transport increment is not finite");
    return result;
  }

  if (waveAccumulator) {
    WaveContribution contribution;
    contribution.turbulenceStateIdentity =
        coefficient.turbulenceStateIdentity;
    contribution.signedStreaming = midpointMu * result.displacementM;
    contribution.intervalS = dtS;
    contribution.displacementM = result.displacementM;
    waveAccumulator->Deposit(contribution);
  }
  result.status = Status::Ok();
  return result;
}

}  // namespace Transport
}  // namespace SEP
