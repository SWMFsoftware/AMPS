#include "sep_focused_transport_mfp_core.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace Transport {
namespace {

double DeterministicMuRate(double mu, double speedMPerS,
                           const FocusedTransportMfpBackground& background) {
  FocusedTransportBackground shared;
  shared.dLnAbsBdsPerM = background.dLnAbsBdsPerM;
  shared.parallelVelocityGradientPerS =
      background.parallelVelocityGradientPerS;
  shared.velocityDivergencePerS = background.velocityDivergencePerS;
  shared.fieldAlignedStrainPerS = background.fieldAlignedStrainPerS;
  shared.equationMode = background.equationMode;
  return FocusedPitchDriftPerS(mu, speedMPerS, shared);
}

double AdvanceMuMidpoint(double mu, double speedMPerS,
                         const FocusedTransportMfpBackground& background,
                         double dtS) {
  const double first = DeterministicMuRate(mu, speedMPerS, background);
  const double midpoint = ReflectPitchAngle(mu + 0.5 * dtS * first, NULL);
  return ReflectPitchAngle(
      mu + dtS * DeterministicMuRate(midpoint, speedMPerS, background), NULL);
}

bool ValidMomentumRange(const MeanFreePathSample& sample,
                        double momentumKgMPerS) {
  if (sample.minimumMomentumKgMPerS == 0.0 &&
      sample.maximumMomentumKgMPerS == 0.0) return true;
  return std::isfinite(sample.minimumMomentumKgMPerS) &&
         std::isfinite(sample.maximumMomentumKgMPerS) &&
         sample.minimumMomentumKgMPerS >= 0.0 &&
         sample.maximumMomentumKgMPerS >= sample.minimumMomentumKgMPerS &&
         momentumKgMPerS >= sample.minimumMomentumKgMPerS &&
         momentumKgMPerS <= sample.maximumMomentumKgMPerS;
}

double LogMomentumRate(double mu,
                       const FocusedTransportMfpBackground& background) {
  FocusedTransportBackground shared;
  shared.velocityDivergencePerS = background.velocityDivergencePerS;
  shared.fieldAlignedStrainPerS = background.fieldAlignedStrainPerS;
  shared.equationMode = background.equationMode;
  return FocusedLogMomentumRatePerS(mu, shared);
}

struct EventRates {
  Status status;
  double plusPerS = 0.0;
  double minusPerS = 0.0;
  double totalPerS = 0.0;
};

EventRates RatesFromSample(const MeanFreePathSample& sample,
                           double speedMPerS) {
  EventRates rates;
  if (sample.hasBranchResolvedRates) {
    if (!std::isfinite(sample.nuPlusPerS) || sample.nuPlusPerS < 0.0 ||
        !std::isfinite(sample.nuMinusPerS) || sample.nuMinusPerS < 0.0) {
      rates.status = Status::Error(StatusCode::InvalidCoefficient,
          "branch-resolved event rates must be finite and nonnegative");
      return rates;
    }
    rates.plusPerS = sample.nuPlusPerS;
    rates.minusPerS = sample.nuMinusPerS;
  }
  else {
    // Compatibility mapping for providers that publish only lambda_parallel.
    // Splitting the total rate equally is explicit and unbiased; production
    // self-consistent turbulence providers set hasBranchResolvedRates=true.
    const double total = std::isinf(sample.lambdaParallelM)
        ? 0.0 : speedMPerS / sample.lambdaParallelM;
    rates.plusPerS = 0.5 * total;
    rates.minusPerS = 0.5 * total;
  }
  rates.totalPerS = rates.plusPerS + rates.minusPerS;
  rates.status = std::isfinite(rates.totalPerS)
      ? Status::Ok()
      : Status::Error(StatusCode::InvalidCoefficient,
                      "total scattering rate is not finite");
  return rates;
}

double IntegratedLinearHazard(double rate0PerS, double rate1PerS,
                              double fullIntervalS, double tS) {
  if (fullIntervalS == 0.0) return 0.0;
  const double slope = (rate1PerS - rate0PerS) / fullIntervalS;
  return rate0PerS * tS + 0.5 * slope * tS * tS;
}

double LocateHazardEvent(double targetOpticalDepth, double rate0PerS,
                         double rate1PerS, double intervalS,
                         std::size_t* iterations) {
  // Monotone bisection is deliberately used instead of the quadratic formula:
  // it remains stable as the rate gradient tends to zero and cannot choose the
  // wrong algebraic root for a decreasing but nonnegative rate profile.
  double left = 0.0, right = intervalS;
  for (std::size_t i = 0; i < 64; ++i) {
    const double middle = 0.5 * (left + right);
    if (IntegratedLinearHazard(rate0PerS, rate1PerS, intervalS, middle) <
        targetOpticalDepth)
      left = middle;
    else
      right = middle;
    if (iterations) ++(*iterations);
  }
  return 0.5 * (left + right);
}

}  // namespace

ScalarResult SampleExponentialWaitingTime(double ratePerS,
                                          RandomStream& random) {
  ScalarResult result;
  if (!std::isfinite(ratePerS) || ratePerS < 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "event rate must be finite and nonnegative");
    return result;
  }
  if (ratePerS == 0.0) {
    result.value = std::numeric_limits<double>::infinity();
    result.status = Status::Ok();
    return result;
  }
  const double xi = random.UniformOpen01();
  if (!std::isfinite(xi) || xi <= 0.0 || xi >= 1.0) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
        "random stream violated its open-unit-interval contract");
    return result;
  }
  result.value = -std::log(xi) / ratePerS;
  result.status = std::isfinite(result.value) && result.value > 0.0
      ? Status::Ok()
      : Status::Error(StatusCode::InvalidParticleState,
                      "exponential waiting time is not finite and positive");
  return result;
}

WaveFrameScatterResult ScatterIsotropicallyInWaveFrame(
    double speedMPerS, double mu, double waveSpeedMPerS,
    double speedOfLightMPerS, RandomStream& random) {
  WaveFrameScatterResult result;
  if (!std::isfinite(speedMPerS) || speedMPerS < 0.0 ||
      speedMPerS >= speedOfLightMPerS || !std::isfinite(mu) ||
      mu < -1.0 || mu > 1.0 || !std::isfinite(waveSpeedMPerS) ||
      std::fabs(waveSpeedMPerS) >= speedOfLightMPerS ||
      !std::isfinite(speedOfLightMPerS) || speedOfLightMPerS <= 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "invalid plasma/wave-frame scattering state");
    return result;
  }

  const double betaParallel = speedMPerS * mu / speedOfLightMPerS;
  const double betaPerpendicular = speedMPerS *
      std::sqrt(std::max(0.0, 1.0 - mu * mu)) / speedOfLightMPerS;
  const double betaWave = waveSpeedMPerS / speedOfLightMPerS;
  const double gammaWave = 1.0 / std::sqrt(1.0 - betaWave * betaWave);
  const double toWaveDenominator = 1.0 - betaParallel * betaWave;
  if (!(toWaveDenominator > 0.0) || !std::isfinite(toWaveDenominator)) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "plasma-to-wave Lorentz denominator is invalid");
    return result;
  }

  const double betaParallelWave =
      (betaParallel - betaWave) / toWaveDenominator;
  const double betaPerpendicularWave =
      betaPerpendicular / (gammaWave * toWaveDenominator);
  const double betaWaveFrame = std::hypot(betaParallelWave,
                                          betaPerpendicularWave);
  if (!std::isfinite(betaWaveFrame) || betaWaveFrame >= 1.0) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "wave-frame particle speed is invalid");
    return result;
  }
  result.waveFrameSpeedBeforeMPerS = betaWaveFrame * speedOfLightMPerS;

  // Uniform mu in [-1,1] is isotropic in the wave frame. Endpoints cannot be
  // drawn because the underlying stream is open on (0,1).
  const double scatteredMuWave = -1.0 + 2.0 * random.UniformOpen01();
  const double scatteredParallelWave = betaWaveFrame * scatteredMuWave;
  const double scatteredPerpendicularWave = betaWaveFrame *
      std::sqrt(std::max(0.0, 1.0 - scatteredMuWave * scatteredMuWave));
  result.waveFrameSpeedAfterMPerS = std::hypot(
      scatteredParallelWave, scatteredPerpendicularWave) * speedOfLightMPerS;

  const double toPlasmaDenominator =
      1.0 + scatteredParallelWave * betaWave;
  if (!(toPlasmaDenominator > 0.0) ||
      !std::isfinite(toPlasmaDenominator)) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "wave-to-plasma Lorentz denominator is invalid");
    return result;
  }
  const double scatteredParallel =
      (scatteredParallelWave + betaWave) / toPlasmaDenominator;
  const double scatteredPerpendicular = scatteredPerpendicularWave /
      (gammaWave * toPlasmaDenominator);
  const double betaFinal = std::hypot(scatteredParallel,
                                      scatteredPerpendicular);
  if (!std::isfinite(betaFinal) || betaFinal >= 1.0 || betaFinal == 0.0) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "scattered plasma-frame speed is invalid");
    return result;
  }
  result.speedMPerS = betaFinal * speedOfLightMPerS;
  result.mu = scatteredParallel / betaFinal;
  result.status = Status::Ok();
  return result;
}

FocusedTransportMfpIncrement AdvanceFocusedTransportMfp(
    const FocusedTransportMfpState& initial,
    const FocusedTransportMfpBackground& background,
    double massKg, double speedOfLightMPerS, double dtS,
    double maximumDeterministicIntervalS,
    const MeanFreePathProvider& provider, RandomStream& random,
    ThreadLocalWaveAccumulator* waveAccumulator) {
  FocusedTransportMfpIncrement result;
  result.state = initial;
  if (!std::isfinite(initial.arcLengthM) ||
      !std::isfinite(initial.momentumKgMPerS) ||
      initial.momentumKgMPerS < 0.0 || !std::isfinite(initial.mu) ||
      initial.mu < -1.0 || initial.mu > 1.0 || !std::isfinite(massKg) ||
      massKg <= 0.0 || !std::isfinite(speedOfLightMPerS) ||
      speedOfLightMPerS <= 0.0 || !std::isfinite(dtS) || dtS < 0.0 ||
      !std::isfinite(maximumDeterministicIntervalS) ||
      maximumDeterministicIntervalS <= 0.0 ||
      !std::isfinite(background.dLnAbsBdsPerM) ||
      !std::isfinite(background.plasmaAdvectionMPerS) ||
      !std::isfinite(background.parallelVelocityGradientPerS) ||
      !std::isfinite(background.velocityDivergencePerS) ||
      !std::isfinite(background.fieldAlignedStrainPerS) ||
      !std::isfinite(background.alfvenSpeedMPerS) ||
      std::fabs(background.alfvenSpeedMPerS) >= speedOfLightMPerS ||
      (!std::isnan(initial.remainingOpticalDepth) &&
       (!std::isfinite(initial.remainingOpticalDepth) ||
        initial.remainingOpticalDepth <= 0.0))) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "invalid event-driven transport input");
    return result;
  }

  if (std::isnan(result.state.remainingOpticalDepth)) {
    const double xi = random.UniformOpen01();
    if (!std::isfinite(xi) || xi <= 0.0 || xi >= 1.0) {
      result.status = Status::Error(StatusCode::InvalidParticleState,
          "random stream violated its open-unit-interval contract");
      return result;
    }
    result.state.remainingOpticalDepth = -std::log(xi);
  }

  double elapsedS = 0.0;
  while (elapsedS < dtS) {
    if (result.diagnostics.deterministicIntervals >= 1000000U) {
      result.status = Status::Error(StatusCode::StepUnderflow,
          "event-driven mover exceeded its declared interval budget");
      return result;
    }
    const ScalarResult speed = SpeedFromMomentum(
        result.state.momentumKgMPerS, massKg, speedOfLightMPerS);
    if (!speed.status.ok()) {
      result.status = speed.status;
      return result;
    }
    const MeanFreePathSample sample0 = provider.Evaluate(
        result.state.arcLengthM, result.state.momentumKgMPerS,
        result.state.mu);
    result.diagnostics.providerEvaluations++;
    const bool lambdaValid =
        (std::isfinite(sample0.lambdaParallelM) ||
         std::isinf(sample0.lambdaParallelM)) &&
        sample0.lambdaParallelM > 0.0;
    if (!sample0.status.ok() || !lambdaValid ||
        !ValidMomentumRange(sample0, result.state.momentumKgMPerS) ||
        sample0.provenance.empty() || sample0.turbulenceStateIdentity.empty()) {
      result.status = Status::Error(StatusCode::InvalidCoefficient,
          sample0.status.message.empty()
              ? "invalid mean-free-path provider result"
              : sample0.status.message);
      return result;
    }
    result.coefficientProvenance = sample0.provenance;
    const EventRates rates0 = RatesFromSample(sample0, speed.value);
    if (!rates0.status.ok()) {
      result.status = rates0.status;
      return result;
    }

    const double remainingS = dtS - elapsedS;
    const double trialIntervalS =
        std::min(remainingS, maximumDeterministicIntervalS);
    if (!std::isfinite(trialIntervalS) || trialIntervalS <= 0.0) {
      result.status = Status::Error(StatusCode::StepUnderflow,
                                    "event interval is not finite and positive");
      return result;
    }

    const double muInitial = result.state.mu;
    const double trialMu = AdvanceMuMidpoint(
        muInitial, speed.value, background, trialIntervalS);
    const double logMomentumTrial = LogMomentumRate(
        0.5 * (muInitial + trialMu), background) * trialIntervalS;
    const double trialMomentum = result.state.momentumKgMPerS *
        std::exp(logMomentumTrial);
    const ScalarResult trialSpeed = SpeedFromMomentum(
        trialMomentum, massKg, speedOfLightMPerS);
    if (!trialSpeed.status.ok()) {
      result.status = trialSpeed.status;
      return result;
    }
    const double trialMidpointMu = 0.5 * (muInitial + trialMu);
    const double trialDisplacement =
        (background.plasmaAdvectionMPerS +
         0.5 * (speed.value + trialSpeed.value) * trialMidpointMu) *
        trialIntervalS;

    // A second sample at the predicted end supplies a linear-in-time hazard
    // model.  This is exact for piecewise-linear rates and converges under the
    // same deterministic limiter used for geometry/background variation.
    const MeanFreePathSample sample1 = provider.Evaluate(
        result.state.arcLengthM + trialDisplacement, trialMomentum, trialMu);
    result.diagnostics.providerEvaluations++;
    if (!sample1.status.ok() ||
        !ValidMomentumRange(sample1, trialMomentum) ||
        sample1.turbulenceStateIdentity != sample0.turbulenceStateIdentity) {
      result.status = Status::Error(StatusCode::InvalidCoefficient,
          "end-of-interval coefficient sample is invalid or changed snapshot");
      return result;
    }
    const EventRates rates1 = RatesFromSample(sample1, trialSpeed.value);
    if (!rates1.status.ok()) {
      result.status = rates1.status;
      return result;
    }
    const double trialHazard = IntegratedLinearHazard(
        rates0.totalPerS, rates1.totalPerS, trialIntervalS, trialIntervalS);
    const bool eventOccurs = trialHazard >= result.state.remainingOpticalDepth &&
        trialHazard > 0.0;
    const double intervalS = eventOccurs
        ? LocateHazardEvent(result.state.remainingOpticalDepth,
                            rates0.totalPerS, rates1.totalPerS,
                            trialIntervalS,
                            &result.diagnostics.opticalDepthRootIterations)
        : trialIntervalS;

    const double muFinal = AdvanceMuMidpoint(
        muInitial, speed.value, background, intervalS);
    const double midpointMu = 0.5 * (muInitial + muFinal);
    const double finalMomentum = result.state.momentumKgMPerS * std::exp(
        LogMomentumRate(midpointMu, background) * intervalS);
    const ScalarResult finalSpeed = SpeedFromMomentum(
        finalMomentum, massKg, speedOfLightMPerS);
    if (!finalSpeed.status.ok()) {
      result.status = finalSpeed.status;
      return result;
    }
    const double displacement =
        (background.plasmaAdvectionMPerS +
         0.5 * (speed.value + finalSpeed.value) * midpointMu) * intervalS;
    if (!std::isfinite(displacement)) {
      result.status = Status::Error(StatusCode::InvalidParticleState,
                                    "event interval displacement is invalid");
      return result;
    }
    result.state.arcLengthM += displacement;
    result.state.momentumKgMPerS = finalMomentum;
    result.state.mu = muFinal;
    result.displacementM += displacement;
    result.diagnostics.deterministicIntervals++;
    if (rates0.totalPerS == 0.0 && rates1.totalPerS == 0.0)
      result.diagnostics.ballisticIntervals++;

    // Stage the complete interval record locally.  It is deposited only after
    // the event/no-event branch below succeeds, so a failed Lorentz transform
    // or invalid next optical-depth draw cannot publish work from a particle
    // update that the caller must reject transactionally.
    WaveContribution contribution;
    contribution.turbulenceStateIdentity = sample0.turbulenceStateIdentity;
    contribution.signedStreaming = midpointMu * displacement;
    contribution.intervalS = intervalS;
    contribution.displacementM = displacement;
    contribution.eventIndex = result.state.nextEventIndex;
    contribution.scatteringEventAtEnd = eventOccurs;

    elapsedS += intervalS;
    if (eventOccurs) {
      const double fraction = intervalS / trialIntervalS;
      const double plusAtEvent = rates0.plusPerS +
          fraction * (rates1.plusPerS - rates0.plusPerS);
      const double minusAtEvent = rates0.minusPerS +
          fraction * (rates1.minusPerS - rates0.minusPerS);
      const double totalAtEvent = plusAtEvent + minusAtEvent;
      if (!(totalAtEvent > 0.0)) {
        result.status = Status::Error(StatusCode::InvalidCoefficient,
            "hazard root reached an event with zero branch rate");
        return result;
      }

      // Conditional branch selection uses nu+/nu- at the event.  A branch
      // whose physical rate is zero has exactly zero selection probability.
      const bool plusBranch = random.UniformOpen01() <
          plusAtEvent / totalAtEvent;
      const double branchSign = plusBranch ? 1.0 : -1.0;
      contribution.resonantBranch = plusBranch ? 1 : -1;
      contribution.preWaveMomentumKgMPerS = finalMomentum;
      const WaveFrameScatterResult scattered =
          ScatterIsotropicallyInWaveFrame(
              finalSpeed.value, result.state.mu,
              branchSign * std::fabs(background.alfvenSpeedMPerS),
              speedOfLightMPerS, random);
      if (!scattered.status.ok()) {
        result.status = scattered.status;
        return result;
      }
      const ScalarResult scatteredMomentum = MomentumFromSpeed(
          scattered.speedMPerS, massKg, speedOfLightMPerS);
      if (!scatteredMomentum.status.ok()) {
        result.status = scatteredMomentum.status;
        return result;
      }
      contribution.postWaveMomentumKgMPerS = scatteredMomentum.value;
      result.state.momentumKgMPerS = scatteredMomentum.value;
      result.state.mu = scattered.mu;
      result.diagnostics.scatteringEvents++;
      if (plusBranch)
        result.diagnostics.plusBranchEvents++;
      else
        result.diagnostics.minusBranchEvents++;
      result.state.nextEventIndex++;

      // Only a completed event consumes a new variate.  Geometry, snapshot, or
      // timestep partitioning merely subtracts integrated optical depth.
      const double xi = random.UniformOpen01();
      if (!std::isfinite(xi) || xi <= 0.0 || xi >= 1.0) {
        result.status = Status::Error(StatusCode::InvalidParticleState,
            "random stream violated its open-unit-interval contract");
        return result;
      }
      result.state.remainingOpticalDepth = -std::log(xi);
    }
    else {
      result.state.remainingOpticalDepth -= trialHazard;
      if (!(result.state.remainingOpticalDepth > 0.0) ||
          !std::isfinite(result.state.remainingOpticalDepth)) {
        result.status = Status::Error(StatusCode::InvalidParticleState,
            "residual scattering optical depth became invalid");
        return result;
      }
    }
    if (waveAccumulator) waveAccumulator->Deposit(contribution);
  }
  result.status = Status::Ok();
  return result;
}

}  // namespace Transport
}  // namespace SEP
