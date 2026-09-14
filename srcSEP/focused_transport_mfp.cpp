#include "coefficient_providers.h"
#include "transport_common.h"

#include "util/sep_focused_transport_mfp_core.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

void AbortMfpStatus(const SEP::Transport::Status& status) {
  exit(__LINE__, __FILE__, status.message.c_str());
}

struct SurvivingMfpWaveContribution {
  std::string turbulenceStateIdentity;
  double dtS = 0.0;
  double vParallelMPerS = 0.0;
  double vNormalMPerS = 0.0;
  double startCoordinate = 0.0;
  double finishCoordinate = 0.0;
  double signedPathM = 0.0;
  std::uint64_t eventIndex = 0;
};

double ShockCrossingTimeS(
    const SEP::Transport::PICAdapter::ParticleContext& context,
    double signedStreamingSpeedMPerS) {
  namespace FL = PIC::FieldLine;
  if (!FL::DatumAtVertexShockLocationDistance.is_active() ||
      signedStreamingSpeedMPerS == 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  double* distance0 = context.segment->GetBegin()->GetDatum_ptr(
      FL::DatumAtVertexShockLocationDistance);
  double* distance1 = context.segment->GetEnd()->GetDatum_ptr(
      FL::DatumAtVertexShockLocationDistance);
  if (!distance0 || !distance1 || !std::isfinite(distance0[0]) ||
      !std::isfinite(distance1[0])) {
    return std::numeric_limits<double>::infinity();
  }

  // Shock distance is stored in AU, whereas every transport operator uses SI.
  // Its sign increases with field-line arc length.  A crossing can therefore
  // occur only when the particle velocity points toward distance zero.
  const double fraction = context.state.coordinate -
                          std::floor(context.state.coordinate);
  const double signedDistanceM =
      ((1.0 - fraction) * distance0[0] + fraction * distance1[0]) * _AU_;
  if (signedDistanceM * signedStreamingSpeedMPerS >= 0.0)
    return std::numeric_limits<double>::infinity();
  return std::fabs(signedDistanceM / signedStreamingSpeedMPerS);
}

}  // namespace

int SEP::ParticleMover_FocusedTransport_EventDriven(
    long int ptr, double dtTotal,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  (void)node;
  using namespace SEP::Transport;

  PICAdapter::ParticleContext context;
  Status status = PICAdapter::LoadParticle(ptr, &context);
  if (!status.ok()) AbortMfpStatus(status);
  if (!std::isfinite(dtTotal) || dtTotal < 0.0) {
    AbortMfpStatus(Status::Error(StatusCode::InvalidArgument,
        "event-driven focused timestep must be finite and nonnegative"));
  }
  if (dtTotal == 0.0) {
    status = PICAdapter::CommitAndAttach(context);
    if (!status.ok()) AbortMfpStatus(status);
    return _PARTICLE_MOTION_FINISHED_;
  }

  double speed = std::hypot(context.state.vParallelMPerS,
                            context.state.vNormalMPerS);
  ScalarResult momentum = MomentumFromSpeed(
      speed, context.state.massKg, SpeedOfLight);
  if (!momentum.status.ok()) AbortMfpStatus(momentum.status);
  double mu = speed > 0.0 ? context.state.vParallelMPerS / speed : 0.0;
  mu = std::max(-1.0, std::min(1.0, mu));

  double elapsedS = 0.0;
  std::uint64_t eventIndex = 0;
  StepDiagnostics stepDiagnostics;
  std::vector<SurvivingMfpWaveContribution> survivingContributions;
  while (elapsedS < dtTotal) {
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, dtTotal, &local);
    if (!status.ok()) AbortMfpStatus(status);

    PICAdapter::PICMeanFreePathProvider provider(context, local.identity);
    const MeanFreePathSample coefficient = provider.Evaluate(
        0.0, momentum.value, mu);
    if (!coefficient.status.ok()) AbortMfpStatus(coefficient.status);

    const double remainingS = dtTotal - elapsedS;
    const double segmentLengthM = context.segment->GetLength();
    const double signedStreamingSpeed =
        local.plasmaAdvectionMPerS + speed * mu;
    std::vector<StepLimit> limits;
    if (signedStreamingSpeed != 0.0) {
      limits.push_back(StepLimit(
          "fte-mfp-segment",
          0.25 * segmentLengthM / std::fabs(signedStreamingSpeed)));
    }
    const double focusingRate =
        0.5 * std::max(0.0, 1.0 - mu * mu) * speed *
        std::fabs(local.dLnAbsBdsPerM);
    if (focusingRate > 0.0)
      limits.push_back(StepLimit("fte-mfp-focusing", 0.05 / focusingRate));
    if (SEP::AccountAdiabaticCoolingFlag &&
        std::fabs(local.velocityDivergencePerS) > 0.0) {
      limits.push_back(StepLimit(
          "fte-mfp-cooling",
          0.05 / std::fabs(local.velocityDivergencePerS)));
    }

    // The elapsed portion of this particle step consumes the same immutable
    // snapshot validity budget.  A non-positive residual is a hard error: using
    // new fields would require ending the global PIC read phase first.
    const double snapshotRemainingS =
        local.snapshotSecondsRemaining - elapsedS;
    if (snapshotRemainingS <= 0.0) {
      AbortMfpStatus(Status::Error(StatusCode::StepUnderflow,
          "event-driven mover reached the background validity boundary"));
    }
    limits.push_back(StepLimit("fte-mfp-snapshot", snapshotRemainingS));

    const double shockCrossingS =
        ShockCrossingTimeS(context, signedStreamingSpeed);
    if (std::isfinite(shockCrossingS) && shockCrossingS > 0.0)
      limits.push_back(StepLimit("fte-mfp-shock", shockCrossingS));

    const ScalarResult selected = SelectSubstep(
        remainingS, limits, 1.0e-12, &stepDiagnostics);
    if (!selected.status.ok()) AbortMfpStatus(selected.status);

    const double startCoordinate = context.state.coordinate;
    KeyedRandomStream random(PICAdapter::CampaignRandomSeed,
                             static_cast<std::uint64_t>(ptr), 9, eventIndex);
    ThreadLocalWaveAccumulator identityAccumulator;
    const FocusedTransportMfpIncrement increment =
        AdvanceFocusedTransportMfp(
            FocusedTransportMfpState(0.0, momentum.value, mu),
            FocusedTransportMfpBackground(
                local.dLnAbsBdsPerM, local.plasmaAdvectionMPerS,
                local.parallelVelocityGradientPerS,
                SEP::AccountAdiabaticCoolingFlag
                    ? local.velocityDivergencePerS : 0.0,
                local.alfvenSpeedMPerS),
            context.state.massKg, SpeedOfLight, selected.value,
            selected.value, provider, random,
            SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode
                ? &identityAccumulator : NULL);
    if (!increment.status.ok()) AbortMfpStatus(increment.status);

    status = PICAdapter::AdvanceAlongFieldLine(
        &context, increment.displacementM);
    momentum.value = increment.state.momentumKgMPerS;
    mu = increment.state.mu;
    const ScalarResult updatedSpeed = SpeedFromMomentum(
        momentum.value, context.state.massKg, SpeedOfLight);
    if (!updatedSpeed.status.ok()) AbortMfpStatus(updatedSpeed.status);
    speed = updatedSpeed.value;

    if (status.code == StatusCode::OutOfDomain) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    if (!status.ok()) AbortMfpStatus(status);

    context.state.vParallelMPerS = speed * mu;
    context.state.vNormalMPerS =
        speed * std::sqrt(std::max(0.0, 1.0 - mu * mu));
    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
        SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
        !identityAccumulator.Contributions().empty()) {
      // Defer publication until the particle survives and is reattached.  The
      // deterministic post-step queue then sorts records independently of the
      // OpenMP worker that happened to advance this event stream.
      SurvivingMfpWaveContribution contribution;
      contribution.turbulenceStateIdentity =
          identityAccumulator.Contributions()[0].turbulenceStateIdentity;
      contribution.dtS = selected.value;
      contribution.vParallelMPerS = context.state.vParallelMPerS;
      contribution.vNormalMPerS = context.state.vNormalMPerS;
      contribution.startCoordinate = startCoordinate;
      contribution.finishCoordinate = context.state.coordinate;
      contribution.signedPathM = increment.displacementM;
      contribution.eventIndex = eventIndex;
      survivingContributions.push_back(contribution);
    }

    elapsedS += selected.value;
    eventIndex += 1 + increment.diagnostics.scatteringEvents;
  }

  if (SEP::Offset::MeanFreePath != -1 &&
      SEP::Sampling::MeanFreePath::active_flag) {
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, dtTotal, &local);
    if (!status.ok()) AbortMfpStatus(status);
    PICAdapter::PICMeanFreePathProvider provider(context, local.identity);
    const MeanFreePathSample sample = provider.Evaluate(
        0.0, momentum.value, mu);
    if (!sample.status.ok()) AbortMfpStatus(sample.status);
    *reinterpret_cast<double*>(context.data + SEP::Offset::MeanFreePath) =
        sample.lambdaParallelM;
  }

  status = PICAdapter::CommitAndAttach(context);
  if (!status.ok()) AbortMfpStatus(status);
  for (const SurvivingMfpWaveContribution& contribution :
       survivingContributions) {
    PICAdapter::QueueFocusedWaveContribution(
        contribution.turbulenceStateIdentity,
        context.state.fieldLineId, ptr, contribution.dtS,
        contribution.vParallelMPerS, contribution.vNormalMPerS,
        contribution.startCoordinate, contribution.finishCoordinate,
        contribution.signedPathM, contribution.eventIndex);
  }
  return _PARTICLE_MOTION_FINISHED_;
}
