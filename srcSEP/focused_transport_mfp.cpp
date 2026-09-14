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
  std::uint64_t eventIndex = *reinterpret_cast<std::uint64_t*>(
      context.data + SEP::Offset::MfpEventIndex);
  std::uint64_t shellIndex = 0;
  double residualOpticalDepth = *reinterpret_cast<double*>(
      context.data + SEP::Offset::MfpOpticalDepth);
  if (!(residualOpticalDepth > 0.0) ||
      !std::isfinite(residualOpticalDepth))
    residualOpticalDepth = std::numeric_limits<double>::quiet_NaN();
  StepDiagnostics stepDiagnostics;
  // One keyed stream spans every deterministic shell in this particle step.
  // Segment/snapshot partitioning therefore cannot restart the hazard stream.
  KeyedRandomStream eventRandom(PICAdapter::CampaignRandomSeed,
                                context.stableParticleId, 9, 0);
  while (elapsedS < dtTotal) {
    const NumericalTolerances& tolerances = ActiveNumericalTolerances();
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, &local);
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
          tolerances.geometryFraction * segmentLengthM /
              std::fabs(signedStreamingSpeed)));
    }
    const double focusingRate =
        0.5 * std::max(0.0, 1.0 - mu * mu) * speed *
        std::fabs(local.dLnAbsBdsPerM);
    if (focusingRate > 0.0)
      limits.push_back(StepLimit(
          "fte-mfp-focusing",
          tolerances.focusingPitchChange / focusingRate));
    if (SEP::AccountAdiabaticCoolingFlag &&
        std::fabs(local.velocityDivergencePerS) > 0.0) {
      limits.push_back(StepLimit(
          "fte-mfp-cooling",
          tolerances.coolingLogChange /
              std::fabs(local.velocityDivergencePerS)));
    }

    // The elapsed portion of this particle step consumes the same immutable
    // snapshot validity budget.  A non-positive residual is a hard error: using
    // new fields would require ending the global PIC read phase first.
    const ScalarResult snapshotLimit = ComposeSnapshotValidityLimit(
        context.particleStepEpochS, elapsedS,
        context.particleStepEpochS + local.snapshotSecondsRemaining,
        context.snapshotGeneration, local.generation);
    if (!snapshotLimit.status.ok() || snapshotLimit.value <= 0.0) {
      AbortMfpStatus(Status::Error(StatusCode::StepUnderflow,
          "event-driven mover reached the background validity boundary"));
    }
    if (std::isfinite(snapshotLimit.value))
      limits.push_back(StepLimit("fte-mfp-snapshot", snapshotLimit.value));

    const double shockCrossingS =
        ShockCrossingTimeS(context, signedStreamingSpeed);
    if (std::isfinite(shockCrossingS) && shockCrossingS > 0.0)
      limits.push_back(StepLimit(
          "fte-mfp-shock", tolerances.shockFraction * shockCrossingS));

    const ScalarResult selected = SelectSubstep(
        remainingS, limits, tolerances.minimumStepS, &stepDiagnostics);
    if (!selected.status.ok()) AbortMfpStatus(selected.status);

    const double startCoordinate = context.state.coordinate;
    const double preMomentum = momentum.value;
    const double preParallel = context.state.vParallelMPerS;
    const double preNormal = context.state.vNormalMPerS;
    ThreadLocalWaveAccumulator identityAccumulator;
    FocusedTransportMfpState coreState(0.0, momentum.value, mu);
    coreState.remainingOpticalDepth = residualOpticalDepth;
    coreState.nextEventIndex = eventIndex;
    FocusedTransportMfpBackground focusedBackground(
        local.dLnAbsBdsPerM, local.plasmaAdvectionMPerS,
        local.parallelVelocityGradientPerS,
        SEP::AccountAdiabaticCoolingFlag
            ? local.velocityDivergencePerS : 0.0,
        local.alfvenSpeedMPerS);
    focusedBackground.fieldAlignedStrainPerS =
        SEP::AccountAdiabaticCoolingFlag
            ? local.fieldAlignedStrainPerS : 0.0;
    focusedBackground.equationMode = FocusedEquationMode::FullGyrotropic;
    const FocusedTransportMfpIncrement increment =
        AdvanceFocusedTransportMfp(
            coreState, focusedBackground,
            context.state.massKg, SpeedOfLight, selected.value,
            selected.value, provider, eventRandom,
            SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode
                ? &identityAccumulator : NULL);
    if (!increment.status.ok()) AbortMfpStatus(increment.status);
    RecordAcceptedStep(&stepDiagnostics);

    status = PICAdapter::AdvanceAlongFieldLine(
        &context, increment.displacementM);
    momentum.value = increment.state.momentumKgMPerS;
    mu = increment.state.mu;
    residualOpticalDepth = increment.state.remainingOpticalDepth;
    *reinterpret_cast<double*>(
        context.data + SEP::Offset::MfpOpticalDepth) = residualOpticalDepth;
    *reinterpret_cast<std::uint64_t*>(
        context.data + SEP::Offset::MfpEventIndex) =
            increment.state.nextEventIndex;
    const ScalarResult updatedSpeed = SpeedFromMomentum(
        momentum.value, context.state.massKg, SpeedOfLight);
    if (!updatedSpeed.status.ok()) AbortMfpStatus(updatedSpeed.status);
    speed = updatedSpeed.value;

    context.state.vParallelMPerS = speed * mu;
    context.state.vNormalMPerS =
        speed * std::sqrt(std::max(0.0, 1.0 - mu * mu));
    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
        SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
        !identityAccumulator.Contributions().empty()) {
      double intervalStart = startCoordinate;
      PIC::FieldLine::cFieldLineSegment* intervalSegment =
          PIC::FieldLine::FieldLinesAll[context.state.fieldLineId].GetSegment(
              intervalStart);
      for (std::size_t i = 0;
           i < identityAccumulator.Contributions().size(); ++i) {
        const WaveContribution& emitted =
            identityAccumulator.Contributions()[i];
        double intervalFinish =
            PIC::FieldLine::FieldLinesAll[context.state.fieldLineId].move(
                intervalStart, emitted.displacementM, intervalSegment);
        intervalSegment =
            PIC::FieldLine::FieldLinesAll[context.state.fieldLineId].GetSegment(
                intervalFinish);

        PICAdapter::CouplingRecord contribution;
        contribution.turbulenceStateIdentity =
            emitted.turbulenceStateIdentity;
        contribution.fieldLineId = context.state.fieldLineId;
        contribution.species = context.state.species;
        contribution.statisticalWeight = context.statisticalWeight;
        contribution.stableParticleId = context.stableParticleId;
        contribution.snapshotGeneration = context.snapshotGeneration;
        contribution.dtS = emitted.intervalS;
        contribution.preMomentumKgMPerS = preMomentum;
        contribution.postMomentumKgMPerS = momentum.value;
        contribution.midpointParallelVelocityMPerS =
            0.5 * (preParallel + context.state.vParallelMPerS);
        contribution.midpointNormalVelocityMPerS =
            0.5 * (preNormal + context.state.vNormalMPerS);
        contribution.startCoordinate = intervalStart;
        contribution.finishCoordinate = intervalFinish;
        contribution.signedPathM = emitted.displacementM;
        contribution.eventIndex = emitted.eventIndex;
        contribution.intervalIndex =
            (shellIndex << 32) | static_cast<std::uint64_t>(i);
        contribution.pitchAngleResolved = true;
        contribution.crossedBoundary = intervalSegment == NULL;
        contribution.eventType = contribution.crossedBoundary
            ? PICAdapter::CouplingEventType::BoundaryExit
            : (emitted.scatteringEventAtEnd
                ? PICAdapter::CouplingEventType::ScatteringEvent
                : PICAdapter::CouplingEventType::DeterministicInterval);
        if (contribution.crossedBoundary) {
          double boundaryCoordinate = 0.0, inDomainPathM = 0.0;
          const Status clip = PICAdapter::ClipPathToFieldLineBoundary(
              contribution.fieldLineId, intervalStart,
              emitted.displacementM, &boundaryCoordinate, &inDomainPathM);
          if (!clip.ok()) AbortMfpStatus(clip);
          contribution.finishCoordinate = boundaryCoordinate;
          contribution.dtS *= std::fabs(
              inDomainPathM / emitted.displacementM);
          contribution.signedPathM = inDomainPathM;
        }
        PICAdapter::QueueWaveContribution(contribution);
        intervalStart = intervalFinish;
        if (!intervalSegment) break;
      }
    }

    if (status.code == StatusCode::OutOfDomain) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    if (!status.ok()) AbortMfpStatus(status);

    elapsedS += selected.value;
    eventIndex = increment.state.nextEventIndex;
    ++shellIndex;
  }

  if (SEP::Offset::MeanFreePath != -1 &&
      SEP::Sampling::MeanFreePath::active_flag) {
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, &local);
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
  return _PARTICLE_MOTION_FINISHED_;
}
