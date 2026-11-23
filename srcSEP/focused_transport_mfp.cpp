#include "coefficient_providers.h"
#include "transport_common.h"
#include "reduced_shock_boundary.h"

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

// Distance to the shock, as a fraction of the particle's segment, below which a
// particle is on the shock (ShockCrossingTimeS).  The coordinate resolves about
// 2e-16 x (segment index) of a segment; 1e-9 is far above that for any line and
// far below any physical scale (0.1 m for a 1e8 m segment).
constexpr double kOnShockSegmentFraction = 1.0e-9;

double ShockCrossingTimeS(
    const SEP::Transport::PICAdapter::ParticleContext& context,
    double signedStreamingSpeedMPerS) {
  namespace FL = PIC::FieldLine;
  double providerTime=0;bool supplied=false;
  const auto providerStatus=SEP::Transport::PICAdapter::QueryShockCrossingTime(
      context,signedStreamingSpeedMPerS,&providerTime,&supplied);
  if(!providerStatus.ok())AbortMfpStatus(providerStatus);
  if(supplied)return providerTime;
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
  // A particle closer to the shock than kOnShockSegmentFraction of its segment
  // is on the shock: the crossing is not resolved further (no shock limit).
  // The field-line coordinate (segment index + fraction) cannot resolve such a
  // distance: a substep limited to the crossing time moves the particle less
  // than the coordinate's resolution, so the distance never shrinks.  On
  // 2026-11-18 a particle born on the shock (1.6e-8 m from it) took 2e6
  // substeps of 2.2e-12 s (just above minimumStepS, so the guard on the
  // crossing time below did not apply) and covered 4 us of its 0.88 s step:
  // SEP--Parker_spiral--field_line never finished iteration 400.
  if (std::fabs(signedDistanceM) <=
      kOnShockSegmentFraction * context.segment->GetLength())
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
  status = PICAdapter::PrepareParticleStep(&context, dtTotal, &dtTotal);
  if (!status.ok()) AbortMfpStatus(status);

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
  bool initialReturn=false;
  status=ReducedParticles::Native::AbsorbInitialReturn(context,momentum.value,&initialReturn);
  if(!status.ok())AbortMfpStatus(status);
  if(initialReturn)return _PARTICLE_LEFT_THE_DOMAIN_;
  double mu = speed > 0.0 ? context.state.vParallelMPerS / speed : 0.0;
  mu = std::max(-1.0, std::min(1.0, mu));

  double elapsedS = 0.0;
  std::uint64_t eventIndex = context.mfpEventIndex;
  std::uint64_t shellIndex = 0;
  double residualOpticalDepth = context.mfpOpticalDepth;
  StepDiagnostics stepDiagnostics;
  // This local numerical hint avoids repeating rejected coarse front checks.
  // It survives internal event intervals only; all new intervals are graded
  // before event/noise draws, and complete PIC dt and persistent RNG history
  // retain their original ownership and meaning.
  BoundaryIntervalControl boundaryControl;
  // Carry the actual consumed stream across outer calls as well as deterministic
  // shells. Residual optical depth alone cannot restore event branch/pitch draws.
  KeyedRandomStream eventRandom(PICAdapter::CampaignRandomSeed,
                                context.stableParticleId, 9, 0);
  status = PICAdapter::ResumeRandomStream(
      context, PICAdapter::StochasticOperator::MeanFreePath,
      &eventRandom, &shellIndex);
  if (!status.ok()) AbortMfpStatus(status);
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
    // Geometry limits and both core predictor/corrector paths use velocity
    // relative to the line. Physical U/grad(U) remain available separately;
    // material coordinates contribute zero bulk displacement, not zero plasma.
    const double signedStreamingSpeed =
        local.coordinateAdvectionMPerS + speed * mu;
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
    const ScalarResult snapshotLimit = PICAdapter::ParticleSnapshotRemaining(
        context, elapsedS, local.generation);
    if (!snapshotLimit.status.ok() || snapshotLimit.value <= 0.0) {
      AbortMfpStatus(Status::Error(StatusCode::StepUnderflow,
          "event-driven mover reached the background validity boundary"));
    }
    if (std::isfinite(snapshotLimit.value))
      limits.push_back(StepLimit("fte-mfp-snapshot", snapshotLimit.value));

    const double shockCrossingS =
        ShockCrossingTimeS(context, signedStreamingSpeed);
    // A particle that is already at the shock within roundoff (it was injected
    // there, or a previous substep limited by the shock stopped it a roundoff
    // distance short of it) has a crossing time far below the minimum substep;
    // the limit then stopped the run ("stability limit is below the declared
    // minimum substep"), and flooring it at the minimum substep let a slow
    // particle approach the shock without ever crossing it (its displacement
    // is below the resolution of the field-line coordinate).  Such a particle
    // is on the shock: the crossing is not resolved further.
    if (!ReducedParticles::Native::Enabled()&&std::isfinite(shockCrossingS) && shockCrossingS > 0.0 &&
        tolerances.shockFraction * shockCrossingS >= tolerances.minimumStepS)
      limits.push_back(StepLimit(
          "fte-mfp-shock", tolerances.shockFraction * shockCrossingS));

    const ScalarResult selected = PICAdapter::SelectParticleSubstep(
        context, elapsedS, limits, tolerances.minimumStepS, &stepDiagnostics);
    if (!selected.status.ok()) AbortMfpStatus(selected.status);
    double intervalS=selected.value;
    ReducedParticles::Native::BoundaryInterval boundary;
    status=ReducedParticles::Native::PrepareBoundary(context,elapsedS,&intervalS,&boundary,&boundaryControl);
    if(!status.ok())AbortMfpStatus(status);
    MfpIntervalBoundary firstPassage;
    if(boundary.active)firstPassage=[&](const FocusedTransportMfpState& begin,
        const FocusedTransportMfpState& end,double elapsedBegin,double elapsedEnd) {
      // The core supplies both closed endpoints. Adding its rounded duration
      // back to the beginning can overshoot this prepared interval by an ULP
      // after a scatter; use the authoritative end without a timestamp clamp.
      const auto hit=ReducedParticles::Native::CheckBoundary(boundary,
          begin.arcLengthM,end.arcLengthM,elapsedBegin,elapsedEnd,0,eventRandom);
      MfpBoundaryDecision decision;decision.status=hit.status;
      decision.absorb=hit.hit;decision.fraction=hit.fraction;return decision;
    };

    const double startCoordinate = context.state.coordinate;
    const double preMomentum = momentum.value;
    const double preParallel = context.state.vParallelMPerS;
    const double preNormal = context.state.vNormalMPerS;
    ThreadLocalWaveAccumulator identityAccumulator;
    FocusedTransportMfpState coreState(0.0, momentum.value, mu);
    coreState.remainingOpticalDepth = residualOpticalDepth;
    coreState.nextEventIndex = eventIndex;
    FocusedTransportMfpBackground focusedBackground(
        local.dLnAbsBdsPerM, local.coordinateAdvectionMPerS,
        local.parallelVelocityGradientPerS,
        SEP::AccountAdiabaticCoolingFlag
            ? local.velocityDivergencePerS : 0.0,
        local.alfvenSpeedMPerS);
    focusedBackground.magneticPolarityAlongLine=local.magneticPolarityAlongLine;
    focusedBackground.fieldAlignedStrainPerS =
        SEP::AccountAdiabaticCoolingFlag
            ? local.fieldAlignedStrainPerS : 0.0;
    focusedBackground.equationMode = FocusedEquationMode::FullGyrotropic;
    const FocusedTransportMfpIncrement increment =
        AdvanceFocusedTransportMfp(
            coreState, focusedBackground,
            context.state.massKg, SpeedOfLight, intervalS,
            intervalS, provider, eventRandom,
            (SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode||ReducedParticles::Native::Enabled())
                ? &identityAccumulator : NULL,firstPassage);
    if (!increment.status.ok()) AbortMfpStatus(increment.status);
    // Prescribed Alfven scatterers can exchange energy with particles even
    // when wave feedback is off. Keep that external work separate from the
    // deterministic flow/strain work; neither is an injected source energy.
    double waveWorkJ=0;
    const auto kinetic=[&](double p) {
      const double pc=p*SpeedOfLight,rest=context.state.massKg*SpeedOfLight*SpeedOfLight;
      return pc*(pc/(std::hypot(pc,rest)+rest));
    };
    for(const auto& event:identityAccumulator.Contributions())if(event.scatteringEventAtEnd)
      waveWorkJ+=kinetic(event.postWaveMomentumKgMPerS)-kinetic(event.preWaveMomentumKgMPerS);
    ReducedParticles::Native::TransportWork(context,preMomentum,increment.state.momentumKgMPerS,waveWorkJ);
    if(increment.absorbed) {
      ReducedParticles::Native::RecordTerminal(context,ReducedParticles::Native::Terminal::ShockReturn,increment.state.momentumKgMPerS);
      PIC::ParticleBuffer::DeleteParticle(ptr);return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    const std::uint64_t intervals = static_cast<std::uint64_t>(
        increment.diagnostics.deterministicIntervals);
    if (intervals > std::numeric_limits<std::uint64_t>::max() - shellIndex)
      AbortMfpStatus(Status::Error(StatusCode::StepUnderflow,
          "MFP accepted-interval identity exhausted"));
    const std::uint64_t nextShellIndex = shellIndex + intervals;
    RecordAcceptedStep(&stepDiagnostics);

    status = PICAdapter::AdvanceAlongFieldLine(
        &context, increment.displacementM);
    momentum.value = increment.state.momentumKgMPerS;
    mu = increment.state.mu;
    residualOpticalDepth = increment.state.remainingOpticalDepth;
    context.mfpOpticalDepth = residualOpticalDepth;
    context.mfpEventIndex = increment.state.nextEventIndex;
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
        // A discrete event must carry the momentum pair immediately before
        // and after its own wave-frame Lorentz scatter.  Shell-wide endpoints
        // include focusing/cooling and possibly several events, so using them
        // here would deposit the wrong particle-energy change into a branch.
        // Deterministic intervals retain the legacy shell endpoints because
        // they contribute streaming diagnostics but no branch-resolved jump.
        contribution.preMomentumKgMPerS = emitted.scatteringEventAtEnd
            ? emitted.preWaveMomentumKgMPerS : preMomentum;
        contribution.postMomentumKgMPerS = emitted.scatteringEventAtEnd
            ? emitted.postWaveMomentumKgMPerS : momentum.value;
        contribution.midpointParallelVelocityMPerS =
            0.5 * (preParallel + context.state.vParallelMPerS);
        contribution.midpointNormalVelocityMPerS =
            0.5 * (preNormal + context.state.vNormalMPerS);
        contribution.startCoordinate = intervalStart;
        contribution.finishCoordinate = intervalFinish;
        contribution.signedPathM = emitted.displacementM;
        contribution.eventIndex = emitted.eventIndex;
        contribution.resonantBranch = emitted.resonantBranch;
        contribution.intervalIndex = shellIndex + static_cast<std::uint64_t>(i);
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
      ReducedParticles::Native::RecordTerminal(context,ReducedParticles::Native::Terminal::DomainExit,momentum.value);
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    if (!status.ok()) AbortMfpStatus(status);

    elapsedS += intervalS;
    eventIndex = increment.state.nextEventIndex;
    shellIndex = nextShellIndex;
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

  status = PICAdapter::StageRandomStream(
      &context, PICAdapter::StochasticOperator::MeanFreePath,
      eventRandom, shellIndex);
  if (!status.ok()) AbortMfpStatus(status);
  status = PICAdapter::CommitAndAttach(context);
  if (!status.ok()) AbortMfpStatus(status);
  return _PARTICLE_MOTION_FINISHED_;
}
