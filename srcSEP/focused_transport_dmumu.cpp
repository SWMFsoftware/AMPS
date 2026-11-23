#include "transport_common.h"
#include "reduced_shock_boundary.h"
#include "coefficient_providers.h"

#include "util/sep_focused_transport_core.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

void AbortFocusedStatus(const SEP::Transport::Status& status) {
  exit(__LINE__, __FILE__, status.message.c_str());
}

}  // namespace

int SEP::ParticleMover_FocusedTransport_Dmumu(
    long int ptr, double dtTotal,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  (void)node;
  using namespace SEP::Transport;

  PICAdapter::ParticleContext context;
  Status status = PICAdapter::LoadParticle(ptr, &context);
  if (!status.ok()) AbortFocusedStatus(status);
  if (!std::isfinite(dtTotal) || dtTotal < 0.0) {
    AbortFocusedStatus(Status::Error(StatusCode::InvalidArgument,
        "focused-transport timestep must be finite and nonnegative"));
  }

  status = PICAdapter::PrepareParticleStep(&context, dtTotal, &dtTotal);
  if (!status.ok()) AbortFocusedStatus(status);

  double speed = std::hypot(context.state.vParallelMPerS,
                            context.state.vNormalMPerS);
  ScalarResult momentum = MomentumFromSpeed(
      speed, context.state.massKg, SpeedOfLight);
  if (!momentum.status.ok()) AbortFocusedStatus(momentum.status);
  bool initialReturn=false;
  status=ReducedParticles::Native::AbsorbInitialReturn(context,momentum.value,&initialReturn);
  if(!status.ok())AbortFocusedStatus(status);
  if(initialReturn)return _PARTICLE_LEFT_THE_DOMAIN_;
  double mu = speed > 0.0 ? context.state.vParallelMPerS / speed : 0.0;
  mu = std::max(-1.0, std::min(1.0, mu));

  double elapsedS = 0.0;
  std::uint64_t eventIndex = 0;
  KeyedRandomStream random(PICAdapter::CampaignRandomSeed,
                           context.stableParticleId, 8, 0);
  status = PICAdapter::ResumeRandomStream(
      context, PICAdapter::StochasticOperator::Dmumu, &random, &eventIndex);
  if (!status.ok()) AbortFocusedStatus(status);
  StepDiagnostics diagnostics;
  // Keep only the last checked numerical interval for this particle move.
  // New geometry is always rechecked before scattering noise; physical dt,
  // its final endpoint and the sampling exposure remain caller-owned.
  BoundaryIntervalControl boundaryControl;
  while (elapsedS < dtTotal) {
    if (eventIndex == std::numeric_limits<std::uint64_t>::max())
      AbortFocusedStatus(Status::Error(StatusCode::StepUnderflow,
          "Dmumu accepted-interval identity exhausted"));
    const NumericalTolerances& tolerances = ActiveNumericalTolerances();
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, &local);
    if (!status.ok()) AbortFocusedStatus(status);
    PICAdapter::PICPitchAngleDiffusionProvider provider(
        context, local.identity);
    const PitchAngleDiffusionSample coefficient = provider.Evaluate(
        0.0, momentum.value, mu);
    if (!coefficient.status.ok()) AbortFocusedStatus(coefficient.status);

    const double segmentLengthM = context.segment->GetLength();
    std::vector<StepLimit> limits;
    // Use the same geometry-relative speed as the actual core displacement.
    // Adding U here alone would impose an unrelated Eulerian crossing limit
    // even after removing the double advection from material particle motion.
    const double streamingSpeed =
        std::fabs(local.coordinateAdvectionMPerS + speed * mu);
    if (streamingSpeed > 0.0) {
      limits.push_back(StepLimit(
          "fte-streaming", tolerances.geometryFraction * segmentLengthM /
                               streamingSpeed));
    }
    const double focusingRate = 0.5 * std::max(0.0, 1.0 - mu * mu) *
        speed * std::fabs(local.dLnAbsBdsPerM);
    if (focusingRate > 0.0)
      limits.push_back(StepLimit(
          "fte-focusing", tolerances.focusingPitchChange / focusingRate));
    if (coefficient.dMuMuPerS > 0.0)
      limits.push_back(StepLimit(
          "fte-diffusion", tolerances.stochasticPitchRms *
                               tolerances.stochasticPitchRms /
                               (2.0 * coefficient.dMuMuPerS)));
    if (std::fabs(coefficient.dDmuMuDmuPerS) > 0.0)
      limits.push_back(StepLimit(
          "fte-diffusion-drift",
          tolerances.focusingPitchChange /
              std::fabs(coefficient.dDmuMuDmuPerS)));
    if (SEP::AccountAdiabaticCoolingFlag &&
        std::fabs(local.velocityDivergencePerS) > 0.0)
      limits.push_back(StepLimit(
          "fte-cooling", tolerances.coolingLogChange /
                             std::fabs(local.velocityDivergencePerS)));
    const ScalarResult snapshotLimit = PICAdapter::ParticleSnapshotRemaining(
        context, elapsedS, local.generation);
    if (!snapshotLimit.status.ok() || snapshotLimit.value <= 0.0)
      AbortFocusedStatus(snapshotLimit.status.ok()
          ? Status::Error(StatusCode::StepUnderflow,
                          "Dmumu mover reached the snapshot boundary")
          : snapshotLimit.status);
    if (std::isfinite(snapshotLimit.value))
      limits.push_back(StepLimit("fte-dmumu-snapshot", snapshotLimit.value));

    const ScalarResult selected = PICAdapter::SelectParticleSubstep(
        context, elapsedS, limits, tolerances.minimumStepS, &diagnostics);
    if (!selected.status.ok()) AbortFocusedStatus(selected.status);
    double intervalS=selected.value;
    ReducedParticles::Native::BoundaryInterval boundary;
    status=ReducedParticles::Native::PrepareBoundary(context,elapsedS,&intervalS,&boundary,&boundaryControl);
    if(!status.ok())AbortFocusedStatus(status);

    const double startCoordinate = context.state.coordinate;
    const double preMomentum = momentum.value;
    const double preParallel = context.state.vParallelMPerS;
    const double preNormal = context.state.vNormalMPerS;
    ThreadLocalWaveAccumulator identityAccumulator;
    FocusedTransportBackground focusedBackground(
        local.dLnAbsBdsPerM, local.coordinateAdvectionMPerS,
        local.parallelVelocityGradientPerS,
        SEP::AccountAdiabaticCoolingFlag
            ? local.velocityDivergencePerS : 0.0);
    // WP07 uses the full gyrotropic coefficients and the independently
    // published bb:grad(U).  Disabling cooling zeros both momentum-driving
    // invariants while retaining magnetic focusing and scattering.
    focusedBackground.fieldAlignedStrainPerS =
        SEP::AccountAdiabaticCoolingFlag
            ? local.fieldAlignedStrainPerS : 0.0;
    focusedBackground.equationMode = FocusedEquationMode::FullGyrotropic;
    // The compatibility cooling switch acts only on the cooling operator.
    // Focusing, velocity-gradient drift, streaming, and scattering still read
    // the same immutable background when cooling is disabled.
    const FocusedTransportIncrement increment =
        AdvanceFocusedTransportDmumu(
            FocusedTransportState(0.0, momentum.value, mu),
            focusedBackground,
            context.state.massKg, SpeedOfLight, intervalS,
            provider, random,
            SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode
                ? &identityAccumulator : NULL);
    if (!increment.status.ok()) AbortFocusedStatus(increment.status);
    RecordAcceptedStep(&diagnostics);
    // This focused numerical interval has a linear accepted streaming path;
    // resolve its first moving-front hit before publishing phase space. Pitch
    // diffusion has already supplied this interval's endpoint, so absorption
    // never retries or redraws the stochastic angular increment.
    const auto hit=ReducedParticles::Native::CheckBoundary(boundary,0,
        increment.displacementM,0,intervalS,0,random);
    if(!hit.status.ok())AbortFocusedStatus(hit.status);
    if(hit.hit) {
      const double terminalP=preMomentum*std::exp(hit.fraction*std::log(increment.state.momentumKgMPerS/preMomentum));
      ReducedParticles::Native::TransportWork(context,preMomentum,terminalP);
      ReducedParticles::Native::RecordTerminal(context,ReducedParticles::Native::Terminal::ShockReturn,terminalP);
      PIC::ParticleBuffer::DeleteParticle(ptr);return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    ReducedParticles::Native::TransportWork(context,preMomentum,increment.state.momentumKgMPerS);

    status = PICAdapter::AdvanceAlongFieldLine(&context,
                                               increment.displacementM);
    momentum.value = increment.state.momentumKgMPerS;
    mu = increment.state.mu;
    const ScalarResult updatedSpeed = SpeedFromMomentum(
        momentum.value, context.state.massKg, SpeedOfLight);
    if (!updatedSpeed.status.ok()) AbortFocusedStatus(updatedSpeed.status);
    speed = updatedSpeed.value;

    context.state.vParallelMPerS = speed * mu;
    context.state.vNormalMPerS =
        speed * std::sqrt(std::max(0.0, 1.0 - mu * mu));
    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
        SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
        !identityAccumulator.Contributions().empty()) {
      PICAdapter::CouplingRecord contribution;
      contribution.turbulenceStateIdentity =
          identityAccumulator.Contributions()[0].turbulenceStateIdentity;
      contribution.fieldLineId = context.state.fieldLineId;
      contribution.species = context.state.species;
      contribution.statisticalWeight = context.statisticalWeight;
      contribution.stableParticleId = context.stableParticleId;
      contribution.snapshotGeneration = context.snapshotGeneration;
      contribution.dtS = intervalS;
      contribution.preMomentumKgMPerS = preMomentum;
      contribution.postMomentumKgMPerS = momentum.value;
      contribution.midpointParallelVelocityMPerS =
          0.5 * (preParallel + context.state.vParallelMPerS);
      contribution.midpointNormalVelocityMPerS =
          0.5 * (preNormal + context.state.vNormalMPerS);
      contribution.startCoordinate = startCoordinate;
      contribution.finishCoordinate = context.state.coordinate;
      contribution.signedPathM = increment.displacementM;
      contribution.eventIndex = eventIndex;
      contribution.intervalIndex = eventIndex;
      contribution.pitchAngleResolved = true;
      contribution.crossedBoundary = status.code == StatusCode::OutOfDomain;
      contribution.eventType = contribution.crossedBoundary
          ? PICAdapter::CouplingEventType::BoundaryExit
          : PICAdapter::CouplingEventType::DeterministicInterval;
      if (contribution.crossedBoundary) {
        double boundaryCoordinate = 0.0, inDomainPathM = 0.0;
        const Status clip = PICAdapter::ClipPathToFieldLineBoundary(
            contribution.fieldLineId, startCoordinate,
            increment.displacementM, &boundaryCoordinate, &inDomainPathM);
        if (!clip.ok()) AbortFocusedStatus(clip);
        contribution.finishCoordinate = boundaryCoordinate;
        contribution.dtS *= std::fabs(
            inDomainPathM / increment.displacementM);
        contribution.signedPathM = inDomainPathM;
      }
      PICAdapter::QueueWaveContribution(contribution);
    }

    if (status.code == StatusCode::OutOfDomain) {
      ReducedParticles::Native::RecordTerminal(context,ReducedParticles::Native::Terminal::DomainExit,momentum.value);
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    if (!status.ok()) AbortFocusedStatus(status);

    elapsedS += intervalS;
    eventIndex++;
  }

  status = PICAdapter::StageRandomStream(
      &context, PICAdapter::StochasticOperator::Dmumu, random, eventIndex);
  if (!status.ok()) AbortFocusedStatus(status);
  status = PICAdapter::CommitAndAttach(context);
  if (!status.ok()) AbortFocusedStatus(status);

  return _PARTICLE_MOTION_FINISHED_;
}
