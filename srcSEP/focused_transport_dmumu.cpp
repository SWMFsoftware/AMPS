#include "transport_common.h"
#include "coefficient_providers.h"

#include "util/sep_focused_transport_core.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

void AbortFocusedStatus(const SEP::Transport::Status& status) {
  exit(__LINE__, __FILE__, status.message.c_str());
}

struct SurvivingWaveContribution {
  std::string turbulenceStateIdentity;
  double dtS;
  double vParallelMPerS;
  double vNormalMPerS;
  double startCoordinate;
  double finishCoordinate;
  double signedPathM;
  std::uint64_t eventIndex;
};

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

  double speed = std::hypot(context.state.vParallelMPerS,
                            context.state.vNormalMPerS);
  ScalarResult momentum = MomentumFromSpeed(
      speed, context.state.massKg, SpeedOfLight);
  if (!momentum.status.ok()) AbortFocusedStatus(momentum.status);
  double mu = speed > 0.0 ? context.state.vParallelMPerS / speed : 0.0;
  mu = std::max(-1.0, std::min(1.0, mu));

  double elapsedS = 0.0;
  std::uint64_t eventIndex = 0;
  StepDiagnostics diagnostics;
  std::vector<SurvivingWaveContribution> survivingWaveContributions;
  while (elapsedS < dtTotal) {
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, dtTotal, &local);
    if (!status.ok()) AbortFocusedStatus(status);
    PICAdapter::PICPitchAngleDiffusionProvider provider(
        context, local.identity);
    const PitchAngleDiffusionSample coefficient = provider.Evaluate(
        0.0, momentum.value, mu);
    if (!coefficient.status.ok()) AbortFocusedStatus(coefficient.status);

    const double segmentLengthM = context.segment->GetLength();
    std::vector<StepLimit> limits;
    const double streamingSpeed =
        std::fabs(local.plasmaAdvectionMPerS + speed * mu);
    if (streamingSpeed > 0.0) {
      limits.push_back(StepLimit(
          "fte-streaming", 0.25 * segmentLengthM / streamingSpeed));
    }
    const double focusingRate = 0.5 * std::max(0.0, 1.0 - mu * mu) *
        speed * std::fabs(local.dLnAbsBdsPerM);
    if (focusingRate > 0.0)
      limits.push_back(StepLimit("fte-focusing", 0.05 / focusingRate));
    if (coefficient.dMuMuPerS > 0.0)
      limits.push_back(StepLimit(
          "fte-diffusion", 0.05 * 0.05 /
                               (2.0 * coefficient.dMuMuPerS)));
    if (std::fabs(coefficient.dDmuMuDmuPerS) > 0.0)
      limits.push_back(StepLimit(
          "fte-diffusion-drift",
          0.05 / std::fabs(coefficient.dDmuMuDmuPerS)));

    const ScalarResult selected = SelectSubstep(
        dtTotal - elapsedS, limits, 1.0e-12, &diagnostics);
    if (!selected.status.ok()) AbortFocusedStatus(selected.status);

    const double startCoordinate = context.state.coordinate;
    KeyedRandomStream random(PICAdapter::CampaignRandomSeed,
                             static_cast<std::uint64_t>(ptr), 8, eventIndex);
    ThreadLocalWaveAccumulator identityAccumulator;
    // The compatibility cooling switch acts only on the cooling operator.
    // Focusing, velocity-gradient drift, streaming, and scattering still read
    // the same immutable background when cooling is disabled.
    const FocusedTransportIncrement increment =
        AdvanceFocusedTransportDmumu(
            FocusedTransportState(0.0, momentum.value, mu),
            FocusedTransportBackground(
                local.dLnAbsBdsPerM, local.plasmaAdvectionMPerS,
                local.parallelVelocityGradientPerS,
                SEP::AccountAdiabaticCoolingFlag
                    ? local.velocityDivergencePerS : 0.0),
            context.state.massKg, SpeedOfLight, selected.value,
            provider, random,
            SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode
                ? &identityAccumulator : NULL);
    if (!increment.status.ok()) AbortFocusedStatus(increment.status);

    status = PICAdapter::AdvanceAlongFieldLine(&context,
                                               increment.displacementM);
    momentum.value = increment.state.momentumKgMPerS;
    mu = increment.state.mu;
    const ScalarResult updatedSpeed = SpeedFromMomentum(
        momentum.value, context.state.massKg, SpeedOfLight);
    if (!updatedSpeed.status.ok()) AbortFocusedStatus(updatedSpeed.status);
    speed = updatedSpeed.value;

    if (status.code == StatusCode::OutOfDomain) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    if (!status.ok()) AbortFocusedStatus(status);

    context.state.vParallelMPerS = speed * mu;
    context.state.vNormalMPerS =
        speed * std::sqrt(std::max(0.0, 1.0 - mu * mu));
    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
        SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
        !identityAccumulator.Contributions().empty()) {
      // The identity recorded by the core came from the same provider sample
      // used for this Dmumu kick.  Keep the record mover-local until the full
      // timestep survives: an absorbing boundary deletes the particle record,
      // so publishing earlier substeps would leave a dangling PIC handle in
      // the post-step wave queue.
      SurvivingWaveContribution contribution;
      contribution.turbulenceStateIdentity =
          identityAccumulator.Contributions()[0].turbulenceStateIdentity;
      contribution.dtS = selected.value;
      contribution.vParallelMPerS = context.state.vParallelMPerS;
      contribution.vNormalMPerS = context.state.vNormalMPerS;
      contribution.startCoordinate = startCoordinate;
      contribution.finishCoordinate = context.state.coordinate;
      contribution.signedPathM = increment.displacementM;
      contribution.eventIndex = eventIndex;
      survivingWaveContributions.push_back(contribution);
    }

    elapsedS += selected.value;
    eventIndex++;
  }

  status = PICAdapter::CommitAndAttach(context);
  if (!status.ok()) AbortFocusedStatus(status);

  // Commit the deferred coupling records only after the particle has been
  // reattached successfully.  FlushWaveContributions() later sorts records
  // from every worker, making accumulation order independent of scheduling.
  for (const SurvivingWaveContribution& contribution :
       survivingWaveContributions) {
    PICAdapter::QueueFocusedWaveContribution(
        contribution.turbulenceStateIdentity,
        context.state.fieldLineId, ptr, contribution.dtS,
        contribution.vParallelMPerS, contribution.vNormalMPerS,
        contribution.startCoordinate, contribution.finishCoordinate,
        contribution.signedPathM, contribution.eventIndex);
  }
  return _PARTICLE_MOTION_FINISHED_;
}
