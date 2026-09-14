#include "transport_common.h"
#include "coefficient_providers.h"

#include "util/sep_parker_core.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

void AbortMoverStatus(const SEP::Transport::Status& status) {
  exit(__LINE__, __FILE__, status.message.c_str());
}

}  // namespace

int SEP::ParticleMover_Parker(
    long int ptr, double dtTotal,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  (void)node;
  using namespace SEP::Transport;

  PICAdapter::ParticleContext context;
  Status status = PICAdapter::LoadParticle(ptr, &context);
  if (!status.ok()) AbortMoverStatus(status);
  if (!std::isfinite(dtTotal) || dtTotal < 0.0) {
    AbortMoverStatus(Status::Error(StatusCode::InvalidArgument,
                                   "Parker timestep must be finite and nonnegative"));
  }

  const double initialCoordinate = context.state.coordinate;
  const double initialSpeed = std::hypot(context.state.vParallelMPerS,
                                         context.state.vNormalMPerS);
  ScalarResult momentum = MomentumFromSpeed(
      initialSpeed, context.state.massKg, SpeedOfLight);
  if (!momentum.status.ok()) AbortMoverStatus(momentum.status);

  double elapsedS = 0.0;
  double currentSpeed = initialSpeed;
  double totalDisplacementM = 0.0;
  std::uint64_t eventIndex = 0;
  StepDiagnostics stepDiagnostics;

  while (elapsedS < dtTotal) {
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, dtTotal, &local);
    if (!status.ok()) AbortMoverStatus(status);

    PICAdapter::PICSpatialDiffusionProvider provider(context, local.identity);
    const SpatialDiffusionSample coefficient =
        provider.Evaluate(0.0, currentSpeed);
    if (!coefficient.status.ok()) AbortMoverStatus(coefficient.status);

    const double segmentLengthM = context.segment->GetLength();
    std::vector<StepLimit> limits;
    const double deterministicSpeed = std::fabs(
        local.plasmaAdvectionMPerS + coefficient.dKappaParallelDsMPerS);
    if (deterministicSpeed > 0.0) {
      limits.push_back(StepLimit(
          "parker-drift", 0.25 * segmentLengthM / deterministicSpeed));
    }
    if (coefficient.kappaParallelM2PerS > 0.0) {
      const double target = 0.25 * segmentLengthM;
      limits.push_back(StepLimit(
          "parker-diffusion",
          target * target / (2.0 * coefficient.kappaParallelM2PerS)));
    }

    const ScalarResult selected = SelectSubstep(
        dtTotal - elapsedS, limits, 1.0e-12, &stepDiagnostics);
    if (!selected.status.ok()) AbortMoverStatus(selected.status);

    // The stream key contains no rank or thread identifier.  For fixed campaign
    // seed, particle handle, mover operator, and substep index, the stochastic
    // displacement is invariant under OpenMP/MPI scheduling.
    KeyedRandomStream random(PICAdapter::CampaignRandomSeed,
                             static_cast<std::uint64_t>(ptr), 7, eventIndex);
    // Preserve the established cooling switch without introducing a second
    // equation: disabling cooling supplies zero divergence to the same exact
    // plasma-frame momentum operator.
    const ParkerIncrement increment = AdvanceParker(
        ParkerState(0.0, momentum.value),
        ParkerBackground(local.plasmaAdvectionMPerS,
                         SEP::AccountAdiabaticCoolingFlag
                             ? local.velocityDivergencePerS : 0.0),
        currentSpeed, selected.value, provider, random);
    if (!increment.status.ok()) AbortMoverStatus(increment.status);

    status = PICAdapter::AdvanceAlongFieldLine(
        &context, increment.displacementM);
    totalDisplacementM += increment.displacementM;
    elapsedS += selected.value;
    eventIndex++;
    momentum.value = increment.state.momentumKgMPerS;

    const ScalarResult updatedSpeed = SpeedFromMomentum(
        momentum.value, context.state.massKg, SpeedOfLight);
    if (!updatedSpeed.status.ok()) AbortMoverStatus(updatedSpeed.status);
    currentSpeed = updatedSpeed.value;

    if (status.code == StatusCode::OutOfDomain) {
      // A deleted particle cannot be deferred because its statistical weight is
      // no longer available at the post-step reduction boundary.  Boundary-flux
      // coupling remains disabled until the accumulator stores weight/species
      // explicitly; transport and absorbing-boundary behavior are unaffected.
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    if (!status.ok()) AbortMoverStatus(status);
  }

  const ScalarResult finalSpeed = SpeedFromMomentum(
      momentum.value, context.state.massKg, SpeedOfLight);
  if (!finalSpeed.status.ok()) AbortMoverStatus(finalSpeed.status);
  if (initialSpeed > 0.0) {
    const double scale = finalSpeed.value / initialSpeed;
    context.state.vParallelMPerS *= scale;
    context.state.vNormalMPerS *= scale;
  }

  if (SEP::Offset::MeanFreePath != -1 &&
      SEP::Sampling::MeanFreePath::active_flag && finalSpeed.value > 0.0) {
    PICAdapter::LocalBackground local;
    status = PICAdapter::EvaluateLocalBackground(context, dtTotal, &local);
    if (!status.ok()) AbortMoverStatus(status);
    PICAdapter::PICSpatialDiffusionProvider provider(context, local.identity);
    const SpatialDiffusionSample coefficient =
        provider.Evaluate(0.0, finalSpeed.value);
    if (!coefficient.status.ok()) AbortMoverStatus(coefficient.status);
    *reinterpret_cast<double*>(context.data + SEP::Offset::MeanFreePath) =
        3.0 * coefficient.kappaParallelM2PerS / finalSpeed.value;
  }

  // Parker feedback is intentionally restricted to the existing isotropic
  // streaming closure.  When that subsystem is inactive, the mover remains
  // one-way coupled and never mutates wave energy directly.
  if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
      SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode && dtTotal > 0.0) {
    PICAdapter::QueueAveragedWaveContribution(
        context.state.fieldLineId, ptr, dtTotal, finalSpeed.value,
        initialCoordinate, context.state.coordinate, totalDisplacementM,
        eventIndex);
  }

  status = PICAdapter::CommitAndAttach(context);
  if (!status.ok()) AbortMoverStatus(status);
  return _PARTICLE_MOTION_FINISHED_;
}
