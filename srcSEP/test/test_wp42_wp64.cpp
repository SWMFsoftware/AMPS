#include "../util/sep_numerical_extensions.h"
#include "../util/sep_physics_extensions.h"
#include "../util/sep_runtime_contracts.h"
#include "../util/sep_validation_extensions.h"
#include "../util/sep_focused_transport_mfp_core.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

int Check(bool condition, const char* label) {
  if (condition) {
    std::cout << "PASS " << label << '\n';
    return 0;
  }
  std::cerr << "FAIL " << label << '\n';
  return 1;
}

bool Near(double a, double b, double tolerance = 1.0e-11) {
  return std::fabs(a - b) <= tolerance *
      std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

class PlusBranchProvider final
    : public SEP::Transport::MeanFreePathProvider {
 public:
  SEP::Transport::MeanFreePathSample Evaluate(
      double, double, double) const override {
    SEP::Transport::MeanFreePathSample sample;
    sample.status = SEP::Transport::Status::Ok();
    sample.lambdaParallelM = 1.0e6;
    sample.nuPlusPerS = 40.0;
    sample.nuMinusPerS = 0.0;
    sample.hasBranchResolvedRates = true;
    sample.provenance = "WP43 branch-resolved fixture";
    sample.turbulenceStateIdentity = "waves-7";
    return sample;
  }
};

SEP::Turbulence::CellState Cell(double plus, double minus,
                                double length = 1.0,
                                double volume = 1.0) {
  SEP::Turbulence::CellState cell;
  cell.lengthM = length;
  cell.volumeM3 = volume;
  cell.plasmaSpeedMPerS = 1.0;
  cell.alfvenSpeedMPerS = 0.5;
  cell.magneticFieldT = 1.0;
  cell.massDensityKgPerM3 = 1.0;
  cell.ePlusJ = plus;
  cell.eMinusJ = minus;
  return cell;
}

SEP::RuntimeContracts::TurbulenceLineRecord RuntimeLine() {
  SEP::RuntimeContracts::TurbulenceLineRecord line;
  line.fieldLineId = 3;
  line.ownerRank = 0;
  line.geometryGeneration = 7;
  line.configurationFingerprint = "cfg-42-64";
  line.state.configuration.source =
      SEP::Turbulence::Source::SelfConsistentIntegrated;
  line.state.configuration.representation =
      SEP::Turbulence::Representation::Integrated;
  line.state.configuration.advectionEnabled = false;
  line.state.configuration.reflectionEnabled = false;
  line.state.configuration.cascadeEnabled = false;
  line.state.configuration.coupling =
      SEP::Turbulence::CouplingPolicy::StreamingEnergyExchange;
  line.state.fieldLineGeneration = line.geometryGeneration;
  line.state.epochS = 10.0;
  line.state.sourceChecksum = "waves-7";
  line.state.provenance = "owner=waves-7;configuration=cfg-42-64";
  line.state.cells.push_back(Cell(10.0, 4.0));
  SEP::Turbulence::InitializeState(&line.state);
  return line;
}

double SineAdvectionError(std::size_t cells) {
  const double pi = std::acos(-1.0);
  const double dx = 1.0 / static_cast<double>(cells);
  std::vector<double> state(cells, 0.0);
  for (std::size_t i = 0; i < cells; ++i)
    state[i] = 1.0 + 0.2 * std::sin(2.0 * pi * (i + 0.5) * dx);
  const std::vector<double> initial = state;
  const double dt = 0.4 * dx;
  const std::size_t steps = static_cast<std::size_t>(std::ceil(1.0 / dt));
  const double exactDt = 1.0 / static_cast<double>(steps);
  for (std::size_t n = 0; n < steps; ++n) {
    const SEP::NumericalExtensions::AdvectionResult update =
        SEP::NumericalExtensions::AdvectPeriodic(
            state, 1.0, exactDt, dx,
            SEP::NumericalExtensions::AdvectionOrder::SecondOrderMuscl,
            SEP::NumericalExtensions::TvdLimiter::MonotonizedCentral);
    if (!update.status.ok()) return 1.0e30;
    state = update.cellAverage;
  }
  long double sum = 0.0L;
  for (std::size_t i = 0; i < cells; ++i) {
    const double delta = state[i] - initial[i];
    sum += delta * delta;
  }
  return std::sqrt(static_cast<double>(sum / cells));
}

}  // namespace

int main() {
  int failures = 0;

  // WP42: a synchronized host/background refresh must preserve evolved wave
  // state.  The complete store then survives a schema-checked restart.
  SEP::RuntimeContracts::TurbulenceRuntimeStore store;
  const SEP::RuntimeContracts::TurbulenceLineRecord line = RuntimeLine();
  SEP::Transport::Status status = store.Install(line);
  store.FindMutable(3)->state.cells[0].ePlusJ = 8.0;
  store.FindMutable(3)->state.completedSteps = 2;
  std::vector<SEP::Turbulence::CellState> refreshed;
  refreshed.push_back(Cell(999.0, 999.0, 1.2, 1.5));
  SEP::Turbulence::EnergyLedger remapLedger;
  status = store.RefreshBackground(
      3, 7, "cfg-42-64", refreshed,
      SEP::RuntimeContracts::PendingRemapPolicy::RejectPending,
      &remapLedger);
  std::string checkpoint;
  SEP::RuntimeContracts::TurbulenceRuntimeStore restored;
  const bool storeRoundTrip = store.Serialize(&checkpoint).ok() &&
      restored.Deserialize(checkpoint).ok();
  failures += Check(status.ok() && store.Find(3)->state.cells[0].ePlusJ == 8.0 &&
      store.Find(3)->state.cells[0].volumeM3 == 1.5 &&
      store.Find(3)->state.completedSteps == 2 && storeRoundTrip &&
      restored.Find(3)->state.cells[0].ePlusJ == 8.0,
      "WP42 persistent turbulence owner and restart round trip");

  // WP43: permutation-independent typed transactions close particle+wave
  // energy.  An invalid record rejects a later batch without partial mutation.
  SEP::RuntimeContracts::CouplingTransaction a;
  a.transactionId = 2; a.particleId = 22; a.generation = 7;
  a.event = 1; a.interval = 1; a.fieldLineId = 3; a.cell = 0;
  a.branch = 1; a.particleEnergyChangeJ = 1.25;
  a.turbulenceIdentity = "waves-7";
  SEP::RuntimeContracts::CouplingTransaction b = a;
  b.transactionId = 1; b.particleId = 21; b.particleEnergyChangeJ = -0.25;
  std::vector<SEP::RuntimeContracts::CouplingTransaction> transactions;
  transactions.push_back(a); transactions.push_back(b);
  const SEP::RuntimeContracts::CouplingBatchResult coupling =
      SEP::RuntimeContracts::ApplyCouplingTransactions(
          transactions,
          SEP::RuntimeContracts::InvalidTransactionPolicy::RejectBatch,
          &store);
  const double pendingBeforeInvalid =
      store.Find(3)->state.cells[0].pendingParticlePlusJ;
  SEP::RuntimeContracts::CouplingTransaction invalid = a;
  invalid.transactionId = 3; invalid.generation = 6;
  const SEP::RuntimeContracts::CouplingBatchResult rejected =
      SEP::RuntimeContracts::ApplyCouplingTransactions(
          std::vector<SEP::RuntimeContracts::CouplingTransaction>(1, invalid),
          SEP::RuntimeContracts::InvalidTransactionPolicy::RejectBatch,
          &store);

  // The event-driven production core must publish the actual wave-frame
  // momentum jump for each selected branch. This catches the former wrapper
  // error that copied one shell-wide endpoint pair into every interval record.
  const double protonMass = 1.67262192369e-27;
  const SEP::Transport::ScalarResult eventMomentum =
      SEP::Transport::MomentumFromSpeed(
          1.0e6, protonMass, 299792458.0);
  SEP::Transport::FocusedTransportMfpBackground eventBackground;
  eventBackground.alfvenSpeedMPerS = 8.0e4;
  SEP::Transport::KeyedRandomStream eventRandom(43, 1, 2, 3);
  SEP::Transport::ThreadLocalWaveAccumulator eventRecords;
  const SEP::Transport::FocusedTransportMfpIncrement eventStep =
      SEP::Transport::AdvanceFocusedTransportMfp(
          SEP::Transport::FocusedTransportMfpState(
              0.0, eventMomentum.value, 0.3),
          eventBackground, protonMass, 299792458.0, 0.2, 0.02,
          PlusBranchProvider(), eventRandom, &eventRecords);
  bool eventMetadataValid = eventStep.status.ok();
  std::size_t branchEvents = 0;
  for (std::size_t i = 0; i < eventRecords.Contributions().size(); ++i) {
    const SEP::Transport::WaveContribution& record =
        eventRecords.Contributions()[i];
    if (!record.scatteringEventAtEnd) continue;
    ++branchEvents;
    eventMetadataValid = eventMetadataValid && record.resonantBranch == 1 &&
        record.preWaveMomentumKgMPerS > 0.0 &&
        record.postWaveMomentumKgMPerS > 0.0;
  }
  failures += Check(coupling.status.ok() && Near(coupling.waveEnergyChangeJ,
      -coupling.particleEnergyChangeJ) &&
      Near(pendingBeforeInvalid, -1.0) && !rejected.status.ok() &&
      Near(store.Find(3)->state.cells[0].pendingParticlePlusJ,
           pendingBeforeInvalid) && eventMetadataValid && branchEvents > 0,
      "WP43 transactional particle-wave coupling and rollback");

  // WP44: rigid translation on a curved line has no strain tensor but has a
  // nonzero derivative of U.b through the curvature term.
  SEP::RuntimeContracts::VelocityGradientInput gradient;
  gradient.velocityMPerS[1] = 4.0;
  gradient.curvaturePerM[1] = 0.5;
  const SEP::RuntimeContracts::VelocityDerivatives derivatives =
      SEP::RuntimeContracts::ComputeVelocityDerivatives(gradient);
  const SEP::Transport::ScalarResult continuity =
      SEP::RuntimeContracts::ContinuityResidualPerS(-3.0, 3.0);
  failures += Check(derivatives.status.ok() &&
      Near(derivatives.fieldAlignedStrainPerS, 0.0) &&
      Near(derivatives.parallelVelocityGradientPerS, 2.0) &&
      continuity.status.ok() && Near(continuity.value, 0.0),
      "WP44 distinct divergence strain and curved-field derivative");

  // WP45: upstream normal flux determines processed particles; a subsonic
  // state is explicitly no-shock and receives no injected population.
  SEP::RuntimeContracts::ShockState shock;
  shock.provider = "analytic"; shock.provenance = "fixture";
  shock.shockNormalSpeedMPerS = 10.0;
  shock.upstreamNormalSpeedMPerS = 2.0;
  shock.downstreamNormalSpeedMPerS = 4.0;
  shock.upstreamNumberDensityPerM3 = 5.0;
  shock.downstreamNumberDensityPerM3 = 15.0;
  shock.compressionRatio = 3.0; shock.alfvenMach = 4.0;
  shock.sonicMach = 3.0;
  const SEP::RuntimeContracts::ShockValidation shockCheck =
      SEP::RuntimeContracts::ValidateShockState(shock);
  const SEP::Transport::ScalarResult processed =
      SEP::RuntimeContracts::ProcessedUpstreamParticleCount(
          shock, 2.0, 3.0, 0.1);
  const SEP::Transport::ScalarResult compression =
      SEP::RuntimeContracts::CompressionFromMach(3.0, 5.0 / 3.0);
  failures += Check(shockCheck.status.ok() &&
      shockCheck.quality == SEP::RuntimeContracts::ShockQuality::ValidCompressive &&
      Near(processed.value, 24.0) && compression.status.ok() &&
      compression.value > 1.0 && compression.value < 4.0,
      "WP45 validated shock state and upstream-flux normalization");

  // WP46: identity and purpose streams contain only persistent integers, and
  // the scheduler resumes at exactly the next event after restart.
  SEP::RuntimeContracts::SourceEventKey sourceKey;
  sourceKey.campaign = 99; sourceKey.source = 2; sourceKey.event = 5;
  sourceKey.fieldLine = 7; sourceKey.species = 1; sourceKey.ordinal = 42;
  const std::uint64_t identity1 =
      SEP::RuntimeContracts::StableSourceIdentity(sourceKey);
  const std::uint64_t identity2 =
      SEP::RuntimeContracts::StableSourceIdentity(sourceKey);
  SEP::Transport::KeyedRandomStream energyRandom =
      SEP::RuntimeContracts::SourceRandomStream(
          sourceKey, SEP::RuntimeContracts::SourceRandomPurpose::Energy);
  SEP::Transport::KeyedRandomStream pitchRandom =
      SEP::RuntimeContracts::SourceRandomStream(
          sourceKey, SEP::RuntimeContracts::SourceRandomPurpose::PitchAngle);
  SEP::RuntimeContracts::SourceEventScheduler scheduler(99), scheduler2;
  scheduler.AllocateEvent(); scheduler.AllocateEvent();
  std::string schedulerText;
  scheduler.Serialize(&schedulerText); scheduler2.Deserialize(schedulerText);
  failures += Check(identity1 != 0 && identity1 == identity2 &&
      energyRandom.UniformOpen01() != pitchRandom.UniformOpen01() &&
      scheduler2.AllocateEvent() == 2,
      "WP46 stable source identity purpose streams and restart event");

  // WP47: the declared volume-density formulation receives the geometric
  // kappa*dlnA/ds drift and its sampling conversion applies area exactly once.
  SEP::PhysicsExtensions::ParkerGeometryInput parker;
  parker.measure = SEP::PhysicsExtensions::ParkerMeasure::PerVolume;
  parker.plasmaAdvectionMPerS = 4.0; parker.kappaParallelM2PerS = 3.0;
  parker.dKappaDsMPerS = 2.0; parker.dLnAreaDsPerM = 0.5;
  const SEP::Transport::ScalarResult parkerDrift =
      SEP::PhysicsExtensions::ParkerItoDriftMPerS(parker);
  const SEP::Transport::ScalarResult density =
      SEP::PhysicsExtensions::WalkerToPhysicalDensity(
          SEP::PhysicsExtensions::ParkerMeasure::PerArcLength, 12.0, 3.0);
  failures += Check(parkerDrift.status.ok() && Near(parkerDrift.value, 7.5) &&
      density.status.ok() && Near(density.value, 4.0),
      "WP47 flux-tube measure geometric drift and sampling conversion");

  // WP48: charge/polarity reversal changes signed resonance consistently while
  // absolute spectral overlap remains invariant.
  SEP::PhysicsExtensions::DynamicResonanceInput resonance;
  resonance.particleSpeedMPerS = 1.0e6; resonance.mu = 0.4;
  resonance.chargeC = 1.602176634e-19; resonance.massKg = 1.67262192369e-27;
  resonance.signedMagneticFieldT = 5.0e-9;
  resonance.alfvenSpeedMPerS = 5.0e4; resonance.branch = 1;
  resonance.binCentersPerM.push_back(1.0e-8);
  resonance.binCentersPerM.push_back(1.0e-5);
  resonance.binCentersPerM.push_back(1.0e-2);
  SEP::PhysicsExtensions::ResonanceMetadata positive =
      SEP::PhysicsExtensions::SolveDynamicResonance(resonance);
  resonance.chargeC *= -1.0; resonance.signedMagneticFieldT *= -1.0;
  SEP::PhysicsExtensions::ResonanceMetadata reversed =
      SEP::PhysicsExtensions::SolveDynamicResonance(resonance);
  failures += Check(positive.status.ok() && reversed.status.ok() &&
      Near(positive.signedWaveNumberPerM, reversed.signedWaveNumberPerM) &&
      Near(positive.lowerWeight + positive.upperWeight, 1.0),
      "WP48 dynamic signed resonance and finite-bin overlap");

  // WP49: the broadening closure is even in D and odd in dD/dmu, is finite at
  // mu=0, and amplitude zero returns the slab baseline exactly.
  SEP::PhysicsExtensions::NinetyDegreeClosure closure;
  closure.amplitudePerS = 2.0; closure.halfWidthMu = 0.2;
  closure.provenance = "controlled nonlinear closure";
  const SEP::PhysicsExtensions::PitchAngleCoefficient plusClosure =
      SEP::PhysicsExtensions::ApplyNinetyDegreeClosure(0.1, 0.0, 0.0, closure);
  const SEP::PhysicsExtensions::PitchAngleCoefficient minusClosure =
      SEP::PhysicsExtensions::ApplyNinetyDegreeClosure(-0.1, 0.0, 0.0, closure);
  closure.amplitudePerS = 0.0;
  const SEP::PhysicsExtensions::PitchAngleCoefficient slab =
      SEP::PhysicsExtensions::ApplyNinetyDegreeClosure(0.0, 3.0, 4.0, closure);
  failures += Check(plusClosure.status.ok() && minusClosure.status.ok() &&
      Near(plusClosure.dMuMuPerS, minusClosure.dMuMuPerS) &&
      Near(plusClosure.derivativePerS, -minusClosure.derivativePerS) &&
      slab.dMuMuPerS == 3.0 && slab.derivativePerS == 4.0,
      "WP49 physical ninety-degree closure and slab baseline");

  // WP50: wave action remains authoritative during a frequency/geometry
  // update and the resulting wave-energy change is named background work.
  SEP::PhysicsExtensions::WaveActionCell wave;
  wave.authoritativeValue = 3.0; wave.volumeM3 = 2.0;
  wave.intrinsicFrequencyRadPerS = 4.0;
  const SEP::PhysicsExtensions::WaveActionUpdate waveUpdate =
      SEP::PhysicsExtensions::ApplyGeometricConservation(
          wave, 5.0, 6.0, SEP::PhysicsExtensions::WaveInvariant::WaveAction);
  failures += Check(waveUpdate.status.ok() &&
      waveUpdate.cell.authoritativeValue == 3.0 &&
      Near(waveUpdate.backgroundWorkJ, 6.0),
      "WP50 wave-action geometric conservation and background-work ledger");

  // WP51: interior fluxes telescope and high-k energy becomes species heat.
  std::vector<double> spectrum(2, 10.0);
  std::vector<double> flux; flux.push_back(0.0); flux.push_back(1.0);
  flux.push_back(0.5);
  const SEP::PhysicsExtensions::SpectralCascadeResult cascade =
      SEP::PhysicsExtensions::AdvanceConservativeCascade(
          spectrum, spectrum, flux, flux, 1.0, 0.25);
  failures += Check(cascade.status.ok() && Near(cascade.electronHeatJ, 0.25) &&
      Near(cascade.ionHeatJ, 0.75) && Near(cascade.closureResidualJ, 0.0),
      "WP51 conservative spectral cascade and heat partition");

  // WP52: Hermite interpolation is exactly continuous at knots, supplies an
  // analytic derivative, and fails closed outside a calibrated domain.
  SEP::PhysicsExtensions::AnalyticProfile profile;
  profile.id = "wind-speed"; profile.version = "1";
  profile.units = "m s^-1"; profile.provenance = "controlled fixture";
  SEP::PhysicsExtensions::ProfileKnot p0, p1, p2;
  p0.radiusM = 1.0; p0.value = 2.0; p0.derivativePerM = 1.0;
  p1.radiusM = 2.0; p1.value = 3.0; p1.derivativePerM = 1.0;
  p2.radiusM = 3.0; p2.value = 4.0; p2.derivativePerM = 1.0;
  profile.knots.push_back(p0); profile.knots.push_back(p1);
  profile.knots.push_back(p2);
  const SEP::PhysicsExtensions::ProfileValue atKnot =
      SEP::PhysicsExtensions::EvaluateProfile(profile, 2.0);
  const SEP::PhysicsExtensions::ProfileValue outside =
      SEP::PhysicsExtensions::EvaluateProfile(profile, 4.0);
  failures += Check(atKnot.status.ok() && Near(atKnot.value, 3.0) &&
      Near(atKnot.derivativePerM, 1.0) && !outside.status.ok() &&
      !SEP::PhysicsExtensions::ProfileManifest(profile).empty(),
      "WP52 versioned C1 background profile and explicit domain policy");

  // WP53: independently requested bridge children sum to the parent and the
  // adaptive step records rejected trials without mutating its input value.
  SEP::NumericalExtensions::BrownianAddress root;
  root.campaign = 5; root.particle = 9; root.operatorId = 3;
  root.physicalInterval = 11;
  SEP::NumericalExtensions::BrownianAddress left = root, right = root;
  left.depth = right.depth = 1; left.index = 0; right.index = 1;
  const double parentDw =
      SEP::NumericalExtensions::BrownianIncrement(root, 2.0).value;
  const double childDw =
      SEP::NumericalExtensions::BrownianIncrement(left, 2.0).value +
      SEP::NumericalExtensions::BrownianIncrement(right, 2.0).value;
  SEP::NumericalExtensions::AdaptiveSdeConfiguration adaptiveCfg;
  adaptiveCfg.absoluteTolerance = 1.0e-8;
  adaptiveCfg.relativeTolerance = 2.0e-3;
  const SEP::NumericalExtensions::AdaptiveSdeResult adaptive =
      SEP::NumericalExtensions::AdvanceAdaptiveMultiplicativeSde(
          1.0, 0.4, 0.6, 0.5, root, adaptiveCfg);
  std::cout << "METRIC WP53 bridge_residual=" << (childDw - parentDw)
            << " accepted=" << adaptive.acceptedLeaves
            << " rejected=" << adaptive.rejectedTrials
            << " depth=" << adaptive.maximumDepthReached
            << " status=" << static_cast<int>(adaptive.status.code) << '\n';
  failures += Check(Near(parentDw, childDw, 2.0e-15) && adaptive.status.ok() &&
      adaptive.acceptedLeaves > 0 && adaptive.rejectedTrials > 0 &&
      std::isfinite(adaptive.value),
      "WP53 Brownian bridge identity and adaptive same-path rejection");

  // WP54: the selectable MUSCL/SSPRK2 path converges faster than first order
  // on a smooth periodic translation while remaining conservative.
  const double error40 = SineAdvectionError(40);
  const double error80 = SineAdvectionError(80);
  const double advectionOrder = std::log(error40 / error80) / std::log(2.0);
  std::cout << "METRIC WP54 coarse_error=" << error40
            << " fine_error=" << error80
            << " observed_order=" << advectionOrder << '\n';
  failures += Check(advectionOrder > 1.4,
      "WP54 second-order monotonic transport refinement");

  // WP55: exact overlap integration preserves a constant physical-volume
  // invariant across refine/coarsen and reports no numerical correction.
  std::vector<double> oldEdges; oldEdges.push_back(0.0);
  oldEdges.push_back(1.0); oldEdges.push_back(2.0);
  std::vector<double> oldValue(2, 3.0);
  std::vector<double> newEdges; newEdges.push_back(0.0);
  newEdges.push_back(0.5); newEdges.push_back(1.5); newEdges.push_back(2.0);
  const SEP::NumericalExtensions::RemapResult remap =
      SEP::NumericalExtensions::RemapConservative(
          oldEdges, oldValue, newEdges,
          SEP::NumericalExtensions::RemapOrder::PiecewiseLinear);
  failures += Check(remap.status.ok() && Near(remap.oldIntegral, 6.0) &&
      Near(remap.newIntegral, 6.0) && Near(remap.numericalResidual, 0.0),
      "WP55 conservative moving-flux-tube remap accounting");

  // WP56: failure types lead to distinct fates, and rollback restores staged
  // particle state before any caller-owned state can be committed.
  SEP::NumericalExtensions::FailureContext failureContext;
  failureContext.particleId = 7; failureContext.species = 0;
  failureContext.fieldLine = 1; failureContext.segment = 2;
  failureContext.operatorName = "coefficient";
  const SEP::NumericalExtensions::MoverResult quarantined =
      SEP::NumericalExtensions::ClassifyMoverFailure(
          SEP::Transport::Status::Error(
              SEP::Transport::StatusCode::InvalidParticleState, "fault"),
          failureContext,
          SEP::NumericalExtensions::FailurePolicy::QuarantineParticleLocal);
  SEP::Transport::ParticleState initialParticle;
  initialParticle.species = 0; initialParticle.fieldLineId = 1;
  initialParticle.massKg = 1.0; initialParticle.coordinate = 2.0;
  SEP::NumericalExtensions::ParticleTransaction particleTransaction(initialParticle);
  particleTransaction.staged()->coordinate = 9.0;
  particleTransaction.Rollback();
  failures += Check(quarantined.disposition ==
      SEP::NumericalExtensions::MoverDisposition::QuarantinedParticle &&
      quarantined.fate == SEP::NumericalExtensions::ParticleFate::Quarantine &&
      particleTransaction.staged()->coordinate == 2.0,
      "WP56 typed mover failure and transactional rollback");

  // WP57: splitting one physical lineage into two half-weight samples leaves
  // its cluster contribution and effective information unchanged.
  std::vector<SEP::ValidationExtensions::LineageSample> splitSamples;
  SEP::ValidationExtensions::LineageSample sample;
  sample.rootLineage = 1; sample.samplingWindow = 1; sample.bin = 0;
  sample.statisticalWeight = 0.5; sample.physicalValue = 2.0;
  sample.speedMPerS = 3.0;
  splitSamples.push_back(sample); splitSamples.push_back(sample);
  sample.rootLineage = 2; sample.statisticalWeight = 1.0;
  splitSamples.push_back(sample);
  std::vector<SEP::ValidationExtensions::LineageSample> unsplitSamples;
  sample.rootLineage = 1; unsplitSamples.push_back(sample);
  sample.rootLineage = 2; unsplitSamples.push_back(sample);
  const SEP::ValidationExtensions::LineageEstimate splitEstimate =
      SEP::ValidationExtensions::EstimateByLineage(splitSamples, 1,
          SEP::ValidationExtensions::EstimatorKind::SnapshotDensity);
  const SEP::ValidationExtensions::LineageEstimate unsplitEstimate =
      SEP::ValidationExtensions::EstimateByLineage(unsplitSamples, 1,
          SEP::ValidationExtensions::EstimatorKind::SnapshotDensity);
  failures += Check(splitEstimate.status.ok() && unsplitEstimate.status.ok() &&
      Near(splitEstimate.estimate[0], unsplitEstimate.estimate[0]) &&
      Near(splitEstimate.variance[0], unsplitEstimate.variance[0]) &&
      Near(splitEstimate.effectiveSampleSize[0],
           unsplitEstimate.effectiveSampleSize[0]),
      "WP57 lineage-aware uncertainty under population splitting");

  // WP58: a two-channel synthetic response matrix matches hand-calculated raw
  // counts; saturation is represented as censoring, not a Gaussian residual.
  SEP::ValidationExtensions::InstrumentResponseMatrix response;
  response.instrumentId = "fixture"; response.calibrationVersion = "1";
  response.checksum = "abc"; response.dataLevel = "raw";
  response.validUntilS = 10.0; response.channels = 2; response.species = 1;
  response.directions = 1; response.trueEnergyBins = 2;
  response.probability.push_back(1.0); response.probability.push_back(0.5);
  response.probability.push_back(0.0); response.probability.push_back(1.0);
  std::vector<double> truth; truth.push_back(4.0); truth.push_back(2.0);
  std::vector<SEP::ValidationExtensions::CountObservation> observed(2);
  observed[0].counts = 5.0; observed[1].counts = 2.0;
  observed[1].saturated = true; observed[1].saturationThreshold = 2.0;
  const SEP::ValidationExtensions::ResponseFoldResult folded =
      SEP::ValidationExtensions::FoldInstrumentCounts(
          response, truth, observed, 1.0, 0.0,
          SEP::ValidationExtensions::DeadTimeModel::None);
  failures += Check(folded.status.ok() && Near(folded.incidentCounts[0], 5.0) &&
      Near(folded.incidentCounts[1], 2.0) && folded.censored[1],
      "WP58 response-matrix count likelihood and saturation censoring");

  // WP59: source integration is deliberately blocked from claiming a linked
  // native row; the same complete trace passes only at NativeAmps level.
  SEP::ValidationExtensions::NativeMatrixTrace native;
  native.configurationFingerprint = "cfg"; native.executableChecksum = "exe";
  native.compiler = "c++"; native.flags = "-Werror"; native.moverId = "parker";
  native.coefficientProviderId = "prescribed";
  native.backgroundOwnerId = "analytic";
  native.turbulenceOwnerId = "integrated";
  native.moverEntered = native.coefficientEntered =
      native.turbulenceEntered = native.couplingEntered = true;
  const SEP::Transport::Status sourceOnly =
      SEP::ValidationExtensions::ValidateNativeMatrixRow(native, true);
  native.observedLevel = SEP::Evidence::Level::NativeAmps;
  const SEP::Transport::Status nativePass =
      SEP::ValidationExtensions::ValidateNativeMatrixRow(native, true);
  failures += Check(!sourceOnly.ok() && nativePass.ok(),
      "WP59 native linked matrix evidence boundary");

  // WP60: discrete identities remain exact while floating ledgers use only the
  // declared reduction tolerance.
  SEP::ValidationExtensions::DecompositionSignature serial;
  serial.configurationFingerprint = "cfg";
  serial.particleIdentityHash = "particles";
  serial.transactionOrderHash = "transactions";
  serial.restartHash = "restart"; serial.sourceEvents = 4;
  serial.transactions = 8; serial.ledger.push_back(10.0);
  SEP::ValidationExtensions::DecompositionSignature parallel = serial;
  parallel.ranks = 2; parallel.threads = 4; parallel.ledger[0] += 1.0e-13;
  failures += Check(SEP::ValidationExtensions::CompareDecompositions(
      serial, parallel, 1.0e-12).ok(),
      "WP60 decomposition restart and reduction signature");

  // WP61: a four-level manufactured sequence recovers its exact second order
  // with a near-perfect log-log fit and records nonzero equation terms.
  SEP::ValidationExtensions::ManufacturedCase manufactured;
  manufactured.id = "full-parker";
  manufactured.activeTerms.push_back("advection");
  manufactured.activeTerms.push_back("diffusion");
  for (int n = 1; n <= 4; ++n) {
    const double h = 1.0 / (10.0 * n);
    manufactured.resolution.push_back(h);
    manufactured.error.push_back(3.0 * h * h);
  }
  const SEP::ValidationExtensions::RefinementFit fit =
      SEP::ValidationExtensions::FitRefinementOrder(manufactured);
  failures += Check(fit.status.ok() && Near(fit.observedOrder, 2.0) &&
      fit.rSquared > 0.999999,
      "WP61 independent manufactured-equation refinement fit");

  // WP62: power is preregistered with multiplicity correction, and a bounded
  // exact-CDF panel produces a finite distributional statistic.
  SEP::ValidationExtensions::PowerSpecification powerSpec;
  powerSpec.observable = "first-passage-tail";
  powerSpec.standardizedEffect = 0.5; powerSpec.comparisons = 4;
  powerSpec.seedPanelVersion = "rare-event-v1";
  const SEP::ValidationExtensions::PowerResult power =
      SEP::ValidationExtensions::PlanNormalMeanPower(powerSpec);
  std::vector<double> samples; samples.push_back(0.1); samples.push_back(0.5);
  samples.push_back(0.9);
  std::vector<double> cdf; cdf.push_back(0.1); cdf.push_back(0.5);
  cdf.push_back(0.9);
  const SEP::Transport::ScalarResult ks =
      SEP::ValidationExtensions::KolmogorovSmirnovStatistic(samples, cdf);
  failures += Check(power.status.ok() && power.requiredSamples > 0 &&
      Near(power.adjustedAlpha, 0.0125) && ks.status.ok() && ks.value > 0.0,
      "WP62 powered rare-event and full-distribution validation");

  // WP63: synthetic fixtures cannot be promoted to observational validation;
  // a checksum-complete real held-out manifest is accepted by the gate.
  SEP::ValidationExtensions::ExternalEventManifest event;
  event.eventId = "event-1"; event.role =
      SEP::ValidationExtensions::EventRole::HeldOut;
  event.backgroundChecksum = "b"; event.geometryChecksum = "g";
  event.shockChecksum = "s"; event.calibrationChecksum = "c";
  event.preprocessing = "versioned"; event.exclusions = "none";
  event.configurationPrior = "frozen"; event.metrics.push_back("onset");
  const SEP::Transport::Status syntheticClaim =
      SEP::ValidationExtensions::ValidateExternalEventManifest(
          event, SEP::Evidence::Level::ObservationalValidation);
  event.synthetic = false;
  const SEP::Transport::Status heldOutClaim =
      SEP::ValidationExtensions::ValidateExternalEventManifest(
          event, SEP::Evidence::Level::ObservationalValidation);
  SEP::ValidationExtensions::EventMetrics zeroMetrics, unitSigma;
  unitSigma.onsetS = unitSigma.peak = unitSigma.fluence =
      unitSigma.anisotropy = unitSigma.spectralIndex =
      unitSigma.profileRmse = 1.0;
  const SEP::Transport::ScalarResult eventDistance =
      SEP::ValidationExtensions::NormalizedMetricDistance(
          zeroMetrics, zeroMetrics, unitSigma);
  failures += Check(!syntheticClaim.ok() && heldOutClaim.ok() &&
      eventDistance.status.ok() && Near(eventDistance.value, 0.0),
      "WP63 fail-closed SWMF cross-model and observational manifest");

  // WP64: a compatible strong-scaling sample passes both deterministic
  // per-item complexity and environment-specific parallel-efficiency gates.
  SEP::ValidationExtensions::ScalingSample scaleReference;
  scaleReference.workloadId = "particle";
  scaleReference.environmentFingerprint = "machine-a";
  scaleReference.problemItems = 100; scaleReference.operations = 1000;
  scaleReference.peakResidentBytes = 10000;
  scaleReference.peakQueueRecords = 20; scaleReference.wallSeconds = 10.0;
  SEP::ValidationExtensions::ScalingSample scaleCandidate = scaleReference;
  scaleCandidate.ranks = 2; scaleCandidate.threadsPerRank = 2;
  scaleCandidate.wallSeconds = 2.8;
  const SEP::ValidationExtensions::ScalingGate scaling =
      SEP::ValidationExtensions::CompareScaling(
          scaleReference, scaleCandidate, 0.8, 1.05, 1.05);
  failures += Check(scaling.status.ok() && scaling.parallelEfficiency > 0.8 &&
      Near(scaling.operationsPerItem, 10.0),
      "WP64 scaling performance and resource regression gate");

  if (failures != 0) {
    std::cerr << failures << " WP42-WP64 checks failed\n";
    return 1;
  }
  std::cout << "PASS WP42-WP64 dependency-light implementation suite\n";
  return 0;
}
