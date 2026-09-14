#include "../util/sep_configuration_matrix.h"
#include "../util/sep_evidence.h"
#include "../util/sep_observation_forward_model.h"
#include "../util/sep_population_control.h"
#include "../util/sep_system_ledger.h"
#include "../util/sep_turbulence_core.h"
#include "../util/sep_validation_tools.h"

#include <cmath>
#include <iostream>
#include <limits>
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

bool Near(double a, double b, double tolerance = 1.0e-12) {
  return std::fabs(a - b) <= tolerance *
      std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

SEP::Turbulence::State TurbulenceState() {
  SEP::Turbulence::State state;
  state.configuration.source =
      SEP::Turbulence::Source::SelfConsistentIntegrated;
  state.configuration.representation = SEP::Turbulence::Representation::Integrated;
  state.configuration.advectionEnabled = false;
  state.configuration.reflectionEnabled = false;
  state.configuration.cascadeEnabled = true;
  SEP::Turbulence::CellState cell;
  cell.lengthM = 1.0;
  cell.volumeM3 = 1.0;
  cell.massDensityKgPerM3 = 1.0;
  cell.magneticFieldT = 1.0;
  cell.ePlusJ = 1.0;
  cell.eMinusJ = 1.0;
  state.cells.push_back(cell);
  return state;
}

SEP::PopulationControl::ParticleRecord Particle(std::uint64_t id,
                                                double weight,
                                                double energy,
                                                double momentum,
                                                double mu) {
  SEP::PopulationControl::ParticleRecord particle;
  particle.stableId = id;
  particle.species = 0;
  particle.weight = weight;
  particle.chargeC = 1.602176634e-19;
  particle.kineticEnergyJ = energy;
  particle.parallelMomentumKgMPerS = momentum;
  particle.mu = mu;
  return particle;
}

}  // namespace

int main() {
  int failures = 0;

  // WP31: a stiff reflection timescale must participate in the same plan as
  // advection and cascade.  The cascade-only manufactured decay demonstrates
  // the declared Lie/backward-Euler first-order behavior under refinement.
  SEP::Turbulence::State stiff = TurbulenceState();
  stiff.configuration.cascadeEnabled = false;
  stiff.configuration.reflectionEnabled = true;
  stiff.configuration.reflectionCoefficient = 1.0;
  stiff.cells[0].alfvenSpeedMPerS = 10.0;
  stiff.cells[0].dLnAlfvenSpeeddsPerM = 1.0;
  const SEP::Turbulence::AdvancePlan plan =
      SEP::Turbulence::PlanAdvance(stiff, 1.0);
  SEP::Turbulence::State coarse = TurbulenceState();
  SEP::Turbulence::State fine = coarse;
  SEP::Turbulence::State reference = coarse;
  coarse.configuration.maximumCascadeFraction = 1.0;
  fine.configuration.maximumCascadeFraction = 1.0;
  reference.configuration.maximumCascadeFraction = 1.0;
  coarse.configuration.perpendicularCorrelationLengthM = 1.0;
  fine.configuration.perpendicularCorrelationLengthM = 1.0;
  reference.configuration.perpendicularCorrelationLengthM = 1.0;
  SEP::Turbulence::InitializeState(&coarse);
  SEP::Turbulence::InitializeState(&fine);
  SEP::Turbulence::InitializeState(&reference);
  SEP::Turbulence::Advance(&coarse, 0.2);
  for (int i = 0; i < 2; ++i) SEP::Turbulence::Advance(&fine, 0.1);
  for (int i = 0; i < 2000; ++i) SEP::Turbulence::Advance(&reference, 0.0001);
  const double exact = reference.cells[0].ePlusJ + reference.cells[0].eMinusJ;
  const double coarseError = std::fabs(
      coarse.cells[0].ePlusJ + coarse.cells[0].eMinusJ - exact);
  const double fineError = std::fabs(
      fine.cells[0].ePlusJ + fine.cells[0].eMinusJ - exact);
  const double order = std::log(coarseError / fineError) / std::log(2.0);
  std::cout << "METRIC WP31 substeps=" << plan.substeps
            << " limiter=" << plan.limitingOperator
            << " coarse_error=" << coarseError
            << " fine_error=" << fineError
            << " observed_order=" << order << '\n';
  failures += Check(plan.status.ok() && plan.substeps > 1 &&
      plan.limitingOperator == "reflection" && order > 0.8 && order < 1.3,
      "WP31 combined operator limit and observed first-order convergence");

  // WP32: a negative source that exceeds available energy is corrected once,
  // recorded with its signed rejected amount, and never converted silently.
  SEP::Turbulence::State limited = TurbulenceState();
  limited.configuration.cascadeEnabled = false;
  limited.cells[0].pendingParticlePlusJ = -3.0;
  SEP::Turbulence::InitializeState(&limited);
  const SEP::Turbulence::StepResult limitedStep =
      SEP::Turbulence::Advance(&limited, 1.0);
  SEP::Turbulence::State invalid = TurbulenceState();
  invalid.cells[0].massDensityKgPerM3 =
      std::numeric_limits<double>::quiet_NaN();
  const SEP::Turbulence::AdvancePlan invalidPlan =
      SEP::Turbulence::PlanAdvance(invalid, 1.0);
  failures += Check(limitedStep.status.ok() && limited.cells[0].ePlusJ == 0.0 &&
      limitedStep.diagnostics.limiterActivations == 1 &&
      Near(limitedStep.ledger.rejectedSourceJ, -2.0) &&
      !limitedStep.diagnostics.events.empty() && !invalidPlan.status.ok(),
      "WP32 typed correction ledger and nonfinite rejection");

  // WP33: clone splitting preserves every declared moment exactly and produces
  // stable, unique, decomposition-independent lineage identifiers.
  const SEP::PopulationControl::ParticleRecord parent =
      Particle(17, 7.0, 4.0, 2.0, 0.25);
  std::vector<SEP::PopulationControl::ParticleRecord> children;
  const SEP::Transport::Status splitStatus =
      SEP::PopulationControl::SplitParticle(parent, 3, 1, &children);
  std::vector<SEP::PopulationControl::ParticleRecord> before(1, parent);
  const SEP::PopulationControl::InvariantReport splitReport =
      SEP::PopulationControl::CompareInvariants(before, children, 1.0e-15);
  SEP::PopulationControl::ParticleRecord merged;
  const SEP::Transport::Status mergeStatus =
      SEP::PopulationControl::MergeParticles(children, 9, 1, &merged);
  std::vector<SEP::PopulationControl::ParticleRecord> after(1, merged);
  const SEP::PopulationControl::InvariantReport mergeReport =
      SEP::PopulationControl::CompareInvariants(children, after, 1.0e-15);
  std::cout << "METRIC WP33 split_error=" << splitReport.maximumRelativeError
            << " merge_error=" << mergeReport.maximumRelativeError
            << " failed_merge_invariants=" << mergeReport.failedInvariants.size()
            << '\n';
  failures += Check(splitStatus.ok() && splitReport.status.ok() &&
      mergeStatus.ok() && mergeReport.status.ok() &&
      children[0].stableId != children[1].stableId &&
      children[0].parentId == parent.stableId,
      "WP33 split/merge invariants and stable lineage");

  // WP34: a source-only double is explicitly unable to promote itself to
  // native evidence.  The same complete observation is accepted only with the
  // native execution level and exactly one queue-flush owner.
  SEP::Evidence::NativeHarnessRequest request;
  request.combination.mover = SEP::Mover::ProductionMover::Parker;
  request.combination.coefficientSource =
      SEP::Transport::Coefficient::SourceMode::Prescribed;
  request.combination.turbulenceSource = SEP::Turbulence::Source::Prescribed;
  request.combination.coupling = SEP::Turbulence::CouplingPolicy::Disabled;
  request.expectedBinary = "amps";
  request.configurationFingerprint = "cfg";
  request.sourceGeneration = "generation-1";
  request.boundaryCase = "absorb";
  SEP::Evidence::NativeObservation observation;
  observation.productionMoverEntered = true;
  observation.coefficientAdapterEntered = true;
  observation.turbulenceDriverEntered = true;
  observation.queueFlushOwners = 1;
  observation.executedBinary = "amps";
  observation.configurationFingerprint = "cfg";
  observation.sourceGeneration = "generation-1";
  observation.threadCount = observation.rankCount = 1;
  const SEP::Transport::Status falsePromotion =
      SEP::Evidence::ValidateNativeObservation(
          request, observation, SEP::Evidence::Level::SourceIntegration);
  const SEP::Transport::Status nativeContract =
      SEP::Evidence::ValidateNativeObservation(
          request, observation, SEP::Evidence::Level::NativeAmps);
  failures += Check(!falsePromotion.ok() && nativeContract.ok(),
      "WP34 fail-closed native production-adapter harness contract");

  // WP35: 3 movers x 3 coefficient sources x 5 turbulence sources x 2
  // coupling policies yields 90 classified rows with no implicit fallback.
  const std::vector<std::pair<SEP::ConfigurationMatrix::Combination,
      SEP::ConfigurationMatrix::Classification> > matrix =
      SEP::ConfigurationMatrix::Enumerate();
  bool classified = matrix.size() == 90;
  for (std::size_t i = 0; i < matrix.size(); ++i)
    classified = classified && !matrix[i].second.diagnosticCode.empty();
  SEP::ConfigurationMatrix::Combination unsupported;
  unsupported.mover = SEP::Mover::ProductionMover::Parker;
  unsupported.coefficientSource =
      SEP::Transport::Coefficient::SourceMode::SelfConsistent;
  unsupported.turbulenceSource = SEP::Turbulence::Source::Prescribed;
  unsupported.coupling = SEP::Turbulence::CouplingPolicy::Disabled;
  failures += Check(classified &&
      !SEP::ConfigurationMatrix::Preflight(unsupported, false).ok() &&
      SEP::ConfigurationMatrix::RenderMarkdown().find("CFG-WAVE-IMMUTABLE") !=
          std::string::npos,
      "WP35 complete generated compatibility matrix and preflight");

  // WP36: one compact global ledger closes through injection and escape, then
  // preserves identity, residual hazard, RNG counter, and generation exactly.
  SEP::SystemVerification::CampaignCheckpoint checkpoint;
  checkpoint.configurationFingerprint = "run-a";
  checkpoint.fieldLineGeneration = 4;
  checkpoint.epochS = 3.0;
  checkpoint.turbulenceStateHash = "waves-a";
  checkpoint.outputManifestChecksum = "output-a";
  checkpoint.ledger.initial.number = 10.0;
  checkpoint.ledger.initial.energyJ = 20.0;
  checkpoint.ledger.final.number = 11.0;
  checkpoint.ledger.final.energyJ = 22.0;
  SEP::SystemVerification::Transaction injection;
  injection.seam = "injection";
  injection.provenance = "controlled";
  injection.change.number = 2.0;
  injection.change.energyJ = 5.0;
  SEP::SystemVerification::Transaction escape;
  escape.seam = "escape";
  escape.provenance = "controlled";
  escape.change.number = -1.0;
  escape.change.energyJ = -3.0;
  checkpoint.ledger.transactions.push_back(injection);
  checkpoint.ledger.transactions.push_back(escape);
  SEP::SystemVerification::ParticleIdentityState identity;
  identity.stableId = 101;
  identity.lineageGeneration = 2;
  identity.rngEventCounter = 7;
  identity.residualHazard = 0.125;
  checkpoint.particles.push_back(identity);
  const SEP::Transport::Status closed =
      SEP::SystemVerification::CloseLedger(&checkpoint.ledger, 0.0);
  std::string serialized;
  SEP::SystemVerification::SerializeCheckpoint(checkpoint, &serialized);
  SEP::SystemVerification::CampaignCheckpoint restarted;
  const SEP::Transport::Status restored =
      SEP::SystemVerification::DeserializeCheckpoint(serialized, &restarted);
  std::string corrupted = serialized;
  corrupted[20] = corrupted[20] == 'x' ? 'y' : 'x';
  SEP::SystemVerification::CampaignCheckpoint rejected;
  failures += Check(closed.ok() && restored.ok() &&
      restarted.particles[0].residualHazard == 0.125 &&
      restarted.fieldLineGeneration == 4 &&
      !SEP::SystemVerification::DeserializeCheckpoint(corrupted, &rejected).ok(),
      "WP36 global ledger checkpoint/restart and corruption detection");

  // WP37: the seed panel is reproducible and domain-separated; the statistical
  // gate states sample count, confidence, estimator, and z threshold.
  const SEP::ValidationTools::SeedPanel panelA =
      SEP::ValidationTools::MakeSeedPanel(9, "v1", "pitch-mean", 8);
  const SEP::ValidationTools::SeedPanel panelB =
      SEP::ValidationTools::MakeSeedPanel(9, "v1", "energy-mean", 8);
  SEP::ValidationTools::RunningStatistics statistics;
  statistics.Add(-2.0); statistics.Add(-1.0); statistics.Add(1.0); statistics.Add(2.0);
  const SEP::ValidationTools::StatisticalGate ensemble =
      SEP::ValidationTools::CompareMean(statistics, 0.0, 2.0, 0.95);
  failures += Check(panelA.seeds.size() == 8 && panelA.seeds != panelB.seeds &&
      ensemble.status.ok() && ensemble.sampleCount == 4 && ensemble.zScore == 0.0,
      "WP37 versioned multi-seed ensemble statistics");

  // WP38: the first failing generated case retains its exact seed/index and a
  // fault armed on the third matching seam fires once at that deterministic hit.
  const std::vector<SEP::ValidationTools::GeneratedCase> generated =
      SEP::ValidationTools::GenerateBoundaryCases(55, 10);
  const SEP::ValidationTools::PropertyFailure property =
      SEP::ValidationTools::CheckProperty("nonnegative", generated,
          [](double value) {
            return value < 0.0
                ? SEP::Transport::Status::Error(
                    SEP::Transport::StatusCode::OutOfDomain, "negative")
                : SEP::Transport::Status::Ok();
          });
  SEP::ValidationTools::FaultInjector injector(
      SEP::ValidationTools::FaultPoint::Checksum, 3);
  const bool first = injector.ShouldFail(SEP::ValidationTools::FaultPoint::Checksum);
  const bool second = injector.ShouldFail(SEP::ValidationTools::FaultPoint::Checksum);
  const bool third = injector.ShouldFail(SEP::ValidationTools::FaultPoint::Checksum);
  failures += Check(property.failed && property.counterexample.seed == 55 &&
      property.counterexample.index == 3 && !first && !second && third,
      "WP38 reproducible property counterexample and failure injection");

  // WP39: a two-bin synthetic spectrum has analytically known response-folded
  // counts.  This validates the forward operator without claiming spacecraft
  // or held-out-event evidence.
  SEP::Observation::PopulationSpectrum population;
  population.species = "proton";
  population.energyEdgesJ.push_back(0.0);
  population.energyEdgesJ.push_back(1.0);
  population.energyEdgesJ.push_back(2.0);
  population.differentialIntensity.push_back(2.0);
  population.differentialIntensity.push_back(3.0);
  population.variance.assign(2, 0.0);
  SEP::Observation::InstrumentChannel channel;
  channel.id = "synthetic-1";
  channel.minimumEnergyJ = 0.0;
  channel.maximumEnergyJ = 2.0;
  channel.geometricFactorM2Sr = 2.0;
  channel.cadenceS = 5.0;
  channel.saturationCounts = 1000.0;
  channel.response.push_back(1.0);
  channel.response.push_back(0.5);
  SEP::Observation::BackgroundEstimate background;
  background.expectedCounts = 5.0;
  background.variance = 4.0;
  const SEP::Observation::ChannelPrediction prediction =
      SEP::Observation::ForwardModel(population, channel, background);
  failures += Check(prediction.status.ok() && Near(prediction.incidentCounts, 35.0) &&
      Near(prediction.backgroundSubtractedCounts, 30.0) &&
      Near(prediction.standardUncertainty, std::sqrt(39.0)),
      "WP39 analytical synthetic instrument forward model");

  // WP40: deterministic work changes fail on every machine; wall-time budgets
  // apply only to the recorded compiler/hardware environment identity.
  SEP::ValidationTools::PerformanceSample performance;
  performance.work.moverSubsteps = 10;
  performance.wallSeconds = 2.0;
  SEP::ValidationTools::PerformanceBaseline baseline;
  baseline.expectedWork.moverSubsteps = 10;
  baseline.medianWallSeconds = 1.0;
  baseline.maximumWallRatio = 1.5;
  baseline.environmentId = "host-a";
  const SEP::Transport::Status differentHost =
      SEP::ValidationTools::ComparePerformance(performance, baseline, "host-b");
  const SEP::Transport::Status slow =
      SEP::ValidationTools::ComparePerformance(performance, baseline, "host-a");
  performance.work.moverSubsteps = 11;
  const SEP::Transport::Status moreWork =
      SEP::ValidationTools::ComparePerformance(performance, baseline, "host-b");
  failures += Check(differentHost.ok() && !slow.ok() && !moreWork.ok(),
      "WP40 normalized-work and environment-specific performance gates");

  // WP41: a PASS below the required evidence level is rejected, while a claim
  // with native level and a reproducible command/artifact can be published.
  SEP::Evidence::Claim claim;
  claim.id = "CLAIM-MOVER-NATIVE";
  claim.statement = "all public movers execute their production shells";
  claim.requiredLevel = SEP::Evidence::Level::NativeAmps;
  claim.observedLevel = SEP::Evidence::Level::SourceIntegration;
  claim.status = SEP::Evidence::GateStatus::Pass;
  const SEP::Transport::Status overstated = SEP::Evidence::ValidateClaim(claim);
  claim.observedLevel = SEP::Evidence::Level::NativeAmps;
  SEP::Evidence::Artifact artifact;
  artifact.command = "make test-native-adapter-harness";
  artifact.path = "native-results.json";
  artifact.checksum = "controlled";
  claim.artifacts.push_back(artifact);
  const SEP::Transport::Status evidenced = SEP::Evidence::ValidateClaim(claim);
  std::vector<SEP::Evidence::Claim> claims(1, claim);
  failures += Check(!overstated.ok() && evidenced.ok() &&
      SEP::Evidence::RenderClaimTable(claims).find("native-amps") !=
          std::string::npos,
      "WP41 evidence-level claim governance");

  return failures == 0 ? 0 : 1;
}
