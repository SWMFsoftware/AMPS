#include "../util/sep_background_snapshot.h"
#include "../util/sep_focused_transport_core.h"
#include "../util/sep_focused_transport_mfp_core.h"
#include "../util/sep_reproducible_reduction.h"
#include "../util/sep_transport_common.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool Near(double a, double b, double tolerance) {
  return std::fabs(a - b) <= tolerance * std::max(1.0, std::fabs(b));
}

class RecordingPitchProvider final
    : public SEP::Transport::PitchAngleDiffusionProvider {
 public:
  mutable double lastLocationM = 0.0;
  SEP::Transport::PitchAngleDiffusionSample Evaluate(
      double sM, double, double) const override {
    lastLocationM = sM;
    SEP::Transport::PitchAngleDiffusionSample sample;
    sample.status = SEP::Transport::Status::Ok();
    sample.provenance = "WP08 recording provider";
    sample.turbulenceStateIdentity = "controlled:g1";
    return sample;
  }
};

class BranchProvider final : public SEP::Transport::MeanFreePathProvider {
 public:
  BranchProvider(double plus, double minus) : plus_(plus), minus_(minus) {}
  SEP::Transport::MeanFreePathSample Evaluate(
      double, double, double) const override {
    SEP::Transport::MeanFreePathSample sample;
    sample.status = SEP::Transport::Status::Ok();
    sample.lambdaParallelM = 1.0e7;
    sample.nuPlusPerS = plus_;
    sample.nuMinusPerS = minus_;
    sample.hasBranchResolvedRates = true;
    sample.provenance = "WP09/WP10 controlled branch rates";
    sample.turbulenceStateIdentity = "controlled:g1";
    return sample;
  }
 private:
  double plus_;
  double minus_;
};

int Check(bool condition, const char* label) {
  if (condition) {
    std::cout << "PASS " << label << '\n';
    return 0;
  }
  std::cerr << "FAIL " << label << '\n';
  return 1;
}

}  // namespace

int main() {
  using namespace SEP::Transport;
  int failures = 0;

  const ScalarResult snapshot = ComposeSnapshotValidityLimit(
      100.0, 2.0, 105.0, 7, 7);
  const ScalarResult stale = ComposeSnapshotValidityLimit(
      100.0, 2.0, 105.0, 7, 8);
  failures += Check(snapshot.status.ok() && snapshot.value == 3.0 &&
                    !stale.status.ok(), "WP04 snapshot limiter/generation");

  const SEP::Background::BackgroundSnapshot epochs(
      SEP::Background::Provider::Analytic,
      SEP::Background::Ownership::ModelOwned, 20.0, 20.0, 30.0, 4,
      "fingerprint", "controlled epochs", 10.0, 20.0);
  failures += Check(epochs.physical_epoch_interval_seconds() == 10.0,
                    "WP05 physical background epochs");

  FocusedTransportBackground full;
  full.dLnAbsBdsPerM = -2.0e-8;
  full.velocityDivergencePerS = 4.0e-4;
  full.fieldAlignedStrainPerS = 1.0e-4;
  full.equationMode = FocusedEquationMode::FullGyrotropic;
  const double mu = 0.25;
  const double speed = 2.0e6;
  const double expectedMuRate = 0.5 * (1.0 - mu * mu) *
      (-speed * full.dLnAbsBdsPerM +
       mu * (full.velocityDivergencePerS -
             3.0 * full.fieldAlignedStrainPerS));
  const double expectedLogP = -0.5 *
      ((1.0 - mu * mu) * full.velocityDivergencePerS +
       (3.0 * mu * mu - 1.0) * full.fieldAlignedStrainPerS);
  failures += Check(
      Near(FocusedPitchDriftPerS(mu, speed, full), expectedMuRate, 1.0e-14) &&
      Near(FocusedLogMomentumRatePerS(mu, full), expectedLogP, 1.0e-14),
      "WP06/WP07 authoritative gyrotropic coefficients");

  RecordingPitchProvider recording;
  KeyedRandomStream locationRandom(1, 2, 3, 4);
  FocusedTransportState locationState(100.0, 1.0e-20, 0.5);
  FocusedTransportBackground uniform;
  uniform.plasmaAdvectionMPerS = 10.0;
  const FocusedTransportIncrement locationStep = AdvanceFocusedTransportDmumu(
      locationState, uniform, 1.67262192369e-27, 299792458.0, 0.01,
      recording, locationRandom, NULL);
  failures += Check(locationStep.status.ok() && recording.lastLocationM > 100.0,
                    "WP08 location-aware midpoint sampling");

  const double mass = 1.67262192369e-27;
  const ScalarResult momentum = MomentumFromSpeed(1.0e6, mass, 299792458.0);
  FocusedTransportMfpBackground mfpBackground;
  mfpBackground.alfvenSpeedMPerS = 0.0;
  BranchProvider balanced(0.8, 0.8);
  KeyedRandomStream oneRandom(91, 77, 9, 0);
  FocusedTransportMfpState oneState(0.0, momentum.value, 0.3);
  const FocusedTransportMfpIncrement one = AdvanceFocusedTransportMfp(
      oneState, mfpBackground, mass, 299792458.0, 1.0, 0.3,
      balanced, oneRandom, NULL);

  KeyedRandomStream splitRandom(91, 77, 9, 0);
  FocusedTransportMfpState splitState(0.0, momentum.value, 0.3);
  const FocusedTransportMfpIncrement first = AdvanceFocusedTransportMfp(
      splitState, mfpBackground, mass, 299792458.0, 0.3, 0.3,
      balanced, splitRandom, NULL);
  const FocusedTransportMfpIncrement second = AdvanceFocusedTransportMfp(
      first.state, mfpBackground, mass, 299792458.0, 0.7, 0.3,
      balanced, splitRandom, NULL);
  failures += Check(one.status.ok() && first.status.ok() && second.status.ok() &&
      Near(one.state.arcLengthM, second.state.arcLengthM, 1.0e-12) &&
      Near(one.state.momentumKgMPerS, second.state.momentumKgMPerS, 1.0e-12) &&
      Near(one.state.mu, second.state.mu, 1.0e-12) &&
      Near(one.state.remainingOpticalDepth,
           second.state.remainingOpticalDepth, 1.0e-12),
      "WP09 carried optical-depth partition invariance");

  BranchProvider plusOnly(8.0, 0.0);
  KeyedRandomStream branchRandom(17, 18, 9, 0);
  const FocusedTransportMfpIncrement branches = AdvanceFocusedTransportMfp(
      FocusedTransportMfpState(0.0, momentum.value, 0.1), mfpBackground,
      mass, 299792458.0, 2.0, 0.1, plusOnly, branchRandom, NULL);
  failures += Check(branches.status.ok() &&
      branches.diagnostics.plusBranchEvents > 0 &&
      branches.diagnostics.minusBranchEvents == 0,
      "WP10 empty branch is never selected");

  SEP::Reproducibility::Contribution duplicate;
  duplicate.key.particle = 5;
  duplicate.waveEnergyJ = 1.0;
  std::vector<std::vector<SEP::Reproducibility::Contribution> > parts(2);
  parts[0].push_back(duplicate);
  parts[1].push_back(duplicate);
  std::vector<SEP::Reproducibility::SegmentAccumulator> reduced;
  failures += Check(!SEP::Reproducibility::CanonicalPartitionReduction(
                        parts, &reduced).ok(),
                    "WP02 duplicate physical keys are rejected");
  return failures == 0 ? 0 : 1;
}
