#include "sep_acceptance_cases.h"

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_background_snapshot.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_registry.h)
#include "sep_focused_transport_core.h"
#include "sep_focused_transport_mfp_core.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace SEP {
namespace Testing {
namespace {

const double kProtonMassKg = 1.67262192369e-27;
const double kSpeedOfLightMPerS = 299792458.0;

Result CheckedResult(bool passed, const std::string& pass_message,
                     const std::string& failure_message) {
  Result result;
  result.status = passed ? Status::Pass : Status::Fail;
  result.message = passed ? pass_message : failure_message;
  // RunOne audits this named counter.  Recording it in every new callback
  // means a future callback refactor cannot accidentally return PASS while an
  // assertion has already failed in a nested diagnostic.
  result.metrics.push_back(
      {"assertion_failures", passed ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

Descriptor MakeDescriptor(const char* id, const char* name, const char* group,
                          const char* description, RuntimeClass runtime,
                          TestCallback callback) {
  Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = group;
  descriptor.description = description;
  descriptor.initialization = InitializationLevel::None;
  descriptor.supportedBuildModes =
      "serial/MPI linked CLI and source-only sanitizer fixture";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = "deterministic keyed stream; fixed seed recorded";
  descriptor.stateIsolation =
      "stack-owned transport state; background store reset before and after fixture";
  descriptor.callback = callback;
  return descriptor;
}

class ZeroPitchAngleDiffusion final
    : public Transport::PitchAngleDiffusionProvider {
 public:
  Transport::PitchAngleDiffusionSample Evaluate(double, double,
                                                  double) const override {
    Transport::PitchAngleDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    sample.dMuMuPerS = 0.0;
    sample.dDmuMuDmuPerS = 0.0;
    sample.provenance = "acceptance:ballistic-dmumu-v1";
    sample.turbulenceStateIdentity = "acceptance:frozen-turbulence:g1";
    return sample;
  }
};

class InfiniteMeanFreePath final : public Transport::MeanFreePathProvider {
 public:
  Transport::MeanFreePathSample Evaluate(double, double,
                                          double) const override {
    Transport::MeanFreePathSample sample;
    sample.status = Transport::Status::Ok();
    sample.lambdaParallelM = std::numeric_limits<double>::infinity();
    sample.provenance = "acceptance:ballistic-mfp-v1";
    sample.turbulenceStateIdentity = "acceptance:frozen-turbulence:g1";
    return sample;
  }
};

bool NearlyEqual(double left, double right, double relative_tolerance) {
  const double scale = std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
  return std::fabs(left - right) <= relative_tolerance * scale;
}

Result RunAnalyticAndSwcmeFixture() {
  using namespace Background;
  const std::string analytic_fingerprint =
      FingerprintConfiguration("provider=analytic;fixture=BG01;version=1");
  const std::string swcme_fingerprint =
      FingerprintConfiguration("provider=swcme;fixture=BG01;version=1");

  // These are intentionally standalone metadata fixtures: they validate the
  // provider/epoch contract without relying on a previously generated AMPS
  // mesh or a mutable SWCME singleton.
  const BackgroundSnapshot analytic(
      Provider::Analytic, Ownership::ModelOwned, 12.0, 10.0, 20.0, 4,
      analytic_fingerprint, "standalone analytic solar-wind fixture");
  const BackgroundSnapshot swcme(
      Provider::Swcme, Ownership::ModelOwned, 12.0, 10.0, 20.0, 4,
      swcme_fingerprint, "standalone SWCME fixture");

  const bool passed = analytic.Covers(10.0) && analytic.Covers(20.0) &&
      swcme.Covers(12.0) && analytic.provider() == Provider::Analytic &&
      swcme.provider() == Provider::Swcme &&
      analytic.ownership() == Ownership::ModelOwned &&
      swcme.ownership() == Ownership::ModelOwned &&
      analytic.configuration_fingerprint() !=
          swcme.configuration_fingerprint();
  Result result = CheckedResult(
      passed,
      "standalone analytic and SWCME fixtures preserve provider, epoch, and ownership",
      "analytic/SWCME background fixture violated provider or epoch metadata");
  result.configuration.push_back("validity_interval_s=[10,20]");
  result.configuration.push_back("field_line_generation=4");
  result.metrics.push_back({"providers_checked", 2.0, 2.0, "==", "count"});
  return result;
}

Result RunMockSwmfFixture() {
  using namespace Background;
  SnapshotStore& store = SnapshotStore::Instance();
  store.ResetForTests();
  bool passed = true;
  std::string failure;
  try {
    const std::string fingerprint =
        FingerprintConfiguration("provider=swmf;fixture=BG02;version=1");
    store.Publish(BackgroundSnapshot(
        Provider::Swmf, Ownership::ImportedReadOnly, 100.0, 100.0, 110.0,
        8, fingerprint, "mock SWMF import; no external coupler required"));
    {
      ParticleReadPhase read = store.BeginParticleRead(105.0);
      passed = passed && read.snapshot().provider() == Provider::Swmf &&
          read.snapshot().ownership() == Ownership::ImportedReadOnly;
    }
    store.PublishHandoff(BackgroundSnapshot(
        Provider::LocalEvolution, Ownership::HandoffCopy, 110.0, 110.0,
        120.0, 9, fingerprint, "explicit copy from mock SWMF generation 8"),
        Provider::Swmf);
    passed = passed && store.Current()->provider() == Provider::LocalEvolution &&
        store.Current()->ownership() == Ownership::HandoffCopy &&
        store.Current()->field_line_generation() == 9;
  }
  catch (const std::exception& exception) {
    passed = false;
    failure = exception.what();
  }
  store.ResetForTests();

  Result result = CheckedResult(
      passed,
      "mock SWMF import remains read-only and transfers through an explicit handoff copy",
      "mock SWMF ownership/handoff fixture failed: " + failure);
  result.configuration.push_back("external_SWMF_dependency=none");
  result.metrics.push_back({"final_generation", passed ? 9.0 : -1.0,
                            9.0, "==", "generation"});
  return result;
}

Result RunBallisticCrossMoverLimit() {
  // Both focused movers consume the same manufactured background and begin
  // from exactly the same SI state.  D_mumu=0 and lambda=+infinity are the two
  // APIs' explicit representations of the same ballistic physical limit.
  const Transport::FocusedTransportState dmumu_initial(50.0, 2.0e-19, 0.35);
  const Transport::FocusedTransportBackground dmumu_background(
      0.0, 4.0e5, 0.0, 0.0);
  ZeroPitchAngleDiffusion dmumu_provider;
  Transport::KeyedRandomStream dmumu_random(1301, 77, 13, 0);
  const Transport::FocusedTransportIncrement dmumu =
      Transport::AdvanceFocusedTransportDmumu(
          dmumu_initial, dmumu_background, kProtonMassKg,
          kSpeedOfLightMPerS, 2.0, dmumu_provider, dmumu_random, NULL);

  const Transport::FocusedTransportMfpState mfp_initial(50.0, 2.0e-19, 0.35);
  const Transport::FocusedTransportMfpBackground mfp_background(
      0.0, 4.0e5, 0.0, 0.0, 5.0e4);
  InfiniteMeanFreePath mfp_provider;
  Transport::KeyedRandomStream mfp_random(1301, 77, 13, 0);
  const Transport::FocusedTransportMfpIncrement mfp =
      Transport::AdvanceFocusedTransportMfp(
          mfp_initial, mfp_background, kProtonMassKg,
          kSpeedOfLightMPerS, 2.0, 2.0, mfp_provider, mfp_random, NULL);

  const bool passed = dmumu.status.ok() && mfp.status.ok() &&
      NearlyEqual(dmumu.state.arcLengthM, mfp.state.arcLengthM, 2.0e-14) &&
      NearlyEqual(dmumu.state.momentumKgMPerS,
                  mfp.state.momentumKgMPerS, 2.0e-14) &&
      NearlyEqual(dmumu.state.mu, mfp.state.mu, 2.0e-14) &&
      mfp.diagnostics.scatteringEvents == 0 &&
      mfp.diagnostics.ballisticIntervals == 1;
  Result result = CheckedResult(
      passed,
      "fte-dmumu and fte-mfp agree in their matched ballistic limit",
      "focused movers diverged for Dmumu=0 and lambda=infinity");
  result.hasSeed = true;
  result.seed = 1301;
  result.configuration.push_back("Dmumu_s^-1=0");
  result.configuration.push_back("lambda_parallel_m=+infinity");
  result.metrics.push_back({"arc_length_difference_m",
      std::fabs(dmumu.state.arcLengthM - mfp.state.arcLengthM), 2.0e-14,
      "relative-scale <=", "m"});
  return result;
}

Result RunMatchedDiffusionClosure() {
  // The comparison is performed at the coefficient boundary, where Parker's
  // kappa, fte-mfp's lambda, and fte-dmumu's D_mumu have an exact common
  // isotropic-scattering closure.  This avoids mistaking finite-time Monte
  // Carlo noise for a disagreement between mover implementations.
  const double speed = 1.2e7;
  const double lambda = 2.4e9;
  const double mu = 0.3;
  const Transport::ScalarResult kappa =
      Transport::Coefficient::KappaFromMeanFreePath(lambda, speed);
  const Transport::ScalarResult lambda_from_kappa =
      Transport::Coefficient::MeanFreePathFromKappa(kappa.value, speed);
  const Transport::ScalarResult dmumu =
      Transport::Coefficient::IsotropicDmumuFromMeanFreePath(
          lambda, speed, mu);
  const Transport::ScalarResult lambda_from_dmumu =
      Transport::Coefficient::MeanFreePathFromIsotropicDmumu(
          dmumu.value, speed, mu);
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  const double relative_kappa_error =
      std::fabs(lambda_from_kappa.value / lambda - 1.0);
  const double relative_dmumu_error =
      std::fabs(lambda_from_dmumu.value / lambda - 1.0);
  const bool passed = kappa.status.ok() && lambda_from_kappa.status.ok() &&
      dmumu.status.ok() && lambda_from_dmumu.status.ok() &&
      relative_kappa_error <= tolerance &&
      relative_dmumu_error <= tolerance;
  Result result = CheckedResult(
      passed,
      "Parker, fte-dmumu, and fte-mfp coefficients share the declared isotropic closure",
      "cross-mover coefficient conversion failed its round-trip tolerance");
  result.configuration.push_back("speed_m_per_s=1.2e7");
  result.configuration.push_back("lambda_parallel_m=2.4e9");
  result.configuration.push_back("mu=0.3");
  result.metrics.push_back({"lambda_roundtrip_from_kappa_relative_error",
      relative_kappa_error, tolerance, "<=", "dimensionless"});
  result.metrics.push_back({"lambda_roundtrip_from_Dmumu_relative_error",
      relative_dmumu_error, tolerance, "<=", "dimensionless"});
  return result;
}

}  // namespace

std::vector<Descriptor> AcceptanceCaseDescriptors() {
  std::vector<Descriptor> descriptors;
  descriptors.push_back(MakeDescriptor(
      "BG01", "Analytic and SWCME background fixtures", "background",
      "Validate standalone analytic and SWCME provider/epoch metadata.",
      RuntimeClass::Routine, RunAnalyticAndSwcmeFixture));
  descriptors.push_back(MakeDescriptor(
      "BG02", "Mock SWMF ownership fixture", "background",
      "Validate read-only SWMF import and the explicit local-evolution handoff.",
      RuntimeClass::Routine, RunMockSwmfFixture));
  descriptors.push_back(MakeDescriptor(
      "CROSS01", "Focused-mover ballistic agreement", "cross-mover",
      "Compare fte-dmumu and fte-mfp under matched no-scattering assumptions.",
      RuntimeClass::Routine, RunBallisticCrossMoverLimit));
  descriptors.push_back(MakeDescriptor(
      "CROSS02", "Three-mover diffusion closure", "cross-mover",
      "Round-trip matched Parker, Dmumu, and mean-free-path coefficients.",
      RuntimeClass::Routine, RunMatchedDiffusionClosure));
  return descriptors;
}

}  // namespace Testing
}  // namespace SEP
