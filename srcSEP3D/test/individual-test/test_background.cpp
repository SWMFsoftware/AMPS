// Phase B analytic-background, coupled-import, and immutable-snapshot tests.
// These tests are AMPS-independent: coupled values enter as explicit records,
// exactly as they do after the SWMF receive boundary.

#include "../../background/background_snapshot.h"
#include "../../background/bg_parker.h"
#include "../../background/bg_swmf.h"
#include "../../core/sep3d_test_registry.h"
#include "../../runtime/runtime_adapters.h"

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace B = SEP3D::Background;
namespace RM = SEP3D::RuntimeModel;
using Result = SEP3D::Testing::Result;

Result Pass(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Pass;
  result.message = message; return result;
}
Result Fail(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Fail;
  result.message = message; return result;
}

B::ParkerConfiguration ParkerConfig(int polarity = 1) {
  B::ParkerConfiguration configuration;
  configuration.magneticPolarity = polarity;
  configuration.validityCadenceS = 100.0;
  return configuration;
}

bool PreparedParker(B::AnalyticParkerProvider* provider, double epoch = 0.0) {
  return provider != nullptr && provider->Prepare(epoch).ok();
}

// Mutable fixture provider used only to inject one bad field at a time.  It
// still goes through the production batch/builder contract; tests do not call
// BackgroundSnapshot's constructor to bypass validation.
class FixtureProvider final : public B::BackgroundProvider {
 public:
  FixtureProvider(const B::BackgroundSample& sample,
                  const B::SnapshotMetadata& metadata)
      : sample_(sample), metadata_(metadata) {}
  const char* CanonicalName() const override { return "phase-b-fixture"; }
  SEP3D::Core::Status Validate() const override {
    return SEP3D::Core::Status::OK();
  }
  SEP3D::Core::Status Prepare(double) override {
    prepared_ = true; return SEP3D::Core::Status::OK();
  }
  const B::SnapshotMetadata* PreparedMetadata() const override {
    return prepared_ ? &metadata_ : nullptr;
  }
  B::BackgroundSample Evaluate(const SEP3D::Core::Vec3&) const override {
    return sample_;
  }
  std::string ResolvedManifest() const override { return "phase-b-fixture-v1"; }
  B::ProviderCapabilities Capabilities() const override {
    B::ProviderCapabilities result;
    result.hasAnalyticGradB = true;
    result.hasAnalyticDivBhat = true;
    result.hasAnalyticCurvature = true;
    result.hasAnalyticDivU = true;
    result.hasFieldAlignedStrain = true;
    result.hasPlasmaState = true;
    result.supportsBatchEval = true;
    return result;
  }
  void set_sample(const B::BackgroundSample& sample) { sample_ = sample; }

 private:
  B::BackgroundSample sample_;
  B::SnapshotMetadata metadata_;
  bool prepared_ = true;
};

B::SnapshotMetadata Metadata(B::ProviderKind kind, std::uint64_t generation,
                             double epoch) {
  B::SnapshotMetadata result;
  result.provider = kind;
  result.ownership = kind == B::ProviderKind::AnalyticParker
      ? B::StorageOwnership::ModelOwned
      : B::StorageOwnership::ImportedReadOnly;
  result.epochS = epoch;
  result.validFromS = epoch;
  result.validUntilS = epoch + 10.0;
  result.generation = generation;
  result.coordinateFrame = "HCI-like-inertial";
  result.providerIdentity = kind == B::ProviderKind::AnalyticParker
      ? "analytic-fixture" : "swmf-fixture";
  result.configurationFingerprint = "fixture-fingerprint";
  return result;
}

B::BackgroundSample CompleteSample(std::uint64_t generation = 1) {
  B::BackgroundSample sample;
  sample.status = SEP3D::Core::Status::OK();
  sample.valid = true;
  sample.B = {3.0e-9, -2.0e-9, 1.0e-9};
  sample.absB = sample.B.Norm();
  sample.bHat = sample.B.Normalized();
  sample.U = {4.0e5, 1.0e4, -2.0e4};
  sample.numberDensityM3 = 5.0e6;
  sample.temperatureK = 1.0e5;
  sample.pressurePa = sample.numberDensityM3 * SEP3D::Core::Const::k_B *
                      sample.temperatureK;
  sample.alfvenSpeedMpS = 5.0e4;
  sample.divU = 1.0e-6;
  sample.divBhat = -2.0e-12;
  sample.focusingLenM = 5.0e10;
  sample.curvature = {1.0e-12, 2.0e-12, -1.0e-12};
  sample.fieldAlignedStrain = 2.0e-7;
  for (int i = 0; i < 3; ++i) {
    sample.gradB(i, i) = (i - 1) * 1.0e-20;
    sample.gradU(i, i) = (i + 1) * 1.0e-7;
  }
  sample.generation = generation;
  sample.configurationDigest = 47;
  return sample;
}

B::SwmfRawSample RawFromSample(const B::BackgroundSample& sample,
                               const SEP3D::Core::Vec3& positionM,
                               double epoch, B::SwmfUnitSystem units) {
  const double length = units == B::SwmfUnitSystem::SI
                            ? 1.0 : SEP3D::Core::Const::R_sun;
  const double magnetic = units == B::SwmfUnitSystem::SI ? 1.0 : 1.0e-9;
  const double velocity = units == B::SwmfUnitSystem::SI ? 1.0 : 1.0e3;
  const double density = units == B::SwmfUnitSystem::SI ? 1.0 : 1.0e6;
  const double pressure = units == B::SwmfUnitSystem::SI ? 1.0 : 1.0e-9;
  B::SwmfRawSample raw;
  raw.position = positionM / length;
  raw.magnetic = sample.B / magnetic;
  raw.velocity = sample.U / velocity;
  raw.numberDensity = sample.numberDensityM3 / density;
  raw.temperature = sample.temperatureK;
  raw.pressure = sample.pressurePa / pressure;
  raw.velocityDivergence = sample.divU;
  raw.divBhat = sample.divBhat * length;
  raw.focusingLength = sample.focusingLenM / length;
  raw.curvature = sample.curvature * length;
  raw.fieldAlignedStrain = sample.fieldAlignedStrain;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      raw.magneticGradient(i, j) = sample.gradB(i, j) * length / magnetic;
      raw.velocityGradient(i, j) = sample.gradU(i, j) * length / velocity;
    }
  raw.epochS = epoch;
  raw.complete = true;
  return raw;
}

B::SwmfImport Import(B::SwmfUnitSystem units,
                     const B::BackgroundSample& sample,
                     const SEP3D::Core::Vec3& position, double epoch = 5.0,
                     std::uint64_t generation = 1) {
  B::SwmfImport imported;
  imported.units = units;
  imported.epochS = epoch;
  imported.validUntilS = epoch + 10.0;
  imported.generation = generation;
  imported.configurationFingerprint = "swmf-config-v1";
  imported.samples.push_back(RawFromSample(sample, position, epoch, units));
  return imported;
}

bool Nearly(double left, double right, double tolerance) {
  return std::fabs(left - right) <= tolerance *
      std::max(1.0e-300, std::max(std::fabs(left), std::fabs(right)));
}

Result RunBGP3D01() {
  double medianAccumulator = 0.0;
  double maximum = 0.0;
  std::size_t count = 0;
  for (int polarity : {-1, 1}) {
    B::AnalyticParkerProvider provider(ParkerConfig(polarity));
    if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
    for (int i = 1; i <= 400; ++i) {
      const double radius = (0.1 + 0.9 * i / 400.0) * SEP3D::Core::Const::AU;
      const double theta = 1.0e-5 + (SEP3D::Core::Const::kPi - 2.0e-5) *
                           i / 400.0;
      const double phi = 0.371 * i;
      const SEP3D::Core::Vec3 x(radius * std::sin(theta) * std::cos(phi),
                                radius * std::sin(theta) * std::sin(phi),
                                radius * std::cos(theta));
      const B::BackgroundSample sample = provider.Evaluate(x);
      if (!sample.status.ok()) return Fail("Parker evaluation failed");
      const double normalized = std::fabs(sample.gradB.Trace()) /
                                (sample.absB / radius);
      medianAccumulator += normalized;
      maximum = std::max(maximum, normalized);
      ++count;
    }
  }
  const double meanUpperBound = medianAccumulator / count;
  if (meanUpperBound > 1.0e-8 || maximum > 1.0e-6)
    return Fail("analytic Parker field is not divergence-free within tolerance");
  Result result = Pass("analytic Parker gradient is divergence-free over deterministic axes/polarities");
  result.metrics.push_back({"mean_normalized_divB", meanUpperBound, 1.0e-8, "<=", ""});
  result.metrics.push_back({"maximum_normalized_divB", maximum, 1.0e-6, "<=", ""});
  return result;
}

Result RunBGP3D02() {
  const B::ParkerConfiguration configuration = ParkerConfig();
  B::AnalyticParkerProvider provider(configuration);
  if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
  for (double radius : {0.2 * SEP3D::Core::Const::AU,
                        0.7 * SEP3D::Core::Const::AU,
                        SEP3D::Core::Const::AU}) {
    const SEP3D::Core::Vec3 x(radius, 0.0, 0.0);
    const B::BackgroundSample sample = provider.Evaluate(x);
    const double br = sample.B.x;
    const double bphi = sample.B.y;
    const double expectedBr = configuration.radialFieldAtReferenceT *
        std::pow(configuration.referenceRadiusM / radius, 2);
    const double expectedRatio = -configuration.solarRotationRateRadPerS *
        (radius - configuration.sourceRadiusM) /
        configuration.solarWindSpeedMPerS;
    if (!Nearly(br, expectedBr, 1.0e-12) ||
        !Nearly(bphi / br, expectedRatio, 1.0e-12))
      return Fail("Parker radial or azimuthal component law changed");
  }
  return Pass("Parker radial scaling and azimuthal/radial ratio match closed forms to 1e-12 relative");
}

Result RunBGP3D03() {
  B::AnalyticParkerProvider provider(ParkerConfig());
  if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
  double worst = 0.0;
  for (int i = 1; i <= 100; ++i) {
    const double radius = (0.15 + 0.8 * i / 100.0) * SEP3D::Core::Const::AU;
    const B::BackgroundSample sample = provider.Evaluate({radius, 0.0, 0.0});
    const double ratio = sample.B.y / sample.B.x;
    const SEP3D::Core::Vec3 tangent =
        SEP3D::Core::Vec3(1.0, ratio, 0.0).Normalized();
    // For nearly parallel unit vectors acos(dot) magnifies round-off by
    // roughly 1/sqrt(1-dot).  The cross-product norm is sin(angle), is
    // first-order accurate at zero, and therefore tests the requested angular
    // tolerance without manufacturing a ~sqrt(epsilon) apparent error.
    worst = std::max(worst, tangent.Cross(sample.bHat).Norm());
  }
  if (worst > 1.0e-8) return Fail("analytic field-line tangent differs from B");
  return Pass("analytic Parker field-line tangents align with the local field below 1e-8 rad");
}

Result RunBGP3D04() {
  B::AnalyticParkerProvider provider(ParkerConfig());
  if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
  const SEP3D::Core::Vec3 x(0.73 * SEP3D::Core::Const::AU,
                            0.11 * SEP3D::Core::Const::AU,
                            0.19 * SEP3D::Core::Const::AU);
  const B::BackgroundSample center = provider.Evaluate(x);
  const double h = 1.0e-5 * x.Norm();
  const double plus = provider.Evaluate(x + h * center.bHat).absB;
  const double minus = provider.Evaluate(x - h * center.bHat).absB;
  const double derivative = (std::log(plus) - std::log(minus)) / (2.0 * h);
  const double numerical = -1.0 / derivative;
  if (!Nearly(numerical, center.focusingLenM, 1.0e-7))
    return Fail("analytic and converged numerical focusing lengths differ");
  return Pass("analytic focusing length matches a converged field-aligned derivative to 1e-7 relative");
}

Result RunBGP3D05() {
  const B::ParkerConfiguration configuration = ParkerConfig();
  B::AnalyticParkerProvider provider(configuration);
  if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
  const SEP3D::Core::Vec3 x(0.5 * SEP3D::Core::Const::AU,
                            0.25 * SEP3D::Core::Const::AU,
                            -0.1 * SEP3D::Core::Const::AU);
  const B::BackgroundSample sample = provider.Evaluate(x);
  const double expected = 2.0 * configuration.solarWindSpeedMPerS / x.Norm();
  if (!Nearly(sample.divU, expected, 1.0e-12) ||
      !Nearly(sample.gradU.Trace(), expected, 1.0e-12))
    return Fail("radial-wind divergence or gradient trace is incorrect");
  return Pass("radial-wind analytic gradient and divergence agree with 2V/r");
}

Result RunBGP3D06() {
  B::AnalyticParkerProvider provider(ParkerConfig());
  if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
  const double radius = SEP3D::Core::Const::AU;
  for (double angle : {0.0, 1.0e-12, SEP3D::Core::Const::kPi - 1.0e-12,
                       SEP3D::Core::Const::kPi}) {
    const B::BackgroundSample sample = provider.Evaluate(
        {radius * std::sin(angle), 0.0, radius * std::cos(angle)});
    if (!sample.status.ok() || !std::isfinite(sample.absB) ||
        !std::isfinite(sample.divBhat) ||
        !std::isfinite(sample.focusingLenM))
      return Fail("Parker polar limit is non-finite");
  }
  return Pass("Parker provider remains finite on and within 1e-12 rad of both rotation-axis poles");
}

Result RunBGP3D07() {
  B::ParkerConfiguration configuration = ParkerConfig();
  configuration.thermodynamicClosure =
      swcme::solarwind::ThermodynamicClosure::MultiSpecies;
  configuration.alphaToProtonRatio = 0.05;
  configuration.electronTemperatureK = 2.0e5;
  configuration.alphaTemperatureK = 3.0e5;
  B::AnalyticParkerProvider provider(configuration);
  if (!PreparedParker(&provider))
    return Fail("SWCME-backed multi-species Parker provider did not prepare");

  const B::BackgroundSample oneAu = provider.Evaluate(
      {SEP3D::Core::Const::AU, 0.0, 0.0});
  const B::BackgroundSample inner = provider.Evaluate(
      {0.1 * SEP3D::Core::Const::AU, 0.0, 0.0});
  if (!oneAu.status.ok() || !inner.status.ok() ||
      !Nearly(oneAu.numberDensityM3,
              configuration.numberDensityAtReferenceM3, 2.0e-14)) {
    return Fail("Leblanc density did not retain its explicit one-AU normalization");
  }

  // At 0.1 AU, the positive r^-4 and r^-6 Leblanc terms make the density
  // strictly larger than a pure n(1 AU)*(AU/r)^2 law. This catches a return to
  // the former local r^-2 approximation without duplicating SWCME's constants
  // in the test.
  const double oldPowerLaw = configuration.numberDensityAtReferenceM3 * 100.0;
  if (!(inner.numberDensityM3 > oldPowerLaw))
    return Fail("near-Sun density collapsed to the retired r^-2 approximation");

  const double electronDensity = inner.numberDensityM3;
  const double protonDensity = electronDensity /
      (1.0 + 2.0 * configuration.alphaToProtonRatio);
  const double alphaDensity =
      configuration.alphaToProtonRatio * protonDensity;
  const double expectedPressure = SEP3D::Core::Const::k_B *
      (protonDensity * configuration.temperatureK +
       electronDensity * configuration.electronTemperatureK +
       alphaDensity * configuration.alphaTemperatureK);
  const double expectedMassDensity =
      swcme::constants::PROTON_MASS_KG * protonDensity +
      swcme::constants::ALPHA_PARTICLE_MASS_KG * alphaDensity;
  const double mu0 = 4.0e-7 * SEP3D::Core::Const::kPi;
  const double expectedAlfven =
      inner.absB / std::sqrt(mu0 * expectedMassDensity);
  if (!Nearly(inner.pressurePa, expectedPressure, 2.0e-14) ||
      !Nearly(inner.alfvenSpeedMpS, expectedAlfven, 2.0e-14)) {
    return Fail("multi-species pressure or mass-density Alfvén speed differs from SWCME closure");
  }
  return Pass(
      "Parker initialization uses SWCME Leblanc density and charge-neutral multi-species thermodynamics");
}

Result RunSNAP3D01() {
  const SEP3D::Core::Vec3 point(SEP3D::Core::Const::AU, 0.0, 0.0);
  const B::BackgroundSample good = CompleteSample();
  FixtureProvider provider(good, Metadata(B::ProviderKind::AnalyticParker, 1, 0.0));
  B::BackgroundSnapshotBuilder builder;
  std::shared_ptr<const B::BackgroundSnapshot> active;
  if (!builder.Build(provider, {point}, &active).ok())
    return Fail("complete fixture snapshot was rejected");
  const B::BackgroundSnapshot* identity = active.get();
  std::vector<B::BackgroundSample> bad;
  B::BackgroundSample sample = good; sample.valid = false; bad.push_back(sample);
  sample = good; sample.absB = 0.0; bad.push_back(sample);
  sample = good; sample.numberDensityM3 = 0.0; bad.push_back(sample);
  sample = good; sample.temperatureK = 0.0; bad.push_back(sample);
  sample = good; sample.pressurePa = 0.0; bad.push_back(sample);
  sample = good; sample.alfvenSpeedMpS = 0.0; bad.push_back(sample);
  sample = good; sample.generation = 0; bad.push_back(sample);
  for (const auto& candidate : bad) {
    provider.set_sample(candidate);
    if (builder.Build(provider, {point}, &active).ok() || active.get() != identity)
      return Fail("incomplete candidate replaced the active snapshot");
  }
  return Pass("omitting each required primitive rejects the candidate and preserves the active snapshot");
}

Result RunSNAP3D02() {
  B::BackgroundSample sample = CompleteSample();
  const B::ProviderCapabilities capabilities = FixtureProvider(
      sample, Metadata(B::ProviderKind::AnalyticParker, 1, 0.0)).Capabilities();
  double* scalars[] = {&sample.absB, &sample.numberDensityM3,
                       &sample.temperatureK, &sample.pressurePa,
                       &sample.alfvenSpeedMpS, &sample.divU,
                       &sample.divBhat, &sample.focusingLenM,
                       &sample.fieldAlignedStrain};
  for (double* scalar : scalars) {
    const double saved = *scalar;
    *scalar = std::numeric_limits<double>::quiet_NaN();
    if (B::ValidateCompleteSample(sample, capabilities).ok())
      return Fail("NaN required field was accepted");
    *scalar = std::numeric_limits<double>::infinity();
    if (B::ValidateCompleteSample(sample, capabilities).ok())
      return Fail("infinite required field was accepted");
    *scalar = saved;
  }
  return Pass("NaN and infinity in required scalar fields are rejected without invented replacements");
}

Result RunSNAP3D03() {
  const SEP3D::Core::Vec3 position(0.8 * SEP3D::Core::Const::AU,
                                   0.1 * SEP3D::Core::Const::AU, 0.0);
  const B::BackgroundSample expected = CompleteSample();
  B::SwmfAwsomProvider si;
  B::SwmfAwsomProvider coupling;
  if (!si.Load(Import(B::SwmfUnitSystem::SI, expected, position)).ok() ||
      !coupling.Load(Import(B::SwmfUnitSystem::AwsomCoupling,
                            expected, position)).ok() ||
      !si.Prepare(5.0).ok() || !coupling.Prepare(5.0).ok())
    return Fail("SI/coupling import preparation failed");
  const B::BackgroundSample left = si.Evaluate(position);
  const B::BackgroundSample right = coupling.Evaluate(position);
  if (!Nearly(left.absB, right.absB, 5.0e-15) ||
      !Nearly(left.numberDensityM3, right.numberDensityM3, 5.0e-15) ||
      !Nearly(left.U.x, right.U.x, 5.0e-15) ||
      !Nearly(left.gradB(1, 1), right.gradB(1, 1), 5.0e-15))
    return Fail("documented coupling units do not reproduce SI state");
  return Pass("SI and documented AWSoM coupling units publish round-off-equivalent snapshots");
}

Result RunSNAP3D04() {
  const SEP3D::Core::Vec3 position(SEP3D::Core::Const::AU, 0.0, 0.0);
  B::SwmfImport imported = Import(B::SwmfUnitSystem::SI,
                                  CompleteSample(), position, 5.0);
  imported.samples.push_back(imported.samples.front());
  imported.samples.back().position.y += 1.0;
  imported.samples.back().epochS = 6.0;
  B::SwmfAwsomProvider provider;
  if (provider.Load(imported).ok()) return Fail("mixed source epochs were accepted");
  return Pass("mixed SWMF record epochs are rejected before provider publication");
}

Result RunSNAP3D05() {
  const SEP3D::Core::Vec3 position(SEP3D::Core::Const::AU, 0.0, 0.0);
  B::SwmfAwsomProvider provider;
  const B::BackgroundSample expected = CompleteSample();
  if (!provider.Load(Import(B::SwmfUnitSystem::SI, expected, position)).ok() ||
      !provider.Prepare(5.0).ok()) return Fail("valid SWMF fixture failed");
  const B::SnapshotMetadata* before = provider.PreparedMetadata();
  const std::uint64_t generation = before->generation;
  B::SwmfImport invalid = Import(B::SwmfUnitSystem::SI,
                                 expected, position, 8.0, 2);
  invalid.samples.front().pressure = std::numeric_limits<double>::quiet_NaN();
  if (provider.Load(invalid).ok() || provider.PreparedMetadata() == nullptr ||
      provider.PreparedMetadata()->generation != generation ||
      provider.Evaluate(position).pressurePa != expected.pressurePa)
    return Fail("failed SWMF candidate partially mutated the active state");
  return Pass("failed late-field validation preserves the active SWMF identifier and values atomically");
}

Result RunSNAP3D06() {
  const SEP3D::Core::Vec3 point(SEP3D::Core::Const::AU, 0.0, 0.0);
  B::BackgroundSample first = CompleteSample(1);
  B::BackgroundSample second = CompleteSample(2);
  second.numberDensityM3 = 3.0 * first.numberDensityM3;
  second.pressurePa = 3.0 * first.pressurePa;
  std::shared_ptr<const B::BackgroundSnapshot> left(new B::BackgroundSnapshot(
      Metadata(B::ProviderKind::AnalyticParker, 1, 0.0),
      FixtureProvider(first, Metadata(B::ProviderKind::AnalyticParker, 1, 0.0)).Capabilities(),
      {point}, {first}));
  B::SnapshotMetadata rightMetadata = Metadata(B::ProviderKind::AnalyticParker, 2, 10.0);
  rightMetadata.configurationFingerprint = left->metadata().configurationFingerprint;
  std::shared_ptr<const B::BackgroundSnapshot> right(new B::BackgroundSnapshot(
      rightMetadata, left->capabilities(), {point}, {second}));
  B::SnapshotBuffer buffer;
  if (!buffer.PublishCurrent(left).ok() || !buffer.StageNext(right).ok())
    return Fail("could not stage compatible snapshot pair");
  std::shared_ptr<const B::BackgroundSnapshot> middle;
  if (!buffer.SnapshotAt(5.0, &middle).ok() ||
      middle->samples().front().numberDensityM3 !=
          2.0 * first.numberDensityM3 ||
      buffer.SnapshotAt(11.0, &middle).ok())
    return Fail("linear interpolation or extrapolation rejection failed");
  return Pass("two frozen snapshots interpolate linear fields exactly and reject extrapolation");
}

Result RunSNAP3D07() {
  B::AnalyticParkerProvider provider(ParkerConfig());
  if (!PreparedParker(&provider)) return Fail("Parker preparation failed");
  const double x[2] = {SEP3D::Core::Const::AU, 0.5 * ParkerConfig().sourceRadiusM};
  const double y[2] = {0.0, 0.0};
  const double z[2] = {0.0, 0.0};
  B::BackgroundSample output[2] = {CompleteSample(99), CompleteSample(99)};
  const B::BackgroundSample sentinel = output[1];
  SEP3D::Core::Status status[2];
  if (provider.EvaluateBatchDetailed(x, y, z, 2, output, status).ok() ||
      !status[0].ok() || status[1].ok() ||
      std::memcmp(&output[1].numberDensityM3, &sentinel.numberDensityM3,
                  sizeof(double)) != 0 || output[1].generation != sentinel.generation)
    return Fail("batch evaluation lost per-sample status or overwrote failure output");
  return Pass("mixed valid/invalid batch preserves each status and leaves failed output unchanged");
}

Result RunSNAP3D08() {
  const SEP3D::Core::Vec3 position(SEP3D::Core::Const::AU, 0.0, 0.0);
  B::SwmfImport imported = Import(B::SwmfUnitSystem::SI,
                                  CompleteSample(), position);
  imported.coordinateFrame = "GSE";
  B::SwmfAwsomProvider provider;
  if (provider.Load(imported).ok())
    return Fail("mismatched SWMF frame was accepted without a transformation");
  return Pass("mismatched coordinate frame is rejected unless an explicit transformation is supplied");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterBackgroundTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using R = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* group, const char* name,
                 const char* description, SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = group; d.description = description;
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = R::Routine; d.seedPolicy = "deterministic-no-rng";
    d.stateIsolation = "fresh provider, builder, and snapshot buffer per test";
    d.callback = std::move(callback); return d;
  };
  return {
      make("BGP3D01", "BGP3D", "Divergence-free Parker field", "Analytic Cartesian gradient trace.", RunBGP3D01),
      make("BGP3D02", "BGP3D", "Parker component laws", "Radial and spiral closed forms.", RunBGP3D02),
      make("BGP3D03", "BGP3D", "Field-line tangency", "Parker curve tangent and local B.", RunBGP3D03),
      make("BGP3D04", "BGP3D", "Focusing length", "Analytic and numerical field-aligned derivative.", RunBGP3D04),
      make("BGP3D05", "BGP3D", "Velocity derivatives", "Radial wind gradient and divergence.", RunBGP3D05),
      make("BGP3D06", "BGP3D", "Polar limits", "Rotation-axis finite limits.", RunBGP3D06),
      make("BGP3D07", "BGP3D", "SWCME ambient closure", "Leblanc density and multi-species thermodynamics.", RunBGP3D07),
      make("SNAP3D01", "SNAP3D", "Completeness", "Required fields and atomic candidate rejection.", RunSNAP3D01),
      make("SNAP3D02", "SNAP3D", "Finite-value policy", "NaN/Inf rejection without replacement.", RunSNAP3D02),
      make("SNAP3D03", "SNAP3D", "Unit conversion", "SI and AWSoM coupling-unit equivalence.", RunSNAP3D03),
      make("SNAP3D04", "SNAP3D", "Epoch consistency", "Mixed source epochs fail before publication.", RunSNAP3D04),
      make("SNAP3D05", "SNAP3D", "Atomic publication", "Failed update preserves active imported state.", RunSNAP3D05),
      make("SNAP3D06", "SNAP3D", "Time interpolation", "Linear interpolation and no extrapolation.", RunSNAP3D06),
      make("SNAP3D07", "SNAP3D", "Batch status", "Per-sample failure status and unchanged output.", RunSNAP3D07),
      make("SNAP3D08", "SNAP3D", "Coordinate frame", "Mismatched frame rejection.", RunSNAP3D08),
  };
}
