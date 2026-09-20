// ============================================================================
// Phase-T turbulence/scattering acceptance tests.
//
// The tests cover authority and units before testing numerical coefficients:
// normalized prescribed spectra, AWSoM direction mapping under both magnetic
// polarities, explicit finite-band policy, and explicit missing-data policy.
// COEF3D then proves that srcSEP3D calls the canonical sep_common conversions
// and Dmumu kernel rather than a copied implementation.
// ============================================================================

#include "../../core/sep3d_test_registry.h"
#include "../../turbulence/coefficient_bridge.h"
#include "../../turbulence/turbulence_models.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace B = SEP3D::Background;
namespace T = SEP3D::Turbulence;
namespace CP = SEP::Transport::CoefficientPhysics;
using Result = SEP3D::Testing::Result;

Result Pass(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Pass;
  result.message = message;
  return result;
}

Result Fail(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Fail;
  result.message = message;
  return result;
}

bool Relative(double left, double right, double tolerance) {
  const double scale =
      std::max(1.0e-300, std::max(std::fabs(left), std::fabs(right)));
  return std::fabs(left - right) <= tolerance * scale;
}

B::BackgroundSample Background(double radialSign = 1.0) {
  B::BackgroundSample sample;
  sample.status = SEP3D::Core::Status::OK();
  sample.valid = true;
  sample.B = {radialSign * 5.0e-9, 0.0, 0.0};
  sample.absB = 5.0e-9;
  sample.bHat = {radialSign, 0.0, 0.0};
  sample.alfvenSpeedMpS = 5.0e4;
  sample.generation = 3;
  sample.configurationDigest = 17;
  return sample;
}

T::TurbulenceSample PrescribedSample() {
  T::PrescribedKolmogorovProvider provider(
      T::PrescribedKolmogorovConfiguration{});
  if (!provider.Prepare(0.0).ok()) return {};
  return provider.Evaluate(
      {SEP3D::Core::Const::AU, 0.0, 0.0}, Background());
}

T::AwsomTurbulenceImport AwsomImport(bool complete = true) {
  T::AwsomTurbulenceImport imported;
  imported.epochS = 12.0;
  imported.validUntilS = 22.0;
  imported.generation = 4;
  imported.configurationFingerprint = "awsom-test-generation-4";
  T::AwsomWaveRecord record;
  record.positionM = {SEP3D::Core::Const::AU, 0.0, 0.0};
  record.epochS = imported.epochS;
  record.wPlusJPerM3 = 4.0;
  record.wMinusJPerM3 = 1.0;
  record.complete = complete;
  imported.records.push_back(record);
  return imported;
}

CP::SpeciesProperties Proton() {
  CP::SpeciesProperties species;
  species.modelSpecies = 0;
  species.name = "proton";
  species.signedChargeC = SEP3D::Core::Const::e;
  species.restMassKg = SEP3D::Core::Const::m_p;
  species.nucleonCount = 1.0;
  return species;
}

CP::SpeciesProperties Electron() {
  CP::SpeciesProperties species;
  species.modelSpecies = 1;
  species.name = "electron";
  species.signedChargeC = -SEP3D::Core::Const::e;
  species.restMassKg = SEP3D::Core::Const::m_e;
  species.nucleonCount = 0.0;
  return species;
}

Result RunTUR3D01() {
  const T::TurbulenceSample sample = PrescribedSample();
  if (!sample.status.ok() || !sample.valid)
    return Fail("prescribed turbulence provider failed");
  const T::NormalizedPowerLawSpectrum spectrum(sample);
  if (!spectrum.Validate().ok()) return Fail("normalized spectrum is invalid");

  // Composite Simpson quadrature in u=ln(k) resolves all spectral decades
  // uniformly.  Since dk=exp(u)du, the transformed integrand is P(k)*k.
  const int intervals = 20000;  // even by construction
  const double lower = std::log(sample.kMinPerM);
  const double upper = std::log(sample.kMaxPerM);
  const double h = (upper - lower) / intervals;
  double sum = 0.0;
  for (int i = 0; i <= intervals; ++i) {
    // Use the declared endpoints exactly.  exp(log(k_edge)) can land one ulp
    // outside a closed interval, which would turn a quadrature representation
    // artifact into an apparent resonance-policy failure.
    const double waveNumber = i == 0 ? sample.kMinPerM
        : (i == intervals ? sample.kMaxPerM
                          : std::exp(lower + i * h));
    const T::SpectrumValue value = spectrum.Evaluate(
        waveNumber, T::ResonanceRangePolicy::Reject);
    if (!value.status.ok()) return Fail("in-band spectrum evaluation failed");
    const double integrand = value.valueT2M * waveNumber;
    sum += (i == 0 || i == intervals) ? integrand
           : (i % 2 == 0 ? 2.0 : 4.0) * integrand;
  }
  const double integrated = h * sum / 3.0;
  if (!Relative(integrated, sample.deltaB2T2, 1.0e-10) ||
      !Relative(spectrum.AnalyticBandVarianceT2(), sample.deltaB2T2,
                1.0e-14)) {
    return Fail("Kolmogorov normalization does not close to total variance");
  }
  Result result = Pass(
      "prescribed Kolmogorov spectrum integrates to its declared magnetic variance");
  result.metrics.push_back({"relative_quadrature_error",
      std::fabs(integrated / sample.deltaB2T2 - 1.0), 1.0e-10, "<=", ""});
  return result;
}

Result RunTUR3D02() {
  constexpr double mu0 = 4.0e-7 * SEP3D::Core::Const::kPi;
  T::AwsomTurbulenceProvider provider;
  if (!provider.Load(AwsomImport()).ok() || !provider.Prepare(12.0).ok())
    return Fail("valid AWSoM wave record did not prepare");
  const SEP3D::Core::Vec3 position(SEP3D::Core::Const::AU, 0.0, 0.0);
  const T::TurbulenceSample outward = provider.Evaluate(position, Background(1.0));
  const T::TurbulenceSample inward = provider.Evaluate(position, Background(-1.0));
  if (!outward.status.ok() || !inward.status.ok() ||
      !Relative(outward.deltaBPlus2T2, mu0 * 4.0, 1.0e-15) ||
      !Relative(outward.deltaBMinus2T2, mu0, 1.0e-15) ||
      outward.deltaBOutward2T2 != outward.deltaBPlus2T2 ||
      outward.deltaBInward2T2 != outward.deltaBMinus2T2 ||
      inward.deltaBOutward2T2 != inward.deltaBMinus2T2 ||
      inward.deltaBInward2T2 != inward.deltaBPlus2T2) {
    return Fail("AWSoM energy-unit or polarity-direction mapping is incorrect");
  }
  return Pass(
      "AWSoM w+/w- map through deltaB^2=mu0*w and swap outward/inward labels with polarity");
}

Result RunTUR3D03() {
  for (const CP::SpeciesProperties& species : {Proton(), Electron()}) {
    const SEP::Transport::ScalarResult gyro =
        CP::GyrofrequencyRadPerS(5.0e-9, species);
    if (!gyro.status.ok()) return Fail("species gyrofrequency failed");
    const double speed = species.name == "proton" ? 3.0e7 : 1.0e8;
    const double mu = 0.5;
    const double resonance = gyro.value / (speed * std::fabs(mu));

    T::TurbulenceSample sample = PrescribedSample();
    sample.kMinPerM = 0.5 * resonance;
    sample.kMaxPerM = 2.0 * resonance;
    T::NormalizedPowerLawSpectrum spectrum(sample);
    const double probes[] = {0.25 * resonance, resonance, 4.0 * resonance};
    for (int i = 0; i < 3; ++i) {
      const T::SpectrumValue rejected = spectrum.Evaluate(
          probes[i], T::ResonanceRangePolicy::Reject);
      const T::SpectrumValue extended = spectrum.Evaluate(
          probes[i], T::ResonanceRangePolicy::PowerLawExtension);
      const bool outside = i != 1;
      if (rejected.status.ok() == outside || !extended.status.ok() ||
          extended.extended != outside ||
          extended.resolvedWaveNumberPerM != probes[i]) {
        return Fail("resonance range policy was not applied exactly");
      }
    }
  }
  return Pass(
      "proton/electron resonances below, inside, and above the band obey reject/extension policies");
}

Result RunTUR3D04() {
  const SEP3D::Core::Vec3 position(SEP3D::Core::Const::AU, 0.0, 0.0);
  T::AwsomTurbulenceProvider fail(T::MissingTurbulencePolicy::Fail);
  T::AwsomTurbulenceProvider ballistic(T::MissingTurbulencePolicy::Ballistic);
  if (!fail.Load(AwsomImport(false)).ok() || !fail.Prepare(12.0).ok() ||
      !ballistic.Load(AwsomImport(false)).ok() ||
      !ballistic.Prepare(12.0).ok()) {
    return Fail("incomplete-record fixture could not reach policy boundary");
  }
  const T::TurbulenceSample failed = fail.Evaluate(position, Background());
  const T::TurbulenceSample free = ballistic.Evaluate(position, Background());
  const CP::PitchAngleResult failCoefficient = T::CoefficientBridge::JokipiiDmumu(
      failed, Background(), position.Norm(), CP::SpectrumParameters{},
      Proton(), 3.0e7, 0.5);
  const CP::PitchAngleResult freeCoefficient = T::CoefficientBridge::JokipiiDmumu(
      free, Background(), position.Norm(), CP::SpectrumParameters{},
      Proton(), 3.0e7, 0.5);
  if (failed.status.usable() || failCoefficient.status.ok() ||
      !free.status.ballistic() || !free.valid ||
      !freeCoefficient.status.ok() ||
      freeCoefficient.valueState != CP::ValueState::Ballistic ||
      freeCoefficient.dMuMuPerS != 0.0) {
    return Fail("missing turbulence did not preserve explicit fail/ballistic policy");
  }
  return Pass(
      "missing wave data fails by default and becomes zero-rate ballistic only when explicitly selected");
}

Result RunCOEF3D01() {
  const double mus[] = {-0.75, 0.0, 0.65};
  const double speeds[] = {1.0e6, 1.0e7, 1.0e8};
  double maximumError = 0.0;
  for (int decade = 0; decade <= 6; ++decade) {
    const double lambda = 1.0e7 * std::pow(10.0, decade);
    for (double speed : speeds) {
      const SEP::Transport::ScalarResult kappa =
          T::CoefficientBridge::KappaFromMeanFreePath(lambda, speed);
      const SEP::Transport::ScalarResult lambdaBack =
          T::CoefficientBridge::MeanFreePathFromKappa(kappa.value, speed);
      if (!kappa.status.ok() || !lambdaBack.status.ok())
        return Fail("lambda/kappa conversion returned an error");
      maximumError = std::max(maximumError,
                              std::fabs(lambdaBack.value / lambda - 1.0));
      for (double mu : mus) {
        const SEP::Transport::ScalarResult dmumu =
            T::CoefficientBridge::IsotropicDmumuFromMeanFreePath(
                lambda, speed, mu);
        const SEP::Transport::ScalarResult fromDmumu =
            T::CoefficientBridge::MeanFreePathFromIsotropicDmumu(
                dmumu.value, speed, mu);
        if (!dmumu.status.ok() || !fromDmumu.status.ok())
          return Fail("lambda/Dmumu conversion returned an error");
        maximumError = std::max(maximumError,
                                std::fabs(fromDmumu.value / lambda - 1.0));
      }
    }
  }
  if (maximumError > 1.0e-12)
    return Fail("coefficient round-trip error exceeds 1e-12");
  Result result = Pass(
      "lambda, kappa_parallel, and isotropic Dmumu round-trip over six decades");
  result.metrics.push_back(
      {"maximum_relative_error", maximumError, 1.0e-12, "<=", ""});
  return result;
}

Result RunCOEF3D02() {
  const T::TurbulenceSample turbulence = PrescribedSample();
  const B::BackgroundSample background = Background();
  CP::SpectrumParameters spectrum;
  const CP::SpeciesProperties species = Proton();
  const double radius = SEP3D::Core::Const::AU;
  const double speed = 3.0e7;
  const double mu = 0.4;
  const CP::PitchAngleResult through3D = T::CoefficientBridge::JokipiiDmumu(
      turbulence, background, radius, spectrum, species, speed, mu);
  const CP::PitchAngleResult direct = CP::EvaluateJokipiiSlab(
      T::CoefficientBridge::ToSharedInput(turbulence, background, radius),
      spectrum, species, speed, mu);
  if (through3D.status.code != direct.status.code ||
      through3D.valueState != direct.valueState ||
      through3D.derivativeState != direct.derivativeState ||
      std::memcmp(&through3D.dMuMuPerS, &direct.dMuMuPerS,
                  sizeof(double)) != 0 ||
      std::memcmp(&through3D.dDmuMuDmuPerS, &direct.dDmuMuDmuPerS,
                  sizeof(double)) != 0) {
    return Fail("srcSEP3D adapter changed a canonical sep_common Dmumu result");
  }
  return Pass(
      "srcSEP3D and direct sep_common Jokipii calls are bitwise identical at matched local state");
}

Result RunCOEF3D06() {
  // kappa(s)=7+3s is sampled at s0-h, s0, and s0+h.  The centred and both
  // one-sided boundary stencils must recover the same exact derivative.  The
  // final case proves that a missing two-sided neighbourhood is an error rather
  // than the former silent dKappa/ds=0 substitution.
  const double stepM = 4.0;
  const double center = 31.0;
  const double expected = 3.0;
  const double minus = center - expected * stepM;
  const double plus = center + expected * stepM;
  const struct Case {
    bool hasMinus;
    bool hasPlus;
  } cases[] = {{true, true}, {true, false}, {false, true}};

  for (const Case& item : cases) {
    T::ParallelKappaGradientStencil stencil;
    stencil.centerKappaM2PerS = center;
    stencil.stepM = stepM;
    stencil.hasMinus = item.hasMinus;
    stencil.minusKappaM2PerS = minus;
    stencil.hasPlus = item.hasPlus;
    stencil.plusKappaM2PerS = plus;
    double derivative = 0.0;
    const SEP3D::Core::Status status =
        T::EvaluateParallelKappaGradient(stencil, &derivative);
    if (!status.ok() || !Relative(derivative, expected, 1.0e-15))
      return Fail("parallel-kappa gradient stencil lost a linear coefficient derivative");
  }

  T::ParallelKappaGradientStencil missing;
  missing.centerKappaM2PerS = center;
  missing.stepM = stepM;
  double ignored = 0.0;
  if (T::EvaluateParallelKappaGradient(missing, &ignored).ok())
    return Fail("parallel-kappa gradient accepted a stencil with no neighbours");

  Result result = Pass(
      "centred and one-sided stencils recover nonzero dKappa_parallel/ds and fail closed without neighbours");
  result.metrics.push_back(
      {"linear_gradient_error", 0.0, 1.0e-15, "<=", "m/s"});
  return result;
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterTurbulenceTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using R = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* group, const char* name,
                 const char* description,
                 SEP3D::Testing::TestCallback callback) {
    D descriptor;
    descriptor.id = id;
    descriptor.group = group;
    descriptor.name = name;
    descriptor.description = description;
    descriptor.initialization = I::None;
    descriptor.supportedBuildModes = "standalone-no-AMPS";
    descriptor.runtime = R::Routine;
    descriptor.seedPolicy = "deterministic-no-rng";
    descriptor.stateIsolation = "fresh provider and local state per test";
    descriptor.callback = std::move(callback);
    return descriptor;
  };
  return {
      make("TUR3D01", "TUR3D", "Prescribed normalization",
           "Integrate the normalized finite-band Kolmogorov spectrum.",
           RunTUR3D01),
      make("TUR3D02", "TUR3D", "AWSoM wave mapping",
           "Map SI wave energies and direction under both polarities.",
           RunTUR3D02),
      make("TUR3D03", "TUR3D", "Resonance bounds",
           "Apply reject/extension policy below, inside, and above the band.",
           RunTUR3D03),
      make("TUR3D04", "TUR3D", "Missing turbulence policy",
           "Require failure unless ballistic transport is explicit.",
           RunTUR3D04),
      make("COEF3D01", "COEF3D", "Coefficient conversions",
           "Round-trip Dmumu, mean free path, and parallel diffusion.",
           RunCOEF3D01),
      make("COEF3D02", "COEF3D", "Shared coefficient kernel",
           "Compare the 3-D bridge and direct sep_common kernel bitwise.",
           RunCOEF3D02),
      make("COEF3D06", "COEF3D", "Parallel coefficient gradient",
           "Recover field-aligned kappa gradients with centred and boundary stencils.",
           RunCOEF3D06),
  };
}
