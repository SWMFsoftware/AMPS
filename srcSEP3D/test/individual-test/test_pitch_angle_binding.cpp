// ============================================================================
// srcSEP3D binding of the pitch-angle diffusion library (PAD3D01-PAD3D03).
//
// The library itself is verified by src/models/sep_common/pitch_angle_diffusion
// /run_all_tests.py.  These tests check what srcSEP3D adds around it:
//
//   PAD3D01  the coefficient bridge (Turbulence::EvaluateLocalScattering with
//            pitch_angle_diffusion_model = pitch-angle-library) returns exactly
//            the library's D, dD/dmu and scattering step limit for the AMPS
//            species/momentum, with the configured solar origin subtracted
//            from the position (heliocentric lambda law);
//   PAD3D02  the focused-diffusion core's stochastic pitch kick
//            (Transport::AdvanceFocused) equals the library's Milstein and
//            Euler-Maruyama kicks (Eqs. 157, 159) for the same normal draw,
//            including mirror folding at mu = +-1;
//   PAD3D03  the production focused-transport composer (AdvanceFocused:
//            deterministic half-steps, library coefficient, Milstein kick):
//            (a) without scattering, the exact constant-focusing solution
//                mu(t) = tanh(atanh(mu0) + v t/(2L)) at second order;
//            (b) with PA-ISO scattering, the stationary mean
//                <mu> = coth(xi) - 1/xi (Eq. 176) within the roadmap budget
//                5e-3, using the predeclared criterion
//                |err| + z SE <= 5e-3 and z SE <= 5e-3/4 (z = 2.576).
// Expected values are exact formulas or the library's own documented output
// for the identical call (PAD3D01/02 test the *binding*, not the physics).
// ============================================================================

#include "sep3d_test_registry.h"

#include "focused_transport.h"
#include "keyed_random.h"
#include "pitch_angle_diffusion/pitch_angle_diffusion.h"
#include "pitch_angle_diffusion/pitch_angle_kernels.h"
#include "turbulence_models.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

using SEP3D::Testing::Result;
namespace PA = SEP::PitchAngleDiffusion;
namespace PK = SEP::PitchAngleDiffusion::Kernels;
namespace T = SEP3D::Transport;
namespace TM = SEP3D::Turbulence;
namespace RM = SEP3D::RuntimeModel;
namespace C = SEP3D::Core;

Result Pass(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Pass;
  result.message = message; return result;
}
Result Fail(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Fail;
  result.message = message; return result;
}

std::string Num(double value) {
  std::ostringstream out; out.precision(17); out << value; return out.str();
}

// A complete PA-EPS configuration (lambda = 0.3 AU (R/1 GV)^(1/3) (r/1 AU)^(2/3)).
PA::ModelConfiguration EpsConfiguration(PA::StepRule rule) {
  PA::ModelConfiguration c;
  c.model = PA::ModelId::EpsilonForm;
  c.shape.epsilon = 0.048;
  c.lambdaLaw.lambdaReferenceM = 0.3 * C::Const::AU;
  c.lambdaLaw.referenceRigidityV = 1.0e9;
  c.lambdaLaw.rigidityExponent = 1.0 / 3.0;
  c.lambdaLaw.referenceRadiusM = C::Const::AU;
  c.lambdaLaw.radialExponent = 2.0 / 3.0;
  c.step.rule = rule;
  c.step.rateFraction = 0.003;
  c.step.ell = 1.0 / 64.0;
  c.step.deltaMuMax = 0.05;
  return c;
}

SEP3D::Background::BackgroundSample Background() {
  SEP3D::Background::BackgroundSample sample;
  sample.status = C::Status::OK();
  sample.valid = true;
  sample.B = {5.0e-9, 0.0, 0.0};
  sample.absB = 5.0e-9;
  sample.bHat = {1.0, 0.0, 0.0};
  sample.alfvenSpeedMpS = 5.0e4;
  sample.generation = 7;
  return sample;
}

Result RunPAD3D01() {
  const double momentum = 1.0e-19;  // ~10 MeV proton, kg m/s
  const C::Vec3 origin(1.0e9, -2.0e9, 3.0e8);
  const C::Vec3 position = origin + C::Vec3(0.6 * C::Const::AU, 0.8 * C::Const::AU, 0.0);
  for (PA::StepRule rule : {PA::StepRule::RateFraction, PA::StepRule::VanDenBerg2026,
                            PA::StepRule::Sandroos2014}) {
    const PA::ModelConfiguration cfg = EpsConfiguration(rule);
    if (!PA::SetActiveConfiguration(cfg).ok())
      return Fail("library configuration was not installed");
    TM::CoefficientSelection selection;
    selection.pitchAngle = RM::PitchAngleDiffusionModel::PitchAngleLibrary;
    selection.requireSpatialDiffusion = false;
    selection.requireMeanFreePath = false;
    selection.requirePitchAngleDiffusion = true;
    selection.solarOriginM = origin;
    selection.deterministicPitchRatePerS = 0.37;  // host dmu/dt
    TM::TurbulenceSample turbulence;  // unused by the library path
    turbulence.status = C::Status::OK();
    turbulence.valid = true;
    for (double mu : {-0.8, -0.05, 0.0, 0.31, 1.0}) {
      const TM::LocalScatteringCoefficients bridged = TM::EvaluateLocalScattering(
          turbulence, Background(), position, 0, C::Const::m_p, C::Const::e,
          momentum, mu, selection);
      if (!bridged.status.ok())
        return Fail("bridge failed at mu=" + Num(mu) + ": " + bridged.status.message);
      PA::ParticleState particle{C::Const::m_p, C::Const::e, momentum};
      PA::LocalState local;
      local.heliocentricPositionM = {{0.6 * C::Const::AU, 0.8 * C::Const::AU, 0.0}};
      const PA::CoefficientResult direct = PA::Evaluate(mu, particle, local, cfg);
      const PA::StepLimitResult limit =
          PA::ScatteringStepLimit(mu, 0.37, particle, local, cfg);
      if (!direct.status.ok() || !direct.dDmuMuDmuPerS || !limit.status.ok())
        return Fail("direct library evaluation failed");
      // The heliocentric position differs from the mesh position only by the
      // origin, so the lambda law (r = 1 AU here) must agree to roundoff.
      const double tol = 1e-13;
      if (std::fabs(bridged.dMuMuPerS - direct.dMuMuPerS) > tol * direct.amplitudeD0PerS ||
          std::fabs(bridged.dDmuMuDmuPerS - *direct.dDmuMuDmuPerS) >
              tol * direct.amplitudeD0PerS ||
          std::fabs(bridged.pitchAngleStepLimitS - limit.maximumStepS) >
              tol * limit.maximumStepS)
        return Fail("bridge differs from the library at mu=" + Num(mu) + " rule=" +
                    PA::Name(rule) + ": D " + Num(bridged.dMuMuPerS) + " vs " +
                    Num(direct.dMuMuPerS) + ", limit " +
                    Num(bridged.pitchAngleStepLimitS) + " vs " + Num(limit.maximumStepS));
    }
  }
  // Without the library selector the old models impose no scattering limit.
  TM::CoefficientSelection constant;
  constant.pitchAngle = RM::PitchAngleDiffusionModel::Constant;
  constant.constantDmumuPerS = 0.1;
  constant.requireSpatialDiffusion = false;
  constant.requireMeanFreePath = false;
  TM::TurbulenceSample turbulence;
  turbulence.status = C::Status::OK();
  turbulence.valid = true;
  const TM::LocalScatteringCoefficients old = TM::EvaluateLocalScattering(
      turbulence, Background(), position, 0, C::Const::m_p, C::Const::e, momentum,
      0.2, constant);
  if (!old.status.ok() || std::isfinite(old.pitchAngleStepLimitS))
    return Fail("a non-library model acquired a scattering step limit");
  return Pass("bridge = library (D, dD/dmu, step limit) for 3 rules x 5 mu with a "
              "translated solar origin; old models unchanged");
}

Result RunPAD3D02() {
  // PA-EPS coefficients from the library at a few mu (D0 = 0.4 s^-1), one
  // AdvanceFocused step without focusing/flow (its deterministic half-steps
  // are then identities) against the library kick with the same normal.
  PA::ShapeParameters shape;
  shape.epsilon = 0.048;
  const double massKg = C::Const::m_p, momentum = 1.0e-19;
  for (T::PitchAngleScheme scheme : {T::PitchAngleScheme::ReflectingMilstein,
                                     T::PitchAngleScheme::ReflectingEulerMaruyama}) {
    for (double mu : {-0.97, -0.2, 0.0, 0.45, 0.999, 1.0}) {
      for (double dt : {0.003, 0.4}) {  // 0.4 s forces large kicks and folding
        const PA::CoefficientResult c = PA::EvaluateShape(mu, 0.4, PA::ModelId::EpsilonForm, shape);
        T::FocusedParticleState state;
        state.positionM = {C::Const::AU, 0.0, 0.0};
        state.momentumKgMPerS = momentum;
        state.mu = mu;
        T::FocusedLocalState local;
        local.bHat = {1.0, 0.0, 0.0};
        local.dMuMuPerS = c.dMuMuPerS;
        local.dDmuMuDmuPerS = *c.dDmuMuDmuPerS;
        local.scheme = scheme;
        T::RandomKey key;
        key.campaignSeed = 11;
        key.particleId = 5;
        key.purpose = T::RandomPurpose::FocusedPitch;
        T::KeyedRandomStream stream(key), replay(key);
        T::FocusedRandomStreams streams;
        streams.pitch = &stream;
        const T::FocusedStepResult moved = T::AdvanceFocused(state, local, massKg, dt, streams);
        if (!moved.status.ok()) return Fail("AdvanceFocused failed: " + moved.status.message);
        const double zeta = replay.Normal01();
        const PK::KickResult kick = scheme == T::PitchAngleScheme::ReflectingMilstein
            ? PK::MilsteinKick(mu, c.dMuMuPerS, *c.dDmuMuDmuPerS, 0.0, dt, zeta)
            : PK::EulerMaruyamaKick(mu, c.dMuMuPerS, *c.dDmuMuDmuPerS, 0.0, dt, zeta);
        if (!kick.status.ok() || std::fabs(moved.state.mu - kick.mu) > 4e-16)
          return Fail("kick differs at mu=" + Num(mu) + " dt=" + Num(dt) + ": " +
                      Num(moved.state.mu) + " vs " + Num(kick.mu));
      }
    }
  }
  return Pass("AdvanceFocused pitch kick = library Milstein/Euler-Maruyama kick "
              "(same normal, folding included) at 6 mu x 2 dt");
}

Result RunPAD3D03() {
  // (a) Deterministic focusing: D = 0, divBhat = 1/L.  Exact solution of
  // dmu/dt = v (1 - mu^2)/(2L): mu(t) = tanh(atanh(mu0) + v t/(2L)).
  const double massKg = C::Const::m_p, momentum = 1.0e-19;
  const double v = T::RelativisticSpeed(momentum, massKg);
  const double inverseL = 1.0 / C::Const::AU;
  auto deterministic = [&](double dt, double mu0, double tEnd) {
    T::FocusedParticleState state;
    state.positionM = {C::Const::AU, 0.0, 0.0};
    state.momentumKgMPerS = momentum;
    state.mu = mu0;
    T::FocusedLocalState local;
    local.bHat = {1.0, 0.0, 0.0};
    local.divBhatPerM = inverseL;
    T::FocusedRandomStreams none;
    for (double t = 0.0; t < tEnd - 1e-9 * dt; t += dt)
      state = T::AdvanceFocused(state, local, massKg, dt, none).state;
    return state.mu;
  };
  const double tEnd = 2.0 * C::Const::AU / v;  // v t/(2L) = 1
  const double exact = std::tanh(std::atanh(-0.6) + 1.0);
  const double e1 = std::fabs(deterministic(tEnd / 50.0, -0.6, tEnd) - exact);
  const double e2 = std::fabs(deterministic(tEnd / 100.0, -0.6, tEnd) - exact);
  const double order = std::log2(e1 / e2);
  if (!(order > 1.8) || !(e2 < 1e-4))
    return Fail("deterministic focusing is not second order: errors " + Num(e1) + ", " +
                Num(e2) + " (order " + Num(order) + ")");

  // (b) Focused transport with PA-ISO scattering through the library pointer:
  // nu = 2 D0 set by the lambda law (lambda = v m, so D0 = 1/2 s^-1), and
  // 1/L = xi/v so that xi = v/(2 L D0) = 1.5.  Uniform start, 4/nu burn-in,
  // time average over 21/nu, step nu dt = 0.005.  Batches of particles give
  // the standard error; N T = 2.1e6 particle-relaxation-times makes z SE
  // <= 5e-3/4 (sigma_mu ~ 0.5, integrated autocorrelation 1/nu).
  PA::ModelConfiguration iso;
  iso.model = PA::ModelId::Isotropic;
  iso.lambdaLaw.lambdaReferenceM = v * 1.0;  // [m]: D0 = v/(2 lambda) = 0.5 s^-1
  iso.lambdaLaw.referenceRigidityV = 1.0;
  iso.lambdaLaw.rigidityExponent = 0.0;
  iso.lambdaLaw.referenceRadiusM = 1.0;
  iso.lambdaLaw.radialExponent = 0.0;
  iso.step.rule = PA::StepRule::RateFraction;
  iso.step.rateFraction = 0.005;
  if (!PA::SetActiveConfiguration(iso).ok()) return Fail("PA-ISO not installed");
  const double xi = 1.5;
  const double dt = 0.005, burn = 4.0, window = 21.0;
  const int batches = 20, perBatch = 5000;
  PA::ParticleState particle{massKg, C::Const::e, momentum};
  PA::LocalState heliocentric;
  heliocentric.heliocentricPositionM = {{C::Const::AU, 0.0, 0.0}};
  std::vector<double> batchMeans;
  for (int b = 0; b < batches; ++b) {
    double sum = 0.0;
    long long samples = 0;
    for (int i = 0; i < perBatch; ++i) {
      T::RandomKey key;
      key.campaignSeed = 0x9ad3d03ull;
      key.particleId = static_cast<std::uint64_t>(b) * perBatch + i + 1;
      key.purpose = T::RandomPurpose::FocusedPitch;
      T::KeyedRandomStream initial(key);
      T::FocusedParticleState state;
      state.positionM = {C::Const::AU, 0.0, 0.0};
      state.momentumKgMPerS = momentum;
      state.mu = 2.0 * initial.UniformOpen01() - 1.0;
      T::FocusedLocalState local;
      local.bHat = {1.0, 0.0, 0.0};
      local.divBhatPerM = xi / v;
      local.scheme = T::PitchAngleScheme::ReflectingMilstein;
      const int steps = static_cast<int>((burn + window) / dt + 0.5);
      const int burnSteps = static_cast<int>(burn / dt + 0.5);
      for (int n = 0; n < steps; ++n) {
        const PA::CoefficientResult c =
            PA::ActivePitchAngleDiffusion(state.mu, particle, heliocentric);
        if (!c.status.ok() || !c.dDmuMuDmuPerS) return Fail("library evaluation failed");
        local.dMuMuPerS = c.dMuMuPerS;
        local.dDmuMuDmuPerS = *c.dDmuMuDmuPerS;
        key.substep = static_cast<std::uint64_t>(n) + 1;
        T::KeyedRandomStream stream(key);
        T::FocusedRandomStreams streams;
        streams.pitch = &stream;
        state = T::AdvanceFocused(state, local, massKg, dt, streams).state;
        if (n >= burnSteps) { sum += state.mu; ++samples; }
      }
    }
    batchMeans.push_back(sum / samples);
  }
  double mean = 0.0, var = 0.0;
  for (double m : batchMeans) mean += m;
  mean /= batches;
  for (double m : batchMeans) var += (m - mean) * (m - mean);
  const double se = std::sqrt(var / (batches - 1) / batches);
  const double reference = 1.0 / std::tanh(xi) - 1.0 / xi;  // 0.438124726315845
  const double z = 2.5758293035489;  // two-sided 99%, one comparison
  if (z * se > 0.25 * 5e-3)
    return Fail("design check: z*SE = " + Num(z * se) + " > 1.25e-3");
  const double used = std::fabs(mean - reference) + z * se;
  if (used > 5e-3)
    return Fail("stationary <mu> " + Num(mean) + " +- " + Num(se) + " vs " + Num(reference));
  return Pass("composer: deterministic focusing order " + Num(order) +
              "; stationary <mu> = " + Num(mean) + " +- " + Num(se) + " vs coth(1.5)-1/1.5 = " +
              Num(reference) + " (|err|+zSE = " + Num(used) + " <= 5e-3)");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterPitchAngleBindingTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* name, const char* description,
                 RC runtime, SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = "PAD3D";
    d.description = description;
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = runtime; d.seedPolicy = "deterministic keyed streams";
    d.stateIsolation = "installs the library's active model (serial)";
    d.callback = std::move(callback); return d;
  };
  return {
      make("PAD3D01", "Pitch-angle library coefficient bridge",
           "EvaluateLocalScattering returns the library D, dD/dmu and step limit.",
           RC::Routine, RunPAD3D01),
      make("PAD3D02", "Focused kick equals library kernel",
           "AdvanceFocused stochastic kick = library Milstein/EM kick, same normal.",
           RC::Routine, RunPAD3D02),
      make("PAD3D03", "Focused composer with library coefficient",
           "Exact deterministic focusing; stationary <mu> = coth(xi) - 1/xi.",
           RC::Extended, RunPAD3D03),
  };
}
