// ============================================================================
// Deterministic unit tests of the pitch-angle diffusion library (PAD01-PAD12).
//
// Test ownership (roadmap Section 26.1): every expected value is either a
// fixture value printed in PITCH_ANGLE_SCATTERING_MODEL.md Section 18 / 5.6
// (computed there at 40 digits by an independent route), an exact algebraic
// identity, or an independent computation in this file (tanh-sinh
// quadrature, a different algorithm from the library's Gauss-Kronrod; finite
// differences; polynomial reconstruction).  No expected value is produced by
// the function under test.
//
// Tolerances (Section 18; roadmap Section 26.5): closed forms 1e-12 relative,
// quadrature-class fixtures 1e-9, values printed to 7 digits 5e-7 relative,
// smooth quadrature round trips 1e-8.
//
// Output: one line "[ID] PASS|FAIL message" per test, read by
// run_all_tests.py; exit status 1 if any test fails.
// ============================================================================

#include "pitch_angle_diffusion.h"
#include "pitch_angle_kernels.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace PA = SEP::PitchAngleDiffusion;
namespace K = SEP::PitchAngleDiffusion::Kernels;

namespace {

// CODATA 2018 / exact SI constants used by the fixtures (Section 18).
constexpr double kC = 299792458.0;
constexpr double kE = 1.602176634e-19;
constexpr double kProtonMass = 1.67262192369e-27;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kAU = 149597870700.0;

struct Outcome {
  bool pass = true;
  std::ostringstream message;
  void Check(bool condition, const std::string& what) {
    if (!condition) {
      if (pass) message << what; else message << "; " << what;
      pass = false;
    }
  }
};

bool Near(double actual, double expected, double relative,
          double absolute = 0.0) {
  return std::fabs(actual - expected) <=
         absolute + relative * std::fabs(expected);
}

std::string Show(double value) {
  std::ostringstream s;
  s.precision(16);
  s << value;
  return s.str();
}

// Independent quadrature: tanh-sinh (double exponential) on [a, b].  It
// tolerates integrable endpoint singularities, which is why the library's
// substitutions are not needed here.
double TanhSinh(const std::function<double(double)>& f, double a, double b) {
  const double half = 0.5 * (b - a);
  const double pi2 = 0.5 * 3.14159265358979323846;
  double h = 1.0, previous = 0.0, sum = 0.0;
  for (int level = 0; level < 12; ++level) {
    sum = 0.0;
    for (int k = -static_cast<int>(6.0 / h); k <= static_cast<int>(6.0 / h); ++k) {
      const double t = k * h;
      const double u = pi2 * std::sinh(t);
      const double x = std::tanh(u);
      const double w = pi2 * std::cosh(t) / (std::cosh(u) * std::cosh(u));
      // Distance to the endpoint computed without cancellation.
      const double complement =
          1.0 / (std::exp(std::fabs(u)) * std::cosh(u));  // 1 - |x|
      if (complement <= 0.0 || w < 1e-300) continue;
      const double point = x >= 0.0 ? b - half * complement
                                    : a + half * complement;
      if (!(point > a && point < b)) continue;
      sum += w * f(point);
    }
    sum *= h * half;
    if (level > 3 && std::fabs(sum - previous) <= 1e-14 * std::fabs(sum)) break;
    previous = sum;
    h *= 0.5;
  }
  return sum;
}

// A complete λ-law configuration for `model`.
PA::ModelConfiguration Config(PA::ModelId model) {
  PA::ModelConfiguration c;
  c.model = model;
  c.shape.q = 5.0 / 3.0;
  c.shape.gapH = 0.05;
  c.shape.epsilon = 0.048;
  c.lambdaLaw.lambdaReferenceM = 0.3 * kAU;
  c.lambdaLaw.referenceRigidityV = 1.0e9;
  c.lambdaLaw.rigidityExponent = 1.0 / 3.0;
  c.lambdaLaw.referenceRadiusM = kAU;
  c.lambdaLaw.radialExponent = 2.0 / 3.0;
  c.step.rule = PA::StepRule::RateFraction;
  c.step.rateFraction = 0.001;
  return c;
}

PA::ParticleState ProtonWithKineticEnergyMeV(double megaElectronVolt) {
  const double t = megaElectronVolt * 1.0e6 * kE;
  const double mc2 = kProtonMass * kC * kC;
  PA::ParticleState p;
  p.massKg = kProtonMass;
  p.chargeC = kE;
  p.momentumKgMPerS = std::sqrt(t * (t + 2.0 * mc2)) / kC;
  return p;
}

// ---------------------------------------------------------------------------
Outcome PAD01() {
  // Kinematics, Section 1.6; F-PUSH-19 (100 MeV proton rigidity).
  Outcome o;
  PA::Kinematics k;
  const PA::ParticleState proton = ProtonWithKineticEnergyMeV(100.0);
  o.Check(PA::ComputeKinematics(proton, &k).ok(), "kinematics failed");
  o.Check(Near(k.rigidityV, 444583420.329749, 1e-12),
          "rigidity " + Show(k.rigidityV) + " != 444583420.329749 V");
  o.Check(Near(k.kineticEnergyJ, 100.0e6 * kE, 1e-13), "kinetic energy round trip");
  const double gamma = 1.0 + 100.0e6 * kE / (kProtonMass * kC * kC);
  o.Check(Near(k.gamma, gamma, 1e-14), "gamma");
  o.Check(Near(k.speedMPerS, kC * std::sqrt(1.0 - 1.0 / (gamma * gamma)), 1e-13),
          "speed");
  // Low-momentum stability: T = p^2/(2m) (1 - p^2/(4 m^2 c^2) + ...).
  PA::ParticleState slow = proton;
  slow.momentumKgMPerS = 1.0e-7 * kProtonMass * kC;
  o.Check(PA::ComputeKinematics(slow, &k).ok(), "slow kinematics failed");
  const double p = slow.momentumKgMPerS;
  const double tSeries = p * p / (2.0 * kProtonMass) *
                         (1.0 - p * p / (4.0 * kProtonMass * kProtonMass * kC * kC));
  o.Check(Near(k.kineticEnergyJ, tSeries, 1e-14), "slow kinetic energy cancels");
  // Signed charge: rigidity uses |q| (an electron of the same p has the same R).
  PA::ParticleState electron;
  electron.massKg = kElectronMass;
  electron.chargeC = -kE;
  electron.momentumKgMPerS = proton.momentumKgMPerS;
  o.Check(PA::ComputeKinematics(electron, &k).ok() &&
              Near(k.rigidityV, 444583420.329749, 1e-12),
          "electron rigidity");
  // Refusals: p = 0, q = 0, m <= 0.
  PA::ParticleState bad = proton;
  bad.momentumKgMPerS = 0.0;
  o.Check(PA::ComputeKinematics(bad, &k).code == PA::StatusCode::InvalidInput,
          "p = 0 accepted");
  bad = proton; bad.chargeC = 0.0;
  o.Check(PA::ComputeKinematics(bad, &k).code == PA::StatusCode::InvalidInput,
          "q = 0 accepted");
  bad = proton; bad.massKg = -1.0;
  o.Check(PA::ComputeKinematics(bad, &k).code == PA::StatusCode::InvalidInput,
          "m < 0 accepted");
  if (o.pass) o.message << "F-PUSH-19 rigidity, stable T, signed charge, refusals";
  return o;
}

Outcome PAD02() {
  // Normalization constants, Section 5.6, F-EMP-01/02/03/04/08, F-GCR-PA-04.
  Outcome o;
  struct QCase { double q, h, expected, tolerance; const char* source; };
  const QCase qcases[] = {
      {5.0 / 3.0, 0.0, 1.92857142857143, 1e-12, "F-EMP-01"},
      {5.0 / 3.0, 0.05, 1.26978038591672, 1e-12, "F-EMP-01"},
      {5.0 / 3.0, 0.2, 0.839883195924871, 1e-12, "F-EMP-01"},
      {1.67, 0.05, 1.27586869391809, 1e-12, "F-EMP-01"},
      {1.41, 0.0, 0.981611151102677, 1e-12, "F-EMP-01"},
      {5.0 / 3.0, 0.01, 1.601995, 5e-7, "Sec. 5.6"},
      {5.0 / 3.0, 0.1, 1.067701, 5e-7, "Sec. 5.6"},
      {5.0 / 3.0, 0.5, 0.5374346, 5e-7, "Sec. 5.6"},
      {1.67, 0.5, 0.5383346, 5e-7, "Sec. 5.6"},
      {1.67, 0.0, 1.950839, 5e-7, "Sec. 5.6"},
      {1.66, 0.0, 1.885370, 5e-7, "Sec. 5.6"},
      {1.63, 0.0, 1.710571, 5e-7, "Sec. 5.6"},
      {1.5, 0.0, 1.2, 1e-12, "Sec. 5.6"},
  };
  for (const QCase& c : qcases) {
    double closed = 0.0, quadrature = 0.0;
    std::string method;
    o.Check(PA::QFormIntegralClosedForm(c.q, c.h, &closed, &method).ok(),
            "closed form failed");
    o.Check(PA::QFormIntegralQuadrature(c.q, c.h, &quadrature).ok(),
            "quadrature failed");
    o.Check(Near(0.75 * closed, c.expected, c.tolerance),
            std::string(c.source) + " q=" + Show(c.q) + " H=" + Show(c.h) +
                " closed D0 lambda/v " + Show(0.75 * closed));
    o.Check(Near(0.75 * quadrature, c.expected, std::max(c.tolerance, 1e-12)),
            std::string(c.source) + " quadrature " + Show(0.75 * quadrature));
  }
  // F-EMP-02 and F-GCR-PA-04 (I itself).
  double v = 0.0;
  std::string method;
  PA::QFormIntegralClosedForm(5.0 / 3.0, 0.05, &v, &method);
  o.Check(Near(v, 1.69304051455562, 1e-12), "I(5/3,0.05)");
  PA::QFormIntegralClosedForm(5.0 / 3.0, 0.2, &v, &method);
  o.Check(Near(v, 1.11984426123316, 1e-12), "I(5/3,0.2)");
  PA::QFormIntegralClosedForm(5.0 / 3.0, 0.01, &v, &method);
  o.Check(Near(v, 2.13599282779915, 1e-12), "I(5/3,0.01) (F-GCR-PA-04)");
  PA::QFormIntegralClosedForm(5.0 / 3.0, 1e-12, &v, &method);
  o.Check(Near((v - 18.0 / 7.0) / 1e-6, -4.71238538038469, 1e-6),
          "small-H slope (F-EMP-02) " + Show((v - 18.0 / 7.0) / 1e-6));
  // Series path away from q = 5/3 at large and moderate H agrees with
  // quadrature (Eq. 37 second form, rigorous tail bound).
  for (double h : {0.003, 0.3, 4.0, 50.0}) {
    double closed = 0.0, quadrature = 0.0;
    PA::QFormIntegralClosedForm(1.55, h, &closed, &method);
    PA::QFormIntegralQuadrature(1.55, h, &quadrature);
    o.Check(Near(closed, quadrature, 1e-12),
            "series vs quadrature at q=1.55 H=" + Show(h));
  }
  // F-EMP-03, F-EMP-04, D_A kappa/v^2 = phi/6.
  struct ECase { double eps, phi; };
  const ECase ecases[] = {{0.01, 7.06953301508509}, {0.045, 4.69919053319457},
                          {0.048, 4.597261773672}, {0.1, 3.45695783731876}};
  for (const ECase& c : ecases) {
    double closed = 0.0, quadrature = 0.0;
    PA::EpsilonPhiClosedForm(c.eps, &closed);
    PA::EpsilonPhiQuadrature(c.eps, &quadrature);
    o.Check(Near(closed, c.phi, 1e-12), "phi(" + Show(c.eps) + ") " + Show(closed));
    o.Check(Near(quadrature, c.phi, 1e-12), "phi quadrature " + Show(quadrature));
  }
  PA::EpsilonPhiClosedForm(0.045, &v);
  o.Check(Near(v / 6.0, 0.783198422199096, 1e-12), "D_A kappa/v^2 (F-EMP-03)");
  PA::EpsilonPhiClosedForm(0.01, &v);
  o.Check(Near(v * 0.51, 3.6054618376934, 1e-12), "nu(1) lambda/v (F-EMP-04)");
  o.Check(PA::EpsilonPhiClosedForm(0.0, &v).code == PA::StatusCode::OutsideDomain,
          "phi(0) not refused");
  // Shape constants of each model (exact rationals and F-EMP-08 factor 2).
  PA::ShapeConstants k;
  o.Check(PA::ValidateConfiguration(Config(PA::ModelId::Isotropic), &k).ok() &&
              k.normalization == 0.5 && k.lambdaEquation23OverLaw == 1.0,
          "PA-ISO N = 1/2");
  o.Check(PA::ValidateConfiguration(Config(PA::ModelId::Eprem), &k).ok() &&
              k.normalization == 0.25 && k.lambdaEquation23OverLaw == 2.0,
          "PA-EPREM N = 1/4, lambda factor 2 (F-EMP-08)");
  if (o.pass) o.message << "F-EMP-01/02/03/04/08, F-GCR-PA-04, Sec. 5.6 table, "
                           "closed form = quadrature";
  return o;
}

Outcome PAD03() {
  // Shape values and derivatives: F-PUSH-01/02/03/11/12, F-EMP-09/12, and an
  // independent central-difference check of every analytic derivative.
  Outcome o;
  PA::ShapeParameters q;
  q.q = 5.0 / 3.0;
  q.gapH = 0.05;
  PA::ShapeParameters q0 = q;
  q0.gapH = 0.0;
  PA::ShapeParameters e;
  e.epsilon = 0.048;
  auto dprime = [](double mu, PA::ModelId m, const PA::ShapeParameters& s) {
    const PA::CoefficientResult r = PA::EvaluateShape(mu, 1.0, m, s);
    return r.dDmuMuDmuPerS ? *r.dDmuMuDmuPerS
                           : std::numeric_limits<double>::quiet_NaN();
  };
  auto dvalue = [](double mu, PA::ModelId m, const PA::ShapeParameters& s) {
    return PA::EvaluateShape(mu, 1.0, m, s).dMuMuPerS;
  };
  o.Check(Near(dprime(0.1, PA::ModelId::QForm, q), 1.36883820162041, 1e-12),
          "F-PUSH-01 D'(0.1)");
  o.Check(Near(dprime(0.01, PA::ModelId::QForm, q0), 3.09315479871956, 1e-12),
          "F-PUSH-01 D'(0.01), H=0");
  o.Check(Near(dprime(0.01, PA::ModelId::QForm, q), 3.09215479871956, 1e-12),
          "Sec. 16.5 D'(0.01), H=0.05");
  o.Check(Near(dprime(1e-12, PA::ModelId::EpsilonForm, e), 0.999999999997904,
               1e-12), "F-PUSH-02 D'(1e-12)");
  o.Check(Near(dprime(0.3, PA::ModelId::EpsilonForm, e), 0.3712, 1e-12),
          "F-PUSH-02 D'(0.3)");
  // F-EMP-12: one-sided limits at 0 equal -nu0/2 and +nu0/2 (D0 = nu0/2).
  const PA::CoefficientResult at0 =
      PA::EvaluateShape(0.0, 0.5, PA::ModelId::EpsilonForm, e);
  o.Check(at0.derivativeMinusPerS && at0.derivativePlusPerS &&
              *at0.derivativeMinusPerS == -0.5 && *at0.derivativePlusPerS == 0.5,
          "F-EMP-12 one-sided limits");
  o.Check(at0.zero == PA::ZeroClass::DriftJump, "PA-EPS zero class");
  // PA-QFORM at mu = 0: no finite derivative (never zero).
  const PA::CoefficientResult q00 = PA::EvaluateShape(0.0, 1.0, PA::ModelId::QForm, q);
  o.Check(!q00.dDmuMuDmuPerS && !q00.derivativeMinusPerS &&
              q00.zero == PA::ZeroClass::UnboundedDrift,
          "PA-QFORM drift at 0 not reported as absent/unbounded");
  // F-PUSH-11: Kruells-Achterberg scale 2 D/D'^2.
  auto ka = [&](double mu, PA::ModelId m, const PA::ShapeParameters& s) {
    const double d = dprime(mu, m, s);
    return 2.0 * dvalue(mu, m, s) / (d * d);
  };
  o.Check(Near(ka(0.001, PA::ModelId::QForm, q), 0.00270009990282157, 1e-12),
          "F-PUSH-11 H=0.05");
  o.Check(Near(ka(0.001, PA::ModelId::QForm, q0), 0.000450003150018, 1e-12),
          "F-PUSH-11 H=0");
  o.Check(Near(ka(0.001, PA::ModelId::EpsilonForm, e), 0.0984100064104579, 1e-12),
          "F-PUSH-11 eps");
  // F-PUSH-03: b b' with b = sqrt(2D) equals D' (independent finite
  // difference of b at mu = 0.3, q = 5/3, H = 0.05).
  const double hstep = 1e-5;
  auto b = [&](double mu) { return std::sqrt(2.0 * dvalue(mu, PA::ModelId::QForm, q)); };
  const double bprime = (b(0.3 + hstep) - b(0.3 - hstep)) / (2.0 * hstep);
  o.Check(Near(b(0.3) * bprime, 0.607355341732575, 1e-8), "F-PUSH-03 b b'");
  o.Check(Near(dprime(0.3, PA::ModelId::QForm, q), 0.607355341732575, 1e-12),
          "F-PUSH-03 D'");
  // F-EMP-09 / F-PUSH-12: eps = 0 hemispheric map D_eta_eta = 4 D_mumu with
  // nu0 = 1 (D0 = 1/2).
  PA::ShapeParameters e0;
  e0.epsilon = 0.0;
  o.Check(Near(4.0 * PA::EvaluateShape(0.2, 0.5, PA::ModelId::EpsilonForm, e0).dMuMuPerS,
               0.32, 1e-14), "F-EMP-09 mu=0.2");
  o.Check(Near(4.0 * PA::EvaluateShape(-0.75, 0.5, PA::ModelId::EpsilonForm, e0).dMuMuPerS,
               0.375, 1e-14), "F-EMP-09 mu=-0.75");
  o.Check(Near(4.0 * PA::EvaluateShape(0.7, 0.5, PA::ModelId::EpsilonForm, e0).dMuMuPerS,
               0.42, 1e-14), "F-PUSH-12 mu=0.7");
  o.Check(PA::EvaluateShape(0.1, 0.5, PA::ModelId::EpsilonForm, e0).zero ==
              PA::ZeroClass::InternalBarrier, "eps=0 zero class");
  // Independent central differences (Richardson) at smooth points.
  struct M { PA::ModelId model; PA::ShapeParameters shape; };
  const M models[] = {{PA::ModelId::Isotropic, {}}, {PA::ModelId::QForm, q},
                      {PA::ModelId::QForm, q0}, {PA::ModelId::EpsilonForm, e},
                      {PA::ModelId::Eprem, {}}};
  for (const M& m : models)
    for (double mu : {-0.93, -0.41, -0.07, 0.03, 0.27, 0.66, 0.995}) {
      auto D = [&](double x) { return dvalue(x, m.model, m.shape); };
      const double h1 = 1e-4, h2 = 5e-5;
      const double c1 = (D(mu + h1) - D(mu - h1)) / (2 * h1);
      const double c2 = (D(mu + h2) - D(mu - h2)) / (2 * h2);
      const double richardson = (4.0 * c2 - c1) / 3.0;
      o.Check(Near(dprime(mu, m.model, m.shape), richardson, 1e-7, 1e-9),
              std::string(PA::ModelName(m.model)) + " D' vs finite difference at " +
                  Show(mu));
    }
  // Endpoint values: D(+-1) = 0 for every shape (zero flux, Section 3.1).
  for (const M& m : models) {
    o.Check(dvalue(1.0, m.model, m.shape) == 0.0 &&
                dvalue(-1.0, m.model, m.shape) == 0.0,
            "D(+-1) != 0");
  }
  if (o.pass) o.message << "F-PUSH-01/02/03/11/12, F-EMP-09/12, finite differences, "
                           "D(+-1)=0";
  return o;
}

Outcome PAD04() {
  // Amplitude -> Eq. (23) -> lambda round trip with an independent
  // (tanh-sinh) quadrature of Int (1-mu^2)^2/D, through the active pointer;
  // and the lambda law itself.
  Outcome o;
  const PA::ParticleState proton = ProtonWithKineticEnergyMeV(10.0);
  PA::LocalState local;
  local.heliocentricPositionM = {{0.3 * kAU, 0.4 * kAU, 0.0}};  // r = 0.5 AU
  PA::Kinematics k;
  PA::ComputeKinematics(proton, &k);
  struct Case { PA::ModelId model; double q, h, eps; };
  const Case cases[] = {
      {PA::ModelId::Isotropic, 0, 0, 0}, {PA::ModelId::Eprem, 0, 0, 0},
      {PA::ModelId::QForm, 5.0 / 3.0, 0.05, 0}, {PA::ModelId::QForm, 5.0 / 3.0, 0.0, 0},
      {PA::ModelId::QForm, 1.41, 0.2, 0}, {PA::ModelId::EpsilonForm, 0, 0, 0.048},
      {PA::ModelId::EpsilonForm, 0, 0, 0.01}, {PA::ModelId::EpsilonForm, 0, 0, 1.0}};
  for (const Case& c : cases) {
    PA::ModelConfiguration cfg = Config(c.model);
    cfg.shape.q = c.q > 0 ? c.q : cfg.shape.q;
    cfg.shape.gapH = c.model == PA::ModelId::QForm ? c.h : cfg.shape.gapH;
    cfg.shape.epsilon = c.eps > 0 ? c.eps : cfg.shape.epsilon;
    if (!PA::SetActiveConfiguration(cfg).ok()) {
      o.Check(false, std::string("cannot activate ") + PA::ModelName(c.model));
      continue;
    }
    const PA::CoefficientResult r =
        PA::ActivePitchAngleDiffusion(0.3, proton, local);
    o.Check(r.status.ok(), "evaluation failed");
    // Independent lambda law.
    const double lawExpected = cfg.lambdaLaw.lambdaReferenceM *
        std::pow(k.rigidityV / cfg.lambdaLaw.referenceRigidityV,
                 cfg.lambdaLaw.rigidityExponent) *
        std::pow(0.5, cfg.lambdaLaw.radialExponent);
    o.Check(Near(r.lambdaLawM, lawExpected, 1e-13), "lambda law value");
    // Eq. (23) by tanh-sinh, split at 0 (integrable |mu|^(1-q) at H = 0).
    auto integrand = [&](double mu) {
      const PA::CoefficientResult x = PA::ActivePitchAngleDiffusion(mu, proton, local);
      return (1.0 - mu * mu) * (1.0 - mu * mu) / x.dMuMuPerS;
    };
    const double integral = TanhSinh(integrand, -1.0, 0.0) +
                            TanhSinh(integrand, 0.0, 1.0);
    const double lambda23 = 3.0 * k.speedMPerS / 8.0 * integral;
    o.Check(Near(lambda23, r.lambdaParallelUnfocusedM, 1e-8),
            std::string(PA::ModelName(c.model)) + " Eq.(23) " + Show(lambda23) +
                " vs " + Show(r.lambdaParallelUnfocusedM));
    o.Check(r.lambdaStatus == PA::LambdaStatus::Finite, "lambda status");
    // Provenance carries the fingerprint and the specification digest.
    o.Check(std::string(r.provenance.configurationFingerprint) ==
                    PA::ConfigurationFingerprint(cfg) &&
                std::string(r.provenance.specificationSha256) ==
                    PA::SpecificationSha256(),
            "provenance");
  }
  // Isotropic identities, Eq. (24): nu = 2 D0, lambda = v/nu.
  PA::SetActiveConfiguration(Config(PA::ModelId::Isotropic));
  const PA::CoefficientResult iso = PA::ActivePitchAngleDiffusion(0.0, proton, local);
  o.Check(Near(iso.localScatteringFrequencyPerS * iso.lambdaLawM, k.speedMPerS, 1e-14),
          "nu lambda = v");
  // r = 0 is refused (InvalidField); mu outside [-1, 1] is refused.
  PA::LocalState origin;
  o.Check(PA::ActivePitchAngleDiffusion(0.2, proton, origin).status.code ==
              PA::StatusCode::InvalidField, "r = 0 accepted");
  o.Check(PA::ActivePitchAngleDiffusion(1.5, proton, local).status.code ==
              PA::StatusCode::InvalidInput, "mu = 1.5 accepted");
  if (o.pass) o.message << "Eq.(23) round trip (tanh-sinh) for 8 shapes, lambda law, "
                           "Eq.(24), refusals";
  return o;
}

PA::Status Parse(const std::vector<std::string>& body, PA::ParsedSection* out) {
  std::vector<PA::SectionLine> lines;
  std::size_t n = 10;
  for (const std::string& text : body) lines.push_back({n++, text});
  return PA::ParseSection(lines, out);
}

Outcome PAD05() {
  // Section parser and parser-neutral bridge.
  Outcome o;
  const std::vector<std::string> base = {
      "lambda_reference_m = 4.48793612e10   ! 0.3 AU",
      "lambda_reference_rigidity_v = 1e9",
      "lambda_rigidity_exponent = 1/3",
      "lambda_reference_radius_m = 1.495978707e11",
      "lambda_radial_exponent = 2/3",
      "step_rule = rate-fraction",
      "rate_fraction = 0.001"};
  auto with = [&](std::vector<std::string> extra) {
    std::vector<std::string> body = base;
    body.insert(body.end(), extra.begin(), extra.end());
    return body;
  };
  PA::ParsedSection parsed;
  PA::Status s = Parse(with({"model = pa-qform", "q = 5/3", "gap_h = 0.05"}), &parsed);
  o.Check(s.ok(), "valid PA-QFORM section rejected: " + s.detail);
  o.Check(parsed.configuration.model == PA::ModelId::QForm &&
              parsed.configuration.shape.q == 5.0 / 3.0 &&
              parsed.configuration.lambdaLaw.rigidityExponent == 1.0 / 3.0,
          "rational values not parsed exactly");
  o.Check(parsed.configurationFingerprint.size() == 64, "fingerprint length");
  s = Parse(with({"model = PA-EPS", "epsilon = 0.048"}), &parsed);
  o.Check(s.ok(), "PA-EPS rejected: " + s.detail);
  s = Parse(with({"model = PA-ISO"}), &parsed);
  o.Check(s.ok(), "PA-ISO rejected: " + s.detail);
  s = Parse(with({"model = pa-eprem"}), &parsed);
  o.Check(s.ok(), "PA-EPREM rejected: " + s.detail);
  // Refusals with their typed codes and line attribution.
  s = Parse(with({"model = QLT-SLAB"}), &parsed);
  o.Check(s.code == PA::StatusCode::UnsupportedModel &&
              s.detail.find("6.1") != std::string::npos &&
              s.detail.find("line 17") != std::string::npos,
          "QLT-SLAB: " + s.detail);
  s = Parse(with({"model = NL-NLD-SH08"}), &parsed);
  o.Check(s.code == PA::StatusCode::UnsupportedModel &&
              s.detail.find("DYN-SLS08") != std::string::npos,
          "alias not resolved: " + s.detail);
  s = Parse(with({"model = PA-FOO"}), &parsed);
  o.Check(s.code == PA::StatusCode::InvalidConfiguration, "unknown ID: " + s.detail);
  s = Parse(with({"model = PA-ISO", "gap_h = 0.05"}), &parsed);
  o.Check(s.code == PA::StatusCode::InvalidConfiguration &&
              s.detail.find("line 18") != std::string::npos,
          "inactive key: " + s.detail);
  s = Parse(with({"model = PA-QFORM", "q = 5/3"}), &parsed);
  o.Check(s.code == PA::StatusCode::MissingInput, "missing gap_h: " + s.detail);
  s = Parse(with({"model = PA-QFORM", "q = 2", "gap_h = 0"}), &parsed);
  o.Check(s.code == PA::StatusCode::OutsideDomain && s.detail.find("line 18") !=
              std::string::npos, "q = 2: " + s.detail);
  s = Parse(with({"model = PA-EPS", "epsilon = 0"}), &parsed);
  o.Check(s.code == PA::StatusCode::OutsideDomain &&
              s.detail.find("divergent") != std::string::npos,
          "eps = 0: " + s.detail);
  s = Parse(with({"model = PA-ISO", "model = PA-EPS"}), &parsed);
  o.Check(s.code == PA::StatusCode::InvalidConfiguration &&
              s.detail.find("duplicate") != std::string::npos, "duplicate: " + s.detail);
  s = Parse(with({"model = PA-ISO", "ell = 0.015625"}), &parsed);
  o.Check(s.code == PA::StatusCode::InvalidConfiguration, "ell with rate-fraction");
  s = Parse(base, &parsed);
  o.Check(s.code == PA::StatusCode::MissingInput, "missing model");
  s = Parse(with({"model PA-ISO"}), &parsed);
  o.Check(s.code == PA::StatusCode::InvalidConfiguration, "missing '='");
  // Rationals only for dimensionless keys.
  std::vector<std::string> bad = with({"model = PA-ISO"});
  bad[0] = "lambda_reference_m = 3/2";
  s = Parse(bad, &parsed);
  o.Check(s.code == PA::StatusCode::InvalidConfiguration, "rational length accepted");
  // Step rules and their parameter sets.
  std::vector<std::string> vdb = with({"model = PA-ISO", "ell = 0.015625"});
  vdb[5] = "step_rule = van-den-berg-2026";
  o.Check(Parse(vdb, &parsed).ok(), "van-den-berg-2026 rejected");
  std::vector<std::string> sand = with({"model = PA-ISO", "delta_mu_max = 0.05"});
  sand[5] = "step_rule = sandroos-2014";
  sand.erase(sand.begin() + 6);  // rate_fraction is not part of sandroos-2014
  o.Check(Parse(sand, &parsed).ok(), "sandroos-2014 rejected");
  sand.push_back("rate_fraction = 0.001");
  o.Check(Parse(sand, &parsed).code == PA::StatusCode::InvalidConfiguration,
          "rate_fraction accepted with sandroos-2014");
  // Fingerprint: identical input -> identical, changed value -> different.
  PA::ParsedSection a, b;
  Parse(with({"model = PA-EPS", "epsilon = 0.048"}), &a);
  Parse(with({"model = PA-EPS", "epsilon = 0.048"}), &b);
  o.Check(a.configurationFingerprint == b.configurationFingerprint, "fingerprint stable");
  Parse(with({"model = PA-EPS", "epsilon = 0.045"}), &b);
  o.Check(a.configurationFingerprint != b.configurationFingerprint, "fingerprint change");
  if (o.pass) o.message << "grammar, typed refusals, line attribution, rationals, "
                           "step-rule key sets, fingerprint";
  return o;
}

Outcome PAD06() {
  // Registry: the Section 4.2 inventory with the Section 4.1 counts.
  Outcome o;
  const auto& rows = PA::Registry();
  o.Check(rows.size() == 172, "registry size " + std::to_string(rows.size()));
  int counts[7] = {0, 0, 0, 0, 0, 0, 0};
  for (const auto& r : rows) ++counts[static_cast<int>(r.state)];
  o.Check(counts[0] == 86 && counts[1] == 6 && counts[2] == 2 && counts[3] == 31 &&
              counts[4] == 35 && counts[5] == 12 && counts[6] == 0,
          "per-state counts differ from Section 4.1");
  for (const char* id : {"PA-ISO", "PA-QFORM", "PA-EPS", "PA-EPREM"}) {
    const PA::RegistryEntry* e = PA::FindRegistryEntry(id);
    o.Check(e && e->state == PA::RuntimeState::ReadyExplicitInputs,
            std::string(id) + " not READY_EXPLICIT_INPUTS");
    PA::ModelId m;
    o.Check(PA::ParseModelId(id, &m).ok(), std::string(id) + " not implemented");
  }
  o.Check(PA::FindRegistryEntry("pa-amps-iii") != nullptr, "case-insensitive lookup");
  o.Check(std::string(PA::SpecificationSha256()).size() == 64, "spec digest");
  if (o.pass) o.message << "172 IDs, Section 4.1 counts, Release A IDs implemented";
  return o;
}

Outcome PAD07() {
  // Active configuration contract (transactional install, bound pointer).
  Outcome o;
  // A failed configuration leaves the active model untouched.
  PA::SetActiveConfiguration(Config(PA::ModelId::EpsilonForm));
  const auto before = PA::ActivePitchAngleDiffusion;
  PA::ModelConfiguration bad = Config(PA::ModelId::QForm);
  bad.shape.q = 2.5;
  o.Check(!PA::SetActiveConfiguration(bad).ok(), "invalid configuration installed");
  o.Check(PA::ActivePitchAngleDiffusion == before &&
              PA::GetActiveConfiguration().model == PA::ModelId::EpsilonForm,
          "failed install changed the active model");
  o.Check(PA::ActivePitchAngleDiffusion ==
              PA::BoundFunctionForModel(PA::ModelId::EpsilonForm), "bound pointer");
  o.Check(PA::HasActiveConfiguration(), "HasActiveConfiguration");
  if (o.pass) o.message << "transactional install, bound evaluator";
  return o;
}

Outcome PAD08() {
  // Coefficient-kernel readiness (Section 16.11).
  Outcome o;
  PA::ModelConfiguration q = Config(PA::ModelId::QForm);
  for (PA::Kernel kernel : {PA::Kernel::EulerMaruyama, PA::Kernel::Milstein}) {
    const PA::Status s = PA::CheckKernelPair(q, kernel);
    o.Check(s.code == PA::StatusCode::IncompatibleKernel &&
                s.detail.find("unbounded") != std::string::npos,
            "PA-QFORM pair not refused: " + s.detail);
  }
  // Every certified entry is accepted at its limit and refused just above.
  for (const PA::PairCertification& entry : PA::CertifiedPairs()) {
    PA::ModelConfiguration c = Config(entry.model);
    if (entry.model == PA::ModelId::EpsilonForm) c.shape.epsilon = entry.epsilonMinimum;
    c.step.rule = PA::StepRule::RateFraction;
    c.step.rateFraction = entry.maximumScaledStep;
    o.Check(PA::CheckKernelPair(c, entry.kernel).ok(),
            std::string("certified pair refused: ") + PA::ModelName(entry.model));
    c.step.rateFraction = entry.maximumScaledStep * 1.01;
    o.Check(PA::CheckKernelPair(c, entry.kernel).code ==
                PA::StatusCode::IncompatibleKernel, "step above limit accepted");
  }
  // A pair absent from the matrix is refused.
  PA::ModelConfiguration iso = Config(PA::ModelId::Isotropic);
  const PA::Status fv = PA::CheckKernelPair(iso, PA::Kernel::FiniteVolumeChain);
  o.Check(fv.code == PA::StatusCode::IncompatibleKernel,
          "finite-volume chain released as a pusher kernel");
  if (o.pass) o.message << "PA-QFORM refused with EM/Milstein; " <<
      PA::CertifiedPairs().size() << " certified pairs checked at their limits";
  return o;
}

Outcome PAD09() {
  // Step rules (Section 16.6; F-PUSH-19).
  Outcome o;
  // Choose lambda_ref so that D0 = 1 s^-1 for the chosen particle: then the
  // rule entries can be compared with the dimensionless F-PUSH-19 values.
  const PA::ParticleState proton = ProtonWithKineticEnergyMeV(100.0);
  PA::Kinematics k;
  PA::ComputeKinematics(proton, &k);
  PA::LocalState local;
  local.heliocentricPositionM = {{kAU, 0.0, 0.0}};
  PA::ModelConfiguration eps = Config(PA::ModelId::EpsilonForm);
  eps.shape.epsilon = 0.045;
  PA::ShapeConstants constants;
  PA::ValidateConfiguration(eps, &constants);
  eps.lambdaLaw.rigidityExponent = 0.0;
  eps.lambdaLaw.radialExponent = 0.0;
  eps.lambdaLaw.lambdaReferenceM = constants.normalization * k.speedMPerS;  // D0 = 1
  eps.step.rule = PA::StepRule::VanDenBerg2026;
  eps.step.ell = 1.0 / 64.0;
  eps.step.rateFraction = 1e-6;
  PA::StepLimitResult s = PA::ScatteringStepLimit(1.0, 0.0, proton, local, eps);
  o.Check(s.status.ok(), "van den Berg at mu = 1: " + s.status.detail);
  o.Check(s.driftEntryS == 0.0, "F-PUSH-19 drift entry at mu = 1 is not 0");
  o.Check(Near(s.noiseEntryS, 0.000223982224770642, 1e-12),
          "F-PUSH-19 noise entry " + Show(s.noiseEntryS));
  o.Check(s.limitingEntry == "rate-fraction-floor" &&
              Near(s.maximumStepS, s.rateFractionS, 1e-15),
          "floor not applied at mu = 1");
  const double nuMax = 2.0 * 1.0 * (0.5 + 0.045);
  o.Check(Near(s.rateFractionS * nuMax, 1e-6, 1e-12), "f/nu_max");
  // Interior point: min of the two entries when above the floor.
  s = PA::ScatteringStepLimit(0.3, 0.0, proton, local, eps);
  o.Check(s.status.ok() && Near(s.maximumStepS, std::min(s.driftEntryS, s.noiseEntryS),
                                1e-15), "van den Berg interior min");
  // The host's deterministic rate enters the drift entry.
  const PA::StepLimitResult withRate =
      PA::ScatteringStepLimit(0.3, 0.5, proton, local, eps);
  const PA::CoefficientResult c = PA::Evaluate(0.3, proton, local, eps);
  o.Check(Near(withRate.driftEntryS,
               std::sqrt(1 - 0.09) / 64.0 / std::fabs(*c.dDmuMuDmuPerS + 0.5), 1e-13),
          "drift entry with dmu/dt");
  // Rate fraction: nu_max dt = f exactly.
  PA::ModelConfiguration iso = Config(PA::ModelId::Isotropic);
  iso.step.rateFraction = 0.003;
  s = PA::ScatteringStepLimit(0.2, 0.0, proton, local, iso);
  const PA::CoefficientResult ci = PA::Evaluate(0.2, proton, local, iso);
  o.Check(s.status.ok() && Near(s.maximumStepS * ci.maximumScatteringFrequencyPerS,
                                0.003, 1e-14), "rate fraction");
  // Sandroos: isotropic a_max = 2 D0, b_max^2 = 2 D0 (independent algebra).
  iso.step.rule = PA::StepRule::Sandroos2014;
  iso.step.deltaMuMax = 0.05;
  s = PA::ScatteringStepLimit(0.2, 0.0, proton, local, iso);
  const double d0 = ci.amplitudeD0PerS;
  o.Check(s.status.ok() && Near(s.maximumStepS,
                                std::min(0.05 / (2 * d0), 0.0025 / (2 * d0)), 1e-14),
          "Sandroos isotropic");
  // Sandroos for PA-EPS: a_max = (1 + 2 eps) D0 (derived, Eq. 164).
  eps.step.rule = PA::StepRule::Sandroos2014;
  eps.step.deltaMuMax = 0.1;
  s = PA::ScatteringStepLimit(0.2, 0.0, proton, local, eps);
  o.Check(s.status.ok() && Near(s.driftEntryS, 0.1 / (1.0 + 2 * 0.045), 1e-12),
          "Sandroos eps drift entry " + Show(s.driftEntryS));
  // sup g for PA-EPS: independent dense scan of (1 - mu^2)(mu/(1+mu) + eps).
  double scan = 0.0;
  for (int i = 0; i <= 200000; ++i) {
    const double mu = i / 200000.0;
    scan = std::max(scan, (1 - mu * mu) * (mu / (1 + mu) + 0.045));
  }
  o.Check(Near(s.noiseEntryS, 0.01 / (2.0 * scan), 1e-9), "Sandroos eps noise entry");
  // PA-QFORM: no Ito step rule.
  PA::ModelConfiguration q = Config(PA::ModelId::QForm);
  o.Check(PA::ScatteringStepLimit(0.3, 0.0, proton, local, q).status.code ==
              PA::StatusCode::IncompatibleKernel, "PA-QFORM step rule");
  if (o.pass) o.message << "F-PUSH-19 entries, floor, dmu/dt, rate fraction, "
                           "Sandroos sup values";
  return o;
}

Outcome PAD10() {
  // Kernels: one-step moments.  An Ito kick (without folding) is the
  // polynomial a + b zeta + c (zeta^2 - 1) in the normal variate; it is
  // reconstructed from three evaluations, which gives the exact moments
  // E = a and E[dmu^2] = a^2 + b^2 + 2 c^2 (F-PUSH-04a: nu = 1, dt = 0.01,
  // mu = 0.3, isotropic D = (1/2)(1 - mu^2)).
  Outcome o;
  const double mu = 0.3, d = 0.5 * (1 - mu * mu), dp = -mu, dt = 0.01;
  auto moments = [&](bool milstein, double* mean, double* second) {
    auto kick = [&](double z) {
      return (milstein ? K::MilsteinKick(mu, d, dp, 0.0, dt, z)
                       : K::EulerMaruyamaKick(mu, d, dp, 0.0, dt, z)).mu - mu;
    };
    const double k0 = kick(0.0), kp = kick(1.0), km = kick(-1.0);
    const double b = 0.5 * (kp - km);
    const double c = 0.5 * (kp + km) - k0;
    const double a = k0 + c;
    *mean = a;
    *second = a * a + b * b + 2 * c * c;
  };
  double mean = 0.0, second = 0.0;
  moments(false, &mean, &second);
  o.Check(Near(mean, -0.003, 1e-12) && Near(second, 0.009109, 1e-12),
          "F-PUSH-04a Euler-Maruyama " + Show(mean) + " " + Show(second));
  moments(true, &mean, &second);
  o.Check(Near(mean, -0.003, 1e-12) && Near(second, 0.0091135, 1e-12),
          "F-PUSH-04a Milstein " + Show(second));
  // Random rotation: exact moments (Eq. 169) by quadrature over the
  // exponential theta^2 and the uniform azimuth (F-PUSH-05: nu dt = 0.01,
  // mu = 0.3).
  auto rotationMoment = [&](int power, double m0, double nudt) {
    // E over uTheta (tanh-sinh on (0,1]) and uPhi (periodic trapezoid).
    auto inner = [&](double uTheta) {
      const int nphi = 64;
      double sum = 0.0;
      for (int j = 0; j < nphi; ++j) {
        const double uPhi = (j + 0.5) / nphi;
        const double x = K::RandomRotationKick(m0, nudt, 1.0, uTheta, uPhi).mu - m0;
        sum += power == 1 ? x : x * x;
      }
      return sum / nphi;
    };
    return TanhSinh(inner, 0.0, 1.0);
  };
  o.Check(Near(rotationMoment(1, 0.3, 0.01), -0.00299001997146029, 1e-11),
          "F-PUSH-05 E[dmu] " + Show(rotationMoment(1, 0.3, 0.01)));
  o.Check(Near(rotationMoment(2, 0.3, 0.01), 0.00899745288638967, 1e-11),
          "F-PUSH-05 E[dmu^2] " + Show(rotationMoment(2, 0.3, 0.01)));
  // C1(m = 0.02) = <cos theta> (F-PUSH-05): mu = 1 gives mu' = cos theta.
  o.Check(Near(1.0 + rotationMoment(1, 1.0, 0.01), 0.990033266761799, 1e-11),
          "F-PUSH-05 C1(0.02)");
  // Refusals: zero step (Nonprogress), non-finite drift (IncompatibleKernel).
  o.Check(K::EulerMaruyamaKick(0.1, 0.1, 0.0, 0.0, 0.0, 0.3).status.code ==
              PA::StatusCode::Nonprogress, "zero-time kick accepted");
  o.Check(K::MilsteinKick(0.1, 0.1, std::numeric_limits<double>::infinity(), 0.0,
                          0.01, 0.3).status.code == PA::StatusCode::IncompatibleKernel,
          "infinite drift accepted");
  // The rotation keeps |mu| <= 1 without reflection for every variate.
  bool inside = true;
  for (int i = 1; i <= 50; ++i)
    for (int j = 0; j < 50; ++j) {
      const auto r = K::RandomRotationKick(0.99, 3.0, 1.0, i / 50.0, j / 50.0);
      inside = inside && std::fabs(r.mu) <= 1.0;
    }
  o.Check(inside, "rotation left [-1, 1]");
  if (o.pass) o.message << "F-PUSH-04a, F-PUSH-05 (Eq. 169), refusals";
  return o;
}

Outcome PAD11() {
  // Mirror folding (Eq. 161 applied until inside, Section 16.4).
  Outcome o;
  auto iterate = [](double x) {
    unsigned n = 0;
    while (x > 1.0 || x < -1.0) {
      x = x > 1.0 ? 2.0 - x : -2.0 - x;
      ++n;
    }
    return std::make_pair(x, n);
  };
  o.Check(Near(K::FoldPitchCosine(0.9 + 0.15), 0.95, 1e-15),
          "Merten & Aerdker example 0.9 + 0.15 -> 0.95");
  for (double x : {1.0, -1.0, 1.3, -1.7, 3.5, -4.25, 7.9, 12.01, -19.6}) {
    unsigned closed = 0;
    const double folded = K::FoldPitchCosine(x, &closed);
    const auto reference = iterate(x);
    o.Check(Near(folded, reference.first, 1e-14, 1e-14) && closed == reference.second,
            "fold " + Show(x));
  }
  if (o.pass) o.message << "closed-form folding equals repeated Eq. (161)";
  return o;
}

std::vector<double> UniformEdges(int n) {
  std::vector<double> e(n + 1);
  for (int i = 0; i <= n; ++i) e[i] = -1.0 + 2.0 * i / n;
  e.front() = -1.0;
  e.back() = 1.0;
  return e;
}

// Nonuniform grid clustered at mu = +-1 and 0 (exercises irregular widths).
std::vector<double> ClusteredEdges(int n) {
  std::vector<double> e(n + 1);
  for (int i = 0; i <= n; ++i) {
    const double s = -1.0 + 2.0 * i / n;
    e[i] = std::sin(0.5 * 3.14159265358979323846 * s) * (0.75 + 0.25 * s * s);
  }
  e.front() = -1.0;
  e.back() = 1.0;
  return e;
}

Outcome PAD12() {
  // Finite-volume reference (Eq. 178; AP-FV-01 deterministic part; PS08
  // reference part).
  Outcome o;
  PA::ShapeParameters q0;
  q0.q = 5.0 / 3.0;
  q0.gapH = 0.0;
  K::FvGenerator g;
  // F-PUSH-21: T between mu = -0.05 and 0.05, q = 5/3, H = 0, D0 = 1.
  o.Check(K::BuildFvGenerator(PA::ModelId::QForm, q0, 1.0, {-1.0, -0.1, 0.0, 0.1, 1.0},
                              &g).ok(), "FV build failed");
  o.Check(Near(g.conductance[1], 0.452241203587312, 1e-9),
          "F-PUSH-21 T " + Show(g.conductance[1]));
  // Generator invariants on an irregular grid: zero row sums, detailed
  // balance m_i Q_ij = m_j Q_ji, positive off-diagonals.
  PA::ShapeParameters e;
  e.epsilon = 0.048;
  const std::vector<double> clustered = ClusteredEdges(101);
  o.Check(K::BuildFvGenerator(PA::ModelId::EpsilonForm, e, 0.5, clustered, &g).ok(),
          "FV eps build");
  double worstRow = 0.0, worstBalance = 0.0;
  for (std::size_t i = 0; i < g.centres.size(); ++i) {
    const double scale = std::fabs(g.rateDiagonal[i]);
    worstRow = std::max(worstRow, std::fabs(g.rateLeft[i] + g.rateRight[i] +
                                            g.rateDiagonal[i]) / scale);
    if (i + 1 < g.centres.size())
      worstBalance = std::max(worstBalance,
          std::fabs(g.stationaryMass[i] * g.rateRight[i] -
                    g.stationaryMass[i + 1] * g.rateLeft[i + 1]) /
              (g.stationaryMass[i] * g.rateRight[i]));
  }
  o.Check(worstRow <= 1e-12 && worstBalance <= 1e-12,
          "row sum " + Show(worstRow) + ", balance " + Show(worstBalance));
  // eps = 0: disconnected hemispheres (Section 16.5).
  PA::ShapeParameters e0;
  e0.epsilon = 0.0;
  K::BuildFvGenerator(PA::ModelId::EpsilonForm, e0, 0.5, UniformEdges(100), &g);
  o.Check(g.internalExtension.find("disconnected") != std::string::npos &&
              g.conductance[49] == 0.0, "eps = 0 not disconnected");
  // Isotropic Legendre decay, Eq. (177) with constant nu = 1 (D0 = 1/2):
  // smooth initial density f0 = 1/2 + sum_l c_l P_l, exact
  // <P_l>(t) = <P_l>(0) exp(-l(l+1) t/2); three refinements of an irregular
  // grid must converge, and the finest must be within a quarter of the 5e-3
  // Monte Carlo budget (reference error, Section 26.5).
  const double coefficients[5] = {0.5, 0.3, 0.2, 0.1, 0.05};
  double previousError = 1.0;
  for (int n : {100, 200, 400}) {
    const std::vector<double> edges = ClusteredEdges(n);
    K::FvGenerator iso;
    K::BuildFvGenerator(PA::ModelId::Isotropic, {}, 0.5, edges, &iso);
    K::FvPropagator p;
    o.Check(K::PrepareFvPropagator(iso, &p).ok(), "propagator");
    std::vector<double> m0(n), m1;
    for (int i = 0; i < n; ++i) {
      double mass = 0.5 * (edges[i + 1] - edges[i]);
      for (int l = 1; l <= 4; ++l)
        mass += coefficients[l] * K::LegendreCellAverage(l, edges[i], edges[i + 1]) *
                (edges[i + 1] - edges[i]);
      m0[i] = mass;
    }
    double worst = 0.0;
    for (double t : {0.1, 1.0, 5.0}) {
      K::PropagateFv(p, m0, t, &m1);
      for (int l = 1; l <= 4; ++l) {
        double moment = 0.0;
        for (int i = 0; i < n; ++i)
          moment += m1[i] * K::LegendreCellAverage(l, edges[i], edges[i + 1]);
        // <P_l>(0) = c_l * 2/(2l+1) for f0 normalized to total mass 1.
        const double exact = coefficients[l] * 2.0 / (2.0 * l + 1.0) *
                             std::exp(-0.5 * l * (l + 1) * t);
        worst = std::max(worst, std::fabs(moment - exact));
      }
    }
    o.Check(worst < previousError, "FV Legendre error did not decrease at n=" +
                                       std::to_string(n));
    previousError = worst;
    if (n == 400)
      o.Check(worst <= 1.25e-3, "FV Legendre reference error " + Show(worst));
  }
  // Asymptotic variance, isotropic: <mu,(-G)^-1 mu> = 1/(3 nu) (Eq. 175).
  K::FvGenerator iso;
  K::BuildFvGenerator(PA::ModelId::Isotropic, {}, 0.5, ClusteredEdges(400), &iso);
  std::vector<double> muCell(iso.centres.size());
  for (std::size_t i = 0; i < muCell.size(); ++i)
    muCell[i] = 0.5 * (iso.edges[i] + iso.edges[i + 1]);
  double variance = 0.0;
  K::FvAsymptoticVariance(iso, muCell, &variance);
  o.Check(Near(variance, 1.0 / 3.0, 1e-4), "isotropic kappa/v^2 " + Show(variance));
  // Constant focusing with isotropic scattering (Section 16.12 second
  // reference): masses m_i = Int rho_+, T^rho = [Int (D rho_+)^-1]^-1 with
  // rho_+ ~ exp(xi mu), D = (1/2)(1 - mu^2) (v = nu = 1, xi = v/(2 L D0)).
  // F-FOC-01 / F-PUSH-08 stationary means and F-FOC-04 kappa_var.
  struct F { double xi, mean, kappa; };
  const F cases[] = {{1.0, 1.0 / std::tanh(1.0) - 1.0, 0.248076740609470},
                     {1.5, 0.438124726315845, 0.178275780337345},
                     {-1.5, -0.438124726315845, 0.178275780337345},
                     {5.0, 0.800090803982019, 0.0103248294848381}};
  for (const F& f : cases) {
    const int n = 1600;
    const std::vector<double> edges = ClusteredEdges(n);
    // Scaled weight exp(xi (mu - sgn(xi))) avoids overflow at large |xi|.
    const double shift = f.xi > 0 ? 1.0 : -1.0;
    auto rho = [&](double mu) { return std::exp(f.xi * (mu - shift)); };
    std::vector<double> masses(n), conductance(n - 1), meanMu(n);
    double total = 0.0, first = 0.0;
    for (int i = 0; i < n; ++i) {
      masses[i] = TanhSinh(rho, edges[i], edges[i + 1]);
      const double moment =
          TanhSinh([&](double mu) { return mu * rho(mu); }, edges[i], edges[i + 1]);
      meanMu[i] = moment / masses[i];
      total += masses[i];
      first += moment;
    }
    for (int i = 0; i + 1 < n; ++i) {
      const double a = 0.5 * (edges[i] + edges[i + 1]);
      const double b = 0.5 * (edges[i + 1] + edges[i + 2]);
      conductance[i] = 1.0 / TanhSinh(
          [&](double mu) { return 1.0 / (0.5 * (1.0 - mu * mu) * rho(mu)); }, a, b);
    }
    o.Check(Near(first / total, f.mean, 1e-12), "stationary mean xi=" + Show(f.xi));
    K::FvGenerator foc;
    o.Check(K::BuildFvGeneratorFromConductances(edges, masses, conductance, &foc).ok(),
            "focusing generator");
    double kappa = 0.0;
    K::FvAsymptoticVariance(foc, meanMu, &kappa);
    o.Check(Near(kappa, f.kappa, 1e-4), "F-FOC-04 kappa_var xi=" + Show(f.xi) + " " +
                                            Show(kappa));
  }
  if (o.pass) o.message << "F-PUSH-21, generator invariants, eps=0 barrier, Legendre "
                           "convergence, Eq.(175) limit, F-FOC-01/04, F-PUSH-08";
  return o;
}

}  // namespace

int main() {
  struct Test { const char* id; const char* title; Outcome (*run)(); };
  const Test tests[] = {
      {"PAD01", "kinematics", PAD01},
      {"PAD02", "lambda normalization constants", PAD02},
      {"PAD03", "shapes and derivatives", PAD03},
      {"PAD04", "Eq. (23) round trip and lambda law", PAD04},
      {"PAD05", "section parser", PAD05},
      {"PAD06", "model registry", PAD06},
      {"PAD07", "active configuration", PAD07},
      {"PAD08", "coefficient-kernel readiness", PAD08},
      {"PAD09", "scattering step rules", PAD09},
      {"PAD10", "kernel one-step moments", PAD10},
      {"PAD11", "mirror folding", PAD11},
      {"PAD12", "finite-volume reference", PAD12},
  };
  int failures = 0;
  for (const Test& test : tests) {
    const Outcome o = test.run();
    std::printf("[%s] %s %s: %s\n", test.id, o.pass ? "PASS" : "FAIL", test.title,
                o.message.str().c_str());
    if (!o.pass) ++failures;
  }
  std::printf("pitch-angle unit summary: %d passed, %d failed\n",
              static_cast<int>(sizeof(tests) / sizeof(tests[0])) - failures, failures);
  return failures ? 1 : 0;
}
