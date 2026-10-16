// ============================================================================
// Pitch-angle diffusion coefficient library, Release A (see the header for
// the scientific contract).  Organization of this file:
//
//   1. status/enum names and the generated model registry (Section 4.2);
//   2. kinematics (Section 1.6);
//   3. adaptive Gauss-Kronrod quadrature (numerical building block);
//   4. shape normalization integrals I(q,H) and phi(eps) (Section 5.6) by a
//      closed form and by an independent quadrature;
//   5. the shapes g(mu) = (1 - mu^2) h(mu), their derivatives and suprema;
//   6. configuration validation, evaluation and step limits;
//   7. (the coefficient-kernel certification matrix is in
//      pitch_angle_certification.cpp);
//   8. the parser-neutral bridge, the section parser and the fingerprint;
//   9. the active configuration and the mover-facing dispatch pointers.
//
// Every formula is quoted with its equation label in
// PITCH_ANGLE_SCATTERING_MODEL.md v1.2; "derived" marks a consequence
// computed for this implementation and verified in test_pitch_angle_*.cpp.
// ============================================================================

#include "pitch_angle_diffusion.h"
#include "pitch_angle_quadrature.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace SEP {
namespace PitchAngleDiffusion {

// ---------------------------------------------------------------------------
// 1. Names and registry.
// ---------------------------------------------------------------------------

Status Status::Error(StatusCode code, const std::string& detail) {
  Status status;
  status.code = code;
  status.detail = detail;
  return status;
}

const char* Name(StatusCode code) {
  switch (code) {
    case StatusCode::Success: return "Success";
    case StatusCode::InvalidInput: return "InvalidInput";
    case StatusCode::InvalidField: return "InvalidField";
    case StatusCode::OutsideDomain: return "OutsideDomain";
    case StatusCode::UnresolvedSource: return "UnresolvedSource";
    case StatusCode::Inadmissible: return "Inadmissible";
    case StatusCode::IncompatibleKernel: return "IncompatibleKernel";
    case StatusCode::Nonconverged: return "Nonconverged";
    case StatusCode::Nonprogress: return "Nonprogress";
    case StatusCode::InvalidConfiguration: return "InvalidConfiguration";
    case StatusCode::MissingInput: return "MissingInput";
    case StatusCode::UnsupportedModel: return "UnsupportedModel";
  }
  return "Unknown";
}

const char* Name(RuntimeState state) {
  switch (state) {
    case RuntimeState::ReadyExplicitInputs: return "READY_EXPLICIT_INPUTS";
    case RuntimeState::ReadyPublishedVariant: return "READY_PUBLISHED_VARIANT";
    case RuntimeState::RequiresExternalInput: return "REQUIRES_EXTERNAL_INPUT";
    case RuntimeState::RequiresUserDecision: return "REQUIRES_USER_DECISION";
    case RuntimeState::RequiresSourceOrCodeAudit:
      return "REQUIRES_SOURCE_OR_CODE_AUDIT";
    case RuntimeState::ReferenceDataOnly: return "REFERENCE_DATA_ONLY";
    case RuntimeState::BlockedDependency: return "BLOCKED_DEPENDENCY";
  }
  return "UNKNOWN";
}

const char* Name(OperatorConvention convention) {
  switch (convention) {
    case OperatorConvention::Standard: return "standard";
    case OperatorConvention::EpremHalfD: return "eprem-half-d";
  }
  return "unknown";
}

const char* Name(ZeroClass zero) {
  switch (zero) {
    case ZeroClass::Regular: return "regular";
    case ZeroClass::DriftJump: return "drift-jump";
    case ZeroClass::UnboundedDrift: return "unbounded-drift";
    case ZeroClass::InternalBarrier: return "internal-barrier";
  }
  return "unknown";
}

const char* Name(Kernel kernel) {
  switch (kernel) {
    case Kernel::EulerMaruyama: return "euler-maruyama";
    case Kernel::Milstein: return "milstein";
    case Kernel::RandomRotation: return "random-rotation";
    case Kernel::FiniteVolumeChain: return "finite-volume-chain";
  }
  return "unknown";
}

const char* Name(StepRule rule) {
  switch (rule) {
    case StepRule::RateFraction: return "rate-fraction";
    case StepRule::VanDenBerg2026: return "van-den-berg-2026";
    case StepRule::Sandroos2014: return "sandroos-2014";
  }
  return "unknown";
}

const char* Name(LambdaStatus status) {
  switch (status) {
    case LambdaStatus::Finite: return "Finite";
    case LambdaStatus::Divergent: return "Divergent";
    case LambdaStatus::Unavailable: return "Unavailable";
    case LambdaStatus::Nonconverged: return "Nonconverged";
  }
  return "Unknown";
}

namespace {

// The generated include supplies PITCH_ANGLE_SPECIFICATION_SHA256 and the 172
// rows {id, aliases, section, equations, state}.
const std::vector<RegistryEntry>& RegistryRows() {
  static const std::vector<RegistryEntry> rows = {
#include "pitch_angle_registry.inc"
  };
  return rows;
}

std::string Upper(std::string text) {
  for (char& c : text)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return text;
}

std::string Lower(std::string text) {
  for (char& c : text)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return text;
}

std::string Trim(const std::string& text) {
  const char* blanks = " \t\r\n\f\v";
  const std::size_t first = text.find_first_not_of(blanks);
  if (first == std::string::npos) return std::string();
  const std::size_t last = text.find_last_not_of(blanks);
  return text.substr(first, last - first + 1);
}

}  // namespace

const std::vector<RegistryEntry>& Registry() { return RegistryRows(); }

const char* SpecificationSha256() { return PITCH_ANGLE_SPECIFICATION_SHA256; }

const RegistryEntry* FindRegistryEntry(const std::string& id) {
  const std::string wanted = Upper(Trim(id));
  if (wanted.empty()) return nullptr;
  for (const RegistryEntry& entry : RegistryRows()) {
    if (Upper(entry.id) == wanted) return &entry;
    // Aliases are listed comma separated (Section 4.2 "alias" notes); a
    // merged ID is implemented once under its primary ID.
    std::stringstream aliases(entry.aliases);
    std::string alias;
    while (std::getline(aliases, alias, ','))
      if (!Trim(alias).empty() && Upper(Trim(alias)) == wanted) return &entry;
  }
  return nullptr;
}

const char* ModelName(ModelId model) {
  switch (model) {
    case ModelId::Isotropic: return "PA-ISO";
    case ModelId::QForm: return "PA-QFORM";
    case ModelId::EpsilonForm: return "PA-EPS";
    case ModelId::Eprem: return "PA-EPREM";
  }
  return "unknown";
}

Status ParseModelId(const std::string& text, ModelId* model) {
  if (model == nullptr)
    return Status::Error(StatusCode::InvalidConfiguration, "null model output");
  const RegistryEntry* entry = FindRegistryEntry(text);
  if (entry == nullptr)
    return Status::Error(StatusCode::InvalidConfiguration,
        "model '" + Trim(text) + "' is not in the specification inventory "
        "(Section 4.2); it is not implementable under this specification");
  const std::string id = entry->id;
  if (id == "PA-ISO") { *model = ModelId::Isotropic; return Status::Success(); }
  if (id == "PA-QFORM") { *model = ModelId::QForm; return Status::Success(); }
  if (id == "PA-EPS") { *model = ModelId::EpsilonForm; return Status::Success(); }
  if (id == "PA-EPREM") { *model = ModelId::Eprem; return Status::Success(); }
  return Status::Error(StatusCode::UnsupportedModel,
      "model '" + id + "' (Section " + entry->section + ", runtime state " +
      Name(entry->state) + ") is registered but outside Release A "
      "(implemented: PA-ISO, PA-QFORM, PA-EPS, PA-EPREM)");
}

// ---------------------------------------------------------------------------
// 2. Kinematics (Section 1.6, Eq. N-set of Section 16.12 "Frames, clock and
//    state").
// ---------------------------------------------------------------------------

Status ComputeKinematics(const ParticleState& particle, Kinematics* out) {
  if (out == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null kinematics output");
  if (!std::isfinite(particle.massKg) || particle.massKg <= 0.0)
    return Status::Error(StatusCode::InvalidInput,
                         "particle mass must be finite and positive");
  if (!std::isfinite(particle.chargeC) || particle.chargeC == 0.0)
    return Status::Error(StatusCode::InvalidInput,
        "particle charge must be finite and nonzero (rigidity pc/|q|)");
  if (!std::isfinite(particle.momentumKgMPerS) ||
      particle.momentumKgMPerS <= 0.0)
    return Status::Error(StatusCode::InvalidInput,
        "particle momentum must be finite and positive (p = 0 has no "
        "pitch-angle direction, Section 16.12)");
  const double c = SpeedOfLightMPerS;
  const double mc = particle.massKg * c;
  // gamma = sqrt(1 + (p/mc)^2) via hypot avoids overflow of p^2.
  const double pOverMc = particle.momentumKgMPerS / mc;
  const double gamma = std::hypot(1.0, pOverMc);
  Kinematics k;
  k.gamma = gamma;
  // v = p c^2 / E = c (p/mc)/gamma, exact and stable for both limits.
  k.speedMPerS = c * pOverMc / gamma;
  k.rigidityV = particle.momentumKgMPerS * c / std::fabs(particle.chargeC);
  k.totalEnergyJ = gamma * mc * c;
  // T = p^2 c^2/(E + m c^2): no cancellation at p << mc (Section 16.12).
  const double pc = particle.momentumKgMPerS * c;
  k.kineticEnergyJ = pc * pc / (k.totalEnergyJ + mc * c);
  if (!std::isfinite(k.speedMPerS) || !(k.speedMPerS > 0.0) ||
      !std::isfinite(k.rigidityV) || !std::isfinite(k.kineticEnergyJ))
    return Status::Error(StatusCode::InvalidInput,
                         "particle kinematics are not representable");
  *out = k;
  return Status::Success();
}

// ---------------------------------------------------------------------------
// 3. Adaptive Gauss-Kronrod quadrature: pitch_angle_quadrature.h.
// ---------------------------------------------------------------------------
using Quadrature::AdaptiveIntegrate;

// ---------------------------------------------------------------------------
// 4. Normalization integrals (Section 5.6).
// ---------------------------------------------------------------------------
namespace {

// 2F1(1,1;c;z) = sum_n n!/(c)_n z^n for c > 1, 0 < z < 1.  The term ratio
// (n+1) z/(c+n) is below z for c > 1, so the tail after the last included
// term t_N is bounded by t_{N+1}/(1 - z) (geometric majorant).  The loop stops
// when that bound is below 2e-17 of the partial sum (derived).
bool Hypergeometric11(double c, double z, long double* sum, long long* terms) {
  if (!(c > 1.0) || !(z > 0.0) || !(z < 1.0)) return false;
  const long double oneMinusZ = 1.0L - static_cast<long double>(z);
  long double term = 1.0L, total = 1.0L;
  const long long maxTerms = 200000000LL;
  for (long long n = 0; n < maxTerms; ++n) {
    term *= static_cast<long double>(n + 1) * z /
            (static_cast<long double>(c) + static_cast<long double>(n));
    total += term;
    if (term / oneMinusZ <= 2.0e-17L * total) {
      *sum = total;
      *terms = n + 1;
      return true;
    }
  }
  return false;
}

bool NearlyFiveThirds(double q) {
  return std::fabs(q - 5.0 / 3.0) <= 4.0 * std::numeric_limits<double>::epsilon();
}

}  // namespace

Status QFormIntegralClosedForm(double q, double gapH, double* value,
                               std::string* method) {
  if (value == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null value output");
  if (!std::isfinite(q) || !(q > 1.0) || !(q < 2.0) || !std::isfinite(gapH) ||
      gapH < 0.0)
    return Status::Error(StatusCode::OutsideDomain,
                         "I(q,H) requires 1 < q < 2 and H >= 0");
  if (gapH == 0.0) {
    // Eq. (27): I(q,0) = 2/((2-q)(4-q)), exact.
    *value = 2.0 / ((2.0 - q) * (4.0 - q));
    if (method) *method = "Eq. (27) exact, H = 0";
    return Status::Success();
  }
  if (NearlyFiveThirds(q) && gapH <= 10.0) {
    // Eq. (38), substitution mu = t^3:
    // I(5/3,H) = 3 [6/7 + H/5 - H^2/3 + H^3 - (H + H^4) atan(1/sqrt H)/sqrt H].
    // Its double-precision error is <= 2.3e-16 for H <= 0.5 and grows by
    // cancellation to 1.5e-12 at H = 10 (Section 5.6); larger H uses the
    // series, which converges fast there.
    const double s = std::sqrt(gapH);
    const double h2 = gapH * gapH;
    *value = 3.0 * (6.0 / 7.0 + gapH / 5.0 - h2 / 3.0 + h2 * gapH -
                    (gapH + h2 * h2) * std::atan(1.0 / s) / s);
    if (method) *method = "Eq. (38), q = 5/3";
    return Status::Success();
  }
  // Second form of Eq. (37): argument z = 1/(1+H) in (0,1),
  // I = [2F1(1,1;q/(q-1);z) - (1/3) 2F1(1,1;(q+2)/(q-1);z)]/(1+H).
  const double z = 1.0 / (1.0 + gapH);
  long double first = 0.0L, second = 0.0L;
  long long n1 = 0, n2 = 0;
  if (!Hypergeometric11(q / (q - 1.0), z, &first, &n1) ||
      !Hypergeometric11((q + 2.0) / (q - 1.0), z, &second, &n2))
    return Status::Error(StatusCode::Nonconverged,
        "Eq. (37) series for I(q,H) needs more than 2e8 terms at this H; "
        "no certified closed-form value");
  *value = static_cast<double>((first - second / 3.0L) /
                               (1.0L + static_cast<long double>(gapH)));
  if (method)
    *method = "Eq. (37) second form, 2F1 series (" + std::to_string(n1) +
              "/" + std::to_string(n2) + " terms)";
  return Status::Success();
}

Status QFormIntegralQuadrature(double q, double gapH, double* value) {
  if (value == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null value output");
  if (!std::isfinite(q) || !(q > 1.0) || !(q < 2.0) || !std::isfinite(gapH) ||
      gapH < 0.0)
    return Status::Error(StatusCode::OutsideDomain,
                         "I(q,H) requires 1 < q < 2 and H >= 0");
  // mu = t^m, m = 1/(2-q), alpha = q-1 (derived): dmu = m t^(m-1) dt and
  // m - 1 = m alpha, so (1-mu^2)/(mu^alpha + H) dmu
  //   = m (1 - t^(2m)) t^(m alpha)/(t^(m alpha) + H) dt,
  // which is bounded on [0,1] for every H >= 0 (= m (1 - t^(2m)) at H = 0).
  const long double m = 1.0L / (2.0L - static_cast<long double>(q));
  const long double ma = m * (static_cast<long double>(q) - 1.0L);
  const long double h = gapH;
  auto integrand = [m, ma, h](long double t) {
    const long double tp = std::pow(t, ma);
    const long double shape = (h == 0.0L) ? 1.0L : tp / (tp + h);
    return m * (1.0L - std::pow(t, 2.0L * m)) * shape;
  };
  std::vector<long double> points = {0.0L, 1.0L};
  // The integrand turns over where t^(m alpha) ~ H; a breakpoint there lets
  // the bisection resolve the knee for small H.
  if (h > 0.0L) {
    const long double knee = std::pow(h, 1.0L / ma);
    if (knee > 0.0L && knee < 1.0L) points = {0.0L, knee, 1.0L};
  }
  long double result = 0.0L, error = 0.0L;
  if (!AdaptiveIntegrate(integrand, points, 1.0e-14L, 1.0e-300L, 4000,
                         &result, &error))
    return Status::Error(StatusCode::Nonconverged,
                         "adaptive quadrature of I(q,H) did not converge");
  *value = static_cast<double>(result);
  return Status::Success();
}

Status EpsilonPhiClosedForm(double epsilon, double* value) {
  if (value == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null value output");
  if (!std::isfinite(epsilon) || !(epsilon > 0.0))
    return Status::Error(StatusCode::OutsideDomain,
        "phi(eps) is finite only for eps > 0 (lambda diverges as "
        "(3/2) ln(1/eps) for eps -> 0, Section 3.4)");
  // Eq. (39), a = 1 + eps:
  // phi = 3/(2 a^4) [1/6 + 2 eps + 5/2 eps^2 + 2/3 eps^3
  //                  + (1 + 2 eps) ln((1 + 2 eps)/eps)].
  // All bracket terms are positive: no cancellation for any eps > 0.  The
  // logarithm is split as log1p(2 eps) - log(eps) for accuracy at small eps.
  const long double e = epsilon;
  const long double a = 1.0L + e;
  const long double a2 = a * a;
  const long double bracket = 1.0L / 6.0L + 2.0L * e + 2.5L * e * e +
      (2.0L / 3.0L) * e * e * e +
      (1.0L + 2.0L * e) * (std::log1p(2.0L * e) - std::log(e));
  *value = static_cast<double>(1.5L * bracket / (a2 * a2));
  return Status::Success();
}

Status EpsilonPhiQuadrature(double epsilon, double* value) {
  if (value == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null value output");
  if (!std::isfinite(epsilon) || !(epsilon > 0.0))
    return Status::Error(StatusCode::OutsideDomain,
                         "phi(eps) is finite only for eps > 0");
  // phi = (3/2) Int_0^1 (1-mu)(1+mu)^2/(a mu + eps) dmu (Eq. 39, first form).
  // mu = (eps/a)(e^u - 1) maps the near-pole at mu = -eps/a to a smooth
  // integrand: dmu/(a mu + eps) = du/a, u in [0, ln(1 + a/eps)] (derived).
  const long double e = epsilon;
  const long double a = 1.0L + e;
  const long double upper = std::log1p(a / e);
  auto integrand = [e, a](long double u) {
    const long double mu = (e / a) * std::expm1(u);
    return 1.5L * (1.0L - mu) * (1.0L + mu) * (1.0L + mu) / a;
  };
  long double result = 0.0L, error = 0.0L;
  if (!AdaptiveIntegrate(integrand, {0.0L, upper}, 1.0e-14L, 1.0e-300L, 4000,
                         &result, &error))
    return Status::Error(StatusCode::Nonconverged,
                         "adaptive quadrature of phi(eps) did not converge");
  *value = static_cast<double>(result);
  return Status::Success();
}

// ---------------------------------------------------------------------------
// 5. Shapes.  D(mu) = D0 g(mu), g = (1 - mu^2) h(mu):
//    PA-ISO, PA-EPREM  h = 1,                     g' = -2 mu;
//    PA-QFORM          h = |mu|^(q-1) + H,        Eq. (163);
//    PA-EPS            h = |mu|/(1+|mu|) + eps,   Eq. (164).
// ---------------------------------------------------------------------------
namespace {

double Sign(double x) { return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); }

double ShapeH(ModelId model, const ShapeParameters& shape, double mu) {
  const double a = std::fabs(mu);
  switch (model) {
    case ModelId::Isotropic:
    case ModelId::Eprem:
      return 1.0;
    case ModelId::QForm:
      return std::pow(a, shape.q - 1.0) + shape.gapH;
    case ModelId::EpsilonForm:
      return a / (1.0 + a) + shape.epsilon;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// g'(mu) away from mu = 0, and at mu = 0 the printed formula with sgn(0)=0
// where it is finite.  Returns false where no finite value exists (PA-QFORM
// at mu = 0, Eq. 163: the term (q-1)(1-mu^2)|mu|^(q-2) sgn(mu) diverges).
bool ShapeGPrime(ModelId model, const ShapeParameters& shape, double mu,
                 double* value) {
  const double a = std::fabs(mu);
  const double h = ShapeH(model, shape, mu);
  switch (model) {
    case ModelId::Isotropic:
    case ModelId::Eprem:
      *value = -2.0 * mu;
      return true;
    case ModelId::QForm:
      if (a == 0.0) return false;
      // Eq. (163) with D0 factored out.
      *value = -2.0 * mu * h +
               (shape.q - 1.0) * (1.0 - mu * mu) * std::pow(a, shape.q - 2.0) *
                   Sign(mu);
      return true;
    case ModelId::EpsilonForm:
      // Eq. (164) with D0 factored out; at mu = 0 the printed sgn(0) = 0
      // gives 0, the mean of the one-sided limits -1 and +1.
      *value = -2.0 * mu * h +
               (1.0 - mu * mu) * Sign(mu) / ((1.0 + a) * (1.0 + a));
      return true;
  }
  return false;
}

// Unique maximum of g on [0, 1] for shapes whose g' is strictly decreasing
// there (PA-QFORM, PA-EPS: g'' = -2h - 4 mu h' + (1-mu^2) h'' < 0 because
// h > 0, h' > 0 and h'' <= 0 on (0,1], derived).  g' > 0 near 0+ and
// g'(1) = -2 h(1) < 0, so the maximizer is the root of g', found by
// bisection to machine resolution.
double MaximumOfG(ModelId model, const ShapeParameters& shape) {
  double low = 0.0, high = 1.0;
  for (int i = 0; i < 200; ++i) {
    const double middle = 0.5 * (low + high);
    if (!(middle > low && middle < high)) break;
    double slope = 0.0;
    if (!ShapeGPrime(model, shape, middle, &slope)) { low = middle; continue; }
    if (slope > 0.0) low = middle; else high = middle;
  }
  const double mu = 0.5 * (low + high);
  return (1.0 - mu * mu) * ShapeH(model, shape, mu);
}

}  // namespace

// ---------------------------------------------------------------------------
// 6. Validation, evaluation and step limits.
// ---------------------------------------------------------------------------
namespace {

Status CheckLambdaLaw(const LambdaLawParameters& law) {
  if (!std::isfinite(law.lambdaReferenceM) || !(law.lambdaReferenceM > 0.0))
    return Status::Error(StatusCode::OutsideDomain,
                         "lambda_reference_m must be finite and positive");
  if (!std::isfinite(law.referenceRigidityV) || !(law.referenceRigidityV > 0.0))
    return Status::Error(StatusCode::OutsideDomain,
                         "lambda_reference_rigidity_v must be finite and positive");
  if (!std::isfinite(law.referenceRadiusM) || !(law.referenceRadiusM > 0.0))
    return Status::Error(StatusCode::OutsideDomain,
                         "lambda_reference_radius_m must be finite and positive");
  if (!std::isfinite(law.rigidityExponent))
    return Status::Error(StatusCode::OutsideDomain,
                         "lambda_rigidity_exponent must be finite");
  if (!std::isfinite(law.radialExponent))
    return Status::Error(StatusCode::OutsideDomain,
                         "lambda_radial_exponent must be finite");
  return Status::Success();
}

Status CheckStepRule(const StepRuleParameters& step) {
  switch (step.rule) {
    case StepRule::RateFraction:
    case StepRule::VanDenBerg2026:
      if (!std::isfinite(step.rateFraction) || !(step.rateFraction > 0.0))
        return Status::Error(StatusCode::OutsideDomain,
                             "rate_fraction must be finite and positive");
      if (step.rule == StepRule::VanDenBerg2026 &&
          (!std::isfinite(step.ell) || !(step.ell > 0.0)))
        return Status::Error(StatusCode::OutsideDomain,
                             "ell must be finite and positive");
      return Status::Success();
    case StepRule::Sandroos2014:
      // dmu_max is a pitch-cosine increment: the interval [-1,1] has length 2.
      if (!std::isfinite(step.deltaMuMax) || !(step.deltaMuMax > 0.0) ||
          step.deltaMuMax > 2.0)
        return Status::Error(StatusCode::OutsideDomain,
                             "delta_mu_max must satisfy 0 < delta_mu_max <= 2");
      return Status::Success();
  }
  return Status::Error(StatusCode::InvalidConfiguration, "unknown step rule");
}

// Shape constants for an amplitude-free shape (also used by EvaluateShape).
Status ComputeShapeConstants(ModelId model, const ShapeParameters& shape,
                             bool requireFiniteLambda,
                             ShapeConstants* constants) {
  ShapeConstants k;
  const double infinity = std::numeric_limits<double>::infinity();
  switch (model) {
    case ModelId::Isotropic:
    case ModelId::Eprem:
      // PA-ISO: lambda = v/(2 D0), so N = 1/2 (Eq. 24).  PA-EPREM: printed
      // D^E = (1-mu^2) v/(2 lambda) with operator d/dmu(D^E/2 d/dmu), so the
      // standard D0 = v/(4 lambda) (N = 1/4) and Eq. (23) returns 2 lambda
      // (Section 5.4; MFP v1.2 D-2).
      k.normalization = model == ModelId::Isotropic ? 0.5 : 0.25;
      k.integralClosedForm = k.integralQuadrature = 1.0;
      k.closedFormMethod = model == ModelId::Isotropic
          ? "Eq. (24) exact" : "Eq. (34) exact (EPREM half-D convention)";
      k.lambdaEquation23OverLaw = model == ModelId::Isotropic ? 1.0 : 2.0;
      k.convention = model == ModelId::Isotropic
          ? OperatorConvention::Standard : OperatorConvention::EpremHalfD;
      k.supH = k.infH = 1.0;
      k.supAbsGPrime = 2.0;  // |g'| = 2|mu| <= 2
      k.supG = 1.0;          // g = 1 - mu^2 <= 1
      k.zero = ZeroClass::Regular;
      break;
    case ModelId::QForm: {
      if (!std::isfinite(shape.q) || !(shape.q > 1.0) || !(shape.q < 2.0))
        return Status::Error(StatusCode::OutsideDomain,
                             "PA-QFORM requires 1 < q < 2 (Eq. 26)");
      if (!std::isfinite(shape.gapH) || shape.gapH < 0.0)
        return Status::Error(StatusCode::OutsideDomain,
                             "PA-QFORM requires gap_h >= 0 (Eq. 26)");
      double closed = 0.0, quadrature = 0.0;
      std::string method;
      const Status closedStatus =
          QFormIntegralClosedForm(shape.q, shape.gapH, &closed, &method);
      if (!closedStatus.ok()) return closedStatus;
      const Status quadratureStatus =
          QFormIntegralQuadrature(shape.q, shape.gapH, &quadrature);
      if (!quadratureStatus.ok()) return quadratureStatus;
      k.integralClosedForm = closed;
      k.integralQuadrature = quadrature;
      k.closedFormMethod = method;
      k.integralRelativeDiscrepancy = std::fabs(closed - quadrature) / closed;
      // Eq. (27): D0 = (3 v/(4 lambda)) I(q,H).
      k.normalization = 0.75 * closed;
      k.supH = 1.0 + shape.gapH;   // h increases with |mu|
      k.infH = shape.gapH;         // at mu = 0
      k.supAbsGPrime = infinity;   // Eq. (163) diverges at mu = 0
      k.supG = MaximumOfG(model, shape);
      k.zero = ZeroClass::UnboundedDrift;
      break;
    }
    case ModelId::EpsilonForm: {
      if (!std::isfinite(shape.epsilon) || shape.epsilon < 0.0)
        return Status::Error(StatusCode::OutsideDomain,
                             "PA-EPS requires epsilon >= 0 (Eq. 29)");
      if (shape.epsilon == 0.0) {
        if (requireFiniteLambda)
          return Status::Error(StatusCode::OutsideDomain,
              "PA-EPS with epsilon = 0 has a logarithmically divergent "
              "lambda_parallel (Section 3.4) and cannot be normalized to a "
              "lambda law; it disconnects the hemispheres (Section 16.5)");
        k.normalization = infinity;
        k.zero = ZeroClass::InternalBarrier;
      } else {
        double closed = 0.0, quadrature = 0.0;
        const Status closedStatus = EpsilonPhiClosedForm(shape.epsilon, &closed);
        if (!closedStatus.ok()) return closedStatus;
        const Status quadratureStatus =
            EpsilonPhiQuadrature(shape.epsilon, &quadrature);
        if (!quadratureStatus.ok()) return quadratureStatus;
        k.integralClosedForm = closed;
        k.integralQuadrature = quadrature;
        k.closedFormMethod = "Eq. (39)";
        k.integralRelativeDiscrepancy = std::fabs(closed - quadrature) / closed;
        // lambda = (v/nu0) phi and D0 = nu0/2, so D0 = phi v/(2 lambda).
        k.normalization = 0.5 * closed;
        k.zero = ZeroClass::DriftJump;
      }
      k.supH = 0.5 + shape.epsilon;  // h increases with |mu|, h(1) = 1/2 + eps
      k.infH = shape.epsilon;
      // g' is strictly decreasing on (0,1] (see MaximumOfG) from g'(0+) = 1
      // to g'(1) = -(1 + 2 eps): sup |g'| = 1 + 2 eps (derived; Eq. 164 and
      // Section 16.6 "D'(1) = -D0(1 + 2 eps)").
      k.supAbsGPrime = 1.0 + 2.0 * shape.epsilon;
      k.supG = MaximumOfG(model, shape);
      break;
    }
  }
  if (constants) *constants = k;
  return Status::Success();
}

// lambda_law = lambda_ref (R/R_ref)^a (r/r_ref)^b, evaluated in log space.
Status LambdaLaw(const LambdaLawParameters& law, double rigidityV,
                 double radiusM, double* lambdaM) {
  if (!std::isfinite(radiusM) || !(radiusM > 0.0))
    return Status::Error(StatusCode::InvalidField,
        "the lambda law needs a positive heliocentric radius");
  const double logLambda = std::log(law.lambdaReferenceM) +
      law.rigidityExponent * std::log(rigidityV / law.referenceRigidityV) +
      law.radialExponent * std::log(radiusM / law.referenceRadiusM);
  const double lambda = std::exp(logLambda);
  if (!std::isfinite(lambda) || !(lambda > 0.0))
    return Status::Error(StatusCode::OutsideDomain,
        "the lambda law value is not representable at this particle state");
  *lambdaM = lambda;
  return Status::Success();
}

CoefficientResult Failure(StatusCode code, const std::string& detail,
                          const char* modelId) {
  CoefficientResult result;
  result.status = Status::Error(code, detail);
  result.provenance.modelId = modelId;
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.specificationSha256 = SpecificationSha256();
  result.provenance.softwareVersion = SoftwareVersion;
  return result;
}

}  // namespace

bool ShapeDmumu(double mu, double amplitudeD0PerS, ModelId model,
                const ShapeParameters& shape, double* dMuMuPerS,
                double* dDmuMuDmuPerS) {
  *dMuMuPerS = amplitudeD0PerS * (1.0 - mu * mu) * ShapeH(model, shape, mu);
  double gPrime = 0.0;
  if (!ShapeGPrime(model, shape, mu, &gPrime)) return false;
  *dDmuMuDmuPerS = amplitudeD0PerS * gPrime;
  return true;
}

namespace {

// Fill D, D' and the derived frequencies for amplitude D0 at mu.
void FillShapeValues(ModelId model, const ShapeParameters& shape,
                     const ShapeConstants& constants, double d0, double mu,
                     CoefficientResult* result) {
  const double h = ShapeH(model, shape, mu);
  result->amplitudeD0PerS = d0;
  result->localScatteringFrequencyPerS = 2.0 * d0 * h;
  result->maximumScatteringFrequencyPerS = 2.0 * d0 * constants.supH;
  result->zero = constants.zero;
  double d = 0.0, dPrime = 0.0;
  const bool finite = ShapeDmumu(mu, d0, model, shape, &d, &dPrime);
  result->dMuMuPerS = d;
  if (finite) result->dDmuMuDmuPerS = dPrime;
  if (mu != 0.0) {
    if (result->dDmuMuDmuPerS) {
      result->derivativeMinusPerS = *result->dDmuMuDmuPerS;
      result->derivativePlusPerS = *result->dDmuMuDmuPerS;
    }
  } else if (model == ModelId::EpsilonForm) {
    // One-sided limits of Eq. (164) at 0: -D0 and +D0 (F-EMP-12).
    result->derivativeMinusPerS = -d0;
    result->derivativePlusPerS = d0;
  } else if (model == ModelId::Isotropic || model == ModelId::Eprem) {
    result->derivativeMinusPerS = 0.0;
    result->derivativePlusPerS = 0.0;
  }
  // PA-QFORM at mu = 0: both one-sided derivatives are infinite; absent.
}

}  // namespace

Status ValidateConfiguration(const ModelConfiguration& configuration,
                             ShapeConstants* constants) {
  switch (configuration.model) {
    case ModelId::Isotropic:
    case ModelId::QForm:
    case ModelId::EpsilonForm:
    case ModelId::Eprem:
      break;
    default:
      return Status::Error(StatusCode::InvalidConfiguration,
                           "unknown pitch-angle model");
  }
  const Status law = CheckLambdaLaw(configuration.lambdaLaw);
  if (!law.ok()) return law;
  const Status step = CheckStepRule(configuration.step);
  if (!step.ok()) return step;
  ShapeConstants k;
  const Status shape = ComputeShapeConstants(configuration.model,
                                             configuration.shape, true, &k);
  if (!shape.ok()) return shape;
  // Section 5.7: the closed form and an independent quadrature are checked
  // against each other for every configured parameter set; the budget is the
  // roadmap's smooth-quadrature tolerance (Section 26.5).
  if (!(k.integralRelativeDiscrepancy <= 1.0e-8)) {
    std::ostringstream message;
    message << std::setprecision(17) << "normalization cross-check failed: "
            << k.closedFormMethod << " = " << k.integralClosedForm
            << ", quadrature = " << k.integralQuadrature;
    return Status::Error(StatusCode::Nonconverged, message.str());
  }
  if (constants) *constants = k;
  return Status::Success();
}

CoefficientResult EvaluateShape(double mu, double amplitudeD0PerS,
                                ModelId model, const ShapeParameters& shape) {
  const char* id = ModelName(model);
  if (!std::isfinite(mu) || mu < -1.0 || mu > 1.0)
    return Failure(StatusCode::InvalidInput, "mu must lie in [-1, 1]", id);
  if (!std::isfinite(amplitudeD0PerS) || amplitudeD0PerS < 0.0)
    return Failure(StatusCode::InvalidInput,
                   "the amplitude D0 must be finite and nonnegative", id);
  // Light path: only the domain checks, sup h and the zero class are needed
  // to fill the values; the normalization integrals (which EvaluateShape
  // does not use) are not recomputed per call.
  ShapeConstants k;
  switch (model) {
    case ModelId::Isotropic:
    case ModelId::Eprem:
      k.supH = 1.0;
      k.zero = ZeroClass::Regular;
      break;
    case ModelId::QForm:
      if (!std::isfinite(shape.q) || !(shape.q > 1.0) || !(shape.q < 2.0) ||
          !std::isfinite(shape.gapH) || shape.gapH < 0.0)
        return Failure(StatusCode::OutsideDomain,
                       "PA-QFORM requires 1 < q < 2 and gap_h >= 0", id);
      k.supH = 1.0 + shape.gapH;
      k.zero = ZeroClass::UnboundedDrift;
      break;
    case ModelId::EpsilonForm:
      if (!std::isfinite(shape.epsilon) || shape.epsilon < 0.0)
        return Failure(StatusCode::OutsideDomain,
                       "PA-EPS requires epsilon >= 0", id);
      k.supH = 0.5 + shape.epsilon;
      k.zero = shape.epsilon == 0.0 ? ZeroClass::InternalBarrier
                                    : ZeroClass::DriftJump;
      break;
    default:
      return Failure(StatusCode::InvalidConfiguration, "unknown model", id);
  }
  CoefficientResult result;
  result.status = Status::Success();
  FillShapeValues(model, shape, k, amplitudeD0PerS, mu, &result);
  result.provenance.modelId = id;
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.specificationSha256 = SpecificationSha256();
  result.provenance.softwareVersion = SoftwareVersion;
  return result;
}

namespace {

// Evaluation with already validated configuration and constants (the active
// path and Evaluate share it).
CoefficientResult EvaluateValidated(double mu, const ParticleState& particle,
                                    const LocalState& local,
                                    const ModelConfiguration& configuration,
                                    const ShapeConstants& constants,
                                    const std::string& fingerprint) {
  const char* id = ModelName(configuration.model);
  if (!std::isfinite(mu) || mu < -1.0 || mu > 1.0)
    return Failure(StatusCode::InvalidInput, "mu must lie in [-1, 1]", id);
  Kinematics kinematics;
  const Status kin = ComputeKinematics(particle, &kinematics);
  if (!kin.ok()) return Failure(kin.code, kin.detail, id);
  const std::array<double, 3>& x = local.heliocentricPositionM;
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    return Failure(StatusCode::InvalidField,
                   "heliocentric position is not finite", id);
  const double radius = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double lambda = 0.0;
  const Status law = LambdaLaw(configuration.lambdaLaw, kinematics.rigidityV,
                               radius, &lambda);
  if (!law.ok()) return Failure(law.code, law.detail, id);
  // D0 = N_shape v / lambda_law (Section 5.6 via Eq. 23).
  const double d0 = constants.normalization * kinematics.speedMPerS / lambda;
  if (!std::isfinite(d0) || !(d0 > 0.0))
    return Failure(StatusCode::OutsideDomain,
                   "the amplitude D0 is not representable", id);
  CoefficientResult result;
  result.status = Status::Success();
  FillShapeValues(configuration.model, configuration.shape, constants, d0, mu,
                  &result);
  result.lambdaStatus = LambdaStatus::Finite;
  result.lambdaLawM = lambda;
  result.lambdaParallelUnfocusedM = constants.lambdaEquation23OverLaw * lambda;
  result.kinematics = kinematics;
  result.provenance.modelId = id;
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.specificationSha256 = SpecificationSha256();
  result.provenance.softwareVersion = SoftwareVersion;
  const std::size_t length = std::min<std::size_t>(fingerprint.size(), 64);
  fingerprint.copy(result.provenance.configurationFingerprint, length);
  result.provenance.configurationFingerprint[length] = '\0';
  result.provenance.backgroundRevision = local.backgroundRevision;
  result.provenance.turbulenceRevision = local.turbulenceRevision;
  return result;
}

StepLimitResult StepFailure(StatusCode code, const std::string& detail) {
  StepLimitResult result;
  result.status = Status::Error(code, detail);
  return result;
}

StepLimitResult StepLimitValidated(double mu, double deterministicMuRatePerS,
                                   const ParticleState& particle,
                                   const LocalState& local,
                                   const ModelConfiguration& configuration,
                                   const ShapeConstants& constants) {
  if (configuration.model == ModelId::QForm)
    return StepFailure(StatusCode::IncompatibleKernel,
        "PA-QFORM has an unbounded Ito drift at mu = 0 (Eq. 163); no Ito "
        "step rule is released for it (Section 16.11)");
  if (!std::isfinite(deterministicMuRatePerS))
    return StepFailure(StatusCode::InvalidInput,
                       "deterministic pitch-angle rate is not finite");
  const CoefficientResult c = EvaluateValidated(
      mu, particle, local, configuration, constants, std::string());
  if (!c.status.ok()) return StepFailure(c.status.code, c.status.detail);
  const double d0 = c.amplitudeD0PerS;
  const double nuMax = c.maximumScatteringFrequencyPerS;
  const double infinity = std::numeric_limits<double>::infinity();
  StepLimitResult result;
  const StepRuleParameters& step = configuration.step;
  result.driftEntryS = infinity;
  result.noiseEntryS = infinity;
  result.rateFractionS = step.rule == StepRule::Sandroos2014
      ? infinity : step.rateFraction / nuMax;
  switch (step.rule) {
    case StepRule::RateFraction:
      result.maximumStepS = result.rateFractionS;
      result.limitingEntry = "rate-fraction";
      break;
    case StepRule::VanDenBerg2026: {
      // Eq. (165) in forward time: A = dD/dmu + dmu/dt (the Ito drift of
      // Eq. 25; the source prints the backward drift).  The noise entry is
      // written with the factored D = D0 (1-mu^2) h so that it stays finite
      // at |mu| = 1: ell^2 (1-mu^2)/(2D) = ell^2/(2 D0 h).
      if (!c.dDmuMuDmuPerS)
        return StepFailure(StatusCode::IncompatibleKernel,
                           "no finite drift for the van den Berg rule");
      const double drift = *c.dDmuMuDmuPerS + deterministicMuRatePerS;
      const double root = std::sqrt(std::max(0.0, 1.0 - mu * mu));
      const double driftEntry = drift != 0.0
          ? step.ell * root / std::fabs(drift) : infinity;
      const double h = c.localScatteringFrequencyPerS / (2.0 * d0);
      const double noiseEntry = h > 0.0
          ? step.ell * step.ell / (2.0 * d0 * h) : infinity;
      result.driftEntryS = driftEntry;
      result.noiseEntryS = noiseEntry;
      double value = std::min(driftEntry, noiseEntry);
      result.limitingEntry = driftEntry <= noiseEntry ? "drift" : "noise";
      // Floor everywhere (user decision 2026-10-13): Eq. (165) is exactly 0
      // at |mu| = 1 whenever D'(+-1) != 0 (Section 16.6, U-55), which would
      // stall an endpoint beam.
      const double floor = result.rateFractionS;
      if (!(value >= floor)) {
        value = floor;
        result.limitingEntry = "rate-fraction-floor";
      }
      result.maximumStepS = value;
      break;
    }
    case StepRule::Sandroos2014: {
      // Eq. (166) with a = dD/dmu and b^2 = 2D maximized over mu (the
      // printed a_max, b_max) at this particle's amplitude.
      const double aMax = d0 * constants.supAbsGPrime;
      const double b2Max = 2.0 * d0 * constants.supG;
      const double driftEntry = step.deltaMuMax / aMax;
      const double noiseEntry = step.deltaMuMax * step.deltaMuMax / b2Max;
      result.driftEntryS = driftEntry;
      result.noiseEntryS = noiseEntry;
      result.maximumStepS = std::min(driftEntry, noiseEntry);
      result.limitingEntry =
          driftEntry <= noiseEntry ? "sandroos-drift" : "sandroos-noise";
      break;
    }
  }
  if (!std::isfinite(result.maximumStepS) || !(result.maximumStepS > 0.0))
    return StepFailure(StatusCode::Nonprogress,
                       "scattering step rule produced no positive interval");
  result.status = Status::Success();
  return result;
}

}  // namespace

CoefficientResult Evaluate(double mu, const ParticleState& particle,
                           const LocalState& local,
                           const ModelConfiguration& configuration) {
  ShapeConstants constants;
  const Status valid = ValidateConfiguration(configuration, &constants);
  if (!valid.ok())
    return Failure(valid.code, valid.detail, ModelName(configuration.model));
  return EvaluateValidated(mu, particle, local, configuration, constants,
                           ConfigurationFingerprint(configuration));
}

StepLimitResult ScatteringStepLimit(double mu, double deterministicMuRatePerS,
                                    const ParticleState& particle,
                                    const LocalState& local,
                                    const ModelConfiguration& configuration) {
  ShapeConstants constants;
  const Status valid = ValidateConfiguration(configuration, &constants);
  if (!valid.ok()) return StepFailure(valid.code, valid.detail);
  return StepLimitValidated(mu, deterministicMuRatePerS, particle, local,
                            configuration, constants);
}

double WorstCaseScaledStep(const ModelConfiguration& configuration,
                           const ShapeConstants& constants) {
  const double infinity = std::numeric_limits<double>::infinity();
  const StepRuleParameters& step = configuration.step;
  switch (step.rule) {
    case StepRule::RateFraction:
      return step.rateFraction;
    case StepRule::VanDenBerg2026: {
      // nu_max dt <= nu_max ell^2/(2 D0 h(mu)) = ell^2 supH/h(mu) for the
      // noise entry; the drift entry can only shorten the step (a host drift
      // may cancel D', so it is not used for the bound).  The floor adds f.
      const double noise = constants.infH > 0.0
          ? step.ell * step.ell * constants.supH / constants.infH : infinity;
      return std::max(step.rateFraction, noise);
    }
    case StepRule::Sandroos2014: {
      const double drift = std::isfinite(constants.supAbsGPrime)
          ? step.deltaMuMax * 2.0 * constants.supH / constants.supAbsGPrime
          : 0.0;
      const double noise =
          step.deltaMuMax * step.deltaMuMax * constants.supH / constants.supG;
      return std::min(drift > 0.0 ? drift : infinity, noise);
    }
  }
  return infinity;
}

// ---------------------------------------------------------------------------
// 7. Certification matrix (PS05): pitch_angle_certification.cpp.  It is a
//    separate translation unit so that recording new certification evidence
//    does not change the objects (and the digest) of the acceptance binary
//    that produced the evidence.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// 8. Parser-neutral bridge, section parser and fingerprint.
// ---------------------------------------------------------------------------
namespace {

bool ParseDecimal(const std::string& text, double* value) {
  if (text.empty()) return false;
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed))
    return false;
  *value = parsed;
  return true;
}

// Decimal, or (dimensionless keys only) an exact rational "a/b" of two
// integers, b > 0.
bool ParseNumber(const std::string& text, bool allowRational, double* value) {
  const std::size_t slash = text.find('/');
  if (slash == std::string::npos) return ParseDecimal(text, value);
  if (!allowRational) return false;
  const std::string numerator = text.substr(0, slash);
  const std::string denominator = text.substr(slash + 1);
  auto integer = [](const std::string& s, bool allowSign) {
    if (s.empty()) return false;
    std::size_t start = 0;
    if (allowSign && (s[0] == '-' || s[0] == '+')) start = 1;
    if (start == s.size()) return false;
    for (std::size_t i = start; i < s.size(); ++i)
      if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    return true;
  };
  if (!integer(numerator, true) || !integer(denominator, false)) return false;
  const double n = std::strtod(numerator.c_str(), nullptr);
  const double d = std::strtod(denominator.c_str(), nullptr);
  if (!(d > 0.0) || !std::isfinite(n) || !std::isfinite(d)) return false;
  *value = n / d;
  return std::isfinite(*value);
}

class KeyReader {
 public:
  explicit KeyReader(const std::vector<InputParameter>& parameters)
      : parameters_(parameters) {}

  Status Prepare() {
    for (const InputParameter& parameter : parameters_) {
      if (parameter.name.empty())
        return Status::Error(StatusCode::InvalidConfiguration,
                             "pitch-angle parameter name is empty");
      if (!values_.emplace(parameter.name, parameter.value).second)
        return Status::Error(StatusCode::InvalidConfiguration,
            "duplicate pitch-angle parameter '" + parameter.name + "'");
    }
    return Status::Success();
  }

  Status Number(const char* key, bool dimensionless, double* output) {
    const auto found = values_.find(key);
    if (found == values_.end())
      return Status::Error(StatusCode::MissingInput,
                           std::string("missing parameter '") + key + "'");
    if (!ParseNumber(found->second, dimensionless, output))
      return Status::Error(StatusCode::InvalidConfiguration,
          std::string("parameter '") + key + "' is not a finite " +
          (dimensionless ? "number or integer ratio a/b" : "decimal number"));
    values_.erase(found);
    return Status::Success();
  }

  Status Text(const char* key, std::string* output) {
    const auto found = values_.find(key);
    if (found == values_.end())
      return Status::Error(StatusCode::MissingInput,
                           std::string("missing parameter '") + key + "'");
    *output = found->second;
    values_.erase(found);
    return Status::Success();
  }

  Status RejectRemaining(const std::string& context) const {
    if (values_.empty()) return Status::Success();
    return Status::Error(StatusCode::InvalidConfiguration,
        "parameter '" + values_.begin()->first + "' is not part of " + context);
  }

 private:
  const std::vector<InputParameter>& parameters_;
  std::map<std::string, std::string> values_;
};

void FingerprintNumber(std::ostringstream* output, const char* name,
                       double value) {
  *output << ';' << name << '=' << std::setprecision(17) << value;
}

std::uint32_t RotateRight(std::uint32_t value, unsigned shift) {
  return (value >> shift) | (value << (32u - shift));
}

// SHA-256 (FIPS 180-4) of the canonical serialization; same algorithm and
// role as the parallel-diffusion library's configuration identity.
std::string Sha256(const std::string& text) {
  static constexpr std::uint32_t k[64] = {
      0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
      0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
      0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
      0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
      0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
      0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
      0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
      0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u};
  std::vector<unsigned char> message(text.begin(), text.end());
  const std::uint64_t bitLength = static_cast<std::uint64_t>(message.size()) * 8u;
  message.push_back(0x80u);
  while ((message.size() % 64u) != 56u) message.push_back(0u);
  for (int shift = 56; shift >= 0; shift -= 8)
    message.push_back(static_cast<unsigned char>(bitLength >> shift));
  std::uint32_t h[8] = {0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,
                        0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u};
  for (std::size_t offset = 0; offset < message.size(); offset += 64u) {
    std::uint32_t w[64];
    for (int i = 0; i < 16; ++i) {
      const std::size_t p = offset + static_cast<std::size_t>(4 * i);
      w[i] = (static_cast<std::uint32_t>(message[p]) << 24) |
             (static_cast<std::uint32_t>(message[p + 1]) << 16) |
             (static_cast<std::uint32_t>(message[p + 2]) << 8) |
             static_cast<std::uint32_t>(message[p + 3]);
    }
    for (int i = 16; i < 64; ++i) {
      const std::uint32_t s0 = RotateRight(w[i - 15], 7) ^
          RotateRight(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = RotateRight(w[i - 2], 17) ^
          RotateRight(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
    for (int i = 0; i < 64; ++i) {
      const std::uint32_t s1 = RotateRight(e,6)^RotateRight(e,11)^RotateRight(e,25);
      const std::uint32_t ch = (e & f) ^ ((~e) & g);
      const std::uint32_t t1 = hh + s1 + ch + k[i] + w[i];
      const std::uint32_t s0 = RotateRight(a,2)^RotateRight(a,13)^RotateRight(a,22);
      const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = s0 + maj;
      hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d;
    h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
  }
  std::ostringstream result;
  result << std::hex << std::setfill('0');
  for (std::uint32_t value : h) result << std::setw(8) << value;
  return result.str();
}

std::string LinePrefix(std::size_t lineNumber) {
  return lineNumber == 0 ? std::string()
                         : "line " + std::to_string(lineNumber) + ": ";
}

}  // namespace

Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          ModelConfiguration* configuration) {
  if (configuration == nullptr)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null configuration output");
  ModelConfiguration candidate;
  const Status model = ParseModelId(modelId, &candidate.model);
  if (!model.ok()) return model;
  KeyReader keys(parameters);
  Status status = keys.Prepare();
  if (!status.ok()) return status;

  // Shape keys of the selected model only.
  switch (candidate.model) {
    case ModelId::QForm:
      if (!(status = keys.Number("q", true, &candidate.shape.q)).ok()) return status;
      if (!(status = keys.Number("gap_h", true, &candidate.shape.gapH)).ok())
        return status;
      break;
    case ModelId::EpsilonForm:
      if (!(status = keys.Number("epsilon", true, &candidate.shape.epsilon)).ok())
        return status;
      break;
    case ModelId::Isotropic:
    case ModelId::Eprem:
      break;
  }
  // Lambda law (user decision 2026-10-13); every value required.
  LambdaLawParameters& law = candidate.lambdaLaw;
  if (!(status = keys.Number("lambda_reference_m", false,
                             &law.lambdaReferenceM)).ok()) return status;
  if (!(status = keys.Number("lambda_reference_rigidity_v", false,
                             &law.referenceRigidityV)).ok()) return status;
  if (!(status = keys.Number("lambda_rigidity_exponent", true,
                             &law.rigidityExponent)).ok()) return status;
  if (!(status = keys.Number("lambda_reference_radius_m", false,
                             &law.referenceRadiusM)).ok()) return status;
  if (!(status = keys.Number("lambda_radial_exponent", true,
                             &law.radialExponent)).ok()) return status;
  // Step rule and only its own parameters.
  std::string rule;
  if (!(status = keys.Text("step_rule", &rule)).ok()) return status;
  const std::string normalizedRule = Lower(Trim(rule));
  StepRuleParameters& step = candidate.step;
  if (normalizedRule == "rate-fraction") {
    step.rule = StepRule::RateFraction;
    if (!(status = keys.Number("rate_fraction", true, &step.rateFraction)).ok())
      return status;
  } else if (normalizedRule == "van-den-berg-2026") {
    step.rule = StepRule::VanDenBerg2026;
    if (!(status = keys.Number("rate_fraction", true, &step.rateFraction)).ok())
      return status;
    if (!(status = keys.Number("ell", true, &step.ell)).ok()) return status;
  } else if (normalizedRule == "sandroos-2014") {
    step.rule = StepRule::Sandroos2014;
    if (!(status = keys.Number("delta_mu_max", true, &step.deltaMuMax)).ok())
      return status;
  } else {
    return Status::Error(StatusCode::InvalidConfiguration,
        "parameter 'step_rule' must be rate-fraction, van-den-berg-2026 or "
        "sandroos-2014");
  }
  status = keys.RejectRemaining(std::string("model ") +
                                ModelName(candidate.model) + " with step_rule " +
                                Name(step.rule));
  if (!status.ok()) return status;
  status = ValidateConfiguration(candidate, nullptr);
  if (!status.ok()) return status;
  *configuration = candidate;
  return Status::Success();
}

std::string ConfigurationFingerprint(const ModelConfiguration& configuration) {
  // Canonical serialization of exactly the active inputs (17 significant
  // digits round-trip a double), plus the specification digest so a changed
  // contract changes the identity.
  std::ostringstream text;
  text << "pitch_angle_diffusion;spec=" << SpecificationSha256()
       << ";model=" << ModelName(configuration.model);
  switch (configuration.model) {
    case ModelId::QForm:
      FingerprintNumber(&text, "q", configuration.shape.q);
      FingerprintNumber(&text, "gap_h", configuration.shape.gapH);
      break;
    case ModelId::EpsilonForm:
      FingerprintNumber(&text, "epsilon", configuration.shape.epsilon);
      break;
    case ModelId::Isotropic:
    case ModelId::Eprem:
      break;
  }
  const LambdaLawParameters& law = configuration.lambdaLaw;
  FingerprintNumber(&text, "lambda_reference_m", law.lambdaReferenceM);
  FingerprintNumber(&text, "lambda_reference_rigidity_v", law.referenceRigidityV);
  FingerprintNumber(&text, "lambda_rigidity_exponent", law.rigidityExponent);
  FingerprintNumber(&text, "lambda_reference_radius_m", law.referenceRadiusM);
  FingerprintNumber(&text, "lambda_radial_exponent", law.radialExponent);
  const StepRuleParameters& step = configuration.step;
  text << ";step_rule=" << Name(step.rule);
  if (step.rule != StepRule::Sandroos2014)
    FingerprintNumber(&text, "rate_fraction", step.rateFraction);
  if (step.rule == StepRule::VanDenBerg2026)
    FingerprintNumber(&text, "ell", step.ell);
  if (step.rule == StepRule::Sandroos2014)
    FingerprintNumber(&text, "delta_mu_max", step.deltaMuMax);
  return Sha256(text.str());
}

Status ParseSection(const std::vector<SectionLine>& lines,
                    ParsedSection* result) {
  if (result == nullptr)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null parsed-section output");
  ParsedSection candidate;
  bool modelSeen = false;
  std::map<std::string, std::size_t> firstLine;
  for (const SectionLine& line : lines) {
    std::string text = line.text;
    const std::size_t comment = text.find_first_of("#!");
    if (comment != std::string::npos) text.erase(comment);
    text = Trim(text);
    if (text.empty()) continue;
    const std::size_t separator = text.find('=');
    if (separator == std::string::npos)
      return Status::Error(StatusCode::InvalidConfiguration,
          LinePrefix(line.lineNumber) +
          "expected 'key = value' in the pitch-angle-diffusion section");
    const std::string key = Trim(text.substr(0, separator));
    const std::string value = Trim(text.substr(separator + 1));
    if (key.empty() || value.empty())
      return Status::Error(StatusCode::InvalidConfiguration,
          LinePrefix(line.lineNumber) +
          "pitch-angle-diffusion key and value must both be nonempty");
    const auto inserted = firstLine.emplace(key, line.lineNumber);
    if (!inserted.second)
      return Status::Error(StatusCode::InvalidConfiguration,
          LinePrefix(line.lineNumber) + "duplicate pitch-angle-diffusion key '" +
          key + "'" +
          (inserted.first->second == 0 ? std::string()
               : " (first assigned at line " +
                 std::to_string(inserted.first->second) + ")"));
    if (key == "model") {
      modelSeen = true;
      candidate.modelId = value;
      candidate.modelLineNumber = line.lineNumber;
    } else {
      candidate.parameters.push_back({key, value, line.lineNumber});
    }
  }
  if (!modelSeen)
    return Status::Error(StatusCode::MissingInput,
        "the pitch-angle-diffusion section must select a model with "
        "'model = <inventory ID>'");
  std::vector<InputParameter> parameters;
  for (const SectionParameter& parameter : candidate.parameters)
    parameters.push_back({parameter.name, parameter.value});
  const Status built = BuildConfiguration(candidate.modelId, parameters,
                                          &candidate.configuration);
  if (!built.ok()) {
    // Attribute the diagnostic to the record whose quoted key it names, or
    // whose key appears as a whole identifier; otherwise to the model line.
    std::size_t lineNumber = candidate.modelLineNumber;
    if (built.code != StatusCode::UnsupportedModel) {
      bool attributed = false;
      for (const SectionParameter& parameter : candidate.parameters)
        if (!attributed &&
            built.detail.find("'" + parameter.name + "'") != std::string::npos) {
          lineNumber = parameter.lineNumber;
          attributed = true;
        }
      for (const SectionParameter& parameter : candidate.parameters) {
        if (attributed) break;
        auto identifier = [](char c) {
          return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
        };
        // Every occurrence is examined: the key "q" also occurs inside words
        // such as "requires", which are not whole-identifier matches.
        for (std::size_t at = built.detail.find(parameter.name);
             at != std::string::npos && !attributed;
             at = built.detail.find(parameter.name, at + 1)) {
          const std::size_t end = at + parameter.name.size();
          if ((at == 0 || !identifier(built.detail[at - 1])) &&
              (end == built.detail.size() || !identifier(built.detail[end]))) {
            lineNumber = parameter.lineNumber;
            attributed = true;
          }
        }
      }
    }
    return Status::Error(built.code, LinePrefix(lineNumber) + built.detail);
  }
  const Status valid =
      ValidateConfiguration(candidate.configuration, &candidate.constants);
  if (!valid.ok()) return valid;
  candidate.configurationFingerprint =
      ConfigurationFingerprint(candidate.configuration);
  *result = candidate;
  return Status::Success();
}

// ---------------------------------------------------------------------------
// 9. Active configuration and dispatch pointers.  Mutated only during serial
//    initialization; read-only afterwards.
// ---------------------------------------------------------------------------
namespace {

ModelConfiguration gActiveConfiguration;
ShapeConstants gActiveConstants;
std::string gActiveFingerprint;
bool gHasActiveConfiguration = false;

CoefficientResult EvaluateUnconfigured(double, const ParticleState&,
                                       const LocalState&) {
  return Failure(StatusCode::InvalidConfiguration,
                 "no pitch-angle-diffusion model is configured", "none");
}

StepLimitResult StepLimitUnconfigured(double, double, const ParticleState&,
                                      const LocalState&) {
  return StepFailure(StatusCode::InvalidConfiguration,
                     "no pitch-angle-diffusion model is configured");
}

// One evaluator per model identity, so the pointer identifies exactly one
// backend; the parameters are those installed with it.
template <ModelId Model>
CoefficientResult EvaluateBound(double mu, const ParticleState& particle,
                                const LocalState& local) {
  return EvaluateValidated(mu, particle, local, gActiveConfiguration,
                           gActiveConstants, gActiveFingerprint);
}

StepLimitResult StepLimitActive(double mu, double deterministicMuRatePerS,
                                const ParticleState& particle,
                                const LocalState& local) {
  return StepLimitValidated(mu, deterministicMuRatePerS, particle, local,
                            gActiveConfiguration, gActiveConstants);
}

}  // namespace

ActiveCoefficientEvaluator ActivePitchAngleDiffusion = &EvaluateUnconfigured;
ActiveStepLimitEvaluator ActiveScatteringStepLimit = &StepLimitUnconfigured;

ActiveCoefficientEvaluator BoundFunctionForModel(ModelId model) {
  switch (model) {
    case ModelId::Isotropic: return &EvaluateBound<ModelId::Isotropic>;
    case ModelId::QForm: return &EvaluateBound<ModelId::QForm>;
    case ModelId::EpsilonForm: return &EvaluateBound<ModelId::EpsilonForm>;
    case ModelId::Eprem: return &EvaluateBound<ModelId::Eprem>;
  }
  return nullptr;
}

Status SetActiveConfiguration(const ModelConfiguration& configuration) {
  ShapeConstants constants;
  const Status valid = ValidateConfiguration(configuration, &constants);
  if (!valid.ok()) return valid;
  ActiveCoefficientEvaluator bound = BoundFunctionForModel(configuration.model);
  if (bound == nullptr)
    return Status::Error(StatusCode::UnsupportedModel,
                         "selected model has no evaluator");
  // Parameters first, then the pointers that make them reachable.
  gActiveConfiguration = configuration;
  gActiveConstants = constants;
  gActiveFingerprint = ConfigurationFingerprint(configuration);
  ActivePitchAngleDiffusion = bound;
  ActiveScatteringStepLimit = &StepLimitActive;
  gHasActiveConfiguration = true;
  return Status::Success();
}

Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters) {
  ModelConfiguration configuration;
  const Status built = BuildConfiguration(modelId, parameters, &configuration);
  if (!built.ok()) return built;
  return SetActiveConfiguration(configuration);
}

Status ConfigureActiveModelFromSection(const std::vector<SectionLine>& lines,
                                       ParsedSection* result) {
  ParsedSection parsed;
  const Status status = ParseSection(lines, &parsed);
  if (!status.ok()) return status;
  const Status installed = SetActiveConfiguration(parsed.configuration);
  if (!installed.ok()) return installed;
  if (result) *result = parsed;
  return Status::Success();
}

ModelConfiguration GetActiveConfiguration() { return gActiveConfiguration; }
ShapeConstants GetActiveShapeConstants() { return gActiveConstants; }
bool HasActiveConfiguration() { return gHasActiveConfiguration; }

}  // namespace PitchAngleDiffusion
}  // namespace SEP
