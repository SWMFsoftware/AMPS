#include "mean_free_path.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

namespace SEP {
namespace MeanFreePath {
namespace {

constexpr double Pi = 3.141592653589793238462643383279502884;

std::string Lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

bool FinitePositive(double x) { return std::isfinite(x) && x > 0.0; }

Status Missing(const std::string& name) {
  return Status::Error(StatusCode::MissingInput,
                       "missing required input '" + name + "'");
}

Status Invalid(const std::string& text) {
  return Status::Error(StatusCode::InvalidInput, text);
}

bool ParsePlainDouble(const std::string& text, double* value) {
  if (!value || text.empty()) return false;
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno != 0 || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParseDimensionless(const std::string& text, double* value) {
  const std::size_t slash = text.find('/');
  if (slash == std::string::npos) return ParsePlainDouble(text, value);
  if (slash == 0 || slash + 1 == text.size() ||
      text.find('/', slash + 1) != std::string::npos) return false;
  auto integerLiteral = [](const std::string& item) {
    std::size_t i = (!item.empty() && (item[0] == '+' || item[0] == '-')) ? 1 : 0;
    if (i == item.size()) return false;
    for (; i < item.size(); ++i)
      if (!std::isdigit(static_cast<unsigned char>(item[i]))) return false;
    return true;
  };
  if (!integerLiteral(text.substr(0, slash)) ||
      !integerLiteral(text.substr(slash + 1))) return false;
  double numerator = 0.0, denominator = 0.0;
  if (!ParsePlainDouble(text.substr(0, slash), &numerator) ||
      !ParsePlainDouble(text.substr(slash + 1), &denominator) ||
      denominator == 0.0) return false;
  *value = numerator / denominator;
  return std::isfinite(*value);
}

bool IsOneOf(const std::string& id,
             std::initializer_list<const char*> values) {
  const std::string lower = Lower(id);
  for (const char* value : values)
    if (lower == Lower(value)) return true;
  return false;
}

const std::vector<ModelDescriptor> kRegistry = {
  // Empirical SEP prescriptions.  Their numeric values are not defaults in
  // this library: BuildConfiguration requires the selected, source-auditable
  // value explicitly, including values selected from a published range.
  {"SEP-PATH09", RuntimeState::RequiresUserDecision, true,
   "Verkhoglyadova2009", "Eq. (2)", "tagged lambda", "U-11: lambda kind must be explicit"},
  {"SEP-EPREM10", RuntimeState::ReadyExplicitInputs, true,
   "Schwadron2010", "Sect. 3; Eq. (2)", "lambda_parallel", "only r=1 AU is source-audited"},
  {"SEP-EPREM13", RuntimeState::ReadyExplicitInputs, true,
   "Kozarev2013", "Eqs. (2)-(3)", "lambda_parallel", "D-2 operator convention is separate"},
  {"SEP-MFLAMPA19", RuntimeState::ReadyExplicitInputs, true,
   "Borovikov2019", "Eq. (6.7)", "lambda_parallel", "lambda0 is selected from a published range"},
  {"SEP-MFLAMPA25", RuntimeState::ReadyExplicitInputs, true,
   "Liu2025", "Eqs. (14)-(16)", "lambda_parallel", "Eq. (9) power law"},
  {"SEP-SOFIE24", RuntimeState::ReadyExplicitInputs, true,
   "Zhao2024", "Sect. 4.5", "lambda_parallel", "constant upstream value"},
  {"SEP-ZHANG23", RuntimeState::ReadyExplicitInputs, true,
   "Zhang2023", "Eqs. (35),(37)", "lambda_radial_sep", "parallel conversion requires psi"},
  {"SEP-PARASOL25", RuntimeState::ReadyExplicitInputs, true,
   "Afanasiev2025", "Eqs. (41)-(44)", "lambda_parallel", "background law only"},
  {"SEP-SOLPENCO05", RuntimeState::RequiresUserDecision, true,
   "Aran2005", "input variable [4]", "tagged lambda", "U-11: lambda kind must be explicit"},
  {"SEP-SPARX15", RuntimeState::ReadyExplicitInputs, true,
   "Marsh2015", "Sect. 3.2", "lambda_isotropic", "Poisson direction reset"},
  {"SEP-MARSH13", RuntimeState::ReadyExplicitInputs, true,
   "Marsh2013", "Sect. 2.2", "lambda_isotropic", "Poisson direction reset"},
  {"SEP-HE11", RuntimeState::ReadyExplicitInputs, true,
   "He2011", "Eqs. (4)-(6)", "lambda_radial_sep", "source equation was reconstructed"},
  {"SEP-WANGQIN15", RuntimeState::ReadyExplicitInputs, true,
   "WangQin2015", "Eqs. (3)-(4)", "lambda_parallel", "Eq. (9) power law"},
  {"SEP-KUBO15", RuntimeState::ReadyExplicitInputs, true,
   "Kubo2015", "Eqs. (6),(8)", "lambda_radial_sep", "source writes lambda_rr"},
  {"SEP-PARADISE19", RuntimeState::ReadyExplicitInputs, true,
   "Wijsen2019", "Eqs. (8),(11)", "lambda_radial_sep", "epsilon shape is separately configured"},
  {"SEP-STRAUSS17G", RuntimeState::ReadyExplicitInputs, true,
   "Strauss2017GLE", "Eqs. (10)-(11)", "lambda_radial_sep", "lambda0 must be supplied"},
  {"SEP-STRAUSS15", RuntimeState::ReadyExplicitInputs, true,
   "StraussFichtner2015", "Sect. III-IV", "lambda_parallel", "published value selection required"},
  {"SEP-DROGE16P", RuntimeState::ReadyExplicitInputs, true,
   "Droge2016ICRC", "Eqs. (3.3)-(3.4)", "lambda_radial_sep", "event/sector value must be explicit"},
  {"SEP-LAITINEN16", RuntimeState::BlockedDependency, false,
   "Laitinen2016", "Eqs. (2),(4),(6)", "lambda_parallel", "requires a configured PARALLEL QLT provider"},
  {"SEP-LAITINEN18", RuntimeState::BlockedDependency, false,
   "Laitinen2018", "text", "lambda_parallel", "requires a configured PARALLEL QLT provider"},
  {"SEP-MINOSHIMA26", RuntimeState::RequiresSourceOrCodeAudit, false,
   "Minoshima2026", "Eq. (9)", "lambda_parallel", "D-29: Omega_n is undefined"},
  {"SEP-CHEN24", RuntimeState::ReadyPublishedKappaOnly, true,
   "Chen2024", "Eq. (4)", "kappa_parallel", "lambda is derived only for explicit species"},
  {"LEGACY-TENISHEV2005AIAA", RuntimeState::ReadyExplicitInputs, true,
   "AMPS-srcSEP", "coefficient_providers.cpp", "lambda_parallel", "code-derived legacy law; not a literature attribution"},

  // Pitch-angle and focusing layer.
  {"PA-QFORM", RuntimeState::ReadyExplicitInputs, true,
   "AguedaVainio2013", "Eqs. (12)-(13)", "D_mumu and lambda_parallel", "standard (1-mu^2)"},
  {"PA-QFORM-LANG-PRINTED", RuntimeState::ReadyPublishedVariant, true,
   "Lang2024", "Eq. (3) as printed", "D_mumu", "D-3; explicit no-flux treatment at mu=-1"},
  {"PA-EPS", RuntimeState::ReadyExplicitInputs, true,
   "Agueda2010", "Eq. (14)", "D_mumu and lambda_parallel", "exact phi(epsilon)"},
  {"PA-EPS-PACHECO-PRINTED", RuntimeState::ReadyPublishedVariant, true,
   "Pacheco2019", "Eq. (4) as printed", "D_mumu", "D-4 printed variant"},
  {"PA-DROGE-VA", RuntimeState::ReadyExplicitInputs, true,
   "Minoshima2026", "Eq. (16)", "D_mumu", "nominal and integral lambda both reported"},
  {"PA-ISO", RuntimeState::ReadyExplicitInputs, true,
   "standard", "Eq. (11) isotropic limit", "D_mumu and lambda_parallel", "standard operator"},
  {"PA-EPREM", RuntimeState::ReadyPublishedVariant, true,
   "Kozarev2013", "Eq. (2)", "D_mumu", "D-2 HalfD choice required"},
  {"PA-KOLMO", RuntimeState::ReadyExplicitInputs, true,
   "Borovikov2019", "Eq. (17)", "D_mumu", "q=5/3, H=0 shape"},
  {"PA-AMPS-I", RuntimeState::RequiresSourceOrCodeAudit, false,
   "Tenishev2022", "Table 1, Type I", "D_mumu", "OCR ambiguity and D-7 normalization prevent source-exact execution"},
  {"PA-AMPS-II", RuntimeState::RequiresSourceOrCodeAudit, false,
   "Tenishev2022", "Table 1, Type II", "D_mumu", "printed extraction is insufficient for a dimensionally closed evaluator"},
  {"PA-AMPS-III", RuntimeState::BlockedDependency, false,
   "Tenishev2022", "Table 1, Type III", "D_mumu", "requires the pinned normalized spectrum provider"},
  {"PA-AMPS-IV", RuntimeState::RequiresSourceOrCodeAudit, false,
   "Tenishev2022", "Table 1, Type IV", "D_mumu", "momentum normalization and amplitude are not fixed by the table"},
  {"PA-AMPS-V", RuntimeState::BlockedDependency, false,
   "Tenishev2022/QinWang2015", "Table 1 / Eqs. (2)-(3)", "D_mumu", "use the pinned QLT spectrum closure; do not use the dimensionally inconsistent printing D-6"},
  {"MFP-FOCUS-HW13", RuntimeState::ReadyExplicitInputs, true,
   "HeWan2013", "Eq. (18); stable Eq. (43)", "lambda_parallel", "explicit focusing correction"},

  // Turbulence and shock closures.
  {"QLT-TS03-P", RuntimeState::ReadyExplicitInputs, true,
   "TS2003", "Eq. (21)", "lambda_parallel", "named slab inputs required"},
  {"QLT-TS03-E-RS", RuntimeState::ReadyExplicitInputs, true,
   "TS2003", "Eqs. (22)-(23), RS", "lambda_parallel", "named slab/dissipation inputs required"},
  {"QLT-TS03-E-DT", RuntimeState::ReadyPublishedVariant, true,
   "TS2003/EB2013b/Lang2024", "Eqs. (22)-(23)", "lambda_parallel", "D-1 variant required"},
  {"QLT-ZANK98", RuntimeState::ReadyExplicitInputs, true,
   "Zank1998", "Eqs. (24)-(25)", "lambda_parallel", "U-12 variance convention required"},
  {"SHOCK-BOHM", RuntimeState::ReadyExplicitInputs, true,
   "Zhang2023", "Eq. (26)", "lambda_parallel", "shock-side tag required"},
  {"SHOCK-AFANASIEV15", RuntimeState::ReadyExplicitInputs, true,
   "Afanasiev2015", "Eq. (27)", "lambda_parallel", "upstream signed coordinate"},
  {"SHOCK-PARASOL", RuntimeState::RequiresSourceOrCodeAudit, false,
   "Afanasiev2025", "Eqs. (28)-(29)", "lambda_parallel", "U-6: Lambda(E), Delta x(E) absent"},
  {"SHOCK-MFLAMPA", RuntimeState::ReadyExplicitInputs, true,
   "Liu2025", "Eqs. (30)-(31)", "lambda_parallel and kappa floor", "downstream state and deltaB required"},

  // GCR forms.  Every normalization is a VALUE in SI at this API; a UNIT-only
  // number from a paper is rejected by the source scalar parser until selected.
  {"GCR-NWU14", RuntimeState::ReadyExplicitInputs, true,
   "Potgieter2014", "Eq. (33)", "kappa_parallel", "explicit field normalization required"},
  {"GCR-CORTI19", RuntimeState::ReadyExplicitInputs, true,
   "Corti2019", "Eq. (34)", "kappa_parallel", "explicit parameter row required"},
  {"GCR-LUO19", RuntimeState::ReferenceDataOnly, false,
   "Luo2019", "Eq. (14) not legible", "reference only", "formula unavailable"},
  {"GCR-HELMOD17", RuntimeState::ReadyExplicitInputs, true,
   "Boschini2018ASR", "Eqs. (36)-(37)", "kappa_parallel", "D-21 numerical-use convention required"},
  {"GCR-HELMOD19", RuntimeState::ReadyExplicitInputs, true,
   "Boschini2019", "Eq. (36), v4 radial term", "kappa_parallel", "D-21 numerical-use convention required"},
  {"GCR-BOBIK12", RuntimeState::RequiresSourceOrCodeAudit, false,
   "Bobik2012", "Eqs. (6),(13)", "kappa_parallel", "normalization units absent"},
  {"GCR-STRAUSS11", RuntimeState::ReadyExplicitInputs, true,
   "Strauss2011", "Eq. (4)", "lambda_parallel", "broken linear rigidity law"},
  {"GCR-EFFENBERGER12", RuntimeState::ReadyExplicitInputs, true,
   "Effenberger2012", "Eqs. (22)-(24)", "kappa_parallel", "pc rather than rigidity axis"},
  {"GCR-WANG19", RuntimeState::RequiresExternalInput, true,
   "Wang2019", "Eq. (4)", "kappa_parallel", "B_c(t) must be supplied"},
  {"GCR-TOMASSETTI17", RuntimeState::ReadyExplicitInputs, true,
   "Tomassetti2017", "Sect. IV", "kappa_parallel", "phi(t) and station fit required"},
  {"GCR-PERUGIA21", RuntimeState::RequiresExternalInput, true,
   "Fiandrini2021", "Eq. (8)", "kappa_parallel", "time-dependent parameters must be supplied"},
  {"GCR-PERUGIA25", RuntimeState::RequiresExternalInput, true,
   "Tomassetti2025", "Eq. (5)", "kappa_parallel", "time-dependent parameters must be supplied"},
  {"GCR-JIANG23", RuntimeState::RequiresUserDecision, false,
   "Jiang2023", "Eq. (8)", "kappa_parallel", "printed K0 unit unresolved"},
  {"GCR-DUAN25", RuntimeState::ReadyPublishedVariant, true,
   "Duan2025", "Eqs. (3.5)-(3.6) as printed", "kappa_parallel", "D-33 parameters and variant explicit"},
  {"GCR-QINSHEN17", RuntimeState::BlockedDependency, false,
   "QinShen2017", "NLGCE-F", "lambda_parallel", "requires configured PARALLEL NLGCE-F"},
  {"GCR-EB13", RuntimeState::ReadyExplicitInputs, true,
   "EB2013a/EB2014", "Eq. (21)", "lambda_parallel", "uses QLT-TS03-P backend"}
};

const ModelDescriptor* Descriptor(const std::string& id) {
  const std::string lower = Lower(id);
  for (const ModelDescriptor& item : kRegistry)
    if (lower == Lower(item.stableId)) return &item;
  return nullptr;
}

Status GateStatus(const ModelDescriptor& model) {
  StatusCode code = StatusCode::UnsupportedModel;
  switch (model.declaredState) {
    case RuntimeState::RequiresExternalInput: code = StatusCode::RequiresExternalInput; break;
    case RuntimeState::RequiresUserDecision: code = StatusCode::RequiresUserDecision; break;
    case RuntimeState::RequiresSourceOrCodeAudit: code = StatusCode::RequiresSourceOrCodeAudit; break;
    case RuntimeState::ReferenceDataOnly: code = StatusCode::ReferenceDataOnly; break;
    case RuntimeState::BlockedDependency: code = StatusCode::BlockedDependency; break;
    default: code = StatusCode::UnsupportedModel; break;
  }
  return Status::Error(code, std::string(model.stableId) + ": " + model.note);
}

bool IsSepPowerLaw(const std::string& id) {
  return IsOneOf(id, {"SEP-PATH09", "SEP-EPREM10", "SEP-EPREM13",
      "SEP-MFLAMPA19", "SEP-MFLAMPA25", "SEP-SOFIE24", "SEP-ZHANG23",
      "SEP-PARASOL25", "SEP-SOLPENCO05", "SEP-SPARX15", "SEP-MARSH13",
      "SEP-HE11", "SEP-WANGQIN15", "SEP-KUBO15", "SEP-PARADISE19",
      "SEP-STRAUSS17G", "SEP-STRAUSS15", "SEP-DROGE16P"});
}

bool IsPitch(const std::string& id) {
  return IsOneOf(id, {"PA-QFORM", "PA-QFORM-LANG-PRINTED", "PA-EPS",
      "PA-EPS-PACHECO-PRINTED", "PA-DROGE-VA", "PA-ISO", "PA-EPREM",
      "PA-KOLMO"});
}

Status ReadNumber(const Configuration& c, const std::string& key,
                  double* value) {
  const auto found = c.numbers.find(key);
  if (found == c.numbers.end()) return Missing(key);
  *value = found->second;
  return Status::Success();
}

Status ReadPositive(const Configuration& c, const std::string& key,
                    double* value) {
  Status status = ReadNumber(c, key, value);
  if (!status.ok()) return status;
  return FinitePositive(*value) ? Status::Success()
      : Invalid("'" + key + "' must be positive and finite");
}

Status RequireLocal(const std::optional<double>& value, const char* name,
                    double* output, bool positive = true) {
  if (!value) return Missing(name);
  if (!std::isfinite(*value) || (positive && *value <= 0.0))
    return Invalid(std::string("local '") + name +
                   (positive ? "' must be positive and finite" : "' must be finite"));
  *output = *value;
  return Status::Success();
}

LambdaKind ParseLambdaKind(const std::string& text, bool* ok) {
  const std::string value = Lower(text);
  *ok = true;
  if (value == "parallel") return LambdaKind::Parallel;
  if (value == "radial_sep") return LambdaKind::RadialSEP;
  if (value == "radial_tensor") return LambdaKind::RadialTensor;
  if (value == "isotropic_scattering") return LambdaKind::IsotropicScattering;
  if (value == "unspecified") return LambdaKind::Unspecified;
  *ok = false;
  return LambdaKind::Unspecified;
}

MomentumVariable ParseMomentumVariable(const std::string& text, bool* ok) {
  const std::string value = Lower(text);
  *ok = true;
  if (value == "rigidity_v") return MomentumVariable::RigidityV;
  if (value == "momentum_pc_ev") return MomentumVariable::MomentumPcEV;
  if (value == "kinetic_total_ev") return MomentumVariable::KineticTotalEV;
  if (value == "kinetic_per_nucleon_ev") return MomentumVariable::KineticPerNucleonEV;
  *ok = false;
  return MomentumVariable::RigidityV;
}

Status CheckKeys(const std::map<std::string, std::string>& raw,
                 const std::set<std::string>& required,
                 const std::set<std::string>& optional) {
  for (const std::string& key : required)
    if (raw.find(key) == raw.end()) return Missing(key);
  for (const auto& entry : raw)
    if (!required.count(entry.first) && !optional.count(entry.first))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "unknown parameter '" + entry.first + "'");
  return Status::Success();
}

Status ParseNumbers(const std::map<std::string, std::string>& raw,
                    const std::set<std::string>& keys,
                    const std::set<std::string>& dimensionless,
                    Configuration* c) {
  for (const std::string& key : keys) {
    auto found = raw.find(key);
    if (found == raw.end()) continue;
    double value = 0.0;
    const bool ok = dimensionless.count(key)
        ? ParseDimensionless(found->second, &value)
        : ParsePlainDouble(found->second, &value);
    if (!ok)
      return Status::Error(StatusCode::InvalidConfiguration,
                           "parameter '" + key + "' is not valid suffix-free " +
                           (dimensionless.count(key) ? "dimensionless" : "SI") + " text");
    c->numbers[key] = value;
  }
  return Status::Success();
}

double IndependentValue(const Kinematics& k, MomentumVariable variable,
                        Status* status) {
  switch (variable) {
    case MomentumVariable::RigidityV: return k.rigidityV;
    case MomentumVariable::MomentumPcEV: return k.momentumPcEV;
    case MomentumVariable::KineticTotalEV: return k.kineticTotalEV;
    case MomentumVariable::KineticPerNucleonEV:
      if (!k.kineticPerNucleonEV) {
        *status = Missing("particle.nucleonCount");
        return 0.0;
      }
      return *k.kineticPerNucleonEV;
  }
  *status = Invalid("unknown momentum variable");
  return 0.0;
}

double Simpson(const std::function<double(double)>& f, double a, double b,
               double fa, double fm, double fb, double whole,
               double tolerance, int depth, bool* ok) {
  const double m = 0.5 * (a + b);
  const double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
  const double flm = f(lm), frm = f(rm);
  if (!std::isfinite(flm) || !std::isfinite(frm)) { *ok = false; return 0.0; }
  const double left = (m - a) * (fa + 4.0 * flm + fm) / 6.0;
  const double right = (b - m) * (fm + 4.0 * frm + fb) / 6.0;
  const double sum = left + right;
  if (depth <= 0 || std::abs(sum - whole) <= 15.0 * tolerance)
    return sum + (sum - whole) / 15.0;
  return Simpson(f, a, m, fa, flm, fm, left, tolerance * 0.5, depth - 1, ok) +
         Simpson(f, m, b, fm, frm, fb, right, tolerance * 0.5, depth - 1, ok);
}

Status Integrate01(const std::function<double(double)>& f, double* value) {
  // Gauss/Simpson endpoints are evaluated only for functions whose limiting
  // values are finite.  Singular q-form integrands use their closed H=0 form.
  const double fa = f(0.0), fm = f(0.5), fb = f(1.0);
  if (!std::isfinite(fa) || !std::isfinite(fm) || !std::isfinite(fb))
    return Status::Error(StatusCode::NumericalFailure,
                         "pitch-angle normalization integrand is nonfinite");
  const double whole = (fa + 4.0 * fm + fb) / 6.0;
  bool ok = true;
  *value = Simpson(f, 0.0, 1.0, fa, fm, fb, whole, 1.0e-12, 24, &ok);
  if (!ok || !FinitePositive(*value))
    return Status::Error(StatusCode::NumericalFailure,
                         "pitch-angle normalization quadrature failed");
  return Status::Success();
}

Status PitchShapeIntegral(const std::string& id, double q, double gap,
                          double vaOverV, double* fullIntegral) {
  const std::string model = Lower(id);
  double half = 0.0;
  if (model == "pa-qform" || model == "pa-kolmo") {
    if (!(q > 1.0 && q < 2.0) || !std::isfinite(gap) || gap < 0.0)
      return Invalid("q-form requires 1 < q < 2 and H >= 0");
    if (gap == 0.0) half = 2.0 / ((2.0 - q) * (4.0 - q));
    else {
      Status s = Integrate01([&](double mu) {
        return (1.0 - mu * mu) / (std::pow(mu, q - 1.0) + gap);
      }, &half);
      if (!s.ok()) return s;
    }
    *fullIntegral = 2.0 * half;
    return Status::Success();
  }
  if (model == "pa-eps") {
    if (!FinitePositive(gap)) return Invalid("epsilon must be positive");
    Status s = Integrate01([&](double mu) {
      return 2.0 * (1.0 - mu * mu) * (1.0 + mu) /
             (mu + gap * (1.0 + mu));
    }, &half);
    if (!s.ok()) return s;
    *fullIntegral = 2.0 * half;
    return Status::Success();
  }
  if (model == "pa-iso" || model == "pa-eprem") {
    // D=(nu/2)(1-mu^2), hence integral of numerator/(D/nu) is 8/3.
    *fullIntegral = 8.0 / 3.0;
    return Status::Success();
  }
  if (model == "pa-droge-va") {
    if (!(q > 1.0 && q < 2.0) || !std::isfinite(vaOverV) || vaOverV < 0.0)
      return Invalid("Droge form requires 1 < q < 2 and V_A/v >= 0");
    if (vaOverV == 0.0)
      half = 2.0 / ((2.0 - q) * (4.0 - q));
    else {
      Status s = Integrate01([&](double mu) {
        return (1.0 - mu * mu) /
            std::pow(mu * mu + vaOverV * vaOverV, 0.5 * (q - 1.0));
      }, &half);
      if (!s.ok()) return s;
    }
    *fullIntegral = 2.0 * half;
    return Status::Success();
  }
  if (model == "pa-qform-lang-printed") {
    if (!(q > 1.0 && q < 2.0) || gap < 0.0)
      return Invalid("printed q-form requires 1 < q < 2 and H >= 0");
    // The printed shape is asymmetric; integrate the exact Eq. (11) ratio.
    auto f = [&](double mu) {
      const double denominator = std::pow(std::abs(mu), q - 1.0) + gap;
      return ((1.0 + mu) * (1.0 + mu)) / denominator;
    };
    if (gap == 0.0) {
      // Split away from the integrable mu=0 singularity by analytic power
      // integration of (1 +/- x)^2/x^(q-1).
      *fullIntegral = 2.0 / (2.0 - q) + 2.0 / (4.0 - q);
      return Status::Success();
    }
    double positive = 0.0, negative = 0.0;
    Status s = Integrate01(f, &positive);
    if (!s.ok()) return s;
    s = Integrate01([&](double x) { return f(-x); }, &negative);
    if (!s.ok()) return s;
    *fullIntegral = positive + negative;
    return Status::Success();
  }
  if (model == "pa-eps-pacheco-printed") {
    if (!FinitePositive(gap)) return Invalid("epsilon must be positive");
    // D=(nu/2)[|mu|(1+|mu|)+eps(1-mu^2)] after cancellation.
    Status s = Integrate01([&](double mu) {
      const double dOverNu = 0.5 *
          (mu * (1.0 + mu) + gap * (1.0 - mu * mu));
      return (1.0 - mu * mu) * (1.0 - mu * mu) / dOverNu;
    }, &half);
    if (!s.ok()) return s;
    *fullIntegral = 2.0 * half;
    return Status::Success();
  }
  return Status::Error(StatusCode::UnsupportedModel,
                       "model has no pitch-angle shape: " + id);
}

double HypergeometricDT(double p, double argument) {
  // For Eq. (23), 2F1(1,b;b+1;z) = integral_0^1 du /
  // [1-z*u^(1/b)] after t=u^(1/b).  This removes the t^(b-1)
  // endpoint singularity and is stable for the required negative z.
  const double b = 1.0 / (p - 1.0);
  double value = 0.0;
  Status status = Integrate01([&](double u) {
    return 1.0 / (1.0 - argument * std::pow(u, 1.0 / b));
  }, &value);
  return status.ok() ? value : std::numeric_limits<double>::quiet_NaN();
}

double LogAddExp(double a, double b) {
  const double upper = std::max(a, b), lower = std::min(a, b);
  return upper + std::log1p(std::exp(lower - upper));
}

double Softplus(double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

void AttachProvenance(const Configuration& c, const LocalState& local,
                      Result* result) {
  const ModelDescriptor* d = Descriptor(c.modelId);
  result->provenance.requestedModelId = c.modelId;
  result->provenance.evaluatedModelId = d ? d->stableId : c.modelId;
  result->provenance.sourceKey = d ? d->sourceKey : std::string();
  result->provenance.sourceLocation = d ? d->sourceLocation : std::string();
  result->provenance.formula = d
      ? std::string(d->directOutput) + "; " + d->sourceLocation
      : std::string();
  result->provenance.configurationFingerprint = c.fingerprint;
  result->provenance.sampleIdentity = local.sampleIdentity;
  result->provenance.backgroundRevision = local.backgroundRevision;
  result->provenance.turbulenceRevision = local.turbulenceRevision;
  result->provenance.runtimeState = c.runtimeState;
}

Result Failure(const Configuration& c, const LocalState& local, Status status) {
  Result result;
  result.status = std::move(status);
  AttachProvenance(c, local, &result);
  return result;
}

Result FromLambda(const Configuration& c, const LocalState& local,
                  const Kinematics& k, double lambdaM, LambdaKind kind) {
  if (!FinitePositive(lambdaM))
    return Failure(c, local, Invalid("mean free path is not positive and finite"));
  Result result;
  result.status = Status::Success();
  result.lambda = LambdaValue{lambdaM, kind, local.parkerSpiralAngleRad};
  if (kind == LambdaKind::Parallel)
    result.kappaParallelM2PerS = k.speedMPerS * lambdaM / 3.0;
  AttachProvenance(c, local, &result);
  return result;
}

Result EvaluateSepPower(const ParticleState& particle, const LocalState& local,
                        const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double lambda0 = 0.0, x0 = 0.0, a = 0.0, r0 = 0.0, b = 0.0;
  if (!(status = ReadPositive(c, "lambda0_m", &lambda0)).ok() ||
      !(status = ReadPositive(c, "x0_SI", &x0)).ok() ||
      !(status = ReadNumber(c, "momentum_exponent", &a)).ok() ||
      !(status = ReadPositive(c, "radius0_m", &r0)).ok() ||
      !(status = ReadNumber(c, "radial_exponent", &b)).ok())
    return Failure(c, local, status);
  double radius = 0.0;
  if (!(status = RequireLocal(local.radiusM, "radiusM", &radius)).ok())
    return Failure(c, local, status);
  Status variableStatus = Status::Success();
  const double x = IndependentValue(k, c.momentumVariable, &variableStatus);
  if (!variableStatus.ok()) return Failure(c, local, variableStatus);

  // SEP-EPREM10 is only source-audited at 1 AU because the sign of its radial
  // exponent is illegible.  Refuse any non-unit radial factor even if a caller
  // attempts to supply one through the generic Eq. (9) schema.
  if (IsOneOf(c.modelId, {"SEP-EPREM10"}) &&
      (b != 0.0 || std::abs(radius / AstronomicalUnitM - 1.0) > 1.0e-12))
    return Failure(c, local, Status::Error(StatusCode::RequiresSourceOrCodeAudit,
        "SEP-EPREM10 away from 1 AU is blocked: radial exponent sign is unreadable"));

  const double logLambda = std::log(lambda0) + a * std::log(x / x0) +
                                   b * std::log(radius / r0);
  if (!std::isfinite(logLambda) ||
      logLambda > std::log(std::numeric_limits<double>::max()) ||
      logLambda < std::log(std::numeric_limits<double>::min()))
    return Failure(c, local, Status::Error(StatusCode::OutsideModelDomain,
                                           "SEP power law over/underflows"));
  return FromLambda(c, local, k, std::exp(logLambda), c.outputKind);
}

Result EvaluateChen(const ParticleState& particle, const LocalState& local,
                    const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double radius = 0.0;
  if (!(status = RequireLocal(local.radiusM, "radiusM", &radius)).ok())
    return Failure(c, local, status);
  const double radiusAU = radius / AstronomicalUnitM;
  const double energyKeV = k.kineticTotalEV / 1.0e3;
  const bool inDomain = radiusAU >= 0.1 && radiusAU <= 0.8 &&
                        energyKeV >= 100.0 && energyKeV <= 1.0e6;
  if (!inDomain && c.domainPolicy == OutOfDomainPolicy::Error)
    return Failure(c, local, Status::Error(StatusCode::OutsideModelDomain,
        "SEP-CHEN24 requires 0.1<=r/AU<=0.8 and 100<=E/keV<=1e6"));

  // Chen et al. Eq. (4): 5.16e18 cm^2/s = 5.16e14 m^2/s.
  const double kappa = 5.16e14 * std::pow(radiusAU, 1.17) *
                                    std::pow(energyKeV, 0.71);
  Result result;
  result.status = Status::Success();
  result.kappaParallelM2PerS = kappa;
  if (!inDomain) result.diagnostics |= Extrapolated;
  if (c.outputKind == LambdaKind::Parallel) {
    const auto identity = c.choices.find("species_id");
    if (identity == c.choices.end() || particle.speciesId.empty() ||
        identity->second != particle.speciesId)
      return Failure(c, local, Invalid(
          "SEP-CHEN24 lambda output requires matching explicit species_id"));
    result.lambda = LambdaValue{3.0 * kappa / k.speedMPerS,
                                LambdaKind::Parallel,
                                local.parkerSpiralAngleRad};
    result.diagnostics |= DerivedQuantity;
    result.provenance.derived = true;
  }
  AttachProvenance(c, local, &result);
  result.provenance.derived = result.lambda.has_value();
  return result;
}

Result EvaluateLegacyTenishev(const ParticleState& particle,
                              const LocalState& local,
                              const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double lambda0, energy0, alpha, radius0, beta, radius;
  if (!(status = ReadPositive(c, "lambda0_m", &lambda0)).ok() ||
      !(status = ReadPositive(c, "kinetic_energy0_J", &energy0)).ok() ||
      !(status = ReadNumber(c, "energy_exponent", &alpha)).ok() ||
      !(status = ReadPositive(c, "radius0_m", &radius0)).ok() ||
      !(status = ReadNumber(c, "radial_exponent", &beta)).ok() ||
      !(status = RequireLocal(local.radiusM, "radiusM", &radius)).ok())
    return Failure(c, local, status);
  const double energyJ = k.kineticTotalEV * ElectronVoltJ;
  const double lambda = lambda0 * std::pow(energyJ / energy0, alpha) *
                                  std::pow(radius / radius0, beta);
  return FromLambda(c, local, k, lambda, LambdaKind::Parallel);
}

Result EvaluatePitch(const ParticleState& particle, const LocalState& local,
                     const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double lambda = 0.0, q = 5.0 / 3.0, gap = 0.0, vaRatio = 0.0;
  if (!(status = ReadPositive(c, "lambda_parallel_m", &lambda)).ok())
    return Failure(c, local, status);
  auto n = c.numbers.find("q"); if (n != c.numbers.end()) q = n->second;
  n = c.numbers.find("gap_parameter"); if (n != c.numbers.end()) gap = n->second;
  n = c.numbers.find("alfven_to_particle_speed");
  if (n != c.numbers.end()) vaRatio = n->second;
  LambdaValue target{lambda, LambdaKind::Parallel, local.parkerSpiralAngleRad};
  double amplitude = 0.0;
  if (IsOneOf(c.modelId, {"PA-DROGE-VA"})) {
    // Equation (16) defines its prefactor with a *nominal* lambda.  It is not
    // normalized by Eq. (11) once VA/v is nonzero.  Preserve that distinction
    // instead of silently rescaling the published shape.
    amplitude = 3.0 * k.speedMPerS /
        (2.0 * (4.0 - q) * (2.0 - q) * lambda);
  } else if (IsOneOf(c.modelId, {"PA-EPREM"}) &&
             c.choices.at("lambda_interpretation") == "printed_parameter") {
    // Kozarev Eq. (2) prints D=(1-mu^2)v/(2 lambda).  Therefore nu=v/lambda;
    // the HalfD operator transports with twice the printed lambda (D-2).
    amplitude = k.speedMPerS / lambda;
  } else status = PitchAmplitudeFromLambda(c.modelId, q, gap, vaRatio,
                                             c.operatorConvention, target,
                                             k.speedMPerS, &amplitude);
  if (!status.ok()) return Failure(c, local, status);
  double effectiveLambda = lambda;
  if (IsOneOf(c.modelId, {"PA-DROGE-VA", "PA-EPREM"})) {
    LambdaValue recovered;
    status = LambdaFromPitchAmplitude(c.modelId, q, gap, vaRatio,
                                      c.operatorConvention, amplitude,
                                      k.speedMPerS, &recovered);
    if (!status.ok()) return Failure(c, local, status);
    effectiveLambda = recovered.metres;
  }
  Result result = FromLambda(c, local, k, effectiveLambda, LambdaKind::Parallel);
  result.pitchAmplitudePerS = amplitude;
  if (IsOneOf(c.modelId, {"PA-QFORM-LANG-PRINTED", "PA-EPS-PACHECO-PRINTED",
                          "PA-EPREM"}))
    result.diagnostics |= PublishedVariant;
  if (local.pitchAngleCosine) {
    double d = 0.0;
    status = EvaluatePitchAngleDiffusion(*local.pitchAngleCosine, particle,
                                         local, c, &d);
    if (!status.ok()) return Failure(c, local, status);
    result.dMuMuPerS = d;
  }
  if (IsOneOf(c.modelId, {"PA-DROGE-VA"}) && vaRatio > 0.0)
    result.diagnostics |= NominalAndIntegralLambdaDiffer;
  return result;
}

Result EvaluateFocus(const ParticleState& particle, const LocalState& local,
                     const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double lambda0 = 0.0, length = 0.0, ratio = 0.0;
  if (!(status = ReadPositive(c, "unfocused_lambda_parallel_m", &lambda0)).ok() ||
      !(status = ReadPositive(c, "focusing_length_m", &length)).ok() ||
      !(status = HeWanFocusingRatio(lambda0 / length, &ratio)).ok())
    return Failure(c, local, status);
  Result result = FromLambda(c, local, k, lambda0 * ratio,
                             LambdaKind::Parallel);
  result.diagnostics |= DerivedQuantity;
  result.provenance.derived = true;
  return result;
}

Result EvaluateTS(const ParticleState& particle, const LocalState& local,
                  const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double field, variance, kmin, s;
  if (!(status = RequireLocal(local.meanFieldMagnitudeT,
                              "meanFieldMagnitudeT", &field)).ok() ||
      !(status = RequireLocal(local.slabVarianceT2,
                              "slabVarianceT2", &variance)).ok() ||
      !(status = RequireLocal(local.kMinPerM, "kMinPerM", &kmin)).ok() ||
      !(status = ReadNumber(c, "inertial_index", &s)).ok())
    return Failure(c, local, status);
  if (!(s > 1.0 && s < 2.0))
    return Failure(c, local, Invalid("TS2003 requires 1 < inertial_index < 2"));
  const double rL = k.rigidityV / (SpeedOfLightMPerS * field);
  const double R = rL * kmin;
  double bracket = 1.0 + 8.0 /
      ((2.0 - s) * (4.0 - s) * std::pow(R, s));
  if (IsOneOf(c.modelId, {"QLT-TS03-E-RS", "QLT-TS03-E-DT"})) {
    double p, kd, va, alpha;
    if (!(status = ReadNumber(c, "dissipation_index", &p)).ok() ||
        !(status = RequireLocal(local.kDPerM, "kDPerM", &kd)).ok() ||
        !(status = RequireLocal(local.alfvenSpeedMPerS,
                                "alfvenSpeedMPerS", &va)).ok() ||
        !(status = ReadPositive(c, "alpha_D", &alpha)).ok())
      return Failure(c, local, status);
    if (!(p > 2.0) || !(kd > kmin))
      return Failure(c, local, Invalid("electron TS requires p>2 and kD>kMin"));
    const double Q = rL * kd;
    const double dyn = k.speedMPerS / (alpha * va);
    double K = 0.0, qPower = p - s;
    if (IsOneOf(c.modelId, {"QLT-TS03-E-RS"})) {
      K = (std::sqrt(Pi) / std::tgamma(0.5 * p) + 1.0 / (p - 2.0)) *
          std::pow(0.5 * dyn, p - 2.0);
    } else {
      const std::string variant = c.choices.at("dt_variant");
      const double f1 = 2.0 * (p - s) /
          (Pi * (p - 2.0) * (2.0 - s));
      double argument = -dyn / (f1 * Q);
      if (variant == "ts_q3ms") qPower = 3.0 - s;
      else if (variant == "eb13") argument = -dyn /
          (f1 * std::pow(Q, p - 2.0));
      else if (variant != "lang24")
        return Failure(c, local, Invalid("unknown DT variant"));
      const double hyper = HypergeometricDT(p, argument);
      if (!std::isfinite(hyper))
        return Failure(c, local, Status::Error(StatusCode::NumericalFailure,
                                               "DT hypergeometric evaluation failed"));
      K = hyper * dyn / f1;
    }
    bracket += 4.0 * K / (std::pow(Q, qPower) * std::pow(R, s));
  }
  const double lambda = 3.0 * s * rL * rL * kmin /
      (4.0 * Pi * (s - 1.0)) * (field * field / variance) * bracket;
  return FromLambda(c, local, k, lambda, LambdaKind::Parallel);
}

Result EvaluateZank(const ParticleState& particle, const LocalState& local,
                    const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  double field, variance, length;
  if (!(status = RequireLocal(local.meanFieldMagnitudeT,
                              "meanFieldMagnitudeT", &field)).ok() ||
      !(status = RequireLocal(local.slabVarianceT2,
                              "slabVarianceT2", &variance)).ok() ||
      !(status = RequireLocal(local.slabCorrelationLengthM,
                              "slabCorrelationLengthM", &length)).ok())
    return Failure(c, local, status);
  const std::string convention = c.choices.at("variance_convention");
  if (k.rigidityV < 1.0e7 || k.rigidityV > 1.0e10)
    return Failure(c, local, Status::Error(StatusCode::OutsideModelDomain,
        "QLT-ZANK98 published validity is 10 MV <= rigidity <= 10 GV"));
  const double coefficient = convention == "per_component" ? 3.1371 : 6.2742;
  const double rL = k.rigidityV / (SpeedOfLightMPerS * field);
  const double z = 0.746834 * rL / length;
  const double z2 = z * z;
  const double A = std::expm1((5.0 / 6.0) * std::log1p(z2));
  const double denom = z2 - std::expm1((1.0 / 6.0) * std::log1p(z2));
  const double qaux = z2 < 1.0e-14 ? 2.0 : (5.0 * z2 / 3.0) / denom;
  const double correction = 1.0 + (7.0 * A / 9.0) /
      ((qaux + 1.0 / 3.0) * (qaux + 7.0 / 3.0));
  const double lambda = coefficient * std::pow(field, 5.0 / 3.0) / variance *
      std::pow(k.rigidityV / SpeedOfLightMPerS, 1.0 / 3.0) *
      std::pow(length, 2.0 / 3.0) * correction;
  return FromLambda(c, local, k, lambda, LambdaKind::Parallel);
}

Result EvaluateShock(const ParticleState& particle, const LocalState& local,
                     const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  if (local.shockSide == LocalState::ShockSide::Unspecified)
    return Failure(c, local, Missing("shockSide"));
  if (IsOneOf(c.modelId, {"SHOCK-BOHM"})) {
    double field;
    if (!(status = RequireLocal(local.meanFieldMagnitudeT,
                                "meanFieldMagnitudeT", &field)).ok())
      return Failure(c, local, status);
    const double rg = particle.momentumKgMPerS /
                      (std::abs(particle.chargeC) * field);
    return FromLambda(c, local, k, rg, LambdaKind::Parallel);
  }
  if (IsOneOf(c.modelId, {"SHOCK-AFANASIEV15"})) {
    if (local.shockSide != LocalState::ShockSide::Upstream)
      return Failure(c, local, Invalid("SHOCK-AFANASIEV15 is upstream-only"));
    double x, u, va, x0;
    if (!(status = RequireLocal(local.upstreamDistanceM,
                                "upstreamDistanceM", &x, false)).ok() ||
        !(status = RequireLocal(local.upstreamFlowShockFrameMPerS,
                                "upstreamFlowShockFrameMPerS", &u)).ok() ||
        !(status = RequireLocal(local.alfvenSpeedMPerS,
                                "alfvenSpeedMPerS", &va)).ok() ||
        !(status = ReadNumber(c, "x0_m", &x0)).ok())
      return Failure(c, local, status);
    if (x < 0.0 || x0 < 0.0 || !(u > va) || !(x + x0 > 0.0))
      return Failure(c, local, Invalid(
          "Afanasiev Eq. (27) requires x>=0, x0>=0, u1>VA, x+x0>0"));
    return FromLambda(c, local, k, 3.0 * (u - va) * (x + x0) /
                                     k.speedMPerS, LambdaKind::Parallel);
  }
  if (local.shockSide != LocalState::ShockSide::Downstream)
    return Failure(c, local, Invalid("SHOCK-MFLAMPA is downstream-only"));
  double field, variance, radius, lmaxFactor, floor;
  if (!(status = RequireLocal(local.meanFieldMagnitudeT,
                              "meanFieldMagnitudeT", &field)).ok() ||
      !(status = RequireLocal(local.waveVarianceT2,
                              "waveVarianceT2", &variance)).ok() ||
      !(status = RequireLocal(local.radiusM, "radiusM", &radius)).ok() ||
      !(status = ReadPositive(c, "Lmax_over_radius", &lmaxFactor)).ok() ||
      !(status = ReadPositive(c, "kappa_floor_m2_per_s", &floor)).ok())
    return Failure(c, local, status);
  const double oneGeVMomentum = 1.0e9 * ElectronVoltJ / SpeedOfLightMPerS;
  const double rL0 = oneGeVMomentum / (ElementaryChargeC * field);
  const double k0 = 2.0 * Pi / (lmaxFactor * radius);
  const double lambda = 81.0 / (7.0 * Pi) * field * field / variance *
      std::pow(rL0, 1.0 / 3.0) / std::pow(k0, 2.0 / 3.0) *
      std::pow(k.momentumPcEV / 1.0e9, 1.0 / 3.0);
  Result result = FromLambda(c, local, k, lambda, LambdaKind::Parallel);
  if (!result.status.ok()) return result;
  if (*result.kappaParallelM2PerS < floor) {
    result.kappaParallelM2PerS = floor;
    result.lambda->metres = 3.0 * floor / k.speedMPerS;
    result.diagnostics |= NumericalFloorApplied;
  }
  return result;
}

Result FromKappa(const Configuration& c, const LocalState& local,
                 const Kinematics& k, double kappa) {
  if (!FinitePositive(kappa))
    return Failure(c, local, Invalid("parallel coefficient is not positive and finite"));
  Result result;
  result.status = Status::Success();
  result.kappaParallelM2PerS = kappa;
  result.lambda = LambdaValue{3.0 * kappa / k.speedMPerS,
                              LambdaKind::Parallel,
                              local.parkerSpiralAngleRad};
  result.diagnostics |= DerivedQuantity;
  result.provenance.derived = true;
  AttachProvenance(c, local, &result);
  result.provenance.derived = true;
  return result;
}

Result EvaluateGcr(const ParticleState& particle, const LocalState& local,
                   const Configuration& c) {
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return Failure(c, local, status);
  const double beta = k.beta;
  if (IsOneOf(c.modelId, {"GCR-STRAUSS11"})) {
    double lambda0, p0, r0, radius;
    if (!(status = ReadPositive(c, "lambda0_m", &lambda0)).ok() ||
        !(status = ReadPositive(c, "rigidity0_V", &p0)).ok() ||
        !(status = ReadPositive(c, "radius0_m", &r0)).ok() ||
        !(status = RequireLocal(local.radiusM, "radiusM", &radius)).ok())
      return Failure(c, local, status);
    return FromLambda(c, local, k, lambda0 * std::max(1.0, k.rigidityV / p0) *
                                     (1.0 + radius / r0), LambdaKind::Parallel);
  }
  double field = 0.0;
  const bool helmod = IsOneOf(c.modelId, {"GCR-HELMOD17", "GCR-HELMOD19"});
  if (!helmod &&
      !(status = RequireLocal(local.meanFieldMagnitudeT,
                              "meanFieldMagnitudeT", &field)).ok())
    return Failure(c, local, status);
  double kappa = 0.0;
  if (IsOneOf(c.modelId, {"GCR-NWU14"})) {
    double K0, B0, P0, Pk, a, b, smooth;
    if (!(status = ReadPositive(c, "K0_m2_per_s", &K0)).ok() ||
        !(status = ReadPositive(c, "field_reference_T", &B0)).ok() ||
        !(status = ReadPositive(c, "rigidity0_V", &P0)).ok() ||
        !(status = ReadPositive(c, "break_rigidity_V", &Pk)).ok() ||
        !(status = ReadNumber(c, "low_slope", &a)).ok() ||
        !(status = ReadNumber(c, "high_slope", &b)).ok() ||
        !(status = ReadPositive(c, "smoothness", &smooth)).ok())
      return Failure(c, local, status);
    const double x = k.rigidityV / P0, xb = Pk / P0;
    const double logShape = a * std::log(x) + (b - a) / smooth *
        (LogAddExp(smooth * std::log(x), smooth * std::log(xb)) -
         Softplus(smooth * std::log(xb)));
    kappa = K0 * beta * (B0 / field) * std::exp(logShape);
  } else if (IsOneOf(c.modelId, {"GCR-CORTI19"})) {
    double K0, B0, Pk, a, b, smooth;
    if (!(status = ReadPositive(c, "K0_m2_per_s", &K0)).ok() ||
        !(status = ReadPositive(c, "field_reference_T", &B0)).ok() ||
        !(status = ReadPositive(c, "break_rigidity_V", &Pk)).ok() ||
        !(status = ReadNumber(c, "low_slope", &a)).ok() ||
        !(status = ReadNumber(c, "high_slope", &b)).ok() ||
        !(status = ReadPositive(c, "smoothness", &smooth)).ok())
      return Failure(c, local, status);
    const double x = k.rigidityV / Pk;
    const double logShape = a * std::log(x) +
        (b - a) / smooth * Softplus(smooth * std::log(x));
    kappa = K0 * beta * B0 / field * std::exp(logShape);
  } else if (IsOneOf(c.modelId, {"GCR-HELMOD17", "GCR-HELMOD19"})) {
    double K0, glow, radius;
    if (!(status = ReadPositive(c, "K0_AU2_per_s_numeric", &K0)).ok() ||
        !(status = ReadNumber(c, "g_low", &glow)).ok() ||
        !(status = RequireLocal(local.radiusM, "radiusM", &radius)).ok())
      return Failure(c, local, status);
    if (c.choices.at("normalization_convention") != "d21_numeric_use")
      return Failure(c, local, Invalid("HelMod requires d21_numeric_use"));
    double radial = 1.0 + radius / AstronomicalUnitM;
    if (IsOneOf(c.modelId, {"GCR-HELMOD19"})) {
      double rc;
      if (!(status = ReadNumber(c, "R_c", &rc)).ok()) return Failure(c, local, status);
      radial = rc + radius / AstronomicalUnitM;
    }
    kappa = (beta / 3.0) * K0 *
        (k.rigidityV / 1.0e9 + glow) * radial *
        AstronomicalUnitM * AstronomicalUnitM;
  } else if (IsOneOf(c.modelId, {"GCR-EFFENBERGER12"})) {
    double K0, pc0, B0, exponent;
    if (!(status = ReadPositive(c, "K0_m2_per_s", &K0)).ok() ||
        !(status = ReadPositive(c, "momentum_pc0_eV", &pc0)).ok() ||
        !(status = ReadPositive(c, "field_reference_T", &B0)).ok() ||
        !(status = ReadNumber(c, "field_exponent", &exponent)).ok())
      return Failure(c, local, status);
    kappa = K0 * beta * (k.momentumPcEV / pc0) * std::pow(B0 / field, exponent);
  } else if (IsOneOf(c.modelId, {"GCR-WANG19"})) {
    double scale, Bc, BE, exponent;
    if (!(status = ReadPositive(c, "base_scale_m2_per_s", &scale)).ok() ||
        !(status = ReadPositive(c, "B_c_T", &Bc)).ok() ||
        !(status = ReadPositive(c, "mean_earth_field_T", &BE)).ok() ||
        !(status = ReadNumber(c, "activity_exponent", &exponent)).ok())
      return Failure(c, local, status);
    const double factor = std::pow(Bc / BE, exponent) * beta * BE / field;
    const double rigidityGV = k.rigidityV / 1.0e9;
    kappa = scale * factor * (rigidityGV < 0.1 ? 1.0 / 30.0
                                                : rigidityGV / 3.0);
  } else if (IsOneOf(c.modelId, {"GCR-TOMASSETTI17"})) {
    double fitA, fitB, phiMV, B0;
    if (!(status = ReadNumber(c, "fit_a_MV", &fitA)).ok() ||
        !(status = ReadNumber(c, "fit_b", &fitB)).ok() ||
        !(status = ReadPositive(c, "modulation_potential_MV", &phiMV)).ok() ||
        !(status = ReadPositive(c, "field_reference_T", &B0)).ok())
      return Failure(c, local, status);
    const double kappa0 = fitA / phiMV + fitB;
    kappa = kappa0 * 1.0e18 * beta * (k.rigidityV / 1.0e9) /
            (3.0 * field / B0);  // 10^22 cm^2/s = 10^18 m^2/s.
  } else if (IsOneOf(c.modelId, {"GCR-PERUGIA21", "GCR-PERUGIA25"})) {
    double K0, B0, P0, Pk, a, b, h;
    if (!(status = ReadPositive(c, "K0_m2_per_s", &K0)).ok() ||
        !(status = ReadPositive(c, "field_reference_T", &B0)).ok() ||
        !(status = ReadPositive(c, "rigidity0_V", &P0)).ok() ||
        !(status = ReadPositive(c, "break_rigidity_V", &Pk)).ok() ||
        !(status = ReadNumber(c, "low_slope", &a)).ok() ||
        !(status = ReadNumber(c, "high_slope", &b)).ok() ||
        !(status = ReadPositive(c, "smoothness", &h)).ok())
      return Failure(c, local, status);
    const double x = k.rigidityV / P0, xb = Pk / P0;
    const double logTransition = LogAddExp(h * std::log(x), h * std::log(xb)) -
                                 Softplus(h * std::log(xb));
    const double logShape = a * std::log(x) + (b - a) / h * logTransition;
    kappa = K0 * beta / 3.0 * B0 / field * std::exp(logShape);
  } else if (IsOneOf(c.modelId, {"GCR-DUAN25"})) {
    double K0, Beq, Pk, a, b, smooth;
    if (!(status = ReadPositive(c, "K0_m2_per_s", &K0)).ok() ||
        !(status = ReadPositive(c, "equatorial_field_T", &Beq)).ok() ||
        !(status = ReadPositive(c, "break_rigidity_V", &Pk)).ok() ||
        !(status = ReadNumber(c, "a", &a)).ok() ||
        !(status = ReadNumber(c, "b", &b)).ok() ||
        !(status = ReadPositive(c, "c", &smooth)).ok())
      return Failure(c, local, status);
    const double x = k.rigidityV / Pk;
    const double logShape = a * std::log(x) + smooth *
        Softplus((b - a) / smooth * std::log(x));
    kappa = K0 * beta * Beq / field * std::exp(logShape);
  } else {
    return Failure(c, local, Status::Error(StatusCode::UnsupportedModel,
                                           "GCR backend is unavailable"));
  }
  return FromKappa(c, local, k, kappa);
}

Configuration gActiveConfiguration;
Result Unconfigured(const ParticleState&, const LocalState& local,
                    const Configuration& c) {
  return Failure(c, local, Status::Error(StatusCode::InvalidConfiguration,
                                         "no active mean-free-path model"));
}

}  // namespace

Status Status::Success() { return Status{}; }
Status Status::Error(StatusCode code, const std::string& detail) {
  Status status; status.code = code; status.detail = detail; return status;
}

Status ComputeKinematics(const ParticleState& p, Kinematics* out) {
  if (!out) return Invalid("kinematics output is null");
  if (!FinitePositive(p.massKg) || !FinitePositive(p.momentumKgMPerS) ||
      !std::isfinite(p.chargeC) || p.chargeC == 0.0)
    return Status::Error(StatusCode::InvalidParticle,
        "mass and momentum must be positive; charge must be finite and nonzero");
  const double mc = p.massKg * SpeedOfLightMPerS;
  const double u = p.momentumKgMPerS / mc;
  out->gamma = std::hypot(1.0, u);
  out->beta = u / out->gamma;
  out->speedMPerS = out->beta * SpeedOfLightMPerS;
  out->rigidityV = p.momentumKgMPerS * SpeedOfLightMPerS /
                   std::abs(p.chargeC);
  out->momentumPcEV = p.momentumKgMPerS * SpeedOfLightMPerS / ElectronVoltJ;
  // gamma-1 loses precision for non-relativistic particles; u^2/(gamma+1)
  // is algebraically identical and retains the low-energy bits.
  out->kineticTotalEV = p.massKg * SpeedOfLightMPerS * SpeedOfLightMPerS /
                        ElectronVoltJ * u * u / (out->gamma + 1.0);
  if (p.nucleonCount) {
    if (!FinitePositive(*p.nucleonCount))
      return Status::Error(StatusCode::InvalidParticle,
                           "nucleonCount must be positive when supplied");
    out->kineticPerNucleonEV = out->kineticTotalEV / *p.nucleonCount;
  } else out->kineticPerNucleonEV.reset();
  return Status::Success();
}

const std::vector<ModelDescriptor>& ModelRegistry() { return kRegistry; }
const ModelDescriptor* FindModel(const std::string& id) { return Descriptor(id); }

Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          Configuration* output) {
  if (!output) return Invalid("configuration output is null");
  const ModelDescriptor* d = Descriptor(modelId);
  if (!d) return Status::Error(StatusCode::UnsupportedModel,
                               "unknown model '" + modelId + "'");
  if (!d->executable) return GateStatus(*d);

  Configuration c;
  c.modelId = d->stableId;
  c.runtimeState = d->declaredState;
  for (const InputParameter& p : parameters) {
    if (p.name.empty())
      return Status::Error(StatusCode::InvalidConfiguration,
                           "parameter name must not be empty");
    if (!c.raw.emplace(p.name, p.value).second)
      return Status::Error(StatusCode::InvalidConfiguration,
                           "duplicate parameter '" + p.name + "'");
  }

  std::set<std::string> required, optional, numeric, dimensionless;
  if (IsSepPowerLaw(c.modelId)) {
    required = {"lambda0_m", "lambda_kind", "momentum_variable", "x0_SI",
                "momentum_exponent", "radius0_m", "radial_exponent"};
    numeric = {"lambda0_m", "x0_SI", "momentum_exponent", "radius0_m",
               "radial_exponent"};
    dimensionless = {"momentum_exponent", "radial_exponent"};
  } else if (IsOneOf(c.modelId, {"SEP-CHEN24"})) {
    required = {"output_quantity", "domain_policy"};
    optional = {"species_id"};
  } else if (IsOneOf(c.modelId, {"LEGACY-TENISHEV2005AIAA"})) {
    required = {"lambda0_m", "kinetic_energy0_J", "energy_exponent",
                "radius0_m", "radial_exponent"};
    numeric = required;
    dimensionless = {"energy_exponent", "radial_exponent"};
  } else if (IsPitch(c.modelId)) {
    required = {"lambda_parallel_m", "operator_convention"};
    numeric = {"lambda_parallel_m"};
    if (IsOneOf(c.modelId, {"PA-QFORM", "PA-QFORM-LANG-PRINTED", "PA-KOLMO"})) {
      required.insert("q"); required.insert("gap_parameter");
      numeric.insert("q"); numeric.insert("gap_parameter");
      dimensionless.insert("q"); dimensionless.insert("gap_parameter");
    } else if (IsOneOf(c.modelId, {"PA-EPS", "PA-EPS-PACHECO-PRINTED"})) {
      required.insert("gap_parameter"); numeric.insert("gap_parameter");
      dimensionless.insert("gap_parameter");
    } else if (IsOneOf(c.modelId, {"PA-DROGE-VA"})) {
      required.insert("q"); required.insert("alfven_to_particle_speed");
      numeric.insert("q"); numeric.insert("alfven_to_particle_speed");
      dimensionless.insert("q"); dimensionless.insert("alfven_to_particle_speed");
    }
    if (IsOneOf(c.modelId, {"PA-EPREM"})) required.insert("lambda_interpretation");
  } else if (IsOneOf(c.modelId, {"MFP-FOCUS-HW13"})) {
    required = {"unfocused_lambda_parallel_m", "focusing_length_m"}; numeric = required;
  } else if (IsOneOf(c.modelId, {"QLT-TS03-P", "GCR-EB13"})) {
    required = {"inertial_index"}; numeric = required; dimensionless = required;
  } else if (IsOneOf(c.modelId, {"QLT-TS03-E-RS"})) {
    required = {"inertial_index", "dissipation_index", "alpha_D"};
    numeric = required; dimensionless = required;
  } else if (IsOneOf(c.modelId, {"QLT-TS03-E-DT"})) {
    required = {"inertial_index", "dissipation_index", "alpha_D", "dt_variant"};
    numeric = {"inertial_index", "dissipation_index", "alpha_D"};
    dimensionless = numeric;
  } else if (IsOneOf(c.modelId, {"QLT-ZANK98"})) {
    required = {"variance_convention"};
  } else if (IsOneOf(c.modelId, {"SHOCK-BOHM"})) {
  } else if (IsOneOf(c.modelId, {"SHOCK-AFANASIEV15"})) {
    required = {"x0_m"}; numeric = required;
  } else if (IsOneOf(c.modelId, {"SHOCK-MFLAMPA"})) {
    required = {"Lmax_over_radius", "kappa_floor_m2_per_s"}; numeric = required;
    dimensionless = {"Lmax_over_radius"};
  } else if (IsOneOf(c.modelId, {"GCR-NWU14"})) {
    required = {"K0_m2_per_s", "field_reference_T", "field_normalization",
                "rigidity0_V", "break_rigidity_V", "low_slope",
                "high_slope", "smoothness"};
    numeric = {"K0_m2_per_s", "field_reference_T", "rigidity0_V",
               "break_rigidity_V", "low_slope", "high_slope", "smoothness"};
    dimensionless = {"low_slope", "high_slope", "smoothness"};
  } else if (IsOneOf(c.modelId, {"GCR-CORTI19"})) {
    required = {"K0_m2_per_s", "field_reference_T", "field_normalization",
                "break_rigidity_V", "low_slope", "high_slope", "smoothness"};
    numeric = {"K0_m2_per_s", "field_reference_T", "break_rigidity_V",
               "low_slope", "high_slope", "smoothness"};
    dimensionless = {"low_slope", "high_slope", "smoothness"};
  } else if (IsOneOf(c.modelId, {"GCR-HELMOD17"})) {
    required = {"K0_AU2_per_s_numeric", "g_low", "normalization_convention"};
    numeric = {"K0_AU2_per_s_numeric", "g_low"}; dimensionless = numeric;
  } else if (IsOneOf(c.modelId, {"GCR-HELMOD19"})) {
    required = {"K0_AU2_per_s_numeric", "g_low", "R_c",
                "normalization_convention"};
    numeric = {"K0_AU2_per_s_numeric", "g_low", "R_c"}; dimensionless = numeric;
  } else if (IsOneOf(c.modelId, {"GCR-STRAUSS11"})) {
    required = {"lambda0_m", "rigidity0_V", "radius0_m"}; numeric = required;
  } else if (IsOneOf(c.modelId, {"GCR-EFFENBERGER12"})) {
    required = {"K0_m2_per_s", "momentum_pc0_eV", "field_reference_T",
                "field_normalization", "field_exponent"};
    numeric = {"K0_m2_per_s", "momentum_pc0_eV", "field_reference_T",
               "field_exponent"}; dimensionless = {"field_exponent"};
  } else if (IsOneOf(c.modelId, {"GCR-WANG19"})) {
    required = {"base_scale_m2_per_s", "B_c_T", "mean_earth_field_T",
                "field_normalization", "activity_exponent"};
    numeric = {"base_scale_m2_per_s", "B_c_T", "mean_earth_field_T",
               "activity_exponent"}; dimensionless = {"activity_exponent"};
  } else if (IsOneOf(c.modelId, {"GCR-TOMASSETTI17"})) {
    required = {"fit_a_MV", "fit_b", "modulation_potential_MV",
                "field_reference_T", "field_normalization"};
    numeric = {"fit_a_MV", "fit_b", "modulation_potential_MV",
               "field_reference_T"}; dimensionless = {"fit_a_MV", "fit_b",
                                                      "modulation_potential_MV"};
  } else if (IsOneOf(c.modelId, {"GCR-PERUGIA21", "GCR-PERUGIA25"})) {
    required = {"K0_m2_per_s", "field_reference_T", "field_normalization",
                "rigidity0_V", "break_rigidity_V", "low_slope",
                "high_slope", "smoothness"};
    numeric = {"K0_m2_per_s", "field_reference_T", "rigidity0_V",
               "break_rigidity_V", "low_slope", "high_slope", "smoothness"};
    dimensionless = {"low_slope", "high_slope", "smoothness"};
  } else if (IsOneOf(c.modelId, {"GCR-DUAN25"})) {
    required = {"K0_m2_per_s", "equatorial_field_T", "field_normalization",
                "break_rigidity_V", "a", "b", "c", "formula_variant"};
    numeric = {"K0_m2_per_s", "equatorial_field_T", "break_rigidity_V",
               "a", "b", "c"}; dimensionless = {"a", "b", "c"};
  } else return Status::Error(StatusCode::UnsupportedModel,
                              "no parser schema for '" + c.modelId + "'");

  Status status = CheckKeys(c.raw, required, optional);
  if (!status.ok()) return status;
  status = ParseNumbers(c.raw, numeric, dimensionless, &c);
  if (!status.ok()) return status;
  for (const auto& entry : c.raw)
    if (!numeric.count(entry.first)) c.choices[entry.first] = Lower(entry.second);

  if (IsSepPowerLaw(c.modelId)) {
    bool ok = false;
    c.outputKind = ParseLambdaKind(c.choices["lambda_kind"], &ok);
    if (!ok || c.outputKind == LambdaKind::Unspecified)
      return Status::Error(StatusCode::RequiresUserDecision,
                           "lambda_kind must explicitly resolve U-11");
    if (c.outputKind == LambdaKind::RadialTensor)
      return Status::Error(StatusCode::InvalidConfiguration,
          "SEP Eq. (9) does not define radial-tensor lambda; assemble Eq. (3) explicitly");
    c.momentumVariable = ParseMomentumVariable(c.choices["momentum_variable"], &ok);
    if (!ok) return Invalid("unknown momentum_variable");
  }
  if (IsOneOf(c.modelId, {"SEP-CHEN24"})) {
    const std::string outputChoice = c.choices["output_quantity"];
    if (outputChoice == "kappa_parallel") c.outputKind = LambdaKind::Unspecified;
    else if (outputChoice == "lambda_parallel") {
      c.outputKind = LambdaKind::Parallel;
      if (!c.raw.count("species_id") || c.raw.at("species_id").empty())
        return Status::Error(StatusCode::RequiresUserDecision,
                             "lambda output requires explicit species_id (U-2)");
      c.choices["species_id"] = c.raw.at("species_id");
    } else return Invalid("output_quantity must be kappa_parallel or lambda_parallel");
    const std::string policy = c.choices["domain_policy"];
    if (policy == "error") c.domainPolicy = OutOfDomainPolicy::Error;
    else if (policy == "warn_and_evaluate") c.domainPolicy = OutOfDomainPolicy::WarnAndEvaluate;
    else return Invalid("domain_policy must be error or warn_and_evaluate");
  }
  if (IsPitch(c.modelId)) {
    const std::string op = c.choices["operator_convention"];
    if (op == "standard") c.operatorConvention = OperatorConvention::Standard;
    else if (op == "half_d") c.operatorConvention = OperatorConvention::HalfD;
    else return Invalid("operator_convention must be standard or half_d");
    if (IsOneOf(c.modelId, {"PA-EPREM"}) && c.operatorConvention != OperatorConvention::HalfD)
      return Status::Error(StatusCode::RequiresUserDecision,
                           "PA-EPREM requires explicit half_d convention (U-8)");
    if (IsOneOf(c.modelId, {"PA-EPREM"}) &&
        !IsOneOf(c.choices["lambda_interpretation"],
                 {"printed_parameter", "transport_lambda"}))
      return Status::Error(StatusCode::RequiresUserDecision,
          "PA-EPREM lambda_interpretation must be printed_parameter or transport_lambda (U-8)");
    if (IsOneOf(c.modelId, {"PA-KOLMO"}) &&
        (std::abs(c.numbers["q"] - 5.0 / 3.0) > 1.0e-14 ||
         c.numbers["gap_parameter"] != 0.0))
      return Status::Error(StatusCode::InvalidConfiguration,
          "PA-KOLMO is the source-exact q=5/3, H=0 shape; use PA-QFORM for other values");
  }
  if (IsOneOf(c.modelId, {"QLT-TS03-E-DT"}) &&
      !IsOneOf(c.choices["dt_variant"], {"ts_q3ms", "eb13", "lang24"}))
    return Status::Error(StatusCode::RequiresUserDecision,
                         "dt_variant must be ts_q3ms, eb13, or lang24 (U-7)");
  if (IsOneOf(c.modelId, {"QLT-ZANK98"}) &&
      !IsOneOf(c.choices["variance_convention"], {"per_component", "total"}))
    return Status::Error(StatusCode::RequiresUserDecision,
                         "variance_convention must be per_component or total (U-12)");
  if (c.raw.count("field_normalization") &&
      !IsOneOf(c.choices["field_normalization"], {"magnitude", "radial_amplitude"}))
    return Status::Error(StatusCode::RequiresUserDecision,
                         "field_normalization must resolve U-13 explicitly");
  if (IsOneOf(c.modelId, {"GCR-HELMOD17", "GCR-HELMOD19"}) &&
      c.choices["normalization_convention"] != "d21_numeric_use")
    return Status::Error(StatusCode::RequiresUserDecision,
                         "HelMod requires explicit d21_numeric_use convention");
  if (IsOneOf(c.modelId, {"GCR-DUAN25"}) &&
      c.choices["formula_variant"] != "duan25_printed_d33")
    return Status::Error(StatusCode::RequiresUserDecision,
                         "Duan requires duan25_printed_d33 variant");

  // Runtime state is per configured run (specification Section 21.1).  Once
  // every external value or U-decision required by an executable schema has
  // been supplied, the configuration is ready; printed variants and Chen's
  // direct-kappa status retain their more specific labels.
  if (c.runtimeState == RuntimeState::RequiresExternalInput ||
      c.runtimeState == RuntimeState::RequiresUserDecision)
    c.runtimeState = RuntimeState::ReadyExplicitInputs;
  if (IsOneOf(c.modelId, {"SEP-CHEN24"}) &&
      c.outputKind == LambdaKind::Parallel)
    c.runtimeState = RuntimeState::ReadyExplicitInputs;

  c.fingerprint = ConfigurationFingerprint(c);
  status = ValidateConfiguration(c);
  if (!status.ok()) return status;
  *output = std::move(c);
  return Status::Success();
}

Status ValidateConfiguration(const Configuration& c) {
  const ModelDescriptor* d = Descriptor(c.modelId);
  if (!d) return Status::Error(StatusCode::UnsupportedModel,
                               "unknown model in configuration");
  if (!d->executable) return GateStatus(*d);
  for (const auto& item : c.numbers)
    if (!std::isfinite(item.second)) return Invalid("nonfinite parameter '" + item.first + "'");
  if (c.fingerprint.empty())
    return Status::Error(StatusCode::InvalidConfiguration,
                         "configuration fingerprint is absent");
  return Status::Success();
}

std::string ConfigurationFingerprint(const Configuration& c) {
  // FNV-1a is an identity/checkpoint key, not a cryptographic asset digest.
  // The companion archive retains SHA-256 separately through SourceRecord.
  std::ostringstream canonical;
  canonical << Lower(c.modelId) << '\n';
  for (const auto& item : c.raw) canonical << item.first << '=' << item.second << '\n';
  const std::string text = canonical.str();
  std::uint64_t hash = 1469598103934665603ULL;
  for (unsigned char ch : text) { hash ^= ch; hash *= 1099511628211ULL; }
  std::ostringstream out;
  out << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}

ModelFunction ActiveModelFunction = Unconfigured;

ModelFunction FunctionForModel(const std::string& id) {
  if (IsSepPowerLaw(id)) return EvaluateSepPower;
  if (IsOneOf(id, {"SEP-CHEN24"})) return EvaluateChen;
  if (IsOneOf(id, {"LEGACY-TENISHEV2005AIAA"})) return EvaluateLegacyTenishev;
  if (IsPitch(id)) return EvaluatePitch;
  if (IsOneOf(id, {"MFP-FOCUS-HW13"})) return EvaluateFocus;
  if (IsOneOf(id, {"QLT-TS03-P", "QLT-TS03-E-RS", "QLT-TS03-E-DT", "GCR-EB13"}))
    return EvaluateTS;
  if (IsOneOf(id, {"QLT-ZANK98"})) return EvaluateZank;
  if (IsOneOf(id, {"SHOCK-BOHM", "SHOCK-AFANASIEV15", "SHOCK-MFLAMPA"}))
    return EvaluateShock;
  if (IsOneOf(id, {"GCR-NWU14", "GCR-CORTI19", "GCR-HELMOD17",
                   "GCR-HELMOD19", "GCR-STRAUSS11", "GCR-EFFENBERGER12",
                   "GCR-WANG19", "GCR-TOMASSETTI17", "GCR-PERUGIA21",
                   "GCR-PERUGIA25", "GCR-DUAN25"})) return EvaluateGcr;
  return nullptr;
}

Result Evaluate(const ParticleState& particle, const LocalState& local,
                const Configuration& c) {
  Status status = ValidateConfiguration(c);
  if (!status.ok()) return Failure(c, local, status);
  ModelFunction function = FunctionForModel(c.modelId);
  if (!function) return Failure(c, local, Status::Error(StatusCode::UnsupportedModel,
                                                        "model has no evaluator"));
  return function(particle, local, c);
}

Result EvaluateActive(const ParticleState& particle, const LocalState& local) {
  return ActiveModelFunction(particle, local, gActiveConfiguration);
}

Status SetActiveConfiguration(const Configuration& c) {
  Status status = ValidateConfiguration(c);
  if (!status.ok()) return status;
  ModelFunction function = FunctionForModel(c.modelId);
  if (!function) return Status::Error(StatusCode::UnsupportedModel,
                                      "model has no evaluator");
  // Publication is intentionally the last operation.  Configure only during
  // serial startup; concurrent reconfiguration is outside this API contract.
  gActiveConfiguration = c;
  ActiveModelFunction = function;
  return Status::Success();
}

Status ConfigureActiveModel(const std::string& id,
                            const std::vector<InputParameter>& parameters) {
  Configuration candidate;
  Status status = BuildConfiguration(id, parameters, &candidate);
  return status.ok() ? SetActiveConfiguration(candidate) : status;
}

Configuration GetActiveConfiguration() { return gActiveConfiguration; }

Status EvaluateBatch(const std::vector<ParticleState>& particles,
                     const std::vector<LocalState>& states,
                     const Configuration& c, std::vector<Result>* results) {
  if (!results) return Invalid("batch result is null");
  if (particles.size() != states.size())
    return Invalid("particle and local-state batch lengths differ");
  std::vector<Result> candidate;
  candidate.reserve(particles.size());
  for (std::size_t i = 0; i < particles.size(); ++i)
    candidate.push_back(Evaluate(particles[i], states[i], c));
  *results = std::move(candidate);
  return Status::Success();
}

Status ParallelToRadialSEP(const LambdaValue& parallel, double psi,
                           LambdaValue* radial) {
  if (!radial) return Invalid("radial output is null");
  if (parallel.kind != LambdaKind::Parallel || !FinitePositive(parallel.metres) ||
      !std::isfinite(psi)) return Invalid("parallel lambda and finite psi required");
  const double cosine = std::cos(psi);
  *radial = LambdaValue{parallel.metres * cosine * cosine,
                        LambdaKind::RadialSEP, psi};
  return Status::Success();
}

Status RadialSEPToParallel(const LambdaValue& radial, double psi,
                           LambdaValue* parallel) {
  if (!parallel) return Invalid("parallel output is null");
  if (radial.kind != LambdaKind::RadialSEP || !FinitePositive(radial.metres) ||
      !std::isfinite(psi)) return Invalid("radial SEP lambda and finite psi required");
  const double cosine = std::cos(psi), factor = cosine * cosine;
  if (factor <= 16.0 * std::numeric_limits<double>::epsilon())
    return Status::Error(StatusCode::OutsideModelDomain,
                         "radial-to-parallel conversion is singular at cos(psi)=0");
  *parallel = LambdaValue{radial.metres / factor, LambdaKind::Parallel, psi};
  return Status::Success();
}

Status ParallelToRadialTensor(const LambdaValue& parallel, double kperp,
                              double speed, double psi,
                              LambdaValue* radial) {
  if (!radial) return Invalid("radial tensor output is null");
  if (parallel.kind != LambdaKind::Parallel || !FinitePositive(parallel.metres) ||
      !std::isfinite(kperp) || kperp < 0.0 || !FinitePositive(speed) ||
      !std::isfinite(psi)) return Invalid("invalid tensor projection input");
  const double cosine = std::cos(psi), sine = std::sin(psi);
  const double lambdaPerp = 3.0 * kperp / speed;
  *radial = LambdaValue{parallel.metres * cosine * cosine +
                        lambdaPerp * sine * sine,
                        LambdaKind::RadialTensor, psi};
  return Status::Success();
}

Status PitchAmplitudeFromLambda(const std::string& id, double q, double gap,
                                double vaOverV, OperatorConvention op,
                                const LambdaValue& lambda, double speed,
                                double* amplitude) {
  if (!amplitude) return Invalid("pitch amplitude output is null");
  if (lambda.kind != LambdaKind::Parallel || !FinitePositive(lambda.metres) ||
      !FinitePositive(speed)) return Invalid("positive parallel lambda and speed required");
  double integral = 0.0;
  Status status = PitchShapeIntegral(id, q, gap, vaOverV, &integral);
  if (!status.ok()) return status;
  // Eq. (11): lambda=(3v/8)*(integral/amplitude_effective).  Under HalfD the
  // printed D is divided by two in the operator, so printed amplitude is twice
  // the standard-operator amplitude required for the same transport lambda.
  *amplitude = 3.0 * speed * integral / (8.0 * lambda.metres);
  if (op == OperatorConvention::HalfD) *amplitude *= 2.0;
  return FinitePositive(*amplitude) ? Status::Success()
      : Status::Error(StatusCode::NumericalFailure, "pitch amplitude is nonfinite");
}

Status LambdaFromPitchAmplitude(const std::string& id, double q, double gap,
                                double vaOverV, OperatorConvention op,
                                double amplitude, double speed,
                                LambdaValue* lambda) {
  if (!lambda) return Invalid("lambda output is null");
  if (!FinitePositive(amplitude) || !FinitePositive(speed))
    return Invalid("positive amplitude and speed required");
  double integral = 0.0;
  Status status = PitchShapeIntegral(id, q, gap, vaOverV, &integral);
  if (!status.ok()) return status;
  const double effectiveAmplitude = op == OperatorConvention::HalfD
      ? 0.5 * amplitude : amplitude;
  *lambda = LambdaValue{3.0 * speed * integral /
                        (8.0 * effectiveAmplitude), LambdaKind::Parallel,
                        std::nullopt};
  return Status::Success();
}

Status EvaluatePitchAngleDiffusion(double mu, const ParticleState& particle,
                                   const LocalState&, const Configuration& c,
                                   double* d) {
  if (!d) return Invalid("D_mumu output is null");
  if (!std::isfinite(mu) || mu < -1.0 || mu > 1.0)
    return Invalid("pitch-angle cosine must be in [-1,1]");
  Kinematics k;
  Status status = ComputeKinematics(particle, &k);
  if (!status.ok()) return status;
  double lambda, q = 5.0 / 3.0, gap = 0.0, ratio = 0.0;
  if (!(status = ReadPositive(c, "lambda_parallel_m", &lambda)).ok()) return status;
  auto it = c.numbers.find("q"); if (it != c.numbers.end()) q = it->second;
  it = c.numbers.find("gap_parameter"); if (it != c.numbers.end()) gap = it->second;
  it = c.numbers.find("alfven_to_particle_speed"); if (it != c.numbers.end()) ratio = it->second;
  double amplitude = 0.0;
  if (IsOneOf(c.modelId, {"PA-DROGE-VA"}))
    amplitude = 3.0 * k.speedMPerS /
        (2.0 * (4.0 - q) * (2.0 - q) * lambda);
  else if (IsOneOf(c.modelId, {"PA-EPREM"}) &&
           c.choices.at("lambda_interpretation") == "printed_parameter")
    amplitude = k.speedMPerS / lambda;
  else status = PitchAmplitudeFromLambda(c.modelId, q, gap, ratio,
      c.operatorConvention, LambdaValue{lambda, LambdaKind::Parallel, std::nullopt},
      k.speedMPerS, &amplitude);
  if (!status.ok()) return status;
  const double x = std::abs(mu), oneMinusMu2 = 1.0 - mu * mu;
  if (IsOneOf(c.modelId, {"PA-QFORM", "PA-KOLMO"}))
    *d = amplitude * oneMinusMu2 * (std::pow(x, q - 1.0) + gap);
  else if (IsOneOf(c.modelId, {"PA-QFORM-LANG-PRINTED"}))
    *d = amplitude * (1.0 - mu) * (1.0 - mu) *
         (std::pow(x, q - 1.0) + gap);
  else if (IsOneOf(c.modelId, {"PA-EPS"}))
    *d = 0.5 * amplitude * (x / (1.0 + x) + gap) * oneMinusMu2;
  else if (IsOneOf(c.modelId, {"PA-EPS-PACHECO-PRINTED"}))
    *d = 0.5 * amplitude * (x * (1.0 + x) + gap * oneMinusMu2);
  else if (IsOneOf(c.modelId, {"PA-DROGE-VA"}))
    *d = amplitude * std::pow(mu * mu + ratio * ratio,
                             0.5 * (q - 1.0)) * oneMinusMu2;
  else if (IsOneOf(c.modelId, {"PA-ISO", "PA-EPREM"}))
    *d = 0.5 * amplitude * oneMinusMu2;
  else return Status::Error(StatusCode::UnsupportedModel,
                            "selected model does not define D_mumu");
  return std::isfinite(*d) && *d >= 0.0 ? Status::Success()
      : Status::Error(StatusCode::NumericalFailure, "D_mumu is invalid");
}

Status HeWanFocusingRatio(double x, double* ratio) {
  if (!ratio) return Invalid("focusing ratio output is null");
  if (!FinitePositive(x)) return Invalid("lambda0/L must be positive and finite");
  if (x <= 0.05) {
    const double y = x * x;
    *ratio = 1.0 + y * (-2.0 / 5.0 + y * (17.0 / 105.0 +
        y * (-62.0 / 945.0 + y * (1382.0 / 51975.0))));
  } else if (x > 1.0e100) {
    *ratio = 3.0 * (1.0 - std::tanh(x) / x) / (x * x);
  } else {
    *ratio = 3.0 * (x - std::tanh(x)) / (x * x * x);
  }
  return FinitePositive(*ratio) ? Status::Success()
      : Status::Error(StatusCode::NumericalFailure, "focusing ratio is invalid");
}

Status ParsePublishedScalar(const std::string& text, bool dimensionless,
                            PublishedScalar* output) {
  if (!output) return Invalid("published scalar output is null");
  PublishedScalar parsed; parsed.raw = text;
  double value = 0.0;
  if (ParsePlainDouble(text, &value)) {
    parsed.kind = PublishedScalarKind::PlainNumber; parsed.value = value;
    *output = parsed; return Status::Success();
  }
  if (text.find('/') != std::string::npos && ParseDimensionless(text, &value)) {
    if (!dimensionless)
      return Invalid("exact rational literals are allowed only for dimensionless fields");
    parsed.kind = PublishedScalarKind::ExactRational; parsed.value = value;
    *output = parsed; return Status::Success();
  }
  // A number followed by a nonempty unit is retained as a typed value, not
  // converted without a field-specific schema.
  std::istringstream input(text);
  std::string number, unit, extra;
  if ((input >> number >> unit) && !(input >> extra) && ParsePlainDouble(number, &value)) {
    parsed.kind = PublishedScalarKind::NumberWithUnit;
    parsed.value = value; parsed.unit = unit; *output = parsed;
    return Status::Success();
  }
  parsed.kind = PublishedScalarKind::NonScalar;
  *output = parsed;
  return Status::Error(StatusCode::RequiresUserDecision,
      "published text is a range, qualifier, uncertainty, inequality, or unresolved expression");
}

Status LoadSourceRecord(const std::string& bundle, const std::string& relative,
                        const std::string& id, SourceRecord* output) {
  if (!output) return Invalid("source record output is null");
  if (bundle.empty() || relative.empty() || id.empty())
    return Invalid("bundle directory, relative file, and record ID are required");
  if (relative.find("..") != std::string::npos || relative.front() == '/')
    return Invalid("relative source path must remain inside the bundle");
  std::ifstream file(bundle + "/" + relative);
  if (!file) return Missing(relative);
  std::ostringstream contents; contents << file.rdbuf();
  const std::string text = contents.str();
  const std::string needle = "\"id\"";
  std::size_t match = 0, objectBegin = std::string::npos, objectEnd = std::string::npos;
  while ((match = text.find(needle, match)) != std::string::npos) {
    const std::size_t colon = text.find(':', match + needle.size());
    const std::size_t quote = colon == std::string::npos ? colon : text.find('"', colon + 1);
    const std::size_t close = quote == std::string::npos ? quote : text.find('"', quote + 1);
    if (quote != std::string::npos && close != std::string::npos &&
        text.substr(quote + 1, close - quote - 1) == id) {
      objectBegin = text.rfind('{', match); break;
    }
    match += needle.size();
  }
  if (objectBegin == std::string::npos)
    return Missing("source record " + id + " in " + relative);
  bool inString = false, escape = false; int depth = 0;
  for (std::size_t i = objectBegin; i < text.size(); ++i) {
    const char ch = text[i];
    if (inString) {
      if (escape) escape = false;
      else if (ch == '\\') escape = true;
      else if (ch == '"') inString = false;
    } else if (ch == '"') inString = true;
    else if (ch == '{') ++depth;
    else if (ch == '}' && --depth == 0) { objectEnd = i + 1; break; }
  }
  if (objectEnd == std::string::npos)
    return Invalid("unterminated JSON object for source record " + id);

  std::ifstream manifest(bundle + "/SHA256SUMS");
  if (!manifest) return Missing("SHA256SUMS");
  std::string line, digest;
  while (std::getline(manifest, line)) {
    std::istringstream row(line); std::string hash, name;
    if ((row >> hash >> name) && name == relative) { digest = hash; break; }
  }
  if (digest.empty()) return Missing("SHA256 manifest entry for " + relative);
  *output = SourceRecord{relative, id,
                         text.substr(objectBegin, objectEnd - objectBegin), digest};
  return Status::Success();
}

}  // namespace MeanFreePath
}  // namespace SEP
