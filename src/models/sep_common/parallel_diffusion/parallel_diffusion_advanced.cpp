#include "parallel_diffusion_advanced.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace SEP {
namespace ParallelDiffusion {
namespace Internal {
namespace {

constexpr double NonlinearResidualLimit = 1.0e-8;
constexpr const char* SpecificationVersion = "1.4";
constexpr const char* SoftwareVersion = "parallel-diffusion-pd10-v1";
constexpr const char* NlgceCoefficientDigest =
    "parallel:7cdc5ab9cddda0c7295a25bfcaff4ba3422002da12ac1b4660a98dcc64eaa762;"
    "perpendicular:7371b6557f40c5efa7a48460a5f1b2971374a2c43d0eb3e209f0a8e62670356a";
#include "nlgce_f_coefficients.inc"

bool Finite(double x) { return std::isfinite(x); }
bool Positive(double x) { return Finite(x) && x > 0.0; }

double Magnitude(const std::array<double, 3>& x) {
  return std::hypot(x[0], std::hypot(x[1], x[2]));
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

ParallelResult Failure(StatusCode code, const std::string& detail,
                       const ModelConfiguration& configuration,
                       const LocalState& local) {
  ParallelResult result;
  result.status = Status::Error(code, detail);
  result.provenance.requestedModelId = ModelName(configuration.model);
  result.provenance.evaluatedModelId = ModelName(configuration.model);
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.softwareVersion = SoftwareVersion;
  result.provenance.configurationFingerprint =
      ConfigurationFingerprint(configuration);
  result.provenance.backgroundRevision = local.backgroundRevision;
  result.provenance.turbulenceRevision = local.turbulenceRevision;
  return result;
}

ParallelResult Success(const ModelConfiguration& configuration,
                       const LocalState& local,
                       const ParticleKinematics& particle,
                       double lambdaM) {
  if (!Positive(lambdaM))
    return Failure(StatusCode::OutsideModelDomain,
                   "parallel mean free path is not finite and positive",
                   configuration, local);
  const double kappa = particle.speedMPerS * lambdaM / 3.0;
  if (!Positive(kappa))
    return Failure(StatusCode::OutsideModelDomain,
                   "parallel diffusion coefficient is not representable",
                   configuration, local);
  ParallelResult result;
  result.status = Status::Success();
  result.lambdaParallelM = lambdaM;
  result.kappaParallelM2PerS = kappa;
  result.provenance.requestedModelId = ModelName(configuration.model);
  result.provenance.evaluatedModelId = ModelName(configuration.model);
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.softwareVersion = SoftwareVersion;
  result.provenance.configurationFingerprint =
      ConfigurationFingerprint(configuration);
  result.provenance.backgroundRevision = local.backgroundRevision;
  result.provenance.turbulenceRevision = local.turbulenceRevision;
  return result;
}

Status FieldAndTurbulence(const LocalState& local, double* fieldT,
                         const LocalState::TurbulenceState** turbulence) {
  if (!local.meanFieldT.has_value())
    return Status::Error(StatusCode::MissingInput,
                         "selected model requires mean_B_T");
  *fieldT = Magnitude(*local.meanFieldT);
  if (!Positive(*fieldT))
    return Status::Error(StatusCode::InvalidBackground,
                         "selected model requires |mean_B_T|>0");
  if (!local.turbulence.has_value())
    return Status::Error(StatusCode::MissingInput,
                         "selected model requires canonical turbulence state");
  *turbulence = &*local.turbulence;
  return Status::Success();
}

Status SlabState(const LocalState& local, double* fieldT, double* varianceT2,
                 double* lengthM, double* index,
                 bool requireIndexGreaterThanOne = true) {
  const LocalState::TurbulenceState* t = nullptr;
  Status status = FieldAndTurbulence(local, fieldT, &t);
  if (!status.ok()) return status;
  if (!t->slabVarianceT2.has_value() ||
      !t->slabBendoverLengthM.has_value() || !t->inertialIndex.has_value())
    return Status::Error(StatusCode::MissingInput,
        "slab closure requires slab_variance_T2, slab_bendover_length_m, and inertial_index");
  *varianceT2 = *t->slabVarianceT2;
  *lengthM = *t->slabBendoverLengthM;
  *index = *t->inertialIndex;
  if (!Finite(*varianceT2) || *varianceT2 < 0.0 || !Positive(*lengthM) ||
      !Finite(*index) || (requireIndexGreaterThanOne && *index <= 1.0))
    return Status::Error(StatusCode::InvalidBackground,
        requireIndexGreaterThanOne
            ? "slab variance must be nonnegative, length positive, and index greater than one"
            : "slab variance must be nonnegative, length positive, and index finite");
  if (*varianceT2 == 0.0)
    return Status::Error(StatusCode::InfiniteMeanFreePath,
                         "zero slab variance supplies no scattering");
  return Status::Success();
}

double Cnu(double nu) {
  // Equation (27). lgamma avoids forming a large numerator and denominator
  // independently for otherwise admissible indices.
  return std::exp(std::lgamma(nu) - std::lgamma(nu - 0.5)) /
      (2.0 * std::sqrt(std::acos(-1.0)));
}

struct Integral {
  bool ok = false;
  double value = 0.0;
  double error = 0.0;
};

double Simpson(double a, double b, double fa, double fm, double fb) {
  return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

bool AdaptiveSimpsonRec(const std::function<double(double)>& f,
                        double a, double b, double fa, double fm, double fb,
                        double whole, double tolerance, int depth,
                        double* value, double* error) {
  const double m = 0.5 * (a + b);
  const double lm = 0.5 * (a + m);
  const double rm = 0.5 * (m + b);
  const double flm = f(lm);
  const double frm = f(rm);
  if (!Finite(flm) || !Finite(frm)) return false;
  const double left = Simpson(a, m, fa, flm, fm);
  const double right = Simpson(m, b, fm, frm, fb);
  const double delta = left + right - whole;
  if (depth <= 0 || std::fabs(delta) <= 15.0 * tolerance) {
    *value += left + right + delta / 15.0;
    *error += std::fabs(delta) / 15.0;
    return depth > 0 || std::fabs(delta) <= 15.0 * tolerance;
  }
  return AdaptiveSimpsonRec(f, a, m, fa, flm, fm, left,
                            tolerance / 2.0, depth - 1, value, error) &&
         AdaptiveSimpsonRec(f, m, b, fm, frm, fb, right,
                            tolerance / 2.0, depth - 1, value, error);
}

Integral Integrate(const std::function<double(double)>& f, double a, double b,
                   double relativeTolerance, int maximumRefinements) {
  Integral result;
  const double m = 0.5 * (a + b);
  const double fa = f(a), fm = f(m), fb = f(b);
  if (!Finite(fa) || !Finite(fm) || !Finite(fb)) return result;
  const double initial = Simpson(a, b, fa, fm, fb);
  // Relative tolerance must scale with the integral itself.  Using max(1,I)
  // would turn weak-turbulence spectral integrals (I << 1 in normalized
  // units) into loose absolute tests and biases the E/F regression states.
  const double tolerance = std::max(1.0e-300,
      relativeTolerance * std::fabs(initial));
  const bool converged = AdaptiveSimpsonRec(
      f, a, b, fa, fm, fb, initial, tolerance, maximumRefinements,
      &result.value, &result.error);
  // A single very narrow leaf may exhaust the recursion counter after the
  // accumulated a-posteriori estimate is already within the requested total
  // tolerance.  Accept on the estimate, not merely the recursion flag.
  result.ok = converged || result.error <= tolerance;
  return result;
}

Integral IntegratePieces(const std::function<double(double)>& f,
                         std::vector<double> points, double relativeTolerance,
                         int maximumRefinements) {
  Integral total;
  total.ok = true;
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  for (std::size_t i = 1; i < points.size(); ++i) {
    const Integral part = Integrate(f, points[i - 1], points[i],
                                    relativeTolerance, maximumRefinements);
    if (!part.ok) total.ok = false;
    total.value += part.value;
    total.error += part.error;
  }
  return total;
}

double GaussLegendre(const std::function<double(double)>& f,
                     double a, double b, int order, bool* ok) {
  // Nodes and weights are generated from the Legendre roots so the nested
  // broadened-slab test does not depend on a compiler-specific table.  The
  // root iteration is purely numerical infrastructure; comparing two orders
  // below supplies the a-posteriori convergence check.
  const double midpoint = 0.5 * (a + b);
  const double half = 0.5 * (b - a);
  double sum = 0.0;
  const int roots = (order + 1) / 2;
  for (int i = 0; i < roots; ++i) {
    double z = std::cos(std::acos(-1.0) * (i + 0.75) / (order + 0.5));
    double derivative = 0.0;
    for (int iteration = 0; iteration < 30; ++iteration) {
      double p0 = 1.0, p1 = z;
      for (int n = 2; n <= order; ++n) {
        const double p2 = ((2.0 * n - 1.0) * z * p1 - (n - 1.0) * p0) / n;
        p0 = p1; p1 = p2;
      }
      derivative = order * (z * p1 - p0) / (z * z - 1.0);
      const double step = p1 / derivative;
      z -= step;
      if (std::fabs(step) < 2.0e-15) break;
    }
    const double weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
    const double left = f(midpoint - half * z);
    const double right = f(midpoint + half * z);
    if (!Finite(left) || !Finite(right)) { *ok = false; return 0.0; }
    sum += weight * (left + right);
  }
  return half * sum;
}

Status SpectrumValue(const SpectrumParameters& spectrum, double k,
                     double* power) {
  if (!Positive(k) || !power)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "invalid supplied-spectrum query");
  const std::vector<double>& x = spectrum.wavenumberRadPerM;
  const std::vector<double>& y = spectrum.powerT2M;
  if (x.size() < 2 || x.size() != y.size())
    return Status::Error(StatusCode::InvalidConfiguration,
                         "supplied spectrum has invalid array shape");
  auto tail = [&](bool low) -> Status {
    const SpectrumTailPolicy policy = low ? spectrum.lowKPolicy
                                          : spectrum.highKPolicy;
    if (policy == SpectrumTailPolicy::OutOfDomain)
      return Status::Error(StatusCode::MissingInput,
                           "resonance lies outside supplied spectrum coverage");
    if (policy == SpectrumTailPolicy::Zero) {
      *power = 0.0;
      return Status::Success();
    }
    const std::size_t i = low ? 0 : x.size() - 1;
    const double exponent = low ? spectrum.lowKPowerIndex
                                : spectrum.highKPowerIndex;
    *power = y[i] * std::exp(exponent * std::log(k / x[i]));
    return Finite(*power) && *power >= 0.0
        ? Status::Success()
        : Status::Error(StatusCode::OutsideModelDomain,
                        "supplied spectrum tail is not representable");
  };
  if (k < x.front()) return tail(true);
  if (k > x.back()) return tail(false);
  const auto upper = std::lower_bound(x.begin(), x.end(), k);
  if (upper == x.begin()) { *power = y.front(); return Status::Success(); }
  if (upper == x.end()) { *power = y.back(); return Status::Success(); }
  const std::size_t hi = static_cast<std::size_t>(upper - x.begin());
  const std::size_t lo = hi - 1;
  if (y[lo] == 0.0 || y[hi] == 0.0) {
    // A true zero is physical absence of resonant power.  Interpolating its
    // logarithm would manufacture a floor, so the entire touching interval is
    // represented as zero and the caller reports infinite mean free path.
    *power = 0.0;
    return Status::Success();
  }
  const double weight = std::log(k / x[lo]) / std::log(x[hi] / x[lo]);
  *power = std::exp((1.0 - weight) * std::log(y[lo]) +
                    weight * std::log(y[hi]));
  return Status::Success();
}

Status MultirangeSpectrumValue(const SpectrumParameters& spectrum,
                               double variance, double bendoverLength,
                               double inertialIndex, double k,
                               double* power) {
  if (!power || !Positive(variance) || !Positive(bendoverLength) ||
      !Positive(k) || !Finite(inertialIndex))
    return Status::Error(StatusCode::InvalidBackground,
        "multirange spectrum requires positive variance, bend-over length, "
        "and wavenumber plus a finite inertial index");

  // Equation (35) is normalized analytically, rather than by a finite
  // numerical k interval.  x_d must be beyond the bend-over point x=1; this
  // depends on the runtime bend-over length and therefore cannot be checked
  // completely when the parser reads k_d.  expm1 evaluates
  // (x_d^(1-s)-1)/(1-s) accurately near the explicitly specified s=1 limit.
  const double xD = spectrum.dissipationWavenumberRadPerM * bendoverLength;
  if (!(xD > 1.0) || !Finite(xD))
    return Status::Error(StatusCode::InvalidBackground,
        "multirange spectrum requires k_d*ell_s > 1");
  const double logXD = std::log(xD);
  const double oneMinusS = 1.0 - inertialIndex;
  const double middle = oneMinusS == 0.0
      ? logXD : std::expm1(oneMinusS * logXD) / oneMinusS;
  const double xDOneMinusS = std::exp(oneMinusS * logXD);
  const double normalization =
      1.0 / (spectrum.energyRangeIndex + 1.0) + middle +
      xDOneMinusS / (spectrum.dissipationIndex - 1.0);
  if (!Positive(normalization) || !Finite(normalization))
    return Status::Error(StatusCode::OutsideModelDomain,
        "multirange spectrum normalization is not representable");

  const double x = k * bendoverLength;
  double shape = 0.0;
  if (x <= 1.0) {
    shape = std::pow(x, spectrum.energyRangeIndex);
  } else if (x <= xD) {
    shape = std::pow(x, -inertialIndex);
  } else {
    // The x_d^(s_d-s) factor makes the inertial and dissipation branches
    // continuous at x_d.  P_s has units T^2 m and integrates over dk to the
    // runtime total two-component slab variance [T^2].
    shape = std::exp((spectrum.dissipationIndex - inertialIndex) * logXD -
                     spectrum.dissipationIndex * std::log(x));
  }
  *power = variance * bendoverLength * shape / normalization;
  return Finite(*power) && *power >= 0.0
      ? Status::Success()
      : Status::Error(StatusCode::OutsideModelDomain,
                      "multirange spectrum value is not representable");
}

Status ValidateSpectrumNormalization(const SpectrumParameters& spectrum,
                                     double relativeTolerance) {
  // The smooth and multirange forms are normalized by their defining closed
  // forms (Equations (27) and (35)); the latter's runtime x_d constraint is
  // checked when ell_s is available.  Only supplied arrays need an explicit
  // numerical normalization audit here.
  if (spectrum.form == SpectrumForm::SmoothBendover ||
      spectrum.form == SpectrumForm::Multirange)
    return Status::Success();
  const auto& k = spectrum.wavenumberRadPerM;
  const auto& p = spectrum.powerT2M;
  double integral = 0.0;
  if (spectrum.lowKPolicy == SpectrumTailPolicy::PowerLaw) {
    if (!(spectrum.lowKPowerIndex > -1.0))
      return Status::Error(StatusCode::InconsistentSpectrum,
                           "low-k power-law tail has divergent variance");
    integral += p.front() * k.front() /
                (spectrum.lowKPowerIndex + 1.0);
  } else if (spectrum.lowKPolicy == SpectrumTailPolicy::OutOfDomain) {
    return Status::Error(StatusCode::MissingInput,
                         "spectrum normalization needs an explicit low-k tail");
  }
  for (std::size_t i = 1; i < k.size(); ++i) {
    if (p[i - 1] == 0.0 || p[i] == 0.0) continue;
    const double exponent = std::log(p[i] / p[i - 1]) /
                            std::log(k[i] / k[i - 1]);
    const double ratio = k[i] / k[i - 1];
    const double segment = std::fabs(exponent + 1.0) < 1.0e-12
        ? p[i - 1] * k[i - 1] * std::log(ratio)
        : p[i - 1] * k[i - 1] *
          (std::pow(ratio, exponent + 1.0) - 1.0) /
          (exponent + 1.0);
    integral += segment;
  }
  if (spectrum.highKPolicy == SpectrumTailPolicy::PowerLaw) {
    if (!(spectrum.highKPowerIndex < -1.0))
      return Status::Error(StatusCode::InconsistentSpectrum,
                           "high-k power-law tail has divergent variance");
    integral += p.back() * k.back() /
                (-spectrum.highKPowerIndex - 1.0);
  } else if (spectrum.highKPolicy == SpectrumTailPolicy::OutOfDomain) {
    return Status::Error(StatusCode::MissingInput,
                         "spectrum normalization needs an explicit high-k tail");
  }
  if (!Positive(integral) ||
      std::fabs(integral / spectrum.declaredVarianceT2 - 1.0) >
          relativeTolerance)
    return Status::Error(StatusCode::InconsistentSpectrum,
        "canonical supplied spectrum does not integrate to declared_variance_T2");
  return Status::Success();
}

Status QltDmu(double mu, const ParticleState& particle,
              const ParticleKinematics& kin, const LocalState& local,
              const QltSlabParameters& parameters, double* dMuMu) {
  if (!dMuMu || !Finite(mu) || std::fabs(mu) > 1.0)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "pitch angle mu must be finite and in [-1,1]");
  double field = 0.0, variance = 0.0, length = 0.0, index = 0.0;
  Status state;
  if (parameters.spectrum.form == SpectrumForm::SmoothBendover ||
      parameters.spectrum.form == SpectrumForm::Multirange) {
    state = SlabState(local, &field, &variance, &length, &index,
        parameters.spectrum.form != SpectrumForm::Multirange);
    if (!state.ok()) return state;
  } else {
    if (!local.meanFieldT.has_value())
      return Status::Error(StatusCode::MissingInput,
                           "supplied-spectrum QLT requires mean_B_T");
    field = Magnitude(*local.meanFieldT);
    if (!Positive(field))
      return Status::Error(StatusCode::InvalidBackground,
                           "supplied-spectrum QLT requires |mean_B_T|>0");
  }
  const double oneMinusMu2 = 1.0 - mu * mu;
  if (oneMinusMu2 == 0.0) { *dMuMu = 0.0; return Status::Success(); }
  if (mu == 0.0) { *dMuMu = 0.0; return Status::Success(); }
  const double omega = std::fabs(particle.chargeC) * field /
      (kin.gamma * particle.massKg);
  const double rL = particle.momentumKgMPerS /
      (std::fabs(particle.chargeC) * field);
  if (parameters.spectrum.form == SpectrumForm::SmoothBendover) {
    const double nu = 0.5 * index;
    const double rstar = rL / length;
    *dMuMu = std::acos(-1.0) * Cnu(nu) * (variance / (field * field)) *
        (kin.speedMPerS / length) * std::pow(rstar, index - 2.0) *
        oneMinusMu2 * std::pow(std::fabs(mu), index - 1.0) *
        std::pow(1.0 + rstar * rstar * mu * mu, -0.5 * index);
  } else {
    double power = 0.0;
    const double resonantK = 1.0 / (rL * std::fabs(mu));
    state = parameters.spectrum.form == SpectrumForm::Multirange
        ? MultirangeSpectrumValue(parameters.spectrum, variance, length,
                                  index, resonantK, &power)
        : SpectrumValue(parameters.spectrum, resonantK, &power);
    if (!state.ok()) return state;
    *dMuMu = std::acos(-1.0) * omega * omega * oneMinusMu2 * power /
        (4.0 * field * field * kin.speedMPerS * std::fabs(mu));
  }
  if (!Finite(*dMuMu) || *dMuMu < 0.0)
    return Status::Error(StatusCode::OutsideModelDomain,
                         "QLT pitch-angle coefficient is not representable");
  return Status::Success();
}

ParallelResult EvaluateQlt(const ParticleState& particle,
                           const LocalState& local,
                           const ModelConfiguration& configuration,
                           bool inertial) {
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  if (!inertial) {
    status = ValidateSpectrumNormalization(configuration.qltSlab.spectrum,
        configuration.qltSlab.numerical.relativeTolerance);
    if (!status.ok())
      return Failure(status.code, status.detail, configuration, local);
  }
  double field = 0.0, variance = 0.0, length = 0.0, index = 0.0;
  const SpectrumForm spectrumForm = configuration.qltSlab.spectrum.form;
  const bool analytic = inertial || spectrumForm == SpectrumForm::SmoothBendover;
  const bool runtimeSpectrum = analytic || spectrumForm == SpectrumForm::Multirange;
  if (runtimeSpectrum) {
    status = SlabState(local, &field, &variance, &length, &index,
        spectrumForm != SpectrumForm::Multirange);
    if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
    if (analytic && !(index < 2.0))
      return Failure(StatusCode::InfiniteMeanFreePath,
                     "selected QLT spectrum does not give a convergent pitch-angle integral",
                     configuration, local);
    if (spectrumForm == SpectrumForm::Multirange) {
      if (!(configuration.qltSlab.spectrum.dissipationIndex < 2.0))
        return Failure(StatusCode::InfiniteMeanFreePath,
            "multirange dissipation index s_d >= 2 makes Equation (20) divergent",
            configuration, local);
      double normalizationProbe = 0.0;
      status = MultirangeSpectrumValue(configuration.qltSlab.spectrum,
          variance, length, index, 1.0 / length, &normalizationProbe);
      if (!status.ok())
        return Failure(status.code, status.detail, configuration, local);
    }
  } else {
    if (!local.meanFieldT.has_value())
      return Failure(StatusCode::MissingInput,
                     "supplied-spectrum QLT requires mean_B_T",
                     configuration, local);
    field = Magnitude(*local.meanFieldT);
    if (!Positive(field))
      return Failure(StatusCode::InvalidBackground,
                     "supplied-spectrum QLT requires |mean_B_T|>0",
                     configuration, local);
  }
  const double rL = particle.momentumKgMPerS /
      (std::fabs(particle.chargeC) * field);
  const double r = analytic ? rL / length : 0.0;
  const double eps2 = analytic ? variance / (field * field) : 0.0;
  const double c = analytic ? Cnu(0.5 * index) : 0.0;
  double lambda = 0.0;
  std::optional<double> lambdaSlope;
  if (inertial) {
    // Equation (31), retained as a separate approximation identity.  The
    // evaluator reports no claim that r* is small; validity-range policy is a
    // production choice (D07), while the arithmetic remains reproducible.
    lambda = 3.0 * length * std::pow(r, 2.0 - index) /
        (2.0 * std::acos(-1.0) * c * (2.0 - index) *
         (4.0 - index) * eps2);
    lambdaSlope = 2.0 - index;
  } else if (configuration.qltSlab.spectrum.form ==
             SpectrumForm::SmoothBendover) {
    // Equation (37) removes the integrable mu^(1-s) endpoint singularity in
    // Equation (30).  No artificial pitch-angle cutoff is introduced.
    const double inv = 1.0 / (2.0 - index);
    const auto transformed = [=](double z) {
      const double mu = z == 0.0 ? 0.0 : std::pow(z, inv);
      return inv * (1.0 - mu * mu) *
          std::pow(1.0 + r * r * mu * mu, 0.5 * index);
    };
    const Integral integral = Integrate(transformed, 0.0, 1.0,
        configuration.qltSlab.numerical.relativeTolerance,
        configuration.qltSlab.numerical.maximumRefinements);
    if (!integral.ok)
      return Failure(StatusCode::IntegrationFailed,
                     "QLT Equation (30) quadrature did not converge",
                     configuration, local);
    lambda = 3.0 * length * std::pow(r, 2.0 - index) * integral.value /
        (4.0 * std::acos(-1.0) * c * eps2);
    const auto derivativeIntegrand = [=](double z) {
      const double mu = z == 0.0 ? 0.0 : std::pow(z, inv);
      const double r2mu2 = r * r * mu * mu;
      return transformed(z) * index * r2mu2 / (1.0 + r2mu2);
    };
    const Integral derivative = Integrate(derivativeIntegrand, 0.0, 1.0,
        configuration.qltSlab.numerical.relativeTolerance,
        configuration.qltSlab.numerical.maximumRefinements);
    if (!derivative.ok)
      return Failure(StatusCode::IntegrationFailed,
                     "QLT rigidity-derivative quadrature did not converge",
                     configuration, local);
    lambdaSlope = 2.0 - index + derivative.value / integral.value;
  } else {
    Status spectrumStatus = Status::Success();
    bool zeroPower = false;
    const auto integrand = [&](double z) {
      // The transformed endpoint is a mathematical limit.  Querying exactly
      // mu=0 would ask a supplied spectrum for k=infinity and then confuse a
      // vanishing D_mu_mu with a divergent transformed integrand.
      const double ze = z == 0.0 ? 1.0e-12 : z;
      const double mu = ze * ze * ze;
      double d = 0.0;
      const Status q = QltDmu(mu, particle, kin, local,
                              configuration.qltSlab, &d);
      if (!q.ok()) { spectrumStatus = q; return 0.0; }
      if (d == 0.0 && mu != 1.0) { zeroPower = true; return 0.0; }
      return d > 0.0 ? 3.0 * z * z *
          (1.0 - mu * mu) * (1.0 - mu * mu) / d : 0.0;
    };
    // Resonance maps k to mu=1/(r_L k), while the endpoint transform uses
    // mu=z^3.  Supplying every physical spectral break to the piecewise
    // integrator prevents an adaptive panel from straddling a slope or tail-
    // policy discontinuity.  Multirange breaks are x=1 and x=x_d; supplied
    // spectra break at every tabulated knot.  Values outside 0<mu<1 do not
    // intersect the pitch-angle integration domain.
    std::vector<double> pitchBreaks{0.0, 1.0};
    auto addWavenumberBreak = [&](double kBreak) {
      const double muBreak = 1.0 / (rL * kBreak);
      if (muBreak > 0.0 && muBreak < 1.0 && Finite(muBreak))
        pitchBreaks.push_back(std::cbrt(muBreak));
    };
    if (spectrumForm == SpectrumForm::Multirange) {
      addWavenumberBreak(1.0 / length);
      addWavenumberBreak(
          configuration.qltSlab.spectrum.dissipationWavenumberRadPerM);
    } else {
      for (double kBreak :
           configuration.qltSlab.spectrum.wavenumberRadPerM)
        addWavenumberBreak(kBreak);
    }
    const Integral integral = IntegratePieces(integrand, pitchBreaks,
        configuration.qltSlab.numerical.relativeTolerance,
        configuration.qltSlab.numerical.maximumRefinements);
    if (!spectrumStatus.ok())
      return Failure(spectrumStatus.code, spectrumStatus.detail,
                     configuration, local);
    if (zeroPower)
      return Failure(StatusCode::InfiniteMeanFreePath,
                     "zero resonant spectrum power makes Equation (20) divergent",
                     configuration, local);
    if (!integral.ok)
      return Failure(StatusCode::IntegrationFailed,
                     "spectral pitch-angle quadrature did not converge",
                     configuration, local);
    lambda = 3.0 * kin.speedMPerS * integral.value / 4.0;
  }
  ParallelResult result = Success(configuration, local, kin, lambda);
  if (result.status.ok()) {
    if (lambdaSlope.has_value()) {
      result.dLnLambdaDLnRigidity = *lambdaSlope;
      result.dLnKappaDLnRigidity = *lambdaSlope +
          1.0 / (kin.gamma * kin.gamma);
    }
    result.diagnosticMask |= QltWeakPerturbationConcern |
                             DerivativeUnavailable;
  }
  return result;
}

double ShapeIntegral(double q, double h, const NumericalParameters& numerical,
                     bool* ok) {
  if (h == 0.0 && q < 2.0) {
    *ok = true;
    return 4.0 / ((2.0 - q) * (4.0 - q));
  }
  const auto f = [=](double z) {
    const double mu = z * z * z;
    const double shape = mu == 0.0
        ? (q > 1.0 ? 0.0 : (q == 1.0 ? 1.0 :
           std::numeric_limits<double>::infinity()))
        : std::pow(mu, q - 1.0);
    const double denominator = shape + h;
    if (std::isinf(denominator)) return 0.0;
    return 6.0 * z * z * (1.0 - mu * mu) / denominator;
  };
  const Integral value = Integrate(f, 0.0, 1.0,
      numerical.relativeTolerance, numerical.maximumRefinements);
  *ok = value.ok && Positive(value.value);
  return value.value;
}

ParallelResult EvaluatePrescribed(const ParticleState& particle,
                                  const LocalState& local,
                                  const ModelConfiguration& configuration) {
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  const PrescribedLambdaMuParameters& p = configuration.prescribedLambdaMu;
  bool ok = false;
  const double integral = ShapeIntegral(p.qMu, p.hMu, p.numerical, &ok);
  if (!ok)
    return Failure(p.hMu == 0.0 && p.qMu >= 2.0
                       ? StatusCode::InfiniteMeanFreePath
                       : StatusCode::IntegrationFailed,
                   "regularized pitch-angle shape integral did not converge",
                   configuration, local);
  const double d0 = p.amplitudeMode == PitchAngleAmplitudeMode::TargetLambda
      ? 3.0 * kin.speedMPerS * integral / (8.0 * p.targetLambdaM)
      : p.fixedD0PerS;
  const double lambda = 3.0 * kin.speedMPerS * integral / (8.0 * d0);
  ParallelResult result = Success(configuration, local, kin, lambda);
  if (result.status.ok()) result.diagnosticMask |= DerivativeUnavailable;
  return result;
}

struct TurbulenceRatios {
  double r = 0.0;
  double fs = 0.0;
  double eps2 = 0.0;
  double rho = 0.0;
  double fieldT = 0.0;
  double slabLengthM = 0.0;
};

Status Ratios(const ParticleState& particle, const LocalState& local,
              TurbulenceRatios* out) {
  const LocalState::TurbulenceState* t = nullptr;
  Status status = FieldAndTurbulence(local, &out->fieldT, &t);
  if (!status.ok()) return status;
  if (!t->slabVarianceT2.has_value() || !t->twoDVarianceT2.has_value() ||
      !t->slabBendoverLengthM.has_value() ||
      !t->twoDBendoverLengthM.has_value())
    return Status::Error(StatusCode::MissingInput,
        "two-component closure requires both variances and bend-over lengths");
  const double slab = *t->slabVarianceT2;
  const double twoD = *t->twoDVarianceT2;
  out->slabLengthM = *t->slabBendoverLengthM;
  const double l2 = *t->twoDBendoverLengthM;
  if (!Positive(slab) || !Positive(twoD) || !Positive(out->slabLengthM) ||
      !Positive(l2))
    return Status::Error(StatusCode::InvalidBackground,
        "mixed two-component closure requires positive component variances and lengths");
  const double total = slab + twoD;
  out->fs = slab / total;
  out->eps2 = total / (out->fieldT * out->fieldT);
  out->rho = out->slabLengthM / l2;
  const double rL = particle.momentumKgMPerS /
      (std::fabs(particle.chargeC) * out->fieldT);
  out->r = rL / out->slabLengthM;
  if (!Positive(out->r) || !(out->fs > 0.0 && out->fs < 1.0) ||
      !Positive(out->eps2) || !Positive(out->rho))
    return Status::Error(StatusCode::InvalidBackground,
                         "invalid two-component dimensionless ratios");
  return Status::Success();
}

double Ax(const TurbulenceRatios& s) {
  // Corrected multiplicative Equation (45): the first denominator term is
  // (xi/(1+xi))/epsilon, not an exponentiation by 1/epsilon.
  const double epsilon = std::sqrt(s.eps2);
  const double xi = s.r / (Cnu(5.0 / 6.0) * epsilon);
  return 0.5 * std::sqrt(s.fs /
      ((xi / (1.0 + xi)) / epsilon + epsilon / (2.0 * xi)));
}

double APrime2(const TurbulenceRatios& s) {
  return 1.0 / (std::sqrt(1.0 / s.rho) / s.fs +
                (4.0 / 3.0) / (1.0 - s.fs));
}

double Poly(const double d[6][4][4][3], const double x[4]) {
  // Section 11.3 requires Horner nesting l,k,j,i.  The compiled array uses
  // explicit d[i][j][k][l] order generated from the audited CSV indices.
  double ai[6];
  for (int i = 0; i < 6; ++i) {
    double bj[4];
    for (int j = 0; j < 4; ++j) {
      double ck[4];
      for (int k = 0; k < 4; ++k)
        ck[k] = (d[i][j][k][2] * x[3] + d[i][j][k][1]) * x[3] +
                d[i][j][k][0];
      bj[j] = ((ck[3] * x[2] + ck[2]) * x[2] + ck[1]) * x[2] + ck[0];
    }
    ai[i] = ((bj[3] * x[1] + bj[2]) * x[1] + bj[1]) * x[1] + bj[0];
  }
  return ((((ai[5] * x[0] + ai[4]) * x[0] + ai[3]) * x[0] + ai[2]) *
          x[0] + ai[1]) * x[0] + ai[0];
}

double PolyDerivative(const double d[6][4][4][3], const double x[4],
                      int axis) {
  // The derivative is evaluated independently as the exact differentiated
  // finite sum.  This makes the index power explicit and avoids numerical
  // differencing noise in the transport rigidity slope.
  double sum = 0.0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        for (int l = 0; l < 3; ++l) {
          const int powers[4] = {i, j, k, l};
          if (powers[axis] == 0) continue;
          double term = d[i][j][k][l] * powers[axis];
          for (int a = 0; a < 4; ++a)
            term *= std::pow(x[a], powers[a] - (a == axis ? 1 : 0));
          sum += term;
        }
  return sum;
}

bool NlgceSpatialGradient(const LocalState& local,
                          const TurbulenceRatios& s,
                          const double polynomial[6][4][4][3],
                          const double x[4],
                          std::array<double, 3>* gradLogLambda) {
  if (!local.meanFieldT.has_value() || !local.dMeanFieldDxTPerM.has_value() ||
      !local.turbulence.has_value()) return false;
  const auto& t = *local.turbulence;
  if (!t.gradSlabVarianceT2PerM.has_value() ||
      !t.gradTwoDVarianceT2PerM.has_value() ||
      !t.gradSlabLengthMPerM.has_value() ||
      !t.gradTwoDLengthMPerM.has_value()) return false;
  const double slab = *t.slabVarianceT2;
  const double twoD = *t.twoDVarianceT2;
  const double total = slab + twoD;
  const double l2 = *t.twoDBendoverLengthM;
  const auto& field = *local.meanFieldT;
  const auto& jacobian = *local.dMeanFieldDxTPerM;
  for (int j = 0; j < 3; ++j) {
    double gradB = 0.0;
    for (int i = 0; i < 3; ++i) {
      if (!Finite(jacobian[i][j])) return false;
      gradB += field[i] * jacobian[i][j] / s.fieldT;
    }
    const double gradSlab = (*t.gradSlabVarianceT2PerM)[j];
    const double gradTwoD = (*t.gradTwoDVarianceT2PerM)[j];
    const double gradLs = (*t.gradSlabLengthMPerM)[j];
    const double gradL2 = (*t.gradTwoDLengthMPerM)[j];
    if (!Finite(gradSlab) || !Finite(gradTwoD) || !Finite(gradLs) ||
        !Finite(gradL2)) return false;
    const double gradLnB = gradB / s.fieldT;
    const double gradLnSlab = gradSlab / slab;
    const double gradLnTotal = (gradSlab + gradTwoD) / total;
    const double gradLnLs = gradLs / s.slabLengthM;
    const double gradLnL2 = gradL2 / l2;
    const double gradX[4] = {
        -gradLnB - gradLnLs,
        gradLnSlab - gradLnTotal,
        gradLnTotal - 2.0 * gradLnB,
        gradLnLs - gradLnL2};
    (*gradLogLambda)[j] = gradLnLs;
    for (int a = 0; a < 4; ++a)
      (*gradLogLambda)[j] +=
          PolyDerivative(polynomial, x, a) * gradX[a];
  }
  return true;
}

bool NlgceBox(const TurbulenceRatios& s) {
  return s.r >= 1.0e-5 && s.r <= 6.3 && s.fs >= 1.0e-3 && s.fs <= 0.85 &&
      s.eps2 >= 1.0e-4 && s.eps2 <= 1.0e2 &&
      s.rho >= 1.0 && s.rho <= 1.0e3;
}

ParallelResult EvaluateNlgceF(const ParticleState& particle,
                              const LocalState& local,
                              const ModelConfiguration& configuration) {
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  TurbulenceRatios s;
  status = Ratios(particle, local, &s);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  if (!NlgceBox(s))
    return Failure(StatusCode::OutsideModelDomain,
                   "NLGCE-F ratios are outside the published fit box",
                   configuration, local);
  const double x[4] = {std::log(s.r), std::log(s.fs), std::log(s.eps2),
                       std::log(s.rho)};
  const double lambdaParallel = s.slabLengthM * std::exp(Poly(NlgceFParallel, x));
  const double lambdaPerp = s.slabLengthM * std::exp(Poly(NlgceFPerpendicular, x));
  ParallelResult result = Success(configuration, local, kin, lambdaParallel);
  if (result.status.ok() && Positive(lambdaPerp)) {
    result.perpendicular = ParallelResult::PerpendicularPair{
        kin.speedMPerS * lambdaPerp / 3.0, lambdaPerp};
    result.dLnLambdaDLnRigidity = PolyDerivative(NlgceFParallel, x, 0);
    result.dLnKappaDLnRigidity = *result.dLnLambdaDLnRigidity +
        1.0 / (kin.gamma * kin.gamma);
    result.diagnosticMask |= SurrogateErrorUnbounded;
    std::array<double, 3> gradLogLambda;
    if (NlgceSpatialGradient(local, s, NlgceFParallel, x,
                             &gradLogLambda)) {
      std::array<double, 3> gradient;
      for (int j = 0; j < 3; ++j)
        gradient[j] = *result.kappaParallelM2PerS * gradLogLambda[j];
      result.gradKappaParallelMPerS = gradient;
    } else {
      result.diagnosticMask |= DerivativeUnavailable;
    }
    result.provenance.coefficientSetIdentity = NlgceCoefficientDigest;
  }
  return result;
}

double IntegrateLogSpectrum(const std::function<double(double)>& f,
                            double alpha, double beta,
                            const NumericalParameters& numerical,
                            bool resonant, double omega) {
  std::vector<double> points{-40.0, 0.0, 40.0};
  if (alpha > 0.0 && beta > 0.0) {
    const double p = 0.5 * std::log(alpha / beta);
    if (p > -40.0 && p < 40.0) points.push_back(p);
  }
  if (resonant && omega > 0.0 && beta > 0.0) {
    const double p = 0.5 * std::log(omega / beta);
    if (p > -40.0 && p < 40.0) points.push_back(p);
  }
  const Integral result = IntegratePieces(f, points,
      std::min(numerical.relativeTolerance, 1.0e-10),
      numerical.maximumRefinements);
  return result.ok ? result.value : std::numeric_limits<double>::quiet_NaN();
}

struct ClosureResidual {
  TurbulenceRatios s;
  NumericalParameters numerical;
  double kappaXGiven = 0.0;  // dimensionless kappa_x/(v ell_s)
  ModelId model = ModelId::NlgceN;

  bool Evaluate(double yz, double yx, double* rz, double* rx) const {
    if (!Finite(yz) || !Finite(yx) || std::fabs(yz) > 300.0 ||
        std::fabs(yx) > 300.0) return false;
    const double kz = std::exp(yz);
    const double kx = model == ModelId::NlpaGivenPerp
        ? kappaXGiven : std::exp(yx);
    const double omega = 1.0 / s.r;
    const double l2 = 1.0 / s.rho;
    const double slab = s.fs * s.eps2;
    const double twoD = (1.0 - s.fs) * s.eps2;
    const double alpha = 1.0 / (3.0 * kz);
    auto resonantIntegral = [&](double beta) {
      const auto f = [&](double u) {
        const double x = std::exp(u);
        const double A = alpha + beta * x * x;
        return 4.0 * Cnu(5.0 / 6.0) *
            std::pow(1.0 + x * x, -5.0 / 6.0) *
            A / (omega * omega + A * A) * x;
      };
      return IntegrateLogSpectrum(f, alpha, beta, numerical, true, omega);
    };
    auto inverseIntegral = [&](double beta) {
      const auto f = [&](double u) {
        const double x = std::exp(u);
        return 4.0 * Cnu(5.0 / 6.0) *
            std::pow(1.0 + x * x, -5.0 / 6.0) /
            (alpha + beta * x * x) * x;
      };
      return IntegrateLogSpectrum(f, alpha, beta, numerical, false, 0.0);
    };
    const double rs = resonantIntegral(kz);
    const double r2 = resonantIntegral(kx / (l2 * l2));
    if (!Positive(rs) || !Positive(r2)) return false;
    const double rhsZ = 3.0 * Ax(s) * omega * omega *
        (slab * rs + twoD * r2);
    if (!Positive(rhsZ)) return false;
    *rz = yz + std::log(rhsZ);
    if (model == ModelId::NlpaGivenPerp) { *rx = 0.0; return true; }
    const double j2 = inverseIntegral(kx / (l2 * l2));
    if (!Positive(j2)) return false;
    double rhsX = 0.0;
    if (model == ModelId::NlgceN) {
      rhsX = APrime2(s) * twoD * j2 / 6.0;
    } else {
      const double js = inverseIntegral(kz);
      if (!Positive(js)) return false;
      rhsX = (1.0 / 3.0) * (slab * js + twoD * j2) / 6.0;
    }
    if (!Positive(rhsX)) return false;
    *rx = yx - std::log(rhsX);
    return true;
  }
};

bool SolveClosure(const ClosureResidual& residual, double* yz, double* yx,
                  int* iterations, double* maximumResidual) {
  const int dimensions = residual.model == ModelId::NlpaGivenPerp ? 1 : 2;
  for (int iteration = 0; iteration < residual.numerical.maximumIterations;
       ++iteration) {
    double f0 = 0.0, f1 = 0.0;
    if (!residual.Evaluate(*yz, *yx, &f0, &f1)) return false;
    const double norm = dimensions == 1 ? std::fabs(f0)
                                        : std::max(std::fabs(f0), std::fabs(f1));
    if (norm < NonlinearResidualLimit) {
      *iterations = iteration;
      *maximumResidual = norm;
      return true;
    }
    const double h = 1.0e-5;
    double zp0 = 0.0, zp1 = 0.0, xp0 = 0.0, xp1 = 0.0;
    if (!residual.Evaluate(*yz + h, *yx, &zp0, &zp1)) return false;
    double dz = -f0 / ((zp0 - f0) / h);
    double dx = 0.0;
    if (dimensions == 2) {
      if (!residual.Evaluate(*yz, *yx + h, &xp0, &xp1)) return false;
      const double a = (zp0 - f0) / h, b = (xp0 - f0) / h;
      const double c = (zp1 - f1) / h, d = (xp1 - f1) / h;
      const double determinant = a * d - b * c;
      if (!Finite(determinant) || std::fabs(determinant) < 1.0e-14) return false;
      dz = (-f0 * d + b * f1) / determinant;
      dx = (c * f0 - a * f1) / determinant;
    }
    bool accepted = false;
    for (int line = 0; line < 20; ++line) {
      const double scale = std::ldexp(1.0, -line);
      double trial0 = 0.0, trial1 = 0.0;
      if (residual.Evaluate(*yz + scale * dz, *yx + scale * dx,
                            &trial0, &trial1)) {
        const double trialNorm = dimensions == 1 ? std::fabs(trial0)
            : std::max(std::fabs(trial0), std::fabs(trial1));
        if (trialNorm < norm) {
          *yz += scale * dz;
          *yx += scale * dx;
          accepted = true;
          break;
        }
      }
    }
    if (!accepted) return false;
  }
  return false;
}

ParallelResult EvaluateNonlinear(const ParticleState& particle,
                                 const LocalState& local,
                                 const ModelConfiguration& configuration) {
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  TurbulenceRatios s;
  status = Ratios(particle, local, &s);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  ClosureResidual residual;
  residual.s = s;
  residual.numerical = configuration.nonlinear.numerical;
  residual.model = configuration.model;
  if (configuration.model == ModelId::NlpaGivenPerp) {
    if (!local.suppliedKappaPerpendicularM2PerS.has_value() ||
        !Positive(*local.suppliedKappaPerpendicularM2PerS) ||
        local.suppliedPerpendicularModelId.empty())
      return Failure(StatusCode::MissingInput,
          "nlpa_given_perp requires positive supplied kappa_perp and its model identity",
          configuration, local);
    residual.kappaXGiven = *local.suppliedKappaPerpendicularM2PerS /
        (kin.speedMPerS * s.slabLengthM);
  }
  double yz = std::log(1.0 / 3.0), yx = std::log(0.01 / 3.0);
  if (NlgceBox(s)) {
    const double x[4] = {std::log(s.r), std::log(s.fs), std::log(s.eps2),
                         std::log(s.rho)};
    yz = Poly(NlgceFParallel, x) - std::log(3.0);
    yx = Poly(NlgceFPerpendicular, x) - std::log(3.0);
  }
  if (configuration.model == ModelId::NlpaGivenPerp)
    yx = std::log(residual.kappaXGiven);
  int iterations = 0;
  double maxResidual = 0.0;
  if (!SolveClosure(residual, &yz, &yx, &iterations, &maxResidual))
    return Failure(StatusCode::NonlinearSolverFailed,
                   "coupled logarithmic closure did not meet the 1e-8 residual gate",
                   configuration, local);
  const double lambdaParallel = 3.0 * std::exp(yz) * s.slabLengthM;
  const double lambdaPerp = 3.0 * std::exp(yx) * s.slabLengthM;
  ParallelResult result = Success(configuration, local, kin, lambdaParallel);
  if (result.status.ok()) {
    result.perpendicular = ParallelResult::PerpendicularPair{
        kin.speedMPerS * lambdaPerp / 3.0, lambdaPerp};
    result.nonlinearIterations = iterations;
    result.nonlinearMaxLogResidual = maxResidual;
    result.diagnosticMask |= DerivativeUnavailable;
  }
  return result;
}

Status BroadenedDmu(double mu, const ParticleState& particle,
                    const ParticleKinematics& kin, const LocalState& local,
                    const BroadenedSlabParameters& p, double* dMuMu) {
  if (!dMuMu || !Finite(mu) || std::fabs(mu) > 1.0)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "pitch angle mu must be finite and in [-1,1]");
  if (p.width0PerS == 0.0) {
    QltSlabParameters qlt;
    qlt.spectrum = p.spectrum;
    qlt.numerical = p.numerical;
    return QltDmu(mu, particle, kin, local, qlt, dMuMu);
  }
  double field = 0.0, variance = 0.0, length = 0.0, index = 0.0;
  Status status = SlabState(local, &field, &variance, &length, &index,
      p.spectrum.form != SpectrumForm::Multirange);
  if (!status.ok()) return status;
  if (mu == 1.0 || mu == -1.0) { *dMuMu = 0.0; return Status::Success(); }
  const double omega = std::fabs(particle.chargeC) * field /
      (kin.gamma * particle.massKg);
  const double nu = 0.5 * index;
  Status spectrumStatus = Status::Success();
  const auto ku = [&](double u) {
      const double x = std::exp(u);
      const double k = x / length;
      const double width = p.kernel == BroadeningKernel::LorentzianLinear
          ? p.width0PerS + p.decorrelationSpeedMPerS * k : p.width0PerS;
      if (!Positive(width)) return std::numeric_limits<double>::quiet_NaN();
      auto kernel = [&](double mismatch) {
        if (p.kernel == BroadeningKernel::Gaussian)
          return std::sqrt(std::acos(-1.0)) / width *
              std::exp(-mismatch * mismatch / (width * width));
        return width / (mismatch * mismatch + width * width);
      };
      double pDk = 0.0;
      if (p.spectrum.form == SpectrumForm::SmoothBendover) {
        pDk = 4.0 * Cnu(nu) * variance *
            std::pow(1.0 + x * x, -nu) * x;
      } else {
        double power = 0.0;
        const Status spectrum = p.spectrum.form == SpectrumForm::Multirange
            ? MultirangeSpectrumValue(p.spectrum, variance, length, index,
                                      k, &power)
            : SpectrumValue(p.spectrum, k, &power);
        if (!spectrum.ok()) { spectrumStatus = spectrum; return 0.0; }
        pDk = power * x / length;
      }
      return pDk * (kernel(k * kin.speedMPerS * mu - omega) +
                    kernel(k * kin.speedMPerS * mu + omega));
  };
  // The upper log-wavenumber bound is a numerical refinement interval, not a
  // physical cutoff.  mu=z^3 in the outer integral can place a narrow smooth-
  // spectrum resonance far above exp(50)/ell; revision-1.4 fixtures establish
  // convergence on [-50,150].
  std::vector<double> points{-50.0, 0.0, 150.0};
  if (mu != 0.0) {
      const double center = omega * length /
          (kin.speedMPerS * std::fabs(mu));
      if (Positive(center)) {
        // Resolve the resonance in its own frequency-width scale.  For a
        // constant kernel, delta-k/k_res = width/Omega; powers of two cover
        // the narrow core and both algebraic/Gaussian shoulders without a
        // fixed grid that can step over the peak (Section 9.5/PD07).
        const double relativeWidth = p.width0PerS / omega;
        for (int m = 0; m < 40; ++m) {
          const double offset = std::ldexp(relativeWidth, m);
          for (double sign : {-1.0, 1.0}) {
            const double shifted = center * (1.0 + sign * offset);
            if (Positive(shifted)) {
              const double point = std::log(shifted);
              if (point > -50.0 && point < 150.0) points.push_back(point);
            }
          }
          if (offset > 1.0e6) break;
        }
      }
  }
  const Integral integral = IntegratePieces(ku, points,
      std::min(p.numerical.relativeTolerance, 1.0e-9),
      p.numerical.maximumRefinements);
  if (!spectrumStatus.ok()) return spectrumStatus;
  if (!integral.ok || !Positive(integral.value))
    return Status::Error(StatusCode::IntegrationFailed,
                         "broadened wavenumber quadrature did not converge");
  *dMuMu = omega * omega * (1.0 - mu * mu) * integral.value /
      (4.0 * field * field);
  return Positive(*dMuMu)
      ? Status::Success()
      : Status::Error(StatusCode::OutsideModelDomain,
                      "broadened pitch-angle coefficient is not representable");
}

ParallelResult EvaluateBroadened(const ParticleState& particle,
                                 const LocalState& local,
                                 const ModelConfiguration& configuration) {
  const BroadenedSlabParameters& p = configuration.broadenedSlab;
  Status normalization = ValidateSpectrumNormalization(
      p.spectrum, p.numerical.relativeTolerance);
  if (!normalization.ok())
    return Failure(normalization.code, normalization.detail,
                   configuration, local);
  if (p.width0PerS == 0.0) {
    ModelConfiguration qlt = configuration;
    qlt.model = ModelId::QltSlabSpectrum;
    qlt.qltSlab.spectrum = p.spectrum;
    qlt.qltSlab.numerical = p.numerical;
    ParallelResult result = EvaluateQlt(particle, local, qlt, false);
    result.provenance.requestedModelId = "broadened_slab";
    result.provenance.evaluatedModelId = "broadened_slab";
    result.provenance.configurationFingerprint =
        ConfigurationFingerprint(configuration);
    return result;
  }
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  Status integrationStatus = Status::Success();
  const auto pitch = [&](double z) {
    const double mu = z * z * z;
    double d = 0.0;
    const Status evaluated = BroadenedDmu(mu, particle, kin, local, p, &d);
    if (!evaluated.ok()) { integrationStatus = evaluated; return 0.0; }
    return d > 0.0 ? 3.0 * z * z *
        std::pow(1.0 - mu * mu, 2) / d : 0.0;
  };
  std::vector<double> pitchEdges{0.0};
  for (int i = 0; i < 90; ++i)
    pitchEdges.push_back(std::exp(std::log(1.0e-9) * (1.0 - i / 89.0)));
  bool coarseOk = true, fineOk = true;
  double coarse = 0.0, fine = 0.0;
  for (std::size_t i = 1; i < pitchEdges.size(); ++i) {
    coarse += GaussLegendre(pitch, pitchEdges[i - 1], pitchEdges[i],
                            24, &coarseOk);
    fine += GaussLegendre(pitch, pitchEdges[i - 1], pitchEdges[i],
                          32, &fineOk);
    if (!integrationStatus.ok()) break;
  }
  const bool pitchConverged = coarseOk && fineOk && Positive(fine) &&
      std::fabs(fine - coarse) <=
          p.numerical.relativeTolerance * std::fabs(fine);
  if (!integrationStatus.ok() || !pitchConverged)
    return Failure(integrationStatus.ok() ? StatusCode::IntegrationFailed
                                          : integrationStatus.code,
                   integrationStatus.ok()
                       ? "nested broadened-slab quadrature did not converge"
                       : integrationStatus.detail,
                   configuration, local);
  ParallelResult result = Success(configuration, local, kin,
      3.0 * kin.speedMPerS * fine / 4.0);
  if (result.status.ok()) result.diagnosticMask |= DerivativeUnavailable;
  return result;
}

double TableCoordinate(TableAxis axis, const ParticleState& particle,
                       const ParticleKinematics& kin, const LocalState& local,
                       Status* status) {
  switch (axis) {
    case TableAxis::Rigidity: return kin.rigidityV;
    case TableAxis::TotalKineticEnergy: return kin.kineticEnergyJ;
    case TableAxis::EnergyPerNucleon:
      if (!particle.nucleonCount.has_value() || !Positive(*particle.nucleonCount)) {
        *status = Status::Error(StatusCode::MissingInput,
            "energy-per-nucleon table axis requires nucleon_count");
        return 0.0;
      }
      return kin.kineticEnergyJ / *particle.nucleonCount;
    case TableAxis::Speed: return kin.speedMPerS;
    case TableAxis::HeliocentricRadius: return Magnitude(local.positionM);
    case TableAxis::Time: return local.timeS;
    case TableAxis::MeanFieldMagnitude:
      if (!local.meanFieldT.has_value()) {
        *status = Status::Error(StatusCode::MissingInput,
                               "mean-field table axis requires mean_B_T");
        return 0.0;
      }
      return Magnitude(*local.meanFieldT);
  }
  return 0.0;
}

ParallelResult EvaluateTable(const ParticleState& particle,
                             const LocalState& local,
                             const ModelConfiguration& configuration) {
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
  const TabulatedParallelParameters& table = configuration.table;
  const std::size_t n = table.axes.size();
  std::vector<std::size_t> lower(n);
  std::vector<double> fraction(n);
  for (std::size_t a = 0; a < n; ++a) {
    const double value = TableCoordinate(table.axes[a], particle, kin, local, &status);
    if (!status.ok()) return Failure(status.code, status.detail, configuration, local);
    const std::vector<double>& axis = table.axisSI[a];
    if (value < axis.front() || value > axis.back())
      return Failure(StatusCode::OutsideModelDomain,
                     "table coordinate is outside its closed data domain",
                     configuration, local);
    auto upper = std::upper_bound(axis.begin(), axis.end(), value);
    std::size_t lo = upper == axis.begin() ? 0 :
        static_cast<std::size_t>(upper - axis.begin() - 1);
    if (lo + 1 >= axis.size()) lo = axis.size() - 2;
    lower[a] = lo;
    if (table.axes[a] == TableAxis::Time)
      fraction[a] = (value - axis[lo]) / (axis[lo + 1] - axis[lo]);
    else
      fraction[a] = std::log(value / axis[lo]) /
          std::log(axis[lo + 1] / axis[lo]);
  }
  std::size_t timeAxis = n;
  for (std::size_t a = 0; a < n; ++a)
    if (table.axes[a] == TableAxis::Time) timeAxis = a;
  const auto interpolatePositiveAxes = [&](std::size_t selectedTime) {
    const std::size_t dimensions = n - (timeAxis < n ? 1 : 0);
    const std::size_t corners = std::size_t(1) << dimensions;
    double logCoefficient = 0.0;
    for (std::size_t corner = 0; corner < corners; ++corner) {
      double weight = 1.0;
      std::size_t flat = 0;
      std::size_t bit = 0;
      for (std::size_t a = 0; a < n; ++a) {
        std::size_t index = selectedTime;
        if (a != timeAxis) {
          const bool high = (corner & (std::size_t(1) << bit)) != 0;
          weight *= high ? fraction[a] : 1.0 - fraction[a];
          index = lower[a] + (high ? 1 : 0);
          ++bit;
        }
        flat = flat * table.axisSI[a].size() + index;
      }
      logCoefficient += weight * std::log(table.coefficientSI[flat]);
    }
    return std::exp(logCoefficient);
  };
  double coefficient = 0.0;
  if (timeAxis == n) {
    coefficient = interpolatePositiveAxes(0);
  } else {
    const double lowerValue = interpolatePositiveAxes(lower[timeAxis]);
    if (*table.timeInterpolation == TimeInterpolation::StepPrevious) {
      coefficient = lowerValue;
    } else {
      const double upperValue = interpolatePositiveAxes(lower[timeAxis] + 1);
      coefficient = (1.0 - fraction[timeAxis]) * lowerValue +
                    fraction[timeAxis] * upperValue;
    }
  }
  const double lambda = table.storedCoefficient == StoredCoefficient::LambdaParallel
      ? coefficient : 3.0 * coefficient / kin.speedMPerS;
  ParallelResult result = Success(configuration, local, kin, lambda);
  if (result.status.ok()) {
    result.diagnosticMask |= DerivativeUnavailable;
    result.provenance.coefficientSetIdentity = table.generationIdentity;
  }
  return result;
}

bool ParseDouble(const std::string& text, double* value) {
  if (!value || text.empty()) return false;
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' || !Finite(parsed))
    return false;
  *value = parsed;
  return true;
}

Status Parameters(const std::vector<InputParameter>& input,
                  std::map<std::string, std::string>* values) {
  for (const InputParameter& parameter : input) {
    if (parameter.name.empty() ||
        !values->emplace(parameter.name, parameter.value).second)
      return Status::Error(StatusCode::InvalidConfiguration,
          "empty or duplicate parallel-diffusion parameter '" + parameter.name + "'");
  }
  return Status::Success();
}

Status Take(std::map<std::string, std::string>* values, const std::string& key,
            bool required, std::string* result) {
  const auto it = values->find(key);
  if (it == values->end())
    return required ? Status::Error(StatusCode::MissingInput,
        "missing parameter '" + key + "'") : Status::Success();
  *result = it->second;
  values->erase(it);
  return Status::Success();
}

Status TakeDouble(std::map<std::string, std::string>* values,
                  const std::string& key, bool required, double* result) {
  std::string text;
  Status status = Take(values, key, required, &text);
  if (!status.ok() || (!required && text.empty())) return status;
  return ParseDouble(text, result) ? Status::Success()
      : Status::Error(StatusCode::InvalidConfiguration,
                      "parameter '" + key + "' is not a finite number");
}

Status ParseList(const std::string& text, std::vector<double>* values) {
  std::istringstream stream(text);
  std::string item;
  while (std::getline(stream, item, ',')) {
    double value = 0.0;
    if (!ParseDouble(item, &value))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "numeric list contains an invalid item");
    values->push_back(value);
  }
  return values->empty()
      ? Status::Error(StatusCode::InvalidConfiguration, "numeric list is empty")
      : Status::Success();
}

Status ParseNumerical(std::map<std::string, std::string>* values,
                      NumericalParameters* numerical) {
  Status status = TakeDouble(values, "relative_tolerance", false,
                             &numerical->relativeTolerance);
  if (!status.ok()) return status;
  std::string text;
  status = Take(values, "maximum_refinements", false, &text);
  if (!status.ok()) return status;
  if (!text.empty()) {
    double number = 0.0;
    if (!ParseDouble(text, &number) || number != std::floor(number))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "maximum_refinements must be an integer");
    numerical->maximumRefinements = static_cast<int>(number);
  }
  status = Take(values, "maximum_iterations", false, &text);
  if (!status.ok()) return status;
  if (!text.empty()) {
    double number = 0.0;
    if (!ParseDouble(text, &number) || number != std::floor(number))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "maximum_iterations must be an integer");
    numerical->maximumIterations = static_cast<int>(number);
  }
  return Status::Success();
}

Status ParseSpectrum(std::map<std::string, std::string>* values,
                     SpectrumParameters* spectrum) {
  std::string form;
  Status status = Take(values, "spectrum_form", true, &form);
  if (!status.ok()) return status;
  form = Lower(form);
  if (form == "smooth_bendover") {
    spectrum->form = SpectrumForm::SmoothBendover;
    return Status::Success();
  }
  if (form == "multirange") {
    spectrum->form = SpectrumForm::Multirange;
    status = TakeDouble(values, "energy_range_index", true,
                        &spectrum->energyRangeIndex);
    if (status.ok())
      status = TakeDouble(values, "dissipation_index", true,
                          &spectrum->dissipationIndex);
    if (status.ok())
      status = TakeDouble(values, "dissipation_wavenumber_rad_per_m", true,
                          &spectrum->dissipationWavenumberRadPerM);
    return status;
  }
  if (form != "supplied_log_log")
    return Status::Error(StatusCode::InvalidConfiguration,
        "spectrum_form must be smooth_bendover, multirange, or supplied_log_log");
  spectrum->form = SpectrumForm::SuppliedLogLog;
  std::string k, p;
  status = Take(values, "spectrum_k_rad_per_m", true, &k);
  if (!status.ok()) return status;
  status = Take(values, "spectrum_power_T2_m", true, &p);
  if (!status.ok()) return status;
  status = Take(values, "spectrum_source_identity", true,
                &spectrum->sourceIdentity);
  if (!status.ok()) return status;
  status = TakeDouble(values, "spectrum_declared_variance_T2", true,
                      &spectrum->declaredVarianceT2);
  if (!status.ok()) return status;
  status = ParseList(k, &spectrum->wavenumberRadPerM);
  if (!status.ok()) return status;
  status = ParseList(p, &spectrum->powerT2M);
  if (!status.ok()) return status;
  auto policy = [&](const char* key, SpectrumTailPolicy* target,
                    double* exponent) -> Status {
    std::string value;
    Status s = Take(values, key, true, &value);
    if (!s.ok()) return s;
    value = Lower(value);
    if (value == "out_of_domain") *target = SpectrumTailPolicy::OutOfDomain;
    else if (value == "zero") *target = SpectrumTailPolicy::Zero;
    else if (value == "power_law") {
      *target = SpectrumTailPolicy::PowerLaw;
      s = TakeDouble(values, std::string(key) + "_index", true, exponent);
      if (!s.ok()) return s;
    } else return Status::Error(StatusCode::InvalidConfiguration,
        std::string(key) + " must be out_of_domain, zero, or power_law");
    return Status::Success();
  };
  status = policy("low_k_policy", &spectrum->lowKPolicy,
                  &spectrum->lowKPowerIndex);
  if (!status.ok()) return status;
  return policy("high_k_policy", &spectrum->highKPolicy,
                &spectrum->highKPowerIndex);
}

Status ParseBroadened(std::map<std::string, std::string>* values,
                      BroadenedSlabParameters* p) {
  Status status = ParseSpectrum(values, &p->spectrum);
  if (!status.ok()) return status;
  std::string kernel;
  status = Take(values, "kernel", true, &kernel);
  if (!status.ok()) return status;
  kernel = Lower(kernel);
  if (kernel == "lorentzian_constant")
    p->kernel = BroadeningKernel::LorentzianConstant;
  else if (kernel == "lorentzian_linear")
    p->kernel = BroadeningKernel::LorentzianLinear;
  else if (kernel == "gaussian")
    p->kernel = BroadeningKernel::Gaussian;
  else
    return Status::Error(StatusCode::InvalidConfiguration,
                         "unknown broadening kernel");
  status = TakeDouble(values, "width0_per_s", true, &p->width0PerS);
  if (status.ok() && p->kernel == BroadeningKernel::LorentzianLinear)
    status = TakeDouble(values, "decorrelation_speed_m_per_s", true,
                        &p->decorrelationSpeedMPerS);
  if (status.ok()) status = ParseNumerical(values, &p->numerical);
  return status;
}

Status ParseClosure(const std::string& value, AdapterClosure* closure) {
  const std::string name = Lower(value);
  if (name == "qlt_slab_spectrum") *closure = AdapterClosure::QltSlabSpectrum;
  else if (name == "qlt_slab_inertial") *closure = AdapterClosure::QltSlabInertial;
  else if (name == "broadened_slab") *closure = AdapterClosure::BroadenedSlab;
  else return Status::Error(StatusCode::InvalidConfiguration,
      "underlying_closure must be qlt_slab_spectrum, qlt_slab_inertial, or broadened_slab");
  return Status::Success();
}

Status ParseTableAxis(const std::string& value, TableAxis* axis) {
  const std::string name = Lower(value);
  if (name == "rigidity") *axis = TableAxis::Rigidity;
  else if (name == "total_kinetic_energy") *axis = TableAxis::TotalKineticEnergy;
  else if (name == "energy_per_nucleon") *axis = TableAxis::EnergyPerNucleon;
  else if (name == "speed") *axis = TableAxis::Speed;
  else if (name == "heliocentric_radius") *axis = TableAxis::HeliocentricRadius;
  else if (name == "time") *axis = TableAxis::Time;
  else if (name == "mean_field_magnitude") *axis = TableAxis::MeanFieldMagnitude;
  else return Status::Error(StatusCode::InvalidConfiguration,
                            "unknown table axis '" + value + "'");
  return Status::Success();
}

}  // namespace

bool IsAdvancedModel(ModelId model) {
  return model != ModelId::ConstantLambda && model != ModelId::ConstantKappa &&
      model != ModelId::PowerLawLambda &&
      model != ModelId::BrokenRigidityKappa && model != ModelId::Bohm;
}

Status ValidateAdvancedConfiguration(const ModelConfiguration& c) {
  auto numerical = [](const NumericalParameters& p) {
    return Positive(p.relativeTolerance) && p.relativeTolerance <= 1.0e-2 &&
           p.maximumRefinements >= 8 && p.maximumRefinements <= 30 &&
           p.maximumIterations > 0;
  };
  auto spectrum = [](const SpectrumParameters& p) {
    if (p.form == SpectrumForm::SmoothBendover) return true;
    if (p.form == SpectrumForm::Multirange)
      return Finite(p.energyRangeIndex) && p.energyRangeIndex > -1.0 &&
             Positive(p.dissipationIndex) && p.dissipationIndex > 1.0 &&
             Positive(p.dissipationWavenumberRadPerM);
    if (p.wavenumberRadPerM.size() < 2 ||
        p.wavenumberRadPerM.size() != p.powerT2M.size() ||
        p.sourceIdentity.empty() || !Positive(p.declaredVarianceT2)) return false;
    for (std::size_t i = 0; i < p.wavenumberRadPerM.size(); ++i) {
      if (!Positive(p.wavenumberRadPerM[i]) || !Finite(p.powerT2M[i]) ||
          p.powerT2M[i] < 0.0 ||
          (i && !(p.wavenumberRadPerM[i] > p.wavenumberRadPerM[i - 1])))
        return false;
    }
    return true;
  };
  switch (c.model) {
    case ModelId::QltSlabSpectrum:
      return spectrum(c.qltSlab.spectrum) && numerical(c.qltSlab.numerical)
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "invalid qlt_slab_spectrum configuration");
    case ModelId::QltSlabInertial: return Status::Success();
    case ModelId::PrescribedLambdaMuShape: {
      const auto& p = c.prescribedLambdaMu;
      const bool amplitude = p.amplitudeMode == PitchAngleAmplitudeMode::TargetLambda
          ? Positive(p.targetLambdaM) : Positive(p.fixedD0PerS);
      return Finite(p.qMu) && Finite(p.hMu) && p.hMu >= 0.0 && amplitude &&
             numerical(p.numerical)
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "invalid prescribed_lambda_mu_shape configuration");
    }
    case ModelId::BroadenedSlab: {
      const auto& p = c.broadenedSlab;
      const bool width = Finite(p.width0PerS) && p.width0PerS >= 0.0 &&
          Finite(p.decorrelationSpeedMPerS) &&
          p.decorrelationSpeedMPerS >= 0.0;
      return spectrum(p.spectrum) && numerical(p.numerical) && width
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "invalid broadened_slab configuration");
    }
    case ModelId::NlpaGivenPerp:
    case ModelId::NlgcE:
    case ModelId::NlgceN:
      return numerical(c.nonlinear.numerical) &&
             c.nonlinear.numerical.relativeTolerance <= 1.0e-8
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "nonlinear closure requires relative_tolerance<=1e-8");
    case ModelId::NlgceF2014:
      return c.nlgceF.coefficientSet == "Qin_Zhang_2014_Tables_3_4"
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "nlgce_f_2014 requires the fixed Qin_Zhang_2014_Tables_3_4 coefficient set");
    case ModelId::TurbulenceAdapter:
      return Positive(c.turbulenceAdapter.vacuumPermeabilityHPerM) &&
             (c.turbulenceAdapter.momentSource == MomentSource::LocalState ||
              (c.turbulenceAdapter.providerMomentM2PerS2.has_value() &&
               c.turbulenceAdapter.residualEnergy.has_value() &&
               c.turbulenceAdapter.slabFraction.has_value() &&
               Positive(*c.turbulenceAdapter.providerMomentM2PerS2) &&
               Finite(*c.turbulenceAdapter.residualEnergy) &&
               *c.turbulenceAdapter.residualEnergy >= -1.0 &&
               *c.turbulenceAdapter.residualEnergy <= 1.0 &&
               *c.turbulenceAdapter.slabFraction >= 0.0 &&
               *c.turbulenceAdapter.slabFraction <= 1.0)) &&
             (c.turbulenceAdapter.closure == AdapterClosure::QltSlabInertial ||
              (c.turbulenceAdapter.closure == AdapterClosure::QltSlabSpectrum &&
               spectrum(c.turbulenceAdapter.qlt.spectrum) &&
               numerical(c.turbulenceAdapter.qlt.numerical)) ||
              (c.turbulenceAdapter.closure == AdapterClosure::BroadenedSlab &&
               spectrum(c.turbulenceAdapter.broadened.spectrum) &&
               numerical(c.turbulenceAdapter.broadened.numerical) &&
               Finite(c.turbulenceAdapter.broadened.width0PerS) &&
               c.turbulenceAdapter.broadened.width0PerS >= 0.0))
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "invalid turbulence_adapter moment or slab fraction");
    case ModelId::WaveSpectrumAdapter:
      return c.waveSpectrumAdapter.propagation == "balanced" &&
             c.waveSpectrumAdapter.polarization == "transverse_axisymmetric" &&
             c.waveSpectrumAdapter.frame == "plasma" &&
             (c.waveSpectrumAdapter.closure == AdapterClosure::QltSlabInertial ||
              (c.waveSpectrumAdapter.closure == AdapterClosure::QltSlabSpectrum &&
               spectrum(c.waveSpectrumAdapter.qlt.spectrum) &&
               numerical(c.waveSpectrumAdapter.qlt.numerical)) ||
              (c.waveSpectrumAdapter.closure == AdapterClosure::BroadenedSlab &&
               spectrum(c.waveSpectrumAdapter.broadened.spectrum) &&
               numerical(c.waveSpectrumAdapter.broadened.numerical) &&
               Finite(c.waveSpectrumAdapter.broadened.width0PerS) &&
               c.waveSpectrumAdapter.broadened.width0PerS >= 0.0))
          ? Status::Success() : Status::Error(StatusCode::InvalidConfiguration,
              "revision 1.4 wave adapter supports only balanced transverse-axisymmetric plasma-frame input");
    case ModelId::TabulatedParallel: {
      const auto& t = c.table;
      if (t.axes.empty() || t.axes.size() != t.axisSI.size() ||
          t.axes.size() > 8 || t.generationIdentity.empty())
        return Status::Error(StatusCode::InvalidConfiguration,
                             "invalid table axes or generation identity");
      std::size_t count = 1;
      for (std::size_t a = 0; a < t.axisSI.size(); ++a) {
        const auto& axis = t.axisSI[a];
        if (axis.size() < 2) return Status::Error(StatusCode::InvalidConfiguration,
                                                  "each table axis needs at least two points");
        for (std::size_t i = 0; i < axis.size(); ++i)
          if (!Finite(axis[i]) ||
              (t.axes[a] != TableAxis::Time && axis[i] <= 0.0) ||
              (i && !(axis[i] > axis[i - 1])))
            return Status::Error(StatusCode::InvalidConfiguration,
                                 "table axes must be finite and strictly increasing");
        for (std::size_t b = a + 1; b < t.axes.size(); ++b)
          if (t.axes[a] == t.axes[b])
            return Status::Error(StatusCode::InvalidConfiguration,
                                 "table axes must be unique");
        count *= axis.size();
      }
      const bool hasTime = std::find(t.axes.begin(), t.axes.end(),
                                     TableAxis::Time) != t.axes.end();
      if (hasTime != t.timeInterpolation.has_value())
        return Status::Error(StatusCode::InvalidConfiguration,
            "a time table axis requires exactly one explicit time rule");
      if (count != t.coefficientSI.size())
        return Status::Error(StatusCode::InvalidConfiguration,
                             "table value count does not match axis shape");
      for (double value : t.coefficientSI)
        if (!Positive(value)) return Status::Error(StatusCode::InvalidConfiguration,
                                                    "table coefficients must be positive");
      return Status::Success();
    }
    default:
      return Status::Error(StatusCode::UnsupportedModel,
                           "model is not an advanced backend");
  }
}

ParallelResult EvaluateAdvanced(const ParticleState& particle,
                                const LocalState& local,
                                const ModelConfiguration& configuration) {
  switch (configuration.model) {
    case ModelId::QltSlabSpectrum:
      return EvaluateQlt(particle, local, configuration, false);
    case ModelId::QltSlabInertial:
      return EvaluateQlt(particle, local, configuration, true);
    case ModelId::PrescribedLambdaMuShape:
      return EvaluatePrescribed(particle, local, configuration);
    case ModelId::BroadenedSlab:
      return EvaluateBroadened(particle, local, configuration);
    case ModelId::NlpaGivenPerp:
    case ModelId::NlgcE:
    case ModelId::NlgceN:
      return EvaluateNonlinear(particle, local, configuration);
    case ModelId::NlgceF2014:
      return EvaluateNlgceF(particle, local, configuration);
    case ModelId::TabulatedParallel:
      return EvaluateTable(particle, local, configuration);
    case ModelId::TurbulenceAdapter: {
      if (!local.turbulence.has_value() ||
          !local.turbulence->densityKgPerM3.has_value())
        return Failure(StatusCode::MissingInput,
                       "turbulence_adapter requires provider density",
                       configuration, local);
      double total = 0.0;
      const auto& a = configuration.turbulenceAdapter;
      const std::optional<double>& moment =
          a.momentSource == MomentSource::InputParameters
              ? a.providerMomentM2PerS2
              : local.turbulence->providerMomentM2PerS2;
      const std::optional<double>& residual =
          a.momentSource == MomentSource::InputParameters
              ? a.residualEnergy : local.turbulence->residualEnergy;
      const std::optional<double>& slabFraction =
          a.momentSource == MomentSource::InputParameters
              ? a.slabFraction : local.turbulence->slabFraction;
      if (!moment.has_value() || !residual.has_value() ||
          !slabFraction.has_value())
        return Failure(StatusCode::MissingInput,
            "turbulence_adapter moment source lacks moment, residual energy, or slab fraction",
            configuration, local);
      if (!Finite(*slabFraction) || *slabFraction < 0.0 ||
          *slabFraction > 1.0)
        return Failure(StatusCode::InvalidBackground,
                       "turbulence_adapter slab fraction must be in [0,1]",
                       configuration, local);
      Status converted = ConvertTurbulenceMoment(*moment,
          *local.turbulence->densityKgPerM3, *residual,
          a.vacuumPermeabilityHPerM,
          a.energyConvention, a.residualConvention, &total);
      if (!converted.ok())
        return Failure(converted.code, converted.detail, configuration, local);
      LocalState adapted = local;
      adapted.turbulence->slabVarianceT2 = *slabFraction * total;
      adapted.turbulence->twoDVarianceT2 = (1.0 - *slabFraction) * total;
      ModelConfiguration child = configuration;
      if (a.closure == AdapterClosure::QltSlabSpectrum) {
        child.model = ModelId::QltSlabSpectrum; child.qltSlab = a.qlt;
      } else if (a.closure == AdapterClosure::QltSlabInertial) {
        child.model = ModelId::QltSlabInertial;
      } else {
        child.model = ModelId::BroadenedSlab; child.broadenedSlab = a.broadened;
      }
      ParallelResult result = EvaluateAdvanced(particle, adapted, child);
      result.provenance.requestedModelId = "turbulence_adapter";
      result.provenance.configurationFingerprint =
          ConfigurationFingerprint(configuration);
      result.provenance.inputMomentConvention =
          std::to_string(static_cast<int>(a.energyConvention)) + ":" +
          std::to_string(static_cast<int>(a.residualConvention));
      return result;
    }
    case ModelId::WaveSpectrumAdapter: {
      const auto& a = configuration.waveSpectrumAdapter;
      ModelConfiguration child = configuration;
      if (a.closure == AdapterClosure::QltSlabSpectrum) {
        child.model = ModelId::QltSlabSpectrum; child.qltSlab = a.qlt;
      } else if (a.closure == AdapterClosure::QltSlabInertial) {
        child.model = ModelId::QltSlabInertial;
      } else {
        child.model = ModelId::BroadenedSlab; child.broadenedSlab = a.broadened;
      }
      ParallelResult result = EvaluateAdvanced(particle, local, child);
      result.provenance.requestedModelId = "wave_spectrum_adapter";
      result.provenance.configurationFingerprint =
          ConfigurationFingerprint(configuration);
      result.provenance.inputMomentConvention =
          a.propagation + ":" + a.polarization + ":" + a.frame;
      return result;
    }
    default:
      return Failure(StatusCode::UnsupportedModel,
                     "advanced evaluator received a non-advanced model",
                     configuration, local);
  }
}

Status EvaluateAdvancedPitchAngle(double mu, const ParticleState& particle,
                                  const LocalState& local,
                                  const ModelConfiguration& configuration,
                                  double* dMuMuPerS) {
  ParticleKinematics kin;
  Status status = ComputeParticleKinematics(particle, &kin);
  if (!status.ok()) return status;
  if (configuration.model == ModelId::QltSlabSpectrum) {
    status = ValidateSpectrumNormalization(configuration.qltSlab.spectrum,
        configuration.qltSlab.numerical.relativeTolerance);
    if (!status.ok()) return status;
    return QltDmu(mu, particle, kin, local, configuration.qltSlab, dMuMuPerS);
  }
  if (configuration.model == ModelId::BroadenedSlab) {
    status = ValidateSpectrumNormalization(
        configuration.broadenedSlab.spectrum,
        configuration.broadenedSlab.numerical.relativeTolerance);
    if (!status.ok()) return status;
    return BroadenedDmu(mu, particle, kin, local,
                        configuration.broadenedSlab, dMuMuPerS);
  }
  if (configuration.model == ModelId::TurbulenceAdapter) {
    if (!local.turbulence.has_value() ||
        !local.turbulence->densityKgPerM3.has_value())
      return Status::Error(StatusCode::MissingInput,
                           "turbulence_adapter requires provider density");
    const auto& a = configuration.turbulenceAdapter;
    const std::optional<double>& moment =
        a.momentSource == MomentSource::InputParameters
            ? a.providerMomentM2PerS2
            : local.turbulence->providerMomentM2PerS2;
    const std::optional<double>& residual =
        a.momentSource == MomentSource::InputParameters
            ? a.residualEnergy : local.turbulence->residualEnergy;
    const std::optional<double>& slabFraction =
        a.momentSource == MomentSource::InputParameters
            ? a.slabFraction : local.turbulence->slabFraction;
    if (!moment.has_value() || !residual.has_value() ||
        !slabFraction.has_value())
      return Status::Error(StatusCode::MissingInput,
          "turbulence_adapter moment source lacks moment, residual energy, or slab fraction");
    if (!Finite(*slabFraction) || *slabFraction < 0.0 ||
        *slabFraction > 1.0)
      return Status::Error(StatusCode::InvalidBackground,
                           "turbulence_adapter slab fraction must be in [0,1]");
    double total = 0.0;
    Status converted = ConvertTurbulenceMoment(*moment,
        *local.turbulence->densityKgPerM3, *residual,
        a.vacuumPermeabilityHPerM,
        a.energyConvention, a.residualConvention, &total);
    if (!converted.ok()) return converted;
    LocalState adapted = local;
    adapted.turbulence->slabVarianceT2 = *slabFraction * total;
    adapted.turbulence->twoDVarianceT2 = (1.0 - *slabFraction) * total;
    ModelConfiguration child = configuration;
    if (a.closure == AdapterClosure::QltSlabSpectrum) {
      child.model = ModelId::QltSlabSpectrum; child.qltSlab = a.qlt;
    } else if (a.closure == AdapterClosure::BroadenedSlab) {
      child.model = ModelId::BroadenedSlab; child.broadenedSlab = a.broadened;
    } else {
      return Status::Error(StatusCode::UnsupportedModel,
                           "inertial QLT has no D_mu_mu evaluator");
    }
    return EvaluateAdvancedPitchAngle(mu, particle, adapted, child, dMuMuPerS);
  }
  if (configuration.model == ModelId::WaveSpectrumAdapter) {
    const auto& a = configuration.waveSpectrumAdapter;
    ModelConfiguration child = configuration;
    if (a.closure == AdapterClosure::QltSlabSpectrum) {
      child.model = ModelId::QltSlabSpectrum; child.qltSlab = a.qlt;
    } else if (a.closure == AdapterClosure::BroadenedSlab) {
      child.model = ModelId::BroadenedSlab; child.broadenedSlab = a.broadened;
    } else {
      return Status::Error(StatusCode::UnsupportedModel,
                           "inertial QLT has no D_mu_mu evaluator");
    }
    return EvaluateAdvancedPitchAngle(mu, particle, local, child, dMuMuPerS);
  }
  if (configuration.model == ModelId::PrescribedLambdaMuShape) {
    if (!dMuMuPerS || !Finite(mu) || std::fabs(mu) > 1.0)
      return Status::Error(StatusCode::InvalidConfiguration,
                           "pitch angle mu must be finite and in [-1,1]");
    const auto& p = configuration.prescribedLambdaMu;
    bool ok = false;
    const double integral = ShapeIntegral(p.qMu, p.hMu, p.numerical, &ok);
    if (!ok) return Status::Error(StatusCode::IntegrationFailed,
                                  "shape normalization integral failed");
    const double d0 = p.amplitudeMode == PitchAngleAmplitudeMode::TargetLambda
        ? 3.0 * kin.speedMPerS * integral / (8.0 * p.targetLambdaM)
        : p.fixedD0PerS;
    *dMuMuPerS = d0 * (1.0 - mu * mu) *
        (std::pow(std::fabs(mu), p.qMu - 1.0) + p.hMu);
    return Finite(*dMuMuPerS) && *dMuMuPerS >= 0.0
        ? Status::Success()
        : Status::Error(StatusCode::OutsideModelDomain,
                        "pitch-angle coefficient is not representable");
  }
  return Status::Error(StatusCode::UnsupportedModel,
                       "selected model has no D_mu_mu interface");
}

Status BuildAdvancedConfiguration(const std::string& modelId,
                                  const std::vector<InputParameter>& parameters,
                                  ModelConfiguration* configuration) {
  ModelId model;
  if (!configuration || !ParseModelId(modelId, &model) || !IsAdvancedModel(model))
    return Status::Error(StatusCode::UnsupportedModel,
                         "unknown advanced parallel-diffusion model '" + modelId + "'");
  ModelConfiguration candidate;
  candidate.model = model;
  std::map<std::string, std::string> values;
  Status status = Parameters(parameters, &values);
  if (!status.ok()) return status;

  // Each branch below is the selected model's semantic reader.  Shared
  // sub-readers handle only genuinely common schemas (numerical controls and
  // canonical spectra); a key is consumed solely by the selected branch.
  // The remaining-map check therefore rejects parameters belonging to a
  // different model instead of letting a broad manager silently ignore them.
  if (model == ModelId::QltSlabSpectrum) {
    status = ParseSpectrum(&values, &candidate.qltSlab.spectrum);
    if (status.ok()) status = ParseNumerical(&values, &candidate.qltSlab.numerical);
  } else if (model == ModelId::QltSlabInertial || model == ModelId::NlpaGivenPerp ||
             model == ModelId::NlgcE || model == ModelId::NlgceN) {
    if (model != ModelId::QltSlabInertial)
      status = ParseNumerical(&values, &candidate.nonlinear.numerical);
  } else if (model == ModelId::PrescribedLambdaMuShape) {
    auto& p = candidate.prescribedLambdaMu;
    std::string mode;
    status = Take(&values, "amplitude_mode", true, &mode);
    if (!status.ok()) return status;
    mode = Lower(mode);
    if (mode == "target_lambda") p.amplitudeMode = PitchAngleAmplitudeMode::TargetLambda;
    else if (mode == "fixed_amplitude") p.amplitudeMode = PitchAngleAmplitudeMode::FixedAmplitude;
    else return Status::Error(StatusCode::InvalidConfiguration,
                              "amplitude_mode must be target_lambda or fixed_amplitude");
    status = TakeDouble(&values, "q_mu", true, &p.qMu);
    if (status.ok()) status = TakeDouble(&values, "h_mu", true, &p.hMu);
    if (status.ok()) status = TakeDouble(&values,
        p.amplitudeMode == PitchAngleAmplitudeMode::TargetLambda
            ? "target_lambda_m" : "D0_per_s", true,
        p.amplitudeMode == PitchAngleAmplitudeMode::TargetLambda
            ? &p.targetLambdaM : &p.fixedD0PerS);
    if (status.ok()) status = ParseNumerical(&values, &p.numerical);
  } else if (model == ModelId::BroadenedSlab) {
    auto& p = candidate.broadenedSlab;
    status = ParseBroadened(&values, &p);
  } else if (model == ModelId::NlgceF2014) {
    std::string set;
    status = Take(&values, "coefficient_set", false, &set);
    if (status.ok() && !set.empty()) candidate.nlgceF.coefficientSet = set;
  } else if (model == ModelId::TurbulenceAdapter) {
    auto& p = candidate.turbulenceAdapter;
    std::string energy, residual, closure, source;
    status = Take(&values, "energy_convention", true, &energy);
    if (status.ok()) status = Take(&values, "residual_energy_convention", true, &residual);
    if (status.ok()) status = Take(&values, "moment_source", true, &source);
    if (status.ok()) status = Take(&values, "underlying_closure", true, &closure);
    if (!status.ok()) return status;
    energy = Lower(energy); residual = Lower(residual);
    if (energy == "kinetic_plus_magnetic_variance") p.energyConvention = MomentEnergyConvention::KineticPlusMagneticVariance;
    else if (energy == "half_elsasser_sum") p.energyConvention = MomentEnergyConvention::HalfElsasserSum;
    else if (energy == "elsasser_sum") p.energyConvention = MomentEnergyConvention::ElsasserSum;
    else if (energy == "specific_total_fluctuation_energy") p.energyConvention = MomentEnergyConvention::SpecificTotalFluctuationEnergy;
    else return Status::Error(StatusCode::InvalidConfiguration, "unknown energy_convention");
    if (residual == "kinetic_minus_magnetic") p.residualConvention = ResidualEnergyConvention::KineticMinusMagnetic;
    else if (residual == "magnetic_minus_kinetic") p.residualConvention = ResidualEnergyConvention::MagneticMinusKinetic;
    else return Status::Error(StatusCode::InvalidConfiguration, "unknown residual_energy_convention");
    source = Lower(source);
    if (source == "input_parameters") p.momentSource = MomentSource::InputParameters;
    else if (source == "local_state") p.momentSource = MomentSource::LocalState;
    else return Status::Error(StatusCode::InvalidConfiguration,
                              "moment_source must be input_parameters or local_state");
    status = ParseClosure(closure, &p.closure);
    double value = 0.0;
    if (status.ok())
      status = TakeDouble(&values, "vacuum_permeability_H_per_m", true,
                          &p.vacuumPermeabilityHPerM);
    if (status.ok() && p.momentSource == MomentSource::InputParameters) {
      status = TakeDouble(&values, "provider_moment_m2_per_s2", true, &value);
      if (status.ok()) p.providerMomentM2PerS2 = value;
      if (status.ok()) status = TakeDouble(&values, "residual_energy", true, &value);
      if (status.ok()) p.residualEnergy = value;
      if (status.ok()) status = TakeDouble(&values, "slab_fraction", true, &value);
      if (status.ok()) p.slabFraction = value;
    }
    if (status.ok() && p.closure == AdapterClosure::QltSlabSpectrum)
      status = ParseSpectrum(&values, &p.qlt.spectrum);
    if (status.ok() && p.closure == AdapterClosure::QltSlabSpectrum)
      status = ParseNumerical(&values, &p.qlt.numerical);
    if (status.ok() && p.closure == AdapterClosure::BroadenedSlab)
      status = ParseBroadened(&values, &p.broadened);
  } else if (model == ModelId::WaveSpectrumAdapter) {
    auto& p = candidate.waveSpectrumAdapter;
    std::string closure;
    status = Take(&values, "propagation", true, &p.propagation);
    if (status.ok()) status = Take(&values, "polarization", true, &p.polarization);
    if (status.ok()) status = Take(&values, "frame", true, &p.frame);
    if (status.ok()) status = Take(&values, "underlying_closure", true, &closure);
    if (!status.ok()) return status;
    p.propagation = Lower(p.propagation); p.polarization = Lower(p.polarization);
    p.frame = Lower(p.frame);
    status = ParseClosure(closure, &p.closure);
    if (status.ok() && p.closure == AdapterClosure::QltSlabSpectrum)
      status = ParseSpectrum(&values, &p.qlt.spectrum);
    if (status.ok() && p.closure == AdapterClosure::QltSlabSpectrum)
      status = ParseNumerical(&values, &p.qlt.numerical);
    if (status.ok() && p.closure == AdapterClosure::BroadenedSlab)
      status = ParseBroadened(&values, &p.broadened);
  } else if (model == ModelId::TabulatedParallel) {
    auto& p = candidate.table;
    std::string quantity, axes, coefficients;
    status = Take(&values, "stored_quantity", true, &quantity);
    if (status.ok()) status = Take(&values, "axes", true, &axes);
    if (status.ok()) status = Take(&values, "coefficient_values_SI", true, &coefficients);
    if (status.ok()) status = Take(&values, "generation_identity", true, &p.generationIdentity);
    if (!status.ok()) return status;
    quantity = Lower(quantity);
    if (quantity == "lambda_parallel") p.storedCoefficient = StoredCoefficient::LambdaParallel;
    else if (quantity == "kappa_parallel") p.storedCoefficient = StoredCoefficient::KappaParallel;
    else return Status::Error(StatusCode::InvalidConfiguration,
                              "stored_quantity must be lambda_parallel or kappa_parallel");
    std::istringstream stream(axes);
    std::string name;
    std::size_t index = 0;
    while (std::getline(stream, name, ',')) {
      TableAxis axis;
      status = ParseTableAxis(name, &axis);
      if (!status.ok()) return status;
      p.axes.push_back(axis);
      std::string axisText;
      status = Take(&values, "axis_" + std::to_string(index) + "_values_SI", true, &axisText);
      if (!status.ok()) return status;
      p.axisSI.emplace_back();
      status = ParseList(axisText, &p.axisSI.back());
      if (!status.ok()) return status;
      ++index;
    }
    if (std::find(p.axes.begin(), p.axes.end(), TableAxis::Time) !=
        p.axes.end()) {
      std::string rule;
      status = Take(&values, "time_rule", true, &rule);
      if (!status.ok()) return status;
      rule = Lower(rule);
      if (rule == "linear")
        p.timeInterpolation = TimeInterpolation::Linear;
      else if (rule == "step_previous")
        p.timeInterpolation = TimeInterpolation::StepPrevious;
      else
        return Status::Error(StatusCode::InvalidConfiguration,
            "time_rule must be linear or step_previous");
    }
    status = ParseList(coefficients, &p.coefficientSI);
  }
  if (!status.ok()) return status;
  if (!values.empty())
    return Status::Error(StatusCode::InvalidConfiguration,
                         "unknown parameter '" + values.begin()->first + "'");
  status = ValidateAdvancedConfiguration(candidate);
  if (!status.ok()) return status;
  *configuration = std::move(candidate);
  return Status::Success();
}

}  // namespace Internal

namespace {

Status SpectrumPowerSum(const std::initializer_list<double>& terms,
                        double multiplier, double* canonicalPowerT2M) {
  if (!canonicalPowerT2M || !std::isfinite(multiplier) || multiplier <= 0.0)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "invalid canonical-spectrum conversion output");
  double sum = 0.0;
  for (double term : terms) {
    if (!std::isfinite(term) || term < 0.0)
      return Status::Error(StatusCode::InvalidConfiguration,
                           "spectrum power must be finite and nonnegative");
    sum += term;
  }
  const double converted = multiplier * sum;
  if (!std::isfinite(converted))
    return Status::Error(StatusCode::OutsideModelDomain,
                         "converted spectrum power is not representable");
  *canonicalPowerT2M = converted;
  return Status::Success();
}

}  // namespace

Status ConvertOneSidedComponentsToCanonical(
    double pXxT2M, double pYyT2M, double* canonicalPowerT2M) {
  // Section 8.7: each one-sided component integrates to its component
  // variance, so their sum integrates to the total transverse variance.
  return SpectrumPowerSum({pXxT2M, pYyT2M}, 1.0, canonicalPowerT2M);
}

Status ConvertTwoSidedComponentsToCanonical(
    double sXxPositiveT2M, double sXxNegativeT2M,
    double sYyPositiveT2M, double sYyNegativeT2M,
    double* canonicalPowerT2M) {
  // No evenness or transverse-axisymmetry is presumed: all four signed-k
  // component samples appear explicitly in the Section 8.7 sum.
  return SpectrumPowerSum({sXxPositiveT2M, sXxNegativeT2M,
                           sYyPositiveT2M, sYyNegativeT2M},
                          1.0, canonicalPowerT2M);
}

Status ConvertEvenTwoSidedTotalToCanonical(
    double signedTotalPowerT2M, double* canonicalPowerT2M) {
  // For an even signed total spectrum, folding negative k onto k>0 doubles
  // the density while leaving its integrated variance unchanged.
  return SpectrumPowerSum({signedTotalPowerT2M}, 2.0,
                          canonicalPowerT2M);
}

Status ConvertOneSidedCyclesPerMToCanonical(
    double wavenumberCyclesPerM, double powerT2PerCyclePerM,
    double* wavenumberRadPerM, double* canonicalPowerT2M) {
  if (!wavenumberRadPerM || !canonicalPowerT2M ||
      !std::isfinite(wavenumberCyclesPerM) || wavenumberCyclesPerM <= 0.0 ||
      !std::isfinite(powerT2PerCyclePerM) || powerT2PerCyclePerM < 0.0)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "invalid cycles-per-metre spectrum sample");
  const double twoPi = 2.0 * std::acos(-1.0);
  // k=2*pi*k_c and P_s(k) dk=P^(c)(k_c) dk_c preserve the variance.
  const double convertedK = twoPi * wavenumberCyclesPerM;
  const double convertedPower = powerT2PerCyclePerM / twoPi;
  if (!std::isfinite(convertedK) || !std::isfinite(convertedPower))
    return Status::Error(StatusCode::OutsideModelDomain,
                         "cycles-per-metre conversion is not representable");
  *wavenumberRadPerM = convertedK;
  *canonicalPowerT2M = convertedPower;
  return Status::Success();
}

Status ConvertFrozenFlowFrequencyToCanonical(
    double frequencyHz, double oneSidedFrequencyPsdT2PerHz,
    const FrozenFlowMapping& mapping, double* wavenumberRadPerM,
    double* canonicalPowerT2M) {
  if (!wavenumberRadPerM || !canonicalPowerT2M ||
      !std::isfinite(frequencyHz) || frequencyHz <= 0.0 ||
      !std::isfinite(oneSidedFrequencyPsdT2PerHz) ||
      oneSidedFrequencyPsdT2PerHz < 0.0 ||
      !std::isfinite(mapping.samplingVelocityProjectionMPerS) ||
      mapping.samplingVelocityProjectionMPerS <= 0.0 ||
      mapping.assumptionIdentity.empty())
    return Status::Error(StatusCode::InvalidConfiguration,
        "frozen-flow conversion requires positive frequency and sampling "
        "projection, nonnegative PSD, and an assumption identity");
  const double twoPi = 2.0 * std::acos(-1.0);
  // Equation (36).  The Jacobian U_T/(2*pi) ensures P_s(k)dk=S_f(f)df.
  const double convertedK =
      twoPi * frequencyHz / mapping.samplingVelocityProjectionMPerS;
  const double convertedPower = mapping.samplingVelocityProjectionMPerS *
      oneSidedFrequencyPsdT2PerHz / twoPi;
  if (!std::isfinite(convertedK) || !std::isfinite(convertedPower))
    return Status::Error(StatusCode::OutsideModelDomain,
                         "frozen-flow conversion is not representable");
  *wavenumberRadPerM = convertedK;
  *canonicalPowerT2M = convertedPower;
  return Status::Success();
}

Status ConvertQinZhangSlabComponentToCanonical(
    double sourcePowerT2M, double* canonicalPowerT2M) {
  // Qin and Zhang's signed, single-component S'_xx^slab integrates to half
  // the total slab variance; folding signed k and summing components gives 4.
  return SpectrumPowerSum({sourcePowerT2M}, 4.0, canonicalPowerT2M);
}

Status ConvertQinZhangTwoDComponentToCanonical(
    double sourcePowerT2M, double* canonicalPowerT2M) {
  // The same factor applies to their reduced radial 2D density under the
  // delta-function convention stated in specification Section 8.7.
  return SpectrumPowerSum({sourcePowerT2M}, 4.0, canonicalPowerT2M);
}

Status ConvertTurbulenceMoment(double providerMomentM2PerS2,
                               double densityKgPerM3,
                               double residualEnergy,
                               double vacuumPermeabilityHPerM,
                               MomentEnergyConvention energyConvention,
                               ResidualEnergyConvention residualConvention,
                               double* totalMagneticVarianceT2) {
  if (!totalMagneticVarianceT2 || !std::isfinite(providerMomentM2PerS2) ||
      providerMomentM2PerS2 < 0.0 || !std::isfinite(densityKgPerM3) ||
      densityKgPerM3 <= 0.0 || !std::isfinite(residualEnergy) ||
      residualEnergy < -1.0 || residualEnergy > 1.0 ||
      !std::isfinite(vacuumPermeabilityHPerM) ||
      vacuumPermeabilityHPerM <= 0.0)
    return Status::Error(StatusCode::InvalidBackground,
                         "invalid turbulence moment, density, or residual energy");
  double z2 = providerMomentM2PerS2;
  if (energyConvention == MomentEnergyConvention::ElsasserSum) z2 *= 0.5;
  else if (energyConvention == MomentEnergyConvention::SpecificTotalFluctuationEnergy)
    z2 *= 2.0;
  double sigma = residualEnergy;
  if (residualConvention == ResidualEnergyConvention::MagneticMinusKinetic)
    sigma = -sigma;
  // Equation (57b) with the named convention already converted to canonical
  // Z^2.  Mu0 is SI, yielding T^2 from kg m^-3 times m^2 s^-2.
  *totalMagneticVarianceT2 = vacuumPermeabilityHPerM * densityKgPerM3 *
      (1.0 - sigma) * z2 / 2.0;
  return std::isfinite(*totalMagneticVarianceT2) &&
         *totalMagneticVarianceT2 >= 0.0
      ? Status::Success()
      : Status::Error(StatusCode::OutsideModelDomain,
                      "converted magnetic variance is not representable");
}

}  // namespace ParallelDiffusion
}  // namespace SEP
