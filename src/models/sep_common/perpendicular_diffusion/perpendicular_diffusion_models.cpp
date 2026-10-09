#include "perpendicular_diffusion_internal.h"

#include "../parallel_diffusion/parallel_diffusion.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <utility>

namespace SEP {
namespace PerpendicularDiffusion {
namespace {

constexpr double Pi = 3.141592653589793238462643383279502884;
constexpr double SqrtPi = 1.77245385090551602729816748334114518;

bool Finite(double x) { return std::isfinite(x); }
bool Positive(double x) { return Finite(x) && x > 0.0; }
bool Nonnegative(double x) { return Finite(x) && x >= 0.0; }

double Parameter(const ModelConfiguration& c, const char* name) {
  const auto found = c.numeric.find(name);
  return found == c.numeric.end() ? std::numeric_limits<double>::quiet_NaN()
                                 : found->second;
}

double OptionalParameter(const ModelConfiguration& c, const char* name,
                         double fallback) {
  const auto found = c.numeric.find(name);
  return found == c.numeric.end() ? fallback : found->second;
}

std::string TextParameter(const ModelConfiguration& c, const char* name) {
  const auto found = c.text.find(name);
  return found == c.text.end() ? std::string() : found->second;
}

double VectorNorm(const std::array<double, 3>& x) {
  return std::hypot(std::hypot(x[0], x[1]), x[2]);
}

double Dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> Cross(const std::array<double, 3>& a,
                            const std::array<double, 3>& b) {
  return {{a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}};
}

Status Normalize(const std::array<double, 3>& input,
                 std::array<double, 3>* output) {
  const double norm = VectorNorm(input);
  if (!output || !Positive(norm))
    return Status::Error(StatusCode::InvalidInput,
                         "a finite nonzero vector is required");
  for (int i = 0; i < 3; ++i) (*output)[i] = input[i] / norm;
  return Status::Success();
}

ModelResult BaseResult(const LocalState& local,
                       const ModelConfiguration& configuration) {
  ModelResult result;
  result.status = Status::Success();
  result.provenance.requestedModel = ModelName(configuration.model);
  result.provenance.usedModel = ModelName(configuration.model);
  result.provenance.equationVersion = "perpendicular-spec-2.1";
  result.provenance.sourceVersion = "PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md:2.1";
  result.provenance.configurationFingerprint =
      ConfigurationFingerprint(configuration);
  result.provenance.sampleFingerprint = local.turbulence.sampleFingerprint;
  return result;
}

ModelResult Failure(const LocalState& local,
                    const ModelConfiguration& configuration,
                    StatusCode code, const std::string& detail) {
  ModelResult result = BaseResult(local, configuration);
  result.status = Status::Error(code, detail);
  return result;
}

Status MeanField(const LocalState& local, double* magnitude,
                 std::array<double, 3>* direction) {
  if (!local.meanFieldT.has_value())
    return Status::Error(StatusCode::MissingInput,
                         "selected model requires mean_field_T");
  const double value = VectorNorm(*local.meanFieldT);
  if (!Positive(value))
    return Status::Error(StatusCode::OutsideModelDomain,
                         "selected ordered-field model requires B0>0");
  if (magnitude) *magnitude = value;
  if (direction) {
    for (int i = 0; i < 3; ++i) (*direction)[i] = (*local.meanFieldT)[i] / value;
  }
  return Status::Success();
}

Status Parallel(const LocalState& local, double* kappa) {
  if (!local.parallelDependency.has_value())
    return Status::Error(StatusCode::MissingInput,
                         "selected model requires a supplied parallel coefficient");
  const ParallelInput& input = *local.parallelDependency;
  if (!Nonnegative(input.kappaM2PerS) || input.modelId.empty() ||
      input.equationVersion.empty() || input.sampleFingerprint.empty())
    return Status::Error(StatusCode::InvalidInput,
        "parallel dependency requires a nonnegative coefficient and complete provenance");
  if (!local.turbulence.sampleFingerprint.empty() &&
      input.sampleFingerprint != local.turbulence.sampleFingerprint)
    return Status::Error(StatusCode::InconsistentPair,
        "parallel and perpendicular inputs do not identify the same sample");
  *kappa = input.kappaM2PerS;
  return Status::Success();
}

Status Kinematics(const ParticleState& p, ParticleKinematics* k) {
  return ComputeParticleKinematics(p, k);
}

Status OrderedFrame(const LocalState& local, bool unequal,
                    PhysicalFrame* frame) {
  if (!frame) return Status::Error(StatusCode::InvalidInput, "null frame output");
  std::array<double, 3> b;
  Status status = MeanField(local, nullptr, &b);
  if (!status.ok()) return status;
  std::array<double, 3> e1;
  if (unequal) {
    if (!local.perpendicularAxis1.has_value() ||
        !local.perpendicularAxis2.has_value())
      return Status::Error(StatusCode::MissingInput,
          "unequal perpendicular eigenvalues require two supplied physical axes");
    status = Normalize(*local.perpendicularAxis1, &e1);
    if (!status.ok()) return status;
    std::array<double, 3> e2;
    status = Normalize(*local.perpendicularAxis2, &e2);
    if (!status.ok()) return status;
    if (std::fabs(Dot(b, e1)) > 1.0e-10 ||
        std::fabs(Dot(b, e2)) > 1.0e-10 ||
        std::fabs(Dot(e1, e2)) > 1.0e-10 ||
        Dot(Cross(b, e1), e2) < 1.0 - 1.0e-10)
      return Status::Error(StatusCode::InvalidInput,
          "supplied perpendicular axes must make a right-handed orthonormal frame");
    frame->kind = FrameKind::OrientedUnequal;
    frame->b = b; frame->e1 = e1; frame->e2 = e2;
    frame->perpendicularAxisDefinition = "caller_supplied_physical_axes";
    return Status::Success();
  }

  // Equal transverse eigenvalues are invariant under rotations about b.  A
  // deterministic Cartesian reference is therefore only a tensor factorization,
  // not an added physical orientation.  Select the least-aligned coordinate
  // direction to avoid loss of significance in the cross product.
  std::array<double, 3> reference{{1.0, 0.0, 0.0}};
  if (std::fabs(b[1]) <= std::fabs(b[0]) && std::fabs(b[1]) <= std::fabs(b[2]))
    reference = {{0.0, 1.0, 0.0}};
  else if (std::fabs(b[2]) <= std::fabs(b[0]))
    reference = {{0.0, 0.0, 1.0}};
  e1 = Cross(reference, b);
  status = Normalize(e1, &e1);
  if (!status.ok()) return status;
  const std::array<double, 3> e2 = Cross(b, e1);
  frame->kind = FrameKind::AxisymmetricOrdered;
  frame->b = b; frame->e1 = e1; frame->e2 = e2;
  frame->perpendicularAxisDefinition =
      "degenerate_equal_eigenvalue_cartesian_factorization";
  return Status::Success();
}

ModelResult EqualCoefficient(const LocalState& local,
                             const ModelConfiguration& configuration,
                             double kappa, Quality quality,
                             DomainState domain = DomainState::NotApplicable) {
  if (!Nonnegative(kappa))
    return Failure(local, configuration, StatusCode::InvalidInput,
                   "model produced a nonfinite or negative coefficient");
  ModelResult result = BaseResult(local, configuration);
  result.observable = Observable::SymmetricCoefficient;
  result.estimator = Estimator::Asymptotic;
  result.quality = quality;
  result.domain = domain;
  result.diffusionRegime = DiffusionRegime::NormalDiffusion;
  SpatialCoefficients coefficients;
  Status frame = OrderedFrame(local, false, &coefficients.frame);
  if (!frame.ok()) return Failure(local, configuration, frame.code, frame.detail);
  coefficients.perpendicular1M2PerS = kappa;
  coefficients.perpendicular2M2PerS = kappa;
  if (local.parallelDependency.has_value()) {
    double parallel = 0.0;
    const Status supplied = Parallel(local, &parallel);
    if (!supplied.ok())
      return Failure(local, configuration, supplied.code, supplied.detail);
    coefficients.parallelM2PerS = parallel;
    coefficients.parallelOwner = DependencyOwner::SuppliedParallel;
    result.provenance.dependencyOwner = DependencyOwner::SuppliedParallel;
    result.provenance.parallelModelId = local.parallelDependency->modelId;
  }
  result.coefficients = coefficients;
  return result;
}

ModelResult UnequalCoefficient(const LocalState& local,
                               const ModelConfiguration& configuration,
                               double parallel, double first, double second,
                               DependencyOwner owner, Quality quality,
                               DomainState domain) {
  if (!Nonnegative(parallel) || !Nonnegative(first) || !Nonnegative(second))
    return Failure(local, configuration, StatusCode::InvalidInput,
                   "model produced a nonfinite or negative coefficient");
  ModelResult result = BaseResult(local, configuration);
  result.observable = Observable::PairedCoefficients;
  result.estimator = Estimator::Asymptotic;
  result.quality = quality;
  result.domain = domain;
  result.diffusionRegime = DiffusionRegime::NormalDiffusion;
  SpatialCoefficients coefficients;
  const Status frame = OrderedFrame(local, first != second, &coefficients.frame);
  if (!frame.ok()) return Failure(local, configuration, frame.code, frame.detail);
  coefficients.parallelM2PerS = parallel;
  coefficients.perpendicular1M2PerS = first;
  coefficients.perpendicular2M2PerS = second;
  coefficients.parallelOwner = owner;
  result.coefficients = coefficients;
  result.provenance.dependencyOwner = owner;
  if (owner == DependencyOwner::SuppliedParallel && local.parallelDependency)
    result.provenance.parallelModelId = local.parallelDependency->modelId;
  return result;
}

ModelResult IsotropicCoefficient(const LocalState& local,
                                 const ModelConfiguration& configuration,
                                 double kappa, Quality quality,
                                 DomainState domain) {
  if (!Nonnegative(kappa))
    return Failure(local, configuration, StatusCode::InvalidInput,
                   "model produced a nonfinite or negative isotropic coefficient");
  ModelResult result = BaseResult(local, configuration);
  result.observable = Observable::SymmetricCoefficient;
  result.estimator = Estimator::Asymptotic;
  result.quality = quality;
  result.domain = domain;
  result.diffusionRegime = DiffusionRegime::NormalDiffusion;
  SpatialCoefficients coefficients;
  coefficients.frame.kind = FrameKind::Isotropic;
  coefficients.frame.perpendicularAxisDefinition = "statistical_isotropy";
  coefficients.isotropicM2PerS = kappa;
  result.coefficients = coefficients;
  return result;
}

struct CompositeState {
  double b0 = 0.0;
  double slabVariance = 0.0;
  double twoDVariance = 0.0;
  double slabLength = 0.0;
  double twoDLength = 0.0;
  double s = 0.0;
  double q = 0.0;
};

Status Composite(const LocalState& local, CompositeState* state,
                 bool needSlab, bool needTwoD, bool requireQGreaterOne = false) {
  if (!state) return Status::Error(StatusCode::InvalidInput, "null turbulence output");
  Status status = MeanField(local, &state->b0, nullptr);
  if (!status.ok()) return status;
  const TurbulenceSample& t = local.turbulence;
  if (needSlab) {
    if (!t.slabVarianceT2)
      return Status::Error(StatusCode::MissingInput,
                           "selected model requires slab variance");
    state->slabVariance = *t.slabVarianceT2;
    if (!Nonnegative(state->slabVariance))
      return Status::Error(StatusCode::InvalidInput,
                           "slab variance must be nonnegative");
    if (state->slabVariance>0.0) {
      if(!t.slabBendoverLengthM)
        return Status::Error(StatusCode::MissingInput,
                             "nonzero slab power requires its bend-over length");
      state->slabLength=*t.slabBendoverLengthM;
      if(!Positive(state->slabLength))
        return Status::Error(StatusCode::InvalidInput,
                             "slab bend-over length must be positive");
    }
  }
  if (needTwoD) {
    if (!t.twoDVarianceT2)
      return Status::Error(StatusCode::MissingInput,
                           "selected model requires 2D variance");
    state->twoDVariance = *t.twoDVarianceT2;
    if (!Nonnegative(state->twoDVariance))
      return Status::Error(StatusCode::InvalidInput,
                           "2D variance must be nonnegative");
    if(state->twoDVariance>0.0) {
      if(!t.twoDBendoverLengthM)
        return Status::Error(StatusCode::MissingInput,
                             "nonzero 2D power requires its bend-over length");
      state->twoDLength=*t.twoDBendoverLengthM;
      if(!Positive(state->twoDLength))
        return Status::Error(StatusCode::InvalidInput,
                             "2D bend-over length must be positive");
    }
  }
  if((state->slabVariance>0.0 || state->twoDVariance>0.0)) {
    if(!t.inertialIndex)
      return Status::Error(StatusCode::MissingInput,
                           "nonzero smooth-spectrum power requires inertial_index");
    state->s=*t.inertialIndex;
    if(!(state->s>1.0))
      return Status::Error(StatusCode::OutsideModelDomain,
                           "smooth spectrum requires s>1");
  }
  if(state->twoDVariance>0.0) {
    if(!t.energyRangeIndex)
      return Status::Error(StatusCode::MissingInput,
                           "nonzero 2D power requires energy_index");
    state->q=*t.energyRangeIndex;
    if(!(state->q>-1.0))
      return Status::Error(StatusCode::OutsideModelDomain,
                           "smooth 2D spectrum requires q>-1");
    if(requireQGreaterOne && !(state->q>1.0))
      return Status::Error(StatusCode::DivergentMoment,
                           "selected closure requires finite I_-2, hence q>1");
  }
  return Status::Success();
}

double SlabFieldLine(const CompositeState& s) {
  return Pi * SpectrumC(s.s) * s.slabLength * s.slabVariance / (s.b0 * s.b0);
}

double TwoDFieldLine(const CompositeState& s) {
  return std::sqrt((s.s - 1.0) / (2.0 * (s.q - 1.0)) *
      s.twoDLength * s.twoDLength * s.twoDVariance / (s.b0 * s.b0));
}

double CompositeFieldLine(const CompositeState& s) {
  const double slab = SlabFieldLine(s);
  const double twoD = TwoDFieldLine(s);
  return 0.5 * (slab + std::hypot(slab, 2.0 * twoD));
}

Status SelectedFieldLine(const LocalState& local, double* value) {
  if (!local.fieldLineCoefficientM.has_value())
    return Status::Error(StatusCode::MissingInput,
                         "selected model requires a field-line coefficient");
  if (!Nonnegative(*local.fieldLineCoefficientM) ||
      local.fieldLineModelId.empty() || local.fieldLineSampleFingerprint.empty())
    return Status::Error(StatusCode::InvalidInput,
        "field-line dependency requires nonnegative value and complete identity");
  if (!local.turbulence.sampleFingerprint.empty() &&
      local.fieldLineSampleFingerprint != local.turbulence.sampleFingerprint)
    return Status::Error(StatusCode::InconsistentPair,
                         "field-line and turbulence samples are incoherent");
  *value = *local.fieldLineCoefficientM;
  return Status::Success();
}

Status IsotropicInputs(const ParticleState& particle, const LocalState& local,
                       ParticleKinematics* kin, double* b0, double* variance,
                       double* totalField);

// Adaptive Simpson integration is applied only on the dimensionless log-k
// coordinate.  The finite log interval is expanded until successive results
// agree; it is an algorithmic tail control and is never called a dissipation
// cutoff.  This keeps integrals positive while spanning many length scales.
double Simpson(double a, double b, double fa, double fm, double fb) {
  return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

bool AdaptiveSimpson(const std::function<double(double)>& f, double a, double b,
                     double fa, double fm, double fb, double whole,
                     double tolerance, int depth, double* value,
                     double* error) {
  const double middle = 0.5 * (a + b);
  const double leftMiddle = 0.5 * (a + middle);
  const double rightMiddle = 0.5 * (middle + b);
  const double fl = f(leftMiddle), fr = f(rightMiddle);
  if (!Finite(fl) || !Finite(fr)) return false;
  const double left = Simpson(a, middle, fa, fl, fm);
  const double right = Simpson(middle, b, fm, fr, fb);
  const double delta = left + right - whole;
  if (depth <= 0 || std::fabs(delta) <= 15.0 * tolerance) {
    *value = left + right + delta / 15.0;
    *error = std::fabs(delta / 15.0);
    return depth > 0 || std::fabs(delta) <= 150.0 * tolerance;
  }
  double lv = 0.0, le = 0.0, rv = 0.0, re = 0.0;
  if (!AdaptiveSimpson(f, a, middle, fa, fl, fm, left,
                       tolerance / 2.0, depth - 1, &lv, &le) ||
      !AdaptiveSimpson(f, middle, b, fm, fr, fb, right,
                       tolerance / 2.0, depth - 1, &rv, &re)) return false;
  *value = lv + rv; *error = le + re;
  return true;
}

Status IntegrateLog(const std::function<double(double)>& f,
                    const ModelConfiguration& configuration,
                    double* value, double* relativeError = nullptr) {
  const double rel = OptionalParameter(configuration, "relative_tolerance", 1.0e-9);
  const int refinements = static_cast<int>(OptionalParameter(
      configuration, "maximum_refinements", 18.0));
  double previous = -1.0;
  for (double bound : {24.0, 32.0, 40.0, 48.0}) {
    // Splitting the wide transformed interval prevents a narrow O(1) peak
    // around y=0 from being judged solely from two remote tail endpoints.
    // The split locations remain numerical controls, never physical cutoffs.
    double current = 0.0, error = 0.0;
    constexpr int Pieces = 16;
    for (int piece=0; piece<Pieces; ++piece) {
      const double a=-bound+2.0*bound*piece/Pieces;
      const double b=-bound+2.0*bound*(piece+1)/Pieces;
      const double middle=0.5*(a+b);
      const double fa=f(a),fm=f(middle),fb=f(b);
      if(!Finite(fa)||!Finite(fm)||!Finite(fb))
        return Status::Error(StatusCode::IntegrationFailed,
                             "nonfinite transformed spectral integrand");
      const double whole=Simpson(a,b,fa,fm,fb);
      double part=0.0,partError=0.0;
      if(!AdaptiveSimpson(f,a,b,fa,fm,fb,whole,rel/Pieces,
                          refinements,&part,&partError))
        return Status::Error(StatusCode::IntegrationFailed,
                             "adaptive spectral quadrature did not converge");
      current+=part; error+=partError;
    }
    if (current < 0.0 || !Finite(current))
      return Status::Error(StatusCode::IntegrationFailed,
                           "spectral quadrature produced an invalid integral");
    if (previous >= 0.0 && std::fabs(current - previous) <=
        rel * std::max({1.0, std::fabs(current), std::fabs(previous)})) {
      *value = current;
      if (relativeError)
        *relativeError = std::max(error, std::fabs(current - previous)) /
                         std::max(current, std::numeric_limits<double>::min());
      return Status::Success();
    }
    previous = current;
  }
  return Status::Error(StatusCode::IntegrationFailed,
                       "log-wavenumber tail expansion did not converge");
}

double LogOnePlusSquare(double y) {
  const double a = 2.0 * y;
  return a > 0.0 ? a + std::log1p(std::exp(-a)) : std::log1p(std::exp(a));
}

double SmoothShapeLog(double y, double s, double q) {
  return q * y - 0.5 * (s + q) * LogOnePlusSquare(y);
}

Status PositiveRoot(const std::function<double(double)>& logResidual,
                    double upper, const ModelConfiguration& configuration,
                    double* root, std::size_t* iterations = nullptr,
                    double* residual = nullptr) {
  if (!Positive(upper))
    return Status::Error(StatusCode::NonlinearSolverFailed,
                         "implicit closure has no positive upper bound");
  const int maximum = static_cast<int>(OptionalParameter(
      configuration, "maximum_iterations", 100.0));
  const double tolerance = OptionalParameter(configuration,
      "relative_tolerance", 1.0e-9);
  double high = std::log(upper);
  double low = high - 80.0;
  double fl = logResidual(low), fh = logResidual(high);
  if (!Finite(fl) || !Finite(fh))
    return Status::Error(StatusCode::NonlinearSolverFailed,
                         "implicit residual is nonfinite at its physical bracket");
  // A root exactly on the analytic upper bound is valid. Otherwise the
  // positive fixed-point kernels have a positive residual at sufficiently
  // small kappa and a nonpositive residual at their bound.
  if (std::fabs(fh) <= tolerance) {
    *root = upper;
    if (iterations) *iterations = 0;
    if (residual) *residual = std::fabs(fh);
    return Status::Success();
  }
  if (!(fl > 0.0 && fh < 0.0))
    return Status::Error(StatusCode::NonlinearSolverFailed,
                         "implicit residual did not bracket one positive root");
  double middle = 0.0, fm = 0.0;
  int i = 0;
  for (; i < maximum; ++i) {
    middle = 0.5 * (low + high);
    fm = logResidual(middle);
    if (!Finite(fm)) break;
    if (std::fabs(fm) <= tolerance || high - low <= tolerance) {
      *root = std::exp(middle);
      if (iterations) *iterations = static_cast<std::size_t>(i + 1);
      if (residual) *residual = std::fabs(fm);
      return Status::Success();
    }
    if (fm > 0.0) low = middle; else high = middle;
  }
  return Status::Error(StatusCode::NonlinearSolverFailed,
                       "positive logarithmic bisection did not converge");
}

DomainState DomainPolicy(bool inside, const ModelConfiguration& c,
                         Status* status, const std::string& description) {
  if (inside) return DomainState::InsideDeclaredDomain;
  if (TextParameter(c, "domain_policy") == "tag_extrapolated")
    return DomainState::Extrapolated;
  *status = Status::Error(StatusCode::OutsideModelDomain, description);
  return DomainState::NotApplicable;
}

struct FitRow { double g, nParallel, rhoParallel, nPerp, aPerp, nA; };

bool CandiaRow(const std::string& name, FitRow* row) {
  if (name == "kraichnan") *row = {1.5, 2.0, 0.22, 0.019, 1.37, 17.6};
  else if (name == "kolmogorov") *row = {5.0/3.0, 1.7, 0.20, 0.025, 1.36, 14.9};
  else if (name == "bykov_toptygin") *row = {2.0, 1.4, 0.16, 0.020, 1.38, 14.2};
  else return false;
  return true;
}

double Gc(double rigidity, double breakRigidity, double low, double high,
          double smoothing) {
  const double x = rigidity / breakRigidity;
  return std::pow(x, low) *
      std::pow(1.0 + std::pow(x, smoothing), (high - low) / smoothing);
}

double InterpolateLinear(const std::vector<double>& x,
                         const std::vector<double>& y, double value) {
  const auto high = std::upper_bound(x.begin(), x.end(), value);
  const std::size_t i = static_cast<std::size_t>(high - x.begin() - 1);
  const double fraction = (value - x[i]) / (x[i + 1] - x[i]);
  return y[i] + fraction * (y[i + 1] - y[i]);
}

}  // namespace

ModelResult EvaluateGcrPrescription(const ParticleState& particle,
                                    const LocalState& local,
                                    const ModelConfiguration& c) {
  // G3--G7 describe eigenvalues in a physical radial/polar transverse frame,
  // not spherical coordinate-tensor entries. The caller therefore supplies
  // the axes and colatitude in radians. NWU/HelMod scale one coherent external
  // parallel coefficient; Corti constructs its complete pair from G2/G5/G6.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double kp=0.0,b0=0.0;
  if(c.model==ModelId::NwuRatioPolar || c.model==ModelId::HelmodRatio) {
    status=Parallel(local,&kp);
    if(!status.ok()) return Failure(local,c,status.code,status.detail);
  }
  if(!local.heliocentricColatitudeRad)
    return Failure(local,c,StatusCode::MissingInput,
                   "latitude-dependent model requires colatitude_rad");
  const double theta=*local.heliocentricColatitudeRad;
  if(!Finite(theta)||theta<0.0||theta>Pi)
    return Failure(local,c,StatusCode::InvalidInput,"colatitude must lie in [0,pi]");
  if(c.model==ModelId::NwuRatioPolar) {
    if(TextParameter(c,"angular_convention")!="folded_radian")
      return Failure(local,c,StatusCode::InvalidConfiguration,
                     "NWU parameterized family requires angular_convention=folded_radian");
    const double thetaA=std::min(theta,Pi-theta);
    const double d=Parameter(c,"polar_enhancement");
    const double f=0.5*(d+1.0)-0.5*(d-1.0)*std::tanh(
        Parameter(c,"width_per_rad")*(thetaA-Pi/2.0+Parameter(c,"theta_F_rad")));
    return UnequalCoefficient(local,c,kp,Parameter(c,"eta_r")*kp,
        Parameter(c,"eta_theta")*kp*f,DependencyOwner::SuppliedParallel,
        Quality::ParameterizedClosure,DomainState::NotApplicable);
  }
  if(c.model==ModelId::HelmodRatio) {
    if(theta<c.tableAxisSI.front()||theta>c.tableAxisSI.back())
      return Failure(local,c,StatusCode::OutsideModelDomain,
                     "HelMod polar table does not cover this colatitude");
    const double factor=theta==c.tableAxisSI.back()?c.tableValuesSI.back():
        InterpolateLinear(c.tableAxisSI,c.tableValuesSI,theta);
    ModelResult result=UnequalCoefficient(local,c,kp,Parameter(c,"rho_H")*kp,
        Parameter(c,"rho_H")*kp*factor,DependencyOwner::SuppliedParallel,
        Quality::ParameterizedClosure,DomainState::InsideDeclaredDomain);
    result.provenance.conventions["polar_function_identity"]=
        TextParameter(c,"polar_function_identity");
    result.provenance.conventions["polar_interpolation"]=
        TextParameter(c,"polar_interpolation");
    return result;
  }
  status=MeanField(local,&b0,nullptr);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(TextParameter(c,"angular_convention")!="parameterized_folded_radian")
    return Failure(local,c,StatusCode::InvalidConfiguration,
        "Corti parameterized family requires angular_convention=parameterized_folded_radian");
  const double base=Parameter(c,"K0_m2_per_s")*kin.beta*
                    Parameter(c,"field_reference_T")/b0;
  const double parallel=base*Gc(kin.rigidityV,Parameter(c,"break_rigidity_V"),
      Parameter(c,"parallel_low_slope"),Parameter(c,"parallel_high_slope"),
      Parameter(c,"smoothness"));
  const double common=base*Gc(kin.rigidityV,Parameter(c,"break_rigidity_V"),
      Parameter(c,"perp_low_slope"),Parameter(c,"perp_high_slope"),
      Parameter(c,"smoothness"));
  const double u=1.5+0.5*std::tanh(Parameter(c,"width_per_rad")*
      (std::fabs(theta-Pi/2.0)-35.0*Pi/180.0));
  ModelResult result=UnequalCoefficient(local,c,parallel,0.02*common,0.01*u*common,
      DependencyOwner::PairedBackend,Quality::ParameterizedClosure,
      DomainState::NotApplicable);
  result.provenance.calibrationId=TextParameter(c,"calibration_id");
  return result;
}

ModelResult EvaluateTabulated(const ParticleState& particle,
                              const LocalState& local,
                              const ModelConfiguration& c) {
  // This bounded first backend is deliberately one-dimensional and linear.
  // It refuses extrapolation because no application validation domain or
  // derivative error budget has been supplied. Table checksums and generation
  // identities are provenance, not evidence that the values are physical.
  if(TextParameter(c,"boundary_policy")!="reject")
    return Failure(local,c,StatusCode::InvalidConfiguration,
                   "this release supports only boundary_policy=reject");
  const std::string axis=TextParameter(c,"axis");
  double coordinate=0.0;
  if(axis=="rigidity") {
    ParticleKinematics kin; Status status=Kinematics(particle,&kin);
    if(!status.ok()) return Failure(local,c,status.code,status.detail);
    coordinate=kin.rigidityV;
  } else if(axis=="time") coordinate=local.timeS;
  else if(axis=="heliocentric_radius") coordinate=VectorNorm(local.positionM);
  else if(axis=="mean_field_magnitude") {
    Status status=MeanField(local,&coordinate,nullptr);
    if(!status.ok()) return Failure(local,c,status.code,status.detail);
  } else return Failure(local,c,StatusCode::InvalidConfiguration,"unknown table axis");
  if(coordinate<c.tableAxisSI.front()||coordinate>c.tableAxisSI.back())
    return Failure(local,c,StatusCode::OutsideModelDomain,
                   "tabulated coefficient axis lies outside supplied coverage");
  double value=c.tableValuesSI.back();
  if(coordinate!=c.tableAxisSI.back()) {
    if(TextParameter(c,"interpolation")=="linear")
      value=InterpolateLinear(c.tableAxisSI,c.tableValuesSI,coordinate);
    else {
      std::vector<double> logAxis(c.tableAxisSI.size());
      std::vector<double> logValue(c.tableValuesSI.size());
      for(std::size_t i=0;i<logAxis.size();++i) {
        logAxis[i]=std::log(c.tableAxisSI[i]);
        logValue[i]=std::log(c.tableValuesSI[i]);
      }
      value=std::exp(InterpolateLinear(logAxis,logValue,std::log(coordinate)));
    }
  }
  ModelResult result=EqualCoefficient(local,c,value,Quality::ParameterizedClosure,
                                      DomainState::InsideDeclaredDomain);
  result.provenance.tableChecksum=TextParameter(c,"table_checksum");
  result.provenance.conventions["generation_identity"]=
      TextParameter(c,"generation_identity");
  result.provenance.conventions["interpolation"]=TextParameter(c,"interpolation");
  result.provenance.conventions["zero_policy"]=TextParameter(c,"zero_policy");
  return result;
}

ModelResult EvaluateScaling(const ParticleState& particle,
                            const LocalState& local,
                            const ModelConfiguration& c) {
  // I16 is retained as an explicitly normalized order-of-magnitude scaling;
  // the classical large-x branch is a diagnostic consequence of I14. Neither
  // is used as a fallback or ceiling for turbulence closures.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double kp=0.0,b0=0.0;
  status=Parallel(local,&kp); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(c.model==ModelId::YanLazarianMa4) {
    const double lambda=3.0*kp/kin.speedMPerS;
    if(!(lambda<Parameter(c,"injection_scale_m")))
      return Failure(local,c,StatusCode::OutsideModelDomain,
                     "M_A^4 branch requires lambda_parallel<L");
    return EqualCoefficient(local,c,Parameter(c,"c_M4")*kp*
        std::pow(Parameter(c,"alfven_mach"),4.0),Quality::ScalingEstimate);
  }
  status=MeanField(local,&b0,nullptr);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double radius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*b0);
  return EqualCoefficient(local,c,kin.speedMPerS*radius*radius/
      (3.0*(3.0*kp/kin.speedMPerS)),Quality::ScalingEstimate);
}

ModelResult EvaluateKuhlen(const ParticleState& particle,
                           const LocalState& local,
                           const ModelConfiguration& c) {
  // I8--I13 use the total-rms-field gyroradius, source correlation length,
  // arc-length field-line fit, and a calibrated decorrelation condition. The
  // published table does not determine gamma_K or L_c,perp, so both and the
  // root rule enter through the strict model schema. Returned running moments
  // retain the continuous post-freeze MSD rather than restarting at tau_c.
  ParticleKinematics kin; double b0=0.0,variance=0.0,brms=0.0;
  Status status=IsotropicInputs(particle,local,&kin,&b0,&variance,&brms);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(!(b0>0.0)) return Failure(local,c,StatusCode::OutsideModelDomain,
                              "Kuhlen field-line fit I11 requires B0>0");
  if(!local.turbulence.correlationLengthM ||
     !Positive(*local.turbulence.correlationLengthM))
    return Failure(local,c,StatusCode::MissingInput,
                   "Kuhlen model requires a tagged positive correlationLengthM");
  const double lc=*local.turbulence.correlationLengthM;
  const double radius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*brms);
  const double x=radius/lc;
  const double parallel=kin.speedMPerS*lc*Parameter(c,"A")*std::pow(x,1.0/3.0)*
      std::pow(1.0+std::pow(x/Parameter(c,"rho_star"),
                           5.0/(3.0*Parameter(c,"s_kappa"))),
               Parameter(c,"s_kappa"));
  const double tau=3.0*parallel/(kin.speedMPerS*kin.speedMPerS);
  const double gamma=Parameter(c,"gamma_K");
  const double z1=Parameter(c,"z1_m"),z2=Parameter(c,"z2_m");
  const double ck=Parameter(c,"C_K");
  const double bxRatio=variance/(3.0*b0*b0);
  auto dLine=[&](double z) {
    if(z==0.0) return 0.0;
    return ck*z*bxRatio*
      std::pow(1.0+std::pow(z/z1,(1.0-gamma)/1.5),-1.5)*
      std::pow(1.0+std::pow(z/z2,-gamma/0.2),0.2);
  };
  auto attached=[&](double time) {
    if(!(time>0.0)) return 0.0;
    const double denominator=std::hypot(time,tau);
    const double dParallel=parallel*time/denominator;
    const double mz=2.0*parallel*time*time/(denominator+tau);
    const double z=std::sqrt(mz);
    return z>0.0?dLine(z)/z*dParallel:0.0;
  };
  const double transverse=Parameter(c,"transverse_correlation_length_m");
  auto crossing=[&](double logTime) {
    const double time=std::exp(logTime);
    return 2.0*time*attached(time)-transverse*transverse;
  };
  // Scan logarithmically for every negative-to-positive crossing. The first
  // upward crossing is selected only because that exact rule and calibration
  // identity are mandatory configuration fields; it is never an implicit default.
  const double center=std::log(tau);
  std::vector<std::pair<double,double> > brackets;
  double previousLog=center-40.0,previous=crossing(previousLog);
  for(int i=1;i<=800;++i) {
    const double currentLog=center-40.0+80.0*i/800.0;
    const double current=crossing(currentLog);
    if(Finite(previous)&&Finite(current)&&previous<=0.0&&current>0.0)
      brackets.push_back({previousLog,currentLog});
    previousLog=currentLog; previous=current;
  }
  if(brackets.empty()) return Failure(local,c,StatusCode::NonlinearSolverFailed,
                                     "Kuhlen decorrelation condition has no upward crossing in the searched positive-time range");
  double low=brackets.front().first,high=brackets.front().second;
  const double tolerance=OptionalParameter(c,"relative_tolerance",1.0e-9);
  const int maximum=static_cast<int>(OptionalParameter(c,"maximum_iterations",100.0));
  int iterations=0;
  for(;iterations<maximum;++iterations) {
    const double mid=0.5*(low+high),f=crossing(mid);
    if(std::fabs(f)<=tolerance*transverse*transverse||high-low<=tolerance) {low=high=mid;break;}
    if(f>0.0) high=mid; else low=mid;
  }
  if(iterations==maximum) return Failure(local,c,StatusCode::NonlinearSolverFailed,
                                        "Kuhlen decorrelation root did not converge");
  const double tc=std::exp(0.5*(low+high));
  const double perpendicular=attached(tc);
  ModelResult result=UnequalCoefficient(local,c,parallel,perpendicular,perpendicular,
      DependencyOwner::PairedBackend,Quality::ParameterizedClosure,
      DomainState::NotApplicable);
  if(!result.status.ok()) return result;
  result.provenance.calibrationId=TextParameter(c,"calibration_id");
  result.provenance.conventions["field_line_coordinate"]="arc_length_sigma_small_turbulence_z_approximation";
  result.provenance.conventions["root_selection"]=TextParameter(c,"root_selection");
  result.numerical.rootMethod="log_time_scan_then_bisection_first_upward_crossing";
  result.numerical.iterations=static_cast<std::size_t>(iterations+1);
  result.numerical.controls["candidate_upward_crossings"]=static_cast<double>(brackets.size());
  result.numerical.controls["decorrelation_time_s"]=tc;
  if(c.numeric.count("age_s")) {
    const double age=Parameter(c,"age_s");
    if(!Nonnegative(age)) return Failure(local,c,StatusCode::InvalidConfiguration,
                                        "Kuhlen age_s must be nonnegative");
    const double used=std::min(age,tc);
    const double denominator=std::hypot(used,tau);
    const double mz=used==0.0?0.0:2.0*parallel*used*used/(denominator+tau);
    const double z=std::sqrt(mz);
    auto integral=[&](double end) {
      if(end==0.0) return 0.0;
      const int n=400;
      const double h=end/n;
      double sum=dLine(0.0)+dLine(end);
      for(int i=1;i<n;++i) sum+=(i%2?4.0:2.0)*dLine(i*h);
      return sum*h/3.0;
    };
    double mx=2.0*integral(z);
    if(age>tc) mx+=2.0*perpendicular*(age-tc);
    ParticleMoments moments; moments.ageS=age;
    moments.rawMsdM2={{mx,mx,age==0.0?0.0:
        2.0*parallel*age*age/(std::hypot(age,tau)+tau)}};
    moments.derivativeM2PerS=std::array<double,3>{{
        age<=tc?attached(age):perpendicular,
        age<=tc?attached(age):perpendicular,
        age==0.0?0.0:parallel*age/std::hypot(age,tau)}};
    result.particleMoments=moments;
  }
  return result;
}

namespace {

ModelResult EvaluateRbd(const ParticleState& particle, const LocalState& local,
                        const ModelConfiguration& c) {
  // B3 is evaluated for the normalized S12 pure-2D area spectrum. V_x is the
  // physical velocity variance [m^2 s^-2] from B2, and erfc implements the
  // backtracking correction B9. Replacing it with erfcx would be a distinct
  // uncorrected model. This is an explicit quadrature, never a root solve.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double kp=0.0,b0=0.0;
  status=Parallel(local,&kp); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  status=MeanField(local,&b0,nullptr); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(local.turbulence.geometry!=GeometryKind::Pure2D ||
     !local.turbulence.twoDVarianceT2 || !local.turbulence.twoDBendoverLengthM)
    return Failure(local,c,StatusCode::IncompatibleGeometry,
                   "RBD B3 implementation requires a pure-2D S12 area spectrum");
  const double variance=*local.turbulence.twoDVarianceT2;
  const double length=*local.turbulence.twoDBendoverLengthM;
  if(!Nonnegative(variance)||!Positive(length))
    return Failure(local,c,StatusCode::InvalidInput,
                   "RBD requires nonnegative variance and positive break scale");
  if(variance==0.0 || kp==0.0) return EqualCoefficient(local,c,0.0,Quality::ParameterizedClosure);
  const double a2=Parameter(c,"a_squared");
  const double ratio=variance/(b0*b0);
  if(a2*ratio>1.0)
    return Failure(local,c,StatusCode::OutsideModelDomain,
                   "RBD velocity variance B2 is negative for a_squared*variance/B0^2>1");
  const double vx=a2*kin.speedMPerS*kin.speedMPerS*ratio/6.0;
  if(vx==0.0) return EqualCoefficient(local,c,0.0,Quality::ParameterizedClosure);
  const double qRate=kin.speedMPerS*kin.speedMPerS/(3.0*kp);
  const double p=Parameter(c,"area_energy_index");
  const double sA=Parameter(c,"area_inertial_index");
  const double c2=(sA-1.0)*(p+2.0)/(2.0*Pi*(p+sA+1.0));
  auto integrand=[&](double y) {
    const double x=std::exp(y);
    const double shape=std::exp((y<=0.0?p:-sA-1.0)*y);
    return shape*std::erfc(qRate*length/(x*std::sqrt(2.0*vx)))*x;
  };
  double integral=0.0,error=0.0;
  status=IntegrateLog(integrand,c,&integral,&error);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  // B3: S2=C2*variance*lambda^2*shape, dk=dx/lambda.
  const double kappa=a2*kin.speedMPerS*kin.speedMPerS/(6.0*b0*b0)*
      std::sqrt(Pi/2.0)*2.0*Pi*c2*variance*length/std::sqrt(vx)*integral;
  ModelResult result=EqualCoefficient(local,c,kappa,Quality::ParameterizedClosure);
  result.numerical.backend="RBD_backtracking_corrected_B3";
  result.numerical.integrationMethod="adaptive_simpson_log_wavenumber_with_tail_refinement";
  result.numerical.relativeErrorEstimate=error;
  result.provenance.conventions["temporal_kernel"]="erfc_backtracking_corrected";
  return result;
}

ModelResult EvaluateCandia(const ParticleState& particle,
                           const LocalState& local,
                           const ModelConfiguration& c, bool drift) {
  // I4/I5 and G10 use the ordered-field gyroradius and sigma^2=deltaB^2/B0^2.
  // The relativistic source fit retains c rather than substituting particle
  // speed. Pair and Hall domains are separately configured and reported.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double b0=0.0; status=MeanField(local,&b0,nullptr);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(local.turbulence.geometry!=GeometryKind::Isotropic3D)
    return Failure(local,c,StatusCode::IncompatibleGeometry,
                   "Candia-Roulet requires isotropic_3d turbulence geometry");
  if(!local.turbulence.totalVarianceT2)
    return Failure(local,c,StatusCode::MissingInput,
                   "Candia-Roulet requires isotropic totalVarianceT2");
  const double variance=*local.turbulence.totalVarianceT2;
  if(!Positive(variance)) return Failure(local,c,StatusCode::OutsideModelDomain,
                                        "Candia-Roulet requires positive turbulence variance");
  FitRow row;
  if(!CandiaRow(TextParameter(c,"spectrum"),&row))
    return Failure(local,c,StatusCode::InvalidConfiguration,
                   "unknown Candia-Roulet spectrum table row");
  const double length=Parameter(c,"outer_scale_m");
  const double radius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*b0);
  const double rho=radius/length;
  const double sigma2=variance/(b0*b0);
  if(drift) {
    const bool inside=rho>=Parameter(c,"rho_min")&&rho<=Parameter(c,"rho_max")&&
        sigma2>=Parameter(c,"sigma2_min")&&sigma2<=Parameter(c,"sigma2_max");
    status=Status::Success();
    const DomainState domain=DomainPolicy(inside,c,&status,"Hall-fit input lies outside its declared supplied domain");
    if(!status.ok()) return Failure(local,c,status.code,status.detail);
    const double sigma0=row.nA*(rho<=0.2?std::pow(rho,0.3):1.9*std::pow(rho,0.7));
    ModelResult result=BaseResult(local,c);
    result.observable=Observable::SignedHallCoefficient;
    result.estimator=Estimator::Asymptotic;
    result.quality=Quality::SourceFit;
    result.domain=domain;
    result.signedHallM2PerS=std::copysign(1.0,particle.chargeC)*
        SpeedOfLightMPerS*radius/3.0/std::sqrt(1.0+std::pow(sigma2/sigma0,2.0));
    result.provenance.conventions["domain_identity"]=TextParameter(c,"domain_identity");
    return result;
  }
  const bool inside=rho>=0.03&&rho<=1.0&&sigma2>=1.0&&sigma2<=10.0;
  status=Status::Success();
  const DomainState domain=DomainPolicy(inside,c,&status,
      "Candia-Roulet pair lies outside 0.03<=rho<=1 and 1<=sigma^2<=10");
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double parallel=SpeedOfLightMPerS*length*rho*row.nParallel/sigma2*
      std::sqrt(std::pow(rho/row.rhoParallel,2.0*(1.0-row.g))+
                std::pow(rho/row.rhoParallel,2.0));
  double ratio=row.nPerp*std::pow(sigma2,row.aPerp);
  if(rho>0.2) ratio*=std::pow(rho/0.2,-2.0);
  return UnequalCoefficient(local,c,parallel,parallel*ratio,parallel*ratio,
                            DependencyOwner::PairedBackend,Quality::SourceFit,domain);
}

Status IsotropicInputs(const ParticleState& particle, const LocalState& local,
                       ParticleKinematics* kin, double* b0, double* variance,
                       double* totalField) {
  Status status=Kinematics(particle,kin); if(!status.ok()) return status;
  if(local.turbulence.geometry!=GeometryKind::Isotropic3D)
    return Status::Error(StatusCode::IncompatibleGeometry,
                         "selected isotropic fit requires isotropic_3d geometry");
  if(!local.meanFieldT || !local.turbulence.totalVarianceT2)
    return Status::Error(StatusCode::MissingInput,
                         "isotropic fit requires meanFieldT and totalVarianceT2");
  *b0=VectorNorm(*local.meanFieldT); *variance=*local.turbulence.totalVarianceT2;
  if(!Nonnegative(*b0)||!Positive(*variance))
    return Status::Error(StatusCode::InvalidInput,
                         "isotropic fit requires finite B0>=0 and variance>0");
  *totalField=std::sqrt((*b0)*(*b0)+*variance);
  return Status::Success();
}

bool SnodinRow(const std::string& spectrum,double* a1,double* a2) {
  if(spectrum=="sharp_kolmogorov") {*a1=0.0031;*a2=0.74;}
  else if(spectrum=="sharp_kraichnan") {*a1=0.0019;*a2=0.76;}
  else if(spectrum=="peaked_kolmogorov") {*a1=0.0017;*a2=0.75;}
  else return false;
  return true;
}

ModelResult EvaluateSnodin(const ParticleState& particle,
                           const LocalState& local,
                           const ModelConfiguration& c) {
  // I6/I7 instead use the total-rms-field radius and the source particle
  // speed. Only the sharp-cutoff Kolmogorov row owns the ordered-field pair;
  // other rows are restricted to the zero-mean isotropic kappa0 observable.
  ParticleKinematics kin; double b0=0.0,variance=0.0,brms=0.0;
  Status status=IsotropicInputs(particle,local,&kin,&b0,&variance,&brms);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double length=Parameter(c,"outer_scale_m");
  const double radius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*brms);
  const double x=radius/length;
  const double eta=variance/(b0*b0+variance);
  double a1=0.0,a2=0.0;
  const std::string spectrum=c.model==ModelId::IsoFitSnodin2016 ?
      "sharp_kolmogorov":TextParameter(c,"spectrum");
  if(!SnodinRow(spectrum,&a1,&a2))
    return Failure(local,c,StatusCode::InvalidConfiguration,"unknown Snodin spectrum row");
  const double kappa0=kin.speedMPerS*length*(a1+a2*x);
  const bool inside=x>=0.004&&x<=0.05&&eta>=0.1&&eta<=1.0;
  status=Status::Success();
  const DomainState domain=DomainPolicy(inside,c,&status,
      "Snodin fit lies outside 0.004<=x<=0.05 and 0.1<=eta_B<=1");
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(c.model==ModelId::IsoKappa0Snodin2016 && b0!=0.0)
    return Failure(local,c,StatusCode::IncompatibleGeometry,
                   "iso_kappa0_snodin_2016 is the zero-mean-field isotropic backend");
  if(c.model==ModelId::IsoKappa0Snodin2016 || eta==1.0)
    return IsotropicCoefficient(local,c,kappa0,Quality::SourceFit,domain);
  const double parallel=kin.speedMPerS*length*(a1+a2*x+
      std::pow(x,1.0/3.0)*(1.0-eta)/(3.0*eta));
  const double perp=kappa0/(1.0+2.35*(1.0-eta)/eta);
  return UnequalCoefficient(local,c,parallel,perp,perp,
                            DependencyOwner::PairedBackend,Quality::SourceFit,domain);
}

ModelResult EvaluateClassical(const ParticleState& particle,
                              const LocalState& local,
                              const ModelConfiguration& c, bool hallOnly) {
  // I14 uses x=lambda_parallel/r_L for one isotropic scattering time. The
  // signed Hall coefficient shares the charge-sign convention of G8, while
  // field polarity remains in the returned ordered frame direction.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double b0=0.0,kp=0.0;
  status=MeanField(local,&b0,nullptr); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  status=Parallel(local,&kp); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double radius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*b0);
  const double lambda=3.0*kp/kin.speedMPerS;
  const double x=lambda/radius;
  if(hallOnly) {
    ModelResult result=BaseResult(local,c);
    result.observable=Observable::SignedHallCoefficient;
    result.quality=Quality::ParameterizedClosure;
    result.signedHallM2PerS=std::copysign(1.0,particle.chargeC)*kp*x/(1.0+x*x);
    result.provenance.dependencyOwner=DependencyOwner::SuppliedParallel;
    result.provenance.parallelModelId=local.parallelDependency->modelId;
    return result;
  }
  ModelResult result=UnequalCoefficient(local,c,kp,kp/(1.0+x*x),kp/(1.0+x*x),
      DependencyOwner::SuppliedParallel,Quality::ParameterizedClosure,
      DomainState::NotApplicable);
  if(result.status.ok())
    result.signedHallM2PerS=std::copysign(1.0,particle.chargeC)*kp*x/(1.0+x*x);
  return result;
}

ModelResult EvaluateDrift(const ParticleState& particle,
                          const LocalState& local,
                          const ModelConfiguration& c) {
  // G8/G9 return the antisymmetric scalar only. Curl construction, current
  // sheets, and the gradient term in G11 remain transport-operator work; the
  // library must not turn this scalar into a drift trajectory by itself.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double b0=0.0; status=MeanField(local,&b0,nullptr);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double radius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*b0);
  double hall=std::copysign(1.0,particle.chargeC)*kin.speedMPerS*radius/3.0;
  if(c.model==ModelId::DriftRigidityReduction) {
    const double x=kin.rigidityV/Parameter(c,"rigidity_A_V");
    hall*=Parameter(c,"K_A0")*x*x/(1.0+x*x);
  }
  ModelResult result=BaseResult(local,c);
  result.observable=Observable::SignedHallCoefficient;
  result.estimator=Estimator::Asymptotic;
  result.quality=Quality::ParameterizedClosure;
  result.signedHallM2PerS=hall;
  return result;
}

}  // namespace

namespace {

ModelResult EvaluateEnlgc(const ParticleState& particle,
                          const LocalState& local,
                          const ModelConfiguration& c) {
  // N5 is solved in the dimensionless coefficient ratio eta. For q=0 its
  // reduced integral is O(1), avoiding dimensional underflow from Q at solar
  // particle speeds. The analytic upper bound brackets the positive root.
  ParticleKinematics kin;
  Status status = Kinematics(particle, &kin);
  if (!status.ok()) return Failure(local,c,status.code,status.detail);
  double kp = 0.0;
  status = Parallel(local, &kp);
  if (!status.ok()) return Failure(local,c,status.code,status.detail);
  CompositeState t;
  status = Composite(local, &t, false, true, false);
  if (!status.ok()) return Failure(local,c,status.code,status.detail);
  if (t.q != 0.0)
    return Failure(local,c,StatusCode::IncompatibleGeometry,
                   "ENLGC N5 requires the exactly flat q=0 2D spectrum");
  if (kp == 0.0 || t.twoDVariance == 0.0)
    return EqualCoefficient(local,c,0.0,Quality::ParameterizedClosure);
  const double upperRatio = t.twoDVariance / (2.0*t.b0*t.b0);
  const double lambdaParallel=3.0*kp/kin.speedMPerS;
  double integrationError = 0.0;
  auto residual = [&](double logRatio) {
    const double ratio=std::exp(logRatio);
    const double a=lambdaParallel*lambdaParallel*ratio/
                   (3.0*t.twoDLength*t.twoDLength);
    auto integrand=[&](double y) {
      const double x=std::exp(y);
      return std::exp(y-0.5*t.s*LogOnePlusSquare(y))/(1.0+a*x*x);
    };
    double integral=0.0,error=0.0;
    const Status integrated=IntegrateLog(integrand,c,&integral,&error);
    if (!integrated.ok()) return std::numeric_limits<double>::quiet_NaN();
    integrationError=std::max(integrationError,error);
    const double rhs=2.0*SpectrumC(t.s)*t.twoDVariance/(t.b0*t.b0)*integral;
    return Positive(rhs) ? std::log(rhs)-logRatio
                         : std::numeric_limits<double>::quiet_NaN();
  };
  double ratio=0.0,residualNorm=0.0; std::size_t iterations=0;
  status=PositiveRoot(residual,upperRatio,c,&ratio,&iterations,&residualNorm);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  ModelResult result=EqualCoefficient(local,c,ratio*kp,Quality::ParameterizedClosure);
  result.numerical.backend="ENLGC_N5";
  result.numerical.integrationMethod="adaptive_simpson_log_wavenumber_with_tail_refinement";
  result.numerical.rootMethod="bracketed_logarithmic_bisection";
  result.numerical.relativeErrorEstimate=integrationError;
  result.numerical.residualNorm=residualNorm;
  result.numerical.iterations=iterations;
  return result;
}

ModelResult EvaluateNlgcePair(const ParticleState& particle,
                              const LocalState& local,
                              const ModelConfiguration& c) {
  // NLGCE is a pair, never a perpendicular closure attached to an unrelated
  // parallel value. The map below performs no repartitioning or scale
  // conversion and delegates E1--E9 to the parallel library's existing owner.
  if (!local.meanFieldT || !local.turbulence.slabVarianceT2 ||
      !local.turbulence.twoDVarianceT2 ||
      !local.turbulence.slabBendoverLengthM ||
      !local.turbulence.twoDBendoverLengthM)
    return Failure(local,c,StatusCode::MissingInput,
        "NLGCE requires mean field, both variances, and both bend-over lengths");

  // The pair has one numerical owner: reuse the already qualified shared
  // parallel library backend rather than copying E1--E9 into a second solver.
  // This mapping preserves total two-component variances and bend-over scales.
  ParallelDiffusion::ParticleState pp;
  pp.massKg=particle.massKg; pp.chargeC=particle.chargeC;
  pp.momentumKgMPerS=particle.momentumKgMPerS;
  ParallelDiffusion::LocalState pl;
  pl.timeS=local.timeS; pl.positionM=local.positionM;
  pl.meanFieldT=local.meanFieldT;
  ParallelDiffusion::LocalState::TurbulenceState turbulence;
  turbulence.slabVarianceT2=local.turbulence.slabVarianceT2;
  turbulence.twoDVarianceT2=local.turbulence.twoDVarianceT2;
  turbulence.slabBendoverLengthM=local.turbulence.slabBendoverLengthM;
  turbulence.twoDBendoverLengthM=local.turbulence.twoDBendoverLengthM;
  turbulence.inertialIndex=5.0/3.0;
  pl.turbulence=turbulence;
  ParallelDiffusion::ModelConfiguration pc;
  pc.model=c.model==ModelId::NlgceN ? ParallelDiffusion::ModelId::NlgceN
                                    : ParallelDiffusion::ModelId::NlgceF2014;
  pc.nonlinear.numerical.relativeTolerance=OptionalParameter(c,"relative_tolerance",1.0e-8);
  pc.nonlinear.numerical.maximumRefinements=static_cast<int>(
      OptionalParameter(c,"maximum_refinements",18.0));
  pc.nonlinear.numerical.maximumIterations=static_cast<int>(
      OptionalParameter(c,"maximum_iterations",80.0));
  const ParallelDiffusion::ParallelResult pair=
      ParallelDiffusion::Evaluate(pp,pl,pc);
  if(!pair.status.ok() || !pair.kappaParallelM2PerS || !pair.perpendicular)
    return Failure(local,c,StatusCode::NonlinearSolverFailed,
        "shared NLGCE pair backend failed: "+pair.status.detail);
  LocalState pairOwnedState=local;
  pairOwnedState.parallelDependency.reset();
  ModelResult result=EqualCoefficient(pairOwnedState,c,
      pair.perpendicular->kappaPerpendicularM2PerS,Quality::SourceFit,
      DomainState::InsideDeclaredDomain);
  if(!result.status.ok()) return result;
  result.observable=Observable::PairedCoefficients;
  result.coefficients->parallelM2PerS=*pair.kappaParallelM2PerS;
  result.coefficients->parallelOwner=DependencyOwner::PairedBackend;
  result.provenance.dependencyOwner=DependencyOwner::PairedBackend;
  result.provenance.parallelModelId=ParallelDiffusion::ModelName(pc.model);
  result.provenance.conventions["paired_backend_software_version"]=
      pair.provenance.softwareVersion;
  result.provenance.conventions["paired_coefficient_set_identity"]=
      pair.provenance.coefficientSetIdentity;
  if(pair.nonlinearIterations)
    result.numerical.iterations=static_cast<std::size_t>(*pair.nonlinearIterations);
  result.numerical.residualNorm=pair.nonlinearMaxLogResidual;
  result.numerical.backend="shared_parallel_diffusion_NLGCE_pair_owner";
  return result;
}

ModelResult EvaluateFiniteTime(const ParticleState&, const LocalState& local,
                               const ModelConfiguration& c) {
  // F8, C3--C6, and N6 are finite-time statistics. Their derivative and
  // secant estimators differ by stated factors and remain tagged; an
  // established subdiffusive asymptote is success with no normal diffusion.
  double kp=0.0;
  Status status=Parallel(local,&kp);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double age=Parameter(c,"age_s");
  if(!Positive(age)) return Failure(local,c,StatusCode::InvalidInput,
                                    "finite-time diagnostic requires age_s>0");
  ModelResult result=BaseResult(local,c);
  result.observable=Observable::ParticleMsd;
  result.quality=Quality::Diagnostic;
  result.diffusionRegime=DiffusionRegime::NoNormalDiffusion;
  ParticleMoments moments; moments.ageS=age;
  if(c.model==ModelId::CompoundDiffusiveLines ||
     c.model==ModelId::EnlgcSlabSecantDiagnostic) {
    double fieldLine=0.0;
    if(c.model==ModelId::CompoundDiffusiveLines) status=SelectedFieldLine(local,&fieldLine);
    else {
      CompositeState t; status=Composite(local,&t,true,false,false);
      if(status.ok()) fieldLine=SlabFieldLine(t);
    }
    if(!status.ok()) return Failure(local,c,status.code,status.detail);
    const double msd=4.0*fieldLine*std::sqrt(kp*age/Pi);
    const double derivative=fieldLine*std::sqrt(kp/(Pi*age));
    const double secant=2.0*derivative;
    moments.rawMsdM2={{msd,msd,2.0*kp*age}};
    moments.derivativeM2PerS=std::array<double,3>{{derivative,derivative,kp}};
    moments.secantM2PerS=std::array<double,3>{{secant,secant,kp}};
    result.estimator=c.model==ModelId::EnlgcSlabSecantDiagnostic ?
        Estimator::Secant : Estimator::MomentValue;
    result.particleMoments=moments;
    result.provenance.conventions["field_line_model_id"]=
        c.model==ModelId::CompoundDiffusiveLines ? local.fieldLineModelId
                                                 : "F2_slab";
    return result;
  }

  CompositeState t;
  status=Composite(local,&t,false,true,false);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(t.q!=0.0)
    return Failure(local,c,StatusCode::IncompatibleGeometry,
                   "GCD C3-C6 requires the flat q=0 2D spectrum");
  const double alpha=std::tgamma(7.0/6.0)/std::sqrt(Pi)*
      std::pow(18.0*SpectrumC(t.s)*std::sqrt(Pi/2.0),2.0/3.0);
  const double amplitude=alpha*std::pow(t.twoDVariance/(t.b0*t.b0),2.0/3.0)*
      std::pow(t.twoDLength,2.0/3.0)*std::pow(2.0*kp,2.0/3.0);
  const double msd=amplitude*std::pow(age,2.0/3.0);
  const double derivative=amplitude/3.0*std::pow(age,-1.0/3.0);
  const double secant=amplitude/2.0*std::pow(age,-1.0/3.0);
  moments.rawMsdM2={{msd,msd,2.0*kp*age}};
  moments.derivativeM2PerS=std::array<double,3>{{derivative,derivative,kp}};
  moments.secantM2PerS=std::array<double,3>{{secant,secant,kp}};
  result.estimator=Estimator::MomentValue;
  result.particleMoments=moments;
  return result;
}

ModelResult EvaluatePreDiffusive(const LocalState& local,
                                 const ModelConfiguration& c) {
  // C10 is a conditional median for returning gyrocentres in the injection
  // plane. It is not an ensemble MSD, covariance, source width, or Markov
  // diffusion coefficient, so only conditionalMedianM2 is populated.
  const double age=Parameter(c,"age_s");
  const double a1=Parameter(c,"A1_m2"), period=Parameter(c,"gyroperiod_s");
  const double t1=Parameter(c,"t1_s"),t2=Parameter(c,"t2_s");
  const double alpha=Parameter(c,"alpha_m"),beta=Parameter(c,"beta_m");
  const double value=a1*std::pow(age/period,alpha)*
      (1.0+std::pow(age/t1,beta-alpha))/(1.0+std::pow(age/t2,beta-1.0));
  ModelResult result=BaseResult(local,c);
  result.observable=Observable::ConditionalReturnMedian;
  result.estimator=Estimator::MomentValue;
  result.quality=Quality::SourceFit;
  result.diffusionRegime=DiffusionRegime::NotEstablished;
  result.conditionalMedianM2=value;
  result.provenance.calibrationId=TextParameter(c,"calibration_id");
  result.provenance.conventions["statistic"]=
      "returning_particle_injection_plane_squared_gyrocentre_median";
  return result;
}

ModelResult EvaluateCompositeClosed(const ParticleState& particle,
                                    const LocalState& local,
                                    const ModelConfiguration& c) {
  // B6 is evaluated in its rationalized form to avoid cancellation when the
  // field-line/parallel product is small. The supplied length keeps its source
  // profile tag; no bend-over/integral conversion is performed here.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  double kp=0.0,fieldLine=0.0;
  status=Parallel(local,&kp); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  status=SelectedFieldLine(local,&fieldLine); if(!status.ok()) return Failure(local,c,status.code,status.detail);
  if(kp==0.0 || fieldLine==0.0) return EqualCoefficient(local,c,0.0,Quality::ParameterizedClosure);
  const double lp=3.0*kp/kin.speedMPerS;
  const double length=Parameter(c,"perpendicular_length_m");
  const double radical=std::sqrt(1.0+8.0*fieldLine*lp/(3.0*length*length));
  const double lambda=4.0*fieldLine*fieldLine*lp/(length*length)/
                      ((radical+1.0)*(radical+1.0));
  ModelResult result=EqualCoefficient(local,c,kin.speedMPerS*lambda/3.0,
                                      Quality::ParameterizedClosure);
  result.provenance.conventions["perpendicular_length_profile"]=
      TextParameter(c,"length_profile");
  result.provenance.conventions["field_line_model_id"]=local.fieldLineModelId;
  return result;
}

ModelResult EvaluatePerturbativeFgr(const ParticleState& particle,
                                    const LocalState& local,
                                    const ModelConfiguration& c) {
  // U10 is only a weak finite-gyroradius diagnostic. Positivity is necessary
  // but insufficient, so the caller supplies a stricter maximum correction;
  // breakdown is reported rather than clipped to a nonnegative coefficient.
  ParticleKinematics kin; Status status=Kinematics(particle,&kin);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  CompositeState t; status=Composite(local,&t,false,true,true);
  if(!status.ok()) return Failure(local,c,status.code,status.detail);
  const double gyroradius=particle.momentumKgMPerS/(std::fabs(particle.chargeC)*t.b0);
  const double correction=(t.q-1.0)/(8.0*(t.s-1.0))*
                          std::pow(gyroradius/t.twoDLength,2.0);
  if(!(correction<Parameter(c,"maximum_relative_correction")))
    return Failure(local,c,StatusCode::OutsideModelDomain,
                   "finite-gyroradius correction exceeds configured perturbative domain");
  const double lambda=1.5*TwoDFieldLine(t)*(1.0-correction);
  if(!Positive(lambda)) return Failure(local,c,StatusCode::OutsideModelDomain,
                                      "perturbative U10 bracket is not positive");
  ModelResult result=EqualCoefficient(local,c,kin.speedMPerS*lambda/3.0,
                                      Quality::Diagnostic);
  result.observable=Observable::PerturbativeDiagnostic;
  result.numerical.controls["relative_correction"]=correction;
  return result;
}

}  // namespace

double SpectrumC(double s) {
  if (!(s > 1.0)) return std::numeric_limits<double>::quiet_NaN();
  return std::tgamma(s / 2.0) /
      (2.0 * std::sqrt(Pi) * std::tgamma((s - 1.0) / 2.0));
}

double SpectrumD(double s, double q) {
  if (!(s > 1.0) || !(q > -1.0))
    return std::numeric_limits<double>::quiet_NaN();
  return std::tgamma((s + q) / 2.0) /
      (2.0 * std::tgamma((s - 1.0) / 2.0) *
       std::tgamma((q + 1.0) / 2.0));
}

Status SmoothTwoDSpectralMoment(double s, double q, int m,
                                double variance, double length,
                                double* value) {
  if (!value || !Nonnegative(variance) || !Positive(length) || !(s > 1.0) ||
      !(q > -1.0))
    return Status::Error(StatusCode::InvalidInput,
                         "invalid smooth-spectrum moment input");
  if (!(q + m > -1.0) || !(s > m + 1.0))
    return Status::Error(StatusCode::DivergentMoment,
                         "requested uncut smooth-spectrum moment diverges");
  *value = SpectrumD(s, q) / Pi * variance * std::pow(length, -m) *
      std::tgamma((q + m + 1.0) / 2.0) *
      std::tgamma((s - m - 1.0) / 2.0) / std::tgamma((s + q) / 2.0);
  return Status::Success();
}

Status SmoothSpectrumLengths(double s, double q, double length,
                             double* ultra, double* integral) {
  if (!ultra || !integral || !Positive(length) || !(s > 1.0))
    return Status::Error(StatusCode::InvalidInput,
                         "invalid smooth-spectrum length input");
  if (!(q > 1.0))
    return Status::Error(StatusCode::DivergentMoment,
                         "ultra-scale requires q>1");
  *ultra = std::sqrt((s - 1.0) / (q - 1.0)) * length;
  *integral = 2.0 * std::tgamma(q / 2.0) * std::tgamma(s / 2.0) /
      (std::tgamma((q + 1.0) / 2.0) *
       std::tgamma((s - 1.0) / 2.0)) * length;
  return Status::Success();
}

Status EvaluateImplicitSlabKernel(double xi, bool rational, double* kernel) {
  if (!kernel || !Nonnegative(xi))
    return Status::Error(StatusCode::InvalidInput,
                         "implicit-slab kernel requires xi>=0");
  if (rational) {
    *kernel = 1.0 / (1.0 + 2.0 * xi * xi);
    return Status::Success();
  }
  if (xi == 0.0) { *kernel = 1.0; return Status::Success(); }
  if (xi >= 12.0) {
    const double inverse2 = 1.0 / (xi * xi);
    // U7 is alternating. Terms are retained only while decreasing, avoiding
    // the catastrophic subtraction in 1-sqrt(pi)*xi*erfcx(xi).
    *kernel = 0.5 * inverse2 - 0.75 * inverse2 * inverse2 +
              1.875 * inverse2 * inverse2 * inverse2;
  } else {
    // erfcx(x)=exp(x^2)erfc(x) is safe in this bounded branch.
    *kernel = 1.0 - SqrtPi * xi * std::exp(xi * xi) * std::erfc(xi);
  }
  if (!Positive(*kernel))
    return Status::Error(StatusCode::IntegrationFailed,
                         "implicit-slab kernel lost positive precision");
  return Status::Success();
}

namespace Internal {

namespace {

ModelResult FieldLineResult(const ParticleState&, const LocalState& local,
                            const ModelConfiguration& c) {
  // F2 uses the two-sided slab normalization, F3 the finite I_-2 2D moment,
  // and F4 the nonlinear composite quadratic. Variances are total transverse
  // component variances [T^2], so the per-axis factors are already in F2/F3.
  const bool slabOnly = c.model == ModelId::FieldLineSlab;
  const bool twoDOnly = c.model == ModelId::FieldLine2D;
  const bool zeroSlab=local.turbulence.slabVarianceT2 &&
                      *local.turbulence.slabVarianceT2==0.0;
  const bool zeroTwoD=local.turbulence.twoDVarianceT2 &&
                      *local.turbulence.twoDVarianceT2==0.0;
  if((slabOnly&&zeroSlab)||(twoDOnly&&zeroTwoD)||
     (!slabOnly&&!twoDOnly&&zeroSlab&&zeroTwoD)) {
    double b0=0.0;
    const Status field=MeanField(local,&b0,nullptr);
    if(!field.ok()) return Failure(local,c,field.code,field.detail);
    ModelResult result=BaseResult(local,c);
    result.observable=Observable::FieldLineCoefficient;
    result.estimator=Estimator::Asymptotic;
    result.quality=Quality::ParameterizedClosure;
    result.diffusionRegime=DiffusionRegime::NormalDiffusion;
    result.fieldLineM=0.0;
    result.provenance.conventions["line_coordinate"]="mean_field_z";
    return result;
  }
  CompositeState t;
  const Status status = Composite(local, &t, !twoDOnly, !slabOnly,
                                  !slabOnly);
  if (!status.ok()) return Failure(local, c, status.code, status.detail);
  double coefficient = 0.0;
  if (slabOnly) coefficient = SlabFieldLine(t);
  else if (twoDOnly) coefficient = TwoDFieldLine(t);
  else coefficient = CompositeFieldLine(t);
  ModelResult result = BaseResult(local, c);
  result.observable = Observable::FieldLineCoefficient;
  result.estimator = Estimator::Asymptotic;
  result.quality = Quality::ParameterizedClosure;
  result.diffusionRegime = DiffusionRegime::NormalDiffusion;
  result.fieldLineM = coefficient;
  result.provenance.conventions["line_coordinate"] = "mean_field_z";
  result.provenance.conventions["variance"] =
      "total_two_component_transverse_magnetic_variance";
  return result;
}

ModelResult EvaluatePrescribed(const ParticleState& particle,
                               const LocalState& local,
                               const ModelConfiguration& c) {
  // P1--P7 are evaluated literally in SI. Configuration supplies every
  // normalization/exponent, while runtime geometry supplies radius, B0,
  // spiral cosine and mu. Pitch-angle results are never exposed as an
  // isotropically averaged coefficient unless that mode is selected.
  ParticleKinematics k;
  Status status = Kinematics(particle, &k);
  if (!status.ok()) return Failure(local, c, status.code, status.detail);
  double parallel = 0.0;
  switch (c.model) {
    case ModelId::ConstantKappaPerp:
      return EqualCoefficient(local, c, Parameter(c, "kappa_perp_m2_per_s"),
                              Quality::ParameterizedClosure);
    case ModelId::ConstantLambdaPerp:
      return EqualCoefficient(local, c,
          k.speedMPerS * Parameter(c, "lambda_perp_m") / 3.0,
          Quality::ParameterizedClosure);
    case ModelId::RatioKappa:
      status = Parallel(local, &parallel);
      if (!status.ok()) return Failure(local, c, status.code, status.detail);
      return EqualCoefficient(local, c, Parameter(c, "eta_kappa") * parallel,
                              Quality::ParameterizedClosure);
    case ModelId::PowerLawPerp: {
      double b0 = 0.0;
      status = MeanField(local, &b0, nullptr);
      if (!status.ok()) return Failure(local, c, status.code, status.detail);
      const double radius = VectorNorm(local.positionM);
      if (!Positive(radius))
        return Failure(local, c, StatusCode::OutsideModelDomain,
                       "power-law radius factor requires r>0");
      const double value = Parameter(c, "kappa0_m2_per_s") *
          std::pow(k.beta, Parameter(c, "beta_exponent")) *
          std::pow(k.rigidityV / Parameter(c, "rigidity0_V"),
                   Parameter(c, "rigidity_exponent")) *
          std::pow(Parameter(c, "field_reference_T") / b0,
                   Parameter(c, "field_exponent")) *
          std::pow(radius / Parameter(c, "radius0_m"),
                   Parameter(c, "radius_exponent"));
      return EqualCoefficient(local, c, value, Quality::ParameterizedClosure);
    }
    case ModelId::PitchAnglePerp: {
      if (!particle.mu.has_value() || !Finite(*particle.mu) ||
          std::fabs(*particle.mu) > 1.0)
        return Failure(local, c, StatusCode::MissingInput,
                       "pitch-angle model requires mu in [-1,1]");
      const double mu = *particle.mu;
      double shape = 1.0;
      if (TextParameter(c, "shape") == "abs_mu") shape = 2.0 * std::fabs(mu);
      else if (TextParameter(c, "shape") == "sqrt_one_minus_mu2")
        shape = 4.0 / Pi * std::sqrt(std::max(0.0, 1.0 - mu * mu));
      ModelResult result = BaseResult(local, c);
      result.observable = Observable::PitchAngleCoefficient;
      result.estimator = Estimator::NotApplicable;
      result.quality = Quality::ParameterizedClosure;
      result.pitchAngleM2PerS = Parameter(c, "D0_m2_per_s") * shape;
      return result;
    }
    case ModelId::DrogeLambdaPerp: {
      if (!particle.mu.has_value() || !local.parkerSpiralCosine.has_value())
        return Failure(local, c, StatusCode::MissingInput,
                       "Droge-Dresing requires mu and Parker-spiral cosine");
      if (std::fabs(*particle.mu) > 1.0 ||
          !Finite(*local.parkerSpiralCosine) ||
          *local.parkerSpiralCosine < 0.0 || *local.parkerSpiralCosine > 1.0)
        return Failure(local, c, StatusCode::InvalidInput,
                       "invalid mu or Parker-spiral cosine");
      status = Parallel(local, &parallel);
      if (!status.ok()) return Failure(local, c, status.code, status.detail);
      const double radius = VectorNorm(local.positionM);
      if (!Positive(radius))
        return Failure(local, c, StatusCode::OutsideModelDomain,
                       "Droge-Dresing requires r>0");
      const double lambdaParallel = 3.0 * parallel / k.speedMPerS;
      const double lambda = Parameter(c, "alpha_D") * lambdaParallel *
          std::pow(radius / Parameter(c, "radius0_m"), 2.0) *
          *local.parkerSpiralCosine *
          std::sqrt(std::max(0.0, 1.0 - (*particle.mu) * (*particle.mu)));
      ModelResult result = BaseResult(local, c);
      result.observable = Observable::PitchAngleCoefficient;
      result.quality = Quality::SourceFit;
      result.pitchAngleM2PerS = k.speedMPerS * lambda / 3.0;
      result.provenance.dependencyOwner = DependencyOwner::SuppliedParallel;
      result.provenance.parallelModelId = local.parallelDependency->modelId;
      return result;
    }
    case ModelId::ParadiseAlpha: {
      double b0 = 0.0;
      status = MeanField(local, &b0, nullptr);
      if (!status.ok()) return Failure(local, c, status.code, status.detail);
      status = Parallel(local, &parallel);
      if (!status.ok()) return Failure(local, c, status.code, status.detail);
      return EqualCoefficient(local, c, Pi / 4.0 * Parameter(c, "alpha_P") *
          parallel * Parameter(c, "field_reference_T") / b0,
          Quality::SourceFit);
    }
    case ModelId::FlrwParticle: {
      double fieldLine = 0.0;
      status = SelectedFieldLine(local, &fieldLine);
      if (!status.ok()) return Failure(local, c, status.code, status.detail);
      const double a = Parameter(c, "a_FL");
      if (TextParameter(c, "pitch_angle_mode") == "local_mu") {
        if (!particle.mu.has_value() || std::fabs(*particle.mu) > 1.0)
          return Failure(local, c, StatusCode::MissingInput,
                         "local FLRW mode requires mu in [-1,1]");
        ModelResult result = BaseResult(local, c);
        result.observable = Observable::PitchAngleCoefficient;
        result.quality = Quality::ParameterizedClosure;
        result.pitchAngleM2PerS = a * k.speedMPerS *
                                 std::fabs(*particle.mu) * fieldLine;
        result.provenance.conventions["field_line_model_id"]=
            local.fieldLineModelId;
        return result;
      }
      ModelResult result=EqualCoefficient(local,c,
          0.5*a*k.speedMPerS*fieldLine,Quality::ParameterizedClosure);
      result.provenance.conventions["field_line_model_id"]=local.fieldLineModelId;
      return result;
    }
    default:
      return Failure(local, c, StatusCode::UnsupportedModel,
                     "internal prescribed-model dispatch mismatch");
  }
}

ModelResult EvaluateImplicitTwoD(const ParticleState& particle,
                                 const LocalState& local,
                                 const ModelConfiguration& c) {
  // N2/U2/U5/D2 share a positive dimensionless spectral solve but retain
  // distinct denominators and physics. eta=kappa_perp/kappa_parallel is the
  // unknown, x=k*ell is integrated on the log coordinate, and each analytic
  // Section 16.2 bound supplies the bracket. No asymptotic model switching,
  // numerical spectral cutoff, clipping, or fallback is performed.
  ParticleKinematics kin;
  Status status = Kinematics(particle, &kin);
  if (!status.ok()) return Failure(local, c, status.code, status.detail);
  double kp = 0.0;
  status = Parallel(local, &kp);
  if (!status.ok()) return Failure(local, c, status.code, status.detail);
  CompositeState t;
  const bool nlgc = c.model == ModelId::Nlgc ||
                    c.model == ModelId::NlgcSlabKernelDiagnostic;
  const bool exactSlab = c.model == ModelId::ImplicitSlabExact2016 ||
                         c.model == ModelId::ImplicitSlabRational2016;
  const bool flpd = c.model == ModelId::FlpdComplete;
  if (c.model == ModelId::Nlgc &&
      local.turbulence.geometry == GeometryKind::PureSlab)
    return Failure(local, c, StatusCode::IncompatibleGeometry,
        "production NLGC rejects the spurious pure-slab asymptotic branch");
  if (c.model == ModelId::Unlt &&
      local.turbulence.geometry == GeometryKind::PureSlab)
    return EqualCoefficient(local, c, 0.0, Quality::ParameterizedClosure);
  const bool slabDiagnostic = c.model == ModelId::NlgcSlabKernelDiagnostic;
  const bool needSlab = slabDiagnostic;
  const bool needTwoD = !slabDiagnostic;
  const bool needQ = flpd;
  status = Composite(local, &t, needSlab, needTwoD, needQ);
  if (!status.ok()) return Failure(local, c, status.code, status.detail);
  if (!needSlab && local.turbulence.slabVarianceT2) {
    t.slabVariance=*local.turbulence.slabVarianceT2;
    if(!Nonnegative(t.slabVariance))
      return Failure(local,c,StatusCode::InvalidInput,
                     "optional slab component is invalid");
    if(t.slabVariance>0.0) {
      if(!local.turbulence.slabBendoverLengthM)
        return Failure(local,c,StatusCode::MissingInput,
                       "nonzero optional slab component requires its bend-over length");
      t.slabLength=*local.turbulence.slabBendoverLengthM;
      if(!Positive(t.slabLength))
        return Failure(local,c,StatusCode::InvalidInput,
                       "optional slab bend-over length must be positive");
    }
  }
  if (local.turbulence.geometry==GeometryKind::CompositeSlab2D &&
      !local.turbulence.slabVarianceT2)
    return Failure(local,c,StatusCode::MissingInput,
                   "composite geometry requires its declared slab component");
  if (c.model == ModelId::NlgcSlabKernelDiagnostic &&
      local.turbulence.geometry != GeometryKind::PureSlab)
    return Failure(local, c, StatusCode::IncompatibleGeometry,
        "NLGC slab-kernel diagnostic requires pure_slab geometry");
  if(c.model==ModelId::Nlgc && t.twoDVariance==0.0)
    return Failure(local,c,StatusCode::IncompatibleGeometry,
                   "production NLGC requires nonzero transverse 2D power");
  if (t.twoDVariance == 0.0 || kp == 0.0) {
    // NLGC retains a nonzero slab kernel when kappa_parallel>0. At exactly
    // zero parallel input every listed continuous branch is zero.
    if (kp == 0.0 || !nlgc || t.slabVariance == 0.0)
      return EqualCoefficient(local, c, 0.0,
          c.model == ModelId::NlgcSlabKernelDiagnostic ? Quality::Diagnostic
                                                       : Quality::ParameterizedClosure);
  }
  const double varianceRatio2 = t.twoDVariance / (t.b0 * t.b0);
  const double totalVariance = t.slabVariance + t.twoDVariance;
  const double a2 = (c.model == ModelId::Nlgc ||
                     c.model == ModelId::NlgcSlabKernelDiagnostic ||
                     c.model == ModelId::Unlt) ? Parameter(c, "a_squared") : 1.0;
  double upperRatio = a2 * totalVariance / (2.0 * t.b0 * t.b0);
  if (!nlgc) upperRatio = a2 * varianceRatio2 / 2.0;
  if (upperRatio == 0.0) return EqualCoefficient(local, c, 0.0,
      c.model == ModelId::NlgcSlabKernelDiagnostic ? Quality::Diagnostic
                                                   : Quality::ParameterizedClosure);

  const double lambdaParallel = 3.0 * kp / kin.speedMPerS;
  const double slabFieldLine = t.slabVariance>0.0 ? SlabFieldLine(t) : 0.0;
  const double compositeFieldLine = flpd ? CompositeFieldLine(t) : 0.0;
  double lastIntegralError = 0.0;
  auto rhsRatio = [&](double eta, bool* ok) {
    *ok = false;
    double predicted = 0.0;
    const double lambdaPerp = eta * lambdaParallel;
    if (nlgc && t.slabVariance > 0.0) {
      const double as=lambdaParallel*lambdaParallel/
                      (3.0*t.slabLength*t.slabLength);
      auto slabIntegrand = [&](double y) {
        const double x = std::exp(y);
        return std::exp(y - 0.5 * t.s * LogOnePlusSquare(y)) /
            (1.0 + as*x*x);
      };
      double slabIntegral = 0.0, error = 0.0;
      status = IntegrateLog(slabIntegrand, c, &slabIntegral, &error);
      if (!status.ok()) return 0.0;
      // N4 with q=0: H=4*C(s)*integral for the slab shape.
      predicted += 2.0*a2*SpectrumC(t.s)*
                   t.slabVariance/(t.b0*t.b0)*slabIntegral;
      lastIntegralError = std::max(lastIntegralError, error);
    }
    if (t.twoDVariance > 0.0) {
      const double a=lambdaParallel*lambdaPerp/
                     (3.0*t.twoDLength*t.twoDLength);
      auto twoDIntegrand = [&](double y) {
        const double x = std::exp(y);
        const double shapeJacobian = std::exp(SmoothShapeLog(y,t.s,t.q) + y);
        double denominator = 1.0+a*x*x;
        double kernel = 1.0;
        if (c.model == ModelId::Unlt)
          denominator = 1.0+4.0*a*x*x/3.0;
        else if (exactSlab) {
          const double xi = slabFieldLine * lambdaParallel * x * x /
              (std::sqrt(3.0 * Pi) * t.twoDLength * t.twoDLength *
               std::sqrt(1.0+a*x*x));
          if (!EvaluateImplicitSlabKernel(xi,
                  c.model == ModelId::ImplicitSlabRational2016, &kernel).ok())
            return std::numeric_limits<double>::quiet_NaN();
        } else if (flpd) {
          denominator += compositeFieldLine/t.twoDLength/std::sqrt(eta)*x*
                         std::sqrt(1.0+a*x*x);
        }
        return shapeJacobian * kernel / denominator;
      };
      double twoDIntegral = 0.0, error = 0.0;
      status = IntegrateLog(twoDIntegrand, c, &twoDIntegral, &error);
      if (!status.ok()) return 0.0;
      predicted += 2.0*a2*SpectrumD(t.s,t.q)*varianceRatio2*twoDIntegral;
      lastIntegralError = std::max(lastIntegralError, error);
    }
    *ok = true;
    return predicted;
  };
  auto residualFunction = [&](double logRatio) {
    bool ok = false;
    const double eta = std::exp(logRatio);
    const double value = rhsRatio(eta, &ok);
    return ok && Positive(value) ? std::log(value) - logRatio
                                 : std::numeric_limits<double>::quiet_NaN();
  };
  double ratio = 0.0, residual = 0.0;
  std::size_t iterations = 0;
  status = PositiveRoot(residualFunction, upperRatio, c, &ratio, &iterations, &residual);
  if (!status.ok()) return Failure(local, c, status.code, status.detail);
  ModelResult result = EqualCoefficient(local, c, ratio*kp,
      c.model == ModelId::NlgcSlabKernelDiagnostic ? Quality::Diagnostic
                                                   : Quality::ParameterizedClosure);
  if (result.status.ok()) {
    result.numerical.backend = "positive_logarithmic_implicit_closure";
    result.numerical.integrationMethod = "adaptive_simpson_log_wavenumber_with_tail_refinement";
    result.numerical.rootMethod = "bracketed_logarithmic_bisection";
    result.numerical.relativeErrorEstimate = lastIntegralError;
    result.numerical.residualNorm = residual;
    result.numerical.iterations = iterations;
    if(flpd)
      result.provenance.conventions["field_line_closure"]="F4_dd_composite";
  }
  return result;
}

}  // namespace

ModelResult EvaluateModel(const ParticleState& particle,
                          const LocalState& local,
                          const ModelConfiguration& configuration) {
  // A single exhaustive dispatcher makes registry coverage reviewable. The
  // five incomplete identities remain explicit typed failures even though
  // normal Evaluate() rejects them during configuration validation first.
  switch(configuration.model) {
    case ModelId::ConstantKappaPerp:
    case ModelId::ConstantLambdaPerp:
    case ModelId::RatioKappa:
    case ModelId::PowerLawPerp:
    case ModelId::PitchAnglePerp:
    case ModelId::DrogeLambdaPerp:
    case ModelId::ParadiseAlpha:
    case ModelId::FlrwParticle:
      return EvaluatePrescribed(particle,local,configuration);
    case ModelId::FieldLineSlab:
    case ModelId::FieldLine2D:
    case ModelId::FieldLineComposite:
      return FieldLineResult(particle,local,configuration);
    case ModelId::Nlgc:
    case ModelId::NlgcSlabKernelDiagnostic:
    case ModelId::Unlt:
    case ModelId::ImplicitSlabExact2016:
    case ModelId::ImplicitSlabRational2016:
    case ModelId::FlpdComplete:
      return EvaluateImplicitTwoD(particle,local,configuration);
    case ModelId::Enlgc2D:
      return EvaluateEnlgc(particle,local,configuration);
    case ModelId::NlgceN:
    case ModelId::NlgceF2014:
      return EvaluateNlgcePair(particle,local,configuration);
    case ModelId::RbdBc:
      return EvaluateRbd(particle,local,configuration);
    case ModelId::CompositeClosed2019:
      return EvaluateCompositeClosed(particle,local,configuration);
    case ModelId::CompoundDiffusiveLines:
    case ModelId::GcdCompound:
    case ModelId::EnlgcSlabSecantDiagnostic:
      return EvaluateFiniteTime(particle,local,configuration);
    case ModelId::PrediffusiveFit:
      return EvaluatePreDiffusive(local,configuration);
    case ModelId::IsoFitCandiaRoulet2004:
      return EvaluateCandia(particle,local,configuration,false);
    case ModelId::IsoFitSnodin2016:
    case ModelId::IsoKappa0Snodin2016:
      return EvaluateSnodin(particle,local,configuration);
    case ModelId::IsoFitKuhlen2025:
      return EvaluateKuhlen(particle,local,configuration);
    case ModelId::ClassicalScattering:
      return EvaluateClassical(particle,local,configuration,false);
    case ModelId::YanLazarianMa4:
      return EvaluateScaling(particle,local,configuration);
    case ModelId::NwuRatioPolar:
    case ModelId::CortiAms02:
    case ModelId::HelmodRatio:
      return EvaluateGcrPrescription(particle,local,configuration);
    case ModelId::DriftWeakScattering:
    case ModelId::DriftRigidityReduction:
      return EvaluateDrift(particle,local,configuration);
    case ModelId::DriftClassical:
      return EvaluateClassical(particle,local,configuration,true);
    case ModelId::DriftCandiaRoulet2004:
      return EvaluateCandia(particle,local,configuration,true);
    case ModelId::TabulatedPerp:
      return EvaluateTabulated(particle,local,configuration);
    case ModelId::UnltFgrPerturbativeDiagnostic:
      return EvaluatePerturbativeFgr(particle,local,configuration);
    case ModelId::FrozenFieldline:
    case ModelId::QltPerp:
    case ModelId::UnltFgr:
    case ModelId::IsoFitCasse2002:
    case ModelId::RestrictedScattering:
      return Failure(local,configuration,StatusCode::SourceGate,
                     "selected registry identity remains source gated");
  }
  return Failure(local,configuration,StatusCode::UnsupportedModel,
                 "invalid perpendicular diffusion model identifier");
}

}  // namespace Internal
}  // namespace PerpendicularDiffusion
}  // namespace SEP
