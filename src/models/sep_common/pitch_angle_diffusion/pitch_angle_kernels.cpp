// Pure scattering kernels and the finite-volume reference operator; see
// pitch_angle_kernels.h for the contract and equation references.

#include "pitch_angle_kernels.h"
#include "pitch_angle_quadrature.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace PitchAngleDiffusion {
namespace Kernels {

double FoldPitchCosine(double mu, unsigned* reflections) {
  if (reflections) *reflections = 0;
  if (!std::isfinite(mu) || (mu >= -1.0 && mu <= 1.0)) return mu;
  // Repeated application of Eq. (161) (mu -> 2 - mu above 1, -2 - mu below
  // -1) is the triangular wave of period 4 through (-1,-1), (1,1), (3,-1).
  // Folding in closed form keeps the cost O(1) for an arbitrarily large kick.
  const double cell = std::floor((mu + 1.0) / 2.0);
  if (reflections) *reflections = static_cast<unsigned>(std::fabs(cell));
  double folded = std::fmod(mu + 1.0, 4.0);
  if (folded < 0.0) folded += 4.0;
  return folded <= 2.0 ? folded - 1.0 : 3.0 - folded;
}

namespace {

// Validation of an Ito kick.  The Status (with its message) is built only on
// failure; the success path costs a few comparisons.
inline bool KickInputValid(double mu, double d, double dPrime, double rate,
                           double dt, double zeta) {
  return std::isfinite(mu) && mu >= -1.0 && mu <= 1.0 && std::isfinite(d) &&
         d >= 0.0 && std::isfinite(dPrime) && std::isfinite(rate) &&
         std::isfinite(dt) && dt > 0.0 && std::isfinite(zeta);
}

Status CheckKickInput(double mu, double d, double dPrime, double rate,
                      double dt, double zeta) {
  if (!std::isfinite(mu) || mu < -1.0 || mu > 1.0)
    return Status::Error(StatusCode::InvalidInput, "mu must lie in [-1, 1]");
  if (!std::isfinite(d) || d < 0.0)
    return Status::Error(StatusCode::InvalidInput,
                         "D_mumu must be finite and nonnegative");
  if (!std::isfinite(dPrime))
    return Status::Error(StatusCode::IncompatibleKernel,
        "an Ito kick needs a finite dD/dmu (Section 16.11)");
  if (!std::isfinite(rate))
    return Status::Error(StatusCode::InvalidInput,
                         "deterministic pitch-angle rate is not finite");
  if (!std::isfinite(dt) || !(dt > 0.0))
    return Status::Error(StatusCode::Nonprogress,
        "a scattering kick needs a positive interval; zero-time success is "
        "forbidden (Section 17)");
  if (!std::isfinite(zeta))
    return Status::Error(StatusCode::InvalidInput, "normal variate not finite");
  return Status::Success();
}

}  // namespace

KickResult EulerMaruyamaKick(double mu, double dMuMuPerS,
                             double dDmuMuDmuPerS,
                             double deterministicMuRatePerS, double dtS,
                             double standardNormal) {
  KickResult result;
  if (!KickInputValid(mu, dMuMuPerS, dDmuMuDmuPerS, deterministicMuRatePerS,
                      dtS, standardNormal)) {
    result.status = CheckKickInput(mu, dMuMuPerS, dDmuMuDmuPerS,
                                   deterministicMuRatePerS, dtS, standardNormal);
    return result;
  }
  result.driftIncrement = (deterministicMuRatePerS + dDmuMuDmuPerS) * dtS;
  result.noiseIncrement = std::sqrt(2.0 * dMuMuPerS * dtS) * standardNormal;
  result.mu = FoldPitchCosine(mu + result.driftIncrement + result.noiseIncrement,
                              &result.reflections);
  return result;
}

KickResult MilsteinKick(double mu, double dMuMuPerS, double dDmuMuDmuPerS,
                        double deterministicMuRatePerS, double dtS,
                        double standardNormal) {
  KickResult result;
  if (!KickInputValid(mu, dMuMuPerS, dDmuMuDmuPerS, deterministicMuRatePerS,
                      dtS, standardNormal)) {
    result.status = CheckKickInput(mu, dMuMuPerS, dDmuMuDmuPerS,
                                   deterministicMuRatePerS, dtS, standardNormal);
    return result;
  }
  result.driftIncrement = (deterministicMuRatePerS + dDmuMuDmuPerS) * dtS;
  result.noiseIncrement = std::sqrt(2.0 * dMuMuPerS * dtS) * standardNormal;
  // The Milstein correction (1/2) b b' dt (zeta^2 - 1) with b b' = D'.  It
  // is applied only where D > 0, where the identity holds.
  if (dMuMuPerS > 0.0)
    result.noiseIncrement +=
        0.5 * dDmuMuDmuPerS * dtS * (standardNormal * standardNormal - 1.0);
  result.mu = FoldPitchCosine(mu + result.driftIncrement + result.noiseIncrement,
                              &result.reflections);
  return result;
}

KickResult RandomRotationKick(double mu, double nuPerS, double dtS,
                              double uniformTheta, double uniformPhi) {
  KickResult result;
  if (!std::isfinite(mu) || mu < -1.0 || mu > 1.0) {
    result.status = Status::Error(StatusCode::InvalidInput,
                                  "mu must lie in [-1, 1]");
    return result;
  }
  if (!std::isfinite(nuPerS) || nuPerS < 0.0) {
    result.status = Status::Error(StatusCode::InvalidInput,
                                  "nu must be finite and nonnegative");
    return result;
  }
  if (!std::isfinite(dtS) || !(dtS > 0.0)) {
    result.status = Status::Error(StatusCode::Nonprogress,
                                  "a rotation kick needs a positive interval");
    return result;
  }
  if (!(uniformTheta > 0.0) || !(uniformTheta <= 1.0) ||
      !(uniformPhi >= 0.0) || !(uniformPhi < 1.0)) {
    result.status = Status::Error(StatusCode::InvalidInput,
        "rotation variates must satisfy 0 < uTheta <= 1, 0 <= uPhi < 1");
    return result;
  }
  // Eq. (167): theta^2 ~ Exp(mean 2 nu dt).
  const double theta = std::sqrt(-2.0 * nuPerS * dtS * std::log(uniformTheta));
  const double phi = 2.0 * 3.14159265358979323846 * uniformPhi;
  const double perpendicular = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  // Eq. (168); the field-aligned component of a rotated unit vector.
  const double rotated = mu * std::cos(theta) +
                         perpendicular * std::sin(theta) * std::cos(phi);
  result.mu = FoldPitchCosine(rotated, &result.reflections);  // roundoff only
  result.noiseIncrement = result.mu - mu;
  result.status = Status::Success();
  return result;
}

// ---------------------------------------------------------------------------
// Finite-volume reference.
// ---------------------------------------------------------------------------
namespace {

// Transmissibility integrals use the shared global-adaptive quadrature
// (pitch_angle_quadrature.h) with a relative tolerance of 1e-13 of the whole
// integral and a 4000-segment budget.
template <class F>
bool Integrate(const F& f, long double a, long double b, long double* value,
               long double* error) {
  return Quadrature::AdaptiveIntegrate(f, {a, b}, 1.0e-13L, 1.0e-300L, 4000,
                                       value, error);
}

double ShapeHValue(ModelId model, const ShapeParameters& shape, double mu) {
  const double a = std::fabs(mu);
  switch (model) {
    case ModelId::Isotropic:
    case ModelId::Eprem: return 1.0;
    case ModelId::QForm: return std::pow(a, shape.q - 1.0) + shape.gapH;
    case ModelId::EpsilonForm: return a / (1.0 + a) + shape.epsilon;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// R(a,b) = Int_a^b dmu/((1-mu^2) h(mu)) for -1 < a < b < 1 (the D0 factor is
// applied by the caller).  Returns +inf for a disconnected (eps = 0) crossing.
bool ReciprocalIntegral(ModelId model, const ShapeParameters& shape, double a,
                        double b, double* value, double* relativeError) {
  *relativeError = 0.0;
  if (model == ModelId::Isotropic || model == ModelId::Eprem) {
    // Int dmu/(1-mu^2) = atanh(b) - atanh(a), exact.
    *value = std::atanh(b) - std::atanh(a);
    return true;
  }
  if (model == ModelId::EpsilonForm) {
    if (shape.epsilon == 0.0 && a <= 0.0 && b >= 0.0) {
      // h ~ |mu| near 0: Int dmu/|mu| diverges (alpha = 1, Section 16.5).
      *value = std::numeric_limits<double>::infinity();
      return true;
    }
    auto f = [&](long double mu) {
      return 1.0L / ((1.0L - mu * mu) *
                     static_cast<long double>(ShapeHValue(model, shape,
                                                          static_cast<double>(mu))));
    };
    long double total = 0.0L, totalError = 0.0L;
    std::vector<long double> points = {a, b};
    if (a < 0.0 && b > 0.0) points = {a, 0.0L, b};  // kink of h at 0
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
      long double piece = 0.0L, pieceError = 0.0L;
      if (!Integrate(f, points[i], points[i + 1], &piece, &pieceError))
        return false;
      total += piece;
      totalError += pieceError;
    }
    *value = static_cast<double>(total);
    *relativeError = static_cast<double>(totalError / std::fabs(total));
    return true;
  }
  // PA-QFORM: h = |mu|^alpha + H, alpha = q - 1 in (0,1).  On each side of 0
  // use mu = t^m, m = 1/(1 - alpha): dmu/mu^alpha = m dt (derived), so the
  // integrand m t^(m alpha)/((1 - t^(2m))(t^(m alpha) + H)) is bounded at
  // t = 0 for every H >= 0 (an integrable |mu|^-alpha resistance stays
  // finite: the hemispheres communicate, Section 16.5).
  const long double alpha = shape.q - 1.0;
  const long double m = 1.0L / (1.0L - alpha);
  const long double ma = m * alpha;
  const long double gap = shape.gapH;
  auto g = [m, ma, gap](long double t) {
    const long double tp = std::pow(t, ma);
    return m * tp / ((1.0L - std::pow(t, 2.0L * m)) * (tp + gap));
  };
  auto onePositive = [&](long double lo, long double hi, long double* out,
                         long double* err) {
    // Int_lo^hi over mu in [0,1) in the t variable.
    const long double tlo = std::pow(lo, 1.0L / m), thi = std::pow(hi, 1.0L / m);
    return Integrate(g, tlo, thi, out, err);
  };
  long double total = 0.0L, totalError = 0.0L, piece = 0.0L, pieceError = 0.0L;
  if (a < 0.0 && b > 0.0) {
    if (!onePositive(0.0L, -static_cast<long double>(a), &piece, &pieceError))
      return false;
    total += piece; totalError += pieceError;
    if (!onePositive(0.0L, b, &piece, &pieceError)) return false;
    total += piece; totalError += pieceError;
  } else {
    const long double lo = std::min(std::fabs(a), std::fabs(b));
    const long double hi = std::max(std::fabs(a), std::fabs(b));
    if (!onePositive(lo, hi, &piece, &pieceError)) return false;
    total = piece; totalError = pieceError;
  }
  *value = static_cast<double>(total);
  *relativeError = static_cast<double>(totalError / std::fabs(total));
  return true;
}

Status CheckEdges(const std::vector<double>& edges) {
  if (edges.size() < 3)
    return Status::Error(StatusCode::InvalidInput,
                         "the finite-volume grid needs at least two cells");
  if (edges.front() != -1.0 || edges.back() != 1.0)
    return Status::Error(StatusCode::InvalidInput,
                         "the finite-volume grid must span exactly [-1, 1]");
  for (std::size_t i = 0; i + 1 < edges.size(); ++i)
    if (!std::isfinite(edges[i]) || !(edges[i + 1] > edges[i]))
      return Status::Error(StatusCode::InvalidInput,
                           "finite-volume edges must be strictly increasing");
  return Status::Success();
}

}  // namespace

Status BuildFvGeneratorFromConductances(const std::vector<double>& edges,
                                        const std::vector<double>& masses,
                                        const std::vector<double>& conductance,
                                        FvGenerator* out) {
  if (out == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null generator output");
  const Status grid = CheckEdges(edges);
  if (!grid.ok()) return grid;
  const std::size_t n = edges.size() - 1;
  if (masses.size() != n || conductance.size() != n - 1)
    return Status::Error(StatusCode::InvalidInput,
                         "masses/conductances do not match the grid");
  FvGenerator g;
  g.edges = edges;
  g.stationaryMass = masses;
  g.conductance = conductance;
  g.centres.resize(n);
  g.rateLeft.assign(n, 0.0);
  g.rateRight.assign(n, 0.0);
  g.rateDiagonal.assign(n, 0.0);
  bool disconnected = false;
  for (std::size_t i = 0; i < n; ++i) {
    g.centres[i] = 0.5 * (edges[i] + edges[i + 1]);
    if (!std::isfinite(masses[i]) || !(masses[i] > 0.0))
      return Status::Error(StatusCode::InvalidInput,
                           "stationary masses must be positive");
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    if (!std::isfinite(conductance[i]) || conductance[i] < 0.0)
      return Status::Error(StatusCode::InvalidInput,
                           "conductances must be finite and nonnegative");
    if (conductance[i] == 0.0) disconnected = true;
    g.rateRight[i] = conductance[i] / masses[i];
    g.rateLeft[i + 1] = conductance[i] / masses[i + 1];
  }
  for (std::size_t i = 0; i < n; ++i)
    g.rateDiagonal[i] = -(g.rateLeft[i] + g.rateRight[i]);
  g.internalExtension = disconnected ? "disconnected (zero conductance)"
                                     : "connected";
  *out = g;
  return Status::Success();
}

Status BuildFvGenerator(ModelId model, const ShapeParameters& shape,
                        double amplitudeD0PerS,
                        const std::vector<double>& edges, FvGenerator* out) {
  if (out == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null generator output");
  const Status grid = CheckEdges(edges);
  if (!grid.ok()) return grid;
  if (!std::isfinite(amplitudeD0PerS) || !(amplitudeD0PerS > 0.0))
    return Status::Error(StatusCode::InvalidInput,
                         "the amplitude D0 must be positive");
  // Validate the shape through the coefficient layer (H = 0, eps = 0 allowed).
  const CoefficientResult probe =
      EvaluateShape(0.5, amplitudeD0PerS, model, shape);
  if (!probe.status.ok()) return probe.status;
  const std::size_t n = edges.size() - 1;
  std::vector<double> masses(n), conductance(n - 1);
  double worst = 0.0;
  for (std::size_t i = 0; i < n; ++i) masses[i] = edges[i + 1] - edges[i];
  for (std::size_t i = 0; i + 1 < n; ++i) {
    const double a = 0.5 * (edges[i] + edges[i + 1]);
    const double b = 0.5 * (edges[i + 1] + edges[i + 2]);
    double resistance = 0.0, error = 0.0;
    if (!ReciprocalIntegral(model, shape, a, b, &resistance, &error))
      return Status::Error(StatusCode::Nonconverged,
                           "transmissibility quadrature did not converge");
    worst = std::max(worst, error);
    conductance[i] = std::isinf(resistance)
        ? 0.0 : amplitudeD0PerS / resistance;  // T = [Int dmu/D]^-1
  }
  const Status built =
      BuildFvGeneratorFromConductances(edges, masses, conductance, out);
  if (!built.ok()) return built;
  out->quadratureRelativeError = worst;
  if (out->internalExtension != "connected")
    out->internalExtension = "disconnected at mu = 0 (Int dmu/D diverges)";
  return Status::Success();
}

namespace {

// Implicit QL with Wilkinson shifts for a symmetric tridiagonal matrix
// (diagonal d, subdiagonal e[1..n-1]); z receives the eigenvectors (columns),
// starting from the identity.  Standard algorithm (Golub & Van Loan 8.3).
bool TridiagonalEigen(std::vector<double>& d, std::vector<double>& e,
                      std::vector<double>& z, std::size_t n) {
  for (std::size_t i = 1; i < n; ++i) e[i - 1] = e[i];
  e[n - 1] = 0.0;
  for (std::size_t l = 0; l < n; ++l) {
    int iterations = 0;
    std::size_t m;
    do {
      for (m = l; m + 1 < n; ++m) {
        const double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
        if (std::fabs(e[m]) <= std::numeric_limits<double>::epsilon() * dd) break;
      }
      if (m != l) {
        if (++iterations > 60) return false;
        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = std::hypot(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? std::fabs(r) : -std::fabs(r)));
        double s = 1.0, c = 1.0, p = 0.0;
        std::size_t i;
        bool underflow = false;
        for (i = m; i-- > l;) {
          double f = s * e[i];
          const double b = c * e[i];
          r = std::hypot(f, g);
          e[i + 1] = r;
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            underflow = true;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          for (std::size_t k = 0; k < n; ++k) {
            f = z[k * n + i + 1];
            z[k * n + i + 1] = s * z[k * n + i] + c * f;
            z[k * n + i] = c * z[k * n + i] - s * f;
          }
        }
        if (underflow) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
  return true;
}

}  // namespace

Status PrepareFvPropagator(const FvGenerator& generator, FvPropagator* out) {
  if (out == nullptr)
    return Status::Error(StatusCode::InvalidInput, "null propagator output");
  const std::size_t n = generator.centres.size();
  if (n < 2) return Status::Error(StatusCode::InvalidInput, "empty generator");
  FvPropagator p;
  p.n = n;
  p.sqrtMass.resize(n);
  for (std::size_t i = 0; i < n; ++i)
    p.sqrtMass[i] = std::sqrt(generator.stationaryMass[i]);
  // S = M^-1/2 K M^-1/2 with K = M Q: S_ii = Q_ii, S_{i,i+1} = T/sqrt(m_i m_j).
  std::vector<double> diagonal(n), sub(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) diagonal[i] = generator.rateDiagonal[i];
  for (std::size_t i = 1; i < n; ++i)
    sub[i] = generator.conductance[i - 1] / (p.sqrtMass[i - 1] * p.sqrtMass[i]);
  p.eigenvectors.assign(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) p.eigenvectors[i * n + i] = 1.0;
  if (!TridiagonalEigen(diagonal, sub, p.eigenvectors, n))
    return Status::Error(StatusCode::Nonconverged,
                         "tridiagonal eigen-decomposition did not converge");
  p.eigenvalues = diagonal;
  *out = p;
  return Status::Success();
}

Status PropagateFv(const FvPropagator& propagator,
                   const std::vector<double>& initialMass, double timeS,
                   std::vector<double>* mass) {
  const std::size_t n = propagator.n;
  if (mass == nullptr || initialMass.size() != n)
    return Status::Error(StatusCode::InvalidInput, "mass vector size mismatch");
  if (!std::isfinite(timeS) || timeS < 0.0)
    return Status::Error(StatusCode::InvalidInput, "time must be nonnegative");
  // p(t) = p0 M^-1/2 V exp(Lambda t) V^T M^1/2 (row-vector convention).
  std::vector<long double> coefficient(n, 0.0L);
  for (std::size_t k = 0; k < n; ++k) {
    long double sum = 0.0L;
    for (std::size_t i = 0; i < n; ++i)
      sum += static_cast<long double>(initialMass[i]) / propagator.sqrtMass[i] *
             propagator.eigenvectors[i * n + k];
    // Eigenvalues are <= 0; a roundoff-positive zero mode is exactly 0.
    coefficient[k] = sum * std::exp(std::min(0.0, propagator.eigenvalues[k]) *
                                    timeS);
  }
  mass->assign(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    long double sum = 0.0L;
    for (std::size_t k = 0; k < n; ++k)
      sum += propagator.eigenvectors[j * n + k] * coefficient[k];
    (*mass)[j] = static_cast<double>(sum * propagator.sqrtMass[j]);
  }
  return Status::Success();
}

Status FvAsymptoticVariance(const FvGenerator& generator,
                            const std::vector<double>& cellValues,
                            double* variance) {
  const std::size_t n = generator.centres.size();
  if (variance == nullptr || cellValues.size() != n)
    return Status::Error(StatusCode::InvalidInput, "cell-value size mismatch");
  long double massTotal = 0.0L, mean = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    massTotal += generator.stationaryMass[i];
    mean += generator.stationaryMass[i] * static_cast<long double>(cellValues[i]);
  }
  mean /= massTotal;
  // -K phi = m o f with K tridiagonal of conductances: the flux through face
  // i+1/2 equals the cumulative source, T_{i,i+1}(phi_i - phi_{i+1}) =
  // sum_{k<=i} m_k (f_k - fbar) (derived; zero flux at mu = +-1).
  std::vector<long double> phi(n, 0.0L);
  long double cumulative = 0.0L;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    cumulative += generator.stationaryMass[i] *
                  (static_cast<long double>(cellValues[i]) - mean);
    if (generator.conductance[i] == 0.0)
      return Status::Error(StatusCode::OutsideDomain,
          "the chain is disconnected; the asymptotic variance is undefined");
    phi[i + 1] = phi[i] - cumulative / generator.conductance[i];
  }
  long double sum = 0.0L;
  for (std::size_t i = 0; i < n; ++i)
    sum += generator.stationaryMass[i] *
           (static_cast<long double>(cellValues[i]) - mean) * phi[i];
  *variance = static_cast<double>(sum / massTotal);
  return Status::Success();
}

double Legendre(int l, double x) {
  if (l == 0) return 1.0;
  double previous = 1.0, current = x;
  for (int k = 1; k < l; ++k) {
    const double next = ((2.0 * k + 1.0) * x * current - k * previous) / (k + 1.0);
    previous = current;
    current = next;
  }
  return current;
}

double LegendreCellAverage(int l, double left, double right) {
  if (l == 0) return 1.0;
  // Int P_l = (P_{l+1} - P_{l-1})/(2l+1).
  const double upper = Legendre(l + 1, right) - Legendre(l - 1, right);
  const double lower = Legendre(l + 1, left) - Legendre(l - 1, left);
  return (upper - lower) / ((2.0 * l + 1.0) * (right - left));
}

}  // namespace Kernels
}  // namespace PitchAngleDiffusion
}  // namespace SEP
