// ============================================================================
// Statistical acceptance campaign of the pitch-angle library (roadmap
// Section 26.5 acceptance cases; stages PS04, PS04a, PS05, PS06, PS08).
//
// Every case integrates an ensemble of pseudo-particles with the library's
// kernels (pitch_angle_kernels.h) driven by the library's coefficient
// (EvaluateShape) and compares observables at common physical times with an
// independent reference:
//   * exact solutions (Eqs. 174, 175, 177, F-FOC-04 values, coth xi - 1/xi);
//   * the finite-volume reference operator (Eq. 178) propagated exactly by
//     eigen-decomposition, with its own discretization error measured by a
//     grid refinement (coarse grid bisected).
//
// Acceptance criterion (Section 26.5), fixed before any run:
//     |Y_hat - Y_ref| + z SE + eps_ref <= eps_Y,   and   z SE <= eps_Y / 4,
// with SE from independent particle batches (128 batches, independent seeded
// streams), z the two-sided normal quantile at family-wise confidence 99%
// with a Bonferroni correction over the K comparisons of the case
// (z = Phi^-1(1 - 0.01/(2K))), and eps_Y the budget of the observable:
//   Legendre moments P1..P4 and hemisphere fractions   5e-3 (absolute)
//   probability in |mu| < 0.01                          5e-4 (absolute)
//   finite-time displacement variance, focusing variance 1%  (relative)
//   focusing stationary mean                            5e-3 (absolute)
// The second inequality is the design check: if the particle count is too
// small for the budget, the case FAILS (it is never passed on a noisy run).
// Seeds and particle counts are fixed in this file.  Each particle count is
// N = (sigma z / (0.81 eps_Y / 4))^2 with sigma the largest standard
// deviation of the observable and z its case's quantile; the factor 0.81
// leaves a 3-sigma margin for the batch estimate of SE (relative standard
// deviation 1/sqrt(2*127) = 6.3% with 128 batches).
//
// Usage:  test_pitch_angle_montecarlo [--case PREFIX] [--threads N]
//                                      [--quick] [--list]
// --quick divides every particle count by 50 for development; its results
// are labelled QUICK and are not acceptance evidence.
// Output: "[ID] PASS|FAIL|SKIP message" per case (read by run_all_tests.py).
// ============================================================================

#include "pitch_angle_diffusion.h"
#include "pitch_angle_kernels.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace PA = SEP::PitchAngleDiffusion;
namespace K = SEP::PitchAngleDiffusion::Kernels;

namespace {

int gThreads = 16;
double gScale = 1.0;  // particle-count scale (1 = acceptance, 1/50 = quick)
bool gQuick = false;
constexpr int kBatches = 128;
constexpr double kPi = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Random numbers: xoshiro256** seeded by splitmix64 (Blackman & Vigna).  One
// independent stream per (case seed, batch).  Normal variates by the
// Marsaglia polar method (exact, no table).
struct Rng {
  std::uint64_t s[4];
  bool hasSpare = false;
  double spare = 0.0;
  explicit Rng(std::uint64_t seed) {
    for (int i = 0; i < 4; ++i) {
      seed += 0x9e3779b97f4a7c15ull;
      std::uint64_t z = seed;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
      s[i] = z ^ (z >> 31);
    }
  }
  static std::uint64_t Rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  std::uint64_t Next() {
    const std::uint64_t result = Rotl(s[1] * 5, 7) * 9;
    const std::uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = Rotl(s[3], 45);
    return result;
  }
  // Uniform on (0, 1]: (k + 1) 2^-53.
  double UniformOpenLeft() { return ((Next() >> 11) + 1) * 0x1.0p-53; }
  // Uniform on [0, 1).
  double Uniform() { return (Next() >> 11) * 0x1.0p-53; }
  double Normal() {
    if (hasSpare) { hasSpare = false; return spare; }
    double u, v, w;
    do {
      u = 2.0 * Uniform() - 1.0;
      v = 2.0 * Uniform() - 1.0;
      w = u * u + v * v;
    } while (w >= 1.0 || w == 0.0);
    const double factor = std::sqrt(-2.0 * std::log(w) / w);
    spare = v * factor;
    hasSpare = true;
    return u * factor;
  }
};

// Run kBatches batch functions on gThreads threads; each batch fills its own
// accumulator vector (no shared mutable state).
void RunBatches(int observables,
                const std::function<void(int, std::vector<double>&)>& batch,
                std::vector<std::vector<double>>* results) {
  results->assign(kBatches, std::vector<double>(observables, 0.0));
  std::atomic<int> next(0);
  auto worker = [&]() {
    for (int b = next++; b < kBatches; b = next++) batch(b, (*results)[b]);
  };
  std::vector<std::thread> pool;
  for (int t = 0; t < std::min(gThreads, kBatches); ++t) pool.emplace_back(worker);
  for (std::thread& t : pool) t.join();
}

struct Estimate { double mean = 0.0, se = 0.0; };
Estimate Combine(const std::vector<std::vector<double>>& batches, int index) {
  double sum = 0.0, sum2 = 0.0;
  for (const auto& b : batches) { sum += b[index]; sum2 += b[index] * b[index]; }
  Estimate e;
  e.mean = sum / kBatches;
  const double var = (sum2 - kBatches * e.mean * e.mean) / (kBatches - 1);
  e.se = std::sqrt(std::max(0.0, var) / kBatches);
  return e;
}

// Inverse standard normal CDF (Acklam's rational approximation refined by one
// Halley step with erfc; relative error ~1e-15).
double InverseNormal(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
      -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,
      2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
      -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
      -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00,
      2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
      2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) {
    const double q = std::sqrt(-2 * std::log(p));
    x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p > 1 - 0.02425) {
    const double q = std::sqrt(-2 * std::log(1 - p));
    x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
        (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
  const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
  const double u = e * std::sqrt(2 * kPi) * std::exp(x * x / 2);
  return x - u / (1 + x * u / 2);
}

struct Comparison {
  std::string name;
  Estimate value;
  double reference = 0.0;
  double referenceError = 0.0;
  double budget = 0.0;
  bool relative = false;
};

struct CaseResult {
  std::string id;
  std::string status;  // PASS / FAIL / SKIP
  std::string message;
};

// Apply the predeclared criterion to a family of comparisons.
CaseResult Judge(const std::string& id, const std::vector<Comparison>& family,
                 double seconds) {
  const int k = static_cast<int>(family.size());
  const double z = InverseNormal(1.0 - 0.01 / (2.0 * k));
  int failures = 0, designFailures = 0;
  double worstUse = 0.0;
  std::string worstName, firstFailure;
  for (const Comparison& c : family) {
    const double scale = c.relative ? std::fabs(c.reference) : 1.0;
    const double budget = c.budget * scale;
    const double error = std::fabs(c.value.mean - c.reference);
    const double used = error + z * c.value.se + c.referenceError;
    if (z * c.value.se > 0.25 * budget) {
      ++designFailures;
      if (firstFailure.empty()) {
        std::ostringstream m;
        m << c.name << ": z*SE=" << z * c.value.se << " exceeds budget/4="
          << 0.25 * budget << " (particle count too small)";
        firstFailure = m.str();
      }
    }
    if (used > budget) {
      ++failures;
      if (firstFailure.empty()) {
        std::ostringstream m;
        m.precision(6);
        m << c.name << ": estimate " << c.value.mean << " +- " << c.value.se
          << " vs reference " << c.reference << " (ref. error "
          << c.referenceError << "), |err|+zSE+eps_ref=" << used
          << " > budget " << budget;
        firstFailure = m.str();
      }
    }
    if (used / budget > worstUse) { worstUse = used / budget; worstName = c.name; }
  }
  CaseResult r;
  r.id = id;
  std::ostringstream m;
  m.precision(3);
  if (failures + designFailures == 0) {
    r.status = "PASS";
    m << k << " comparisons, z=" << z << ", worst budget use " << worstUse * 100
      << "% (" << worstName << "), " << seconds << " s";
  } else {
    r.status = "FAIL";
    m << failures << " budget and " << designFailures << " design failures of " << k
      << "; first: " << firstFailure;
  }
  if (gQuick) {
    r.status = "QUICK-" + r.status;
    m << " [quick: not acceptance evidence]";
  }
  r.message = m.str();
  return r;
}

long long Count(double nominal) {
  return std::max<long long>(kBatches, static_cast<long long>(nominal * gScale));
}

double Seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// ---------------------------------------------------------------------------
// One scattering substep of a released kernel, with the coefficient from the
// library (EvaluateShape) and the host's deterministic rate added once.
enum class KernelChoice { EulerMaruyama, Milstein, Rotation };
const char* Tag(KernelChoice k) {
  return k == KernelChoice::EulerMaruyama ? "em"
       : k == KernelChoice::Milstein ? "milstein" : "rotation";
}

struct ShapeSpec {
  PA::ModelId model = PA::ModelId::Isotropic;
  PA::ShapeParameters shape;
  double d0 = 0.5;  // amplitude [s^-1]
};

inline double Kick(KernelChoice kernel, const ShapeSpec& spec, double mu,
                   double rate, double dt, Rng& rng) {
  if (kernel == KernelChoice::Rotation) {
    // Isotropic only: nu = 2 D0.  The deterministic rate (if any) is applied
    // as an explicit Euler step before the rotation (operator splitting).
    double m = mu;
    if (rate != 0.0) m = K::FoldPitchCosine(m + rate * dt);
    const double u1 = rng.UniformOpenLeft(), u2 = rng.Uniform();
    return K::RandomRotationKick(m, 2.0 * spec.d0, dt, u1, u2).mu;
  }
  // The library's allocation-free coefficient core (the same function every
  // evaluator, including the active pointer, uses).
  double d = 0.0, dPrime = 0.0;
  PA::ShapeDmumu(mu, spec.d0, spec.model, spec.shape, &d, &dPrime);
  const double zeta = rng.Normal();
  const K::KickResult r = kernel == KernelChoice::EulerMaruyama
      ? K::EulerMaruyamaKick(mu, d, dPrime, rate, dt, zeta)
      : K::MilsteinKick(mu, d, dPrime, rate, dt, zeta);
  return r.mu;
}

std::string StepTag(double f) {
  std::ostringstream s;
  s << "f" << f;
  return s.str();
}

// ---------------------------------------------------------------------------
// AP-ISO-01: isotropic Legendre decay, Eq. (174)/(177).
CaseResult ApIso01(KernelChoice kernel, double f) {
  const auto start = std::chrono::steady_clock::now();
  const std::string id = std::string("AP-ISO-01.") + Tag(kernel) + "." + StepTag(f);
  ShapeSpec spec;  // PA-ISO, D0 = 1/2: nu = 1
  const double nu = 1.0, dt = f / nu;
  const double times[3] = {0.1, 1.0, 5.0};
  // Starts: beams at 0.5, 0, +1, -1 and the isotropic distribution (code 9).
  const double starts[5] = {0.5, 0.0, 1.0, -1.0, 9.0};
  // K = 60, z = 3.77, sigma(P1, isotropic) = 1/sqrt(3): N = 4.8e6.
  const long long perBatch = Count(4.8e6) / kBatches;
  const int nObs = 5 * 3 * 4;
  std::vector<std::vector<double>> results;
  RunBatches(nObs, [&](int b, std::vector<double>& out) {
    Rng rng(0x5150a0011ull * 1000003ull + static_cast<std::uint64_t>(b) +
            static_cast<std::uint64_t>(f * 1e7) * 7919ull +
            static_cast<std::uint64_t>(kernel) * 104729ull);
    std::vector<double> mu(perBatch);
    for (int s = 0; s < 5; ++s) {
      for (long long i = 0; i < perBatch; ++i)
        mu[i] = starts[s] > 2.0 ? 2.0 * rng.Uniform() - 1.0 : starts[s];
      double t = 0.0;
      for (int ti = 0; ti < 3; ++ti) {
        // Advance to the common output time with the fixed step (the last
        // substep is shortened to land exactly on the output time).
        while (t < times[ti] - 1e-12) {
          const double h = std::min(dt, times[ti] - t);
          for (long long i = 0; i < perBatch; ++i)
            mu[i] = Kick(kernel, spec, mu[i], 0.0, h, rng);
          t += h;
        }
        for (int l = 1; l <= 4; ++l) {
          double sum = 0.0;
          for (long long i = 0; i < perBatch; ++i) sum += K::Legendre(l, mu[i]);
          out[(s * 3 + ti) * 4 + (l - 1)] = sum / perBatch;
        }
      }
    }
  }, &results);
  std::vector<Comparison> family;
  for (int s = 0; s < 5; ++s)
    for (int ti = 0; ti < 3; ++ti)
      for (int l = 1; l <= 4; ++l) {
        Comparison c;
        std::ostringstream name;
        name << "P" << l << " start=" << (starts[s] > 2 ? std::string("iso")
                                                        : std::to_string(starts[s]))
             << " nu*t=" << times[ti];
        c.name = name.str();
        c.value = Combine(results, (s * 3 + ti) * 4 + (l - 1));
        const double initial = starts[s] > 2.0 ? 0.0 : K::Legendre(l, starts[s]);
        c.reference = initial * std::exp(-0.5 * l * (l + 1) * nu * times[ti]);
        c.budget = 5e-3;
        family.push_back(c);
      }
  return Judge(id, family, Seconds(start));
}

// AP-ISO-02: finite-time displacement from an isotropic stationary start,
// Eq. (175), with streaming along the subcycled trajectory (trapezoid
// ds = v (mu_n + mu_{n+1})/2 dt, Section 16.12) and the ballistic support
// |dz| <= v t.
CaseResult ApIso02(KernelChoice kernel, double f, int nTimes) {
  const auto start = std::chrono::steady_clock::now();
  std::ostringstream idStream;
  idStream << "AP-ISO-02." << Tag(kernel) << "." << StepTag(f) << ".t"
           << (nTimes == 4 ? 50 : 10);
  const std::string id = idStream.str();
  ShapeSpec spec;
  const double nu = 1.0, v = 1.0, dt = f / nu;
  const double times[4] = {0.1, 1.0, 10.0, 50.0};
  // Relative SE of a variance ~ sqrt(2/N) (Gaussian-like displacement);
  // K <= 4, z = 3.02: N = 2 (3.02/(0.81*0.0025))^2 = 4.5e6.
  const long long perBatch = Count(4.5e6) / kBatches;
  std::vector<std::vector<double>> results;
  // Observables per time: mean dz, variance/(2t), support violation flag.
  RunBatches(12, [&](int b, std::vector<double>& out) {
    Rng rng(0x150a0022ull * 1000003ull + static_cast<std::uint64_t>(b) +
            static_cast<std::uint64_t>(f * 1e7) * 7919ull +
            static_cast<std::uint64_t>(kernel) * 104729ull);
    std::vector<double> mu(perBatch), z(perBatch, 0.0);
    for (long long i = 0; i < perBatch; ++i) mu[i] = 2.0 * rng.Uniform() - 1.0;
    double t = 0.0;
    for (int ti = 0; ti < nTimes; ++ti) {
      while (t < times[ti] - 1e-12) {
        const double h = std::min(dt, times[ti] - t);
        for (long long i = 0; i < perBatch; ++i) {
          const double next = Kick(kernel, spec, mu[i], 0.0, h, rng);
          z[i] += v * 0.5 * (mu[i] + next) * h;
          mu[i] = next;
        }
        t += h;
      }
      double sum = 0.0, sum2 = 0.0, violations = 0.0;
      for (long long i = 0; i < perBatch; ++i) {
        sum += z[i];
        sum2 += z[i] * z[i];
        if (std::fabs(z[i]) > v * times[ti] * (1.0 + 1e-12)) violations += 1.0;
      }
      out[ti * 3 + 0] = sum / perBatch;
      out[ti * 3 + 1] = (sum2 / perBatch) / (2.0 * times[ti]);
      out[ti * 3 + 2] = violations;
    }
  }, &results);
  std::vector<Comparison> family;
  for (int ti = 0; ti < nTimes; ++ti) {
    const double t = times[ti];
    Comparison c;
    c.name = "kappa(t) nu*t=" + std::to_string(t);
    c.value = Combine(results, ti * 3 + 1);
    // Eq. (175): kappa(t) = v^2/(3 nu) [1 - (1 - e^{-nu t})/(nu t)].
    c.reference = v * v / (3.0 * nu) * (1.0 - (1.0 - std::exp(-nu * t)) / (nu * t));
    c.budget = 0.01;
    c.relative = true;
    family.push_back(c);
    double violations = 0.0;
    for (const auto& b : results) violations += b[ti * 3 + 2];
    if (violations > 0.0) {
      CaseResult r{id, "FAIL", "ballistic support |dz| <= v t violated"};
      return r;
    }
  }
  return Judge(id, family, Seconds(start));
}

// AP-TIME-01: a rate update inside a host interval (nu = 0.4 for 0.3 s then
// 2.1 for 0.7 s); optical time 1.59 and Eq. (177).
CaseResult ApTime01(KernelChoice kernel) {
  const auto start = std::chrono::steady_clock::now();
  const std::string id = std::string("AP-TIME-01.") + Tag(kernel);
  const double f = 0.001, mu0 = 0.5;
  // K = 4, z = 3.02, sigma <= 1/sqrt(3): N = 3.1e6 (rounded up).
  const long long perBatch = Count(3.2e6) / kBatches;
  std::vector<std::vector<double>> results;
  RunBatches(4, [&](int b, std::vector<double>& out) {
    Rng rng(0x7130001ull * 1000003ull + static_cast<std::uint64_t>(b) +
            static_cast<std::uint64_t>(kernel) * 104729ull);
    std::vector<double> mu(perBatch, mu0);
    // The host splits its interval at the published update time (Section
    // 16.12, "the host subdivides an interval at a published update").
    const double segments[2][2] = {{0.4, 0.3}, {2.1, 0.7}};
    for (const auto& segment : segments) {
      ShapeSpec spec;
      spec.d0 = 0.5 * segment[0];
      const double dt = f / segment[0];
      double t = 0.0;
      while (t < segment[1] - 1e-12) {
        const double h = std::min(dt, segment[1] - t);
        for (long long i = 0; i < perBatch; ++i)
          mu[i] = Kick(kernel, spec, mu[i], 0.0, h, rng);
        t += h;
      }
    }
    for (int l = 1; l <= 4; ++l) {
      double sum = 0.0;
      for (long long i = 0; i < perBatch; ++i) sum += K::Legendre(l, mu[i]);
      out[l - 1] = sum / perBatch;
    }
  }, &results);
  std::vector<Comparison> family;
  const double optical = 0.4 * 0.3 + 2.1 * 0.7;  // 1.59
  for (int l = 1; l <= 4; ++l) {
    Comparison c;
    c.name = "P" + std::to_string(l);
    c.value = Combine(results, l - 1);
    c.reference = K::Legendre(l, mu0) * std::exp(-0.5 * l * (l + 1) * optical);
    c.budget = 5e-3;
    family.push_back(c);
  }
  return Judge(id, family, Seconds(start));
}

// ---------------------------------------------------------------------------
// Finite-volume reference grids for PA-EPS: edges forced at -0.01, 0, 0.01;
// cells clustered toward mu = +-1 and toward the occupation band.  The fine
// grid bisects every coarse cell, so a coarse-cell initial condition is
// represented exactly on both.
std::vector<double> EpsCoarseEdges() {
  const int inner = 20, outer = 380;  // per hemisphere: n = 800
  std::vector<double> positive;
  for (int i = 0; i <= inner; ++i) positive.push_back(0.01 * i / inner);
  for (int i = 1; i <= outer; ++i) {
    const double s = static_cast<double>(i) / outer;  // (0,1]
    // Clustered near both ends of [0.01, 1].
    const double w = 0.5 - 0.5 * std::cos(kPi * s);
    positive.push_back(0.01 + 0.99 * w);
  }
  positive.back() = 1.0;
  std::vector<double> edges;
  for (std::size_t i = positive.size(); i-- > 1;) edges.push_back(-positive[i]);
  for (double p : positive) edges.push_back(p);
  edges.front() = -1.0;
  return edges;
}

std::vector<double> Bisect(const std::vector<double>& coarse) {
  std::vector<double> fine;
  for (std::size_t i = 0; i + 1 < coarse.size(); ++i) {
    fine.push_back(coarse[i]);
    fine.push_back(0.5 * (coarse[i] + coarse[i + 1]));
  }
  fine.push_back(coarse.back());
  return fine;
}

struct FvObservables { double legendre[4], occupation, negative; };

FvObservables FvMeasure(const std::vector<double>& edges,
                        const std::vector<double>& mass) {
  FvObservables o{};
  for (std::size_t i = 0; i < mass.size(); ++i) {
    for (int l = 1; l <= 4; ++l)
      o.legendre[l - 1] += mass[i] * K::LegendreCellAverage(l, edges[i], edges[i + 1]);
    if (edges[i] >= -0.01 - 1e-15 && edges[i + 1] <= 0.01 + 1e-15) o.occupation += mass[i];
    if (edges[i + 1] <= 1e-15) o.negative += mass[i];
  }
  return o;
}

// Initial mass on `edges` of a uniform distribution on [a, b] (a, b edges).
std::vector<double> InitialMass(const std::vector<double>& edges, double a, double b) {
  std::vector<double> m(edges.size() - 1, 0.0);
  for (std::size_t i = 0; i + 1 < edges.size(); ++i)
    if (edges[i] >= a - 1e-15 && edges[i + 1] <= b + 1e-15)
      m[i] = (edges[i + 1] - edges[i]) / (b - a);
  return m;
}

// AP-ZERO-02 / PS05 for PA-EPS with eps > 0.
CaseResult ApZero02(double epsilon, KernelChoice kernel, double f) {
  const auto start = std::chrono::steady_clock::now();
  std::ostringstream idStream;
  idStream << "AP-ZERO-02.eps" << epsilon << "." << Tag(kernel) << "." << StepTag(f);
  const std::string id = idStream.str();
  ShapeSpec spec;
  spec.model = PA::ModelId::EpsilonForm;
  spec.shape.epsilon = epsilon;
  spec.d0 = 0.5;  // nu0 = 1
  const double nuMax = 2.0 * spec.d0 * (0.5 + epsilon);
  const double dt = f / nuMax;
  const double times[2] = {0.2, 2.0};
  // Starts (uniform on a coarse-grid interval): uniform on [-1,1]; zero-
  // centred (the two coarse cells adjacent to 0); near zero (the coarse cell
  // containing 0.005); beam (cell containing 0.95); endpoint beam (last cell).
  const std::vector<double> coarse = EpsCoarseEdges();
  const std::vector<double> fine = Bisect(coarse);
  auto cellOf = [&](double mu) {
    std::size_t i = 0;
    while (i + 2 < coarse.size() && !(mu >= coarse[i] && mu < coarse[i + 1])) ++i;
    return i;
  };
  const std::size_t zeroEdge = static_cast<std::size_t>(
      std::find(coarse.begin(), coarse.end(), 0.0) - coarse.begin());
  const std::size_t near = cellOf(0.005), beam = cellOf(0.95);
  struct Start { const char* name; double a, b; bool heavy; };
  const Start starts[5] = {
      {"uniform", -1.0, 1.0, true},
      {"zero", coarse[zeroEdge - 1], coarse[zeroEdge + 1], true},
      {"near-zero", coarse[near], coarse[near + 1], true},
      {"beam-0.95", coarse[beam], coarse[beam + 1], false},
      {"endpoint", coarse[coarse.size() - 2], 1.0, false}};
  // Reference: exact propagation on the coarse and the fine grid; the fine
  // value is the reference and |fine - coarse| bounds its error.  The
  // eigen-decompositions depend only on eps and are cached across the kernel
  // and step cases of the same eps.
  struct CachedReference { double epsilon; K::FvPropagator coarse, fine; };
  static std::vector<CachedReference> cache;
  const CachedReference* cached = nullptr;
  for (const CachedReference& c : cache)
    if (c.epsilon == epsilon) cached = &c;
  if (cached == nullptr) {
    CachedReference entry;
    entry.epsilon = epsilon;
    K::FvGenerator gc, gf;
    if (!K::BuildFvGenerator(spec.model, spec.shape, spec.d0, coarse, &gc).ok() ||
        !K::BuildFvGenerator(spec.model, spec.shape, spec.d0, fine, &gf).ok() ||
        !K::PrepareFvPropagator(gc, &entry.coarse).ok() ||
        !K::PrepareFvPropagator(gf, &entry.fine).ok())
      return CaseResult{id, "FAIL", "finite-volume reference could not be built"};
    cache.push_back(entry);
    cached = &cache.back();
  }
  const K::FvPropagator& pc = cached->coarse;
  const K::FvPropagator& pf = cached->fine;
  FvObservables refFine[5][2], refCoarse[5][2];
  for (int s = 0; s < 5; ++s)
    for (int ti = 0; ti < 2; ++ti) {
      std::vector<double> m;
      K::PropagateFv(pf, InitialMass(fine, starts[s].a, starts[s].b), times[ti], &m);
      refFine[s][ti] = FvMeasure(fine, m);
      K::PropagateFv(pc, InitialMass(coarse, starts[s].a, starts[s].b), times[ti], &m);
      refCoarse[s][ti] = FvMeasure(coarse, m);
    }
  // K = 60, z = 3.77.  P(|mu|<0.01) ~ 0.01 (sigma 0.0995) needs
  // N = (0.0995*3.77/(0.81*1.25e-4))^2 = 1.37e7 for the starts that populate
  // the band; Legendre/hemisphere (sigma <= 0.58) need 4.6e6.
  // Known design defect (first campaign, kept frozen and reported): for the
  // concentrated "zero" and "near-zero" starts the band probability is far
  // from 0.01 at nu0 t = 0.2 (sigma up to ~0.5), so N = 1.37e7 does not meet
  // z SE <= budget/4 there and those comparisons FAIL as design failures.
  // The PA-EPS verdict does not rest on them: eps = 0.01 and 0.048 also have
  // budget failures of the uniform start, whose sizing is correct.
  const long long heavy = Count(1.37e7) / kBatches, light = Count(4.6e6) / kBatches;
  const int perStart = 2 * 6;
  std::vector<std::vector<double>> results;
  RunBatches(5 * perStart, [&](int b, std::vector<double>& out) {
    Rng rng(0x2e400002ull * 1000003ull + static_cast<std::uint64_t>(b) +
            static_cast<std::uint64_t>(f * 1e7) * 7919ull +
            static_cast<std::uint64_t>(epsilon * 1e6) * 15485863ull +
            static_cast<std::uint64_t>(kernel) * 104729ull);
    for (int s = 0; s < 5; ++s) {
      const long long n = starts[s].heavy ? heavy : light;
      std::vector<double> mu(n);
      for (long long i = 0; i < n; ++i)
        mu[i] = starts[s].a + (starts[s].b - starts[s].a) * rng.Uniform();
      double t = 0.0;
      for (int ti = 0; ti < 2; ++ti) {
        while (t < times[ti] - 1e-12) {
          const double h = std::min(dt, times[ti] - t);
          for (long long i = 0; i < n; ++i) mu[i] = Kick(kernel, spec, mu[i], 0.0, h, rng);
          t += h;
        }
        double sums[6] = {0, 0, 0, 0, 0, 0};
        for (long long i = 0; i < n; ++i) {
          for (int l = 1; l <= 4; ++l) sums[l - 1] += K::Legendre(l, mu[i]);
          if (std::fabs(mu[i]) < 0.01) sums[4] += 1.0;
          if (mu[i] < 0.0) sums[5] += 1.0;
        }
        for (int k = 0; k < 6; ++k) out[s * perStart + ti * 6 + k] = sums[k] / n;
      }
    }
  }, &results);
  std::vector<Comparison> family;
  for (int s = 0; s < 5; ++s)
    for (int ti = 0; ti < 2; ++ti) {
      const FvObservables& rf = refFine[s][ti];
      const FvObservables& rc = refCoarse[s][ti];
      for (int k = 0; k < 6; ++k) {
        Comparison c;
        std::ostringstream name;
        const char* what[6] = {"P1", "P2", "P3", "P4", "P(|mu|<0.01)", "P(mu<0)"};
        name << what[k] << " " << starts[s].name << " nu0*t=" << times[ti];
        c.name = name.str();
        c.value = Combine(results, s * perStart + ti * 6 + k);
        const double fineValue = k < 4 ? rf.legendre[k] : (k == 4 ? rf.occupation : rf.negative);
        const double coarseValue = k < 4 ? rc.legendre[k] : (k == 4 ? rc.occupation : rc.negative);
        c.reference = fineValue;
        c.referenceError = std::fabs(fineValue - coarseValue);
        c.budget = k == 4 ? 5e-4 : 5e-3;
        family.push_back(c);
      }
    }
  return Judge(id, family, Seconds(start));
}

// AP-ZERO-02 pure-scattering disconnection (eps = 0): the reference operator
// is disconnected (Eq. 178, Section 16.5), so the hemisphere-crossing
// fraction from a positive start must vanish; the discrete kernels cross
// only at O(dt) rate.  The case checks that the measured crossing fraction
// decreases with every step refinement and that the finite-volume reference
// is exactly disconnected; it is a diagnostic of the barrier, not a release
// test (eps = 0 is not normalizable to a lambda law).
CaseResult ApZero02Disconnected(KernelChoice kernel) {
  const auto start = std::chrono::steady_clock::now();
  const std::string id = std::string("AP-ZERO-02.eps0.") + Tag(kernel);
  ShapeSpec spec;
  spec.model = PA::ModelId::EpsilonForm;
  spec.shape.epsilon = 0.0;
  spec.d0 = 0.5;
  K::FvGenerator g;
  K::BuildFvGenerator(spec.model, spec.shape, spec.d0, EpsCoarseEdges(), &g);
  if (g.internalExtension.find("disconnected") == std::string::npos)
    return CaseResult{id, "FAIL", "reference operator is not disconnected"};
  const double steps[3] = {0.01, 0.003, 0.001};
  double crossed[3];
  const long long perBatch = Count(3.2e6) / kBatches;
  for (int k = 0; k < 3; ++k) {
    std::vector<std::vector<double>> results;
    const double dt = steps[k] / (2.0 * spec.d0 * 0.5);
    RunBatches(1, [&](int b, std::vector<double>& out) {
      Rng rng(0x0d15c0ull * 1000003ull + static_cast<std::uint64_t>(b) + k * 7919ull +
              static_cast<std::uint64_t>(kernel) * 104729ull);
      std::vector<double> mu(perBatch);
      for (long long i = 0; i < perBatch; ++i) mu[i] = rng.Uniform();  // (0,1)
      double t = 0.0;
      while (t < 1.0 - 1e-12) {
        const double h = std::min(dt, 1.0 - t);
        for (long long i = 0; i < perBatch; ++i) mu[i] = Kick(kernel, spec, mu[i], 0.0, h, rng);
        t += h;
      }
      double sum = 0.0;
      for (long long i = 0; i < perBatch; ++i) sum += mu[i] < 0.0 ? 1.0 : 0.0;
      out[0] = sum / perBatch;
    }, &results);
    crossed[k] = Combine(results, 0).mean;
  }
  std::ostringstream m;
  m.precision(3);
  m << "crossing fraction at nu_max*dt = 0.01/0.003/0.001: " << crossed[0] << "/"
    << crossed[1] << "/" << crossed[2] << " (reference: 0, disconnected), "
    << Seconds(start) << " s";
  // Non-increasing under refinement (an exactly zero sequence, which the
  // Milstein correction produces for D ~ |mu| near 0, is the exact result).
  const bool decreasing = crossed[1] <= crossed[0] && crossed[2] <= crossed[1];
  CaseResult r{id, decreasing ? "PASS" : "FAIL", m.str()};
  if (gQuick) r.status = "QUICK-" + r.status;
  return r;
}

// AP-FOC-01: constant focusing with isotropic scattering (v = nu = 1,
// D0 = 1/2, 1/L = xi).  Forward Ito drift dmu/dt = v (1 - mu^2)/(2L) added
// once by this composer (Section 16.12 static-tube reference).  Observables:
// stationary <mu> = coth xi - 1/xi (Eq. 176) and the long-time centred
// variance slope kappa_var (Eq. 115; F-FOC-04, or the finite-volume Poisson
// reference of PAD12 for xi = 0.1).  The slope (Var z(t2) - Var z(t1)) /
// (2 (t2 - t1)) removes the constant transient offset.
CaseResult ApFoc01(KernelChoice kernel, double xi, double f) {
  const auto start = std::chrono::steady_clock::now();
  std::ostringstream idStream;
  idStream << "AP-FOC-01." << Tag(kernel) << ".xi" << xi << "." << StepTag(f);
  const std::string id = idStream.str();
  ShapeSpec spec;  // isotropic, nu = 1
  const double v = 1.0, inverseL = xi, dt = f;
  // Burn-in 5/nu (uniform start relaxes to rho_+), then the variance slope
  // between 5 and 20 after the burn-in.  Relative SE of the slope
  // ~1.82/sqrt(N) (sampling of Var Z and of Cov(z1, Z)); K = 2, z = 2.81:
  // N = (1.82*2.81/(0.81*0.0025))^2 = 6.4e6.
  // For |xi| = 5 the measured relative noise of the slope is 1.33 times
  // larger (pilot run, quick mode): N = 1.33^2 * 6.4e6 = 1.13e7.
  const double burn = 5.0, t1 = 5.0, t2 = 20.0;
  const long long perBatch =
      Count(std::fabs(xi) >= 5.0 ? 1.13e7 : 6.4e6) / kBatches;
  std::vector<std::vector<double>> results;
  RunBatches(2, [&](int b, std::vector<double>& out) {
    Rng rng(0xf0c0001ull * 1000003ull + static_cast<std::uint64_t>(b) +
            static_cast<std::uint64_t>((xi + 10.0) * 1e4) * 7919ull +
            static_cast<std::uint64_t>(f * 1e7) * 15485863ull +
            static_cast<std::uint64_t>(kernel) * 104729ull);
    std::vector<double> mu(perBatch), z(perBatch, 0.0);
    for (long long i = 0; i < perBatch; ++i) mu[i] = 2.0 * rng.Uniform() - 1.0;
    double t = 0.0, meanAccumulator = 0.0, var1 = 0.0;
    long long meanSamples = 0;
    auto advance = [&](double until, bool stream, bool average) {
      while (t < until - 1e-12) {
        const double h = std::min(dt, until - t);
        for (long long i = 0; i < perBatch; ++i) {
          const double rate = v * (1.0 - mu[i] * mu[i]) * 0.5 * inverseL;
          const double next = Kick(kernel, spec, mu[i], rate, h, rng);
          if (stream) z[i] += v * 0.5 * (mu[i] + next) * h;
          mu[i] = next;
        }
        t += h;
        if (average) {
          double s = 0.0;
          for (long long i = 0; i < perBatch; ++i) s += mu[i];
          meanAccumulator += s / perBatch;
          ++meanSamples;
        }
      }
    };
    advance(burn, false, false);
    const double origin = t;
    advance(origin + t1, true, true);
    double s = 0.0, s2 = 0.0;
    for (long long i = 0; i < perBatch; ++i) { s += z[i]; s2 += z[i] * z[i]; }
    var1 = s2 / perBatch - (s / perBatch) * (s / perBatch);
    advance(origin + t2, true, true);
    s = s2 = 0.0;
    for (long long i = 0; i < perBatch; ++i) { s += z[i]; s2 += z[i] * z[i]; }
    const double var2 = s2 / perBatch - (s / perBatch) * (s / perBatch);
    out[0] = meanAccumulator / meanSamples;
    out[1] = (var2 - var1) / (2.0 * (t2 - t1));
  }, &results);
  // References.
  double meanReference, kappaReference, kappaError = 0.0;
  if (xi == 0.0) {
    meanReference = 0.0;
    kappaReference = 1.0 / 3.0;  // Eq. (175) limit v^2/(3 nu)
  } else {
    const double a = std::fabs(xi);
    // coth a - 1/a, series for small a (Section 26.3 PS08).
    double m = a < 1e-3 ? a / 3.0 - a * a * a / 45.0 : 1.0 / std::tanh(a) - 1.0 / a;
    meanReference = xi > 0 ? m : -m;
    if (a == 1.5) kappaReference = 0.178275780337345;        // F-FOC-04
    else if (a == 5.0) kappaReference = 0.0103248294848381;  // F-FOC-04
    else if (a == 1.0) kappaReference = 0.248076740609470;   // F-FOC-04
    else {
      // Finite-volume Poisson reference (validated against F-FOC-04 in
      // PAD12) on two grids; the difference bounds its error.
      double values[2];
      for (int g = 0; g < 2; ++g) {
        const int n = g == 0 ? 800 : 1600;
        std::vector<double> edges(n + 1);
        for (int i = 0; i <= n; ++i) edges[i] = -std::cos(kPi * i / n);
        edges.front() = -1.0; edges.back() = 1.0;
        const double shift = xi > 0 ? 1.0 : -1.0;
        std::vector<double> masses(n), conductance(n - 1), meanMu(n);
        auto rho = [&](double x) { return std::exp(xi * (x - shift)); };
        for (int i = 0; i < n; ++i) {
          // Closed-form cell integrals of exp(xi (mu - shift)) and mu times it.
          const double ea = rho(edges[i]), eb = rho(edges[i + 1]);
          masses[i] = (eb - ea) / xi;
          meanMu[i] = ((edges[i + 1] * eb - edges[i] * ea) / xi - (eb - ea) / (xi * xi)) /
                      masses[i];
        }
        for (int i = 0; i + 1 < n; ++i) {
          // T = [Int dmu / (D rho)]^-1 by 64-point composite midpoint rule
          // (smooth integrand away from the endpoints).
          const double a0 = 0.5 * (edges[i] + edges[i + 1]);
          const double b0 = 0.5 * (edges[i + 1] + edges[i + 2]);
          double integral = 0.0;
          const int m2 = 64;
          for (int j = 0; j < m2; ++j) {
            const double x = a0 + (b0 - a0) * (j + 0.5) / m2;
            integral += 1.0 / (0.5 * (1.0 - x * x) * rho(x));
          }
          conductance[i] = 1.0 / (integral * (b0 - a0) / m2);
        }
        K::FvGenerator foc;
        K::BuildFvGeneratorFromConductances(edges, masses, conductance, &foc);
        K::FvAsymptoticVariance(foc, meanMu, &values[g]);
      }
      kappaReference = values[1];
      kappaError = std::fabs(values[1] - values[0]);
    }
  }
  std::vector<Comparison> family(2);
  family[0].name = "stationary <mu>";
  family[0].value = Combine(results, 0);
  family[0].reference = meanReference;
  family[0].budget = 5e-3;
  family[1].name = "kappa_var";
  family[1].value = Combine(results, 1);
  family[1].reference = kappaReference;
  family[1].referenceError = kappaError;
  family[1].budget = 0.01;
  family[1].relative = true;
  return Judge(id, family, Seconds(start));
}

CaseResult Skip(const std::string& id, const std::string& reason) {
  return CaseResult{id, "SKIP", reason};
}

}  // namespace

int main(int argc, char** argv) {
  std::string prefix;
  bool list = false;
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--case") && i + 1 < argc) prefix = argv[++i];
    else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) gThreads = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--quick")) { gQuick = true; gScale = 1.0 / 50.0; }
    else if (!std::strcmp(argv[i], "--list")) list = true;
    else { std::fprintf(stderr, "unknown argument %s\n", argv[i]); return 2; }
  }
  gThreads = std::max(1, gThreads);
  using Factory = std::function<CaseResult()>;
  std::vector<std::pair<std::string, Factory>> cases;
  const KernelChoice ito[2] = {KernelChoice::EulerMaruyama, KernelChoice::Milstein};
  const KernelChoice all[3] = {KernelChoice::EulerMaruyama, KernelChoice::Milstein,
                               KernelChoice::Rotation};
  // Three step levels for every kernel (Section 26.5: at least three levels).
  for (KernelChoice k : all)
    for (double f : {0.01, 0.003, 0.001})
      cases.push_back({std::string("AP-ISO-01.") + Tag(k) + "." + StepTag(f),
                       [k, f]() { return ApIso01(k, f); }});
  // Campaign scope after the pilot run (IMPLEMENTATION_STATUS.md):
  // Euler-Maruyama with mirror folding failed AP-ISO-01 at every tested step
  // (endpoint-beam bias, Section 16.4), so it cannot be certified; it is kept
  // only where it documents that outcome (AP-ISO-01, AP-TIME-01, the eps = 0
  // barrier).  Milstein is the certification candidate at nu_max dt = 0.003.
  //
  // Finite-time streaming: Milstein at the candidate step to nu t = 50 and at
  // f = 0.001 to nu t = 10 (convergence); rotation at 0.003 to nu t = 10.
  cases.push_back({"AP-ISO-02.milstein.f0.003.t50",
                   []() { return ApIso02(KernelChoice::Milstein, 0.003, 4); }});
  cases.push_back({"AP-ISO-02.milstein.f0.001.t10",
                   []() { return ApIso02(KernelChoice::Milstein, 0.001, 3); }});
  cases.push_back({"AP-ISO-02.rotation.f0.003.t10",
                   []() { return ApIso02(KernelChoice::Rotation, 0.003, 3); }});
  // Second campaign (scope frozen after the first full campaign, before any
  // of these cases ran): AP-FOC-01.milstein.xi5.f0.003 exceeded its budget
  // (101%), so f = 0.003 cannot be certified for Milstein.  The roadmap's next
  // refinement level f = 0.001 is therefore taken through every Milstein case
  // it lacked: streaming to nu t = 50 here and AP-FOC-01 below.  No budget,
  // seed or count of an existing case was changed; the failed cases stay in
  // the record.
  cases.push_back({"AP-ISO-02.milstein.f0.001.t50",
                   []() { return ApIso02(KernelChoice::Milstein, 0.001, 4); }});
  for (KernelChoice k : all)
    cases.push_back({std::string("AP-TIME-01.") + Tag(k), [k]() { return ApTime01(k); }});
  for (double eps : {0.01, 0.048, 1.0})
    for (KernelChoice k : {KernelChoice::Milstein})
      for (double f : {0.01, 0.003, 0.001}) {
        std::ostringstream id;
        id << "AP-ZERO-02.eps" << eps << "." << Tag(k) << "." << StepTag(f);
        cases.push_back({id.str(), [eps, k, f]() { return ApZero02(eps, k, f); }});
      }
  for (KernelChoice k : ito)
    cases.push_back({std::string("AP-ZERO-02.eps0.") + Tag(k),
                     [k]() { return ApZero02Disconnected(k); }});
  for (double xi : {0.0, 0.1, -0.1, 1.5, -1.5, 5.0, -5.0}) {
    std::ostringstream id;
    id << "AP-FOC-01.milstein.xi" << xi << "." << StepTag(0.003);
    cases.push_back({id.str(), [xi]() {
                       return ApFoc01(KernelChoice::Milstein, xi, 0.003); }});
  }
  // Second campaign: the same focusing cases at the refined step f = 0.001.
  for (double xi : {0.0, 0.1, -0.1, 1.5, -1.5, 5.0, -5.0}) {
    std::ostringstream id;
    id << "AP-FOC-01.milstein.xi" << xi << "." << StepTag(0.001);
    cases.push_back({id.str(), [xi]() {
                       return ApFoc01(KernelChoice::Milstein, xi, 0.001); }});
  }
  // Acceptance cases of Section 26.5 that this release does not exercise,
  // with the reason (never counted as passed).
  const std::pair<const char*, const char*> skipped[] = {
      {"AP-ZERO-01", "PA-QFORM is refused with Ito kernels (user decision "
                     "2026-10-13); no q-form pair is released"},
      {"AP-RNG-01", "no adaptive rejection: step rules select dt before the "
                    "noise is drawn, so no Brownian bridge (Eq. N4) is used"},
      {"AP-CLOCK-01", "no clock factor is applied (user decision: host time, "
                      "order 0 in U/c)"},
      {"AP-FOC-02", "no backward pusher in Release A"},
      {"AP-TENSOR-01", "no (mu,p) tensor capability in Release A"},
      {"AP-FRAME-01", "no frame transformation in the library (plasma-frame "
                      "p, mu kept by the host)"},
      {"AP-SPEC-01", "spectral models are Release B"},
      {"AP-WAVE-01", "wave populations are Release C"},
      {"AP-ENERGY-01", "wave-frame energy exchange is Release C"},
      {"AP-ENERGY-02", "particle-wave coupling is Release C"},
      {"AP-GEOM-01..04", "host adapter lifecycle: tested in srcSEP3D, not in the "
                         "standalone library"}};
  for (const auto& s : skipped) {
    const std::string id = s.first, reason = s.second;
    cases.push_back({id, [id, reason]() { return Skip(id, reason); }});
  }
  if (list) {
    for (const auto& c : cases) std::printf("%s\n", c.first.c_str());
    return 0;
  }
  int failures = 0, passes = 0, skips = 0;
  for (const auto& c : cases) {
    if (!prefix.empty() && c.first.rfind(prefix, 0) != 0) continue;
    const CaseResult r = c.second();
    std::printf("[%s] %s %s\n", r.id.c_str(), r.status.c_str(), r.message.c_str());
    std::fflush(stdout);
    if (r.status.find("FAIL") != std::string::npos) ++failures;
    else if (r.status == "SKIP") ++skips;
    else ++passes;
  }
  std::printf("pitch-angle acceptance summary: %d passed, %d failed, %d skipped%s\n",
              passes, failures, skips, gQuick ? " (quick, not acceptance)" : "");
  return failures ? 1 : 0;
}
