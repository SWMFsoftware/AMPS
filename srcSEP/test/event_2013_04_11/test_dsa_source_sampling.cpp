// SRC01-SRC03: importance sampling of the field-line DSA shock source.
//
// Production path under test: srcSEP/field_line.cpp,
// GetMomentum_Tenishev2005AIAA, branch LogUniformMomentumSampling (CLI
// --shock-energy-sampling log-uniform).  That branch draws the dimensionless
// momentum x = p/pmin log-uniformly on [1, X], X = pmax/pmin, and gives each
// macroparticle the weight w(x) = SEP::Injection::
// ImportanceWeightForLogUniformProposal(target, x) for the target number
// density dN/dp ~ p^(2-q) (q = 3 s / (s - 1), the DSA phase-space index).
//
// Independent references (closed forms, not calls of the helper's own
// normalization):
//   * unbiasedness: Int_0^1 w(X^u) du = 1 exactly, for every q;
//   * the weighted proposal reproduces the target CDF
//       F(x_c) = (x_c^(3-q) - 1) / (X^(3-q) - 1)        (q != 3);
//   * every weight stays finite for weak shocks (q up to 63, s = 1.05),
//     which is why the production branch works in x rather than SI momentum.
// The integrals use deterministic composite Simpson quadrature in u, so the
// test has no stochastic tolerance.  A regression of the proposal density
// (e.g. log-uniform in E instead of p without the matching Jacobian), of the
// target index (q instead of q-2), or of the scaling would fail SRC01/SRC02.
#include "sep_injection_spectrum.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void Require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

// Proton momentum [kg m/s] for kinetic energy T [J] (exact relativistic).
double ProtonMomentum(double kineticJ) {
  const double c = 299792458.0;
  const double restJ = 1.67262192595e-27 * c * c;
  return std::sqrt(kineticJ * kineticJ + 2.0 * kineticJ * restJ) / c;
}

double Weight(const SEP::Injection::Spectrum& scaled, double x) {
  const SEP::Transport::ScalarResult w =
      SEP::Injection::ImportanceWeightForLogUniformProposal(scaled, x);
  Require(w.status.ok() && std::isfinite(w.value) && w.value >= 0.0,
          "weight is not finite/non-negative");
  return w.value;
}

// Composite Simpson rule of f(u) on [a,b] with an even number n of panels.
template <class F>
double Simpson(F f, double a, double b, int n) {
  const double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; ++i) sum += f(a + i * h) * ((i % 2) ? 4.0 : 2.0);
  return sum * h / 3.0;
}

}  // namespace

int main() {
  const double eV = 1.602176634e-19;
  // The 2013-04-11 event source bounds: 10 keV - 200 MeV protons.
  const double pmin = ProtonMomentum(10.0e3 * eV);
  const double pmax = ProtonMomentum(200.0e6 * eV);
  const double X = pmax / pmin;
  const double logX = std::log(X);
  const double qs[] = {4.0, 4.5, 6.0, 30.0, 63.0};

  for (double q : qs) {
    SEP::Injection::Spectrum scaled;
    scaled.measure = SEP::Injection::Measure::Momentum;
    scaled.minimum = 1.0;
    scaled.maximum = X;
    scaled.powerIndex = q - 2.0;  // dN/dp ~ p^(2-q)
    // The integrand is sharply peaked at u=0 for weak shocks (width ~
    // 1/((q-3) ln X)); 400000 panels resolve it with large margin.
    const int panels = 400000;
    const double mean = Simpson(
        [&](double u) { return Weight(scaled, std::exp(u * logX)); },
        0.0, 1.0, panels);
    Require(std::fabs(mean - 1.0) < 1.0e-9,
            "SRC01 mean importance weight differs from 1 for q=" +
                std::to_string(q));

    // SRC02: weighted proposal CDF at three cut points versus closed form.
    for (double fraction : {0.05, 0.3, 0.7}) {
      const double xc = std::exp(fraction * logX);
      const double weighted = Simpson(
          [&](double u) { return Weight(scaled, std::exp(u * logX)); },
          0.0, fraction, panels);
      const double expected =
          (std::pow(xc, 3.0 - q) - 1.0) / (std::pow(X, 3.0 - q) - 1.0);
      Require(std::fabs(weighted - expected) <= 1.0e-9 * std::max(1.0, expected),
              "SRC02 weighted CDF differs from the DSA target for q=" +
                  std::to_string(q));
    }
  }
  std::cout << "PASS SRC01: log-uniform importance weights are unbiased "
               "(q = 4 .. 63)\n";
  std::cout << "PASS SRC02: weighted proposal reproduces the DSA momentum CDF\n";

  // SRC03: the weight is a dimensionless density ratio, so evaluating it in
  // the scaled coordinate x = p/pmin (production) or in SI momentum p must
  // give the same value wherever the SI form is representable (strong shock,
  // q = 4.5).  This verifies that the production rescaling changes only the
  // floating-point range, not the physical weight given to a macroparticle.
  {
    SEP::Injection::Spectrum si;
    si.measure = SEP::Injection::Measure::Momentum;
    si.minimum = pmin;
    si.maximum = pmax;
    si.powerIndex = 4.5 - 2.0;
    SEP::Injection::Spectrum scaled = si;
    scaled.minimum = 1.0;
    scaled.maximum = X;
    for (double x : {1.0, 2.0, 17.0, 0.999 * X}) {
      const double wSi = Weight(si, pmin * x);
      const double wScaled = Weight(scaled, x);
      Require(std::fabs(wSi / wScaled - 1.0) < 1.0e-12,
              "SRC03 scaled and SI importance weights differ at x=" +
                  std::to_string(x));
    }
  }
  std::cout << "PASS SRC03: scaled-coordinate weights equal SI-momentum weights\n";
  return 0;
}
