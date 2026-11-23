#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/pfss_harmonics.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

// Observation-side batch adapter for the existing coronal PFSS evaluator.
// Coordinates are Carrington Cartesian metres. This program never initializes
// the particle mover or installs the derived map as a Stage 13 wind input.
namespace {
using SEP::CoronalCME::HarmonicCoefficient;
using SEP::CoronalCME::MapSample;
using SEP::CoronalCME::MonopolePolicy;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::PfssHarmonics;
using SEP::CoronalCME::Vec3;
using SEP::CoronalCME::Dot;
constexpr double pi = SEP::CoronalCME::Constants::kPi;
constexpr double solarRadiusM = 695700000.0; // IAU nominal photospheric radius.

enum class End { Photosphere, SourceSurface, Unresolved };

struct Endpoint {
  End boundary = End::Unresolved;
  Vec3 position{};
  double outwardPolarityT = 0.0;
};

Endpoint Trace(const PfssHarmonics& field, Vec3 x, double sign,
               double stepM) {
  const double outerM = field.SourceSurfaceRadiusM();
  const double toleranceM = 1.0e-6 * solarRadiusM;
  // A length cap bounds computation near nulls or exceptionally long tubes.
  // Such lines stay unresolved. Boundary steps stay inside the sphere/shell
  // by the triangle inequality, independent of the direction of curvature.
  const int maximum = static_cast<int>(12.0 * solarRadiusM / stepM) + 128;
  auto direction = [&](Vec3 position, Vec3& result) {
    const auto b = field.EvaluateCartesian(position);
    const double magnitude = Norm(b.value);
    if (!b.ok() || !std::isfinite(magnitude) || magnitude < 1.0e-30) {
      return false;
    }
    result = (sign / magnitude) * b.value;
    return true;
  };
  for (int n = 0; n < maximum; ++n) {
    const double radiusM = Norm(x);
    const double gapInner = radiusM - solarRadiusM;
    const double gapOuter = outerM - radiusM;
    Vec3 k1{}, k2{}, k3{}, k4{};
    if (!direction(x, k1)) break;
    const double incidence = Dot(k1, x) / radiusM;
    if (gapInner <= toleranceM && incidence < 0.0) {
      return {End::Photosphere, (solarRadiusM / radiusM) * x, 0.0};
    }
    if (gapOuter <= toleranceM && incidence > 0.0) {
      const auto b = field.EvaluateCartesian(x);
      return {End::SourceSurface, (outerM / radiusM) * x,
              Dot(b.value, x) / radiusM};
    }
    const double h = std::min(stepM, 0.8 * std::min(gapInner, gapOuter));
    if (!(h > 0.0) || !direction(x + (h / 2.0) * k1, k2)
        || !direction(x + (h / 2.0) * k2, k3)
        || !direction(x + h * k3, k4)) break;
    x = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }
  return {};
}

struct TraceResult {
  const char* topology;
  Endpoint outer;
};

TraceResult Classify(const PfssHarmonics& field, Vec3 seed, double stepM) {
  // An eligible open tube must actually connect the photosphere to the
  // source surface. A partially traced branch does not establish that link.
  const Endpoint forward = Trace(field, seed, 1.0, stepM);
  const Endpoint backward = Trace(field, seed, -1.0, stepM);
  if (forward.boundary == End::Photosphere) {
    if (backward.boundary == End::SourceSurface) return {"open", backward};
    if (backward.boundary == End::Photosphere) return {"closed", {}};
  }
  if (backward.boundary == End::Photosphere && forward.boundary == End::SourceSurface) {
    return {"open", forward};
  }
  return {"unresolved", {}};
}
} // namespace

int main() {
  try {
    int degree = 0;
    int numberOfSeeds = 0;
    double stepSolarRadii = 0.0;
    double outerSolarRadii = 0.0;
    if (!(std::cin >> degree >> stepSolarRadii >> outerSolarRadii >> numberOfSeeds)
        || degree < 1 || degree > 12 || numberOfSeeds < 1
        || !(stepSolarRadii >= 0.001 && stepSolarRadii <= 0.04)
        || !(outerSolarRadii >= 2.0 && outerSolarRadii <= 3.0)) {
      throw std::runtime_error("invalid PFSS diagnostic request");
    }
    // Kitt Peak's published low-resolution f map gives longitude in unit
    // degree steps and equal-area latitude bins in sin(latitude), south first.
    // Keep the observed sign and amplitude, converting gauss to tesla. A
    // monopole is incompatible with a globally divergence-free solar field.
    std::vector<MapSample> map;
    map.reserve(360 * 180);
    for (int row = 0; row < 180; ++row) {
      const double mu = -1.0 + (row + 0.5) / 90.0;
      for (int column = 0; column < 360; ++column) {
        double fieldGauss;
        if (!(std::cin >> fieldGauss) || !std::isfinite(fieldGauss)) {
          throw std::runtime_error("missing or invalid radial magnetogram pixel");
        }
        map.push_back({std::acos(mu), pi * column / 180.0,
                       fieldGauss * 1.0e-4, 4.0 * pi / (360 * 180)});
      }
    }
    std::vector<Vec3> seeds(numberOfSeeds);
    for (Vec3& seed : seeds) {
      if (!(std::cin >> seed.x >> seed.y >> seed.z)
          || !(Norm(seed) > solarRadiusM && Norm(seed) < outerSolarRadii * solarRadiusM)) {
        throw std::runtime_error("invalid diagnostic seed");
      }
    }
    SEP::CoronalCME::ReconstructionMetrics metrics;
    const auto coefficients = PfssHarmonics::ProjectMap(
        map, degree, MonopolePolicy::Remove, 0.0, &metrics);
    if (!coefficients.ok()) throw std::runtime_error(coefficients.status.message);
    const auto field = PfssHarmonics::Create(
        solarRadiusM, outerSolarRadii * solarRadiusM, coefficients.value);
    if (!field.ok()) throw std::runtime_error(field.status.message);
    std::vector<TraceResult> results(numberOfSeeds);
    // Every ray is independent. Store results by input index so thread count
    // and scheduling cannot change the numerical answer or output order.
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numberOfSeeds; ++i) {
      results[i] = Classify(field.value, seeds[i], stepSolarRadii * solarRadiusM);
    }
    std::cout << std::setprecision(17) << metrics.removedMonopoleT << '\n';
    for (const auto& result : results) {
      const auto& endpoint = result.outer;
      std::cout << result.topology << ' ' << endpoint.position.x << ' '
                << endpoint.position.y << ' ' << endpoint.position.z << ' '
                << endpoint.outwardPolarityT << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
