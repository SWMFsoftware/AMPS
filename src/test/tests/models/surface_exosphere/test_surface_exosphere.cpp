#include "test_harness.h"

#include "Exosphere_DefaultSurfaceDataStructure.h"
#include "Exosphere_SurfaceSource.h"
#include "SurfaceElementGeometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

void exit(long int line, const char* file, const char* message) {
  std::string text = AMPS::Test::Location(file, static_cast<int>(line));
  text += message == NULL ? "AMPS surface-data fatal error" : message;
  throw AMPS::Test::Failure(text);
}

namespace {
const double kPi = 3.141592653589793238462643383279502884;

template <class Function>
void ExpectInvalid(Function function) {
  bool rejected = false;
  try { function(); }
  catch (const std::invalid_argument&) { rejected = true; }
  AMPS_EXPECT_TRUE(rejected);
}

template <class Function>
void ExpectRange(Function function) {
  bool rejected = false;
  try { function(); }
  catch (const std::out_of_range&) { rejected = true; }
  AMPS_EXPECT_TRUE(rejected);
}

double Dot(const std::array<double, 3>& a,
           const std::array<double, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Norm(const std::array<double, 3>& value) {
  return std::sqrt(Dot(value, value));
}

/**
 * Adapter required by cInternalSphericalData_UserDefined::Allocate().  Both
 * its area and normal come from the production dependency-light spherical
 * geometry class, so Allocate() is tested against real geometry rather than a
 * constant-area mock.
 */
class TestSphere {
 public:
  explicit TestSphere(const Surface::Geometry::SphericalElementGrid& grid)
      : grid_(grid) {}

  double GetSurfaceElementArea(int element) const {
    return grid_.area(static_cast<std::size_t>(element));
  }

  void GetSurfaceElementNormal(double* normal, int element) const {
    const std::array<double, 3> value =
        grid_.centre_normal(static_cast<std::size_t>(element));
    std::copy(value.begin(), value.end(), normal);
  }

 private:
  const Surface::Geometry::SphericalElementGrid& grid_;
};

// The legacy default surface buffer has no destructor.  Explicit cleanup in
// this test prevents leak reports from masking arithmetic failures; it does
// not alter ownership rules in production applications.
void Release(cInternalSphericalData_UserDefined& data) {
  delete[] data.SurfaceElementDesorptionFluxUP[0];
  delete[] data.SurfaceElementDesorptionFluxUP;
  delete[] data.SurfaceElementAdsorptionFluxDOWN[0];
  delete[] data.SurfaceElementAdsorptionFluxDOWN;
  delete[] data.SurfaceElementPopulation[0];
  delete[] data.SurfaceElementPopulation;
  delete[] data.SurfaceElementExternalNormal;
  delete[] data.SurfaceElementArea;
  delete[] data.ElementSourceRate[0];
  delete[] data.ElementSourceRate;
  delete[] data.SolarWindSurfaceFlux;
  delete[] data.SampleSpeciesSurfaceSourceRate[0][0];
  delete[] data.SampleSpeciesSurfaceSourceRate[0];
  delete[] data.SampleSpeciesSurfaceSourceRate;
  delete[] data.SampleSpeciesSurfaceAreaDensity[0];
  delete[] data.SampleSpeciesSurfaceAreaDensity;
  delete[] data.SampleSpeciesSurfaceReturnFlux[0];
  delete[] data.SampleSpeciesSurfaceReturnFlux;
  delete[] data.SampleSpeciesSurfaceInjectionFlux[0];
  delete[] data.SampleSpeciesSurfaceInjectionFlux;
  delete[] data.SampleReturnFluxBulkSpeed[0];
  delete[] data.SampleReturnFluxBulkSpeed;
  delete[] data.SampleInjectedFluxBulkSpeed[0];
  delete[] data.SampleInjectedFluxBulkSpeed;
}
}  // namespace

AMPS_TEST_CASE("SURFACE-SPHERE-GEOMETRY-001",
               "spherical indexing is bijective and element areas close to 4*pi*R^2") {
  const double radius_m = 2.5;
  const Surface::Geometry::ZenithSpacing modes[] = {
      Surface::Geometry::ZenithSpacing::EqualCosine,
      Surface::Geometry::ZenithSpacing::EqualAngle};

  for (std::size_t mode = 0; mode < 2; ++mode) {
    const Surface::Geometry::SphericalElementGrid grid(7, 11, radius_m,
                                                        modes[mode]);
    double summed_area = 0.0;
    for (std::size_t azimuth = 0; azimuth < grid.azimuth_count(); ++azimuth) {
      for (std::size_t zenith = 0; zenith < grid.zenith_count(); ++zenith) {
        const std::size_t flat = grid.flatten(zenith, azimuth);
        const Surface::Geometry::ElementIndex restored = grid.unflatten(flat);
        AMPS_EXPECT_EQ(restored.zenith, zenith);
        AMPS_EXPECT_EQ(restored.azimuth, azimuth);
        AMPS_EXPECT_TRUE(grid.area(flat) > 0.0);
        AMPS_EXPECT_NEAR(Norm(grid.centre_normal(flat)), 1.0, 2.0e-15);
        summed_area += grid.area(flat);
      }
    }
    // This oracle is the analytic area of a sphere, independent of the cell
    // edge calculation used by the production grid.
    AMPS_EXPECT_NEAR(summed_area, 4.0 * kPi * radius_m * radius_m, 2.0e-13);
  }

  ExpectInvalid([]() {
    Surface::Geometry::SphericalElementGrid bad(
        0, 4, 1.0, Surface::Geometry::ZenithSpacing::EqualCosine);
    (void)bad;
  });
  ExpectInvalid([]() {
    Surface::Geometry::SphericalElementGrid bad(
        4, 4, -1.0, Surface::Geometry::ZenithSpacing::EqualCosine);
    (void)bad;
  });
  const Surface::Geometry::SphericalElementGrid grid(
      2, 3, 1.0, Surface::Geometry::ZenithSpacing::EqualCosine);
  ExpectRange([&grid]() { (void)grid.unflatten(grid.size()); });
  ExpectRange([&grid]() { (void)grid.flatten(2, 0); });
}

AMPS_TEST_CASE("SURFACE-STATE-FLUX-001",
               "default exosphere surface buffers preserve inventory fluxes and reset sampled fluxes") {
  const int species_count = 2;
  const int source_max_id = 2;
  const Surface::Geometry::SphericalElementGrid grid(
      3, 4, 10.0, Surface::Geometry::ZenithSpacing::EqualCosine);
  TestSphere sphere(grid);
  cInternalSphericalData_UserDefined data;
  data.Allocate(species_count, static_cast<int>(grid.size()), source_max_id,
                &sphere);

  for (std::size_t element = 0; element < grid.size(); ++element) {
    AMPS_EXPECT_NEAR(data.SurfaceElementArea[element], grid.area(element),
                     2.0e-14 * grid.area(element));
    const std::array<double, 3> cached = {{
        data.SurfaceElementExternalNormal[element].norm[0],
        data.SurfaceElementExternalNormal[element].norm[1],
        data.SurfaceElementExternalNormal[element].norm[2]}};
    AMPS_EXPECT_NEAR(Norm(cached), 1.0, 2.0e-15);
    AMPS_EXPECT_NEAR(data.SolarWindSurfaceFlux[element], -1.0, 0.0);
  }

  // The physical inventory accumulators are intentionally not sampling-window
  // buffers: adsorption/desorption/population survive flush().  All sampled
  // rates, counts, and speed moments must start a new window at zero.
  data.SurfaceElementAdsorptionFluxDOWN[1][5] = 7.5;
  data.SurfaceElementDesorptionFluxUP[1][5] = 2.0;
  data.SurfaceElementPopulation[1][5] = 40.0;
  data.SampleSpeciesSurfaceSourceRate[1][5][2] = 11.0;
  data.SampleSpeciesSurfaceAreaDensity[1][5] = 12.0;
  data.SampleSpeciesSurfaceReturnFlux[1][5] = 13.0;
  data.SampleSpeciesSurfaceInjectionFlux[1][5] = 14.0;
  data.SampleReturnFluxBulkSpeed[1][5] = 15.0;
  data.SampleInjectedFluxBulkSpeed[1][5] = 16.0;
  data.flush();

  AMPS_EXPECT_NEAR(data.SurfaceElementAdsorptionFluxDOWN[1][5], 7.5, 0.0);
  AMPS_EXPECT_NEAR(data.SurfaceElementDesorptionFluxUP[1][5], 2.0, 0.0);
  AMPS_EXPECT_NEAR(data.SurfaceElementPopulation[1][5], 40.0, 0.0);
  AMPS_EXPECT_NEAR(data.SampleSpeciesSurfaceSourceRate[1][5][2], 0.0, 0.0);
  AMPS_EXPECT_NEAR(data.SampleSpeciesSurfaceAreaDensity[1][5], 0.0, 0.0);
  AMPS_EXPECT_NEAR(data.SampleSpeciesSurfaceReturnFlux[1][5], 0.0, 0.0);
  AMPS_EXPECT_NEAR(data.SampleSpeciesSurfaceInjectionFlux[1][5], 0.0, 0.0);
  AMPS_EXPECT_NEAR(data.SampleReturnFluxBulkSpeed[1][5], 0.0, 0.0);
  AMPS_EXPECT_NEAR(data.SampleInjectedFluxBulkSpeed[1][5], 0.0, 0.0);
  // Net surface inventory for this deterministic fixture is the independent
  // molecule-count balance N_final=N_initial+A_down-D_up=45.5.
  AMPS_EXPECT_NEAR(data.SurfaceElementPopulation[1][5] +
                       data.SurfaceElementAdsorptionFluxDOWN[1][5] -
                       data.SurfaceElementDesorptionFluxUP[1][5],
                   45.5, 1.0e-14);

  Release(data);
}

AMPS_TEST_CASE("EXOSPHERE-SOURCE-SHAPE-001",
               "surface source shapes integrate to the requested total physical rate") {
  using Exosphere::SurfaceSource::NormalizeShape;
  const std::vector<double> area = {1.0, 2.0, 3.0, 4.0};
  const std::vector<double> uniform_shape(4, 1.0);
  const std::vector<double> uniform = NormalizeShape(uniform_shape, area, 100.0);
  AMPS_EXPECT_NEAR(uniform[0], 10.0, 2.0e-14);
  AMPS_EXPECT_NEAR(uniform[1], 20.0, 2.0e-14);
  AMPS_EXPECT_NEAR(uniform[2], 30.0, 2.0e-14);
  AMPS_EXPECT_NEAR(uniform[3], 40.0, 2.0e-14);

  const std::vector<double> shaped = NormalizeShape(
      std::vector<double>{0.0, 1.0, 2.0, 0.5}, area, 37.0);
  double total = 0.0;
  for (std::size_t i = 0; i < shaped.size(); ++i) total += shaped[i];
  AMPS_EXPECT_NEAR(total, 37.0, 1.0e-14);
  AMPS_EXPECT_NEAR(shaped[0], 0.0, 0.0);

  ExpectInvalid([&area]() {
    (void)NormalizeShape(std::vector<double>{0.0, 0.0, 0.0, 0.0}, area, 1.0);
  });
  ExpectInvalid([&area]() {
    (void)NormalizeShape(std::vector<double>{1.0, -1.0, 1.0, 1.0}, area, 1.0);
  });
  ExpectInvalid([&uniform_shape]() {
    (void)NormalizeShape(uniform_shape,
                         std::vector<double>{1.0, 0.0, 1.0, 1.0}, 1.0);
  });
  ExpectInvalid([&uniform_shape, &area]() {
    (void)NormalizeShape(uniform_shape, area,
                         std::numeric_limits<double>::infinity());
  });
}

AMPS_TEST_CASE("EXOSPHERE-ANGULAR-UNIFORM-001",
               "uniform-hemisphere quadrature reproduces analytic angular moments") {
  using Exosphere::SurfaceSource::AngularLaw;
  using Exosphere::SurfaceSource::DirectionFromUnitSamples;
  const std::array<double, 3> normal = {{2.0, -1.0, 3.0}};
  const double inverse_norm = 1.0 / Norm(normal);
  const std::array<double, 3> unit_normal = {{
      normal[0] * inverse_norm, normal[1] * inverse_norm,
      normal[2] * inverse_norm}};
  const int n_mu = 160;
  const int n_phi = 128;
  double mean_mu = 0.0;
  double mean_mu2 = 0.0;
  std::array<double, 3> mean_direction = {{0.0, 0.0, 0.0}};

  for (int i = 0; i < n_mu; ++i) {
    for (int j = 0; j < n_phi; ++j) {
      const std::array<double, 3> direction = DirectionFromUnitSamples(
          normal, (i + 0.5) / n_mu, (j + 0.5) / n_phi,
          AngularLaw::UniformSolidAngle);
      const double mu = Dot(direction, unit_normal);
      AMPS_EXPECT_NEAR(Norm(direction), 1.0, 8.0e-15);
      AMPS_EXPECT_TRUE(mu >= -2.0e-15);
      mean_mu += mu;
      mean_mu2 += mu * mu;
      for (int d = 0; d < 3; ++d) mean_direction[d] += direction[d];
    }
  }
  const double count = static_cast<double>(n_mu * n_phi);
  mean_mu /= count;
  mean_mu2 /= count;
  for (int d = 0; d < 3; ++d) mean_direction[d] /= count;
  // Analytic hemisphere moments for p(mu)=1 are <mu>=1/2 and
  // <mu^2>=1/3.  Midpoint quadrature integrates <mu> exactly.
  AMPS_EXPECT_NEAR(mean_mu, 0.5, 3.0e-14);
  AMPS_EXPECT_NEAR(mean_mu2, 1.0 / 3.0, 1.4e-5);
  for (int d = 0; d < 3; ++d)
    AMPS_EXPECT_NEAR(mean_direction[d], 0.5 * unit_normal[d], 2.0e-4);
}

AMPS_TEST_CASE("EXOSPHERE-ANGULAR-COSINE-001",
               "Knudsen cosine quadrature closes moments and handles the anti-z pole") {
  using Exosphere::SurfaceSource::AngularLaw;
  using Exosphere::SurfaceSource::DirectionFromUnitSamples;
  const std::array<double, 3> south = {{0.0, 0.0, -4.0}};
  const std::array<double, 3> unit_south = {{0.0, 0.0, -1.0}};
  const int n_mu = 400;
  const int n_phi = 64;
  double mean_mu = 0.0;
  double mean_mu2 = 0.0;
  double max_tangent_mean = 0.0;
  std::array<double, 3> mean_direction = {{0.0, 0.0, 0.0}};

  for (int i = 0; i < n_mu; ++i) {
    for (int j = 0; j < n_phi; ++j) {
      const std::array<double, 3> direction = DirectionFromUnitSamples(
          south, (i + 0.5) / n_mu, (j + 0.5) / n_phi,
          AngularLaw::KnudsenCosine);
      const double mu = Dot(direction, unit_south);
      AMPS_EXPECT_TRUE(mu >= 0.0);
      mean_mu += mu;
      mean_mu2 += mu * mu;
      for (int d = 0; d < 3; ++d) mean_direction[d] += direction[d];
    }
  }
  const double count = static_cast<double>(n_mu * n_phi);
  mean_mu /= count;
  mean_mu2 /= count;
  for (int d = 0; d < 3; ++d) mean_direction[d] /= count;
  max_tangent_mean = std::max(std::fabs(mean_direction[0]),
                              std::fabs(mean_direction[1]));

  // For p(mu)=2mu, <mu>=2/3 and <mu^2>=1/2.  The midpoint sum of sqrt(u)
  // converges algebraically, so use a deterministic quadrature tolerance.
  AMPS_EXPECT_NEAR(mean_mu, 2.0 / 3.0, 4.0e-5);
  AMPS_EXPECT_NEAR(mean_mu2, 0.5, 3.0e-14);
  AMPS_EXPECT_TRUE(max_tangent_mean < 2.0e-14);
  AMPS_EXPECT_TRUE(mean_direction[2] < 0.0);

  ExpectInvalid([]() {
    (void)Exosphere::SurfaceSource::DirectionFromUnitSamples(
        std::array<double, 3>{{0.0, 0.0, 0.0}}, 0.5, 0.5,
        Exosphere::SurfaceSource::AngularLaw::KnudsenCosine);
  });
  ExpectInvalid([&south]() {
    (void)Exosphere::SurfaceSource::DirectionFromUnitSamples(
        south, 1.0, 0.5,
        Exosphere::SurfaceSource::AngularLaw::KnudsenCosine);
  });
}

AMPS_TEST_CASE("EXOSPHERE-POISSON-SCHEDULER-001",
               "production surface events use the stationary Poisson waiting-time law") {
  using Exosphere::SurfaceSource::PoissonWaitingTime;
  const double rate_per_s = 4.0;
  const double expected_intervals_s[] = {0.0, 0.125, 0.5, 2.0};
  for (std::size_t i = 0;
       i < sizeof(expected_intervals_s) / sizeof(expected_intervals_s[0]);
       ++i) {
    // For F(t)=1-exp(-lambda*t), supplying u=F(t) must invert exactly to t.
    // The CDF is the independent analytic oracle used by the production
    // injection loop; no random sample or production-generated fixture is
    // embedded in this comparison.
    const double u = 1.0 - std::exp(-rate_per_s * expected_intervals_s[i]);
    AMPS_EXPECT_NEAR(PoissonWaitingTime(rate_per_s, u),
                     expected_intervals_s[i],
                     64.0 * std::numeric_limits<double>::epsilon() *
                         std::fmax(1.0, expected_intervals_s[i]));
  }
  AMPS_EXPECT_TRUE(std::isinf(PoissonWaitingTime(0.0, 0.25)));
  ExpectInvalid([]() {
    (void)Exosphere::SurfaceSource::PoissonWaitingTime(-1.0, 0.25);
  });
  ExpectInvalid([]() {
    (void)Exosphere::SurfaceSource::PoissonWaitingTime(1.0, 1.0);
  });
}

AMPS_TEST_CASE("EXOSPHERE-COMPILED-SOURCE-SELECTION-001",
               "compiled source-process selection follows physical-rate fractions") {
  using Exosphere::SurfaceSource::SelectSourceProcess;
  const std::vector<double> rates = {0.0, 2.0, 3.0, 5.0, 0.0};
  const std::size_t expected[] = {1, 1, 2, 2, 3, 3};
  const double samples[] = {0.0, 0.199999, 0.2, 0.499999, 0.5, 0.999999};
  for (std::size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); ++i)
    AMPS_EXPECT_EQ(SelectSourceProcess(rates, 10.0, samples[i]), expected[i]);

  // Midpoint stratification gives an exact discrete count oracle because the
  // requested probabilities 0.2, 0.3, and 0.5 align with 10,000 strata.
  std::size_t counts[5] = {0, 0, 0, 0, 0};
  for (std::size_t i = 0; i < 10000; ++i)
    ++counts[SelectSourceProcess(rates, 10.0, (i + 0.5) / 10000.0)];
  AMPS_EXPECT_EQ(counts[0], static_cast<std::size_t>(0));
  AMPS_EXPECT_EQ(counts[1], static_cast<std::size_t>(2000));
  AMPS_EXPECT_EQ(counts[2], static_cast<std::size_t>(3000));
  AMPS_EXPECT_EQ(counts[3], static_cast<std::size_t>(5000));
  AMPS_EXPECT_EQ(counts[4], static_cast<std::size_t>(0));

  // Source rates are dimensional s^-1 and may legitimately be far below
  // unity.  The selector must use a relative closure tolerance rather than an
  // absolute max(1,Q) scale that would accept order-one disagreement here.
  const std::vector<double> tiny_rates = {2.0e-30, 3.0e-30};
  AMPS_EXPECT_EQ(SelectSourceProcess(tiny_rates, 5.0e-30, 0.1),
                 static_cast<std::size_t>(0));
  AMPS_EXPECT_EQ(SelectSourceProcess(tiny_rates, 5.0e-30, 0.9),
                 static_cast<std::size_t>(1));
  ExpectInvalid([&tiny_rates]() {
    (void)SelectSourceProcess(tiny_rates, 6.0e-30, 0.5);
  });

  ExpectInvalid([&rates]() { (void)SelectSourceProcess(rates, 9.0, 0.5); });
  ExpectInvalid([]() {
    (void)Exosphere::SurfaceSource::SelectSourceProcess(
        std::vector<double>{0.0, 0.0}, 1.0, 0.5);
  });
}

AMPS_TEST_CASE("EXOSPHERE-FLUX-CLOSURE-001",
               "source, impact, inventory, and sampled duplicate ledgers close") {
  Exosphere::SurfaceSource::FluxLedger ledger;
  // Event 1 represents 12 particles over 3 s and depletes the surface.
  // Event 2 is a nondepleting external source representing 5 particles over
  // 2 s.  The independent expected diagnostic is 12/3+5/2=6.5 s^-1.
  ledger.RecordInjection(12.0, 3.0, true);
  ledger.RecordInjection(5.0, 2.0, false);
  ledger.RecordImpact(9.0, 7.0);
  ledger.RecordImpact(4.0, 1.5);

  AMPS_EXPECT_NEAR(ledger.physical_source_rate_per_s(), 6.5, 0.0);
  AMPS_EXPECT_NEAR(ledger.sampled_source_rate_per_s(), 6.5, 0.0);
  AMPS_EXPECT_NEAR(ledger.sampled_injected_particles(), 17.0, 0.0);
  AMPS_EXPECT_NEAR(ledger.desorbed_particles(), 12.0, 0.0);
  AMPS_EXPECT_NEAR(ledger.physical_returned_particles(), 13.0, 0.0);
  AMPS_EXPECT_NEAR(ledger.sampled_returned_particles(), 13.0, 0.0);
  AMPS_EXPECT_NEAR(ledger.adsorbed_particles(), 8.5, 0.0);
  AMPS_EXPECT_NEAR(ledger.inventory_change_particles(), -3.5, 0.0);
  AMPS_EXPECT_TRUE(ledger.DuplicateLedgersClose(0.0));

  ExpectInvalid([&ledger]() { ledger.RecordInjection(-1.0, 1.0, true); });
  ExpectInvalid([&ledger]() { ledger.RecordInjection(1.0, 0.0, true); });
  ExpectInvalid([&ledger]() { ledger.RecordImpact(1.0, 1.1); });
  ExpectInvalid([&ledger]() { (void)ledger.DuplicateLedgersClose(-1.0); });
}
