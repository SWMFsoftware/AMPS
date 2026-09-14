#include "../../util/sep_flux_tube_geometry_core.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

int failures = 0;

void Check(bool condition, const std::string& id, const std::string& detail) {
  if (!condition) {
    std::cerr << "FAIL " << id << ": " << detail << '\n';
    ++failures;
  } else {
    std::cout << "PASS " << id << ": " << detail << '\n';
  }
}

bool Near(double actual, double expected, double relative_tolerance) {
  return std::fabs(actual - expected) <=
         relative_tolerance * std::max(std::fabs(expected), 1.0);
}

}  // namespace

int main() {
  using SEP::FieldLine::FluxTubeGeometryCore::AreaFromMagneticFlux;
  using SEP::FieldLine::FluxTubeGeometryCore::InjectedPhysicalParticleCount;
  using SEP::FieldLine::FluxTubeGeometryCore::IntegrateLinearArea;
  using SEP::FieldLine::FluxTubeGeometryCore::IntegrateSimpsonArea;
  using SEP::FieldLine::FluxTubeGeometryCore::MacroparticleWeight;
  using SEP::FieldLine::FluxTubeGeometryCore::SpectralNumberDensityPerJ;
  using SEP::FieldLine::FluxTubeGeometryCore::SweptVolume;
  using namespace SEP::Units;

  // GEOA01: all sampled areas must carry exactly the same magnetic flux.
  const AreaM2 reference_area_m2(12.5);
  const MagneticFieldT reference_B_T(8.0e-8);
  const double B_samples_T[] = {2.0e-8, 4.0e-8, 8.0e-8, 1.6e-7};
  bool flux_conserved = true;
  for (std::size_t i = 0; i < sizeof(B_samples_T) / sizeof(B_samples_T[0]); ++i) {
    const AreaM2 area_m2 = AreaFromMagneticFlux(
        reference_area_m2, reference_B_T, MagneticFieldT(B_samples_T[i]));
    flux_conserved = flux_conserved &&
        Near(area_m2.Value() * B_samples_T[i],
             reference_area_m2.Value() * reference_B_T.Value(), 1.0e-14);
  }
  Check(flux_conserved, "GEOA01", "A|B| is invariant");

  // GEOA02: the production Simpson area integration must converge at fourth
  // order for a smooth area profile.  Each cell samples the same begin/mid/end
  // area contract used by the PIC adapter.
  const double exact_volume_m3 = std::exp(1.0) - 1.0;
  double previous_error = 0.0;
  bool converges = true;
  for (int level = 0; level < 4; ++level) {
    const int cells = 8 << level;
    double volume_m3 = 0.0;
    for (int i = 0; i < cells; ++i) {
      const double x0 = static_cast<double>(i) / cells;
      const double x1 = static_cast<double>(i + 1) / cells;
      volume_m3 += IntegrateSimpsonArea(
          AreaM2(std::exp(x0)), AreaM2(std::exp(0.5*(x0+x1))),
          AreaM2(std::exp(x1)), LengthM(1.0 / cells)).Value();
    }
    const double error = std::fabs(volume_m3 - exact_volume_m3);
    if (level > 0) converges = converges && previous_error / error > 15.0;
    previous_error = error;
  }
  Check(converges, "GEOA02", "segment-volume error converges at fourth order");

  // SRC01: 4 m2 * 500 km/s * 2 s = 4,000,000 m3.  The typed arguments make
  // an accidental radius, distance, or energy substitution a compile error.
  const VolumeM3 swept_m3 =
      SweptVolume(AreaM2(4.0), SpeedMPerS(500000.0), TimeS(2.0));
  Check(Near(swept_m3.Value(), 4.0e6, 1.0e-14), "SRC01",
        "swept volume is A v dt in SI");

  // SRC02: provider identity must not change a source assembled from the same
  // density, geometry, time, speed, and efficiency.  Dividing into a different
  // number of computational particles only changes each particle's weight.
  const double physical_count = InjectedPhysicalParticleCount(
      NumberDensityPerM3(5.0e6), swept_m3, 3.4e-4);
  const double analytic_weight = MacroparticleWeight(physical_count, 100);
  const double swcme_weight = MacroparticleWeight(physical_count, 250);
  const double swmf_weight = MacroparticleWeight(physical_count, 400);
  Check(Near(100.0 * analytic_weight, physical_count, 1.0e-14) &&
            Near(250.0 * swcme_weight, physical_count, 1.0e-14) &&
            Near(400.0 * swmf_weight, physical_count, 1.0e-14),
        "SRC02", "provider-equivalent sources have equal physical weight");

  // SRC03: a normalized piecewise-constant probability density must integrate
  // back to the physical number density used by sampling and comparisons.
  const double bin_width_J = 2.0e-13;
  const int bins = 5;
  const double probability_per_J = 1.0 / (bins * bin_width_J);
  const NumberDensityPerM3 density_per_m3(7.5e5);
  double reconstructed_density_per_m3 = 0.0;
  for (int i = 0; i < bins; ++i) {
    reconstructed_density_per_m3 +=
        SpectralNumberDensityPerJ(density_per_m3, probability_per_J) *
        bin_width_J;
  }
  Check(Near(reconstructed_density_per_m3, density_per_m3.Value(), 1.0e-14),
        "SRC03", "spectral bins reconstruct number density");

  const EnergyJ one_MeV_J = SEP::Units::EnergyFromMeV(1.0);
  Check(Near(one_MeV_J.Value(), 1.602176634e-13, 1.0e-14), "SRC03-UNITS",
        "MeV-to-joule conversion is explicit and exact");

  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
