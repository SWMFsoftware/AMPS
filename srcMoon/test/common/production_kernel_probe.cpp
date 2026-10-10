// Linked verification probe for pure srcMoon production kernels.
//
// This executable intentionally calls the functions compiled into libAMPS.a;
// it is not a second implementation of the model.  Closed-form calculations
// below are independent acceptance oracles for deliberately simple fixtures.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

#include "pic.h"
#include "Moon.h"
#include "LunarSurface.h"
#include "MoonInput.h"

double ArgonStickingProbability(double &reemission_fraction, double temperature_k);
double SodiumStickingProbability(double &reemission_fraction, double temperature_k);
double localResolution(double *position_m);
double localSphericalSurfaceResolution(double *position_m);
double SodiumRadiationPressureAcceleration__Combi_1997_icarus(
    double heliocentric_velocity_m_s, double heliocentric_distance_m);

namespace {

bool close(double actual, double expected, double relative_tolerance,
    double absolute_tolerance = 0.0) {
  // Combine relative and absolute tolerances so zero-valued invariants can be
  // tested without dividing by zero.  Individual contracts choose tolerances;
  // this helper never relaxes them after observing a result.
  return std::abs(actual - expected) <=
      std::max(absolute_tolerance, relative_tolerance * std::abs(expected));
}

int report(const char *name, double actual, double expected,
    double relative_tolerance, double absolute_tolerance = 0.0) {
  const bool passed = close(actual, expected, relative_tolerance,
      absolute_tolerance);
  std::printf("%s %s actual=%.17e expected=%.17e\n",
      passed ? "PASS" : "FAIL", name, actual, expected);
  return passed ? 0 : 1;
}

int sticking() {
  // Exercise compiled production lookup/piecewise functions at frozen table
  // or branch control points.  These checks do not sample the stochastic
  // Bernoulli collision outcome or establish a residence-time model.
  int failures = 0;
  double reemission_fraction = -1.0;

  failures += report("sodium_sticking_100K",
      SodiumStickingProbability(reemission_fraction, 100.0), 0.99830,
      1.0e-14);
  failures += report("sodium_reemission_fraction", reemission_fraction,
      1.0, 0.0);

  failures += report("argon_sticking_88K",
      ArgonStickingProbability(reemission_fraction, 88.0),
      std::pow(10.0, -0.72), 1.0e-14);
  failures += report("argon_sticking_110K",
      ArgonStickingProbability(reemission_fraction, 110.0),
      std::pow(10.0, -2.2), 1.0e-14);
  failures += report("argon_sticking_158K",
      ArgonStickingProbability(reemission_fraction, 158.0),
      1.0e-4, 1.0e-14);
  failures += report("argon_reemission_fraction", reemission_fraction,
      1.0, 0.0);
  return failures == 0 ? 0 : 1;
}

int radiation_shadow() {
  // Coordinates are metres in the Moon-centred solar-orbital frame used by
  // EarthShadowCheck().  The fixture isolates axial, radial, and sunward
  // branches of the production cylindrical-shadow predicate.
  const double earth_radius = _RADIUS_(_EARTH_);
  double behind_earth[3] = {-12.0 * earth_radius, 0.0, 0.0};
  double off_axis[3] = {-12.0 * earth_radius, 2.0 * earth_radius, 0.0};
  double sunward_of_earth[3] = {0.0, 0.0, 0.0};

  Exosphere::xEarth_SO[0] = -10.0 * earth_radius;
  Exosphere::xEarth_SO[1] = 0.0;
  Exosphere::xEarth_SO[2] = 0.0;
  Exosphere::xSun_SO[0] = 1000.0 * earth_radius;
  Exosphere::xSun_SO[1] = 0.0;
  Exosphere::xSun_SO[2] = 0.0;

  int failures = 0;
  failures += report("earth_umbra_axis",
      Moon::EarthShadowCheck(behind_earth) ? 1.0 : 0.0, 1.0, 0.0);
  failures += report("earth_umbra_off_axis",
      Moon::EarthShadowCheck(off_axis) ? 1.0 : 0.0, 0.0, 0.0);
  failures += report("earth_umbra_sunward",
      Moon::EarthShadowCheck(sunward_of_earth) ? 1.0 : 0.0, 0.0, 0.0);

  const double acceleration_1au =
      SodiumRadiationPressureAcceleration__Combi_1997_icarus(0.0, _AU_);
  const double acceleration_2au =
      SodiumRadiationPressureAcceleration__Combi_1997_icarus(0.0, 2.0 * _AU_);
  failures += report("radiation_pressure_inverse_square",
      acceleration_1au / acceleration_2au, 4.0, 1.0e-14);
  return failures == 0 ? 0 : 1;
}

int sodium_sources() {
  // The rate is the configured one-AU production value.  The PDF expectation
  // is evaluated independently at 1 eV and compared with the production
  // function compiled into libAMPS.a; this is not a total source-budget test.
  Exosphere::xObjectRadial = _AU_;
  int species = _NA_SPEC_;
  const double energy_j = 1.0 * eV2J;
  const double exponent = 0.7;
  const double binding_energy_ev = 0.052;
  const double expected_psd = exponent * (1.0 + exponent) *
      std::pow(binding_energy_ev, exponent) /
      std::pow(1.0 + binding_energy_ev, 2.0 + exponent);

  int failures = 0;
  failures += report("impact_vaporization_rate_s-1",
      Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(
          species, _INTERNAL_BOUNDARY_TYPE_SPHERE_, nullptr),
      1.69e22, 1.0e-14);
  failures += report("psd_energy_pdf_at_1eV",
      Exosphere::SourceProcesses::PhotonStimulatedDesorption::
          EnergyDistributionFunction(energy_j, &species),
      expected_psd, 1.0e-14);
  return failures == 0 ? 0 : 1;
}

int temperature() {
  // Verify closed-form points of the active analytic temperature law.  Earth
  // is placed away from the solar ray so the fixture tests temperature rather
  // than the separate terrestrial-shadow branch.
  int failures = 0;
  double surface_position[3] = {_RADIUS_(_MOON_), 0.0, 0.0};

  // Keep the Earth far from the ray used by this fixture. Coordinates are in
  // the Moon-centred solar-orbital frame and use metres.
  Exosphere::xSun_SO[0] = 1.0e11;
  Exosphere::xSun_SO[1] = 0.0;
  Exosphere::xSun_SO[2] = 0.0;
  Exosphere::xEarth_SO[0] = 0.0;
  Exosphere::xEarth_SO[1] = 1.0e9;
  Exosphere::xEarth_SO[2] = 0.0;

  failures += report("nightside_temperature",
      Exosphere::GetSurfaceTemperature(-0.25, surface_position), 100.0, 0.0);
  failures += report("terminator_temperature",
      Exosphere::GetSurfaceTemperature(0.0, surface_position), 100.0, 0.0);
  failures += report("subsolar_temperature",
      Exosphere::GetSurfaceTemperature(1.0, surface_position), 380.0, 0.0);
  failures += report("quarter_cosine_temperature",
      Exosphere::GetSurfaceTemperature(0.0625, surface_position), 240.0,
      1.0e-14);
  return failures == 0 ? 0 : 1;
}

int gravity() {
  // At x=2 lunar radii with zero velocity and orbit terms compiled out, the
  // production acceleration must reduce to the point-mass identity -GM/r^2.
  double position_m[3] = {2.0 * _RADIUS_(_MOON_), 0.0, 0.0};
  double velocity_m_s[3] = {0.0, 0.0, 0.0};
  double acceleration_m_s2[3] = {0.0, 0.0, 0.0};

  Moon::TotalParticleAcceleration(acceleration_m_s2, _NA_SPEC_, -1,
      position_m, velocity_m_s, nullptr);

  const double radius_m = position_m[0];
  const double expected_x =
      -GravityConstant * _MASS_(_MOON_) / (radius_m * radius_m);

  int failures = 0;
  failures += report("gravity_x", acceleration_m_s2[0], expected_x, 1.0e-14);
  failures += report("gravity_y", acceleration_m_s2[1], 0.0, 0.0);
  failures += report("gravity_z", acceleration_m_s2[2], 0.0, 0.0);
  return failures == 0 ? 0 : 1;
}

int photochemistry() {
  // Current orbit-off Na mode uses the legacy constant lifetime.  This local
  // check is not evidence of daughter-ion production or a D04-driven network.
  double position_m[3] = {2.0 * _RADIUS_(_MOON_), 0.0, 0.0};
  bool allowed = false;
  const double actual = Moon::ExospherePhotoionizationLifeTime(position_m,
      _NA_SPEC_, -1, allowed, nullptr);
  const double expected = 3600.0 * 5.8 / std::pow(0.4, 2.0);

  int failures = 0;
  failures += report("sodium_photo_lifetime_s", actual, expected, 1.0e-14);
  failures += report("sodium_photo_allowed", allowed ? 1.0 : 0.0, 1.0,
      0.0);
  return failures == 0 ? 0 : 1;
}

int mesh_resolution() {
  // Freeze the effective callbacks in the current build.  Their constant
  // values document baseline wiring but do not establish mesh convergence.
  double position_m[3] = {_RADIUS_(_MOON_), 0.0, 0.0};
  const double expected_surface =
      _RADIUS_(_MOON_) * 4.0 * 4.0 / 100.0;
  const double expected_volume = _RADIUS_(_MOON_) * 4.0 * 2.0;

  int failures = 0;
  failures += report("surface_resolution_m",
      localSphericalSurfaceResolution(position_m), expected_surface, 1.0e-14);
  failures += report("volume_resolution_m", localResolution(position_m),
      expected_volume, 1.0e-14);
  return failures == 0 ? 0 : 1;
}

int lola_geometry() {
  // U08 is intentionally coarse (500 km) so decoding and topology checks are
  // fast.  The linked smoke run may reuse this receipt, but the local probe
  // itself does not initialize AMR or move particles.
  namespace fs = std::filesystem;
  const char *output_environment=std::getenv("MOON_U08_OUTPUT");
  const fs::path output_directory=output_environment==nullptr ?
      fs::path("test_output/srcMoon/unit/U08/surface_probe") :
      fs::path(output_environment);
  fs::create_directories(output_directory);

  const fs::path input_file=output_directory/"u08.in";
  const fs::path cea_file=output_directory/"lola_surface.cea";
  const fs::path tecplot_file=output_directory/"lola_surface.dat";
  {
    // Include an unrelated section to prove the production parser respects
    // application ownership.  Absolute raw-data paths identify the exact D01
    // files; generated products remain under the selected test output root.
    std::ofstream input(input_file);
    input << "#section begin: unrelated\nignored = yes\n#section end\n"
          << "#section begin: moon\n"
          << "spice_path = /home/vtenishe/SPICE\n"
          << "surface_geometry = lola\n"
          << "surface_mesh_resolution_m = 500000\n"
          << "lola_product_id = LDEM_4\n"
          << "lola_image_file = "
          << "/data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG\n"
          << "lola_label_file = "
          << "/data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL\n"
          << "surface_cea_file = " << cea_file.string() << "\n"
          << "surface_tecplot_file = " << tecplot_file.string() << "\n"
          << "#section end\n";
  }

  Moon::Runtime::Configuration configuration;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(input_file.string(),
                                             &configuration,&error)) {
    std::fprintf(stderr,"FAIL moon_input_parser %s\n",error.c_str());
    return 1;
  }

  Moon::Surface::LolaDem dem;
  if (!dem.Load(configuration,&error)) {
    std::fprintf(stderr,"FAIL lola_load %s\n",error.c_str());
    return 1;
  }

  int failures=0;
  failures+=report("label_lines",dem.metadata().lines,720.0,0.0);
  failures+=report("label_samples",dem.metadata().lineSamples,1440.0,0.0);
  failures+=report("label_reference_radius_m",
      dem.metadata().referenceRadiusM,1737400.0,0.0);
  // These four pixel-centre values were decoded independently and frozen
  // before running the model.  They span both poles, the prime-meridian seam,
  // the equator, and the antimeridian, detecting endian/axis/wrap mistakes.
  failures+=report("native_pixel_north_west_m",
      dem.ElevationM(0.125,89.875),-119.5,0.0);
  failures+=report("native_pixel_equator_west_m",
      dem.ElevationM(0.125,-0.125),-721.5,0.0);
  failures+=report("native_pixel_equator_antimeridian_m",
      dem.ElevationM(180.125,-0.125),2836.5,0.0);
  failures+=report("native_pixel_south_east_m",
      dem.ElevationM(359.875,-89.875),91.0,0.0);

  Moon::Surface::Triangulation mesh;
  if (!Moon::Surface::BuildIcosphere(dem.metadata().referenceRadiusM,
          configuration.surfaceMeshResolutionM,&mesh,&error) ||
      !Moon::Surface::ApplyLolaTopography(dem,&mesh,&error)) {
    std::fprintf(stderr,"FAIL lola_triangulation %s\n",error.c_str());
    return 1;
  }
  failures+=report("icosphere_vertices",mesh.verticesM.size(),642.0,0.0);
  failures+=report("icosphere_faces",mesh.faces.size(),1280.0,0.0);
  if (!(mesh.maximumEdgeLengthM<=configuration.surfaceMeshResolutionM)) {
    std::fprintf(stderr,"FAIL maximum_edge actual=%.17e limit=%.17e\n",
        mesh.maximumEdgeLengthM,configuration.surfaceMeshResolutionM);
    ++failures;
  }

  // Check geometry invariants independently from the writer/loader: every
  // face must point away from the origin and global chord lengths must retain
  // the declared approximately uniform icosphere character.
  double minimum_edge=1.0e100,maximum_edge=0.0;
  for (const auto& face : mesh.faces) {
    const auto& a=mesh.verticesM[face[0]];
    const auto& b=mesh.verticesM[face[1]];
    const auto& c=mesh.verticesM[face[2]];
    const double ab[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
    const double ac[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
    const double normal[3]={ab[1]*ac[2]-ab[2]*ac[1],
        ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]};
    const double center[3]={a[0]+b[0]+c[0],a[1]+b[1]+c[1],
        a[2]+b[2]+c[2]};
    if (normal[0]*center[0]+normal[1]*center[1]+normal[2]*center[2]<=0.0) {
      ++failures;
    }
    const std::array<const std::array<double,3>*,3> point={&a,&b,&c};
    for (int edge=0;edge<3;edge++) {
      const auto& p=*point[edge];
      const auto& q=*point[(edge+1)%3];
      const double length=std::sqrt((p[0]-q[0])*(p[0]-q[0])+
          (p[1]-q[1])*(p[1]-q[1])+(p[2]-q[2])*(p[2]-q[2]));
      minimum_edge=std::min(minimum_edge,length);
      maximum_edge=std::max(maximum_edge,length);
    }
  }
  if (minimum_edge/maximum_edge<0.75) {
    std::fprintf(stderr,"FAIL edge_uniformity min_max_ratio=%.17e\n",
        minimum_edge/maximum_edge);
    ++failures;
  }
  else {
    std::printf("PASS edge_uniformity min_max_ratio=%.17e\n",
        minimum_edge/maximum_edge);
  }

  if (!Moon::Surface::WriteCeaSurface(mesh,cea_file.string(),&error) ||
      !Moon::Surface::WriteTecplotSurface(mesh,tecplot_file.string(),&error)) {
    std::fprintf(stderr,"FAIL surface_writers %s\n",error.c_str());
    return 1;
  }
  // Re-open both artifacts instead of trusting write success alone.  U08
  // checks record counts and Tecplot topology declaration; the linked smoke
  // run separately demonstrates that AMPS can consume the CEA product.
  std::ifstream cea(cea_file);
  std::size_t written_vertices=0,written_faces=0;
  cea >> written_vertices >> written_faces;
  failures+=report("cea_vertices",written_vertices,mesh.verticesM.size(),0.0);
  failures+=report("cea_faces",written_faces,mesh.faces.size(),0.0);
  std::ifstream tecplot(tecplot_file);
  std::string variables,zone;
  std::getline(tecplot,variables);
  std::getline(tecplot,zone);
  if (variables.find("VARIABLES=")!=0 || zone.find("ZONETYPE=FETRIANGLE")==
      std::string::npos) {
    std::fprintf(stderr,"FAIL tecplot_header\n");
    ++failures;
  }
  else std::printf("PASS tecplot_header\n");

  return failures==0 ? 0 : 1;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr,
        "usage: production_kernel_probe <gravity|mesh-resolution|"
        "photochemistry|radiation-shadow|sodium-sources|sticking|"
        "temperature|lola-geometry>\n");
    return 2;
  }

  const std::string test(argv[1]);
  if (test == "gravity") return gravity();
  if (test == "mesh-resolution") return mesh_resolution();
  if (test == "lola-geometry") return lola_geometry();
  if (test == "photochemistry") return photochemistry();
  if (test == "radiation-shadow") return radiation_shadow();
  if (test == "sodium-sources") return sodium_sources();
  if (test == "sticking") return sticking();
  if (test == "temperature") return temperature();

  std::fprintf(stderr, "unknown test: %s\n", argv[1]);
  return 2;
}
