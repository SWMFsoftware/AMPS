// Linked verification probe for pure srcMoon production kernels.
//
// This executable intentionally calls the functions compiled into libAMPS.a;
// it is not a second implementation of the model.  Closed-form calculations
// below are independent acceptance oracles for deliberately simple fixtures.

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "pic.h"
#include "Moon.h"
#include "LunarSurface.h"
#include "MoonInput.h"
#include "MoonSurfaceEnvironment.h"
#include "MoonSurfacePhysics.h"
#include "MoonSurfaceTemperature.h"
#include "MoonDivinerObservation.h"
#include "MoonThermalModel.h"
#include "species/Argon.h"

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

int report_string(const char *name, const std::string& actual,
                  const char *expected) {
  const bool passed=actual==expected;
  std::printf("%s %s actual=%s expected=%s\n",
      passed ? "PASS" : "FAIL",name,actual.c_str(),expected);
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

int radiation_pressure_values(int count, char **velocity_tokens) {
  // This mode is deliberately only a production-kernel adapter.  The
  // authoritative Combi-curve values and their digitization uncertainty live
  // in U05's reference package and are evaluated by the Python runner.  Keeping
  // expected values out of this executable prevents the embedded production
  // table from becoming its own oracle.
  if (count <= 0) {
    std::fprintf(stderr,
        "radiation-pressure-values requires at least one velocity in m/s\n");
    return 2;
  }

  for (int i = 0; i < count; ++i) {
    errno = 0;
    char *end = nullptr;
    const double velocity_m_s = std::strtod(velocity_tokens[i], &end);
    if (errno != 0 || end == velocity_tokens[i] || *end != '\0' ||
        !std::isfinite(velocity_m_s)) {
      std::fprintf(stderr, "invalid heliocentric velocity: %s\n",
          velocity_tokens[i]);
      return 2;
    }

    // The publication curve and the compiled table are both normalized to
    // one astronomical unit.  The separate radiation_shadow() fixture tests
    // the production kernel's inverse-square scaling away from 1 AU.
    const double acceleration_m_s2 =
        SodiumRadiationPressureAcceleration__Combi_1997_icarus(
            velocity_m_s, _AU_);
    std::printf(
        "VALUE sodium_radiation_pressure velocity_m_s=%.17e "
        "acceleration_m_s2=%.17e\n",
        velocity_m_s, acceleration_m_s2);
  }

  return 0;
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

int argon_source() {
  // Exercise the species-specific production kernel independently of the
  // unavailable global potassium geography.  Areas deliberately differ so
  // this detects an equal-facet implementation that would bias latitude on
  // the logical sphere.  Conservation is the independent oracle: every
  // element receives A_i/sum(A) of the declared global rate.
  Moon::Species::Argon::SourceParameters parameters;
  parameters.enabled=true;
  parameters.spatialModel="uniform-global";
  parameters.timeModel="continuous";
  parameters.totalRateAtomsS=3.8e20;
  parameters.reference=
      "Killen_2002_doi_10.1111/j.1945-5100.2002.tb00891.x_lower_bound";
  std::string error;
  int failures=0;
  if (!Moon::Species::Argon::ValidateSourceParameters(parameters,&error)) {
    std::fprintf(stderr,"FAIL argon_source_parameters %s\n",error.c_str());
    return 1;
  }

  const std::vector<double> areasM2{1.0,2.0,3.0};
  std::vector<double> ratesAtomsS;
  if (!Moon::Species::Argon::BuildUniformElementRates(
          areasM2,parameters.totalRateAtomsS,&ratesAtomsS,&error)) {
    std::fprintf(stderr,"FAIL argon_source_partition %s\n",error.c_str());
    return 1;
  }
  failures+=report("argon_rate_area_1",ratesAtomsS[0],
      parameters.totalRateAtomsS/6.0,1.0e-15);
  failures+=report("argon_rate_area_2",ratesAtomsS[1],
      parameters.totalRateAtomsS/3.0,1.0e-15);
  failures+=report("argon_rate_area_3",ratesAtomsS[2],
      parameters.totalRateAtomsS/2.0,1.0e-15);
  failures+=report("argon_rate_conservation",
      std::accumulate(ratesAtomsS.begin(),ratesAtomsS.end(),0.0),
      parameters.totalRateAtomsS,1.0e-15);

  // Unsupported geography and physically invalid areas must fail closed.
  parameters.spatialModel="potassium-weighted";
  failures+=report("argon_reject_regional_k_map",
      Moon::Species::Argon::ValidateSourceParameters(parameters,&error)?1.0:0.0,
      0.0,0.0);
  std::vector<double> invalidAreasM2{1.0,-1.0};
  failures+=report("argon_reject_negative_area",
      Moon::Species::Argon::BuildUniformElementRates(invalidAreasM2,3.8e20,
          &ratesAtomsS,&error)?1.0:0.0,0.0,0.0);

  // Parse the maintained M2 receipt through the production parser.  The
  // standard build remains Na-only and therefore records the candidate source
  // parameters while explicitly disabling injection.
  Moon::Runtime::Configuration configuration;
  const std::filesystem::path receipt=
      std::filesystem::path("srcMoon")/"examples"/"lola_surface.in";
  if (!Moon::Runtime::ParseApplicationInput(receipt.string(),&configuration,
                                            &error)) {
    std::fprintf(stderr,"FAIL argon_source_receipt %s\n",error.c_str());
    return 1;
  }
  failures+=report("argon_source_enabled",
      configuration.argonSourceEnabled?1.0:0.0,0.0,0.0);
  failures+=report_string("argon_source_spatial_model",
      configuration.argonSourceSpatialModel,"uniform-global");
  failures+=report("argon_source_rate_atoms_s",
      configuration.argonSourceRateAtomsS,3.8e20,0.0);
  failures+=report_string("argon_source_time_model",
      configuration.argonSourceTimeModel,"continuous");
  failures+=report_string("argon_source_map_file",
      configuration.argonSourceMapFile,"none");
  return failures==0?0:1;
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

int terrain_shadow() {
  // Construct a deterministic body-fixed fixture with one horizontal source
  // facet and a vertical ridge east of it.  Expected visibility follows from
  // elementary line/plane geometry and does not invoke a second ray tracer.
  Moon::Surface::Triangulation mesh;
  mesh.verticesM={{{-1.0,-1.0,0.0}},{{1.0,-1.0,0.0}},{{0.0,1.0,0.0}},
      {{1.0,-1.0,0.0}},{{1.0,1.0,0.0}},{{1.0,1.0,2.0}},
      {{1.0,-1.0,2.0}}};
  mesh.faces={{{0,1,2}},{{3,4,5}},{{3,5,6}}};
  Moon::SurfaceEnvironment::TerrainRayTracer tracer;
  std::string error;
  if (!tracer.Initialize(mesh,&error)) {
    std::fprintf(stderr,"FAIL terrain_initialize %s\n",error.c_str());
    return 1;
  }
  const std::array<double,3> point{{0.0,0.0,0.0}};
  const std::array<double,3> normal{{0.0,0.0,1.0}};
  bool visible=false;
  int failures=0;
  if (!tracer.IsSunVisible(point,normal,{{0.0,0.0,1.0}},1.0e-6,0,
                           &visible,&error)) {
    std::fprintf(stderr,"FAIL overhead_visibility %s\n",error.c_str());
    return 1;
  }
  failures+=report("overhead_sun_visible",visible?1.0:0.0,1.0,0.0);
  if (!tracer.IsSunVisible(point,normal,{{1.0,0.0,0.5}},1.0e-6,0,
                           &visible,&error)) {
    std::fprintf(stderr,"FAIL ridge_visibility %s\n",error.c_str());
    return 1;
  }
  failures+=report("ridge_blocks_low_sun",visible?1.0:0.0,0.0,0.0);
  if (!tracer.IsSunVisible(point,normal,{{0.0,0.0,-1.0}},1.0e-6,0,
                           &visible,&error)) {
    std::fprintf(stderr,"FAIL backface_visibility %s\n",error.c_str());
    return 1;
  }
  failures+=report("backface_self_shadow",visible?1.0:0.0,0.0,0.0);
  Moon::SurfaceEnvironment::RayHit hit;
  if (!tracer.IntersectRay(point,{{1.0,0.0,0.5}},1.0e-9,100.0,0,&hit,
                           &error)) {
    std::fprintf(stderr,"FAIL ridge_intersection %s\n",error.c_str());
    return 1;
  }
  failures+=report("ridge_hit_distance_m",hit.distanceM,
      std::sqrt(1.25),1.0e-14);

  // Independently verify the analytical overlap branches used by production
  // eclipse forcing.  These fixtures are dimensionless circle geometry, so
  // they do not depend on SPICE ephemerides or an observed eclipse.
  const double angularRadius=0.01;
  failures+=report("solar_disc_no_overlap",
      Moon::SurfaceTemperature::UniformSolarDiscVisibleFraction(
          angularRadius,angularRadius,2.1*angularRadius),1.0,0.0);
  failures+=report("solar_disc_totality",
      Moon::SurfaceTemperature::UniformSolarDiscVisibleFraction(
          angularRadius,2.0*angularRadius,0.0),0.0,0.0);
  failures+=report("solar_disc_concentric_small_occulter",
      Moon::SurfaceTemperature::UniformSolarDiscVisibleFraction(
          angularRadius,0.5*angularRadius,0.0),0.75,1.0e-14);
  const double pi=std::acos(-1.0);
  const double expectedEqualDiscVisible=1.0-
      (2.0*pi/3.0-std::sqrt(3.0)/2.0)/pi;
  failures+=report("solar_disc_equal_radius_partial",
      Moon::SurfaceTemperature::UniformSolarDiscVisibleFraction(
          angularRadius,angularRadius,angularRadius),
      expectedEqualDiscVisible,1.0e-14);
  return failures==0?0:1;
}

int diviner_temperature() {
  // The fixture is a complete 2x2x2 grid with an affine temperature field.
  // Trilinear interpolation must reproduce that field exactly, while a
  // non-GOOD corner must make the enclosing cell unavailable.
  const char* outputRoot=std::getenv("MOON_U10_OUTPUT");
  if (outputRoot==nullptr) {
    std::fprintf(stderr,"MOON_U10_OUTPUT is required\n");
    return 2;
  }
  const std::filesystem::path directory(outputRoot);
  std::filesystem::create_directories(directory);
  const auto tablePath=directory/"diviner_fixture.csv";
  {
    std::ofstream output(tablePath);
    output << "latitude_deg,longitude_east_deg,local_time_hour,temperature_k,"
              "uncertainty_k,quality_flag,product_id,source_record\n";
    for (double latitude:{-10.0,10.0}) for (double longitude:{0.0,20.0})
      for (double localTime:{0.0,4.0}) {
        const double temperature=200.0+latitude+0.5*longitude+2.0*localTime;
        output << latitude << ',' << longitude << ',' << localTime << ','
               << temperature << ",1,GOOD,D03_FIXTURE,row-"
               << latitude << '-' << longitude << '-' << localTime << '\n';
      }
  }
  Moon::SurfaceTemperature::DivinerTable table;
  std::string error;
  if (!table.LoadCsv(tablePath.string(),&error)) {
    std::fprintf(stderr,"FAIL diviner_load %s\n",error.c_str());
    return 1;
  }
  double temperature=0.0,uncertainty=0.0;
  if (!table.Interpolate(0.0,10.0,2.0,&temperature,&uncertainty,&error)) {
    std::fprintf(stderr,"FAIL diviner_interpolate %s\n",error.c_str());
    return 1;
  }
  int failures=0;
  failures+=report("diviner_affine_midpoint_k",temperature,209.0,1.0e-14);
  failures+=report("diviner_midpoint_uncertainty_k",uncertainty,
      std::sqrt(8.0)/8.0,1.0e-14);
  if (!table.Interpolate(-10.0,0.0,0.0,&temperature,&uncertainty,&error)) {
    std::fprintf(stderr,"FAIL diviner_native_node %s\n",error.c_str());
    return 1;
  }
  failures+=report("diviner_native_node_k",temperature,190.0,0.0);

  // Verify the production observation operator separately from table
  // interpolation.  This fixture uses no observations: its oracle is the
  // Stefan--Boltzmann identity for a grey surface and the monotonic inverse of
  // Planck radiance.  The deliberately sub-threshold C3--C8 temperatures also
  // expose the exact Paige et al. longest-wavelength substitution chain.
  constexpr double kineticTemperatureK=240.0;
  constexpr double emissivity=0.81;
  const double expectedBolometricK=std::pow(emissivity,0.25)*
      kineticTemperatureK;
  failures+=report("diviner_grey_bolometric_temperature_k",
      Moon::DivinerObservation::GreyBolometricBrightnessTemperatureK(
          kineticTemperatureK,emissivity),expectedBolometricK,1.0e-14);
  failures+=report("diviner_grey_exitance_w_m2",
      Moon::DivinerObservation::GreySurfaceExitanceWm2(
          kineticTemperatureK,emissivity),
      emissivity*5.670374419e-8*std::pow(kineticTemperatureK,4),1.0e-14);

  for (int channel=6;channel<=9;++channel) {
    const double radiance=
        Moon::DivinerObservation::RectangularBandRadianceWm2Sr(
            channel,kineticTemperatureK,1.0);
    double recoveredTemperatureK=0.0;
    if (!Moon::DivinerObservation::RectangularBandBrightnessTemperatureK(
            channel,radiance,&recoveredTemperatureK,&error)) {
      std::fprintf(stderr,"FAIL diviner_channel_inverse_c%d %s\n",channel,
          error.c_str());
      ++failures;
      continue;
    }
    const std::string metric="diviner_channel_inverse_c"+
        std::to_string(channel)+"_k";
    failures+=report(metric.c_str(),recoveredTemperatureK,
        kineticTemperatureK,1.0e-12);
  }

  const std::array<double,Moon::DivinerObservation::InfraredChannelCount>
      channelTemperatureK{{100.0,110.0,120.0,80.0,50.0,30.0,55.0}};
  std::array<double,Moon::DivinerObservation::InfraredChannelCount>
      temperaturesUsedK{};
  double tbolK=0.0;
  if (!Moon::DivinerObservation::BolometricBrightnessTemperatureK(
          channelTemperatureK,&tbolK,&temperaturesUsedK,&error)) {
    std::fprintf(stderr,"FAIL diviner_tbol_operator %s\n",error.c_str());
    return 1;
  }
  for (std::size_t index=0;index<temperaturesUsedK.size();++index) {
    // Every C3--C8 fixture value is below its channel-specific threshold, so
    // the publication rule must carry the valid C9 value backward through the
    // entire chain.  This checks rule ordering, not observational agreement.
    failures+=report(("diviner_sensitivity_substitution_c"+
        std::to_string(static_cast<int>(index)+3)+"_k").c_str(),
        temperaturesUsedK[index],55.0,0.0);
  }
  // Independent high-resolution Python quadrature of the published 0--1000
  // um partitions gives 54.988557963911575 K.  It is slightly below 55 K
  // because Paige et al. deliberately terminate C9 at 1000 um rather than
  // assigning the unmeasured >1000-um Planck tail to Diviner.
  failures+=report("diviner_tbol_substitution_uniform_k",tbolK,
      54.988557963911575,0.0,2.0e-6);
  return failures==0?0:1;
}

Moon::Thermal::ColumnConfiguration thermal_fixture_configuration() {
  Moon::Thermal::ColumnConfiguration configuration;
  configuration.solarConstantWm2=1361.0;
  configuration.infraredEmissivity=0.95;
  configuration.normalBondAlbedo=0.12;
  configuration.albedoCoefficientA=0.06;
  configuration.albedoCoefficientB=0.25;
  configuration.surfaceDensityKgM3=1100.0;
  configuration.deepDensityKgM3=1800.0;
  configuration.densityScaleHeightM=0.06;
  configuration.surfaceConductivityWmK=7.4e-4;
  configuration.deepConductivityWmK=3.4e-3;
  configuration.radiativeConductivityParameter=2.7;
  configuration.heatCapacityPolynomial={-3.6125,2.7431,2.3616e-3,
                                         -1.2340e-5,8.9093e-9};
  configuration.geothermalFluxWm2=0.018;
  configuration.bottomDepthM=1.0;
  configuration.layerCount=40;
  configuration.nonlinearToleranceK=1.0e-10;
  configuration.nonlinearMaximumIterations=100;
  return configuration;
}

int thermal_inertia() {
  Moon::Thermal::Column column;
  auto configuration=thermal_fixture_configuration();
  std::string error;
  if (!column.Initialize(configuration,250.0,&error)) {
    std::fprintf(stderr,"FAIL thermal_initialize %s\n",error.c_str());
    return 1;
  }
  int failures=0;
  failures+=report("thermal_density_surface",
      column.DensityKgM3(0.0),1100.0,0.0);
  failures+=report("thermal_density_one_scale_height",
      column.DensityKgM3(0.06),1800.0-700.0/std::exp(1.0),1.0e-14);
  failures+=report("thermal_conductivity_surface_350k",
      column.ConductivityWmK(0.0,350.0),7.4e-4*(1.0+2.7),1.0e-14);
  const double before=column.SurfaceTemperatureK();
  if (!column.Advance(60.0,0.0,0.0,1.0,&error)) {
    std::fprintf(stderr,"FAIL thermal_eclipse_step %s\n",error.c_str());
    return 1;
  }
  if (!(column.SurfaceTemperatureK()<before)) {
    std::fprintf(stderr,"FAIL thermal_eclipse_cooling before=%.17e after=%.17e\n",
        before,column.SurfaceTemperatureK());
    ++failures;
  }
  else std::printf("PASS thermal_eclipse_cooling actual=%.17e expected_less_than=%.17e\n",
      column.SurfaceTemperatureK(),before);
  const double boundaryEnergy=60.0*(configuration.geothermalFluxWm2-
      configuration.infraredEmissivity*5.670374419e-8*
      std::pow(column.SurfaceTemperatureK(),4));
  failures+=report("thermal_energy_relative_residual",
      column.LastEnergyResidualJm2()/std::max(1.0,std::abs(boundaryEnergy)),
      0.0,0.0,2.0e-5);
  return failures==0?0:1;
}

Moon::SurfacePhysics::SpeciesParameters argon_parameters() {
  Moon::SurfacePhysics::SpeciesParameters parameters;
  parameters.species=Moon::SurfacePhysics::Species::Argon;
  parameters.residenceLaw=
      Moon::SurfacePhysics::ResidenceLaw::ModifiedArrhenius;
  parameters.activationEnergyJmol=28000.0;
  parameters.stickingTemperatureK={10.0,500.0};
  parameters.stickingProbability={1.0,1.0};
  parameters.citation="Kegerreis et al. 2017, doi:10.1002/2017JE005352";
  return parameters;
}

int surface_physics() {
  const auto argon=argon_parameters();
  const double temperature=100.0;
  const double expected=6.62607015e-34/(1.380649e-23*temperature)*
      std::exp(28000.0/(8.31446261815324*temperature));
  int failures=0;
  failures+=report("argon_modified_arrhenius_residence_s",
      Moon::SurfacePhysics::ResidenceTimeS(argon,temperature),expected,
      1.0e-14);
  failures+=report("argon_sticking_explicit",
      Moon::SurfacePhysics::StickingProbability(argon,temperature),1.0,0.0);
  failures+=report("argon_desorption_probability",
      Moon::SurfacePhysics::DesorptionProbability(argon,temperature,5.0),
      -std::expm1(-5.0/expected),1.0e-14);

  Moon::SurfacePhysics::SpeciesParameters sodium;
  sodium.species=Moon::SurfacePhysics::Species::Sodium;
  sodium.residenceLaw=Moon::SurfacePhysics::ResidenceLaw::FirstOrderArrhenius;
  sodium.attemptFrequencyHz=1.0e13;
  sodium.activationEnergyJ=1.85*eV2J;
  sodium.stickingTemperatureK={100.0,250.0,500.0};
  sodium.stickingProbability={1.0,0.5,0.2};
  failures+=report("sodium_sticking_250k_explicit",
      Moon::SurfacePhysics::StickingProbability(sodium,250.0),0.5,0.0);
  failures+=report("sodium_sticking_375k_linear",
      Moon::SurfacePhysics::StickingProbability(sodium,375.0),0.35,1.0e-14);

  // Exercise the same generic scalar kernel called after a production
  // non-sticking collision. The incident and thermal speeds differ so the
  // historical linear-speed blend cannot accidentally satisfy the published
  // Tenishev et al. (2013) kinetic-energy identity.
  constexpr double incidentSpeed=2000.0;
  constexpr double thermalSpeed=1000.0;
  constexpr double omega=0.2;
  failures+=report("accommodation_linear_speed_m_s",
      Exosphere::SurfaceInteraction::EvaluateAccommodationSpeed(
          incidentSpeed,thermalSpeed,omega,
          Exosphere::SurfaceInteraction::AccommodationSpeedLaw::LinearSpeed),
      1800.0,1.0e-14);
  failures+=report("accommodation_tenishev2013_energy_m_s",
      Exosphere::SurfaceInteraction::EvaluateAccommodationSpeed(
          incidentSpeed,thermalSpeed,omega,
          Exosphere::SurfaceInteraction::AccommodationSpeedLaw::KineticEnergy),
      std::sqrt(omega*thermalSpeed*thermalSpeed+
          (1.0-omega)*incidentSpeed*incidentSpeed),1.0e-14);
  failures+=report("accommodation_zero_limit_m_s",
      Exosphere::SurfaceInteraction::EvaluateAccommodationSpeed(
          incidentSpeed,thermalSpeed,0.0,
          Exosphere::SurfaceInteraction::AccommodationSpeedLaw::KineticEnergy),
      incidentSpeed,0.0);
  failures+=report("accommodation_unity_limit_m_s",
      Exosphere::SurfaceInteraction::EvaluateAccommodationSpeed(
          incidentSpeed,thermalSpeed,1.0,
          Exosphere::SurfaceInteraction::AccommodationSpeedLaw::KineticEnergy),
      thermalSpeed,0.0);

  // Parse the committed M2 example through the production input parser.  The
  // scalar checks above verify the accommodation mathematics; these checks
  // separately prove that the advertised species selections actually reach
  // the runtime configuration instead of remaining documentation-only keys.
  // __FILE__ is compiled from an absolute path by moon_testlib, so resolving
  // relative to this source file is independent of the runner's working
  // directory.
  const std::filesystem::path m2Input=
      std::filesystem::path(__FILE__).parent_path()/".."/".."/"examples"/
      "lola_surface.in";
  Moon::Runtime::Configuration configuration;
  std::string inputError;
  if (!Moon::Runtime::ParseApplicationInput(m2Input.lexically_normal().string(),
                                             &configuration,&inputError)) {
    std::fprintf(stderr,"FAIL m2_surface_input %s\n",inputError.c_str());
    return 1;
  }
  failures+=report_string("sodium_sticking_model",
      configuration.sodiumStickingModel,"temperature-table");
  failures+=report_string("sodium_accommodation_speed_model",
      configuration.sodiumAccommodationSpeedModel,"kinetic-energy");
  failures+=report("sodium_accommodation_coefficient",
      configuration.sodiumAccommodationCoefficient,0.2,0.0);
  failures+=report("sodium_sticking_table_nodes",
      static_cast<double>(configuration.sodiumStickingProbability.size()),
      80.0,0.0);
  failures+=report_string("argon_sticking_model",
      configuration.argonStickingModel,"constant");
  failures+=report_string("argon_accommodation_speed_model",
      configuration.argonAccommodationSpeedModel,"kinetic-energy");
  failures+=report("argon_accommodation_coefficient",
      configuration.argonAccommodationCoefficient,1.0,0.0);
  return failures==0?0:1;
}

int cold_trapping() {
  const auto argon=argon_parameters();
  Moon::SurfacePhysics::Inventory inventory;
  Moon::SurfacePhysics::Budget budget;
  std::string error;
  if (!inventory.Initialize(1,1000.0,&error)) {
    std::fprintf(stderr,"FAIL inventory_initialize %s\n",error.c_str());
    return 1;
  }
  budget.initialInventory=1000.0;
  double coldRelease=0.0;
  if (!inventory.AdvanceCell(0,argon,40.0,3600.0,86400.0,&coldRelease,
                             &budget,&error)) {
    std::fprintf(stderr,"FAIL cold_inventory_step %s\n",error.c_str());
    return 1;
  }
  const double beforeWarm=inventory.cells()[0].adsorbedAtoms;
  double warmRelease=0.0;
  if (!inventory.AdvanceCell(0,argon,150.0,3600.0,86400.0,&warmRelease,
                             &budget,&error)) {
    std::fprintf(stderr,"FAIL warm_inventory_step %s\n",error.c_str());
    return 1;
  }
  budget.finalInventory=inventory.cells()[0].adsorbedAtoms;
  budget.finalExosphere=budget.desorbed;
  int failures=0;
  if (!(warmRelease>coldRelease && warmRelease>0.0 &&
        inventory.cells()[0].adsorbedAtoms<beforeWarm)) {
    std::fprintf(stderr,"FAIL reversible_cold_trap cold=%.17e warm=%.17e\n",
        coldRelease,warmRelease);
    ++failures;
  }
  else std::printf("PASS reversible_cold_trap cold=%.17e warm=%.17e\n",
      coldRelease,warmRelease);
  failures+=report("inventory_budget_closure_atoms",budget.ClosureResidual(),
      0.0,0.0,1.0e-10);
  return failures==0?0:1;
}

int gravity() {
  // At x=2 lunar radii the isolated production lunar-gravity kernel must
  // reduce to the point-mass identity -GM/r^2.  Testing the isolated kernel
  // keeps U03 valid in both orbit-off and orbit-on production builds.
  double position_m[3] = {2.0 * _RADIUS_(_MOON_), 0.0, 0.0};
  double acceleration_m_s2[3] = {0.0, 0.0, 0.0};

  Moon::OrbitalDynamics::AddLunarPointMassAcceleration(acceleration_m_s2,
      position_m);

  const double radius_m = position_m[0];
  const double expected_x =
      -GravityConstant * _MASS_(_MOON_) / (radius_m * radius_m);

  int failures = 0;
  failures += report("gravity_x", acceleration_m_s2[0], expected_x, 1.0e-14);
  failures += report("gravity_y", acceleration_m_s2[1], 0.0, 0.0);
  failures += report("gravity_z", acceleration_m_s2[2], 0.0, 0.0);
  return failures == 0 ? 0 : 1;
}

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
void independent_cross(const double *left, const double *right, double *out) {
  // This test-only vector primitive is the independent oracle used below; it
  // does not implement or replace a production force law.
  out[0]=left[1]*right[2]-left[2]*right[1];
  out[1]=left[2]*right[0]-left[0]*right[2];
  out[2]=left[0]*right[1]-left[1]*right[0];
}

int rotating_frame() {
  // Load the exact kernel list compiled into Exosphere.h and evaluate the
  // frozen M0 epoch.  Missing kernels/frames are a CSPICE runtime error, not a
  // condition that this scientific test is allowed to convert into PASS.
  Exosphere::Init_SPICE();
  SpiceDouble epoch;
  utc2et_c("2009-01-24T00:00:00",&epoch);
  int failures=0;
  failures+=report("inertial_frame_is_J2000",
      std::strcmp(Moon::Frames::Inertial,"J2000")==0 ? 1.0 : 0.0,1.0,0.0);
  failures+=report("body_frame_is_MOON_ME_DE421",
      std::strcmp(Exosphere::IAU_FRAME,"MOON_ME_DE421")==0 ? 1.0 : 0.0,
      1.0,0.0);
  failures+=report("solar_orbital_frame_is_LSO",
      std::strcmp(Exosphere::SO_FRAME,"LSO")==0 ? 1.0 : 0.0,1.0,0.0);

  // Differential gravity uses a non-axial fixture so every component and the
  // subtraction of the lunar-origin acceleration are observable.
  const double particle_m[3]={2.1e6,-3.2e5,4.7e5};
  const double attractor_m[3]={1.23e11,-4.56e10,7.89e9};
  const double attractor_mass_kg=_MASS_(_SUN_);
  double actual_differential[3]={0.0,0.0,0.0};
  double expected_differential[3]={0.0,0.0,0.0};
  Moon::OrbitalDynamics::AddDifferentialPointMassAcceleration(
      actual_differential,particle_m,attractor_m,attractor_mass_kg);
  double origin_distance2=0.0,particle_distance2=0.0;
  for (int i=0;i<3;i++) {
    origin_distance2+=attractor_m[i]*attractor_m[i];
    const double delta=attractor_m[i]-particle_m[i];
    particle_distance2+=delta*delta;
  }
  const double origin_distance3=origin_distance2*std::sqrt(origin_distance2);
  const double particle_distance3=
      particle_distance2*std::sqrt(particle_distance2);
  for (int i=0;i<3;i++) {
    expected_differential[i]=GravityConstant*attractor_mass_kg*(
        (attractor_m[i]-particle_m[i])/particle_distance3-
        attractor_m[i]/origin_distance3);
    const std::string name="differential_gravity_"+std::to_string(i);
    failures+=report(name.c_str(),actual_differential[i],
        expected_differential[i],1.0e-12,1.0e-18);
  }

  // Independently form -omega x (omega x r) - 2 omega x v.  The arbitrary
  // non-coplanar fixture detects sign, component-order, and omitted-term bugs.
  const double position_m[3]={1.9e6,-2.4e5,8.1e5};
  const double velocity_m_s[3]={730.0,-410.0,95.0};
  const double omega_rad_s[3]={-3.1e-9,7.2e-8,2.0e-7};
  double actual_rotating[3]={0.0,0.0,0.0};
  double omega_cross_r[3],omega_cross_omega_cross_r[3],omega_cross_v[3];
  double expected_rotating[3];
  Moon::OrbitalDynamics::AddRotatingFrameAcceleration(actual_rotating,
      position_m,velocity_m_s,omega_rad_s);
  independent_cross(omega_rad_s,position_m,omega_cross_r);
  independent_cross(omega_rad_s,omega_cross_r,
      omega_cross_omega_cross_r);
  independent_cross(omega_rad_s,velocity_m_s,omega_cross_v);
  for (int i=0;i<3;i++) {
    expected_rotating[i]=-omega_cross_omega_cross_r[i]-2.0*omega_cross_v[i];
    const std::string name="rotating_acceleration_"+std::to_string(i);
    failures+=report(name.c_str(),actual_rotating[i],expected_rotating[i],
        1.0e-12,1.0e-18);
  }

  // The production angular velocity comes from xf2rav.  The independent
  // reference below extracts dR/dt directly from the lower-left block of the
  // LSO->J2000 state transform and evaluates W=(dR/dt)R^T.  These are distinct
  // algorithms operating on the same authoritative SPICE frame definition.
  double actual_omega_lso[3];
  Moon::OrbitalDynamics::GetSolarOrbitalAngularVelocityLSO(epoch,
      actual_omega_lso);
  SpiceDouble state_transform[6][6],rotation_j2000_to_lso[3][3];
  sxform_c("LSO","J2000",epoch,state_transform);
  pxform_c("J2000","LSO",epoch,rotation_j2000_to_lso);
  double skew_j2000[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
    skew_j2000[i][j]+=state_transform[i+3][k]*state_transform[j][k];
  }
  const double omega_j2000[3]={skew_j2000[2][1],skew_j2000[0][2],
      skew_j2000[1][0]};
  double expected_omega_lso[3];
  mxv_c(rotation_j2000_to_lso,omega_j2000,expected_omega_lso);
  for (int i=0;i<3;i++) {
    const std::string name="spice_angular_velocity_"+std::to_string(i);
    failures+=report(name.c_str(),actual_omega_lso[i],expected_omega_lso[i],
        1.0e-12,1.0e-20);
  }

  // Orthogonality is evaluated independently on the position rotation block.
  // This detects a malformed or incorrectly indexed state transform.
  double maximum_orthogonality_error=0.0;
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
    double dot=0.0;
    for (int k=0;k<3;k++) {
      dot+=rotation_j2000_to_lso[i][k]*rotation_j2000_to_lso[j][k];
    }
    maximum_orthogonality_error=std::max(maximum_orthogonality_error,
        std::abs(dot-(i==j ? 1.0 : 0.0)));
  }
  failures+=report("rotation_orthogonality",maximum_orthogonality_error,0.0,
      0.0,1.0e-13);
  return failures==0 ? 0 : 1;
}
#endif

int photochemistry() {
  // This local check is not evidence of daughter-ion production or a
  // D04-driven network. In an orbit-enabled build initialize the declared
  // SPICE state so Earth-shadow gating uses physical Sun/Earth positions.
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  Exosphere::Init_SPICE();
  // Init_SPICE establishes the inertial ephemeris receipt; the normal
  // amps_time_step hook subsequently fills the Moon-centred LSO force state.
  // This thin fixture performs that same state query without advancing PIC,
  // so U14 isolates illumination gating rather than running a particle step.
  SpiceDouble sun_state_lso[6],earth_state_lso[6],light_time;
  spkezr_c("SUN",Exosphere::OrbitalMotion::et,Moon::Frames::SolarOrbital,
      Moon::Frames::ForceAberrationCorrection,Moon::ObjectName,
      sun_state_lso,&light_time);
  spkezr_c("EARTH",Exosphere::OrbitalMotion::et,Moon::Frames::SolarOrbital,
      Moon::Frames::ForceAberrationCorrection,Moon::ObjectName,
      earth_state_lso,&light_time);
  for (int i=0;i<3;i++) {
    Moon::xSun_SO[i]=1.0e3*sun_state_lso[i];
    Moon::xEarth_SO[i]=1.0e3*earth_state_lso[i];
  }
#endif
  double sunlit_position_m[3] = {2.0 * _RADIUS_(_MOON_), 0.0, 0.0};
  double shadow_position_m[3] = {-2.0 * _RADIUS_(_MOON_), 0.0, 0.0};
  bool sunlit_allowed = false;
  const double actual = Moon::ExospherePhotoionizationLifeTime(
      sunlit_position_m, _NA_SPEC_, -1, sunlit_allowed, nullptr);
  const double expected = 3600.0 * 5.8 / std::pow(0.4, 2.0);

  int failures = 0;
  failures += report("sodium_photo_lifetime_s", actual, expected, 1.0e-14);
  failures += report("sodium_photo_sunlit_allowed",
      sunlit_allowed ? 1.0 : 0.0, 1.0, 0.0);
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  bool shadow_allowed = true;
  const double shadow_lifetime = Moon::ExospherePhotoionizationLifeTime(
      shadow_position_m, _NA_SPEC_, -1, shadow_allowed, nullptr);
  failures += report("sodium_photo_shadow_lifetime_s", shadow_lifetime,
      -1.0, 0.0);
  failures += report("sodium_photo_shadow_allowed",
      shadow_allowed ? 1.0 : 0.0, 0.0, 0.0);
#endif
  return failures == 0 ? 0 : 1;
}

int mesh_resolution() {
  // U07 calls the same input-driven production kernels used by the AMR
  // callbacks.  Frame transformation and realized-cell placement remain the
  // linked I12 responsibility; this local layer independently fixes the
  // piecewise interpolation, terminator override, polar minimum, and radial
  // transition arithmetic before a costly mesh is built.
  Moon::Runtime::Configuration configuration;
  configuration.meshNearSurfaceCellSizeM=10000.0;
  configuration.meshNearSurfaceExtentM=200000.0;
  configuration.meshFarFieldCellSizeM=40000.0;
  configuration.meshFarFieldExtentM=10000000.0;
  configuration.meshOuterCellSizeM=17371000.0;
  configuration.meshAngularSzaDeg={0.0,30.0,60.0,90.0,120.0,180.0};
  configuration.meshAngularResolutionM={10000.0,10000.0,5000.0,5000.0,
                                         10000.0,10000.0};
  configuration.meshTerminatorHalfWidthDeg=15.0;
  configuration.meshPolarCapStartLatitudeDeg=75.0;
  configuration.meshPolarResolutionM=5000.0;

  double position_m[3] = {_RADIUS_(_MOON_), 0.0, 0.0};
  const double expected_surface =
      _RADIUS_(_MOON_) * 4.0 * 4.0 / 100.0;
  const double expected_volume = _RADIUS_(_MOON_) * 4.0 * 2.0;

  int failures = 0;
  failures += report("surface_resolution_m",
      localSphericalSurfaceResolution(position_m), expected_surface, 1.0e-14);
  failures += report("volume_resolution_m", localResolution(position_m),
      expected_volume, 1.0e-14);
  const double szaDeg[6]={0.0,30.0,60.0,90.0,120.0,180.0};
  const double expectedAngularM[6]={10000.0,10000.0,5000.0,5000.0,
                                    10000.0,10000.0};
  for (int i=0;i<6;i++) {
    const std::string metric="configured_sza_"+std::to_string(i)+
        "_resolution_m";
    failures+=report(metric.c_str(),
        Moon::Runtime::EvaluateAngularMeshResolutionM(
            configuration,szaDeg[i],0.0),expectedAngularM[i],1.0e-14);
  }
  // SZA=45 linearly interpolates to 7.5 km.  SZA=80 lies inside the approved
  // +/-15-degree terminator band and is therefore capped at the explicit
  // 90-degree value of 5 km.  A polar point independently receives 5 km.
  failures+=report("configured_sza_45_resolution_m",
      Moon::Runtime::EvaluateAngularMeshResolutionM(
          configuration,45.0,0.0),7500.0,1.0e-14);
  failures+=report("configured_terminator_band_resolution_m",
      Moon::Runtime::EvaluateAngularMeshResolutionM(
          configuration,80.0,0.0),5000.0,1.0e-14);
  failures+=report("configured_polar_resolution_m",
      Moon::Runtime::EvaluateAngularMeshResolutionM(
          configuration,0.0,80.0),5000.0,1.0e-14);

  const double referenceRadiusM=_RADIUS_(_MOON_);
  const double outerRadiusM=400.0*referenceRadiusM;
  failures+=report("configured_volume_surface_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(
          configuration,referenceRadiusM,referenceRadiusM,outerRadiusM),
      10000.0,1.0e-14);
  failures+=report("configured_volume_extent_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(configuration,
          referenceRadiusM+200000.0,referenceRadiusM,outerRadiusM),
      10000.0,1.0e-14);
  failures+=report("configured_volume_far_extent_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(configuration,
          referenceRadiusM+10000000.0,referenceRadiusM,outerRadiusM),
      40000.0,1.0e-14);
  const double nearTransitionRadiusM=referenceRadiusM+200000.0;
  const double farTransitionRadiusM=referenceRadiusM+10000000.0;
  failures+=report("configured_volume_science_log_midpoint_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(configuration,
          std::sqrt(nearTransitionRadiusM*farTransitionRadiusM),
          referenceRadiusM,outerRadiusM),25000.0,1.0e-14);
  failures+=report("configured_volume_outer_log_midpoint_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(configuration,
          std::sqrt(farTransitionRadiusM*outerRadiusM),
          referenceRadiusM,outerRadiusM),8705500.0,1.0e-14);
  failures+=report("configured_volume_outer_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(
          configuration,outerRadiusM,referenceRadiusM,outerRadiusM),
      17371000.0,1.0e-14);

  // QUICK changes resolution only: physical extents and angular-region
  // boundaries are identical to FULL, while every length target is tripled.
  // Factor two was measured to retain one excessive dyadic AMR level, so the
  // explicit factor-three receipt is now part of the production-kernel
  // contract rather than an untested performance assumption.
  // Exercise the production selectors with that independently constructed
  // receipt so the fast profile cannot accidentally retain one FULL-scale
  // branch and create an unexpectedly expensive mixed mesh.
  Moon::Runtime::Configuration quick=configuration;
  quick.meshNearSurfaceCellSizeM=30000.0;
  quick.meshFarFieldCellSizeM=120000.0;
  quick.meshOuterCellSizeM=52113000.0;
  quick.meshAngularResolutionM={30000.0,30000.0,15000.0,15000.0,
                                30000.0,30000.0};
  quick.meshPolarResolutionM=15000.0;
  failures+=report("quick_terminator_resolution_m",
      Moon::Runtime::EvaluateAngularMeshResolutionM(quick,80.0,0.0),
      15000.0,1.0e-14);
  failures+=report("quick_polar_resolution_m",
      Moon::Runtime::EvaluateAngularMeshResolutionM(quick,0.0,80.0),
      15000.0,1.0e-14);
  failures+=report("quick_volume_surface_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(
          quick,referenceRadiusM,referenceRadiusM,outerRadiusM),
      30000.0,1.0e-14);
  failures+=report("quick_volume_far_extent_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(quick,
          referenceRadiusM+10000000.0,referenceRadiusM,outerRadiusM),
      120000.0,1.0e-14);
  failures+=report("quick_volume_outer_resolution_m",
      Moon::Runtime::EvaluateVolumeMeshResolutionM(
          quick,outerRadiusM,referenceRadiusM,outerRadiusM),
      52113000.0,1.0e-14);
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
    // Exercise the same generic ``application: moon`` block consumed by the
    // production executable. Absolute raw-data paths identify the exact D01
    // files; generated products remain under the selected test output root.
    // The shared grammar itself is independently tested in
    // src/general/input/test_amps_input.cpp; this probe verifies Moon's owner
    // dispatch and its existing physical/domain validation.
    std::ofstream input(input_file);
    input << "#begin application: moon\n"
          << "spice_path = /home/vtenishe/SPICE\n"
          << "surface_geometry = lola\n"
          << "surface_mesh_resolution_m = 500000\n"
          << "mesh_profile = custom\n"
          << "mesh_near_surface_cell_size_m = 10000\n"
          << "mesh_near_surface_extent_m = 200000\n"
          << "mesh_far_field_cell_size_m = 40000\n"
          << "mesh_far_field_extent_m = 10000000\n"
          << "mesh_outer_cell_size_m = 17371000\n"
          << "mesh_output_tecplot = false\n"
          << "mesh_angular_resolution_m_by_sza = "
          << "0:10000,30:10000,60:5000,90:5000,120:10000,180:10000\n"
          << "mesh_terminator_half_width_deg = 15\n"
          << "mesh_polar_cap_start_latitude_deg = 75\n"
          << "mesh_polar_resolution_m = 5000\n"
          << "lola_product_id = LDEM_4\n"
          << "lola_image_file = "
          << "/data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG\n"
          << "lola_label_file = "
          << "/data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL\n"
          << "surface_cea_file = " << cea_file.string() << "\n"
          << "surface_tecplot_file = " << tecplot_file.string() << "\n"
          << "m2_environment = disabled\n"
          << "simulation_start_utc = 2009-01-24T00:00:00\n"
          << "simulation_environment_mode = frozen\n"
          << "random_seed = 731993\n"
          << "profile_altitude_near_maximum_m = 200000\n"
          << "profile_altitude_near_step_m = 5000\n"
          << "profile_altitude_middle_maximum_m = 1000000\n"
          << "profile_altitude_middle_step_m = 25000\n"
          << "profile_altitude_far_maximum_m = 10000000\n"
          << "profile_altitude_far_points_per_decade = 8\n"
          << "profile_solar_zenith_angles_deg = 0,30,60,90,120,180\n"
          << "profile_include_body_fixed_poles = true\n"
          << "profile_include_shell_average = true\n"
          << "profile_shell_angular_weighting = solid-angle\n"
          << "limb_observer = earth\n"
          << "limb_plane = sun-moon-earth\n"
          << "limb_sides = both\n"
          << "limb_branch_reference = projected-sun\n"
          << "limb_tangent_maximum_altitude_m = 1000000\n"
          << "tail_center_direction = anti-lunar\n"
          << "tail_map_coordinate_system = j2000-ra-dec\n"
          << "tail_map_display_orientation = west-right\n"
          << "tail_map_half_width_deg = 7.5\n"
          << "tail_beam_diameter_deg = 1.0\n"
          << "tail_exposure_duration_s = 120\n"
          << "tail_beam_response = unresolved\n"
          << "tail_temporal_operator = exposure-average\n"
          << "#end application\n";
  }

  Moon::Runtime::Configuration configuration;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(input_file.string(),
                                             &configuration,&error)) {
    std::fprintf(stderr,"FAIL moon_input_parser %s\n",error.c_str());
    return 1;
  }

  // The migration is intentionally strict: the retired #section spelling
  // must fail instead of being accepted through an undocumented compatibility
  // path. This prevents two Moon grammars from drifting independently.
  const fs::path legacy_input_file=output_directory/"u08-legacy-section.in";
  std::ifstream generic_input(input_file);
  std::string legacy_text((std::istreambuf_iterator<char>(generic_input)),
                           std::istreambuf_iterator<char>());
  const std::size_t generic_begin=legacy_text.find("#begin application: moon");
  const std::size_t generic_end=legacy_text.find("#end application");
  if (generic_begin==std::string::npos || generic_end==std::string::npos) {
    std::fprintf(stderr,"FAIL legacy_input_fixture missing generic boundary\n");
    return 1;
  }
  legacy_text.replace(generic_end,std::string("#end application").size(),
                      "#section end");
  legacy_text.replace(generic_begin,
                      std::string("#begin application: moon").size(),
                      "#section begin: moon");
  {
    std::ofstream legacy_input(legacy_input_file);
    legacy_input << legacy_text;
  }
  Moon::Runtime::Configuration rejected_configuration;
  if (Moon::Runtime::ParseApplicationInput(legacy_input_file.string(),
                                            &rejected_configuration,&error) ||
      error.find("unknown directive '#section'")==std::string::npos) {
    std::fprintf(stderr,
        "FAIL legacy #section input was not rejected by generic grammar: %s\n",
        error.c_str());
    return 1;
  }

  // Parse a second, independently named receipt with the other D14-backed
  // coordinate/display choices and one selected limb branch.  This proves the
  // selectors are genuinely adjustable through .in rather than constants
  // that happen to be serialized in the primary fixture.  It does not run an
  // observation operator or turn these parser checks into U26/I33 evidence.
  const fs::path alternate_input_file=output_directory/"u08-alternate.in";
  std::ifstream primary_input(input_file);
  std::string alternate_text((std::istreambuf_iterator<char>(primary_input)),
                              std::istreambuf_iterator<char>());
  const auto replace_once=[&alternate_text](const std::string& from,
                                             const std::string& to) {
    const std::size_t location=alternate_text.find(from);
    if (location==std::string::npos) return false;
    alternate_text.replace(location,from.size(),to);
    return true;
  };
  if (!replace_once("limb_sides = both",
                    "limb_sides = away-from-projected-sun") ||
      !replace_once("tail_map_coordinate_system = j2000-ra-dec",
                    "tail_map_coordinate_system = eclipj2000-lon-lat") ||
      !replace_once("tail_map_display_orientation = west-right",
                    "tail_map_display_orientation = east-right")) {
    std::fprintf(stderr,"FAIL alternate_input_fixture replacement failed\n");
    return 1;
  }
  {
    std::ofstream alternate_input(alternate_input_file);
    alternate_input << alternate_text;
  }
  Moon::Runtime::Configuration alternate_configuration;
  if (!Moon::Runtime::ParseApplicationInput(alternate_input_file.string(),
                                             &alternate_configuration,&error)) {
    std::fprintf(stderr,"FAIL alternate_moon_input_parser %s\n",error.c_str());
    return 1;
  }

  // The steady/transient choice is part of the physical run contract, not a
  // free-form provenance label. Prove that an unknown value fails before AMR
  // initialization instead of silently reverting to the legacy evolving
  // behavior.
  const fs::path invalid_environment_file=
      output_directory/"u08-invalid-environment.in";
  std::string invalid_environment_text=alternate_text;
  const std::string valid_environment="simulation_environment_mode = frozen";
  const std::size_t environment_location=
      invalid_environment_text.find(valid_environment);
  if (environment_location==std::string::npos) {
    std::fprintf(stderr,"FAIL invalid_environment_fixture missing key\n");
    return 1;
  }
  invalid_environment_text.replace(environment_location,
      valid_environment.size(),"simulation_environment_mode = assumed");
  {
    std::ofstream invalid_environment(invalid_environment_file);
    invalid_environment << invalid_environment_text;
  }
  Moon::Runtime::Configuration invalid_configuration;
  if (Moon::Runtime::ParseApplicationInput(invalid_environment_file.string(),
                                            &invalid_configuration,&error)) {
    std::fprintf(stderr,
        "FAIL invalid simulation_environment_mode was accepted\n");
    return 1;
  }

  // The visualization switch must be a real Boolean.  Rejecting a misspelled
  // value prevents an expensive run from silently enabling the historical
  // full-volume text dump and exhausting the output filesystem.
  const fs::path invalid_mesh_output_file=
      output_directory/"u08-invalid-mesh-output.in";
  std::string invalid_mesh_output_text=alternate_text;
  const std::string valid_mesh_output="mesh_output_tecplot = false";
  const std::size_t mesh_output_location=
      invalid_mesh_output_text.find(valid_mesh_output);
  if (mesh_output_location==std::string::npos) {
    std::fprintf(stderr,"FAIL invalid_mesh_output_fixture missing key\n");
    return 1;
  }
  invalid_mesh_output_text.replace(mesh_output_location,
      valid_mesh_output.size(),"mesh_output_tecplot = automatic");
  {
    std::ofstream invalid_mesh_output(invalid_mesh_output_file);
    invalid_mesh_output << invalid_mesh_output_text;
  }
  if (Moon::Runtime::ParseApplicationInput(invalid_mesh_output_file.string(),
                                            &invalid_configuration,&error)) {
    std::fprintf(stderr,"FAIL invalid mesh_output_tecplot was accepted\n");
    return 1;
  }

  // QUICK, FULL, and CUSTOM are the only declared profile identities. The
  // profile is provenance rather than a source of hidden resolution values,
  // but accepting an old or misspelled label would make run manifests
  // ambiguous about which committed receipt was intended.
  const fs::path invalid_profile_file=
      output_directory/"u08-invalid-profile.in";
  std::string invalid_profile_text=alternate_text;
  const std::string valid_profile="mesh_profile = custom";
  const std::size_t profile_location=invalid_profile_text.find(valid_profile);
  if (profile_location==std::string::npos) {
    std::fprintf(stderr,"FAIL invalid_profile_fixture missing key\n");
    return 1;
  }
  invalid_profile_text.replace(profile_location,valid_profile.size(),
                               "mesh_profile = production");
  {
    std::ofstream invalid_profile(invalid_profile_file);
    invalid_profile << invalid_profile_text;
  }
  if (Moon::Runtime::ParseApplicationInput(invalid_profile_file.string(),
                                            &invalid_configuration,&error)) {
    std::fprintf(stderr,"FAIL invalid mesh_profile was accepted\n");
    return 1;
  }

  Moon::Surface::LolaDem dem;
  if (!dem.Load(configuration,&error)) {
    std::fprintf(stderr,"FAIL lola_load %s\n",error.c_str());
    return 1;
  }

  int failures=0;
  // The agreed M1 grid has 41 near-surface points (including zero), 32
  // additional middle-altitude points and 8 logarithmic points from above
  // 1000 km through 10000 km.  These checks exercise the production parser's
  // endpoint rules without treating configuration plumbing as an observation-
  // operator PASS.
  failures+=report("profile_altitude_grid_points",
      configuration.profileAltitudeGridM.size(),81.0,0.0);
  failures+=report("profile_altitude_grid_first_m",
      configuration.profileAltitudeGridM.front(),0.0,0.0);
  failures+=report("profile_altitude_grid_near_boundary_m",
      configuration.profileAltitudeGridM[40],200000.0,0.0);
  failures+=report("profile_altitude_grid_middle_boundary_m",
      configuration.profileAltitudeGridM[72],1000000.0,0.0);
  failures+=report("profile_altitude_grid_last_m",
      configuration.profileAltitudeGridM.back(),10000000.0,1.0e-14);
  failures+=report("profile_sza_count",
      configuration.profileSolarZenithAnglesDeg.size(),6.0,0.0);
  const double expected_sza_deg[6]={0.0,30.0,60.0,90.0,120.0,180.0};
  for (int i=0;i<6;i++) {
    const std::string metric="profile_sza_"+std::to_string(i)+"_deg";
    failures+=report(metric.c_str(),
        configuration.profileSolarZenithAnglesDeg[i],expected_sza_deg[i],0.0);
  }
  failures+=report("profile_include_body_fixed_poles",
      configuration.profileIncludeBodyFixedPoles ? 1.0 : 0.0,1.0,0.0);
  failures+=report("profile_include_shell_average",
      configuration.profileIncludeShellAverage ? 1.0 : 0.0,1.0,0.0);
  failures+=report_string("profile_shell_angular_weighting",
      configuration.profileShellAngularWeighting,"solid-angle");
  failures+=report_string("simulation_start_utc",
      configuration.simulationStartUtc,"2009-01-24T00:00:00");
  failures+=report_string("simulation_environment_mode",
      configuration.simulationEnvironmentMode,"frozen");
  failures+=report("random_seed",configuration.randomSeed,731993.0,0.0);
  failures+=report_string("mesh_profile",configuration.meshProfile,"custom");
  failures+=report("mesh_far_field_extent_m",
      configuration.meshFarFieldExtentM,10000000.0,0.0);
  failures+=report("mesh_outer_cell_size_m",
      configuration.meshOuterCellSizeM,17371000.0,0.0);
  failures+=report("mesh_output_tecplot",
      configuration.meshOutputTecplot ? 1.0 : 0.0,0.0,0.0);
  failures+=report_string("limb_observer",
      configuration.limbObserver,"earth");
  failures+=report_string("limb_plane",
      configuration.limbPlane,"sun-moon-earth");
  failures+=report_string("limb_sides",configuration.limbSides,"both");
  failures+=report_string("limb_branch_reference",
      configuration.limbBranchReference,"projected-sun");
  failures+=report("limb_tangent_maximum_altitude_m",
      configuration.limbTangentMaximumAltitudeM,1000000.0,0.0);
  failures+=report_string("tail_center_direction",
      configuration.tailCenterDirection,"anti-lunar");
  failures+=report_string("tail_map_coordinate_system",
      configuration.tailMapCoordinateSystem,"j2000-ra-dec");
  failures+=report_string("tail_map_display_orientation",
      configuration.tailMapDisplayOrientation,"west-right");
  failures+=report("tail_map_half_width_deg",
      configuration.tailMapHalfWidthDeg,7.5,0.0);
  failures+=report("tail_beam_diameter_deg",
      configuration.tailBeamDiameterDeg,1.0,0.0);
  failures+=report("tail_exposure_duration_s",
      configuration.tailExposureDurationS,120.0,0.0);
  failures+=report_string("tail_beam_response",
      configuration.tailBeamResponse,"unresolved");
  failures+=report_string("tail_temporal_operator",
      configuration.tailTemporalOperator,"exposure-average");
  failures+=report_string("alternate_limb_sides",
      alternate_configuration.limbSides,"away-from-projected-sun");
  failures+=report_string("alternate_tail_map_coordinate_system",
      alternate_configuration.tailMapCoordinateSystem,
      "eclipj2000-lon-lat");
  failures+=report_string("alternate_tail_map_display_orientation",
      alternate_configuration.tailMapDisplayOrientation,"east-right");
  // Exercise the production epoch override API used by amps_init().  CSPICE
  // UTC-to-ET interpretation remains a linked orbit-enabled responsibility;
  // this local check verifies that the immutable parser value reaches the
  // exosphere module without a dangling c_str() pointer or truncation.
  Exosphere::SetSimulationStartTimeString(
      configuration.simulationStartUtc.c_str());
  failures+=report_string("simulation_start_override",
      Exosphere::GetSimulationStartTimeString(),"2009-01-24T00:00:00");
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

int lola_surface_demo(int argument_count,char** arguments) {
  // This is a thin presentation adapter around the same parser, D01 decoder,
  // icosphere builder, topographic displacement, and production writers used
  // by U08 and normal srcMoon initialization.  It deliberately supplies no
  // independent oracle and therefore cannot return an I14 validation result;
  // the Python demonstration runner performs the independent D01 comparison
  // and records the formal V&V status separately.
  if (argument_count!=2) {
    std::fprintf(stderr,
        "usage: production_kernel_probe lola-surface-demo "
        "<application.in> <vertex-output.csv>\n");
    return 2;
  }
  const std::string input_path(arguments[0]);
  const std::string vertex_path(arguments[1]);
  Moon::Runtime::Configuration configuration;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(input_path,&configuration,&error)) {
    std::fprintf(stderr,"ERROR lola_demo_input %s\n",error.c_str());
    return 2;
  }
  if (configuration.surfaceGeometry!=Moon::Runtime::SurfaceGeometry::Lola) {
    std::fprintf(stderr,"ERROR lola_demo_requires_surface_geometry_lola\n");
    return 2;
  }

  Moon::Surface::LolaDem dem;
  Moon::Surface::Triangulation mesh;
  if (!dem.Load(configuration,&error) ||
      !Moon::Surface::BuildIcosphere(dem.metadata().referenceRadiusM,
          configuration.surfaceMeshResolutionM,&mesh,&error) ||
      !Moon::Surface::ApplyLolaTopography(dem,&mesh,&error) ||
      !Moon::Surface::WriteCeaSurface(
          mesh,configuration.surfaceCeaFile,&error) ||
      !Moon::Surface::WriteTecplotSurface(
          mesh,configuration.surfaceTecplotFile,&error)) {
    std::fprintf(stderr,"ERROR lola_demo_production_kernel %s\n",error.c_str());
    return 2;
  }

  // Preserve the complete vertex population rather than selecting attractive
  // locations for a meeting figure.  Coordinates remain label-native
  // planetocentric latitude and east-positive longitude, while radii and
  // elevations are SI metres relative to the LDEM_4 reference sphere.
  std::ofstream vertices(vertex_path);
  if (!vertices) {
    std::fprintf(stderr,"ERROR cannot_write_demo_vertices %s\n",
                 vertex_path.c_str());
    return 2;
  }
  vertices << std::setprecision(17)
           << "vertex_index,x_m,y_m,z_m,latitude_deg,longitude_east_deg,"
              "radius_m,model_elevation_m\n";
  std::size_t inward_faces=0,degenerate_faces=0;
  const double pi=std::acos(-1.0);
  for (std::size_t index=0;index<mesh.verticesM.size();++index) {
    const auto& vertex=mesh.verticesM[index];
    const double radius=std::sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+
                                  vertex[2]*vertex[2]);
    const double latitude=std::asin(vertex[2]/radius)*180.0/pi;
    double longitude=std::atan2(vertex[1],vertex[0])*180.0/pi;
    if (longitude<0.0) longitude+=360.0;
    vertices << index << ',' << vertex[0] << ',' << vertex[1] << ','
             << vertex[2] << ',' << latitude << ',' << longitude << ','
             << radius << ','
             << radius-dem.metadata().referenceRadiusM << '\n';
  }
  for (const auto& face:mesh.faces) {
    const auto& a=mesh.verticesM[face[0]];
    const auto& b=mesh.verticesM[face[1]];
    const auto& c=mesh.verticesM[face[2]];
    const double ab[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
    const double ac[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
    const double normal[3]={ab[1]*ac[2]-ab[2]*ac[1],
        ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]};
    const double area2=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+
                                 normal[2]*normal[2]);
    const double center[3]={a[0]+b[0]+c[0],a[1]+b[1]+c[1],
                            a[2]+b[2]+c[2]};
    if (!(area2>0.0) || !std::isfinite(area2)) ++degenerate_faces;
    if (normal[0]*center[0]+normal[1]*center[1]+normal[2]*center[2]<=0.0)
      ++inward_faces;
  }
  if (!vertices) {
    std::fprintf(stderr,"ERROR incomplete_demo_vertex_output\n");
    return 2;
  }
  std::printf("surface_resolution_m=%.17g\n",
              configuration.surfaceMeshResolutionM);
  std::printf("reference_radius_m=%.17g\n",dem.metadata().referenceRadiusM);
  std::printf("subdivision_level=%u\n",mesh.subdivisionLevel);
  std::printf("vertex_count=%zu\n",mesh.verticesM.size());
  std::printf("face_count=%zu\n",mesh.faces.size());
  std::printf("maximum_edge_length_m=%.17g\n",mesh.maximumEdgeLengthM);
  std::printf("minimum_radius_m=%.17g\n",mesh.minimumRadiusM);
  std::printf("maximum_radius_m=%.17g\n",mesh.maximumRadiusM);
  std::printf("degenerate_face_count=%zu\n",degenerate_faces);
  std::printf("inward_face_count=%zu\n",inward_faces);
  return (degenerate_faces==0 && inward_faces==0) ? 0 : 1;
}

int surface_environment_demo(int argument_count,char** arguments) {
  // Export production illumination and thermal state without duplicating
  // either physical kernel in Python.  The adapter constructs the same D01
  // terrain/BVH and logical 60x100 surface used by normal srcMoon, initializes
  // the normal thermal columns, and asks the read-only production diagnostic
  // service for every cell.  It deliberately performs no observational score.
  if (argument_count!=5) {
    std::fprintf(stderr,
        "usage: production_kernel_probe surface-environment-demo "
        "<application.in> <surface-output.csv> <profile-output.csv> "
        "<snapshot-utc> <input-initial|periodic-spinup>\n");
    return 2;
  }
#ifdef _NO_SPICE_CALLS_
  std::fprintf(stderr,
      "ERROR surface_environment_demo_requires_spice_enabled_build\n");
  return 2;
#else
  const std::string input_path(arguments[0]);
  const std::string surface_path(arguments[1]);
  const std::string profile_path(arguments[2]);
  const std::string snapshot_utc(arguments[3]);
  const std::string initialization(arguments[4]);
  if (initialization!="input-initial" && initialization!="periodic-spinup") {
    std::fprintf(stderr,"ERROR unsupported_thermal_initialization %s\n",
                 initialization.c_str());
    return 2;
  }

  Moon::Runtime::Configuration configuration;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(input_path,&configuration,&error)) {
    std::fprintf(stderr,"ERROR environment_demo_input %s\n",error.c_str());
    return 2;
  }
  if (configuration.surfaceGeometry!=Moon::Runtime::SurfaceGeometry::Lola ||
      !configuration.m2EnvironmentEnabled ||
      configuration.surfaceTemperatureMode!="thermal-model") {
    std::fprintf(stderr,
        "ERROR environment_demo_requires_lola_m2_thermal_model\n");
    return 2;
  }

  Moon::Surface::LolaDem dem;
  Moon::Surface::Triangulation mesh;
  if (!dem.Load(configuration,&error) ||
      !Moon::Surface::BuildIcosphere(dem.metadata().referenceRadiusM,
          configuration.surfaceMeshResolutionM,&mesh,&error) ||
      !Moon::Surface::ApplyLolaTopography(dem,&mesh,&error) ||
      !Moon::SurfaceEnvironment::InitializeProductionTerrain(
          mesh,configuration,&error) ||
      !Moon::Surface::WriteCeaSurface(
          mesh,configuration.surfaceCeaFile,&error) ||
      !Moon::Surface::WriteTecplotSurface(
          mesh,configuration.surfaceTecplotFile,&error)) {
    std::fprintf(stderr,"ERROR environment_demo_surface %s\n",error.c_str());
    return 2;
  }

  // This logical equal-area grid is the production inventory/thermal grid,
  // not the D01 triangle population.  Each radial logical cell is intersected
  // with the terrain BVH before its actual position, normal and horizon are
  // evaluated.  The fixed 60x100 dimensions therefore reproduce normal
  // srcMoon behavior rather than inventing a presentation-only raster.
  cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);
  auto* logical_surface=new cInternalSphericalData;
  double origin_m[3]={0.0,0.0,0.0};
  logical_surface->SetSphereGeometricalParameters(origin_m,_RADIUS_(_MOON_));
  const int element_count=logical_surface->GetTotalSurfaceElementsNumber();
  logical_surface->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,
      element_count,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,logical_surface);

  // Init_SPICE furnishes the exact kernels compiled by Config.pl.  Applying
  // the input UTC before initialization mirrors amps_init(); no identity frame
  // or approximate lunar rotation is permitted as a fallback.
  Exosphere::SetSimulationStartTimeString(
      configuration.simulationStartUtc.c_str());
  Exosphere::Init_SPICE();
  SpiceDouble initial_et=0.0,snapshot_et=0.0;
  str2et_c(configuration.simulationStartUtc.c_str(),&initial_et);
  str2et_c(snapshot_utc.c_str(),&snapshot_et);
  if (!std::isfinite(initial_et) || !std::isfinite(snapshot_et)) {
    std::fprintf(stderr,"ERROR non_finite_demo_epoch\n");
    return 2;
  }
  if (!Moon::SurfaceTemperature::InitializeProductionState(
          configuration,logical_surface,&error)) {
    std::fprintf(stderr,"ERROR environment_demo_thermal_init %s\n",
                 error.c_str());
    return 2;
  }

  Moon::SurfaceTemperature::SpinUpMetrics spinup;
  double state_epoch=initial_et;
  if (initialization=="periodic-spinup") {
    if (!Moon::SurfaceTemperature::SpinUpProductionState(&spinup,&error)) {
      std::fprintf(stderr,"ERROR environment_demo_spinup %s\n",error.c_str());
      return 2;
    }
    state_epoch=spinup.cycleEndEphemerisTimeS;
  }
  if (snapshot_et<state_epoch) {
    std::fprintf(stderr,
        "ERROR snapshot_precedes_thermal_state snapshot_et=%.17g "
        "state_et=%.17g; choose a later UTC or input-initial mode\n",
        snapshot_et,state_epoch);
    return 2;
  }
  if (!Moon::SurfaceTemperature::AdvanceProductionState(
          state_epoch,snapshot_et,&error)) {
    std::fprintf(stderr,"ERROR environment_demo_thermal_advance %s\n",
                 error.c_str());
    return 2;
  }
  std::vector<Moon::SurfaceTemperature::SurfaceEnvironmentSample> samples;
  if (!Moon::SurfaceTemperature::SampleProductionSurfaceEnvironment(
          snapshot_et,&samples,&error)) {
    std::fprintf(stderr,"ERROR environment_demo_snapshot %s\n",error.c_str());
    return 2;
  }

  std::ofstream surface_output(surface_path);
  if (!surface_output) {
    std::fprintf(stderr,"ERROR cannot_write_environment_csv %s\n",
                 surface_path.c_str());
    return 2;
  }
  surface_output << std::setprecision(17)
      << "surface_element,zenith_index,azimuth_index,latitude_deg,"
         "longitude_east_deg,local_solar_time_hour,x_body_fixed_m,"
         "y_body_fixed_m,z_body_fixed_m,radius_m,elevation_m,"
         "normal_x_body_fixed,normal_y_body_fixed,normal_z_body_fixed,"
         "solar_incidence_cosine,sun_above_local_horizon,"
         "terrain_sun_visible,terrain_occluded,"
         "solar_disc_visible_fraction,earth_eclipse_fraction,"
         "direct_beam_fraction,absorbed_solar_flux_w_m2,"
         "surface_temperature_k,last_energy_residual_j_m2,"
         "last_nonlinear_iterations\n";
  for (const auto& sample:samples) {
    const double radius=std::sqrt(
        sample.surfacePositionBodyFixedM[0]*
            sample.surfacePositionBodyFixedM[0]+
        sample.surfacePositionBodyFixedM[1]*
            sample.surfacePositionBodyFixedM[1]+
        sample.surfacePositionBodyFixedM[2]*
            sample.surfacePositionBodyFixedM[2]);
    surface_output << sample.surfaceElement << ',' << sample.zenithIndex << ','
        << sample.azimuthIndex << ',' << sample.latitudeDeg << ','
        << sample.longitudeEastDeg << ',' << sample.localSolarTimeHour << ','
        << sample.surfacePositionBodyFixedM[0] << ','
        << sample.surfacePositionBodyFixedM[1] << ','
        << sample.surfacePositionBodyFixedM[2] << ',' << radius << ','
        << radius-dem.metadata().referenceRadiusM << ','
        << sample.terrainNormalBodyFixed[0] << ','
        << sample.terrainNormalBodyFixed[1] << ','
        << sample.terrainNormalBodyFixed[2] << ','
        << sample.solarIncidenceCosine << ','
        << (sample.sunAboveLocalHorizon?1:0) << ','
        << (sample.terrainSunVisible?1:0) << ','
        << (sample.terrainOccluded?1:0) << ','
        << sample.solarDiscVisibleFraction << ','
        << 1.0-sample.solarDiscVisibleFraction << ','
        << sample.directBeamFraction << ','
        << sample.absorbedSolarFluxWm2 << ','
        << sample.surfaceTemperatureK << ','
        << sample.lastEnergyResidualJm2 << ','
        << sample.lastNonlinearIterations << '\n';
  }
  surface_output.close();
  if (!surface_output) {
    std::fprintf(stderr,"ERROR incomplete_environment_csv\n");
    return 2;
  }

  // Preserve every layer of every logical column.  Selecting only visually
  // interesting sites here would make the output irreproducible and could
  // hide a deep-column initialization problem.
  std::ofstream profile_output(profile_path);
  if (!profile_output) {
    std::fprintf(stderr,"ERROR cannot_write_thermal_profile_csv %s\n",
                 profile_path.c_str());
    return 2;
  }
  profile_output << std::setprecision(17)
      << "surface_element,latitude_deg,longitude_east_deg,layer_index,"
         "layer_center_depth_m,layer_temperature_k\n";
  std::size_t profile_rows=0;
  for (const auto& sample:samples) {
    std::vector<double> depth_m,temperature_k;
    if (!Moon::SurfaceTemperature::ProductionThermalProfile(
            sample.surfaceElement,&depth_m,&temperature_k,&error)) {
      std::fprintf(stderr,"ERROR environment_demo_profile %s\n",error.c_str());
      return 2;
    }
    for (std::size_t layer=0;layer<depth_m.size();++layer) {
      profile_output << sample.surfaceElement << ',' << sample.latitudeDeg
          << ',' << sample.longitudeEastDeg << ',' << layer << ','
          << depth_m[layer] << ',' << temperature_k[layer] << '\n';
      ++profile_rows;
    }
  }
  profile_output.close();
  if (!profile_output) {
    std::fprintf(stderr,"ERROR incomplete_thermal_profile_csv\n");
    return 2;
  }

  const auto count_true=[&samples](auto member) {
    return std::count_if(samples.begin(),samples.end(),
        [member](const auto& sample) { return sample.*member; });
  };
  const auto temperature_range=std::minmax_element(samples.begin(),samples.end(),
      [](const auto& first,const auto& second) {
        return first.surfaceTemperatureK<second.surfaceTemperatureK;
      });
  std::printf("surface_element_count=%zu\n",samples.size());
  std::printf("thermal_profile_row_count=%zu\n",profile_rows);
  std::printf("terrain_occluded_count=%zu\n",count_true(
      &Moon::SurfaceTemperature::SurfaceEnvironmentSample::terrainOccluded));
  std::printf("terrain_visible_count=%zu\n",count_true(
      &Moon::SurfaceTemperature::SurfaceEnvironmentSample::terrainSunVisible));
  std::printf("minimum_surface_temperature_k=%.17g\n",
      temperature_range.first->surfaceTemperatureK);
  std::printf("maximum_surface_temperature_k=%.17g\n",
      temperature_range.second->surfaceTemperatureK);
  std::printf("thermal_initialization=%s\n",initialization.c_str());
  std::printf("snapshot_et_s=%.17g\n",snapshot_et);
  std::printf("advanced_duration_s=%.17g\n",snapshot_et-state_epoch);
  std::printf("spinup_completed_cycles=%lld\n",spinup.completedCycles);
  std::printf("spinup_maximum_node_change_k=%.17g\n",
      spinup.maximumNodeTemperatureChangeK);
  std::printf("spinup_maximum_enthalpy_change_j_m2=%.17g\n",
      spinup.maximumColumnEnthalpyChangeJm2);
  return 0;
#endif
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
        "usage: production_kernel_probe <gravity|mesh-resolution|"
        "photochemistry|radiation-shadow|radiation-pressure-values|"
        "rotating-frame|sodium-sources|sticking|temperature|lola-geometry|"
        "lola-surface-demo|surface-environment-demo|"
        "terrain-shadow|diviner-temperature|thermal-inertia|surface-physics|"
        "cold-trapping|argon-source> "
        "[mode arguments]\n");
    return 2;
  }

  const std::string test(argv[1]);
  if (test == "radiation-pressure-values") {
    return radiation_pressure_values(argc - 2, argv + 2);
  }
  if (test == "lola-surface-demo") {
    return lola_surface_demo(argc - 2,argv + 2);
  }
  if (test == "surface-environment-demo") {
    return surface_environment_demo(argc - 2,argv + 2);
  }
  if (argc != 2) {
    std::fprintf(stderr, "mode %s does not accept additional arguments\n",
        argv[1]);
    return 2;
  }
  if (test == "gravity") return gravity();
  if (test == "mesh-resolution") return mesh_resolution();
  if (test == "lola-geometry") return lola_geometry();
  if (test == "photochemistry") return photochemistry();
  if (test == "radiation-shadow") return radiation_shadow();
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  if (test == "rotating-frame") return rotating_frame();
#endif
  if (test == "sodium-sources") return sodium_sources();
  if (test == "sticking") return sticking();
  if (test == "temperature") return temperature();
  if (test == "terrain-shadow") return terrain_shadow();
  if (test == "diviner-temperature") return diviner_temperature();
  if (test == "thermal-inertia") return thermal_inertia();
  if (test == "surface-physics") return surface_physics();
  if (test == "cold-trapping") return cold_trapping();
  if (test == "argon-source") return argon_source();

  std::fprintf(stderr, "unknown test: %s\n", argv[1]);
  return 2;
}
