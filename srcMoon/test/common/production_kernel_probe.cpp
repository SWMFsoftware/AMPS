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

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
        "usage: production_kernel_probe <gravity|mesh-resolution|"
        "photochemistry|radiation-shadow|radiation-pressure-values|"
        "rotating-frame|sodium-sources|sticking|temperature|lola-geometry> "
        "[mode arguments]\n");
    return 2;
  }

  const std::string test(argv[1]);
  if (test == "radiation-pressure-values") {
    return radiation_pressure_values(argc - 2, argv + 2);
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

  std::fprintf(stderr, "unknown test: %s\n", argv[1]);
  return 2;
}
