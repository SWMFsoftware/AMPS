#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

// PST07 is an integration-equivalence test, so copied adapter fields are
// required to be exactly equal to the direct production value.  A relative
// tolerance is used only for an independently recomputed spectrum expression,
// where a different but algebraically equivalent operation order is expected.
void expect_exact(swcme_test::Context& context,const std::string& label,
                  double adapter_value,double direct_value) {
  context.expect_true(adapter_value==direct_value,
                      label+" is copied without conversion or recomputation");
}

void expect_relative(swcme_test::Context& context,const std::string& label,
                     double actual,double expected,double tolerance=2.0e-13) {
  const double scale=std::max({1.0e-300,std::abs(actual),std::abs(expected)});
  context.expect_near(actual,expected,tolerance*scale,label);
}

std::array<double,3> normalized(const std::array<double,3>& value) {
  const double magnitude=std::hypot(value[0],std::hypot(value[1],value[2]));
  return {{value[0]/magnitude,value[1]/magnitude,value[2]/magnitude}};
}

// Configure a fast spherical shock and an explicitly normalized proton source.
// The same physical values are used in both dimensions; only the 3-D geometry
// controls are additional.  This prevents parameter retuning from concealing
// an adapter conversion or default mismatch.
void make_fixture(swcme1d::Params& one,swcme3d::Params& three,
                  swcme::sep::SpectrumConfig& spectrum) {
  one.V_sw_kms=415.0;
  one.n1AU_cm3=6.25;
  one.B1AU_nT=5.75;
  one.T_K=1.4e5;
  one.gamma_ad=5.0/3.0;
  one.sin_theta=1.0;
  one.kinematics_mode=swcme::kinematics::Mode::DBM;
  one.r0_Rs=20.0;
  one.V0_sh_kms=1500.0;
  one.Gamma_kmInv=4.0e-8;
  one.region_mode=swcme::regions::Mode::ShockOnly;
  one.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  one.relative_source_weight_per_area=1.75;

  three.V_sw_kms=one.V_sw_kms;
  three.n1AU_cm3=one.n1AU_cm3;
  three.B1AU_nT=one.B1AU_nT;
  three.T_K=one.T_K;
  three.gamma_ad=one.gamma_ad;
  three.sin_theta=one.sin_theta;
  three.kinematics_mode=one.kinematics_mode;
  three.r0_Rs=one.r0_Rs;
  three.V0_sh_kms=one.V0_sh_kms;
  three.Gamma_kmInv=one.Gamma_kmInv;
  three.region_mode=one.region_mode;
  three.shock_acceleration_mode=one.shock_acceleration_mode;
  three.relative_source_weight_per_area=one.relative_source_weight_per_area;
  three.shape=swcme3d::ShockShape::Sphere;
  three.cme_dir[0]=1.0; three.cme_dir[1]=0.0; three.cme_dir[2]=0.0;
  three.solar_rotation_axis[0]=0.0;
  three.solar_rotation_axis[1]=0.0;
  three.solar_rotation_axis[2]=1.0;

  spectrum.kinetic_energy_min_MeV=1.0;
  spectrum.kinetic_energy_max_MeV=1000.0;
  spectrum.reference_energy_MeV=10.0;
  spectrum.normalization=
      swcme::sep::NormalizationMode::ReferenceDifferentialIntensity;
  spectrum.reference_differential_intensity_SI=
      swcme::sep::differential_intensity_common_to_SI(12.5);
}

// Compare the common direct ShockAccelerationState with the stable
// AMPS-facing SEPSourceState.  This explicitly covers every shared scalar and
// vector rather than relying on the serialized record, which could hide an
// equal formatting error on both sides.
void compare_source_mapping(
    swcme_test::Context& context,const std::string& prefix,
    const swcme::acceleration::ShockAccelerationState& direct,
    double direct_pressure_Pa,double direct_focusing_m,
    const swcme::sep::SEPSourceState& adapter) {
  context.expect_true(adapter.status.ok(),prefix+" adapter source status is OK");
  context.expect_true(adapter.active==direct.source_enabled &&
                          adapter.acceleration_mode==direct.mode,
                      prefix+" activation and acceleration mode are identical");
  expect_exact(context,prefix+" time",adapter.time_s,direct.time_s);
  for (std::size_t component=0;component<3;++component) {
    expect_exact(context,prefix+" source position",adapter.position_m[component],
                 direct.position_m[component]);
    expect_exact(context,prefix+" shock normal",adapter.normal[component],
                 direct.normal[component]);
  }
  expect_exact(context,prefix+" compression",adapter.compression,direct.compression);
  expect_exact(context,prefix+" theta_Bn",adapter.theta_Bn_rad,direct.theta_Bn_rad);
  expect_exact(context,prefix+" fast Mach",adapter.fast_mach,direct.fast_mach);
  expect_exact(context,prefix+" normal shock speed",adapter.normal_speed_m_s,
               direct.normal_speed_m_s);
  expect_exact(context,prefix+" upstream density",adapter.upstream_density_m3,
               direct.upstream_density_m3);
  expect_exact(context,prefix+" upstream magnetic magnitude",adapter.upstream_B_T,
               direct.upstream_B_T);
  expect_exact(context,prefix+" upstream pressure",adapter.upstream_pressure_Pa,
               direct_pressure_Pa);
  expect_exact(context,prefix+" focusing length",adapter.focusing_length_m,
               direct_focusing_m);
  expect_exact(context,prefix+" DSA phase-space slope",adapter.q_phase_space,
               direct.dsa_q_phase_space);
  expect_exact(context,prefix+" momentum-intensity index",
               adapter.momentum_intensity_index,
               direct.dsa_q_phase_space-2.0);
  expect_exact(context,prefix+" nonrelativistic energy-intensity index",
               adapter.nonrel_energy_intensity_index,
               0.5*(direct.dsa_q_phase_space-2.0));
  expect_exact(context,prefix+" source weight control",
               adapter.relative_source_weight_per_area,
               direct.relative_source_weight_per_area);
  expect_exact(context,prefix+" point-source relative weight",
               adapter.relative_patch_weight,
               direct.relative_source_weight_per_area);
}

// Recompute J(E)/J(Eref) from the direct shock slope and relativistic momentum.
// The adapter is then exercised at the exact same energy nodes.  This detects
// energy-unit, reference-energy, exponent-sign, and normalization drift.
void compare_energy_grid(swcme_test::Context& context,const std::string& prefix,
                         const swcme::acceleration::ShockAccelerationState& direct,
                         const swcme::sep::SEPSourceState& adapter) {
  const std::vector<double> energy_MeV={1.0,3.0,10.0,100.0,1000.0};
  const double reference_momentum=
      swcme::sep::momentum_kg_m_s_from_kinetic_MeV(
          adapter.spectrum.reference_energy_MeV,
          adapter.spectrum.particle_mass_kg);
  for (double energy : energy_MeV) {
    const double momentum=swcme::sep::momentum_kg_m_s_from_kinetic_MeV(
        energy,adapter.spectrum.particle_mass_kg);
    const double expected_shape=std::pow(
        momentum/reference_momentum,2.0-direct.dsa_q_phase_space);
    double shape=std::numeric_limits<double>::quiet_NaN();
    const swcme::ModelStatus shape_status=
        swcme::sep::relative_intensity_shape(adapter,energy,shape);
    context.expect_true(shape_status.ok(),prefix+" spectrum node status is OK");
    expect_relative(context,prefix+" relative spectrum node",shape,expected_shape);

    double intensity=std::numeric_limits<double>::quiet_NaN();
    const swcme::ModelStatus intensity_status=
        swcme::sep::differential_intensity_SI(adapter,energy,intensity);
    context.expect_true(intensity_status.ok(),prefix+" SI intensity status is OK");
    expect_relative(context,prefix+" SI spectrum normalization",intensity,
                    adapter.spectrum.reference_differential_intensity_SI*
                        expected_shape);
  }
}

void compare_background_1d(swcme_test::Context& context,
                           const swcme::sep::Interface1D& adapter,
                           const swcme1d::StepState& step,double radius_m) {
  swcme::sep::BackgroundState wrapped;
  const swcme::ModelStatus adapter_status=
      adapter.evaluate_background(step,radius_m,wrapped);
  double density=0.0,velocity=0.0,Br=0.0,Bphi=0.0,Bmag=0.0,divV=0.0;
  const swcme::ModelStatus direct_status=
      adapter.model().evaluate_radii_with_B_div_checked(
          step,&radius_m,&density,&velocity,&Br,&Bphi,&Bmag,&divV,1);
  context.expect_true(adapter_status.code==direct_status.code &&
                          adapter_status.ok(),
                      "1-D direct and adapter background statuses are identical");
  context.expect_true(wrapped.owner_model_identity==step.owner_model_identity &&
                          wrapped.owner_model_identity==adapter.model().model_identity(),
                      "1-D adapter preserves the model identity");
  context.expect_true(wrapped.configuration_digest==step.configuration_digest,
                      "1-D adapter preserves the configuration digest");
  expect_exact(context,"1-D reported radius",wrapped.position_m[0],radius_m);
  expect_exact(context,"1-D reported y",wrapped.position_m[1],0.0);
  expect_exact(context,"1-D reported z",wrapped.position_m[2],0.0);
  expect_exact(context,"1-D density",wrapped.density_m3,density);
  expect_exact(context,"1-D pressure",wrapped.pressure_Pa,
               swcme::solarwind::proton_pressure_Pa(
                   step.common.solar_wind,density));
  expect_exact(context,"1-D radial velocity",wrapped.velocity_m_s[0],velocity);
  expect_exact(context,"1-D transverse velocity y",wrapped.velocity_m_s[1],0.0);
  expect_exact(context,"1-D transverse velocity z",wrapped.velocity_m_s[2],0.0);
  expect_exact(context,"1-D Br",wrapped.magnetic_T[0],Br);
  expect_exact(context,"1-D Bphi",wrapped.magnetic_T[1],Bphi);
  expect_exact(context,"1-D |B|",wrapped.magnetic_magnitude_T,Bmag);
  expect_exact(context,"1-D div(V)",wrapped.div_velocity_s_inv,divV);
  expect_exact(context,"1-D focusing",wrapped.focusing_length_m,
               swcme::solarwind::parker_focusing_length_m(
                   step.common.solar_wind,radius_m,
                   step.common.solar_wind.reference_sin_theta));
}

void compare_background_3d(swcme_test::Context& context,
                           const swcme::sep::Interface3D& adapter,
                           const swcme3d::StepState& step,
                           const std::array<double,3>& position_m) {
  swcme::sep::BackgroundState wrapped;
  const swcme::ModelStatus adapter_status=
      adapter.evaluate_background(step,position_m,wrapped);
  const double x=position_m[0],y=position_m[1],z=position_m[2];
  double density=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0,divV=0.0;
  const swcme::ModelStatus direct_status=
      adapter.model().evaluate_cartesian_with_B_div_checked(
          step,&x,&y,&z,&density,&vx,&vy,&vz,&bx,&by,&bz,&divV,1);
  context.expect_true(adapter_status.code==direct_status.code &&
                          adapter_status.ok(),
                      "3-D direct and adapter background statuses are identical");
  context.expect_true(wrapped.owner_model_identity==step.owner_model_identity &&
                          wrapped.owner_model_identity==adapter.model().model_identity(),
                      "3-D adapter preserves the model identity");
  context.expect_true(wrapped.configuration_digest==step.configuration_digest,
                      "3-D adapter preserves the configuration digest");
  for (std::size_t component=0;component<3;++component)
    expect_exact(context,"3-D reported position",wrapped.position_m[component],
                 position_m[component]);
  expect_exact(context,"3-D density",wrapped.density_m3,density);
  expect_exact(context,"3-D pressure",wrapped.pressure_Pa,
               swcme::solarwind::proton_pressure_Pa(
                   step.common.solar_wind,density));
  const std::array<double,3> direct_velocity={{vx,vy,vz}};
  const std::array<double,3> direct_magnetic={{bx,by,bz}};
  for (std::size_t component=0;component<3;++component) {
    expect_exact(context,"3-D velocity component",wrapped.velocity_m_s[component],
                 direct_velocity[component]);
    expect_exact(context,"3-D magnetic component",wrapped.magnetic_T[component],
                 direct_magnetic[component]);
  }
  expect_exact(context,"3-D |B|",wrapped.magnetic_magnitude_T,
               std::hypot(bx,std::hypot(by,bz)));
  expect_exact(context,"3-D div(V)",wrapped.div_velocity_s_inv,divV);

  const double radius_m=std::hypot(x,std::hypot(y,z));
  const std::array<double,3> radial_hat={{x/radius_m,y/radius_m,z/radius_m}};
  const std::array<double,3> axis={{step.solar_axis_hat[0],
                                    step.solar_axis_hat[1],
                                    step.solar_axis_hat[2]}};
  expect_exact(context,"3-D focusing",wrapped.focusing_length_m,
               swcme::solarwind::parker_focusing_length_m(
                   step.common.solar_wind,radius_m,
                   swcme::solarwind::parker_sin_colatitude(axis,radial_hat)));
}

}  // namespace

void test_pst07(swcme_test::Context& context) {
  std::cout << "PST07 AMPS adapter equivalence\n";

  swcme1d::Params one_params;
  swcme3d::Params three_params;
  swcme::sep::SpectrumConfig spectrum;
  make_fixture(one_params,three_params,spectrum);
  swcme::sep::Interface1D one(one_params,spectrum);
  swcme::sep::Interface3D three(three_params,spectrum);

  // Multiple upstream radii catch accidental distance-unit conversion and
  // verify that pressure/focusing use the same prepared solar-wind state as
  // density and B rather than adapter-owned defaults.
  const auto one_step=one.prepare(2.0*3600.0);
  const auto three_step=three.prepare(2.0*3600.0);
  for (double radius_AU : {0.35,0.8,1.0,1.7}) {
    const double radius_m=radius_AU*swcme::constants::AU_M;
    compare_background_1d(context,one,one_step,radius_m);
    const std::array<double,3> direction=normalized({{0.82,0.43,0.17}});
    compare_background_3d(
        context,three,three_step,
        {{radius_m*direction[0],radius_m*direction[1],radius_m*direction[2]}});
  }

  // Failure status is part of the adapter contract too.  The same unsupported
  // radius must be classified identically by the direct and AMPS-facing APIs;
  // a wrapper may not turn it into a finite fallback background.
  const double invalid_radius=swcme::solarwind::MIN_RADIUS_M*0.99;
  double n=17.0,v=18.0;
  const swcme::ModelStatus direct_invalid_1d=
      one.model().evaluate_radii_fast_checked(
          one_step,&invalid_radius,&n,&v,1);
  swcme::sep::BackgroundState invalid_background_1d;
  const swcme::ModelStatus adapter_invalid_1d=
      one.evaluate_background(one_step,invalid_radius,invalid_background_1d);
  context.expect_true(direct_invalid_1d.code==adapter_invalid_1d.code &&
                          adapter_invalid_1d.code==
                              swcme::StatusCode::OutsideModelDomain,
                      "1-D failure status is unchanged by the adapter");

  const std::array<double,3> invalid_position={{invalid_radius,0.0,0.0}};
  const double ix=invalid_position[0],iy=0.0,iz=0.0;
  double nv=1.0,vx=2.0,vy=3.0,vz=4.0;
  const swcme::ModelStatus direct_invalid_3d=
      three.model().evaluate_cartesian_fast_checked(
          three_step,&ix,&iy,&iz,&nv,&vx,&vy,&vz,1);
  swcme::sep::BackgroundState invalid_background_3d;
  const swcme::ModelStatus adapter_invalid_3d=
      three.evaluate_background(
          three_step,invalid_position,invalid_background_3d);
  context.expect_true(direct_invalid_3d.code==adapter_invalid_3d.code &&
                          adapter_invalid_3d.code==
                              swcme::StatusCode::OutsideModelDomain,
                      "3-D failure status is unchanged by the adapter");

  // Direct shock records are compared with adapter source records before any
  // serialization.  The +X spherical fixture also makes the 1-D and 3-D
  // states physically equivalent without weakening their separate identity.
  swcme::acceleration::ShockAccelerationState direct_source_1d;
  context.expect_true(one.model().shock_acceleration_state_checked(
                          one_step,direct_source_1d).ok(),
                      "direct 1-D shock/source state is available");
  swcme::sep::SEPSourceState adapter_source_1d;
  context.expect_true(one.source_at_shock(one_step,adapter_source_1d).ok(),
                      "AMPS-facing 1-D source state is available");
  compare_source_mapping(
      context,"1-D",direct_source_1d,
      one_step.shock_jump.upstream.pressure_Pa,
      swcme::solarwind::parker_focusing_length_m(
          one_step.common.solar_wind,direct_source_1d.radius_m,
          one_step.common.solar_wind.reference_sin_theta),adapter_source_1d);
  compare_energy_grid(context,"1-D",direct_source_1d,adapter_source_1d);

  // The serialized AMPS handoff must name every newly exposed quantity with
  // its unit-bearing field name.  Checking the header separately from values
  // prevents a correct in-memory record from being shifted or omitted in a
  // campaign file consumed outside this process.
  const std::string source_header=swcme::sep::source_csv_header();
  context.expect_true(source_header.find("upstream_pressure_Pa")!=
                          std::string::npos,
                      "source CSV declares upstream pressure units");
  context.expect_true(source_header.find("focusing_length_m")!=
                          std::string::npos,
                      "source CSV declares focusing-length units");
  context.expect_true(source_header.find("field_line_path_length_m")!=
                          std::string::npos,
                      "source CSV declares observer path-length units");

  const double shock_direction[3]={1.0,0.0,0.0};
  swcme3d::LocalShockState direct_shock_3d;
  context.expect_true(three.model().shock_state_direction_checked(
                          three_step,shock_direction,direct_shock_3d).ok(),
                      "direct 3-D local shock state is available");
  swcme::acceleration::ShockAccelerationState direct_source_3d;
  context.expect_true(three.model().shock_acceleration_state_checked(
                          three_step,shock_direction,direct_source_3d).ok(),
                      "direct 3-D acceleration state is available");
  swcme::sep::SEPSourceState adapter_source_3d;
  context.expect_true(three.source_at_direction(
                          three_step,{{1.0,0.0,0.0}},adapter_source_3d).ok(),
                      "AMPS-facing 3-D directional source is available");
  compare_source_mapping(
      context,"3-D",direct_source_3d,direct_shock_3d.upstream.pressure_Pa,
      swcme::solarwind::parker_focusing_length_m(
          three_step.common.solar_wind,direct_source_3d.radius_m,1.0),
      adapter_source_3d);
  compare_energy_grid(context,"3-D",direct_source_3d,adapter_source_3d);

  // Finally compare the operational observer path.  The adapter returns a
  // copy of the complete connectivity state plus a source record; both the
  // selected path length and local focusing must match direct SWCME exactly.
  const auto connection_step=three.prepare(8.0*3600.0);
  const std::array<double,3> observer={{swcme::constants::AU_M,0.0,0.0}};
  const swcme3d::ConnectivityState direct_connection=
      three.model().observer_connectivity(connection_step,observer.data());
  swcme3d::ConnectivityState adapter_connection;
  swcme::sep::SEPSourceState connected_source;
  const swcme::ModelStatus connection_status=
      three.source_at_observer_cobpoint(
          connection_step,observer,connected_source,&adapter_connection);
  context.expect_true(connection_status.ok() && direct_connection.connected &&
                          adapter_connection.connected,
                      "direct and adapter observer paths are connected");
  context.expect_true(adapter_connection.status==direct_connection.status &&
                          adapter_connection.selected_root==
                              direct_connection.selected_root &&
                          adapter_connection.roots.size()==
                              direct_connection.roots.size(),
                      "adapter preserves connectivity status and root selection");
  if (direct_connection.connected && !direct_connection.roots.empty() &&
      direct_connection.selected_root<direct_connection.roots.size()) {
    const auto& direct_root=
        direct_connection.roots[direct_connection.selected_root];
    const auto& adapter_root=
        adapter_connection.roots[adapter_connection.selected_root];
    expect_exact(context,"connectivity path length",adapter_root.path_length_m,
                 direct_root.path_length_m);
    expect_exact(context,"source path length",
                 connected_source.field_line_path_length_m,
                 direct_root.path_length_m);

    // The 1-D shared-core helper is also the implementation now used by the
    // direct 3-D connectivity class.  Comparing it here closes the entire
    // common-core -> direct model -> adapter chain for the selected cobpoint.
    expect_exact(context,"common/direct Parker path length",
                 direct_root.path_length_m,
                 swcme::solarwind::parker_path_length_m(
                     connection_step.common.solar_wind,1.0,
                     direct_root.radius_m,swcme::constants::AU_M));

    const std::array<double,3> source_hat=normalized(connected_source.position_m);
    const std::array<double,3> axis={{connection_step.solar_axis_hat[0],
                                      connection_step.solar_axis_hat[1],
                                      connection_step.solar_axis_hat[2]}};
    const double source_radius=std::hypot(
        connected_source.position_m[0],
        std::hypot(connected_source.position_m[1],
                   connected_source.position_m[2]));
    expect_exact(context,"connected-source focusing",
                 connected_source.focusing_length_m,
                 swcme::solarwind::parker_focusing_length_m(
                     connection_step.common.solar_wind,source_radius,
                     swcme::solarwind::parker_sin_colatitude(axis,source_hat)));
  }
}
