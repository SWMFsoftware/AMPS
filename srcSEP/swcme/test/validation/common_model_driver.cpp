// Shared production-API probe for observational validation cases VP06-VP15.
//
// The case-owned Python runners never reimplement SWCME physics.  They invoke
// this small warning-clean executable to obtain shock, Parker, and SEP records
// from the public production headers, then compare those records with an
// independently written reference_solution.py in the relevant case directory.
// Keeping the probe shared avoids eleven subtly different unit conversions;
// the selected mode and every fixture parameter are written into the CSV so a
// campaign artifact remains auditable without reading this source file.

#include <swcme_constants.hpp>
#include <swcme_sep_interface.hpp>
#include <swcme_shock.hpp>
#include <swcme_solarwind.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void shock_probe() {
  // This finite oblique fixture is intentionally the same physical state at
  // several normal angles.  VP06 consumes the conservation diagnostics while
  // VP07 consumes the independently reconstructible obliquity and downstream
  // geometry, so neither case silently substitutes a gas-dynamic jump.
  swcme::shock::PrimitiveState upstream;
  upstream.rho_kg_m3 = 5.0e6 * swcme::constants::PROTON_MASS_KG;
  upstream.pressure_Pa = 5.0e6 * swcme::constants::BOLTZMANN_J_K * 1.2e5;
  upstream.velocity_m_s = {{4.0e5, 0.0, 0.0}};
  upstream.magnetic_T = {{4.0e-9, 3.0e-9, 1.0e-9}};
  std::cout << "angle_deg,status,compression,theta_bn_deg,fast_mach,"
               "rho2_kg_m3,p2_pa,v2x_m_s,v2y_m_s,v2z_m_s,"
               "b2x_t,b2y_t,b2z_t,mass_residual,normal_b_residual,"
               "electric_residual,momentum_residual,energy_residual,entropy_ratio\n";
  for (double angle_deg : {0.0, 20.0, 40.0}) {
    const double angle = angle_deg * swcme::constants::PI / 180.0;
    const swcme::shock::Vec3 normal{{std::cos(angle), std::sin(angle), 0.0}};
    const auto result = swcme::shock::solve_ideal_mhd_fast_shock(
        upstream, normal, 1.2e6, 5.0 / 3.0);
    std::cout << angle_deg << ',' << swcme::shock::solve_status_name(result.status)
              << ',' << result.compression << ','
              << result.theta_Bn_rad * 180.0 / swcme::constants::PI << ','
              << result.fast_mach << ',' << result.downstream.rho_kg_m3 << ','
              << result.downstream.pressure_Pa << ','
              << result.downstream.velocity_m_s[0] << ','
              << result.downstream.velocity_m_s[1] << ','
              << result.downstream.velocity_m_s[2] << ','
              << result.downstream.magnetic_T[0] << ','
              << result.downstream.magnetic_T[1] << ','
              << result.downstream.magnetic_T[2] << ','
              << result.mass_residual << ',' << result.normal_B_residual << ','
              << result.electric_residual << ',' << result.momentum_residual
              << ',' << result.energy_residual << ',' << result.entropy_ratio
              << '\n';
  }
}

swcme::solarwind::PreparedState parker_state() {
  swcme::solarwind::ConfigSI config;
  config.V_sw_m_s = 4.0e5;
  config.n1AU_m3 = 6.0e6;
  config.B1AU_T = 5.0e-9;
  config.T_K = 1.2e5;
  config.parker_source_radius_m = 2.5 * swcme::constants::SOLAR_RADIUS_M;
  return swcme::solarwind::prepare(config);
}

void parker_probe() {
  const auto state = parker_state();
  const double source = 2.5 * swcme::constants::SOLAR_RADIUS_M;
  std::cout << "sin_theta,radius_au,path_length_m,focusing_length_m,bmag_t\n";
  for (double sin_theta : {0.0, 0.5, 1.0}) {
    for (double radius_au : {0.1, 0.3, 0.5, 1.0}) {
      const double radius = radius_au * swcme::constants::AU_M;
      const auto field = swcme::solarwind::parker_components(
          state, radius, sin_theta);
      std::cout << sin_theta << ',' << radius_au << ','
                << swcme::solarwind::parker_path_length_m(
                       state, sin_theta, source, radius)
                << ',' << swcme::solarwind::parker_focusing_length_m(
                              state, radius, sin_theta)
                << ',' << field.Bmag_T << '\n';
    }
  }
}

void equivalent_parameters(swcme1d::Params& one, swcme3d::Params& three) {
  // Every source-affecting control is copied explicitly.  Defaults in the two
  // dimensional parameter types are deliberately independent and therefore
  // cannot be assumed equal by an AMPS equivalence validation.
  one.V_sw_kms = 410.0;
  one.n1AU_cm3 = 5.5;
  one.B1AU_nT = 5.2;
  one.T_K = 1.25e5;
  one.gamma_ad = 5.0 / 3.0;
  one.sin_theta = 1.0;
  one.kinematics_mode = swcme::kinematics::Mode::DBM;
  one.r0_Rs = 20.0;
  one.V0_sh_kms = 1450.0;
  one.Gamma_kmInv = 5.0e-8;
  one.region_mode = swcme::regions::Mode::ShockOnly;
  one.shock_acceleration_mode = swcme::acceleration::Mode::Source;
  one.relative_source_weight_per_area = 1.0;
  three.V_sw_kms = one.V_sw_kms;
  three.n1AU_cm3 = one.n1AU_cm3;
  three.B1AU_nT = one.B1AU_nT;
  three.T_K = one.T_K;
  three.gamma_ad = one.gamma_ad;
  three.sin_theta = one.sin_theta;
  three.kinematics_mode = one.kinematics_mode;
  three.r0_Rs = one.r0_Rs;
  three.V0_sh_kms = one.V0_sh_kms;
  three.Gamma_kmInv = one.Gamma_kmInv;
  three.region_mode = one.region_mode;
  three.shock_acceleration_mode = one.shock_acceleration_mode;
  three.relative_source_weight_per_area = one.relative_source_weight_per_area;
  three.shape = swcme3d::ShockShape::Sphere;
  three.cme_dir[0] = 1.0;
  three.cme_dir[1] = 0.0;
  three.cme_dir[2] = 0.0;
  three.solar_rotation_axis[0] = 0.0;
  three.solar_rotation_axis[1] = 0.0;
  three.solar_rotation_axis[2] = 1.0;
}

void sep_probe() {
  swcme1d::Params one_parameters;
  swcme3d::Params three_parameters;
  equivalent_parameters(one_parameters, three_parameters);
  swcme::sep::SpectrumConfig spectrum;
  spectrum.kinetic_energy_min_MeV = 1.0;
  spectrum.kinetic_energy_max_MeV = 500.0;
  spectrum.reference_energy_MeV = 20.0;
  const swcme::sep::Interface1D one(one_parameters, spectrum);
  const swcme::sep::Interface3D three(three_parameters, spectrum);
  std::cout << "time_h,dimension,active,radius_m,compression,fast_mach,"
               "q_phase_space,source_weight,focusing_length_m,pressure_pa\n";
  for (double time_h : {2.0, 6.0, 12.0, 24.0, 36.0}) {
    const auto one_step = one.prepare(time_h * 3600.0);
    const auto three_step = three.prepare(time_h * 3600.0);
    swcme::sep::SEPSourceState one_source;
    swcme::sep::SEPSourceState three_source;
    const auto one_status = one.source_at_shock(one_step, one_source);
    const auto three_status = three.source_at_direction(
        three_step, {{1.0, 0.0, 0.0}}, three_source);
    if (!one_status.ok() || !three_status.ok())
      throw std::runtime_error("production SEP source query failed");
    for (const auto& item : {std::pair<const char*, const swcme::sep::SEPSourceState&>(
                                  "1D", one_source),
                              std::pair<const char*, const swcme::sep::SEPSourceState&>(
                                  "3D", three_source)}) {
      const auto& source = item.second;
      const double radius = std::hypot(source.position_m[0],
                                       std::hypot(source.position_m[1], source.position_m[2]));
      std::cout << time_h << ',' << item.first << ',' << (source.active ? 1 : 0)
                << ',' << radius << ',' << source.compression << ','
                << source.fast_mach << ',' << source.q_phase_space << ','
                << source.relative_patch_weight << ',' << source.focusing_length_m
                << ',' << source.upstream_pressure_Pa << '\n';
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("usage: common_model_driver shock|parker|sep");
    std::cout << std::setprecision(17);
    const std::string mode(argv[1]);
    if (mode == "shock") shock_probe();
    else if (mode == "parker") parker_probe();
    else if (mode == "sep") sep_probe();
    else throw std::runtime_error("unknown probe mode: " + mode);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "validation model probe: " << error.what() << '\n';
    return 2;
  }
}
