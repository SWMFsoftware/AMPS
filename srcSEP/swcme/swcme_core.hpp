#pragma once

// ============================================================================
// swcme_core.hpp
// ----------------------------------------------------------------------------
// Shared dimensionality-independent SWCME preparation layer.
//
// Both the 1-D and 3-D public models expose different geometry APIs, but their
// ambient solar-wind normalization and CME/shock apex kinematics are the same
// physics.  This component converts one common public-unit configuration to SI,
// prepares the shared solar-wind cache, and evaluates the common kinematics.
// Dimensional wrappers may mirror selected fields for source compatibility,
// but they no longer recompute these equations independently.
// ============================================================================

#include "swcme_kinematics.hpp"
#include "swcme_solarwind.hpp"
#include "swcme_units.hpp"

#include <vector>

namespace swcme {
namespace core {

struct CommonConfig {
  // Ambient public units.
  double V_sw_kms = 400.0;
  double n1AU_cm3 = 5.0;
  double B1AU_nT = 5.0;
  double T_K = 1.2e5;
  double gamma_ad = 5.0 / 3.0;
  double parker_reference_sin_theta = 1.0;
  double solar_rotation_rate_rad_s =
      swcme::constants::SOLAR_ROTATION_RAD_S;

  // Common apex kinematics public units.
  swcme::kinematics::Mode kinematics_mode = swcme::kinematics::Mode::DBM;
  double r0_Rs = 20.0;
  double V0_sh_kms = 1500.0;
  double Gamma_kmInv = 1.0e-7;
  std::vector<double> data_time_s;
  std::vector<double> data_radius_Rs;
  swcme::kinematics::ExtrapolationPolicy data_extrapolation =
      swcme::kinematics::ExtrapolationPolicy::OutsideTime;
};

struct PreparedState {
  swcme::solarwind::PreparedState solar_wind;
  swcme::kinematics::State apex;
  double r0_m = 0.0;
};

inline PreparedState prepare(const CommonConfig& cfg, double time_s) {
  PreparedState state;
  state.r0_m = swcme::units::solar_radii_to_m(cfg.r0_Rs);

  swcme::solarwind::ConfigSI ambient;
  ambient.V_sw_m_s = swcme::units::km_per_s_to_m_per_s(cfg.V_sw_kms);
  ambient.n1AU_m3 = swcme::units::cm3_to_m3(cfg.n1AU_cm3);
  ambient.B1AU_T = swcme::units::nT_to_T(cfg.B1AU_nT);
  ambient.T_K = cfg.T_K;
  ambient.gamma_ad = cfg.gamma_ad;
  ambient.reference_sin_theta = cfg.parker_reference_sin_theta;
  ambient.solar_rotation_rate_rad_s = cfg.solar_rotation_rate_rad_s;
  state.solar_wind = swcme::solarwind::prepare(ambient);

  swcme::kinematics::Config kin;
  kin.mode = cfg.kinematics_mode;
  kin.r0_m = state.r0_m;
  kin.V0_m_s = swcme::units::km_per_s_to_m_per_s(cfg.V0_sh_kms);
  kin.Vsw_m_s = state.solar_wind.V_sw_m_s;
  kin.Gamma_m_inv = swcme::units::km_inverse_to_m_inverse(cfg.Gamma_kmInv);
  kin.extrapolation = cfg.data_extrapolation;
  kin.data_time_s = cfg.data_time_s;
  kin.data_radius_m.reserve(cfg.data_radius_Rs.size());
  for (double radius_Rs : cfg.data_radius_Rs) {
    kin.data_radius_m.push_back(swcme::units::solar_radii_to_m(radius_Rs));
  }
  state.apex = swcme::kinematics::evaluate(kin, time_s);
  return state;
}

}  // namespace core
}  // namespace swcme
