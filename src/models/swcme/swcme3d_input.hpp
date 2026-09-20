#ifndef SWCME3D_INPUT_HPP
#define SWCME3D_INPUT_HPP

// Canonical, complete standalone-input resolver for the 3-D SWCME model.
//
// Common solar-wind, kinematic, region, and SEP-spectrum assignments are
// deliberately delegated to swcme1d_input.hpp.  The 1-D and 3-D applications
// therefore cannot acquire different unit conversions or validation rules for
// the same physics.  This file owns only the genuinely three-dimensional
// geometry, rotation-axis, and shock-surface discretization fields.

#include "swcme1d_input.hpp"
#include "swcme3d.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace swcme {
namespace input3d {

using Assignment = input1d::Assignment;
using Code = input1d::Code;
using Status = input1d::Status;

struct ResolvedConfiguration {
  swcme3d::Params model;
  sep::SpectrumConfig spectrum;
  double injection_efficiency = 0.0;
  double launch_epoch_s = 0.0;
  double valid_from_s = 0.0;
  double valid_until_s = 0.0;
  std::size_t surface_theta_intervals = 0;
  std::size_t surface_phi_points = 0;
  std::string normalized_manifest;
  std::string fingerprint;
};

struct ResolveResult {
  Status status;
  ResolvedConfiguration configuration;
  bool ok() const { return status.ok(); }
};

namespace detail {

inline Status error(Code code, const Assignment& assignment,
                    const std::string& message) {
  Status result;
  result.code = code;
  result.key = assignment.key;
  result.layer = "srcSEP3D [swcme]";
  result.origin = assignment.origin;
  result.line = assignment.line;
  result.message = message;
  return result;
}

inline bool parse_double(const std::string& text, double* value) {
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (text.empty() || end == text.c_str() || *end != '\0' ||
      errno == ERANGE || !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

inline bool parse_size(const std::string& text, std::size_t* value) {
  if (text.empty() || text[0] == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE ||
      parsed > std::numeric_limits<std::size_t>::max()) return false;
  *value = static_cast<std::size_t>(parsed);
  return true;
}

inline const std::set<std::string>& geometry_keys() {
  static const std::set<std::string> keys = {
      "geometry.shape", "geometry.axis_ratio_y", "geometry.axis_ratio_z",
      "geometry.half_width_rad", "geometry.cme_direction_x",
      "geometry.cme_direction_y", "geometry.cme_direction_z",
      "geometry.solar_rotation_axis_x", "geometry.solar_rotation_axis_y",
      "geometry.solar_rotation_axis_z",
      "parker.solar_rotation_rate_rad_per_s",
      "surface.theta_intervals", "surface.phi_points"};
  return keys;
}

inline const std::set<std::string>& common_required_keys() {
  // These names match the canonical SWCME1D resolver exactly.  All are
  // required even when a field is inactive in the selected mode so a later
  // reviewed mode change cannot revive a hidden preset value.
  static const std::set<std::string> keys = {
      "preset", "ambient.wind_speed", "ambient.density_1au",
      "ambient.magnetic_field_1au", "ambient.proton_temperature",
      "ambient.adiabatic_index", "ambient.alpha_to_proton_ratio",
      "ambient.electron_temperature", "ambient.alpha_temperature",
      "ambient.thermodynamic_closure", "parker.radial_polarity",
      "parker.sin_theta", "parker.source_radius", "cme.kinematics",
      "cme.launch_radius", "cme.launch_speed", "cme.drag_coefficient",
      "cme.extrapolation", "shock.region_mode",
      "shock.acceleration_mode", "shock.relative_source_weight_per_area",
      "geometry.sheath_thickness_1au", "geometry.ejecta_thickness_1au",
      "smoothing.shock_width_1au", "smoothing.leading_edge_width_1au",
      "smoothing.trailing_edge_width_1au", "sheath.ramp_power",
      "sheath.leading_edge_speed_factor", "ejecta.density_factor",
      "ejecta.speed_factor", "event.launch_epoch", "event.valid_from",
      "event.valid_until", "source.particle_mass", "source.charge_number",
      "source.energy_min", "source.energy_max", "source.reference_energy",
      "source.injection_efficiency", "source.normalization",
      "source.reference_intensity_si"};
  return keys;
}

inline void copy_common(const swcme1d::Params& from, swcme3d::Params* to) {
  // Assignment-by-field is intentional: it is an auditable statement of the
  // common physics shared by the dimensional interfaces and avoids any ABI or
  // layout assumption between the two public Params structures.
  to->sin_theta = from.sin_theta;
  to->parker_source_radius_Rs = from.parker_source_radius_Rs;
  to->kinematics_mode = from.kinematics_mode;
  to->r0_Rs = from.r0_Rs;
  to->V0_sh_kms = from.V0_sh_kms;
  to->V_sw_kms = from.V_sw_kms;
  to->Gamma_kmInv = from.Gamma_kmInv;
  to->data_time_s = from.data_time_s;
  to->data_radius_Rs = from.data_radius_Rs;
  to->data_extrapolation = from.data_extrapolation;
  to->n1AU_cm3 = from.n1AU_cm3;
  to->B1AU_nT = from.B1AU_nT;
  to->T_K = from.T_K;
  to->gamma_ad = from.gamma_ad;
  to->thermodynamic_closure = from.thermodynamic_closure;
  to->alpha_to_proton_ratio = from.alpha_to_proton_ratio;
  to->electron_T_K = from.electron_T_K;
  to->alpha_T_K = from.alpha_T_K;
  to->parker_radial_polarity = from.parker_radial_polarity;
  to->region_mode = from.region_mode;
  to->shock_acceleration_mode = from.shock_acceleration_mode;
  to->relative_source_weight_per_area = from.relative_source_weight_per_area;
  to->sheath_thick_AU_at1AU = from.sheath_thick_AU_at1AU;
  to->ejecta_thick_AU_at1AU = from.ejecta_thick_AU_at1AU;
  to->edge_smooth_shock_AU_at1AU = from.edge_smooth_shock_AU_at1AU;
  to->edge_smooth_le_AU_at1AU = from.edge_smooth_le_AU_at1AU;
  to->edge_smooth_te_AU_at1AU = from.edge_smooth_te_AU_at1AU;
  to->V_sheath_LE_factor = from.V_sheath_LE_factor;
  to->V_ME_factor = from.V_ME_factor;
  to->sheath_ramp_power = from.sheath_ramp_power;
  to->f_ME = from.f_ME;
}

inline std::string spectrum_manifest(const ResolvedConfiguration& c) {
  std::ostringstream out;
  out << std::scientific << std::setprecision(17)
      << "event.launch_epoch_s=" << c.launch_epoch_s << '\n'
      << "event.valid_from_s=" << c.valid_from_s << '\n'
      << "event.valid_until_s=" << c.valid_until_s << '\n'
      << "surface.theta_intervals=" << c.surface_theta_intervals << '\n'
      << "surface.phi_points=" << c.surface_phi_points << '\n'
      << "source.particle_mass_kg=" << c.spectrum.particle_mass_kg << '\n'
      << "source.charge_number=" << c.spectrum.charge_number << '\n'
      << "source.energy_min_MeV=" << c.spectrum.kinetic_energy_min_MeV << '\n'
      << "source.energy_max_MeV=" << c.spectrum.kinetic_energy_max_MeV << '\n'
      << "source.reference_energy_MeV=" << c.spectrum.reference_energy_MeV << '\n'
      << "source.normalization="
      << sep::normalization_mode_name(c.spectrum.normalization) << '\n'
      << "source.reference_intensity_si=";
  if (std::isfinite(c.spectrum.reference_differential_intensity_SI))
    out << c.spectrum.reference_differential_intensity_SI;
  else
    out << "NA";
  out << '\n' << "source.injection_efficiency="
      << c.injection_efficiency << '\n';
  return out.str();
}

}  // namespace detail

inline ResolveResult Resolve(const std::vector<Assignment>& assignments) {
  ResolveResult result;
  std::set<std::string> seen;
  std::vector<Assignment> common;
  std::vector<Assignment> geometry;
  std::string kinematics;
  for (const Assignment& raw : assignments) {
    Assignment assignment = raw;
    assignment.key = input1d::detail::canonical_key(assignment.key);
    if (!seen.insert(assignment.key).second) {
      result.status = detail::error(
          Code::DuplicateKey, assignment, "key occurs more than once");
      return result;
    }
    if (assignment.key == "cme.kinematics")
      kinematics = input1d::detail::canonical_key(assignment.value);
    if (detail::geometry_keys().count(assignment.key) != 0)
      geometry.push_back(assignment);
    else
      common.push_back(assignment);
  }
  for (const std::string& key : detail::common_required_keys()) {
    if (seen.count(key) == 0) {
      Assignment missing;
      missing.key = key;
      result.status = detail::error(
          Code::InvalidCombination, missing, "required key is missing");
      return result;
    }
  }
  for (const std::string& key : detail::geometry_keys()) {
    if (seen.count(key) == 0) {
      Assignment missing;
      missing.key = key;
      result.status = detail::error(
          Code::InvalidCombination, missing, "required 3-D key is missing");
      return result;
    }
  }
  const bool hasTimes = seen.count("cme.data_times") != 0;
  const bool hasRadii = seen.count("cme.data_radii") != 0;
  if (kinematics == "data_driven") {
    if (!hasTimes || !hasRadii) {
      Assignment missing;
      missing.key = "cme.data_times/cme.data_radii";
      result.status = detail::error(
          Code::InvalidCombination, missing,
          "data_driven kinematics requires both knot lists");
      return result;
    }
  } else if (hasTimes || hasRadii) {
    Assignment invalid;
    invalid.key = "cme.data_times/cme.data_radii";
    result.status = detail::error(
        Code::InvalidCombination, invalid,
        "data knots are legal only for data_driven kinematics");
    return result;
  }

  input1d::Layer commonLayer;
  commonLayer.name = "srcSEP3D [swcme]";
  commonLayer.assignments = common;
  const input1d::ResolveResult commonResolved =
      input1d::Resolve(input1d::Preset::Fast, {commonLayer});
  if (!commonResolved.ok()) {
    result.status = commonResolved.status;
    return result;
  }
  detail::copy_common(commonResolved.configuration.model,
                      &result.configuration.model);
  result.configuration.spectrum = commonResolved.configuration.source.spectrum;
  result.configuration.injection_efficiency =
      commonResolved.configuration.source.injection_efficiency;
  result.configuration.launch_epoch_s =
      commonResolved.configuration.launch_epoch_s;
  result.configuration.valid_from_s = commonResolved.configuration.valid_from_s;
  result.configuration.valid_until_s =
      commonResolved.configuration.valid_until_s;

  for (const Assignment& assignment : geometry) {
    const std::string& key = assignment.key;
    double value = 0.0;
    if (key == "geometry.shape") {
      const std::string name = input1d::detail::canonical_key(assignment.value);
      if (name == "sphere") result.configuration.model.shape = swcme3d::ShockShape::Sphere;
      else if (name == "ellipsoid") result.configuration.model.shape = swcme3d::ShockShape::Ellipsoid;
      else if (name == "sse") result.configuration.model.shape = swcme3d::ShockShape::SSE;
      else {
        result.status = detail::error(
            Code::InvalidValue, assignment, "expected sphere, ellipsoid, or sse");
        return result;
      }
      continue;
    }
    if (key == "surface.theta_intervals" || key == "surface.phi_points") {
      std::size_t parsed = 0;
      if (!detail::parse_size(assignment.value, &parsed) || parsed < 2) {
        result.status = detail::error(
            Code::InvalidValue, assignment, "surface resolution must be an integer >= 2");
        return result;
      }
      if (key == "surface.theta_intervals")
        result.configuration.surface_theta_intervals = parsed;
      else
        result.configuration.surface_phi_points = parsed;
      continue;
    }
    if (!detail::parse_double(assignment.value, &value)) {
      result.status = detail::error(
          Code::InvalidValue, assignment,
          "3-D geometry values are bare SI/dimensionless scalars named by the key");
      return result;
    }
    if (key == "geometry.axis_ratio_y") result.configuration.model.axis_ratio_y = value;
    else if (key == "geometry.axis_ratio_z") result.configuration.model.axis_ratio_z = value;
    else if (key == "geometry.half_width_rad") result.configuration.model.half_width_rad = value;
    else if (key == "geometry.cme_direction_x") result.configuration.model.cme_dir[0] = value;
    else if (key == "geometry.cme_direction_y") result.configuration.model.cme_dir[1] = value;
    else if (key == "geometry.cme_direction_z") result.configuration.model.cme_dir[2] = value;
    else if (key == "geometry.solar_rotation_axis_x") result.configuration.model.solar_rotation_axis[0] = value;
    else if (key == "geometry.solar_rotation_axis_y") result.configuration.model.solar_rotation_axis[1] = value;
    else if (key == "geometry.solar_rotation_axis_z") result.configuration.model.solar_rotation_axis[2] = value;
    else if (key == "parker.solar_rotation_rate_rad_per_s")
      result.configuration.model.solar_rotation_rate_rad_s = value;
  }

  const config::ValidationResult modelStatus =
      swcme3d::validate_params(result.configuration.model);
  if (!modelStatus.ok()) {
    Assignment invalid;
    invalid.key = "canonical_model";
    result.status = detail::error(
        Code::CanonicalValidationFailure, invalid, modelStatus.summary());
    return result;
  }
  const ModelStatus spectrumStatus =
      sep::validate_spectrum_config(result.configuration.spectrum);
  if (!spectrumStatus.ok()) {
    Assignment invalid;
    invalid.key = "canonical_source";
    result.status = detail::error(
        Code::CanonicalValidationFailure, invalid, spectrumStatus.summary());
    return result;
  }

  result.configuration.normalized_manifest =
      std::string("swcme3d_input_schema=1\n") +
      swcme3d::resolved_configuration_manifest(result.configuration.model) +
      detail::spectrum_manifest(result.configuration);
  result.configuration.fingerprint = input1d::detail::fingerprint(
      result.configuration.normalized_manifest);
  result.status = Status{};
  return result;
}

}  // namespace input3d
}  // namespace swcme

#endif  // SWCME3D_INPUT_HPP
