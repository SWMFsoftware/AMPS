#ifndef SWCME_SEP_SOURCE_HPP
#define SWCME_SEP_SOURCE_HPP

// ============================================================================
// swcme_sep_source.hpp
// ----------------------------------------------------------------------------
// Common SWCME -> SEP source contract.
//
// This header is deliberately independent of the 1-D and 3-D model classes.
// It converts the already validated common ShockAccelerationState into a
// transport-facing SEPSourceState with explicit units and spectrum semantics.
// Both dimensional adapters therefore use exactly the same source algebra.
//
// Spectrum convention
// -------------------
// In SOURCE mode SWCME reports the test-particle DSA phase-space slope
//
//      f(p) proportional to p^{-q},       q = 3 r_c/(r_c-1).
//
// For an isotropic distribution the differential directional intensity obeys
// J(E)=p^2 f(p), because dE/dp=v.  The exact relativistic shape used here is
// therefore
//
//      J(E)/J(E_ref) = [p(E)/p(E_ref)]^{2-q},
//
// with
//
//      p c = sqrt(K (K + 2 m c^2)).
//
// No non-relativistic power-law approximation is required.  The helper still
// exposes the familiar non-relativistic energy index (q-2)/2 as a diagnostic.
//
// Normalization convention
// ------------------------
// RelativeOnly is the canonical controlled-transport baseline.  It carries a
// dimensionless area weight but deliberately does NOT invent an absolute
// injection rate.  ReferenceDifferentialIntensity instead accepts a physical
// J(E_ref) in SI [particles m^-2 s^-1 sr^-1 J^-1].  This explicit choice keeps
// source normalization auditable and prevents a unitless test weight from
// being mistaken for a physical particle flux.
// ============================================================================

#include "swcme_acceleration.hpp"
#include "swcme_constants.hpp"
#include "swcme_status.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace swcme {
namespace sep {

enum class NormalizationMode {
  RelativeOnly,
  ReferenceDifferentialIntensity
};

inline const char* normalization_mode_name(NormalizationMode mode) {
  return mode == NormalizationMode::RelativeOnly
             ? "RELATIVE_ONLY"
             : "REFERENCE_DIFFERENTIAL_INTENSITY";
}

// External source-spectrum configuration.  Kinetic energies are expressed in
// MeV at the public interface, while all momentum/intensity calculations are
// performed in SI.  charge_number is retained for rigidity diagnostics; the
// DSA momentum shape itself depends only on mass and kinetic energy.
struct SpectrumConfig {
  double particle_mass_kg = swcme::constants::PROTON_MASS_KG;
  int charge_number = 1;
  double kinetic_energy_min_MeV = 1.0;
  double kinetic_energy_max_MeV = 1000.0;
  double reference_energy_MeV = 10.0;
  NormalizationMode normalization = NormalizationMode::RelativeOnly;

  // Used only by ReferenceDifferentialIntensity.  SI convention:
  // particles m^-2 s^-1 sr^-1 J^-1.  A NaN is intentional in RelativeOnly
  // mode so downstream code cannot accidentally treat the relative baseline
  // as a calibrated physical flux.
  double reference_differential_intensity_SI =
      std::numeric_limits<double>::quiet_NaN();
};

inline ModelStatus validate_spectrum_config(const SpectrumConfig& c) {
  if (!std::isfinite(c.particle_mass_kg) || !(c.particle_mass_kg > 0.0)) {
    return ModelStatus::make_value(StatusCode::InvalidConfiguration,
                                   "sep spectrum particle_mass_kg",
                                   c.particle_mass_kg);
  }
  if (c.charge_number == 0) {
    return ModelStatus::make(StatusCode::InvalidConfiguration,
                             "sep spectrum charge_number");
  }
  if (!std::isfinite(c.kinetic_energy_min_MeV) ||
      !(c.kinetic_energy_min_MeV > 0.0)) {
    return ModelStatus::make_value(StatusCode::InvalidConfiguration,
                                   "sep spectrum kinetic_energy_min_MeV",
                                   c.kinetic_energy_min_MeV);
  }
  if (!std::isfinite(c.kinetic_energy_max_MeV) ||
      !(c.kinetic_energy_max_MeV > c.kinetic_energy_min_MeV)) {
    return ModelStatus::make_value(StatusCode::InvalidConfiguration,
                                   "sep spectrum kinetic_energy_max_MeV",
                                   c.kinetic_energy_max_MeV);
  }
  if (!std::isfinite(c.reference_energy_MeV) ||
      c.reference_energy_MeV < c.kinetic_energy_min_MeV ||
      c.reference_energy_MeV > c.kinetic_energy_max_MeV) {
    return ModelStatus::make_value(StatusCode::InvalidConfiguration,
                                   "sep spectrum reference_energy_MeV",
                                   c.reference_energy_MeV);
  }
  if (c.normalization == NormalizationMode::ReferenceDifferentialIntensity &&
      (!std::isfinite(c.reference_differential_intensity_SI) ||
       !(c.reference_differential_intensity_SI > 0.0))) {
    return ModelStatus::make_value(
        StatusCode::InvalidConfiguration,
        "sep spectrum reference_differential_intensity_SI",
        c.reference_differential_intensity_SI);
  }
  return ModelStatus::success();
}

// Relativistic momentum from kinetic energy.  The returned SI momentum is
// kg m/s.  No ultrarelativistic/nonrelativistic branch is needed.
inline double momentum_kg_m_s_from_kinetic_MeV(double kinetic_energy_MeV,
                                                double mass_kg) {
  if (!std::isfinite(kinetic_energy_MeV) || kinetic_energy_MeV < 0.0 ||
      !std::isfinite(mass_kg) || !(mass_kg > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double K_J = kinetic_energy_MeV * swcme::constants::MEV_TO_J;
  const double mc2 = mass_kg * swcme::constants::SPEED_OF_LIGHT_M_S *
                     swcme::constants::SPEED_OF_LIGHT_M_S;
  const double pc = std::sqrt(K_J * (K_J + 2.0 * mc2));
  return pc / swcme::constants::SPEED_OF_LIGHT_M_S;
}

// Rigidity in GV for an integer charge number Z.  This is a diagnostic helper
// for SEP interfaces that use rigidity bins; the source spectrum remains
// parameterized by kinetic energy and phase-space momentum slope.
inline double rigidity_GV_from_kinetic_MeV(double kinetic_energy_MeV,
                                           double mass_kg,
                                           int charge_number) {
  if (charge_number == 0) return std::numeric_limits<double>::quiet_NaN();
  const double p = momentum_kg_m_s_from_kinetic_MeV(kinetic_energy_MeV,mass_kg);
  if (!std::isfinite(p)) return std::numeric_limits<double>::quiet_NaN();
  const double pc_J = p * swcme::constants::SPEED_OF_LIGHT_M_S;
  const double volts = pc_J /
      (std::abs(charge_number) * swcme::constants::ELEMENTARY_CHARGE_C);
  return volts * 1.0e-9;
}

// Common heliophysics intensity units <-> SI.
//   common: particles cm^-2 s^-1 sr^-1 MeV^-1
//   SI:     particles m^-2  s^-1 sr^-1 J^-1
inline double differential_intensity_common_to_SI(double value) {
  return value * 1.0e4 / swcme::constants::MEV_TO_J;
}
inline double differential_intensity_SI_to_common(double value) {
  return value * swcme::constants::MEV_TO_J / 1.0e4;
}

// Stable AMPS-facing source record.  `active` means that SOURCE mode is active
// and a physical fast shock exists at this source location.  Connectivity is
// orthogonal: surface-source tables do not evaluate an observer connection,
// whereas a cobpoint record sets connection_evaluated=true.
struct SEPSourceState {
  ModelStatus status;
  bool active = false;
  bool connection_evaluated = false;
  bool connected = false;
  acceleration::Mode acceleration_mode = acceleration::Mode::Source;

  double time_s = 0.0;
  std::size_t source_id = 0;
  std::array<double,3> position_m{{0.0,0.0,0.0}};
  std::array<double,3> normal{{0.0,0.0,0.0}};

  // Optional source-patch geometry.  A point/cobpoint source may use zero area
  // and area_fraction=1.  Mesh-derived source surfaces use the physical cell
  // area and normalize area_fraction over ACTIVE physical-shock cells only.
  double patch_area_m2 = 0.0;
  double active_surface_area_m2 = 0.0;
  double area_fraction = 1.0;
  double relative_patch_weight = 0.0;

  double compression = 1.0;
  double theta_Bn_rad = 0.0;
  double fast_mach = 0.0;
  double normal_speed_m_s = 0.0;
  double upstream_density_m3 = 0.0;
  double upstream_B_T = 0.0;

  // Transport context attached by the dimensional adapter after the common
  // shock/source conversion.  Pressure uses the production proton-only
  // closure, and focusing uses the shared Parker derivative.  Path length is
  // defined only for an observer-connected cobpoint; directional and surface
  // sources retain NaN so zero cannot be mistaken for a colocated observer.
  double upstream_pressure_Pa = 0.0;
  double focusing_length_m = std::numeric_limits<double>::quiet_NaN();
  double field_line_path_length_m = std::numeric_limits<double>::quiet_NaN();

  double q_phase_space = std::numeric_limits<double>::quiet_NaN();
  double momentum_intensity_index = std::numeric_limits<double>::quiet_NaN();
  double nonrel_energy_intensity_index = std::numeric_limits<double>::quiet_NaN();

  SpectrumConfig spectrum;

  // Copy of the source-per-area control from ShockAccelerationState.  The mesh
  // adapter converts it into relative_patch_weight after active area is known.
  double relative_source_weight_per_area = 0.0;
};

inline ModelStatus make_source_state(
    const acceleration::ShockAccelerationState& acceleration_state,
    const SpectrumConfig& spectrum,
    bool connection_evaluated,
    bool connected,
    std::size_t source_id,
    double patch_area_m2,
    SEPSourceState& out) {
  out = SEPSourceState{};
  out.status = validate_spectrum_config(spectrum);
  if (out.status.failure()) return out.status;
  if (!std::isfinite(patch_area_m2) || patch_area_m2 < 0.0) {
    out.status = ModelStatus::make_value(StatusCode::InvalidConfiguration,
                                         "sep source patch_area_m2",
                                         patch_area_m2);
    return out.status;
  }

  out.connection_evaluated = connection_evaluated;
  out.connected = connection_evaluated ? connected : false;
  out.acceleration_mode = acceleration_state.mode;
  out.time_s = acceleration_state.time_s;
  out.source_id = source_id;
  out.position_m = acceleration_state.position_m;
  out.normal = acceleration_state.normal;
  out.patch_area_m2 = patch_area_m2;
  out.compression = acceleration_state.compression;
  out.theta_Bn_rad = acceleration_state.theta_Bn_rad;
  out.fast_mach = acceleration_state.fast_mach;
  out.normal_speed_m_s = acceleration_state.normal_speed_m_s;
  out.upstream_density_m3 = acceleration_state.upstream_density_m3;
  out.upstream_B_T = acceleration_state.upstream_B_T;
  out.spectrum = spectrum;
  out.relative_source_weight_per_area =
      acceleration_state.relative_source_weight_per_area;

  out.active = acceleration_state.source_enabled &&
               acceleration_state.physical_shock &&
               (!connection_evaluated || connected);
  if (out.active) {
    out.q_phase_space = acceleration_state.dsa_q_phase_space;
    if (!std::isfinite(out.q_phase_space) || !(out.q_phase_space > 2.0)) {
      out.status = ModelStatus::make_value(StatusCode::NonFiniteResult,
                                           "sep source q_phase_space",
                                           out.q_phase_space);
      out.active = false;
      return out.status;
    }
    out.momentum_intensity_index = out.q_phase_space - 2.0;
    out.nonrel_energy_intensity_index = 0.5*(out.q_phase_space - 2.0);
    out.relative_patch_weight = acceleration_state.relative_source_weight_per_area;
  } else {
    // Deliberately unavailable when SOURCE is not active.  NaN makes misuse
    // visible in debug/output rather than silently carrying a stale slope.
    out.q_phase_space = std::numeric_limits<double>::quiet_NaN();
    out.momentum_intensity_index = std::numeric_limits<double>::quiet_NaN();
    out.nonrel_energy_intensity_index = std::numeric_limits<double>::quiet_NaN();
    out.relative_patch_weight = 0.0;
  }
  out.status = ModelStatus::success();
  return out.status;
}

// Exact DSA intensity shape normalized to unity at E_ref.  It is valid only
// for an active SOURCE record and within the configured source energy range.
inline ModelStatus relative_intensity_shape(const SEPSourceState& source,
                                            double kinetic_energy_MeV,
                                            double& shape) {
  shape = std::numeric_limits<double>::quiet_NaN();
  if (!source.active) {
    return ModelStatus::make(StatusCode::SourceInactive,
                             "sep relative_intensity_shape");
  }
  if (!std::isfinite(kinetic_energy_MeV) ||
      kinetic_energy_MeV < source.spectrum.kinetic_energy_min_MeV ||
      kinetic_energy_MeV > source.spectrum.kinetic_energy_max_MeV) {
    return ModelStatus::make_value(StatusCode::OutsideModelDomain,
                                   "sep kinetic energy MeV",
                                   kinetic_energy_MeV);
  }
  const double p = momentum_kg_m_s_from_kinetic_MeV(
      kinetic_energy_MeV,source.spectrum.particle_mass_kg);
  const double p_ref = momentum_kg_m_s_from_kinetic_MeV(
      source.spectrum.reference_energy_MeV,source.spectrum.particle_mass_kg);
  if (!std::isfinite(p) || !std::isfinite(p_ref) || !(p > 0.0) || !(p_ref > 0.0)) {
    return ModelStatus::make(StatusCode::NonFiniteResult,
                             "sep momentum conversion");
  }
  shape = std::pow(p/p_ref,2.0-source.q_phase_space);
  if (!std::isfinite(shape) || !(shape > 0.0)) {
    return ModelStatus::make_value(StatusCode::NonFiniteResult,
                                   "sep relative intensity shape",shape);
  }
  return ModelStatus::success();
}

// Physical differential intensity, available only when the source config uses
// ReferenceDifferentialIntensity.  This prevents an uncalibrated relative
// source from being silently exported with physical units.
inline ModelStatus differential_intensity_SI(const SEPSourceState& source,
                                             double kinetic_energy_MeV,
                                             double& intensity_SI) {
  intensity_SI = std::numeric_limits<double>::quiet_NaN();
  if (source.spectrum.normalization !=
      NormalizationMode::ReferenceDifferentialIntensity) {
    return ModelStatus::make(StatusCode::InvalidConfiguration,
                             "sep differential intensity normalization");
  }
  double shape=0.0;
  const ModelStatus status=relative_intensity_shape(source,kinetic_energy_MeV,shape);
  if (!status.ok()) return status;
  intensity_SI=source.spectrum.reference_differential_intensity_SI*shape;
  if (!std::isfinite(intensity_SI) || !(intensity_SI>0.0)) {
    return ModelStatus::make_value(StatusCode::NonFiniteResult,
                                   "sep differential intensity SI",
                                   intensity_SI);
  }
  return ModelStatus::success();
}

inline const char* source_csv_header() {
  return "status,active,connection_evaluated,connected,acceleration_mode,time_s,source_id,"
         "x_m,y_m,z_m,nx,ny,nz,patch_area_m2,active_surface_area_m2,area_fraction,"
         "relative_source_weight_per_area,relative_patch_weight,compression,theta_Bn_rad,"
         "fast_mach,Vsh_n_m_s,upstream_density_m3,upstream_B_T,upstream_pressure_Pa,"
         "focusing_length_m,field_line_path_length_m,q_phase_space,"
         "momentum_intensity_index,nonrel_energy_intensity_index,normalization,"
         "particle_mass_kg,charge_number,Emin_MeV,Emax_MeV,Eref_MeV,Jref_SI";
}

// Byte-stable scientific-notation record used by standalone diagnostics, the
// campaign runner and 1-D/3-D regression checks before AMPS propagates any
// particles.  Unavailable floating-point fields serialize as NA.
inline std::string serialize_source_csv(const SEPSourceState& s) {
  std::ostringstream out;
  out.setf(std::ios::scientific);
  out << std::setprecision(17)
      << status_code_name(s.status.code) << ','
      << (s.active?1:0) << ','
      << (s.connection_evaluated?1:0) << ','
      << (s.connected?1:0) << ','
      << acceleration::mode_name(s.acceleration_mode) << ','
      << s.time_s << ',' << s.source_id << ','
      << s.position_m[0] << ',' << s.position_m[1] << ',' << s.position_m[2] << ','
      << s.normal[0] << ',' << s.normal[1] << ',' << s.normal[2] << ','
      << s.patch_area_m2 << ',' << s.active_surface_area_m2 << ','
      << s.area_fraction << ',' << s.relative_source_weight_per_area << ','
      << s.relative_patch_weight << ','
      << s.compression << ',' << s.theta_Bn_rad << ',' << s.fast_mach << ','
      << s.normal_speed_m_s << ',' << s.upstream_density_m3 << ',' << s.upstream_B_T
      << ',' << s.upstream_pressure_Pa << ',';
  auto emit = [&](double value) {
    if (std::isfinite(value)) out << value;
    else out << "NA";
  };
  emit(s.focusing_length_m); out << ',';
  emit(s.field_line_path_length_m); out << ',';
  emit(s.q_phase_space); out << ',';
  emit(s.momentum_intensity_index); out << ',';
  emit(s.nonrel_energy_intensity_index); out << ',';
  out << normalization_mode_name(s.spectrum.normalization) << ','
      << s.spectrum.particle_mass_kg << ','
      << s.spectrum.charge_number << ','
      << s.spectrum.kinetic_energy_min_MeV << ','
      << s.spectrum.kinetic_energy_max_MeV << ','
      << s.spectrum.reference_energy_MeV << ',';
  emit(s.spectrum.reference_differential_intensity_SI);
  return out.str();
}

}  // namespace sep
}  // namespace swcme

#endif
