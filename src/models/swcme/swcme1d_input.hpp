#ifndef SWCME1D_INPUT_HPP
#define SWCME1D_INPUT_HPP

// ============================================================================
// Canonical typed input configuration for the 1-D SWCME/SEP interface.
//
// This header is the single textual schema for standalone applications. It is
// intentionally owned by src/models/swcme rather than srcSEP: both 1-D and 3-D
// applications can transport key/value assignments to this resolver without
// copying model fields, unit conversions, presets, or validation rules.
//
// Resolution order is deterministic:
//
//   named preset -> input-file layer -> command-line/programmatic layer
//
// A key may appear once per layer; a more authoritative layer may override it.
// Every quantity is converted at this boundary into the units of
// swcme1d::Params (or the explicitly documented SEP source units), then the
// canonical model/source validators run before a Model is constructed.
// ============================================================================

#include "swcme1d.hpp"
#include "swcme_sep_source.hpp"
#include "swcme_units.hpp"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace swcme {
namespace input1d {

enum class Preset { Fast, Slow };

inline const char* preset_name(Preset preset) {
  return preset == Preset::Fast ? "FAST" : "SLOW";
}

enum class Code {
  Ok,
  UnknownKey,
  DuplicateKey,
  InvalidValue,
  InvalidUnit,
  UnsupportedField,
  InvalidCombination,
  CanonicalValidationFailure
};

inline const char* code_name(Code code) {
  switch (code) {
    case Code::Ok: return "OK";
    case Code::UnknownKey: return "UNKNOWN_KEY";
    case Code::DuplicateKey: return "DUPLICATE_KEY";
    case Code::InvalidValue: return "INVALID_VALUE";
    case Code::InvalidUnit: return "INVALID_UNIT";
    case Code::UnsupportedField: return "UNSUPPORTED_FIELD";
    case Code::InvalidCombination: return "INVALID_COMBINATION";
    case Code::CanonicalValidationFailure:
      return "CANONICAL_VALIDATION_FAILURE";
  }
  return "UNKNOWN";
}

struct Status {
  Code code = Code::Ok;
  std::string key;
  std::string layer;
  // Human-readable assignment source (file, API, or CLI spelling). It is
  // diagnostic provenance only and is intentionally excluded from the
  // effective-physics fingerprint.
  std::string origin;
  std::size_t line = 0;
  std::string message;
  bool ok() const { return code == Code::Ok; }
};

// Assignment is deliberately a transport record, not a second parameter
// schema. The key is interpreted only by ApplyAssignment below.
struct Assignment {
  std::string key;
  std::string value;
  std::string origin;
  std::size_t line = 0;
};

struct Layer {
  std::string name;
  std::vector<Assignment> assignments;
};

struct SourceConfiguration {
  sep::SpectrumConfig spectrum;
  // srcSEP multiplies its swept-volume source by this dimensionless fraction.
  // Absolute differential-intensity normalization remains owned by spectrum.
  double injection_efficiency = 3.4e-4;
};

struct ResolvedConfiguration {
  swcme1d::Params model;
  SourceConfiguration source;
  Preset preset = Preset::Fast;

  // External simulation-clock contract. The canonical model itself consumes
  // launch-relative seconds: t_model = simulation_epoch - launch_epoch_s.
  double launch_epoch_s = 0.0;
  double valid_from_s = 0.0;
  double valid_until_s = std::numeric_limits<double>::infinity();

  std::string normalized_manifest;
  std::string fingerprint;
};

struct ResolveResult {
  Status status;
  ResolvedConfiguration configuration;
  bool ok() const { return status.ok(); }
};

namespace detail {

inline std::string trim(const std::string& input) {
  const std::string whitespace = " \t\r\n";
  const std::size_t begin = input.find_first_not_of(whitespace);
  if (begin == std::string::npos) return std::string();
  const std::size_t end = input.find_last_not_of(whitespace);
  return input.substr(begin,end-begin+1);
}

inline std::string lower(std::string value) {
  std::transform(value.begin(),value.end(),value.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

inline std::string canonical_key(std::string key) {
  key=lower(trim(key));
  std::replace(key.begin(),key.end(),'-','_');
  return key;
}

inline Status error(Code code,const Assignment& assignment,
                    const std::string& layer,const std::string& message) {
  Status status;
  status.code=code;
  status.key=assignment.key;
  status.layer=layer;
  status.origin=assignment.origin;
  status.line=assignment.line;
  status.message=message;
  return status;
}

struct Unit {
  const char* name;
  double factor;
};

// Parse one scalar and convert it directly into the destination public unit.
// `units` contains accepted spellings and their destination-unit factors.
inline Status scalar(const Assignment& assignment,const std::string& layer,
                     const Unit* units,std::size_t unit_count,
                     bool unit_required,double* destination) {
  const std::string text=trim(assignment.value);
  errno=0;
  char* end=nullptr;
  const double value=std::strtod(text.c_str(),&end);
  if (errno!=0 || end==text.c_str() || !std::isfinite(value))
    return error(Code::InvalidValue,assignment,layer,
                 "expected one finite numeric value");
  const std::string unit=lower(trim(end ? std::string(end) : std::string()));
  if (unit.empty() && !unit_required) {
    *destination=value;
    return Status{};
  }
  for (std::size_t i=0;i<unit_count;++i) {
    if (unit==units[i].name) {
      const double converted=value*units[i].factor;
      if (!std::isfinite(converted))
        return error(Code::InvalidValue,assignment,layer,
                     "unit conversion overflowed the destination quantity");
      *destination=converted;
      return Status{};
    }
  }
  std::ostringstream message;
  message << "invalid or missing unit; accepted:";
  for (std::size_t i=0;i<unit_count;++i) message << ' ' << units[i].name;
  return error(Code::InvalidUnit,assignment,layer,message.str());
}

inline Status integer(const Assignment& assignment,const std::string& layer,
                      int* destination) {
  const std::string text=trim(assignment.value);
  errno=0;
  char* end=nullptr;
  const long value=std::strtol(text.c_str(),&end,10);
  if (errno!=0 || end==text.c_str() ||
      !trim(end ? std::string(end) : "").empty() ||
      value<std::numeric_limits<int>::min() ||
      value>std::numeric_limits<int>::max())
    return error(Code::InvalidValue,assignment,layer,"expected one integer");
  *destination=static_cast<int>(value);
  return Status{};
}

inline Status scalar_list(const Assignment& assignment,const std::string& layer,
                          const Unit* units,std::size_t unit_count,
                          std::vector<double>* destination) {
  std::vector<double> values;
  std::size_t begin=0;
  while (begin<=assignment.value.size()) {
    const std::size_t comma=assignment.value.find(',',begin);
    Assignment item=assignment;
    item.value=assignment.value.substr(
        begin,comma==std::string::npos ? std::string::npos : comma-begin);
    double parsed=0.0;
    Status status=scalar(item,layer,units,unit_count,true,&parsed);
    if (!status.ok()) return status;
    values.push_back(parsed);
    if (comma==std::string::npos) break;
    begin=comma+1;
  }
  if (values.empty())
    return error(Code::InvalidValue,assignment,layer,"list must not be empty");
  *destination=values;
  return Status{};
}

inline std::uint64_t fnv1a(const std::string& text) {
  std::uint64_t value=UINT64_C(14695981039346656037);
  for (std::size_t i=0;i<text.size();++i) {
    value^=static_cast<unsigned char>(text[i]);
    value*=UINT64_C(1099511628211);
  }
  return value;
}

inline std::string fingerprint(const std::string& manifest) {
  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << fnv1a(manifest);
  return out.str();
}

inline swcme1d::Params preset_params(Preset preset) {
  swcme1d::Params p;
  if (preset==Preset::Fast) {
    p.V_sw_kms=400.0; p.n1AU_cm3=6.0; p.B1AU_nT=5.0; p.T_K=1.2e5;
    p.r0_Rs=1.05; p.V0_sh_kms=1900.0; p.Gamma_kmInv=8.0e-8;
    p.sheath_thick_AU_at1AU=0.12; p.ejecta_thick_AU_at1AU=0.22;
    p.edge_smooth_shock_AU_at1AU=0.010;
    p.edge_smooth_le_AU_at1AU=0.020;
    p.edge_smooth_te_AU_at1AU=0.030;
    p.sheath_comp_floor=1.25; p.sheath_ramp_power=2.0;
    p.V_sheath_LE_factor=1.12; p.f_ME=0.50; p.V_ME_factor=0.80;
  } else {
    p.V_sw_kms=380.0; p.n1AU_cm3=5.0; p.B1AU_nT=4.5; p.T_K=1.0e5;
    p.r0_Rs=1.05; p.V0_sh_kms=950.0; p.Gamma_kmInv=3.0e-8;
    p.sheath_thick_AU_at1AU=0.08; p.ejecta_thick_AU_at1AU=0.18;
    p.edge_smooth_shock_AU_at1AU=0.015;
    p.edge_smooth_le_AU_at1AU=0.030;
    p.edge_smooth_te_AU_at1AU=0.050;
    p.sheath_comp_floor=1.15; p.sheath_ramp_power=1.5;
    p.V_sheath_LE_factor=1.08; p.f_ME=0.60; p.V_ME_factor=0.90;
  }
  return p;
}

inline Status ApplyAssignment(const Assignment& assignment,
                              const std::string& layer,
                              ResolvedConfiguration* c,
                              bool* explicit_valid_from) {
  const std::string key=canonical_key(assignment.key);
  const Unit dimensionless[]={{"",1.0},{"1",1.0}};
  const Unit speed[]={{"km/s",1.0},{"m/s",1.0e-3}};
  const Unit density[]={{"cm^-3",1.0},{"cm-3",1.0},
                        {"m^-3",1.0e-6},{"m-3",1.0e-6}};
  const Unit magnetic[]={{"nt",1.0},{"t",1.0e9}};
  const Unit temperature[]={{"k",1.0}};
  const Unit radius_rs[]={{"rs",1.0},{"r_sun",1.0},
      {"m",1.0/swcme::constants::SOLAR_RADIUS_M},
      {"au",swcme::constants::AU_M/swcme::constants::SOLAR_RADIUS_M}};
  const Unit distance_au[]={{"au",1.0},{"m",1.0/swcme::constants::AU_M}};
  const Unit inverse_km[]={{"1/km",1.0},{"km^-1",1.0},
                           {"1/m",1.0e3},{"m^-1",1.0e3}};
  const Unit time_s[]={{"s",1.0},{"min",60.0},{"h",3600.0},{"hr",3600.0}};
  const Unit mass_kg[]={{"kg",1.0},{"mp",swcme::constants::PROTON_MASS_KG}};
  const Unit energy_mev[]={{"mev",1.0},{"kev",1.0e-3},
                           {"j",1.0/swcme::constants::MEV_TO_J}};

  swcme1d::Params& p=c->model;
  SourceConfiguration& source=c->source;
  Status status;
  // Preset selection is consumed in Resolve's pre-pass because it must happen
  // before any field override. It remains a recognized key here so the normal
  // duplicate/unknown-key machinery still owns its diagnostics.
  if (key=="preset") return Status{};
  if (key=="ambient.wind_speed")
    return scalar(assignment,layer,speed,2,true,&p.V_sw_kms);
  if (key=="ambient.density_1au")
    return scalar(assignment,layer,density,4,true,&p.n1AU_cm3);
  if (key=="ambient.magnetic_field_1au")
    return scalar(assignment,layer,magnetic,2,true,&p.B1AU_nT);
  if (key=="ambient.proton_temperature")
    return scalar(assignment,layer,temperature,1,true,&p.T_K);
  if (key=="ambient.adiabatic_index")
    return scalar(assignment,layer,dimensionless,2,false,&p.gamma_ad);
  if (key=="ambient.alpha_to_proton_ratio")
    return scalar(assignment,layer,dimensionless,2,false,
                  &p.alpha_to_proton_ratio);
  if (key=="ambient.electron_temperature")
    return scalar(assignment,layer,temperature,1,true,&p.electron_T_K);
  if (key=="ambient.alpha_temperature")
    return scalar(assignment,layer,temperature,1,true,&p.alpha_T_K);
  if (key=="ambient.thermodynamic_closure") {
    const std::string value=canonical_key(assignment.value);
    if (value=="proton_only")
      p.thermodynamic_closure=swcme::solarwind::ThermodynamicClosure::ProtonOnly;
    else if (value=="multi_species")
      p.thermodynamic_closure=swcme::solarwind::ThermodynamicClosure::MultiSpecies;
    else return error(Code::InvalidValue,assignment,layer,
                      "expected proton_only or multi_species");
    return Status{};
  }
  if (key=="parker.radial_polarity")
    return integer(assignment,layer,&p.parker_radial_polarity);
  if (key=="parker.sin_theta")
    return scalar(assignment,layer,dimensionless,2,false,&p.sin_theta);
  if (key=="parker.source_radius")
    return scalar(assignment,layer,radius_rs,4,true,&p.parker_source_radius_Rs);
  if (key=="cme.kinematics") {
    const std::string value=canonical_key(assignment.value);
    if (value=="ballistic") p.kinematics_mode=swcme::kinematics::Mode::Ballistic;
    else if (value=="dbm") p.kinematics_mode=swcme::kinematics::Mode::DBM;
    else if (value=="data_driven")
      p.kinematics_mode=swcme::kinematics::Mode::DataDriven;
    else return error(Code::InvalidValue,assignment,layer,
                      "expected ballistic, dbm, or data_driven");
    return Status{};
  }
  if (key=="cme.launch_radius")
    return scalar(assignment,layer,radius_rs,4,true,&p.r0_Rs);
  if (key=="cme.launch_speed")
    return scalar(assignment,layer,speed,2,true,&p.V0_sh_kms);
  if (key=="cme.drag_coefficient")
    return scalar(assignment,layer,inverse_km,4,true,&p.Gamma_kmInv);
  if (key=="cme.data_times")
    return scalar_list(assignment,layer,time_s,4,&p.data_time_s);
  if (key=="cme.data_radii")
    return scalar_list(assignment,layer,radius_rs,4,&p.data_radius_Rs);
  if (key=="cme.extrapolation") {
    const std::string value=canonical_key(assignment.value);
    if (value=="outside_time" || value=="reject")
      p.data_extrapolation=swcme::kinematics::ExtrapolationPolicy::OutsideTime;
    else if (value=="ballistic")
      p.data_extrapolation=swcme::kinematics::ExtrapolationPolicy::Ballistic;
    else return error(Code::InvalidValue,assignment,layer,
                      "expected outside_time or ballistic");
    return Status{};
  }
  if (key=="shock.region_mode") {
    const std::string value=canonical_key(assignment.value);
    if (value=="shock_only") p.region_mode=swcme::regions::Mode::ShockOnly;
    else if (value=="full_icme") p.region_mode=swcme::regions::Mode::FullICME;
    else return error(Code::InvalidValue,assignment,layer,
                      "expected shock_only or full_icme");
    return Status{};
  }
  if (key=="shock.acceleration_mode") {
    const std::string value=canonical_key(assignment.value);
    if (value=="source") p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
    else if (value=="resolved_compression")
      p.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
    else return error(Code::InvalidValue,assignment,layer,
                      "expected source or resolved_compression");
    return Status{};
  }
  if (key=="shock.relative_source_weight_per_area")
    return scalar(assignment,layer,dimensionless,2,false,
                  &p.relative_source_weight_per_area);
  if (key=="geometry.sheath_thickness_1au")
    return scalar(assignment,layer,distance_au,2,true,
                  &p.sheath_thick_AU_at1AU);
  if (key=="geometry.ejecta_thickness_1au")
    return scalar(assignment,layer,distance_au,2,true,
                  &p.ejecta_thick_AU_at1AU);
  if (key=="smoothing.shock_width_1au")
    return scalar(assignment,layer,distance_au,2,true,
                  &p.edge_smooth_shock_AU_at1AU);
  if (key=="smoothing.leading_edge_width_1au")
    return scalar(assignment,layer,distance_au,2,true,
                  &p.edge_smooth_le_AU_at1AU);
  if (key=="smoothing.trailing_edge_width_1au")
    return scalar(assignment,layer,distance_au,2,true,
                  &p.edge_smooth_te_AU_at1AU);
  if (key=="sheath.ramp_power")
    return scalar(assignment,layer,dimensionless,2,false,&p.sheath_ramp_power);
  if (key=="sheath.leading_edge_speed_factor")
    return scalar(assignment,layer,dimensionless,2,false,
                  &p.V_sheath_LE_factor);
  if (key=="ejecta.density_factor")
    return scalar(assignment,layer,dimensionless,2,false,&p.f_ME);
  if (key=="ejecta.speed_factor")
    return scalar(assignment,layer,dimensionless,2,false,&p.V_ME_factor);
  if (key=="event.launch_epoch")
    return scalar(assignment,layer,time_s,4,true,&c->launch_epoch_s);
  if (key=="event.valid_from") {
    *explicit_valid_from=true;
    return scalar(assignment,layer,time_s,4,true,&c->valid_from_s);
  }
  if (key=="event.valid_until")
    return scalar(assignment,layer,time_s,4,true,&c->valid_until_s);
  if (key=="source.particle_mass")
    return scalar(assignment,layer,mass_kg,2,true,
                  &source.spectrum.particle_mass_kg);
  if (key=="source.charge_number")
    return integer(assignment,layer,&source.spectrum.charge_number);
  if (key=="source.energy_min")
    return scalar(assignment,layer,energy_mev,3,true,
                  &source.spectrum.kinetic_energy_min_MeV);
  if (key=="source.energy_max")
    return scalar(assignment,layer,energy_mev,3,true,
                  &source.spectrum.kinetic_energy_max_MeV);
  if (key=="source.reference_energy")
    return scalar(assignment,layer,energy_mev,3,true,
                  &source.spectrum.reference_energy_MeV);
  if (key=="source.injection_efficiency")
    return scalar(assignment,layer,dimensionless,2,false,
                  &source.injection_efficiency);
  if (key=="source.normalization") {
    const std::string value=canonical_key(assignment.value);
    if (value=="relative_only")
      source.spectrum.normalization=sep::NormalizationMode::RelativeOnly;
    else if (value=="reference_differential_intensity")
      source.spectrum.normalization=
          sep::NormalizationMode::ReferenceDifferentialIntensity;
    else return error(Code::InvalidValue,assignment,layer,
                      "expected relative_only or reference_differential_intensity");
    return Status{};
  }
  if (key=="source.reference_intensity_si") {
    const Unit intensity[]={{"si",1.0}};
    return scalar(assignment,layer,intensity,1,true,
                  &source.spectrum.reference_differential_intensity_SI);
  }

  // This compatibility parameter is present in Params but deliberately has no
  // physical effect. Rejecting rather than accepting it prevents an input deck
  // from claiming to control compression when RH physics is authoritative.
  if (key=="sheath.compression_floor")
    return error(Code::UnsupportedField,assignment,layer,
                 "deprecated and ignored by canonical RH physics");

  return error(Code::UnknownKey,assignment,layer,
               "unknown canonical SWCME1D configuration key");
}

inline std::string source_manifest(const SourceConfiguration& source) {
  std::ostringstream out;
  out << std::scientific << std::setprecision(17)
      << "source.particle_mass_kg=" << source.spectrum.particle_mass_kg << '\n'
      << "source.charge_number=" << source.spectrum.charge_number << '\n'
      << "source.energy_min_MeV="
      << source.spectrum.kinetic_energy_min_MeV << '\n'
      << "source.energy_max_MeV="
      << source.spectrum.kinetic_energy_max_MeV << '\n'
      << "source.reference_energy_MeV="
      << source.spectrum.reference_energy_MeV << '\n'
      << "source.normalization="
      << sep::normalization_mode_name(source.spectrum.normalization) << '\n'
      << "source.reference_differential_intensity_SI="
      << source.spectrum.reference_differential_intensity_SI << '\n'
      << "source.injection_efficiency=" << source.injection_efficiency << '\n';
  return out.str();
}

}  // namespace detail

inline ResolveResult Resolve(Preset preset,const std::vector<Layer>& layers) {
  ResolveResult result;
  // Resolve the preset itself with the same layer precedence before applying
  // physical fields. A command-line `preset=...` can therefore override the
  // PARAM choice while every later field still overrides the expanded preset.
  Preset effective_preset=preset;
  for (std::size_t layer_index=0;layer_index<layers.size();++layer_index) {
    for (std::size_t i=0;i<layers[layer_index].assignments.size();++i) {
      const Assignment& assignment=layers[layer_index].assignments[i];
      if (detail::canonical_key(assignment.key)!="preset") continue;
      const std::string value=detail::canonical_key(assignment.value);
      if (value=="fast") effective_preset=Preset::Fast;
      else if (value=="slow") effective_preset=Preset::Slow;
      else {
        result.status=detail::error(
            Code::InvalidValue,assignment,layers[layer_index].name,
            "preset must be fast or slow");
        return result;
      }
    }
  }
  result.configuration.preset=effective_preset;
  result.configuration.model=detail::preset_params(effective_preset);
  // Preserve srcSEP's established energy interval while using the canonical
  // source type and validation. Presets differ only in the SW/CME realization.
  result.configuration.source.spectrum.kinetic_energy_min_MeV=0.1;
  result.configuration.source.spectrum.kinetic_energy_max_MeV=500.0;
  result.configuration.source.spectrum.reference_energy_MeV=10.0;

  bool explicit_valid_from=false;
  for (std::size_t layer_index=0;layer_index<layers.size();++layer_index) {
    const Layer& layer=layers[layer_index];
    std::set<std::string> seen;
    for (std::size_t i=0;i<layer.assignments.size();++i) {
      const Assignment& assignment=layer.assignments[i];
      const std::string key=detail::canonical_key(assignment.key);
      if (!seen.insert(key).second) {
        result.status=detail::error(
            Code::DuplicateKey,assignment,layer.name,
            "key occurs more than once in the same authority layer");
        return result;
      }
      result.status=detail::ApplyAssignment(
          assignment,layer.name,&result.configuration,&explicit_valid_from);
      if (!result.status.ok()) return result;
    }
  }

  if (!explicit_valid_from)
    result.configuration.valid_from_s=result.configuration.launch_epoch_s;
  if (!std::isfinite(result.configuration.launch_epoch_s) ||
      !std::isfinite(result.configuration.valid_from_s) ||
      result.configuration.valid_from_s<result.configuration.launch_epoch_s ||
      std::isnan(result.configuration.valid_until_s) ||
      !(result.configuration.valid_until_s>
        result.configuration.valid_from_s)) {
    Assignment synthetic;
    synthetic.key="event.validity";
    result.status=detail::error(
        Code::InvalidCombination,synthetic,"resolved",
        "require finite launch/valid-from, valid-from >= launch, and valid-until > valid-from");
    return result;
  }
  if (!std::isfinite(result.configuration.source.injection_efficiency) ||
      result.configuration.source.injection_efficiency<0.0 ||
      result.configuration.source.injection_efficiency>1.0) {
    Assignment synthetic;
    synthetic.key="source.injection_efficiency";
    result.status=detail::error(
        Code::InvalidValue,synthetic,"resolved",
        "injection efficiency must be finite and in [0,1]");
    return result;
  }

  const swcme::config::ValidationResult model_validation=
      swcme1d::validate_params(result.configuration.model);
  if (!model_validation.ok()) {
    Assignment synthetic;
    synthetic.key="canonical_model";
    result.status=detail::error(
        Code::CanonicalValidationFailure,synthetic,"resolved",
        model_validation.summary());
    return result;
  }
  const swcme::ModelStatus source_validation=
      swcme::sep::validate_spectrum_config(result.configuration.source.spectrum);
  if (!source_validation.ok()) {
    Assignment synthetic;
    synthetic.key="canonical_source";
    std::ostringstream message;
    message << swcme::status_code_name(source_validation.code)
            << " at " << source_validation.context;
    result.status=detail::error(
        Code::CanonicalValidationFailure,synthetic,"resolved",message.str());
    return result;
  }

  std::ostringstream manifest;
  manifest << "swcme1d_input_schema=1\n"
           << "preset=" << preset_name(effective_preset) << '\n'
           << std::scientific << std::setprecision(17)
           << "event.launch_epoch_s=" << result.configuration.launch_epoch_s << '\n'
           << "event.valid_from_s=" << result.configuration.valid_from_s << '\n'
           << "event.valid_until_s=" << result.configuration.valid_until_s << '\n'
           << "fixed.frame=" << swcme::defaults::FRAME_NAME << '\n'
           << "fixed.solar_rotation_rate_rad_s="
           << swcme::defaults::SOLAR_ROTATION_RATE_RAD_S << '\n'
           << "unsupported.sheath.compression_floor=deprecated-no-effect\n"
           << swcme1d::resolved_configuration_manifest(result.configuration.model)
           << detail::source_manifest(result.configuration.source);
  result.configuration.normalized_manifest=manifest.str();
  result.configuration.fingerprint=
      detail::fingerprint(result.configuration.normalized_manifest);
  result.status=Status{};
  return result;
}

}  // namespace input1d
}  // namespace swcme

#endif  // SWCME1D_INPUT_HPP
