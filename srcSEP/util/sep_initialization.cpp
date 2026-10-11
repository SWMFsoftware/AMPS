#include "sep_initialization.h"

// Generic AMPS runtime-input reader (src/INPUT_FILE.md).  As for the
// SEP-common headers (sep_common_header_path.h), prefer the AMPS-root include
// path that every AMPS compilation context provides; dependency-light test
// builds supply -I<AMPS_ROOT>/src/general/input instead.
#if defined(__has_include)
#  if __has_include(<src/general/input/amps_input.h>)
#    include <src/general/input/amps_input.h>
#  else
#    include "amps_input.h"
#  endif
#else
#  include "amps_input.h"
#endif

// Shared SEP input vocabulary (schema 6) and the canonical SWCME1D resolver,
// both header-only; resolved like the generic reader above.
#include SRCSEP_SEP_COMMON_HEADER(sep_input_layout.h)
#if defined(__has_include)
#  if __has_include(<src/models/swcme/swcme1d_input.hpp>)
#    include <src/models/swcme/swcme1d_input.hpp>
#  else
#    include "swcme1d_input.hpp"
#  endif
#else
#  include "swcme1d_input.hpp"
#endif

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>

namespace SEP {
namespace Initialization {
namespace {

std::unique_ptr<const Configuration> gActive;

Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                  message);
}

std::string Trim(const std::string& value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return std::string();
  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

bool ParseDouble(const std::string& text, double* value) {
  errno = 0;
  char* end = NULL;
  const double parsed = std::strtod(text.c_str(), &end);
  if (text.empty() || end == text.c_str() || *end != '\0' || errno == ERANGE ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParseUnsigned64(const std::string& text, std::uint64_t* value) {
  errno = 0;
  char* end = NULL;
  if (text.empty() || text[0] == '-') return false;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE) return false;
  *value = static_cast<std::uint64_t>(parsed);
  return true;
}

bool ParseBool(const std::string& text, bool* value) {
  const std::string normalized = Lower(text);
  if (normalized == "true" || normalized == "yes" || normalized == "on") {
    *value = true;
    return true;
  }
  if (normalized == "false" || normalized == "no" || normalized == "off") {
    *value = false;
    return true;
  }
  return false;
}

double Dot(const Vec3& left, const Vec3& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}
Vec3 Add(const Vec3& left, const Vec3& right) {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}
Vec3 Subtract(const Vec3& left, const Vec3& right) {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}
Vec3 Scale(double value, const Vec3& vector) {
  return {value * vector.x, value * vector.y, value * vector.z};
}
Vec3 Cross(const Vec3& left, const Vec3& right) {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}
double Norm(const Vec3& vector) { return std::sqrt(Dot(vector, vector)); }
Vec3 Unit(const Vec3& vector) {
  const double length = Norm(vector);
  return length > 0.0 ? Scale(1.0 / length, vector) : Vec3();
}
bool Finite(const Vec3& vector) {
  return std::isfinite(vector.x) && std::isfinite(vector.y) &&
         std::isfinite(vector.z);
}

double Clamp(double value, double lower, double upper) {
  return std::max(lower, std::min(value, upper));
}

double Profile(double coordinate, RefinementProfile profile, double exponent) {
  const double x = Clamp(coordinate, 0.0, 1.0);
  if (profile == RefinementProfile::Linear) return x;
  if (profile == RefinementProfile::PowerLaw) return std::pow(x, exponent);
  const double smooth = x * x * (3.0 - 2.0 * x);
  return std::pow(smooth, exponent);
}

// Complete set of physically effective canonical SWCME1D inputs.  Data knots
// are conditional because an empty list is the only physically correct value
// for ballistic/DBM, whereas data-driven kinematics requires both lists.  The
// deprecated sheath.compression_floor is intentionally absent: canonical RH
// physics rejects it because it has no effect.
const std::set<std::string>& RequiredSwcmeKeys() {
  static const std::set<std::string> keys = {
      "preset",
      "ambient.wind_speed", "ambient.density_1au",
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

// Positive-polarity outward tangent. Magnetic polarity changes B, not this
// geometric curve, so it is deliberately absent from initialization input.
Vec3 ParkerTangent(const Vec3& relativePosition,
                   const Configuration& configuration) {
  const double radius = Norm(relativePosition);
  if (!(radius > 0.0)) return Vec3();
  const Vec3 radial = Scale(1.0 / radius, relativePosition);
  const Vec3 axis = {0.0, 0.0, 1.0};
  const double sourceRadius =
      Norm(Subtract(configuration.parkerInitialPointM,
                    configuration.parkerOriginM));
  const double winding = configuration.solarRotationRateRadPerS *
      std::max(0.0, radius - sourceRadius) /
      configuration.solarWindSpeedMPerS;
  return Unit(Subtract(radial, Scale(winding, Cross(axis, radial))));
}

Vec3 ParkerDirectionAtRadius(double radius,
                             const Configuration& configuration) {
  const Vec3 source = Unit(Subtract(configuration.parkerInitialPointM,
                                    configuration.parkerOriginM));
  const double sourceRadius = configuration.innerRadiusM;
  const double angle = -configuration.solarRotationRateRadPerS *
      std::max(0.0, radius - sourceRadius) /
      configuration.solarWindSpeedMPerS;
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  // Rodrigues rotation about +Z, written explicitly to keep this small layer
  // independent of a general matrix library.
  return {cosine * source.x - sine * source.y,
          sine * source.x + cosine * source.y, source.z};
}

ObserverConfiguration* ObserverForSection(const std::string& section,
                                          Configuration* configuration) {
  const std::string prefix = "observer.";
  if (section.rfind(prefix, 0) != 0 || section.size() == prefix.size())
    return NULL;
  const std::string id = section.substr(prefix.size());
  for (char value : id)
    if (!std::isalnum(static_cast<unsigned char>(value)) && value != '_' &&
        value != '-') return NULL;
  for (ObserverConfiguration& observer : configuration->observers)
    if (observer.id == id) return &observer;
  ObserverConfiguration observer;
  observer.id = id;
  configuration->observers.push_back(observer);
  return &configuration->observers.back();
}

// Parse a value restricted to one accepted spelling.  D3: the shared
// vocabulary may offer values that srcSEP does not implement (srcSEP3D's
// providers, observer kinds, ...); those are typed "not supported" errors
// rather than generic invalid values.
Transport::Status RequireSrcSepValue(const std::string& field,
                                     const std::string& value,
                                     const std::string& accepted,
                                     const std::set<std::string>& shared) {
  const std::string normalized = Lower(value);
  if (normalized == accepted) return Transport::Status::Ok();
  if (shared.count(normalized) != 0)
    return Error(SEP::InputLayout::NotSupportedMessage(
        "'" + field + " = " + value + "'", SEP::InputLayout::kSrcSep) +
        "; srcSEP accepts only '" + accepted + "'");
  return Error("invalid value for '" + field + "': " + value);
}

// Per-walk state that is not part of Configuration itself.
struct ApplyState {
  // source: shock reference spectrum (decision D4): becomes the canonical
  // SWCME source.* assignments once the document is complete.
  SEP::InputLayout::SwcmeSourceSettings spectrum;
};

// Apply one srcSEP key under its qualified name (block path + key).  The
// walker has already checked that the key belongs to the shared vocabulary
// and is implemented by srcSEP.
Transport::Status Apply(const std::string& section, const std::string& key,
                        const std::string& value, Configuration* c,
                        ApplyState* state) {
  const std::string field = section + "." + key;
  auto invalid = [&]() { return Error("invalid value for '" + field + "': " + value); };
  auto boolean = [&](bool* destination) {
    return ParseBool(value, destination) ? Transport::Status::Ok() : invalid();
  };
  auto real = [&](double* destination) {
    return ParseDouble(value, destination) ? Transport::Status::Ok() : invalid();
  };
  auto count = [&](std::uint64_t* destination) {
    return ParseUnsigned64(value, destination) ? Transport::Status::Ok()
                                               : invalid();
  };
  if (field == "application.sep.schema_version") {
    // Decision D1: srcSEP and srcSEP3D share content schema 6.  The srcSEP
    // schemas 1-4 are retired without a compatibility path.
    std::uint64_t parsed = 0;
    if (!ParseUnsigned64(value, &parsed)) return invalid();
    if (parsed != SEP::InputLayout::kSchemaVersion)
      return Error("application.sep.schema_version = " + value + " is "
                   "retired; srcSEP and srcSEP3D read only the shared schema "
                   "6 (convert with tools/convert_amps_input.py "
                   "sep6-from-srcsep)");
  } else if (field == "application.sep.transport") {
    // One shared key selects the transport formulation in both applications.
    // srcSEP maps it onto its three production movers, which implement the
    // same formulations: parker (pitch-angle-averaged Parker transport),
    // focused-diffusion = fte-dmumu (continuous D_mumu pitch-angle diffusion)
    // and focused-scattering = fte-mfp (discrete scattering events at the
    // rate nu = v / lambda_parallel); srcSEP/PRODUCTION_MOVER_API.md.
    const std::string normalized = Lower(value);
    if (normalized == "parker")
      c->particleMover = Mover::ProductionMover::Parker;
    else if (normalized == "focused-diffusion")
      c->particleMover = Mover::ProductionMover::FocusedTransportDiffusion;
    else if (normalized == "focused-scattering")
      c->particleMover = Mover::ProductionMover::FocusedTransportMeanFreePath;
    else
      return Error("invalid value for '" + field + "': " + value +
                   " (accepted: parker, focused-diffusion, focused-scattering)");
    c->particleMoverSpecified = true;
  } else if (field == "run.time_step_s") {
    return real(&c->timeStepS);
  } else if (field == "run.maximum_time_steps") {
    c->maximumTimeStepsSpecified = true;
    return count(&c->maximumTimeSteps);
  } else if (field == "source.shock.enabled") {
    // srcSEP's SWCME shock source is always active; a disabled source is a
    // srcSEP3D propagation-only configuration.
    bool enabled = false;
    if (!ParseBool(value, &enabled)) return invalid();
    if (!enabled)
      return Error(SEP::InputLayout::NotSupportedMessage(
          "'source.shock.enabled = false'", SEP::InputLayout::kSrcSep));
  } else if (field == "source.shock.samples_per_step") {
    return count(&c->macroparticlesPerStep);
  } else if (field == "source.shock.minimum_energy_j") {
    return real(&state->spectrum.minimumEnergyJ);
  } else if (field == "source.shock.maximum_energy_j") {
    return real(&state->spectrum.maximumEnergyJ);
  } else if (field == "source.shock.injection_efficiency") {
    return real(&state->spectrum.injectionEfficiency);
  } else if (field == "source.shock.spectrum_reference_energy_j") {
    return real(&state->spectrum.referenceEnergyJ);
  } else if (field == "source.shock.spectrum_particle_mass_kg") {
    return real(&state->spectrum.particleMassKg);
  } else if (field == "source.shock.spectrum_charge_number") {
    std::uint64_t charge = 0;
    if (!ParseUnsigned64(value, &charge) || charge == 0 ||
        charge > static_cast<std::uint64_t>(std::numeric_limits<long>::max()))
      return invalid();
    state->spectrum.chargeNumber = static_cast<long>(charge);
  } else if (field == "source.shock.spectrum_normalization") {
    // srcSEP's swept-volume source is defined only for the relative SWCME
    // spectrum (SW1DAdapter rejects a calibrated reference intensity), so the
    // shared value reference-differential-intensity is a typed D3 error here
    // rather than a failure deferred to model configuration.
    const std::string normalized = Lower(value);
    if (normalized == "reference-differential-intensity")
      return Error(SEP::InputLayout::NotSupportedMessage(
          "'" + field + " = " + value + "'", SEP::InputLayout::kSrcSep) +
          "; srcSEP accepts only 'relative-only'");
    if (normalized != "relative-only") return invalid();
    state->spectrum.normalization = normalized;
  } else if (field ==
             "source.shock.spectrum_reference_intensity_per_m2_s_sr_j") {
    return real(&state->spectrum.referenceIntensitySi);
  } else if (field == "species.macroparticle_weight") {
    return real(&c->particleWeight);
  } else if (field == "physics.background.provider") {
    // srcSEP's ambient plasma is the standalone SWCME StepState sampled along
    // the Parker field line (Background::Provider::Swcme).
    return RequireSrcSepValue(field, value, "swcme",
        {"analytic-parker", "python-interpolator", "swmf", "runtime-model"});
  } else if (field == "physics.shock.authority") {
    // The particle-facing shock is the canonical SWCME1D front.
    return RequireSrcSepValue(field, value, "swcme", {"none"});
  } else if (field == "output.initialization_mesh_tecplot_file") {
    c->meshTecplotFile = value;
  } else if (field == "output.initialization_field_line_tecplot_file") {
    c->fieldLineTecplotFile = value;
  } else if (field == "output.initialization_data_tecplot_file") {
    c->dataTecplotFile = value;
  } else if (field == "mesh.field_line.origin_x_m") {
    return real(&c->parkerOriginM.x);
  } else if (field == "mesh.field_line.origin_y_m") {
    return real(&c->parkerOriginM.y);
  } else if (field == "mesh.field_line.origin_z_m") {
    return real(&c->parkerOriginM.z);
  } else if (field == "mesh.field_line.initial_x_m") {
    return real(&c->parkerInitialPointM.x);
  } else if (field == "mesh.field_line.initial_y_m") {
    return real(&c->parkerInitialPointM.y);
  } else if (field == "mesh.field_line.initial_z_m") {
    return real(&c->parkerInitialPointM.z);
  } else if (field == "mesh.field_line.length_m") {
    return real(&c->parkerLengthM);
  } else if (field == "mesh.field_line.point_count") {
    return count(&c->parkerPointCount);
  } else if (field == "mesh.domain.inner_radius_m") {
    return real(&c->innerRadiusM);
  } else if (field == "mesh.domain.outer_radius_m") {
    return real(&c->outerRadiusM);
  } else if (field == "mesh.global_cell_size_m") {
    return real(&c->globalCellSizeM);
  } else if (field == "mesh.minimum_cell_size_m") {
    return real(&c->minimumCellSizeM);
  } else if (field == "mesh.maximum_level") {
    std::uint64_t parsed = 0;
    if (!ParseUnsigned64(value, &parsed) ||
        parsed > std::numeric_limits<unsigned>::max()) return invalid();
    c->maximumMeshLevel = static_cast<unsigned>(parsed);
  } else if (field == "mesh.refinement.solar.enabled") {
    return boolean(&c->solarRefinementEnabled);
  } else if (field == "mesh.refinement.solar.surface_cell_size_m") {
    return real(&c->solarSurfaceCellSizeM);
  } else if (field == "mesh.refinement.solar.transition_outer_radius_m") {
    return real(&c->solarTransitionOuterRadiusM);
  } else if (field == "mesh.refinement.solar.profile") {
    const std::string normalized = Lower(value);
    if (normalized == "linear") c->solarProfile = RefinementProfile::Linear;
    else if (normalized == "power-law") c->solarProfile = RefinementProfile::PowerLaw;
    else if (normalized == "smoothstep") c->solarProfile = RefinementProfile::Smoothstep;
    else return invalid();
  } else if (field == "mesh.refinement.solar.exponent") {
    return real(&c->solarExponent);
  } else if (field == "mesh.refinement.tube.enabled") {
    return boolean(&c->tubeRefinementEnabled);
  } else if (field == "mesh.refinement.tube.reference_radius_m") {
    return real(&c->tubeReferenceRadiusM);
  } else if (field == "mesh.refinement.tube.radius_at_reference_m") {
    return real(&c->tubeRadiusAtReferenceM);
  } else if (field == "mesh.refinement.tube.radius_mode") {
    const std::string normalized = Lower(value);
    if (normalized == "physical-constant")
      c->tubeRadiusMode = TubeRadiusMode::PhysicalConstant;
    else if (normalized == "constant-angular-width")
      c->tubeRadiusMode = TubeRadiusMode::ConstantAngularWidth;
    else return invalid();
  } else if (field == "mesh.refinement.tube.center_cell_size_m") {
    return real(&c->tubeCenterCellSizeM);
  } else if (field == "mesh.refinement.tube.transverse_profile") {
    const std::string normalized = Lower(value);
    if (normalized == "linear") c->tubeProfile = RefinementProfile::Linear;
    else if (normalized == "power-law") c->tubeProfile = RefinementProfile::PowerLaw;
    else if (normalized == "smoothstep") c->tubeProfile = RefinementProfile::Smoothstep;
    else return invalid();
  } else if (field == "mesh.refinement.tube.transverse_exponent") {
    return real(&c->tubeExponent);
  } else if (section == "application.sep.turbulence") {
    // srcSEP's content of the only application-specific block (D2): the
    // self-consistent Alfven-wave turbulence configuration.  Spellings are
    // those of the matching command-line options (sep_turbulence_core.cpp).
    Turbulence::Configuration& t = c->turbulence;
    c->turbulenceSpecified = true;
    const std::string normalized = Lower(value);
    if (key == "source") {
      if (!Turbulence::ParseSource(normalized, &t.source)) return invalid();
    } else if (key == "representation") {
      if (normalized == "integrated")
        t.representation = Turbulence::Representation::Integrated;
      else if (normalized == "spectral")
        t.representation = Turbulence::Representation::Spectral;
      else return invalid();
    } else if (key == "coupling") {
      if (!Turbulence::ParseCouplingPolicy(normalized, &t.coupling)) return invalid();
    } else if (key == "advection") {
      return boolean(&t.advectionEnabled);
    } else if (key == "reflection") {
      return boolean(&t.reflectionEnabled);
    } else if (key == "cascade") {
      return boolean(&t.cascadeEnabled);
    } else if (key == "shock_injection") {
      return boolean(&t.shockInjectionEnabled);
    } else if (key == "inner_boundary") {
      if (!Turbulence::ParseBoundaryPolicy(normalized, &t.innerBoundary.policy))
        return invalid();
    } else if (key == "inner_boundary_value_si") {
      return real(&t.innerBoundary.value);
    } else if (key == "outer_boundary") {
      if (!Turbulence::ParseBoundaryPolicy(normalized, &t.outerBoundary.policy))
        return invalid();
    } else if (key == "outer_boundary_value_si") {
      return real(&t.outerBoundary.value);
    } else if (key == "reflection_coefficient") {
      return real(&t.reflectionCoefficient);
    } else if (key == "cascade_coefficient") {
      return real(&t.cascadeCoefficient);
    } else if (key == "perpendicular_correlation_length_m") {
      return real(&t.perpendicularCorrelationLengthM);
    } else if (key == "electron_heating_fraction") {
      return real(&t.electronHeatingFraction);
    } else if (key == "spectral_k_min_per_m") {
      return real(&t.spectralKMinPerM);
    } else if (key == "spectral_k_max_per_m") {
      return real(&t.spectralKMaxPerM);
    } else if (key == "spectral_bins") {
      std::uint64_t bins = 0;
      if (!ParseUnsigned64(value, &bins) || bins == 0) return invalid();
      t.spectralBins = static_cast<std::size_t>(bins);
    } else if (key == "cfl_safety") {
      return real(&t.cflSafety);
    } else if (key == "operator_accuracy_safety") {
      return real(&t.operatorAccuracySafety);
    } else if (key == "maximum_source_fraction") {
      return real(&t.maximumSourceFraction);
    } else if (key == "maximum_cascade_fraction") {
      return real(&t.maximumCascadeFraction);
    } else if (key == "minimum_substep_s") {
      return real(&t.minimumSubstepS);
    } else if (key == "maximum_substeps") {
      if (!ParseUnsigned64(value, &t.maximumSubsteps) || t.maximumSubsteps == 0)
        return invalid();
    } else if (key == "conservation_relative_tolerance") {
      return real(&t.conservationRelativeTolerance);
    } else if (key == "shock_source_efficiency") {
      return real(&c->shockSourceEfficiency);
    } else if (key == "shock_source_plus_fraction") {
      return real(&c->shockSourcePlusFraction);
    } else if (key == "diagnostic_cadence_steps") {
      return count(&c->diagnosticCadenceSteps);
    } else {
      return Error("unknown configuration key '" + field + "'");
    }
  } else {
    ObserverConfiguration* observer = ObserverForSection(section, c);
    if (observer == NULL)
      return Error("unknown configuration key '" + field + "'");
    if (key == "kind") {
      // D5: srcSEP's observers are points on the field line.
      return RequireSrcSepValue(field, value, "field-line-radius",
          {"fixed-cartesian", "fixed-heliographic", "moving-cartesian",
           "spherical-shell", "field-connected"});
    } else if (key == "heliocentric_radius_m") {
      return real(&observer->heliocentricRadiusM);
    } else if (key == "minimum_energy_j") {
      return real(&observer->minimumEnergyJ);
    } else if (key == "maximum_energy_j") {
      return real(&observer->maximumEnergyJ);
    } else if (key == "energy_bins") {
      std::uint64_t parsed = 0;
      if (!ParseUnsigned64(value, &parsed) ||
          parsed > std::numeric_limits<unsigned>::max()) return invalid();
      observer->energyChannels = static_cast<unsigned>(parsed);
    } else if (key == "energy_spacing") {
      const std::string normalized = Lower(value);
      if (normalized == "logarithmic")
        observer->energySpacing = EnergyChannelSpacing::Logarithmic;
      else if (normalized == "linear")
        observer->energySpacing = EnergyChannelSpacing::Linear;
      else return invalid();
    } else if (key == "pitch_angle_bins") {
      std::uint64_t parsed = 0;
      if (!ParseUnsigned64(value, &parsed) ||
          parsed > std::numeric_limits<unsigned>::max()) return invalid();
      observer->pitchAngleBins = static_cast<unsigned>(parsed);
    } else {
      return Error("unknown configuration key '" + field + "'");
    }
  }
  return Transport::Status::Ok();
}

}  // namespace

const char* Name(RefinementProfile profile) {
  if (profile == RefinementProfile::Linear) return "linear";
  if (profile == RefinementProfile::PowerLaw) return "power-law";
  return "smoothstep";
}
const char* Name(TubeRadiusMode mode) {
  return mode == TubeRadiusMode::PhysicalConstant
      ? "physical-constant" : "constant-angular-width";
}
const char* Name(EnergyChannelSpacing spacing) {
  return spacing == EnergyChannelSpacing::Logarithmic
      ? "logarithmic" : "linear";
}

Transport::Status Validate(const Configuration& c) {
  if (!Finite(c.parkerOriginM) || !Finite(c.parkerInitialPointM))
    return Error("Parker origin or initial point is not finite");
  const double values[] = {
      c.parkerOriginM.x, c.parkerOriginM.y, c.parkerOriginM.z,
      c.parkerInitialPointM.x, c.parkerInitialPointM.y,
      c.parkerInitialPointM.z, c.parkerLengthM, c.solarWindSpeedMPerS,
      c.solarRotationRateRadPerS, c.innerRadiusM, c.outerRadiusM,
      c.globalCellSizeM, c.minimumCellSizeM, c.solarSurfaceCellSizeM,
      c.solarTransitionOuterRadiusM, c.solarExponent,
      c.tubeReferenceRadiusM, c.tubeRadiusAtReferenceM,
      c.tubeCenterCellSizeM, c.tubeExponent};
  for (double value : values)
    if (!std::isfinite(value)) return Error("configuration contains a non-finite value");
  if (c.parkerPointCount < 2 ||
      c.parkerPointCount > 10000000ULL || c.parkerLengthM <= 0.0 ||
      c.solarWindSpeedMPerS <= 0.0 || c.solarRotationRateRadPerS < 0.0 ||
      c.innerRadiusM <= 0.0 || c.outerRadiusM <= c.innerRadiusM ||
      c.minimumCellSizeM <= 0.0 || c.globalCellSizeM < c.minimumCellSizeM ||
      c.maximumMeshLevel > 19 || c.solarExponent <= 0.0 ||
      c.tubeExponent <= 0.0) {
    return Error("Parker, domain, or global mesh values are outside their supported range");
  }
  if (!std::isfinite(c.timeStepS) || c.timeStepS <= 0.0 ||
      c.macroparticlesPerStep == 0 ||
      c.macroparticlesPerStep >
          static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
      !std::isfinite(c.particleWeight) || c.particleWeight <= 0.0 ||
      c.meshTecplotFile.empty() || c.fieldLineTecplotFile.empty() ||
      c.dataTecplotFile.empty()) {
    return Error("time step, particle sampling, or Tecplot output "
                 "configuration is invalid");
  }
  if (c.maximumTimeStepsSpecified &&
      (c.maximumTimeSteps == 0 ||
       c.maximumTimeSteps >
           static_cast<std::uint64_t>(std::numeric_limits<int>::max())))
    return Error("run.maximum_time_steps must be in [1, 2147483647]");
  if (c.turbulenceSpecified) {
    const Transport::Status turbulence =
        Turbulence::ValidateConfiguration(c.turbulence);
    if (!turbulence.ok())
      return Error("application.sep.turbulence: " + turbulence.message);
    if (!(c.shockSourceEfficiency >= 0.0 && c.shockSourceEfficiency <= 1.0) ||
        !(c.shockSourcePlusFraction >= 0.0 && c.shockSourcePlusFraction <= 1.0) ||
        c.diagnosticCadenceSteps >
            static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
      return Error("application.sep.turbulence shock fractions must be in "
                   "[0,1] and diagnostic_cadence_steps must fit an int");
  }
  if (c.observers.empty())
    return Error("at least one '#begin observer: ID' block is required");
  std::set<std::string> observerIds;
  for (const ObserverConfiguration& observer : c.observers) {
    if (observer.id.empty() || observer.id.size() > 48 ||
        !observerIds.insert(observer.id).second ||
        !std::isfinite(observer.heliocentricRadiusM) ||
        observer.heliocentricRadiusM < c.innerRadiusM ||
        observer.heliocentricRadiusM > c.outerRadiusM ||
        !std::isfinite(observer.minimumEnergyJ) ||
        !std::isfinite(observer.maximumEnergyJ) ||
        observer.minimumEnergyJ <= 0.0 ||
        observer.maximumEnergyJ <= observer.minimumEnergyJ ||
        observer.energyChannels == 0 || observer.energyChannels > 1000000U ||
        observer.pitchAngleBins == 0 || observer.pitchAngleBins > 1000000U)
      return Error("one-dimensional observer identity, location, or spectrum is invalid");
  }
  const double sourceRadius = Norm(Subtract(c.parkerInitialPointM, c.parkerOriginM));
  if (std::fabs(sourceRadius - c.innerRadiusM) > 1.0e-10 * c.innerRadiusM)
    return Error("Parker initial point must lie on mesh.domain.inner_radius_m");
  if (c.solarRefinementEnabled &&
      (c.solarSurfaceCellSizeM < c.minimumCellSizeM ||
       c.solarSurfaceCellSizeM > c.globalCellSizeM ||
       c.solarTransitionOuterRadiusM <= c.innerRadiusM ||
       c.solarTransitionOuterRadiusM > c.outerRadiusM)) {
    return Error("near-Sun refinement sizes/radii are inconsistent");
  }
  if (c.tubeRefinementEnabled &&
      (c.tubeReferenceRadiusM <= c.innerRadiusM ||
       c.tubeRadiusAtReferenceM <= 0.0 ||
       c.tubeCenterCellSizeM < c.minimumCellSizeM ||
       c.tubeCenterCellSizeM > c.globalCellSizeM)) {
    return Error("magnetic-tube refinement sizes/radii are inconsistent");
  }
  return Transport::Status::Ok();
}

// ---------------------------------------------------------------------------
// Generic AMPS runtime input -> srcSEP configuration
// ---------------------------------------------------------------------------
//
// The srcSEP input is written in the generic AMPS runtime-input format
// (src/INPUT_FILE.md; reader src/general/input/amps_input.h) using the SEP
// content schema 6 shared with srcSEP3D.  The complete block and key
// vocabulary is the table in src/models/sep_common/sep_input_layout.h; the
// walker consults it for every block and key:
//
//   * a block or key that is not in the table is unknown and fatal;
//   * a srcSEP3D-only block or key (mesh.active_region, run.campaign_seed,
//     an observer position, ...) is a typed "not supported by srcSEP" error
//     (decision D3) -- it is never ignored;
//   * everything else reaches Apply() under its qualified name.
//
// Qualified names.  A key's qualified name is the dot-joined path of its
// enclosing blocks, where '#begin kind: name' contributes 'kind.name':
// 'mesh.domain.inner_radius_m', 'mesh.refinement.solar.profile',
// 'application.sep.turbulence.cascade', 'observer.earth.energy_bins'.
//
// Model blocks keep the vocabularies of the libraries that own them:
//   physics: swcme              canonical SWCME1D assignments (dotted keys,
//                               units in the value), forwarded by main.cpp to
//                               the model-owned resolver; its source.* keys
//                               are synthesized from source: shock (D4) and
//                               its 3-D-only keys are typed D3 errors;
//   physics: parallel_diffusion forwarded to SEP::ParallelDiffusion::
//                               ParseSection by main.cpp.
namespace {

namespace Layout = SEP::InputLayout;

struct WalkState {
  Configuration* configuration = nullptr;
  std::set<std::string> sections;  // qualified block paths
  std::set<std::string> assigned;  // qualified application keys
  std::map<std::string, std::size_t> keyLines;  // qualified key -> line
  ApplyState apply;
};

// The 3-D SWCME rotation-rate key is honoured by srcSEP only when it states
// the fixed SWCME1D Carrington rate that srcSEP uses anyway (D3: a setting is
// rejected only when the application cannot honour it).  The value is the
// plain SI number [rad/s] that the SWCME3D resolver reads.
bool IsSwcme1dRotationRate(const std::string& value) {
  double parsed = 0.0;
  if (!ParseDouble(Trim(value), &parsed)) return false;
  const double fixed = swcme::constants::SOLAR_ROTATION_RAD_S;
  return std::fabs(parsed - fixed) <=
      64.0 * std::numeric_limits<double>::epsilon() * std::fabs(fixed);
}

Transport::Status WalkBlocks(const AMPS::Input::Block& parent,
                             const std::string& parentPath, WalkState* state) {
  for (const AMPS::Input::Block& block : parent.blocks) {
    const Layout::BlockRule* rule =
        Layout::FindBlock(parentPath, block.kind, block.name);
    if (rule == NULL)
      return Error(AMPS::Input::ErrorAt(block.begin, "block '" + block.Label() +
          "' is not part of the SEP input" +
          (parentPath.empty() ? std::string(" at the top level")
                              : " inside '" + parentPath + "'") +
          " (allowed: " + Layout::AllowedChildren(parentPath) + ")").message);
    const std::string path =
        Layout::QualifiedPath(parentPath, block.kind, block.name);
    if ((rule->support & Layout::kSrcSep) == 0)
      return Error(AMPS::Input::ErrorAt(block.begin,
          Layout::NotSupportedMessage("block '" + path + "'",
                                      Layout::kSrcSep)).message);
    if (!state->sections.insert(path).second)
      return Error(AMPS::Input::ErrorAt(block.begin, "duplicate block '" +
                                        path + "'").message);

    if (path == "physics.swcme") {
      // SWCME owns its dotted, unit-bearing keys; nested blocks are the same
      // key space.  main.cpp resolves them with the model-owned resolver.
      std::set<std::string> seen;
      for (const AMPS::Input::Assignment& item : AMPS::Input::Flatten(block)) {
        if (!seen.insert(item.key).second)
          return Error(AMPS::Input::ErrorAt(item.where, "duplicate SWCME key '" +
                                            item.key + "'", item.text).message);
        if (const char* owner = Layout::SourceShockKeyForSwcmeKey(item.key))
          return Error(AMPS::Input::ErrorAt(item.where, "SWCME key '" +
              item.key + "' is written in source: shock as '" + owner +
              "' (each quantity is set in one place, src/INPUT_FILE.md "
              "section 5)", item.text).message);
        if (Layout::IsSwcmeThreeDimensionalKey(item.key)) {
          if (item.key == "parker.solar_rotation_rate_rad_per_s" &&
              IsSwcme1dRotationRate(item.value))
            continue;  // states the fixed SWCME1D rate; nothing to forward
          return Error(AMPS::Input::ErrorAt(item.where,
              Layout::NotSupportedMessage("SWCME key '" + item.key + "'",
                                          Layout::kSrcSep) +
              "; srcSEP's SWCME1D model is spherically symmetric about the "
              "field line with the fixed rotation rate " +
              Layout::ExactDecimal(swcme::constants::SOLAR_ROTATION_RAD_S) +
              " rad/s", item.text).message);
        }
        SwcmeAssignment assignment;
        assignment.key = item.key;
        assignment.value = item.value;
        assignment.line = item.where.line;
        state->configuration->swcmeAssignments.push_back(assignment);
      }
      continue;
    }
    if (path == "physics.parallel_diffusion") {
      // The shared library parser owns this block's model selection and key
      // schema; it receives "key = value" records with their line numbers.
      if (!block.blocks.empty())
        return Error(AMPS::Input::ErrorAt(block.blocks.front().begin,
            "physics: parallel_diffusion holds only key = value records").message);
      state->configuration->hasParallelDiffusionSection = true;
      for (const AMPS::Input::Assignment& item : block.assignments) {
        ParallelDiffusionLine record;
        record.line = item.where.line;
        record.text = item.key + " = " + item.value;
        state->configuration->parallelDiffusionLines.push_back(record);
      }
      continue;
    }

    for (const AMPS::Input::Assignment& item : block.assignments) {
      if (item.key.find('.') != std::string::npos)
        return Error(AMPS::Input::ErrorAt(item.where, "key '" + item.key +
            "' contains '.'; SEP keys are plain names and a nested path "
            "is written as nested blocks", item.text).message);
      const std::string qualified = path + "." + item.key;
      const Layout::KeyRule* key = Layout::FindKey(path, item.key);
      if (key == NULL)
        return Error(AMPS::Input::ErrorAt(item.where,
            "unknown configuration key '" + qualified + "'", item.text).message);
      if ((key->support & Layout::kSrcSep) == 0)
        return Error(AMPS::Input::ErrorAt(item.where,
            Layout::NotSupportedMessage("'" + qualified + "'",
                                        Layout::kSrcSep), item.text).message);
      state->assigned.insert(qualified);
      state->keyLines[qualified] = item.where.line;
      const Transport::Status status = Apply(path, item.key, item.value,
          state->configuration, &state->apply);
      if (!status.ok())
        return Error(AMPS::Input::ErrorAt(item.where, status.message,
                                          item.text).message);
      if (qualified == "application.sep.transport")
        state->configuration->particleMoverLine = item.where.line;
    }
    const Transport::Status nested = WalkBlocks(block, path, state);
    if (!nested.ok()) return nested;
  }
  return Transport::Status::Ok();
}

}  // namespace

static Transport::Status ParseDocument(const AMPS::Input::Document& document,
                                       Configuration* result) {
  if (result == NULL) return Error("configuration output is null");
  Configuration candidate;
  WalkState walk;
  walk.configuration = &candidate;
  {
    const Transport::Status status = WalkBlocks(document.root, "", &walk);
    if (!status.ok()) return status;
  }
  const std::set<std::string>& sections = walk.sections;
  const std::set<std::string>& assigned = walk.assigned;
  if (assigned.count("application.sep.schema_version") == 0)
    return Error("missing required key application.sep.schema_version "
                 "(write schema_version = 6 in '#begin application: sep')");

  // ---- Schema-6 completeness ----------------------------------------------
  //
  // Every block srcSEP implements is required, and every srcSEP key of those
  // blocks is required (the vocabulary table is the list), so no value is
  // silently inherited from a production default.  The single optional key
  // is run.maximum_time_steps, whose absence keeps the historical driver
  // step budget.
  const char* requiredSections[] = {
      "run", "application.sep", "application.sep.turbulence", "mesh",
      "mesh.domain", "mesh.field_line", "mesh.refinement.solar",
      "mesh.refinement.tube", "physics.background", "physics.shock",
      "physics.swcme", "source.shock", "species", "output"};
  for (const char* required : requiredSections)
    if (sections.count(required) == 0)
      return Error(std::string("missing required block '") + required + "'");
  for (const Layout::KeyRule& rule : Layout::kKeys) {
    if ((rule.support & Layout::kSrcSep) == 0) continue;
    const std::string block = rule.block;
    const std::string qualified = block + "." + rule.key;
    if (block == "observer.*" || qualified == "run.maximum_time_steps") continue;
    if (assigned.count(qualified) == 0)
      return Error("missing required key '" + qualified + "'");
  }
  if (candidate.observers.empty())
    return Error("at least one '#begin observer: ID' block is required");
  for (const ObserverConfiguration& observer : candidate.observers)
    for (const Layout::KeyRule& rule : Layout::kKeys) {
      if (std::string(rule.block) != "observer.*" ||
          (rule.support & Layout::kSrcSep) == 0) continue;
      const std::string qualified = "observer." + observer.id + "." + rule.key;
      if (assigned.count(qualified) == 0)
        return Error("missing required key '" + qualified + "'");
    }

  // ---- Decision D4: synthesize the SWCME source from source: shock --------
  //
  // The SWCME source.* assignments are generated from the SI source: shock
  // values and carry the line of the key that supplies each of them, so the
  // canonical resolver's diagnostics point at an editable line.
  struct SourceLine { const char* swcme; const char* sourceKey; };
  const SourceLine sourceLines[] = {
      {"source.particle_mass", "spectrum_particle_mass_kg"},
      {"source.charge_number", "spectrum_charge_number"},
      {"source.energy_min", "minimum_energy_j"},
      {"source.energy_max", "maximum_energy_j"},
      {"source.reference_energy", "spectrum_reference_energy_j"},
      {"source.injection_efficiency", "injection_efficiency"},
      {"source.normalization", "spectrum_normalization"},
      {"source.reference_intensity_si",
       "spectrum_reference_intensity_per_m2_s_sr_j"}};
  for (const Layout::SynthesizedAssignment& synthesized :
       Layout::SynthesizeSwcmeSource(walk.apply.spectrum)) {
    SwcmeAssignment assignment;
    assignment.key = synthesized.key;
    assignment.value = synthesized.value;
    for (const SourceLine& line : sourceLines)
      if (synthesized.key == line.swcme)
        assignment.line =
            walk.keyLines[std::string("source.shock.") + line.sourceKey];
    candidate.swcmeAssignments.push_back(assignment);
  }

  // The complete SWCME1D key set: every physically effective key is
  // explicit, and the data knots exist exactly for data-driven kinematics.
  std::set<std::string> swcmeKeys;
  std::string kinematics;
  for (const SwcmeAssignment& assignment : candidate.swcmeAssignments) {
    swcmeKeys.insert(assignment.key);
    if (assignment.key == "cme.kinematics")
      kinematics = Lower(Trim(assignment.value));
    if (RequiredSwcmeKeys().count(assignment.key) == 0 &&
        assignment.key != "cme.data_times" &&
        assignment.key != "cme.data_radii")
      return Error("line " + std::to_string(assignment.line) +
                   ": unknown SWCME configuration key '" + assignment.key + "'");
  }
  for (const std::string& required : RequiredSwcmeKeys())
    if (swcmeKeys.count(required) == 0)
      return Error("missing required SWCME key '" + required + "'");
  const bool hasDataTimes = swcmeKeys.count("cme.data_times") != 0;
  const bool hasDataRadii = swcmeKeys.count("cme.data_radii") != 0;
  if (kinematics == "data_driven") {
    if (!hasDataTimes || !hasDataRadii)
      return Error("data_driven CME kinematics requires cme.data_times and "
                   "cme.data_radii");
  } else if (hasDataTimes || hasDataRadii) {
    return Error("CME data knots are legal only for data_driven kinematics");
  }

  // ---- Decision D4: the Parker winding takes the SWCME ambient ------------
  //
  // Resolve the file's SWCME layer with the canonical SWCME1D resolver (the
  // same resolver main.cpp calls through SW1DAdapter::Configure) before AMPS
  // exists: it validates every model value and unit now, and it supplies the
  // ambient wind speed V_sw [m/s] of the Parker winding angle
  // phi = -Omega (r - r_s) / V_sw used by the mesh refinement tube and the
  // field-line generator.  Omega is the fixed SWCME1D Carrington rate.  A
  // later --swcme-override of the wind is rejected by main.cpp's agreement
  // check against the final SWCME summary.
  {
    swcme::input1d::Layer layer;
    layer.name = "input file";
    for (const SwcmeAssignment& raw : candidate.swcmeAssignments) {
      swcme::input1d::Assignment assignment;
      assignment.key = raw.key;
      assignment.value = raw.value;
      assignment.origin = "srcSEP input";
      assignment.line = raw.line;
      layer.assignments.push_back(assignment);
    }
    const swcme::input1d::ResolveResult resolved =
        swcme::input1d::Resolve(swcme::input1d::Preset::Fast, {layer});
    if (!resolved.ok()) {
      std::ostringstream message;
      message << "invalid physics: swcme configuration key='"
              << resolved.status.key << "'";
      if (resolved.status.line != 0) message << " line=" << resolved.status.line;
      if (!resolved.status.message.empty())
        message << ": " << resolved.status.message;
      return Error(message.str());
    }
    candidate.solarWindSpeedMPerS = resolved.configuration.model.V_sw_kms * 1.0e3;
    candidate.solarRotationRateRadPerS = swcme::constants::SOLAR_ROTATION_RAD_S;
  }

  const Transport::Status valid = Validate(candidate);
  if (!valid.ok()) return valid;
  *result = candidate;
  return Transport::Status::Ok();
}

Transport::Status ParseText(const std::string& text, Configuration* result) {
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed =
      AMPS::Input::ParseText(text, "<srcSEP input text>", ".", &document);
  if (!parsed.ok()) return Error(parsed.message);
  return ParseDocument(document, result);
}

Transport::Status LoadFile(const std::string& path, Configuration* result) {
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed = AMPS::Input::ParseFile(path, &document);
  if (!parsed.ok()) return Error(parsed.message);
  return ParseDocument(document, result);
}

Transport::Status ApplyOutputDirectoryOverride(
    const std::string& directory, Configuration* configuration) {
  if (configuration == NULL) return Error("initialization configuration is null");
  if (directory.empty()) return Error("initialization output directory is empty");

  // Preserve the reviewed product names and remove only their former parent.
  // Both separator spellings are recognized because input decks are commonly
  // copied between workstation and HPC filesystems before a production run.
  const auto leafName = [](const std::string& path) {
    const std::size_t separator = path.find_last_of("/\\");
    return separator == std::string::npos ? path : path.substr(separator + 1);
  };
  const std::string meshLeaf = leafName(configuration->meshTecplotFile);
  const std::string lineLeaf = leafName(configuration->fieldLineTecplotFile);
  const std::string dataLeaf = leafName(configuration->dataTecplotFile);
  if (meshLeaf.empty() || lineLeaf.empty() || dataLeaf.empty() ||
      meshLeaf == "." || meshLeaf == ".." || lineLeaf == "." ||
      lineLeaf == ".." || dataLeaf == "." || dataLeaf == "..") {
    return Error("initialization Tecplot paths must end in file names before "
                 "the output-directory override can be applied");
  }
  const std::string separator = directory.back() == '/' ? "" : "/";
  configuration->meshTecplotFile = directory + separator + meshLeaf;
  configuration->fieldLineTecplotFile = directory + separator + lineLeaf;
  configuration->dataTecplotFile = directory + separator + dataLeaf;
  return Validate(*configuration);
}

Transport::Status Install(const Configuration& configuration) {
  if (gActive.get() != NULL) return Error("initialization configuration is already installed");
  const Transport::Status valid = Validate(configuration);
  if (!valid.ok()) return valid;
  gActive.reset(new Configuration(configuration));
  return Transport::Status::Ok();
}

bool HasActive() { return gActive.get() != NULL; }
const Configuration& Active() { return *gActive; }

std::string Fingerprint(const Configuration& c) {
  // Canonical text of every startup value of the schema-6 record, hashed with
  // FNV-1a below.  Comments, spacing, key order and include structure never
  // enter: the text is built from parsed values, and the SWCME assignments
  // are sorted.  The identity tag names the content schema so a schema-6
  // identity can never collide with a retired schema-1..4 identity.
  std::ostringstream canonical;
  canonical << std::setprecision(17) << std::scientific
      << "srcsep-initialization-v" << SEP::InputLayout::kSchemaVersion
      << ";dt=" << c.timeStepS
      << ";macro_per_step=" << c.macroparticlesPerStep
      << ";particle_weight=" << c.particleWeight
      << ";mesh_output=" << c.meshTecplotFile
      << ";line_output=" << c.fieldLineTecplotFile
      << ";data_output=" << c.dataTecplotFile;
  if (c.maximumTimeStepsSpecified)
    canonical << ";maximum_time_steps=" << c.maximumTimeSteps;
  for (const ObserverConfiguration& observer : c.observers)
    canonical << ";observer=" << observer.id << ','
        << observer.heliocentricRadiusM << ',' << observer.minimumEnergyJ
        << ',' << observer.maximumEnergyJ << ',' << observer.energyChannels
        << ',' << Name(observer.energySpacing) << ','
        << observer.pitchAngleBins;

  canonical << ";origin=" << c.parkerOriginM.x << ','
      << c.parkerOriginM.y << ',' << c.parkerOriginM.z << ";initial="
      << c.parkerInitialPointM.x << ',' << c.parkerInitialPointM.y << ','
      << c.parkerInitialPointM.z << ";length=" << c.parkerLengthM
      << ";points=" << c.parkerPointCount << ";wind=" << c.solarWindSpeedMPerS
      << ";omega=" << c.solarRotationRateRadPerS << ";inner=" << c.innerRadiusM
      << ";outer=" << c.outerRadiusM << ";global=" << c.globalCellSizeM
      << ";minimum=" << c.minimumCellSizeM << ";level=" << c.maximumMeshLevel
      << ";solar=" << c.solarRefinementEnabled << ',' << c.solarSurfaceCellSizeM
      << ',' << c.solarTransitionOuterRadiusM << ',' << Name(c.solarProfile)
      << ',' << c.solarExponent << ";tube=" << c.tubeRefinementEnabled << ','
      << c.tubeReferenceRadiusM << ',' << c.tubeRadiusAtReferenceM << ','
      << Name(c.tubeRadiusMode) << ',' << c.tubeCenterCellSizeM << ','
      << Name(c.tubeProfile) << ',' << c.tubeExponent;

  // Field order has no physical meaning.  Sort a private copy so two
  // complete decks with the same SWCME assignments have the same startup
  // fingerprint.  The canonical SWCME resolver separately normalizes units
  // and produces the model-owned physical fingerprint.
  std::vector<SwcmeAssignment> assignments = c.swcmeAssignments;
  std::sort(assignments.begin(), assignments.end(),
      [](const SwcmeAssignment& left, const SwcmeAssignment& right) {
        if (left.key != right.key) return left.key < right.key;
        return left.value < right.value;
      });
  for (const SwcmeAssignment& assignment : assignments)
    canonical << ";swcme." << assignment.key << '=' << assignment.value;
  if (c.particleMoverSpecified) {
    // The file-selected mover is part of the startup contract.  The
    // effective mover (after any --particle-mover override) is recorded
    // separately in the WP30 run-configuration fingerprint.
    canonical << ";particle_mover=" << Mover::Describe(c.particleMover).canonicalName;
  }
  if (c.turbulenceSpecified) {
    // File turbulence values (before command-line overrides; the effective
    // turbulence configuration enters the WP30 run fingerprint).
    const Turbulence::Configuration& t = c.turbulence;
    canonical << ";turbulence=" << Turbulence::SourceName(t.source) << ','
        << Turbulence::RepresentationName(t.representation) << ','
        << Turbulence::CouplingPolicyName(t.coupling) << ','
        << t.advectionEnabled << t.reflectionEnabled << t.cascadeEnabled
        << t.shockInjectionEnabled << ','
        << Turbulence::BoundaryPolicyName(t.innerBoundary.policy) << ','
        << t.innerBoundary.value << ','
        << Turbulence::BoundaryPolicyName(t.outerBoundary.policy) << ','
        << t.outerBoundary.value << ',' << t.reflectionCoefficient << ','
        << t.cascadeCoefficient << ',' << t.perpendicularCorrelationLengthM
        << ',' << t.electronHeatingFraction << ',' << t.spectralKMinPerM
        << ',' << t.spectralKMaxPerM << ',' << t.spectralBins << ','
        << t.cflSafety << ',' << t.operatorAccuracySafety << ','
        << t.maximumSourceFraction << ',' << t.maximumCascadeFraction << ','
        << t.minimumSubstepS << ',' << t.maximumSubsteps << ','
        << t.conservationRelativeTolerance << ',' << c.shockSourceEfficiency
        << ',' << c.shockSourcePlusFraction << ',' << c.diagnosticCadenceSteps;
  }
  if (c.hasParallelDiffusionSection) {
    // The physical identity of physics: parallel_diffusion is the shared library's
    // SHA-256 configuration fingerprint, copied here by the srcSEP adapter
    // after the library parser accepted the section (main.cpp, before
    // Install()).  It is independent of comments, spacing, and key order.
    // Until then the fields are empty and the identity is marked unresolved.
    canonical << ";parallel_diffusion_model="
        << (c.parallelDiffusionModelId.empty() ? "unresolved"
                                               : c.parallelDiffusionModelId)
        << ";parallel_diffusion_fingerprint="
        << c.parallelDiffusionConfigurationFingerprint;
  }
  // FNV-1a is used as a compact deterministic identity, not as a security
  // primitive. The canonical manifest above retains full scientific meaning.
  std::uint64_t hash = 1469598103934665603ULL;
  const std::string bytes = canonical.str();
  for (unsigned char byte : bytes) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  std::ostringstream result;
  result << std::hex << std::setfill('0') << std::setw(16) << hash;
  return result.str();
}

Transport::Status BuildParkerLine(const Configuration& c,
                                  std::vector<Vec3>* points) {
  if (points == NULL) return Error("Parker line output is null");
  const Transport::Status valid = Validate(c);
  if (!valid.ok()) return valid;
  std::vector<Vec3> candidate;
  candidate.reserve(static_cast<std::size_t>(c.parkerPointCount));
  Vec3 point = c.parkerInitialPointM;
  candidate.push_back(point);
  const double step = c.parkerLengthM /
      static_cast<double>(c.parkerPointCount - 1);
  for (std::uint64_t index = 1; index < c.parkerPointCount; ++index) {
    const Vec3 relative = Subtract(point, c.parkerOriginM);
    const Vec3 first = ParkerTangent(relative, c);
    const Vec3 midpoint = Add(relative, Scale(0.5 * step, first));
    const Vec3 tangent = ParkerTangent(midpoint, c);
    if (!(Norm(first) > 0.0) || !(Norm(tangent) > 0.0))
      return Error("Parker tangent vanished during line generation");
    point = Add(point, Scale(step, tangent));
    candidate.push_back(point);
  }
  points->swap(candidate);
  return Transport::Status::Ok();
}

Transport::Status WriteParkerLineTecplot(const Configuration& c,
                                         const std::string& path) {
  if (path.empty()) return Error("Parker-line Tecplot path is empty");
  std::vector<Vec3> points;
  const Transport::Status built = BuildParkerLine(c, &points);
  if (!built.ok()) return built;

  // Open only after all geometry succeeds.  A stream failure is returned to
  // the caller and is fatal during initialization; visualization is part of
  // the declared startup contract rather than a best-effort diagnostic.
  std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good())
    return Error("cannot open Parker-line Tecplot file '" + path + "'");
  output << "TITLE=\"srcSEP initialized Parker field line\"\n"
         << "VARIABLES=\"arc_length_m\",\"x_m\",\"y_m\",\"z_m\","
            "\"heliocentric_radius_m\",\"requested_cell_size_m\"\n"
         << "ZONE T=\"field-line\", I=" << points.size()
         << ", F=POINT\n" << std::scientific << std::setprecision(17);
  double arcLength = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i != 0) arcLength += Norm(Subtract(points[i], points[i - 1]));
    output << arcLength << ' ' << points[i].x << ' ' << points[i].y << ' '
           << points[i].z << ' '
           << Norm(Subtract(points[i], c.parkerOriginM)) << ' '
           << RequestedCellSizeM(points[i], c) << '\n';
  }
  output.flush();
  if (!output.good())
    return Error("failed while writing Parker-line Tecplot file '" + path + "'");
  return Transport::Status::Ok();
}

double TubeRadiusM(double radius, const Configuration& c) {
  if (!(radius > 0.0) || !(c.tubeReferenceRadiusM > 0.0))
    return std::numeric_limits<double>::quiet_NaN();
  return c.tubeRadiusMode == TubeRadiusMode::PhysicalConstant
      ? c.tubeRadiusAtReferenceM
      : c.tubeRadiusAtReferenceM * radius / c.tubeReferenceRadiusM;
}

double TubeDistanceM(const Vec3& position, const Configuration& c) {
  const Vec3 relative = Subtract(position, c.parkerOriginM);
  const double radius = Norm(relative);
  if (!(radius > 0.0)) return std::numeric_limits<double>::infinity();
  const Vec3 direction = Unit(relative);
  const Vec3 centre = ParkerDirectionAtRadius(radius, c);
  const double sine = Norm(Cross(direction, centre));
  const double cosine = Clamp(Dot(direction, centre), -1.0, 1.0);
  return radius * std::atan2(sine, cosine);
}

double RequestedCellSizeM(const Vec3& position, const Configuration& c) {
  const double radius = Norm(Subtract(position, c.parkerOriginM));
  double requested = c.globalCellSizeM;
  if (c.solarRefinementEnabled) {
    const double fraction = (radius - c.innerRadiusM) /
        (c.solarTransitionOuterRadiusM - c.innerRadiusM);
    requested = std::min(requested, c.solarSurfaceCellSizeM +
        Profile(fraction, c.solarProfile, c.solarExponent) *
        (c.globalCellSizeM - c.solarSurfaceCellSizeM));
  }
  if (c.tubeRefinementEnabled) {
    const double distance = TubeDistanceM(position, c);
    const double fraction = distance / TubeRadiusM(radius, c);
    requested = std::min(requested, c.tubeCenterCellSizeM +
        Profile(fraction, c.tubeProfile, c.tubeExponent) *
        (c.globalCellSizeM - c.tubeCenterCellSizeM));

    // The production AMR builder samples this point function on a Cartesian
    // lattice. The closest lattice point to a curve crossing a cell is within
    // half a cell diagonal, so max(h_tube,2*d) forces progressive capture of
    // the centreline even when the physical tube is initially sub-cell wide.
    constexpr double kStrictRefinement = 2.0 * (1.0 - 1.0e-12);
    const double capture = std::max(
        c.tubeCenterCellSizeM, kStrictRefinement * distance);
    requested = std::min(requested, capture);
  }
  return Clamp(requested, c.minimumCellSizeM, c.globalCellSizeM);
}

}  // namespace Initialization
}  // namespace SEP
