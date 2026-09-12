#pragma once

// ============================================================================
// swcme_defaults.hpp
// ----------------------------------------------------------------------------
// Canonical SWCME baseline parameters and model-scope conventions.
//
// Why this header exists
// ----------------------
// Historically the public 1-D and 3-D Params structures carried independent
// default values.  That made a default-constructed 1-D model and a
// default-constructed 3-D model physically different before any geometry was
// involved (for example, different 1-AU density and ejecta thickness).  It
// also made the intended science scope implicit: the same class could be used
// as an upstream Parker/shock-source model or as a phenomenological ICME model
// without an explicit statement of which assumptions were scientifically
// valid.
//
// This header is now the single source of truth for values shared by the two
// dimensional interfaces.  Public Params fields remain in place for source
// compatibility, but their default initializers reference these constants.
// Event-specific values should always be assigned explicitly by the caller and
// recorded by the run/campaign metadata layer.
//
// Scope convention
// ----------------
// The recommended SEP baseline is CONTROLLED_SEP_PRE_SHOCK:
//   regions      = SHOCK_ONLY
//   acceleration = SOURCE
// In this scope the transport-facing plasma remains the analytical
// Parker/Leblanc background and the shock is used for geometry/connectivity/
// source bookkeeping.  The Parker background is declared scientifically valid
// at an observer only while the modeled shock has not reached that observer.
//
// FULL_ICME_DIAGNOSTIC is the optional phenomenological region model:
//   regions      = FULL_ICME
//   acceleration = RESOLVED_COMPRESSION
// It is useful for controlled diagnostics of sheath/ejecta and resolved
// compression, but it must not be confused with a validated global ICME/MHD
// background.  In particular, the ejecta magnetic field remains Parker-like.
// ============================================================================

#include "swcme_constants.hpp"
#include "swcme_solarwind.hpp"
#include "swcme_kinematics.hpp"
#include "swcme_regions.hpp"
#include "swcme_acceleration.hpp"

#include <cmath>

namespace swcme {
namespace defaults {

// Configuration/schema version associated with the conventions documented in
// README.md and DEFAULTS_SCOPE_FIX_NOTES.md.  The later campaign/AMPS adapter
// can write this exact value into run manifests.
constexpr int CONFIG_VERSION = 2;
constexpr const char* FRAME_NAME = "HCI_like_inertial";

// ---- Canonical ambient solar-wind baseline ---------------------------------
constexpr double V_SW_KMS = 400.0;
constexpr double N1AU_CM3 = 5.0;
constexpr double B1AU_TOTAL_NT = 5.0;
constexpr double T_K = 1.2e5;
constexpr double GAMMA_AD = 5.0 / 3.0;

// The public B1AU value is a POSITIVE total field magnitude at 1 AU at the
// equatorial reference latitude (sin(colatitude)=1).  The prepared common core
// converts this to the radial Br normalization.  Positive Br is the adopted
// outward magnetic-polarity convention.  The 3-D field always uses local
// latitude for Parker winding; this reference quantity is normalization only.
constexpr double PARKER_REFERENCE_SIN_THETA = 1.0;
constexpr int PARKER_RADIAL_POLARITY = +1;
constexpr const char* PARKER_NORMALIZATION_CONVENTION =
    "TOTAL_B_AT_1AU_REFERENCE_LATITUDE";
constexpr double SOLAR_ROTATION_RATE_RAD_S =
    swcme::constants::SOLAR_ROTATION_RAD_S;

// ---- Canonical CME/shock-apex kinematics -----------------------------------
constexpr swcme::kinematics::Mode KINEMATICS_MODE =
    swcme::kinematics::Mode::DBM;
constexpr double DBM_R0_RS = 20.0;
constexpr double V0_SH_KMS = 1500.0;
constexpr double DBM_GAMMA_KM_INV = 1.0e-7;
constexpr swcme::kinematics::ExtrapolationPolicy DATA_EXTRAPOLATION =
    swcme::kinematics::ExtrapolationPolicy::OutsideTime;

// ---- Recommended controlled SEP scope -------------------------------------
constexpr swcme::regions::Mode REGION_MODE = swcme::regions::Mode::ShockOnly;
constexpr swcme::acceleration::Mode ACCELERATION_MODE =
    swcme::acceleration::Mode::Source;
constexpr double RELATIVE_SOURCE_WEIGHT_PER_AREA = 1.0;

// ---- Optional FULL_ICME region parameters ---------------------------------
// These are still default-initialized even though they are inactive in the
// controlled SHOCK_ONLY baseline.  Keeping one common value set ensures that
// changing both modes to FULL_ICME/RESOLVED_COMPRESSION yields identical
// region parameters in 1-D and 3-D unless the caller explicitly overrides one.
constexpr double SHEATH_THICK_AU_AT_1AU = 0.10;
constexpr double EJECTA_THICK_AU_AT_1AU = 0.20;
constexpr double EDGE_SMOOTH_SHOCK_AU_AT_1AU = 0.01;
constexpr double EDGE_SMOOTH_LE_AU_AT_1AU = 0.02;
constexpr double EDGE_SMOOTH_TE_AU_AT_1AU = 0.03;
constexpr double SHEATH_RAMP_POWER = 2.0;
constexpr double V_SHEATH_LE_FACTOR = 1.10;
constexpr double F_ME = 0.50;
constexpr double V_ME_FACTOR = 0.80;

// Deprecated compatibility input.  Physical compression is set only by the
// ideal-MHD RH solver.  The canonical ignored value is therefore neutral 1.0
// instead of two different historical 1-D/3-D pseudo-compression floors.
constexpr double SHEATH_COMP_FLOOR_COMPAT = 1.0;

// ---- 3-D science-geometry defaults -----------------------------------------
constexpr double SSE_HALF_WIDTH_RAD =
    40.0 * swcme::constants::PI / 180.0;
constexpr double ELLIPSOID_AXIS_RATIO_Y = 1.0;
constexpr double ELLIPSOID_AXIS_RATIO_Z = 1.0;

// Stable textual names are used by resolved-configuration manifests.  They
// intentionally live beside the defaults so downstream campaign code does not
// invent a second spelling for the same option.
inline const char* kinematics_mode_name(swcme::kinematics::Mode mode) {
  switch (mode) {
    case swcme::kinematics::Mode::Ballistic: return "BALLISTIC";
    case swcme::kinematics::Mode::DBM: return "DBM";
    case swcme::kinematics::Mode::DataDriven: return "DATA_DRIVEN";
  }
  return "UNKNOWN";
}

inline const char* extrapolation_policy_name(
    swcme::kinematics::ExtrapolationPolicy policy) {
  switch (policy) {
    case swcme::kinematics::ExtrapolationPolicy::OutsideTime:
      return "OUTSIDE_TIME";
    case swcme::kinematics::ExtrapolationPolicy::Ballistic:
      return "BALLISTIC";
  }
  return "UNKNOWN";
}

inline const char* region_mode_name(swcme::regions::Mode mode) {
  switch (mode) {
    case swcme::regions::Mode::ShockOnly: return "SHOCK_ONLY";
    case swcme::regions::Mode::FullICME: return "FULL_ICME";
  }
  return "UNKNOWN";
}

inline const char* acceleration_mode_name(swcme::acceleration::Mode mode) {
  switch (mode) {
    case swcme::acceleration::Mode::Source: return "SOURCE";
    case swcme::acceleration::Mode::ResolvedCompression:
      return "RESOLVED_COMPRESSION";
  }
  return "UNKNOWN";
}

// ModelScope is deliberately derived from the two already-validated mode
// switches rather than stored as a third independent switch that could become
// contradictory.  INVALID is returned only for an unvalidated combination.
enum class ModelScope {
  ControlledSEPPreShock,
  FullICMEDiagnostic,
  Invalid
};

constexpr ModelScope model_scope(swcme::regions::Mode region_mode,
                                 swcme::acceleration::Mode acceleration_mode) {
  return (region_mode == swcme::regions::Mode::ShockOnly &&
          acceleration_mode == swcme::acceleration::Mode::Source)
             ? ModelScope::ControlledSEPPreShock
         : (region_mode == swcme::regions::Mode::FullICME &&
            acceleration_mode == swcme::acceleration::Mode::ResolvedCompression)
             ? ModelScope::FullICMEDiagnostic
             : ModelScope::Invalid;
}

inline const char* model_scope_name(ModelScope scope) {
  switch (scope) {
    case ModelScope::ControlledSEPPreShock: return "CONTROLLED_SEP_PRE_SHOCK";
    case ModelScope::FullICMEDiagnostic: return "FULL_ICME_DIAGNOSTIC";
    case ModelScope::Invalid: return "INVALID";
  }
  return "INVALID";
}

// Observer-local declaration of whether the selected model is being used
// within its documented scope.  A finite-SSE front may not exist on an
// observer ray at all; in that case a SHOCK_ONLY Parker background remains
// valid because the modeled disturbance has not reached that observer.
struct ObserverScopeStatus {
  ModelScope scope = ModelScope::Invalid;
  bool valid_input = false;
  bool shock_surface_on_observer_ray = false;
  bool shock_has_reached_observer = false;
  bool within_declared_scope = false;
};

inline ObserverScopeStatus observer_scope_status(
    swcme::regions::Mode region_mode,
    swcme::acceleration::Mode acceleration_mode,
    bool shock_surface_on_observer_ray,
    double shock_radius_m,
    double observer_radius_m) {
  ObserverScopeStatus out;
  out.scope = model_scope(region_mode, acceleration_mode);
  out.shock_surface_on_observer_ray = shock_surface_on_observer_ray;

  // Scope metadata must never turn a malformed radius into a science-valid
  // result.  Numerical/model evaluators have richer status codes; this compact
  // helper simply marks malformed metadata input as outside declared scope.
  // CFG04 defines observer validity against the same analytical-domain
  // boundary used by field evaluation.  A positive but sub-domain observer
  // must not be advertised as scientifically in scope and then fail only when
  // a later Parker/Leblanc query is attempted.  Equality is supported.
  if (!std::isfinite(observer_radius_m) ||
      observer_radius_m < swcme::solarwind::MIN_RADIUS_M ||
      (shock_surface_on_observer_ray &&
       (!std::isfinite(shock_radius_m) ||
        shock_radius_m < swcme::solarwind::MIN_RADIUS_M))) {
    return out;
  }
  out.valid_input = true;

  if (shock_surface_on_observer_ray) {
    // Equality is counted as arrival.  The controlled Parker-background scope
    // is explicitly PRE-shock, so a shock exactly at the observer is already
    // outside that assumption.
    out.shock_has_reached_observer = shock_radius_m >= observer_radius_m;
  }

  switch (out.scope) {
    case ModelScope::ControlledSEPPreShock:
      out.within_declared_scope = !out.shock_has_reached_observer;
      break;
    case ModelScope::FullICMEDiagnostic:
      // The FULL_ICME model remains mathematically evaluable before/inside/
      // behind the modeled front, but this flag means only "within the declared
      // phenomenological diagnostic scope", not observational validation of a
      // realistic ICME magnetic structure.
      out.within_declared_scope = true;
      break;
    case ModelScope::Invalid:
      out.within_declared_scope = false;
      break;
  }
  return out;
}

}  // namespace defaults
}  // namespace swcme
