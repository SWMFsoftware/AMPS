#ifndef SWCME_CONFIG_HPP
#define SWCME_CONFIG_HPP

// ============================================================================
// swcme_config.hpp
// ----------------------------------------------------------------------------
// Shared SWCME configuration-validation infrastructure.
//
// This layer intentionally contains no 1-D or 3-D geometry implementation.
// Instead, dimensional wrappers populate a CommonConfigView and then append
// only the checks that are genuinely specific to their geometry.  This keeps
// the physical-range rules for ambient plasma, kinematics, region thicknesses,
// and smoothing identical in both interfaces while still allowing the 3-D
// model to validate vector axes and finite-SSE parameters.
//
// Validation philosophy
// ---------------------
// * Invalid configuration is rejected once, before any physics calculation.
// * Unit conversion never repairs a value; validation decides admissibility.
// * Every failure identifies the offending field and the violated rule.
// * No NaN/Inf or silent normalization fallback is considered a valid input.
// ============================================================================

#include "swcme_kinematics.hpp"
#include "swcme_regions.hpp"
#include "swcme_acceleration.hpp"
#include "swcme_defaults.hpp"
#include "swcme_solarwind.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace swcme {
namespace config {

enum class Code {
  NonFinite,
  NonPositive,
  Negative,
  OutOfRange,
  ZeroVector,
  InvalidKinematicsTable,
  IncompatibleOptions
};

inline const char* code_name(Code code) {
  switch (code) {
    case Code::NonFinite: return "NON_FINITE";
    case Code::NonPositive: return "NON_POSITIVE";
    case Code::Negative: return "NEGATIVE";
    case Code::OutOfRange: return "OUT_OF_RANGE";
    case Code::ZeroVector: return "ZERO_VECTOR";
    case Code::InvalidKinematicsTable: return "INVALID_KINEMATICS_TABLE";
    case Code::IncompatibleOptions: return "INCOMPATIBLE_OPTIONS";
  }
  return "UNKNOWN";
}

struct Issue {
  std::string field;
  Code code = Code::NonFinite;
  double value = 0.0;
  std::string requirement;
};

struct ValidationResult {
  std::vector<Issue> issues;

  bool ok() const { return issues.empty(); }

  void add(const std::string& field, Code code, double value,
           const std::string& requirement) {
    issues.push_back(Issue{field, code, value, requirement});
  }

  // Human-readable summary used by prepare_step() exceptions and by CFG01.
  // The complete list is preserved so callers can diagnose more than the
  // first invalid field in one validation pass.
  std::string summary(const std::string& prefix = "SWCME invalid configuration") const {
    if (issues.empty()) return prefix + ": none";
    std::ostringstream out;
    out << prefix << ':';
    for (const Issue& issue : issues) {
      out << " " << issue.field << " [" << code_name(issue.code) << "] "
          << issue.requirement << ";";
    }
    return out.str();
  }
};

inline bool finite(double value) { return std::isfinite(value); }

inline void require_finite(ValidationResult& result, const char* field,
                           double value) {
  if (!finite(value)) {
    result.add(field, Code::NonFinite, value, "must be finite");
  }
}

inline void require_positive(ValidationResult& result, const char* field,
                             double value) {
  if (!finite(value)) {
    result.add(field, Code::NonFinite, value, "must be finite and > 0");
  } else if (!(value > 0.0)) {
    result.add(field, Code::NonPositive, value, "must be > 0");
  }
}

inline void require_nonnegative(ValidationResult& result, const char* field,
                                double value) {
  if (!finite(value)) {
    result.add(field, Code::NonFinite, value, "must be finite and >= 0");
  } else if (value < 0.0) {
    result.add(field, Code::Negative, value, "must be >= 0");
  }
}

inline void require_range(ValidationResult& result, const char* field,
                          double value, double minimum, double maximum,
                          bool minimum_inclusive = true,
                          bool maximum_inclusive = true) {
  if (!finite(value)) {
    result.add(field, Code::NonFinite, value, "must be finite and within range");
    return;
  }
  const bool low_ok = minimum_inclusive ? value >= minimum : value > minimum;
  const bool high_ok = maximum_inclusive ? value <= maximum : value < maximum;
  if (!low_ok || !high_ok) {
    std::ostringstream rule;
    rule << "must be " << (minimum_inclusive ? "[" : "(") << minimum
         << ',' << maximum << (maximum_inclusive ? "]" : ")");
    result.add(field, Code::OutOfRange, value, rule.str());
  }
}

inline void require_nonzero_vector(ValidationResult& result, const char* field,
                                   const double vector[3]) {
  if (!finite(vector[0]) || !finite(vector[1]) || !finite(vector[2])) {
    result.add(field, Code::NonFinite, 0.0,
               "all vector components must be finite");
    return;
  }
  const double norm2 = vector[0]*vector[0] + vector[1]*vector[1] +
                       vector[2]*vector[2];
  if (!(norm2 > 0.0) || !finite(norm2)) {
    result.add(field, Code::ZeroVector, 0.0,
               "vector magnitude must be finite and non-zero");
  }
}

// Common physical/model values expressed in the units of the public parameter
// interface.  This is intentionally a non-owning view: dimensional models can
// validate their existing Params without introducing a second configuration
// object or changing user-facing source compatibility.
struct CommonConfigView {
  double V_sw_kms = 0.0;
  double n1AU_cm3 = 0.0;
  double B1AU_nT = 0.0;
  double T_K = 0.0;
  double gamma_ad = 0.0;
  double sin_theta = 0.0;

  swcme::kinematics::Mode kinematics_mode = swcme::defaults::KINEMATICS_MODE;
  double r0_Rs = 0.0;
  double V0_sh_kms = 0.0;
  double Gamma_kmInv = 0.0;
  const std::vector<double>* data_time_s = nullptr;
  const std::vector<double>* data_radius_Rs = nullptr;

  swcme::regions::Mode region_mode = swcme::defaults::REGION_MODE;
  swcme::acceleration::Mode acceleration_mode =
      swcme::defaults::ACCELERATION_MODE;
  double relative_source_weight_per_area =
      swcme::defaults::RELATIVE_SOURCE_WEIGHT_PER_AREA;
  double sheath_thick_AU_at1AU = 0.0;
  double ejecta_thick_AU_at1AU = 0.0;
  double edge_smooth_shock_AU_at1AU = 0.0;
  double edge_smooth_le_AU_at1AU = 0.0;
  double edge_smooth_te_AU_at1AU = 0.0;
  double sheath_comp_floor = 1.0;
  double sheath_ramp_power = 0.0;
  double V_sheath_LE_factor = 0.0;
  double f_ME = 0.0;
  double V_ME_factor = 0.0;
};

inline ValidationResult validate_common(const CommonConfigView& c) {
  ValidationResult out;

  // A strictly positive ambient speed is required by the Parker spiral and by
  // the DBM reference frame.  Importantly, this check happens after pure unit
  // conversion; zero is not silently replaced by an arbitrary 1 m/s value.
  require_positive(out, "V_sw_kms", c.V_sw_kms);
  require_positive(out, "n1AU_cm3", c.n1AU_cm3);
  require_nonnegative(out, "B1AU_nT", c.B1AU_nT);
  require_positive(out, "T_K", c.T_K);
  if (!finite(c.gamma_ad)) {
    out.add("gamma_ad", Code::NonFinite, c.gamma_ad, "must be finite and > 1");
  } else if (!(c.gamma_ad > 1.0)) {
    out.add("gamma_ad", Code::OutOfRange, c.gamma_ad, "must be > 1");
  }
  require_range(out, "sin_theta", c.sin_theta, 0.0, 1.0, true, true);

  // r0 is not merely a positive mathematical radius: it is the physical
  // handoff/reference point at which the analytical Parker/Leblanc model must
  // already be valid.  Accept the exact 1.05-Rsun boundary and reject the next
  // representable value below it before any kinematic state is prepared.
  require_range(out,"r0_Rs",c.r0_Rs,
                swcme::solarwind::MIN_RADIUS_RS,
                std::numeric_limits<double>::infinity(),true,true);
  require_nonnegative(out, "V0_sh_kms", c.V0_sh_kms);
  require_nonnegative(out, "Gamma_kmInv", c.Gamma_kmInv);

  require_nonnegative(out, "sheath_thick_AU_at1AU", c.sheath_thick_AU_at1AU);
  require_nonnegative(out, "ejecta_thick_AU_at1AU", c.ejecta_thick_AU_at1AU);
  // Region thicknesses are specified as AU at a 1-AU shock and therefore are
  // also the self-similar radial fractions used by swcme_regions.hpp.  Their
  // sum must remain below unity so the trailing edge stays at positive radius
  // and no runtime sorting/clipping is needed to repair inverted layers.
  if (finite(c.sheath_thick_AU_at1AU) && finite(c.ejecta_thick_AU_at1AU) &&
      c.sheath_thick_AU_at1AU + c.ejecta_thick_AU_at1AU >= 1.0) {
    out.add("sheath_thick_AU_at1AU+ejecta_thick_AU_at1AU",
            Code::OutOfRange,
            c.sheath_thick_AU_at1AU + c.ejecta_thick_AU_at1AU,
            "self-similar sheath+ejecta fractions must sum to < 1");
  }
  require_nonnegative(out, "edge_smooth_shock_AU_at1AU",
                      c.edge_smooth_shock_AU_at1AU);
  require_nonnegative(out, "edge_smooth_le_AU_at1AU", c.edge_smooth_le_AU_at1AU);
  require_nonnegative(out, "edge_smooth_te_AU_at1AU", c.edge_smooth_te_AU_at1AU);

  // A successful configuration must never be modified later by the region
  // builder.  The former implementation silently capped these inputs at 90%
  // of an adjacent layer, so the archived Params could disagree with the
  // effective physics.  Validate the exact same self-similar limits here,
  // accepting equality and rejecting even the next representable value above
  // it with the individual public field name.
  if (finite(c.sheath_thick_AU_at1AU) &&
      finite(c.ejecta_thick_AU_at1AU) &&
      c.sheath_thick_AU_at1AU>=0.0 &&
      c.ejecta_thick_AU_at1AU>=0.0) {
    const swcme::regions::SmoothingFractionLimits limits=
        swcme::regions::smoothing_fraction_limits(
            c.sheath_thick_AU_at1AU,c.ejecta_thick_AU_at1AU);
    if (finite(c.edge_smooth_shock_AU_at1AU) &&
        c.edge_smooth_shock_AU_at1AU>=0.0 &&
        c.edge_smooth_shock_AU_at1AU>limits.shock) {
      out.add("edge_smooth_shock_AU_at1AU",Code::OutOfRange,
              c.edge_smooth_shock_AU_at1AU,
              "must be <= 0.90 * sheath_thick_AU_at1AU");
    }
    if (finite(c.edge_smooth_le_AU_at1AU) &&
        c.edge_smooth_le_AU_at1AU>=0.0 &&
        c.edge_smooth_le_AU_at1AU>limits.leading) {
      out.add("edge_smooth_le_AU_at1AU",Code::OutOfRange,
              c.edge_smooth_le_AU_at1AU,
              "must be <= 0.90 * min(sheath_thick_AU_at1AU, "
              "ejecta_thick_AU_at1AU)");
    }
    if (finite(c.edge_smooth_te_AU_at1AU) &&
        c.edge_smooth_te_AU_at1AU>=0.0 &&
        c.edge_smooth_te_AU_at1AU>limits.trailing) {
      out.add("edge_smooth_te_AU_at1AU",Code::OutOfRange,
              c.edge_smooth_te_AU_at1AU,
              "must be <= 0.90 * ejecta_thick_AU_at1AU");
    }
  }
  require_nonnegative(out, "relative_source_weight_per_area",
                      c.relative_source_weight_per_area);

  // Shock acceleration has one and only one representation.  The currently
  // supported combinations are intentionally strict because they make double
  // counting impossible by construction:
  //   SOURCE               <-> SHOCK_ONLY
  //   RESOLVED_COMPRESSION <-> FULL_ICME with a finite shock width.
  // A future region model may relax these pairings only if it can prove that a
  // SOURCE transport field contains no second RH compression accelerator.
  if (c.acceleration_mode==swcme::acceleration::Mode::Source &&
      c.region_mode!=swcme::regions::Mode::ShockOnly) {
    out.add("acceleration_mode/region_mode", Code::IncompatibleOptions, 0.0,
            "SOURCE requires SHOCK_ONLY to prevent resolved-shock double counting");
  }
  if (c.acceleration_mode==swcme::acceleration::Mode::ResolvedCompression &&
      c.region_mode!=swcme::regions::Mode::FullICME) {
    out.add("acceleration_mode/region_mode", Code::IncompatibleOptions, 0.0,
            "RESOLVED_COMPRESSION requires FULL_ICME transport flow");
  }
  if (c.acceleration_mode==swcme::acceleration::Mode::ResolvedCompression &&
      finite(c.edge_smooth_shock_AU_at1AU) &&
      !(c.edge_smooth_shock_AU_at1AU > 0.0)) {
    out.add("edge_smooth_shock_AU_at1AU", Code::OutOfRange,
            c.edge_smooth_shock_AU_at1AU,
            "must be > 0 in RESOLVED_COMPRESSION mode");
  }

  if (!finite(c.sheath_comp_floor)) {
    out.add("sheath_comp_floor", Code::NonFinite, c.sheath_comp_floor,
            "must be finite and >= 1");
  } else if (c.sheath_comp_floor < 1.0) {
    out.add("sheath_comp_floor", Code::OutOfRange, c.sheath_comp_floor,
            "must be >= 1");
  }

  if (!finite(c.sheath_ramp_power)) {
    out.add("sheath_ramp_power", Code::NonFinite, c.sheath_ramp_power,
            "must be finite and >= 1");
  } else if (c.sheath_ramp_power < 1.0) {
    out.add("sheath_ramp_power", Code::OutOfRange, c.sheath_ramp_power,
            "must be >= 1");
  }
  // The repaired region model uses the same physical convention in 1-D and
  // 3-D: a forward-shock sheath relaxes toward, but not below, the ambient
  // radial speed before entering the ejecta.  Values below one are therefore
  // rejected instead of being silently clamped by one dimensional wrapper.
  if (!finite(c.V_sheath_LE_factor)) {
    out.add("V_sheath_LE_factor", Code::NonFinite, c.V_sheath_LE_factor,
            "must be finite and >= 1");
  } else if (c.V_sheath_LE_factor < 1.0) {
    out.add("V_sheath_LE_factor", Code::OutOfRange, c.V_sheath_LE_factor,
            "must be >= 1 for the forward-shock sheath model");
  }
  require_nonnegative(out, "f_ME", c.f_ME);
  require_nonnegative(out, "V_ME_factor", c.V_ME_factor);

  // Validate the DATA_DRIVEN table at the public-unit boundary as well as in
  // the SI kinematics layer.  This lets CFG01 identify the offending public
  // field directly instead of reporting only a generic kinematics failure.
  if (c.kinematics_mode == swcme::kinematics::Mode::DataDriven) {
    if (c.data_time_s == nullptr || c.data_radius_Rs == nullptr ||
        c.data_time_s->size() < 2 ||
        c.data_time_s->size() != c.data_radius_Rs->size()) {
      out.add("data_time_s/data_radius_Rs", Code::InvalidKinematicsTable, 0.0,
              "DATA_DRIVEN mode requires matching arrays with at least two knots");
    } else {
      for (std::size_t i=0; i<c.data_time_s->size(); ++i) {
        const double t = (*c.data_time_s)[i];
        const double r = (*c.data_radius_Rs)[i];
        const std::string time_field=
            "data_time_s["+std::to_string(i)+"]";
        const std::string radius_field=
            "data_radius_Rs["+std::to_string(i)+"]";

        // Report malformed members independently so a mixed-validity table
        // identifies every bad knot in one setup pass.  Indexed field names
        // are part of CFG04's diagnostic contract and prevent users from
        // having to locate a sub-domain value by trial and error.
        if (!finite(t)) {
          out.add(time_field,Code::NonFinite,t,
                  "DATA_DRIVEN knot time must be finite");
        }
        if (!finite(r)) {
          out.add(radius_field,Code::NonFinite,r,
                  "DATA_DRIVEN knot radius must be finite");
        } else if (r<swcme::solarwind::MIN_RADIUS_RS) {
          out.add(radius_field,Code::OutOfRange,r,
                  "DATA_DRIVEN knot radius must be >= 1.05 R_sun");
        }

        // Ordering comparisons are meaningful only for finite neighboring
        // values.  Skipping a derivative diagnostic beside a non-finite knot
        // avoids a misleading second error while preserving all primary knot
        // diagnostics above.
        if (i > 0 && finite(t) && finite((*c.data_time_s)[i-1]) &&
            !(t > (*c.data_time_s)[i-1])) {
          out.add("data_time_s", Code::InvalidKinematicsTable, t,
                  "DATA_DRIVEN times must be strictly increasing");
        }
        if (i > 0 && finite(r) && finite((*c.data_radius_Rs)[i-1]) &&
            r < (*c.data_radius_Rs)[i-1]) {
          out.add("data_radius_Rs", Code::InvalidKinematicsTable, r,
                  "DATA_DRIVEN radius must be nondecreasing");
        }
      }
    }
  }

  return out;
}

inline void append(ValidationResult& destination, const ValidationResult& source) {
  destination.issues.insert(destination.issues.end(), source.issues.begin(),
                            source.issues.end());
}

}  // namespace config
}  // namespace swcme

#endif
