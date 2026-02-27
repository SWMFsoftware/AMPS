#pragma once
/**
 * @file spectrum.h
 * @brief cSpectrum definition parsed from AMPS_PARAM.in and used for particle injection.
 *
 * The website/wizard produces several spectrum types. This header wraps the spectrum
 * type and parameters into a single class with a unified API:
 *
 *   double cSpectrum::GetSpectrum(double E_J) const;
 *
 * where E_J is kinetic energy (per nucleon, if applicable) in Joules.
 *
 * Notes on units:
 * - Many user-facing spectra are specified as differential flux per (MeV/n):
 *     dF/dE_MeV  [ (m^2 s sr MeV)^-1 ]   (or similar)
 * - GetSpectrum(E_J) returns the differential flux per Joule:
 *     dF/dE_J = dF/dE_MeV * dE_MeV/dE_J = dF/dE_MeV / MeV_IN_J
 *
 * If your injector only needs relative weights (PDF), the absolute scaling is irrelevant,
 * but the Joule-vs-MeV conversion is still important for consistent integration.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>
#include <vector>

class cSpectrum {
public:
  enum class Type : uint8_t {
    PowerLaw,
    PowerLawCutoff,
    LisForceField,
    Band,
    Table,
    Unknown
  };

  // ---------- Construction ----------
  cSpectrum() = default;

  static cSpectrum MakePowerLaw(double J0_perMeV, double gamma, double E0_MeV,
                               double Emin_MeV, double Emax_MeV) {
    cSpectrum s;
    s.type_ = Type::PowerLaw;
    s.spec_j0_ = J0_perMeV;
    s.spec_gamma_ = gamma;
    s.spec_e0_MeV_ = E0_MeV;
    s.spec_emin_MeV_ = Emin_MeV;
    s.spec_emax_MeV_ = Emax_MeV;
    s.ValidateOrThrow();
    return s;
  }

  static cSpectrum MakePowerLawCutoff(double J0_perMeV, double gamma, double E0_MeV,
                                     double Ec_MeV, double Emin_MeV, double Emax_MeV) {
    cSpectrum s;
    s.type_ = Type::PowerLawCutoff;
    s.spec_j0_ = J0_perMeV;
    s.spec_gamma_ = gamma;
    s.spec_e0_MeV_ = E0_MeV;
    s.spec_ec_MeV_ = Ec_MeV;
    s.spec_emin_MeV_ = Emin_MeV;
    s.spec_emax_MeV_ = Emax_MeV;
    s.ValidateOrThrow();
    return s;
  }

  /**
   * @brief Force-field modulated LIS power law.
   *
   * Parameters are specified in the same convention as the wizard:
   * - LIS normalization/spec index and pivot energy define the unmodulated LIS.
   * - phi_MV is the modulation potential in MV (i.e., MeV for charge=1).
   *
   * Implementation uses a common "force-field approximation" (Gleeson & Axford):
   *   J_mod(E) = J_LIS(E + phi) * (E (E + 2 m)) / ((E+phi) (E+phi + 2 m))
   * where E, m, phi are in MeV (per nucleon) and m is rest mass energy.
   *
   * Default rest mass is proton mass energy (938.2720813 MeV).
   */
  static cSpectrum MakeLisForceField(double LisJ0_perMeV, double LisGamma,
                                    double E0_MeV, double Phi_MV,
                                    double Emin_MeV, double Emax_MeV,
                                    double rest_mass_MeV = 938.2720813) {
    cSpectrum s;
    s.type_ = Type::LisForceField;
    s.spec_lis_j0_ = LisJ0_perMeV;
    s.spec_lis_gamma_ = LisGamma;
    s.spec_e0_MeV_ = E0_MeV;
    s.spec_phi_MV_ = Phi_MV;
    s.rest_mass_MeV_ = rest_mass_MeV;
    s.spec_emin_MeV_ = Emin_MeV;
    s.spec_emax_MeV_ = Emax_MeV;
    s.ValidateOrThrow();
    return s;
  }

  /**
   * @brief Band-like spectrum (piecewise with smooth transition).
   *
   * We follow the standard "Band function" form used in high-energy astrophysics,
   * adapted to your parameter names:
   * - gamma1 : low-energy spectral index (positive -> falling with energy)
   * - gamma2 : high-energy spectral index (positive -> falling with energy)
   * - E0_MeV : characteristic energy (MeV) controlling the rollover
   *
   * In canonical Band notation, indices are often (alpha, beta) where N(E) ~ E^alpha
   * at low E; here we assume user enters positive slopes so we use -gamma.
   *
   * Break energy: E_break = (gamma2 - gamma1) * E0.
   */
  static cSpectrum MakeBand(double J0_perMeV, double gamma1, double gamma2,
                           double E0_MeV, double Emin_MeV, double Emax_MeV) {
    cSpectrum s;
    s.type_ = Type::Band;
    s.spec_j0_ = J0_perMeV;
    s.spec_gamma1_ = gamma1;
    s.spec_gamma2_ = gamma2;
    s.spec_e0_MeV_ = E0_MeV;
    s.spec_emin_MeV_ = Emin_MeV;
    s.spec_emax_MeV_ = Emax_MeV;
    s.ValidateOrThrow();
    return s;
  }

  /**
   * @brief Table spectrum: energy(MeV) and flux(per MeV) columns.
   * The file may contain comments beginning with #.
   * Interpolation: log-log linear interpolation (positive values required).
   */
  static cSpectrum MakeTable(std::string table_file,
                            double Emin_MeV, double Emax_MeV) {
    cSpectrum s;
    s.type_ = Type::Table;
    s.table_file_ = std::move(table_file);
    s.spec_emin_MeV_ = Emin_MeV;
    s.spec_emax_MeV_ = Emax_MeV;
    s.LoadTableOrThrow();
    s.ValidateOrThrow();
    return s;
  }

  // ---------- Parsing helpers ----------
  /**
   * @brief Build a cSpectrum from key/value pairs (strings), e.g. parsed from AMPS_PARAM.in.
   *
   * Expected keys (case-insensitive):
   *  - SPECTRUM_TYPE
   *  - SPEC_EMIN, SPEC_EMAX
   *  - POWER_LAW: SPEC_J0, SPEC_GAMMA, SPEC_E0
   *  - POWER_LAW_CUTOFF: + SPEC_EC
   *  - LIS_FORCE_FIELD: SPEC_LIS_J0, SPEC_LIS_GAMMA, SPEC_E0, SPEC_PHI
   *  - BAND: SPEC_J0, SPEC_GAMMA1, SPEC_GAMMA2, SPEC_E0
   *  - TABLE: SPEC_TABLE_FILE
   */
  static cSpectrum FromKeyValueMap(const std::unordered_map<std::string, std::string>& kv) {
    const auto get = [&](const std::string& k) -> std::string {
      auto it = kv.find(k);
      if (it != kv.end()) return it->second;
      // try case-insensitive lookup
      std::string kl = ToUpper(k);
      for (auto it2 = kv.begin(); it2 != kv.end(); ++it2) {
        if (ToUpper(it2->first) == kl) return it2->second;
      }
      return {};
    };

    const std::string type_s = ToUpper(Trim(get("SPECTRUM_TYPE")));
    const double Emin = ParseDoubleOrThrow(get("SPEC_EMIN"), "SPEC_EMIN");
    const double Emax = ParseDoubleOrThrow(get("SPEC_EMAX"), "SPEC_EMAX");

    if (type_s == "POWER_LAW") {
      return MakePowerLaw(
          ParseDoubleOrThrow(get("SPEC_J0"), "SPEC_J0"),
          ParseDoubleOrThrow(get("SPEC_GAMMA"), "SPEC_GAMMA"),
          ParseDoubleOrThrow(get("SPEC_E0"), "SPEC_E0"),
          Emin, Emax);
    }

    if (type_s == "POWER_LAW_CUTOFF") {
      return MakePowerLawCutoff(
          ParseDoubleOrThrow(get("SPEC_J0"), "SPEC_J0"),
          ParseDoubleOrThrow(get("SPEC_GAMMA"), "SPEC_GAMMA"),
          ParseDoubleOrThrow(get("SPEC_E0"), "SPEC_E0"),
          ParseDoubleOrThrow(get("SPEC_EC"), "SPEC_EC"),
          Emin, Emax);
    }

    if (type_s == "LIS_FORCE_FIELD") {
      return MakeLisForceField(
          ParseDoubleOrThrow(get("SPEC_LIS_J0"), "SPEC_LIS_J0"),
          ParseDoubleOrThrow(get("SPEC_LIS_GAMMA"), "SPEC_LIS_GAMMA"),
          ParseDoubleOrThrow(get("SPEC_E0"), "SPEC_E0"),
          ParseDoubleOrThrow(get("SPEC_PHI"), "SPEC_PHI"),
          Emin, Emax);
    }

    if (type_s == "BAND") {
      return MakeBand(
          ParseDoubleOrThrow(get("SPEC_J0"), "SPEC_J0"),
          ParseDoubleOrThrow(get("SPEC_GAMMA1"), "SPEC_GAMMA1"),
          ParseDoubleOrThrow(get("SPEC_GAMMA2"), "SPEC_GAMMA2"),
          ParseDoubleOrThrow(get("SPEC_E0"), "SPEC_E0"),
          Emin, Emax);
    }

    if (type_s == "TABLE") {
      std::string f = Trim(get("SPEC_TABLE_FILE"));
      if (f.empty()) throw std::runtime_error("Missing required key SPEC_TABLE_FILE for SPECTRUM_TYPE=TABLE");
      return MakeTable(f, Emin, Emax);
    }

    std::ostringstream oss;
    oss << "Unrecognized SPECTRUM_TYPE='" << type_s
        << "'. Supported types: POWER_LAW, POWER_LAW_CUTOFF, LIS_FORCE_FIELD, BAND, TABLE";
    throw std::runtime_error(oss.str());
  }

  /**
   * @brief Build a cSpectrum from a parsed key/value map (std::map variant).
   *
   * Many existing AMPS utilities store parsed sections as std::map for deterministic ordering.
   * This overload avoids forcing the caller (e.g., the parser) to convert containers.
   */
  static cSpectrum FromKeyValueMap(const std::map<std::string, std::string>& kv) {
    // Convert to an unordered_map and reuse the primary implementation.
    std::unordered_map<std::string, std::string> u;
    u.reserve(kv.size());
    for (const auto& it : kv) u.emplace(it.first, it.second);
    return FromKeyValueMap(u);
  }

  // ---------- Accessors ----------
  Type GetType() const noexcept { return type_; }
  const std::string& GetTableFile() const noexcept { return table_file_; }

  /**
   * @brief Direct access to table points (only meaningful for Type::Table).
   *
   * Exposed for diagnostics and Tecplot writers that must preserve the exact
   * user-provided table values.
   */
  const std::vector<double>& TableEnergy_MeV() const noexcept { return table_E_MeV_; }
  const std::vector<double>& TableFlux_perMeV() const noexcept { return table_J_perMeV_; }

  double Emin_MeV() const noexcept { return spec_emin_MeV_; }
  double Emax_MeV() const noexcept { return spec_emax_MeV_; }

  // ---------- Main API ----------
  /**
   * @brief Differential spectrum at kinetic energy E_J (Joules).
   * @return dF/dE_J  (same physical flux as user inputs, but per Joule).
   */
  double GetSpectrum(double E_J) const {
    if (!(E_J > 0.0)) return 0.0;

    const double E_MeV = E_J / MeV_IN_J;
    if (E_MeV < spec_emin_MeV_ || E_MeV > spec_emax_MeV_) return 0.0;

    const double val_perMeV = GetSpectrumPerMeV_(E_MeV);
    if (!(val_perMeV > 0.0)) return 0.0;

    return val_perMeV / MeV_IN_J;
  }

  /**
   * @brief Convenience: spectrum per MeV at E_MeV.
   */
  double GetSpectrumPerMeV(double E_MeV) const {
    if (!(E_MeV > 0.0)) return 0.0;
    if (E_MeV < spec_emin_MeV_ || E_MeV > spec_emax_MeV_) return 0.0;
    return GetSpectrumPerMeV_(E_MeV);
  }

  // ---------- Utilities ----------
  static constexpr double EV_IN_J = 1.602176634e-19;
  static constexpr double MeV_IN_J = 1.0e6 * EV_IN_J;

  static std::string TypeToString(Type t) {
    switch (t) {
      case Type::PowerLaw: return "POWER_LAW";
      case Type::PowerLawCutoff: return "POWER_LAW_CUTOFF";
      case Type::LisForceField: return "LIS_FORCE_FIELD";
      case Type::Band: return "BAND";
      case Type::Table: return "TABLE";
      default: return "UNKNOWN";
    }
  }

private:
  // ---------- Internal evaluation ----------
  double GetSpectrumPerMeV_(double E_MeV) const {
    switch (type_) {
      case Type::PowerLaw: {
        // J(E) = J0 * (E/E0)^(-gamma)
        const double x = E_MeV / spec_e0_MeV_;
        return spec_j0_ * std::pow(x, -spec_gamma_);
      }
      case Type::PowerLawCutoff: {
        // J(E) = J0 * (E/E0)^(-gamma) * exp(-E/Ec)
        const double x = E_MeV / spec_e0_MeV_;
        return spec_j0_ * std::pow(x, -spec_gamma_) * std::exp(-E_MeV / spec_ec_MeV_);
      }
      case Type::LisForceField: {
        // Force-field modulation. Treat Phi_MV as MeV for Z=1. For other charge states, caller can
        // scale Phi before construction.
        const double phi_MeV = spec_phi_MV_;
        const double E_LIS = E_MeV + phi_MeV;

        const double x = E_LIS / spec_e0_MeV_;
        const double J_LIS = spec_lis_j0_ * std::pow(x, -spec_lis_gamma_);

        const double m = rest_mass_MeV_;
        const double num = E_MeV * (E_MeV + 2.0 * m);
        const double den = E_LIS * (E_LIS + 2.0 * m);
        if (!(den > 0.0)) return 0.0;

        return J_LIS * (num / den);
      }
      case Type::Band: {
        // Standard Band-like form with indices -gamma1, -gamma2 and rollover E0.
        // Define alpha=-gamma1, beta=-gamma2.
        const double alpha = -spec_gamma1_;
        const double beta  = -spec_gamma2_;
        const double E0 = spec_e0_MeV_;

        const double Eb = (alpha - beta) * E0; // = (gamma2-gamma1)*E0
        if (!(Eb > 0.0) || !(E0 > 0.0)) return 0.0;

        if (E_MeV <= Eb) {
          // A * (E/E0)^alpha * exp(-E/E0)
          const double x = E_MeV / E0;
          return spec_j0_ * std::pow(x, alpha) * std::exp(-x);
        } else {
          // A * (Eb/E0)^(alpha-beta) * exp(beta-alpha) * (E/E0)^beta
          const double pref = std::pow(Eb / E0, (alpha - beta)) * std::exp(beta - alpha);
          const double x = E_MeV / E0;
          return spec_j0_ * pref * std::pow(x, beta);
        }
      }
      case Type::Table: {
        return InterpLogLog_(E_MeV);
      }
      default:
        return 0.0;
    }
  }

  // ---------- Validation / loading ----------
  void ValidateOrThrow() const {
    if (type_ == Type::Unknown) {
      throw std::runtime_error("cSpectrum has Unknown type");
    }
    if (!(spec_emin_MeV_ > 0.0) || !(spec_emax_MeV_ > spec_emin_MeV_)) {
      throw std::runtime_error("cSpectrum invalid energy bounds: require 0 < SPEC_EMIN < SPEC_EMAX (MeV)");
    }

    switch (type_) {
      case Type::PowerLaw:
        RequirePositive_("SPEC_J0", spec_j0_);
        RequirePositive_("SPEC_E0", spec_e0_MeV_);
        RequireFinite_("SPEC_GAMMA", spec_gamma_);
        break;
      case Type::PowerLawCutoff:
        RequirePositive_("SPEC_J0", spec_j0_);
        RequirePositive_("SPEC_E0", spec_e0_MeV_);
        RequirePositive_("SPEC_EC", spec_ec_MeV_);
        RequireFinite_("SPEC_GAMMA", spec_gamma_);
        break;
      case Type::LisForceField:
        RequirePositive_("SPEC_LIS_J0", spec_lis_j0_);
        RequirePositive_("SPEC_E0", spec_e0_MeV_);
        RequirePositive_("REST_MASS_MEV", rest_mass_MeV_);
        RequireFinite_("SPEC_LIS_GAMMA", spec_lis_gamma_);
        RequireFinite_("SPEC_PHI", spec_phi_MV_);
        if (spec_phi_MV_ < 0.0) throw std::runtime_error("SPEC_PHI must be >= 0");
        break;
      case Type::Band:
        RequirePositive_("SPEC_J0", spec_j0_);
        RequirePositive_("SPEC_E0", spec_e0_MeV_);
        RequireFinite_("SPEC_GAMMA1", spec_gamma1_);
        RequireFinite_("SPEC_GAMMA2", spec_gamma2_);
        if ((spec_gamma2_ - spec_gamma1_) <= 0.0) {
          throw std::runtime_error("BAND requires SPEC_GAMMA2 > SPEC_GAMMA1 (so break energy is positive)");
        }
        break;
      case Type::Table:
        if (table_E_MeV_.size() < 2) throw std::runtime_error("TABLE spectrum requires >= 2 data points");
        break;
      default:
        break;
    }
  }

  void LoadTableOrThrow() {
    table_E_MeV_.clear();
    table_J_perMeV_.clear();

    std::ifstream in(table_file_);
    if (!in) {
      throw std::runtime_error("Failed to open SPEC_TABLE_FILE='" + table_file_ + "'");
    }

    std::string line;
    size_t lineno = 0;
    while (std::getline(in, line)) {
      ++lineno;
      line = Trim(line);
      if (line.empty() || line[0] == '#') continue;

      std::istringstream iss(line);
      double e = 0.0, j = 0.0;
      if (!(iss >> e >> j)) {
        std::ostringstream oss;
        oss << "Bad TABLE line " << lineno << " in '" << table_file_ << "': expected two columns (E_MeV flux_perMeV)";
        throw std::runtime_error(oss.str());
      }
      if (!(e > 0.0) || !(j > 0.0)) continue; // skip non-positive entries
      table_E_MeV_.push_back(e);
      table_J_perMeV_.push_back(j);
    }

    if (table_E_MeV_.size() < 2) {
      throw std::runtime_error("TABLE file '" + table_file_ + "' did not contain >=2 valid (positive) rows");
    }

    // Ensure monotonic increasing E; if not, sort pairs.
    std::vector<size_t> idx(table_E_MeV_.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return table_E_MeV_[a] < table_E_MeV_[b]; });

    std::vector<double> E2, J2;
    E2.reserve(idx.size()); J2.reserve(idx.size());
    for (size_t k : idx) { E2.push_back(table_E_MeV_[k]); J2.push_back(table_J_perMeV_[k]); }
    table_E_MeV_.swap(E2);
    table_J_perMeV_.swap(J2);

    // Remove duplicates in energy (keep last).
    std::vector<double> E3, J3;
    for (size_t i = 0; i < table_E_MeV_.size(); ++i) {
      if (!E3.empty() && std::fabs(table_E_MeV_[i] - E3.back()) <= 0.0) {
        J3.back() = table_J_perMeV_[i];
      } else {
        E3.push_back(table_E_MeV_[i]);
        J3.push_back(table_J_perMeV_[i]);
      }
    }
    table_E_MeV_.swap(E3);
    table_J_perMeV_.swap(J3);
  }

  double InterpLogLog_(double E_MeV) const {
    if (table_E_MeV_.empty()) return 0.0;
    if (E_MeV <= table_E_MeV_.front()) return table_J_perMeV_.front();
    if (E_MeV >= table_E_MeV_.back())  return table_J_perMeV_.back();

    auto it = std::upper_bound(table_E_MeV_.begin(), table_E_MeV_.end(), E_MeV);
    const size_t i1 = size_t(it - table_E_MeV_.begin());
    const size_t i0 = i1 - 1;

    const double x0 = table_E_MeV_[i0];
    const double x1 = table_E_MeV_[i1];
    const double y0 = table_J_perMeV_[i0];
    const double y1 = table_J_perMeV_[i1];

    if (!(x1 > x0) || !(y0 > 0.0) || !(y1 > 0.0)) return 0.0;

    const double lx0 = std::log(x0), lx1 = std::log(x1), lx = std::log(E_MeV);
    const double ly0 = std::log(y0), ly1 = std::log(y1);

    const double t = (lx - lx0) / (lx1 - lx0);
    const double ly = (1.0 - t) * ly0 + t * ly1;
    return std::exp(ly);
  }

  // ---------- String/parse helpers ----------
  static std::string Trim(std::string s) {
    auto issp = [](unsigned char c){ return std::isspace(c) != 0; };
    while (!s.empty() && issp((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && issp((unsigned char)s.back())) s.pop_back();
    return s;
  }

  static std::string ToUpper(std::string s) {
    for (char& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
  }

  static double ParseDoubleOrThrow(const std::string& v, const char* key) {
    const std::string s = Trim(v);
    if (s.empty()) {
      std::ostringstream oss;
      oss << "Missing required key " << key;
      throw std::runtime_error(oss.str());
    }
    char* end = nullptr;
    const double x = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || !std::isfinite(x)) {
      std::ostringstream oss;
      oss << "Invalid numeric value for " << key << ": '" << s << "'";
      throw std::runtime_error(oss.str());
    }
    return x;
  }

  static void RequirePositive_(const char* key, double v) {
    if (!(v > 0.0) || !std::isfinite(v)) {
      std::ostringstream oss;
      oss << key << " must be > 0";
      throw std::runtime_error(oss.str());
    }
  }

  static void RequireFinite_(const char* key, double v) {
    if (!std::isfinite(v)) {
      std::ostringstream oss;
      oss << key << " must be finite";
      throw std::runtime_error(oss.str());
    }
  }

private:
  // ---------- Stored definition ----------
  Type type_ = Type::Unknown;

  // Bounds (MeV/n)
  double spec_emin_MeV_ = 0.0;
  double spec_emax_MeV_ = 0.0;

  // Common / power law style
  double spec_j0_ = 0.0;        // per MeV
  double spec_gamma_ = 0.0;
  double spec_e0_MeV_ = 0.0;

  // Cutoff
  double spec_ec_MeV_ = 0.0;

  // LIS force-field
  double spec_lis_j0_ = 0.0;    // per MeV

double spec_lis_gamma_ = 0.0;
  double spec_phi_MV_ = 0.0;    // MV ~ MeV for Z=1
  double rest_mass_MeV_ = 938.2720813;

  // Band
  double spec_gamma1_ = 0.0;
  double spec_gamma2_ = 0.0;

  // Table
  std::string table_file_;
  std::vector<double> table_E_MeV_;
  std::vector<double> table_J_perMeV_;
};

//------------------------------------------------------------------------------
// Global spectrum instance
//------------------------------------------------------------------------------
/**
 * @brief Global spectrum instance initialized by the AMPS_PARAM.in parser.
 *
 * Why global?
 *  - Injection routines are often called from multiple places and kept lightweight; passing
 *    a spectrum object through many call chains can be intrusive in a large legacy codebase.
 *  - Keeping a single, validated instance ensures consistent spectrum usage in injection,
 *    diagnostics, and Tecplot outputs.
 *
 * IMPORTANT:
 *  - This object is defined in spectrum.cpp.
 *  - It must be initialized exactly once after parsing #SPECTRUM, before any injection uses it.
 */
extern cSpectrum gSpectrum;

//------------------------------------------------------------------------------
// Global initializer
//------------------------------------------------------------------------------
/**
 * @brief Initialize the global spectrum from a parsed key/value map (unordered_map variant).
 * Throws std::runtime_error on unknown type or missing/invalid parameters.
 */
void InitGlobalSpectrumFromKeyValueMap(
    const std::unordered_map<std::string, std::string>& kv);

/**
 * @brief Initialize the global spectrum from a parsed key/value map (std::map variant).
 * Convenience overload for parsers that store sections as std::map.
 */
void InitGlobalSpectrumFromKeyValueMap(
    const std::map<std::string, std::string>& kv);

//------------------------------------------------------------------------------
// Tecplot writer for spectrum input
//------------------------------------------------------------------------------
/**
 * @brief Write a Tecplot ASCII file with the spectrum that was parsed from the input.
 *
 * Requirements implemented:
 *  1) Energy is written in MeV.
 *  2) If the spectrum is TABLE, the original table points are written.
 *  3) Otherwise (analytic), energies are log-spaced ("exponential splitting") between
 *     SPEC_EMIN and SPEC_EMAX using @p n_analytic points.
 *
 * The output is meant as a faithful record of what was parsed and what will be used
 * by injection/boundary-condition code.
 */
void WriteSpectrumInputTecplot(const std::string& filename,
                               const cSpectrum& s,
                               int n_analytic = 200);



