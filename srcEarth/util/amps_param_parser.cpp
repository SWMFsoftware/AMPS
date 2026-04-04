//======================================================================================
// amps_param_parser.cpp
//======================================================================================
//
// See amps_param_parser.h for the complete input file format specification,
// type conversion rules, and error-handling contract.
//
//======================================================================================
// PARSING ALGORITHM
//======================================================================================
//
// The parser is a single-pass line-by-line state machine. State is a "current section"
// string, updated whenever a '#'-prefixed section header is encountered:
//
//   state = ""                 (top of file, no section active)
//   for each line L:
//     strip leading/trailing whitespace from L
//     strip comment text after '!' from L
//     if L starts with '#': state = ToUpper(L without '#')
//     else if L is "POINTS_BEGIN": start accumulating points
//     else if L is "POINTS_END":   stop accumulating points
//     else if L is blank: continue
//     else: tokenise L into KEY and VALUE; dispatch to section handler
//
// Dispatch (section handlers are static lambdas capturing the AmpsParam& result):
//   "RUN_INFO"             -> ParseRunInfo
//   "CALCULATION_MODE"     -> ParseCalcMode
//   "CUTOFF_RIGIDITY"      -> ParseCutoffRigidity
//   "DENSITY_SPECTRUM"     -> ParseDensitySpectrum
//   "BOUNDARY_ANISOTROPY"  -> ParseBoundaryAnisotropy
//   "PARTICLE_SPECIES"     -> ParseSpecies
//   "BACKGROUND_FIELD"     -> ParseBackgroundField
//   "DOMAIN_BOUNDARY"      -> ParseDomainBoundary
//   "OUTPUT_DOMAIN"        -> ParseOutputDomain
//   "NUMERICAL"            -> ParseNumerical
//   "SPECTRUM"             -> store raw keys in result.spectrum and build
//                              result.particleSpectrum during post-parse init
//   everything else        -> store in result.unknown raw map
//
//======================================================================================
// TOKENISATION
//======================================================================================
//
// Each non-blank, non-comment line within a section is split on whitespace into
// exactly two tokens: KEY and VALUE. If the line has more than two tokens, everything
// after the first whitespace gap is taken as VALUE (allowing values with spaces, e.g.
// epoch strings like "2003-11-20 06:00"). The rule is:
//
//   pos = line.find_first_of(" \t")
//   key   = line.substr(0, pos)
//   value = trim(line.substr(pos+1))
//
// This handles the most common format (key and value separated by one or more spaces
// or a tab) robustly without requiring a specific delimiter character.
//
//======================================================================================
// POINTS_BEGIN / POINTS_END BLOCK
//======================================================================================
//
// Point coordinates are parsed as three space-separated doubles per line:
//   x  y  z
// Only lines with exactly three tokens in floating-point format are accepted.
// Lines with fewer tokens raise std::runtime_error. Extra trailing tokens are ignored.
//
//======================================================================================
// POST-PARSE VALIDATION
//======================================================================================
//
// After the parse loop completes, a validation pass checks:
//
//   (1) If DS_BOUNDARY_MODE = ANISOTROPIC and the #BOUNDARY_ANISOTROPY section
//       was not present: emit a std::runtime_error with a helpful message listing
//       the required keys.
//
//   (2) If CUTOFF_SAMPLING is neither VERTICAL nor ISOTROPIC: warning on stderr,
//       silently fall back to ISOTROPIC.
//
//   (3) If DS_ENERGY_SPACING is neither LOG nor LINEAR: warning on stderr,
//       silently fall back to LOG.
//
//   (4) If FIELD_MODEL is not in {T96, T05, DIPOLE}: std::runtime_error.
//
// Physical range checks (e.g., Emin < Emax, rInner > 0) are delegated to the
// solver startup rather than the parser, so that the parser remains a pure
// structural parser that can be unit-tested independently of physics.
//
//======================================================================================
// DESIGN CHOICES
//======================================================================================
//
//   (A) All string keys are ToUpper'd before comparison so that both
//       "Dst" and "DST" are accepted.
//
//   (B) Comment stripping happens before tokenisation, so:
//         DST  -50.0  ! nT storm main phase
//       correctly extracts value = "-50.0".
//
//   (C) Boolean values (T/F/TRUE/FALSE/1/0) are parsed by the ToBool helper,
//       which is also exported in the header for use by the CLI.
//
//   (D) The #SPECTRUM section is stored in two forms:
//         * result.spectrum         — raw key/value strings for the existing
//                                     boundary/spectrum.cpp initializer
//         * result.particleSpectrum — typed metadata parsed here for early
//                                     validation and clearer initialization
//       This keeps the parser decoupled from the injection implementation while
//       still making the spectrum choice explicit in AmpsParam.
//
//======================================================================================

#include "amps_param_parser.h"
#include "../boundary/spectrum.h"  // cSpectrum + global spectrum init
#include "specfunc.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <vector>
#include <map>
#include <cstdlib>
#include <cerrno>
#include <cmath>
#include <limits>
#include <cstdio>

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

namespace EarthUtil {

//======================================================================================
// String / unit-conversion utilities
// These must be defined before any function that calls them.
//======================================================================================

static inline std::string Trim(const std::string& s) {
  size_t a=0;
  while (a<s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  size_t b=s.size();
  while (b>a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
  return s.substr(a,b-a);
}

std::string ToUpper(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c){return std::toupper(c);});
  return s;
}

static inline std::string ToLower(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c){return std::tolower(c);});
  return s;
}

bool ToBool(const std::string& sIn) {
  std::string s=ToUpper(Trim(sIn));
  if (s=="T"||s=="TRUE"||s=="1"||s=="YES"||s=="Y") return true;
  if (s=="F"||s=="FALSE"||s=="0"||s=="NO"||s=="N") return false;
  { std::ostringstream _exit_msg; _exit_msg << "Cannot parse boolean token: '"+sIn+"'"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  return false; //The code should not come to this point -- it is here just to make the compiler happy
}

//======================================================================================
// Density + spectrum (gridless) parsing helpers
//======================================================================================

static EarthUtil::DensitySpectrumParam::Spacing ParseEnergySpacingToken(const std::string& s) {
  const std::string t = ToUpper(Trim(s));
  if (t=="LOG") return EarthUtil::DensitySpectrumParam::Spacing::LOG;
  if (t=="LINEAR") return EarthUtil::DensitySpectrumParam::Spacing::LINEAR;
  { std::ostringstream _exit_msg; _exit_msg << "DS_ENERGY_SPACING must be LOG or LINEAR (got '" + Trim(s) + "')"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  return EarthUtil::DensitySpectrumParam::Spacing::LINEAR; //The code should not come to this point -- it is here just to make the compiler happy
}

static inline void SplitKV(const std::string& line,std::string& key,std::string& value) {
  std::istringstream iss(line);
  iss >> key;
  std::getline(iss,value);
  value=Trim(value);
}

static inline std::string StripComment(const std::string& line) {
  size_t p=line.find('!');
  if (p==std::string::npos) return line;
  return line.substr(0,p);
}

// Split an input line into the "code" part and the optional comment after '!'.
// We intentionally preserve the comment text because some AMPS_PARAM files carry
// UNITS in the comment, for example:
//   DOMAIN_X_MAX   20     ! Re
// In that case the numeric token is unitless in the code part, and the parser
// should use the unit hint from the comment.
static inline void SplitCodeAndComment(const std::string& line,std::string& code,std::string& comment) {
  size_t p=line.find('!');
  if (p==std::string::npos) {
    code=line;
    comment.clear();
  }
  else {
    code=line.substr(0,p);
    comment=line.substr(p+1);
  }
}

//======================================================================================
// Spectrum parsing
//======================================================================================
// The AMPS_PARAM.in format contains a #SPECTRUM section that defines the energy spectrum
// of the injected flux (typically differential flux in energy).
//
// We parse the spectrum into a strongly-typed object so downstream output code can
// reproduce the spectrum definition exactly (e.g., in Tecplot headers) and we can
// fail fast if the spectrum type is unknown.
//
// Supported spectrum types (as produced by the wizard website):
//   - POWER_LAW
//   - POWER_LAW_CUTOFF
//   - LIS_FORCE_FIELD
//   - BAND
//   - TABLE
//======================================================================================

/**
 * @brief Enumerates supported spectrum functional forms read from #SPECTRUM.
 *
 * These correspond 1:1 to the wizard-generated keywords in AMPS_PARAM.in.
 * Parsing into an enum avoids string comparisons during hot loops and allows us
 * to validate the input once at startup (fail-fast).
 */

// Parsed spectrum object.
// Note: this file was uploaded without amps_param_parser.h, so we keep this struct
// local. In the real codebase, move this to the header and add a field to AmpsParam.
/**
 * @brief Holds a parsed spectrum definition from AMPS_PARAM.in.
 *
 * This is a "data-only" representation produced by the parser.
 * - It is kept separate from any injection logic on purpose: the parser should
 *   not depend on physics modules, and injection should not depend on parsing.
 * - The Tecplot writer can use this object to emit metadata (AUXDATA) so the
 *   output file records exactly which spectrum was used.
 *
 * All energies are stored in MeV (typically MeV/n). Conversion to Joules is done
 * by downstream physics/injection code as needed.
 */

static inline double GetDoubleOrThrow(const std::map<std::string,std::string>& m,
                                      const std::string& k,
                                      const std::string& ctx) {
  auto it=m.find(k);
  if (it==m.end()) { std::ostringstream _exit_msg; _exit_msg << "Missing required spectrum key '"+k+"' in "+ctx; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  return std::stod(it->second);
}

static inline std::string GetStringOrThrow(const std::map<std::string,std::string>& m,
                                           const std::string& k,
                                           const std::string& ctx) {
  auto it=m.find(k);
  if (it==m.end()) { std::ostringstream _exit_msg; _exit_msg << "Missing required spectrum key '"+k+"' in "+ctx; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  return Trim(it->second);
}

static bool TryGetKey(const std::map<std::string,std::string>& kv,
                      const std::vector<std::string>& keys,
                      std::string& valueOut) {
  for (const auto& k : keys) {
    auto it = kv.find(k);
    if (it != kv.end()) {
      valueOut = Trim(it->second);
      return true;
    }
  }
  return false;
}

static bool TryGetDouble(const std::map<std::string,std::string>& kv,
                         const std::vector<std::string>& keys,
                         double& valueOut) {
  std::string s;
  if (!TryGetKey(kv, keys, s)) return false;
  valueOut = std::stod(s);
  return true;
}

static EarthUtil::ParticleSpectrum::Type ParseSpectrumTypeToken(const std::string& s) {
  const std::string t = ToUpper(Trim(s));
  if (t=="POWER_LAW") return EarthUtil::ParticleSpectrum::Type::POWER_LAW;
  if (t=="POWER_LAW_CUTOFF") return EarthUtil::ParticleSpectrum::Type::POWER_LAW_CUTOFF;
  if (t=="LIS_FORCE_FIELD") return EarthUtil::ParticleSpectrum::Type::LIS_FORCE_FIELD;
  if (t=="BAND") return EarthUtil::ParticleSpectrum::Type::BAND;
  if (t=="TABLE") return EarthUtil::ParticleSpectrum::Type::TABLE;
  return EarthUtil::ParticleSpectrum::Type::UNKNOWN;
}

EarthUtil::ParticleSpectrum ParseParticleSpectrum(const std::map<std::string,std::string>& kv) {
  EarthUtil::ParticleSpectrum spec;
  spec.raw = kv;
  if (kv.empty()) return spec;

  const std::string ctx = "#SPECTRUM";
  std::string typeToken;
  if (!TryGetKey(kv, {"SPECTRUM_TYPE", "TYPE"}, typeToken)) {
    exit(__LINE__,__FILE__,"Missing required key SPECTRUM_TYPE in #SPECTRUM");
  }

  spec.typeName = ToUpper(typeToken);
  spec.type = ParseSpectrumTypeToken(spec.typeName);
  spec.parsed = true;

  if (spec.type == EarthUtil::ParticleSpectrum::Type::UNKNOWN) {
    std::ostringstream _exit_msg;
    _exit_msg << "Unsupported SPECTRUM_TYPE='" << typeToken
              << "'. Supported tokens: POWER_LAW | POWER_LAW_CUTOFF | LIS_FORCE_FIELD | BAND | TABLE";
    exit(__LINE__,__FILE__,_exit_msg.str().c_str());
  }

  // Common aliases produced by different wizard / legacy inputs.
  const bool hasJ0    = TryGetDouble(kv, {"SPEC_J0", "J0", "J_REFERENCE", "SPEC_J_REFERENCE"}, spec.J0);
  const bool hasGamma = TryGetDouble(kv, {"SPEC_GAMMA", "GAMMA", "SPECTRAL_INDEX", "SPEC_SPECTRAL_INDEX"}, spec.gamma);
  const bool hasE0    = TryGetDouble(kv, {"SPEC_E0", "E0", "E0_MEV", "E_REFERENCE", "SPEC_E0_MEV", "SPEC_E_REFERENCE"}, spec.E0_MeV);
  TryGetDouble(kv, {"SPEC_EMIN", "EMIN", "ENERGY_MIN", "SPEC_ENERGY_MIN"}, spec.Emin_MeV);
  TryGetDouble(kv, {"SPEC_EMAX", "EMAX", "ENERGY_MAX", "SPEC_ENERGY_MAX"}, spec.Emax_MeV);
  TryGetDouble(kv, {"SPEC_ECUT", "ECUT", "E_CUTOFF", "SPEC_E_CUTOFF", "CUTOFF_ENERGY_MEV"}, spec.cutoffEnergy_MeV);
  TryGetDouble(kv, {"PHI_MV", "MODULATION_POTENTIAL_MV", "FORCE_FIELD_PHI_MV"}, spec.modulationPotential_MV);
  TryGetDouble(kv, {"ALPHA", "BAND_ALPHA"}, spec.alpha);
  TryGetDouble(kv, {"BETA", "BAND_BETA"}, spec.beta);
  TryGetDouble(kv, {"E_BREAK", "BREAK_ENERGY_MEV", "BAND_EBREAK"}, spec.breakEnergy_MeV);
  {
    std::string tableFile;
    if (TryGetKey(kv, {"SPEC_TABLE_FILE", "TABLE_FILE", "SPECTRUM_FILE", "FILE"}, tableFile)) spec.tableFile = tableFile;
  }

  // Strong validation for the currently used/tested spectrum form.
  if (spec.type == EarthUtil::ParticleSpectrum::Type::POWER_LAW) {
    if (!hasJ0 || !(spec.J0 > 0.0)) {
      exit(__LINE__,__FILE__,"POWER_LAW spectrum requires SPEC_J0/J0 > 0");
    }
    if (!hasE0 || !(spec.E0_MeV > 0.0)) {
      exit(__LINE__,__FILE__,"POWER_LAW spectrum requires SPEC_E0/E0 > 0 (MeV)");
    }
    if (!hasGamma) {
      exit(__LINE__,__FILE__,"POWER_LAW spectrum requires SPEC_GAMMA/GAMMA (or SPECTRAL_INDEX)");
    }
    if (spec.Emin_MeV > 0.0 && spec.Emax_MeV > 0.0 && !(spec.Emax_MeV > spec.Emin_MeV)) {
      exit(__LINE__,__FILE__,"POWER_LAW spectrum requires SPEC_EMAX > SPEC_EMIN when both are provided");
    }
  }

  return spec;
}


namespace {

enum class SpectrumTableFileFormat { UNKNOWN, TWO_COLUMN, TIME_DEPENDENT };

static inline bool TryParseScalarToken(const std::string& sIn, double& valueOut) {
  const std::string s = Trim(sIn);
  if (s.empty()) return false;
  char* end = nullptr;
  errno = 0;
  const double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str() || errno == ERANGE) return false;
  while (end && *end && std::isspace(static_cast<unsigned char>(*end))) ++end;
  if (end && *end) return false;
  if (!std::isfinite(v)) return false;
  valueOut = v;
  return true;
}

static inline bool ParseUtcTimestampToUnixSeconds(const std::string& sIn, double& secondsOut) {
  std::string s = Trim(sIn);
  if (s.empty()) return false;
  if (!s.empty() && (s.back()=='Z' || s.back()=='z')) s.pop_back();

  int year=0, month=0, day=0, hour=0, minute=0;
  double second = 0.0;
  int consumed = 0;

  if (std::sscanf(s.c_str(), "%d-%d-%dT%d:%d:%lf%n", &year, &month, &day, &hour, &minute, &second, &consumed) < 6) {
    consumed = 0;
    second = 0.0;
    if (std::sscanf(s.c_str(), "%d-%d-%d %d:%d:%lf%n", &year, &month, &day, &hour, &minute, &second, &consumed) < 6) {
      consumed = 0;
      if (std::sscanf(s.c_str(), "%d-%d-%dT%d:%d%n", &year, &month, &day, &hour, &minute, &consumed) < 5) {
        consumed = 0;
        if (std::sscanf(s.c_str(), "%d-%d-%d %d:%d%n", &year, &month, &day, &hour, &minute, &consumed) < 5) {
          consumed = 0;
          if (std::sscanf(s.c_str(), "%d-%d-%d%n", &year, &month, &day, &consumed) < 3) return false;
          hour = 0; minute = 0; second = 0.0;
        }
        else {
          second = 0.0;
        }
      }
      else {
        second = 0.0;
      }
    }
  }

  while (consumed < static_cast<int>(s.size()) && std::isspace(static_cast<unsigned char>(s[consumed]))) ++consumed;
  if (consumed != static_cast<int>(s.size())) return false;

  if (month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 || hour > 23 || minute < 0 || minute > 59) return false;
  if (!(second >= 0.0) || !(second < 61.0) || !std::isfinite(second)) return false;

  const int y = year - (month <= 2 ? 1 : 0);
  const int era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned mp = static_cast<unsigned>(month + (month > 2 ? -3 : 9));
  const unsigned doy = (153u * mp + 2u) / 5u + static_cast<unsigned>(day - 1);
  const unsigned doe = yoe * 365u + yoe / 4u - yoe / 100u + doy;
  const long long days = static_cast<long long>(era) * 146097LL + static_cast<long long>(doe) - 719468LL;

  secondsOut = static_cast<double>(days) * 86400.0 + static_cast<double>(hour * 3600 + minute * 60) + second;
  return true;
}

static inline std::vector<std::string> SplitWhitespaceTokens(const std::string& line) {
  std::istringstream iss(line);
  std::vector<std::string> out;
  std::string tok;
  while (iss >> tok) out.push_back(tok);
  return out;
}

static inline std::string StripLeadingCommentMarkers(const std::string& line) {
  std::string s = Trim(line);
  while (!s.empty() && (s[0]=='#' || s[0]=='!')) s = Trim(s.substr(1));
  return s;
}

static SpectrumTableFileFormat DetectSpectrumTableFileFormat(const std::string& fileName) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::ostringstream oss;
    oss << "Cannot open spectrum table file: " << fileName;
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  std::string line;
  std::size_t lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string trimmed = Trim(line);
    if (trimmed.empty()) continue;

    if (trimmed[0]=='#' || trimmed[0]=='!') {
      const std::string header = ToUpper(StripLeadingCommentMarkers(trimmed));
      if (header.find("ENERGY_MEV:") == 0 ||
          header.find("COLUMNS: TIMESTAMP") == 0 ||
          header.find("COLUMN 1") == 0 ||
          header.find("N_ENERGY_BINS") == 0) {
        return SpectrumTableFileFormat::TIME_DEPENDENT;
      }
      continue;
    }

    const std::vector<std::string> toks = SplitWhitespaceTokens(trimmed);
    if (toks.empty()) continue;

    double tSeconds = 0.0;
    if (ParseUtcTimestampToUnixSeconds(toks[0], tSeconds)) {
      return SpectrumTableFileFormat::TIME_DEPENDENT;
    }

    if (toks.size() == 2) {
      double e = 0.0, j = 0.0;
      if (TryParseScalarToken(toks[0], e) && TryParseScalarToken(toks[1], j)) {
        return SpectrumTableFileFormat::TWO_COLUMN;
      }
    }

    std::ostringstream oss;
    oss << "Cannot determine spectrum-table format for '" << fileName
        << "' from line " << lineNo << ": '" << trimmed << "'";
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  std::ostringstream oss;
  oss << "Spectrum table file '" << fileName << "' is empty";
  exit(__LINE__,__FILE__,oss.str().c_str());
  return SpectrumTableFileFormat::UNKNOWN;
}

static std::vector<double> ParseTimeDependentSpectrumEnergyGridOrThrow(const std::string& fileName) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::ostringstream oss;
    oss << "Cannot open spectrum table file: " << fileName;
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  std::string line;
  while (std::getline(in, line)) {
    const std::string header = StripLeadingCommentMarkers(line);
    const std::string upper = ToUpper(header);
    if (upper.find("ENERGY_MEV:") != 0) continue;

    const std::string payload = Trim(header.substr(std::string("Energy_MeV:").size()));
    const std::vector<std::string> toks = SplitWhitespaceTokens(payload);
    std::vector<double> energyGrid;
    energyGrid.reserve(toks.size());
    for (const auto& tok : toks) {
      double e = 0.0;
      if (!TryParseScalarToken(tok, e) || !(e > 0.0)) {
        std::ostringstream oss;
        oss << "Invalid energy value '" << tok << "' in time-dependent spectrum header of '" << fileName << "'";
        exit(__LINE__,__FILE__,oss.str().c_str());
      }
      energyGrid.push_back(e);
    }
    if (energyGrid.size() < 2) {
      std::ostringstream oss;
      oss << "Time-dependent spectrum file '" << fileName << "' must define at least two energies in Energy_MeV header";
      exit(__LINE__,__FILE__,oss.str().c_str());
    }
    return energyGrid;
  }

  std::ostringstream oss;
  oss << "Time-dependent spectrum file '" << fileName << "' is missing the Energy_MeV header";
  exit(__LINE__,__FILE__,oss.str().c_str());
  return std::vector<double>();
}

static inline std::string BasenameOnly(const std::string& path) {
  const std::size_t pos = path.find_last_of("/\\");
  return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

static inline std::string SanitizeFileStem(std::string s) {
  for (char& c : s) {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c=='.' || c=='_' || c=='-')) c = '_';
  }
  return s;
}

static std::string SelectSpectrumInitializationEpochUTC(const EarthUtil::AmpsParam& p) {
  const std::string mode = ToUpper(p.output.mode);
  if (mode == "TRAJECTORY") {
    for (const auto& tr : p.output.trajectories) {
      if (!tr.samples.empty()) return tr.samples.front().timeUTC;
    }
  }
  return p.field.epoch;
}

static std::string MaterializeTimeDependentSpectrumSnapshotOrThrow(const std::string& fileName,
                                                                   const std::string& refEpochUTC) {
  double refSeconds = 0.0;
  if (!ParseUtcTimestampToUnixSeconds(refEpochUTC, refSeconds)) {
    std::ostringstream oss;
    oss << "Cannot parse reference epoch for time-dependent spectrum selection: '" << refEpochUTC << "'";
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  const std::vector<double> energyGrid = ParseTimeDependentSpectrumEnergyGridOrThrow(fileName);

  std::ifstream in(fileName.c_str());
  if (!in) {
    std::ostringstream oss;
    oss << "Cannot open spectrum table file: " << fileName;
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  double bestDelta = std::numeric_limits<double>::infinity();
  std::string bestEpochUTC;
  std::vector<double> bestFlux;

  std::string line;
  std::size_t lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string trimmed = Trim(line);
    if (trimmed.empty() || trimmed[0]=='#' || trimmed[0]=='!') continue;

    const std::vector<std::string> toks = SplitWhitespaceTokens(trimmed);
    if (toks.size() < 2) continue;

    double rowSeconds = 0.0;
    if (!ParseUtcTimestampToUnixSeconds(toks[0], rowSeconds)) continue;

    if (toks.size() < 1 + energyGrid.size()) {
      std::ostringstream oss;
      oss << "Time-dependent spectrum row " << lineNo << " in '" << fileName
          << "' has " << (toks.size() - 1) << " flux columns but the Energy_MeV header declares "
          << energyGrid.size();
      exit(__LINE__,__FILE__,oss.str().c_str());
    }

    std::vector<double> flux(energyGrid.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i=0; i<energyGrid.size(); ++i) {
      double v = 0.0;
      if (TryParseScalarToken(toks[1+i], v) && v > 0.0) flux[i] = v;
    }

    const double delta = std::fabs(rowSeconds - refSeconds);
    if (delta < bestDelta) {
      bestDelta = delta;
      bestEpochUTC = toks[0];
      bestFlux.swap(flux);
    }
  }

  if (bestFlux.empty()) {
    std::ostringstream oss;
    oss << "No timestamped spectrum rows were found in time-dependent spectrum file '" << fileName << "'";
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  std::ostringstream name;
  name << "spectrum_snapshot_" << SanitizeFileStem(BasenameOnly(fileName)) << ".dat";
  const std::string outFile = name.str();

  std::ofstream out(outFile.c_str());
  if (!out) {
    std::ostringstream oss;
    oss << "Cannot create snapshot spectrum file: " << outFile;
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  out << "# Extracted from time-dependent spectrum file\n";
  out << "# source_file: " << fileName << "\n";
  out << "# reference_epoch_utc: " << refEpochUTC << "\n";
  out << "# selected_epoch_utc: " << bestEpochUTC << "\n";
  out << "# selected_time_offset_s: " << bestDelta << "\n";
  out << "# columns: Energy_MeV  differential_flux_perMeV\n";

  std::size_t nValid = 0;
  out.setf(std::ios::scientific);
  out.precision(12);
  for (std::size_t i=0; i<energyGrid.size() && i<bestFlux.size(); ++i) {
    if (!(bestFlux[i] > 0.0) || !std::isfinite(bestFlux[i])) continue;
    out << energyGrid[i] << " " << bestFlux[i] << "\n";
    ++nValid;
  }

  if (nValid < 2) {
    std::ostringstream oss;
    oss << "Selected time-dependent spectrum snapshot from '" << fileName
        << "' at epoch '" << bestEpochUTC << "' contains fewer than two valid positive bins";
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  return outFile;
}

static void ResolveTableSpectrumFileForInitialization(EarthUtil::AmpsParam& p) {
  if (p.particleSpectrum.type != EarthUtil::ParticleSpectrum::Type::TABLE) return;

  std::string tableFile = Trim(p.particleSpectrum.tableFile);
  if (tableFile.empty()) {
    std::map<std::string,std::string>::const_iterator it = p.spectrum.find("SPEC_TABLE_FILE");
    if (it != p.spectrum.end()) tableFile = Trim(it->second);
  }
  if (tableFile.empty()) {
    std::map<std::string,std::string>::const_iterator it = p.spectrum.find("TABLE_FILE");
    if (it != p.spectrum.end()) tableFile = Trim(it->second);
  }
  if (tableFile.empty()) return;

  p.spectrum["SPEC_TABLE_FILE"] = tableFile;
  p.particleSpectrum.tableFile = tableFile;

  const SpectrumTableFileFormat fmt = DetectSpectrumTableFileFormat(tableFile);
  if (fmt == SpectrumTableFileFormat::TWO_COLUMN) {
    p.particleSpectrum.tableFormat = EarthUtil::ParticleSpectrum::TableFormat::TWO_COLUMN;
  }
  else if (fmt == SpectrumTableFileFormat::TIME_DEPENDENT) {
    p.particleSpectrum.tableFormat = EarthUtil::ParticleSpectrum::TableFormat::TIME_DEPENDENT;
  }
  else {
    std::ostringstream oss;
    oss << "Unsupported spectrum table format for file '" << tableFile << "'";
    exit(__LINE__,__FILE__,oss.str().c_str());
  }

  const std::string refEpochUTC = SelectSpectrumInitializationEpochUTC(p);
  p.spectrum["SPEC_TABLE_REFERENCE_EPOCH_UTC"] = refEpochUTC;
}


} // anonymous namespace

/**
 * @brief Parse and validate the #SPECTRUM section.
 *
 * @param kv  Key/value table collected by the generic AMPS_PARAM.in reader for the
 *            #SPECTRUM section only (keys already uppercased by the parser).
 *
 * @return SpectrumSpec with type and parameters filled in.
 *
 * Failure policy (important):
 * - Unknown SPECTRUM_TYPE => throw with a clear list of supported types.
 * - Missing required parameters for a known type => throw (fail-fast).
 * - Emin/Emax invalid => throw.
 *
 * Rationale:
 * Spectra define injected particle populations. If the spectrum is misread, the
 * resulting simulation can be scientifically invalid. Therefore we prefer explicit
 * errors over silent defaults.
 */

// Write spectrum definition as Tecplot AUXDATA entries.
// Typical usage: call this while writing the Tecplot header (before ZONE).
/**
 * @brief Emit spectrum parameters into Tecplot header as AUXDATA entries.
 *
 * This is meant to be called by the Tecplot writer (not by the parser itself)
 * so that each output file records the exact spectrum configuration used.
 *
 * Why AUXDATA:
 * - It is preserved with the dataset.
 * - It is easy to parse later for provenance.
 * - It does not interfere with VARIABLES/ZONE data.
 */

// Earth radius used for conversion of "Re" to kilometers in the parser.
// The parser stores all length-like quantities internally in km. Downstream code
// can convert km->Re (for Tsyganenko calls) or km->m (for particle pushers).
#ifndef AMPS_EARTH_RADIUS_KM
static const double kEarthRadiusKm = 6371.2;
#else
static const double kEarthRadiusKm = AMPS_EARTH_RADIUS_KM;
#endif

// Parse a single length token with optional inline unit suffix (e.g. "10Re",
// "2000km"). If no inline suffix is present, optionally use a unit hint taken
// from the trailing comment after '!'. If neither is present, default to km.
//
// IMPORTANT MODIFICATION (comment-embedded units):
//   Some AMPS_PARAM files document units as part of a human-readable comment,
//   e.g.
//      DOMAIN_X_MAX  15.0   ! RE dayside
//      R_INNER        2.0   ! RE inner loss sphere
//   In such cases the unit token is NOT equal to the entire comment string.
//   Therefore this parser now scans the comment text and extracts the first
//   recognized unit token ("km" or "Re", case-insensitive) from anywhere in
//   the comment.
//
// Supported units (case-insensitive): km, Re.
static double ParseLengthToKm(const std::string& token,const std::string& unitHintFromComment="") {
  std::string t=Trim(token);
  if (t.empty()) exit(__LINE__,__FILE__,"Empty length token");

  const char* begin=t.c_str();
  char* end=nullptr;
  errno=0;
  double value=std::strtod(begin,&end);
  if (begin==end) {
    { std::ostringstream _exit_msg; _exit_msg << "Cannot parse numeric length token: '"+t+"'"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }
  if (errno==ERANGE) {
    { std::ostringstream _exit_msg; _exit_msg << "Length token is out of range: '"+t+"'"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }

  std::string suffix=ToLower(Trim(std::string(end)));
  std::string hint=ToLower(Trim(unitHintFromComment));

  // Helper: normalize a potential unit token extracted from comments.
  // We strip common punctuation so comments like "! RE, dayside" or
  // "! (km) inner sphere" are accepted. We intentionally keep only alphabetic
  // characters because supported units are alphabetic (km, Re).
  auto NormalizeUnitToken = [](std::string x) {
    std::string y;
    y.reserve(x.size());
    for (char c : x) {
      unsigned char uc = static_cast<unsigned char>(c);
      if (std::isalpha(uc)) y.push_back(static_cast<char>(std::tolower(uc)));
    }
    return y;
  };

  // If the value token itself has an explicit suffix, it takes precedence over
  // any comment hint.
  if (!suffix.empty()) {
    if (suffix=="km") return value;
    if (suffix=="re") return value*kEarthRadiusKm;
    { std::ostringstream _exit_msg; _exit_msg << "Unsupported inline length unit '"+suffix+"' in token '"+t+"'"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }

  // No inline suffix: use comment hint if present. The comment may be a pure
  // unit token ("Re") OR a longer descriptive phrase containing the unit
  // token ("RE dayside", "km inner loss sphere"). We scan tokens and use the
  // first recognized unit.
  if (!hint.empty()) {
    if (hint=="km") return value;
    if (hint=="re") return value*kEarthRadiusKm;

    // Scan the *entire* comment text token-by-token because the unit may be
    // embedded in a longer phrase, e.g. "RE dayside", "RE nightside",
    // "km inner loss sphere". The first recognized unit wins.
    std::istringstream hss(unitHintFromComment);
    std::string ht;
    while (hss >> ht) {
      std::string hu = NormalizeUnitToken(Trim(ht));
      if (hu=="km") return value;
      if (hu=="re") return value*kEarthRadiusKm;
    }

    // Free-text comment without a recognizable unit token: fall back to the
    // backward-compatible default (km) instead of throwing. This keeps comments
    // like "! rectangular box in GSM" harmless for length fields while still
    // supporting explicit unit hints embedded in prose.
    return value;
  }

  // Backward-compatible default when no units are provided anywhere.
  return value; // km
}

// Parse a sequence of length values from a whitespace-separated string. Units may
// be provided inline (e.g. "1Re 200km") and/or in the comment after '!'.
// Supported comment conventions:
//   ! Re           -> apply one unit to all values that have no inline suffix
//   ! km km Re     -> per-value unit hints in order
// Any value with an inline suffix keeps that suffix and ignores the comment hint.
static std::vector<double> ParseLengthListToKm(const std::string& values,const std::string& comment) {
  std::vector<std::string> valueToks;
  {
    std::istringstream iss(values);
    std::string t;
    while (iss>>t) valueToks.push_back(t);
  }

  // IMPORTANT ROBUSTNESS FIX
  // ------------------------
  // A comment attached to a POINT line often contains *descriptive prose*, e.g.
  //
  //   POINT 50969.6 0 0   ! TC-EQ8: 8 Re equator (GSM x-axis)
  //
  // The earlier implementation scanned the full comment and collected every token that
  // looked like a supported unit.  In the example above it would find the word "Re"
  // in the prose fragment "8 Re equator" and incorrectly interpret that as a unit hint
  // for the coordinate list.  The result was catastrophic: a coordinate already given in
  // km (50969.6) was reinterpreted as 50969.6 Re, i.e. multiplied by Earth's radius.
  //
  // That bug affected any list-style length field parsed with ParseLengthListToKm(), but
  // it was most visible for OUTPUT_MODE=POINTS because the analytic dipole benchmark then
  // saw observation points absurdly far from Earth and returned near-zero/zero reference
  // values.  The numerical solver itself still used the same misparsed coordinates, which
  // is why the output X_km/Y_km/Z_km columns looked obviously wrong (orders of magnitude
  // too large).
  //
  // The corrected rule is intentionally strict:
  //   * a comment-derived unit hint is accepted only if the comment is *purely a unit
  //     specification*, either
  //         ! Re
  //     or  ! km km km
  //     i.e. one global unit token or one token per listed value;
  //   * as soon as the comment contains any non-unit text, it is treated as descriptive
  //     prose and contributes NO unit hint.
  //
  // This keeps explicit unit comments working, while preventing accidental capture of
  // words like "Re" in free-form descriptions.
  std::vector<std::string> rawCommentToks;
  {
    std::istringstream iss(comment);
    std::string t;
    while (iss>>t) rawCommentToks.push_back(t);
  }

  auto NormalizeUnitToken = [](std::string x) {
    std::string y;
    y.reserve(x.size());
    for (char c : x) {
      unsigned char uc = static_cast<unsigned char>(c);
      if (std::isalpha(uc)) y.push_back(static_cast<char>(std::tolower(uc)));
    }
    return y;
  };

  bool commentIsPureUnitList = !rawCommentToks.empty();
  std::vector<std::string> unitHints;
  unitHints.reserve(rawCommentToks.size());
  for (const auto& t : rawCommentToks) {
    const std::string u = NormalizeUnitToken(Trim(t));
    if (u=="km" || u=="re") {
      unitHints.push_back(u);
    }
    else {
      // Any non-unit token means the comment is descriptive prose rather than an
      // unambiguous unit declaration.  Discard *all* comment-derived hints.
      commentIsPureUnitList = false;
      unitHints.clear();
      break;
    }
  }

  std::vector<double> out;
  out.reserve(valueToks.size());

  const bool oneGlobalHint = commentIsPureUnitList && (unitHints.size()==1);
  const bool perValueHints = commentIsPureUnitList && (unitHints.size()==valueToks.size());

  for (size_t i=0;i<valueToks.size();++i) {
    std::string hint;
    if (perValueHints) hint=unitHints[i];
    else if (oneGlobalHint) hint=unitHints[0];
    out.push_back(ParseLengthToKm(valueToks[i],hint));
  }
  return out;
}


static inline bool ParseIsoUtcToEt(const std::string& timeUTC, double& et) {
#ifndef _NO_SPICE_CALLS_
  try {
    SpiceDouble etSpice = 0.0;
    str2et_c(timeUTC.c_str(), &etSpice);
    et = etSpice;
    return true;
  }
  catch (...) {
    return false;
  }
#else
  (void)timeUTC; (void)et;
  return false;
#endif
}

static inline EarthUtil::Vec3 RotateVectorWithSpiceOrThrow(const EarthUtil::Vec3& vin,
                                                           const char* fromFrame,
                                                           const char* toFrame,
                                                           const std::string& timeUTC,
                                                           const std::string& context) {
#ifndef _NO_SPICE_CALLS_
  SpiceDouble et = 0.0;
  if (!ParseIsoUtcToEt(timeUTC, et)) {
    { std::ostringstream _exit_msg; _exit_msg << "Cannot parse trajectory timestamp for SPICE transform in "
                                              << context << ": '" << timeUTC << "'"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }

  SpiceDouble rot[3][3];
  pxform_c(fromFrame, toFrame, et, rot);

  SpiceDouble in[3]  = {vin.x, vin.y, vin.z};
  SpiceDouble out[3] = {0.0,0.0,0.0};
  mxv_c(rot, in, out);
  return {out[0], out[1], out[2]};
#else
  (void)vin; (void)fromFrame; (void)toFrame; (void)timeUTC; (void)context;
  exit(__LINE__,__FILE__,"Trajectory frame conversion requires SPICE, but the code was built with _NO_SPICE_CALLS_");
  return EarthUtil::Vec3{};
#endif
}

static inline EarthUtil::Vec3 GeoLatLonAltToCartesianKm(double lat_deg,double lon_deg,double alt_km) {
  const double r_km = kEarthRadiusKm + alt_km;
  const double pi = std::acos(-1.0);
  const double lat = lat_deg * pi / 180.0;
  const double lon = lon_deg * pi / 180.0;
  const double clat = std::cos(lat);
  return {r_km*clat*std::cos(lon), r_km*clat*std::sin(lon), r_km*std::sin(lat)};
}

static EarthUtil::Vec3 ConvertTrajectorySampleToGsmM(const std::string& frame,
                                                     const std::string& timeUTC,
                                                     double c1,double c2,double c3,
                                                     const std::string& fileName,
                                                     int lineNo) {
  const std::string f = ToUpper(Trim(frame));
  const std::string ctx = std::string("trajectory file '") + fileName + "' line " + std::to_string(lineNo);

  // All branches ultimately produce a GSM Cartesian vector in SI meters.
  if (f=="GSM") {
    // Input: Cartesian GSM in Earth radii.  Convert Re -> m.
    const double Re_m = kEarthRadiusKm * 1000.0;
    return {c1*Re_m, c2*Re_m, c3*Re_m};
  }
  else if (f=="SM") {
    const double Re_m = kEarthRadiusKm * 1000.0;
    const EarthUtil::Vec3 xSM_km{c1*kEarthRadiusKm, c2*kEarthRadiusKm, c3*kEarthRadiusKm};
    const EarthUtil::Vec3 xGSM_km = RotateVectorWithSpiceOrThrow(xSM_km, "SM", "GSM", timeUTC, ctx);
    return {xGSM_km.x*1000.0, xGSM_km.y*1000.0, xGSM_km.z*1000.0};
  }
  else if (f=="GEO") {
    const EarthUtil::Vec3 xGeo_km = GeoLatLonAltToCartesianKm(c1, c2, c3);
    const EarthUtil::Vec3 xGSM_km = RotateVectorWithSpiceOrThrow(xGeo_km, "ITRF93", "GSM", timeUTC, ctx);
    return {xGSM_km.x*1000.0, xGSM_km.y*1000.0, xGSM_km.z*1000.0};
  }

  { std::ostringstream _exit_msg; _exit_msg << "Unsupported TRAJ_FRAME='" << frame
                                            << "' in " << ctx
                                            << ". Supported values: GEO | GSM | SM"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  return EarthUtil::Vec3{};
}

EarthUtil::SpacecraftTrajectory LoadTrajectoryFileAsGsm(const std::string& fileName,
                                                     const std::string& trajFrame) {
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    { std::ostringstream _exit_msg; _exit_msg << "Cannot open trajectory file: " << fileName; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }

  EarthUtil::SpacecraftTrajectory tr;
  tr.name = fileName;
  tr.sourceFrame = ToUpper(Trim(trajFrame));

  std::string line;
  int lineNo = 0;
  while (std::getline(fin, line)) {
    ++lineNo;

    // Strip comments: in the trajectory file both '#' and '!' introduce a comment
    // that extends to the end of the line.  We strip whichever character comes
    // first so that lines like
    //   2017-09-10T00:00:00Z  4.716  -0.684  -0.415  # first point
    //   # full-line comment
    //   ! also a comment
    // are all handled correctly.
    // Note: '#' is NOT treated as a comment in the main AMPS_PARAM parser (it marks
    // section headers there), so this stripping is intentionally local to this
    // function and not applied via the global StripComment helper.
    {
      const std::size_t ph = line.find('#');
      const std::size_t pe = line.find('!');
      const std::size_t pcut = std::min(ph, pe);   // npos if both absent, else the earlier one
      if (pcut != std::string::npos) line = line.substr(0, pcut);
    }

    line = Trim(line);
    if (line.empty()) continue;   // blank or comment-only lines are skipped

    std::istringstream iss(line);
    std::string timeUTC;
    double c1=0.0, c2=0.0, c3=0.0;
    if (!(iss >> timeUTC >> c1 >> c2 >> c3)) {
      continue;  // fewer than 4 tokens on a non-blank line — skip gracefully
    }

    // Require an ISO-8601 'T' separator to distinguish data rows from any
    // residual header lines that slip through (e.g. column-name rows whose
    // first token is a plain string with no 'T').
    if (timeUTC.find('T') == std::string::npos) continue;

    // Position stored in meters (SI) in xGSM_m.
    const EarthUtil::Vec3 xGSM_m = ConvertTrajectorySampleToGsmM(trajFrame, timeUTC, c1, c2, c3, fileName, lineNo);
    tr.AddSample(timeUTC, xGSM_m);
  }

  if (tr.samples.empty()) {
    { std::ostringstream _exit_msg; _exit_msg << "No trajectory samples were loaded from " << fileName; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }

  return tr;
}

static void PackStandalonePointsIntoSyntheticTrajectories(EarthUtil::AmpsParam& p) {
  p.output.trajectories.clear();
  p.output.trajectories.reserve(p.output.points.size());

  for (size_t i=0; i<p.output.points.size(); ++i) {
    EarthUtil::SpacecraftTrajectory tr;
    tr.name = std::string("POINT_") + std::to_string(i);
    tr.sourceFrame = ToUpper(p.output.coords);
    // output.points is in km; trajectory samples are stored in meters.
    const EarthUtil::Vec3& pt_km = p.output.points[i];
    tr.AddSample(p.field.epoch, {pt_km.x*1000.0, pt_km.y*1000.0, pt_km.z*1000.0});
    p.output.trajectories.push_back(tr);
  }
}

void InitializePointLikeOutput(EarthUtil::AmpsParam& p) {
  const std::string mode = ToUpper(Trim(p.output.mode));

  if (mode=="TRAJECTORY") {
    if (Trim(p.output.trajFile).empty()) {
      exit(__LINE__,__FILE__,"OUTPUT_MODE=TRAJECTORY requires TRAJ_FILE in #OUTPUT_DOMAIN");
    }
    p.output.trajectories.clear();
    p.output.trajectories.push_back(LoadTrajectoryFileAsGsm(p.output.trajFile, p.output.trajFrame));
    p.output.RebuildFlattenedPointsFromTrajectories();
    return;
  }

  if (mode=="POINTS") {
    if (p.output.points.empty()) {
      exit(__LINE__,__FILE__,"OUTPUT_MODE=POINTS but no POINT entries were found in POINTS_BEGIN/END block");
    }
    PackStandalonePointsIntoSyntheticTrajectories(p);
    p.output.RebuildFlattenedPointsFromTrajectories();
    return;
  }

  if (mode=="SHELLS") return;

  {
    std::ostringstream _exit_msg;
    _exit_msg << "Unsupported OUTPUT_MODE='" << p.output.mode
              << "'. Supported values: POINTS | TRAJECTORY | SHELLS";
    exit(__LINE__,__FILE__,_exit_msg.str().c_str());
  }
}


//======================================================================================
// TsDriverTable — implementation
//======================================================================================

static inline double Lerp(double a, double b, double t) { return a + t*(b-a); }

EarthUtil::TsDriverRecord EarthUtil::TsDriverTable::Lookup(double et) const {
  if (records_.empty())
    exit(__LINE__,__FILE__,"TsDriverTable::Lookup called on an empty table");

  if (et <= records_.front().et) {
    if (!clampWarnedLow_) {
      std::cerr << "[TsDriverTable] WARNING: time before table start ("
                << records_.front().timeUTC << "); clamping.\n";
      clampWarnedLow_ = true;
    }
    return records_.front();
  }
  if (et >= records_.back().et) {
    if (!clampWarnedHigh_) {
      std::cerr << "[TsDriverTable] WARNING: time after table end ("
                << records_.back().timeUTC << "); clamping.\n";
      clampWarnedHigh_ = true;
    }
    return records_.back();
  }

  // Binary search: records_[lo].et <= et < records_[lo+1].et
  std::size_t lo = 0, hi = records_.size() - 1;
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    if (records_[mid].et <= et) lo = mid; else hi = mid;
  }

  const EarthUtil::TsDriverRecord& r0 = records_[lo];
  const EarthUtil::TsDriverRecord& r1 = records_[hi];
  const double dt = r1.et - r0.et;
  const double t  = (dt > 0.0) ? (et - r0.et) / dt : 0.0;

  EarthUtil::TsDriverRecord out;
  out.et       = et;
  out.timeUTC  = r0.timeUTC;
  out.imfBy_nT = Lerp(r0.imfBy_nT, r1.imfBy_nT, t);
  out.imfBz_nT = Lerp(r0.imfBz_nT, r1.imfBz_nT, t);
  out.swVx_kms = Lerp(r0.swVx_kms, r1.swVx_kms, t);
  out.swN_cm3  = Lerp(r0.swN_cm3,  r1.swN_cm3,  t);
  out.pdyn_nPa = Lerp(r0.pdyn_nPa, r1.pdyn_nPa, t);
  out.dst_nT   = Lerp(r0.dst_nT,   r1.dst_nT,   t);
  for (int i = 0; i < 6; ++i) out.w[i] = Lerp(r0.w[i], r1.w[i], t);
  return out;
}

void EarthUtil::TsDriverTable::ApplyToField(const EarthUtil::TsDriverRecord& rec,
                                             EarthUtil::BackgroundField& field) {
  field.imfBy_nT = rec.imfBy_nT;
  field.imfBz_nT = rec.imfBz_nT;
  field.swVx_kms = rec.swVx_kms;
  field.swN_cm3  = rec.swN_cm3;
  field.pdyn_nPa = rec.pdyn_nPa;
  field.dst_nT   = rec.dst_nT;
  for (int i = 0; i < 6; ++i) field.w[i] = rec.w[i];
}

//======================================================================================
//======================================================================================
// NormalizeTsModelName
//======================================================================================
// PURPOSE
//   Map every known spelling of a Tsyganenko model name to the single canonical
//   string used throughout this codebase.  This centralises all alias handling
//   so that no downstream comparison needs to know about alternate spellings.
//
// WHY MULTIPLE SPELLINGS EXIST
//   The Tsyganenko models have been distributed with slightly different name
//   conventions over the years:
//     • CCMC Runs-on-Request uses "T96", "T05", "TS05" interchangeably depending
//       on the interface version.
//     • The ViRBO / Qin-Denton data portal labels the T05 model "TS05".
//     • Some older AMPS input files were generated with "T04S" (the Fortran
//       function name of the T05 external field is t04_s_).
//     • User-supplied input files may use any of the above, or variants like
//       "T05S", "TS04".
//   Rather than spreading if-chains throughout the code, we normalise once here
//   and everywhere else compares against the canonical form only.
//
// CANONICAL NAMES (matching the strings used in GetB_T / cFieldEvaluator)
//   "T96"   — Tsyganenko (1996)
//   "T05"   — Tsyganenko-Sitnov (2005); Fortran entry t04_s_
//   "T01"   — Tsyganenko (2001)
//   "TA15N" — Tsyganenko-Andreeva (2015), northward IMF variant
//   "TA15B" — Tsyganenko-Andreeva (2015), southward IMF variant
//   "TA16"  — Tsyganenko-Andreeva (2016) (future; Fortran not yet linked)
//   "DIPOLE"— Analytic centred dipole; no external-field call
//
// EXTENSIBILITY
//   To add a new alias, append one line to the table below.
//   To add a new canonical model, add a new entry in GetB_T, RequiredColumnsForModel,
//   and cFieldEvaluator::ReinitGeopack, then optionally add aliases here.
static std::string NormalizeTsModelName(const std::string& raw) {
  const std::string u = ToUpper(Trim(raw));

  // T05 family — all refer to the Tsyganenko-Sitnov (2005) model whose
  // Fortran entry point is t04_s_.  "TS04" / "T04S" come from that name.
  if (u == "TS05" || u == "T05S" || u == "T04S" || u == "TS04") return "T05";

  // T96 family — alternate CCMC/ViRBO spellings.
  if (u == "TS96" || u == "T96S")                                 return "T96";

  // T01 family — Tsyganenko (2001).
  if (u == "TS01" || u == "T01S")                                 return "T01";

  // TA15 — two distinct sub-variants (northward / southward IMF orientation);
  // we keep the sub-variant suffix because they use different PARMOD layouts.
  if (u == "TA15N" || u == "TA15B")                               return u;

  // TA16 — placeholder; no aliases known yet.
  if (u == "TA16")                                                 return "TA16";

  // All other values (T96, T05, T01, DIPOLE, ...) are already canonical.
  return u;
}

//======================================================================================
//======================================================================================
// TsColMap — variable name (uppercase) -> 0-based column index
//======================================================================================
//
// WHY WE NEED THIS
//   Different versions of the Qin-Denton / ViRBO format, and driver files produced
//   by other sources (CCMC, custom scripts), may place the same physical quantity
//   in a different column.  Hardcoding column numbers (e.g. "Dst is always column
//   25") makes the loader fragile and produces silently wrong results whenever the
//   layout differs — which happens, for example, between the 1-minute and 5-minute
//   Qin-Denton products, or between different model driver formats (T96 vs T05 vs
//   TA15).
//
//   We avoid this entirely by parsing the JSON-like metadata block that every
//   ViRBO-format file embeds in its comment header.  That block explicitly maps
//   every variable NAME to its START_COLUMN, so the loader never has to guess.
//
// HEADER FORMAT (excerpt from qd_sep2017.txt)
//   The metadata block is delimited by comment lines:
//     # {
//     ...
//     # } End JSON
//
//   Inside, each variable is described by a JSON-like stanza, e.g.:
//
//     #  "NAME": "Pdyn",
//     #  "TITLE": "Solar Wind Dynamic Pressure",
//     #  "START_COLUMN": 11,
//     #  "UNITS": "nPa"
//
//   Vector (array) variables also carry an ELEMENT_NAMES list:
//
//     #  "NAME": "W",
//     #  "DIMENSION": [ 6 ],
//     #  "START_COLUMN": 32,
//     #  "ELEMENT_NAMES": [ "W1", "W2", "W3", "W4", "W5", "W6" ]
//
//   For vectors we expand the list into individual entries so that the caller can
//   look up any element by name: "W1"->32, "W2"->33, ..., "W6"->37.
//   The base name ("W") also gets an entry pointing at the first element.
//
// HOW PARSEDTSDRIVER HEADER BUILDS THE MAP
//   The function is a simple state machine:
//   • It accumulates (NAME, START_COLUMN, ELEMENT_NAMES) in three local variables.
//   • Whenever it sees a new NAME entry it "commits" the previous variable by
//     writing its entry (or entries, for vectors) into the map.
//   • At the end of the JSON block it commits whatever is still pending.
//   • Lines outside the JSON block and non-# data lines are ignored.
//
// THREAD SAFETY
//   TsColMap is a plain std::map<std::string,int>; it is built once by the parser
//   (single-threaded) and then used read-only by all MPI ranks/threads.
using TsColMap = std::map<std::string, int>;

//======================================================================================
// ParseTsDriverHeader
//======================================================================================
// Reads the JSON metadata block embedded in the comment header of a ViRBO-format
// driver file and returns a TsColMap mapping each variable name (uppercase) to its
// 0-based column index in the data rows.
//
// ARGUMENTS
//   fin      — open input stream positioned at the beginning of the file.
//   fileName — used only in warning messages.
//
// RETURN VALUE
//   A TsColMap with one entry per scalar variable and one entry per vector element
//   (plus one entry for the vector base name).  If no JSON block is found, an empty
//   map is returned and a warning is written to stderr.
//
// STATE MACHINE OVERVIEW
//   State variable: inJsonBlock (false until "# {" is seen)
//   Per-variable accumulator: currentName, currentStartCol, currentElements
//   The accumulator is committed to the map each time a new "NAME" line is seen,
//   ensuring that the final variable in the block is also written when "} End JSON"
//   is encountered.
//
// PARSING RULES
//   • Only '#'-prefixed lines are processed; the first non-# line ends the header.
//   • Three JSON keys are recognised (case-insensitive match on the uppercased line):
//       "NAME"          — begins a new variable entry; commits any pending entry first
//       "START_COLUMN"  — the 0-based column index of the first (or only) element
//       "ELEMENT_NAMES" — optional list of per-element names for array variables
//   • All other JSON keys (TITLE, LABEL, UNITS, DIMENSION, DESCRIPTION, ...) are
//     intentionally ignored; they carry no information needed by the loader.
//
// ROBUSTNESS
//   • The function does NOT require strict JSON: missing commas, extra whitespace,
//     and reordered keys within a variable block are all tolerated.
//   • If a variable block lacks a START_COLUMN it is silently skipped.
//   • If a variable block lacks ELEMENT_NAMES it is treated as a scalar.
static TsColMap ParseTsDriverHeader(std::istream& fin,
                                    const std::string& fileName) {
  TsColMap colMap;

  // ── Per-variable accumulator ──────────────────────────────────────────────
  // These three variables hold the state for the variable block currently being
  // parsed.  They are reset by commitCurrent() and updated incrementally as the
  // corresponding JSON keys are encountered.
  std::string currentName;             // text of the most recent "NAME" value
  int         currentStartCol = -1;    // value of "START_COLUMN" (-1 = not yet seen)
  std::vector<std::string> currentElements; // contents of "ELEMENT_NAMES" (empty = scalar)

  bool inJsonBlock = false; // true after "# {" has been seen
  bool headerDone  = false; // true once we have left the header section

  // ── Lambda: commitCurrent ─────────────────────────────────────────────────
  // Write the accumulated variable state into colMap and reset the accumulator.
  // Called each time a new "NAME" line is seen (before updating currentName) and
  // at the end of the JSON block.
  //
  // Scalar variable: one map entry   NAME -> startCol
  // Vector variable: one entry per   ELEMENT_NAME[k] -> (startCol + k)
  //                  plus a base     NAME -> startCol   (points at element 0)
  //
  // If currentName is empty or currentStartCol is -1 (incomplete block),
  // the function does nothing — this handles the first call before any NAME
  // has been seen, and any block whose START_COLUMN was missing.
  auto commitCurrent = [&]() {
    if (currentName.empty() || currentStartCol < 0) return;
    const std::string uName = ToUpper(currentName);
    if (currentElements.empty()) {
      // Scalar: single column entry.
      colMap[uName] = currentStartCol;
    } else {
      // Vector: one entry per named element at consecutive columns.
      for (int k = 0; k < static_cast<int>(currentElements.size()); ++k) {
        colMap[ToUpper(currentElements[k])] = currentStartCol + k;
      }
      // Base name also maps to the first element (column startCol + 0).
      colMap[uName] = currentStartCol;
    }
    // Reset for the next variable.
    currentName.clear();
    currentStartCol = -1;
    currentElements.clear();
  };

  // ── Lambda: extractQuotedValue ────────────────────────────────────────────
  // Extract the text between the first pair of double-quotes that follows the
  // colon on a line such as:
  //   #  "NAME": "ByIMF",
  // Returns "" if the expected pattern is not found.
  auto extractQuotedValue = [](const std::string& line) -> std::string {
    const std::size_t colon = line.find(':');
    if (colon == std::string::npos) return "";
    const std::string after = line.substr(colon + 1);
    const std::size_t q1 = after.find('"');
    if (q1 == std::string::npos) return "";
    const std::size_t q2 = after.find('"', q1 + 1);
    if (q2 == std::string::npos) return "";
    return after.substr(q1 + 1, q2 - q1 - 1);
  };

  // ── Lambda: extractIntValue ───────────────────────────────────────────────
  // Extract the integer that follows the colon on a line such as:
  //   #  "START_COLUMN": 7,
  // The trailing comma is stripped before parsing.  Returns -1 on failure.
  auto extractIntValue = [](const std::string& line) -> int {
    const std::size_t colon = line.find(':');
    if (colon == std::string::npos) return -1;
    std::string after = line.substr(colon + 1);
    // Replace the first comma with a space so stoi/>> does not stumble on it.
    for (char& c : after) if (c == ',') { c = ' '; break; }
    std::istringstream iss(after);
    int v = -1;
    iss >> v;
    return v;
  };

  // ── Lambda: extractStringArray ────────────────────────────────────────────
  // Extract a JSON array of quoted strings from a line such as:
  //   #  "ELEMENT_NAMES": [ "W1", "W2", "W3", "W4", "W5", "W6" ],
  // Returns an empty vector if the '[' ... ']' delimiters are absent.
  // Scanning is character-by-character within the bracket content so that
  // element names containing digits (e.g. "G1") or underscores are preserved.
  auto extractStringArray = [](const std::string& line) -> std::vector<std::string> {
    std::vector<std::string> result;
    const std::size_t lb = line.find('[');
    const std::size_t rb = line.find(']');
    if (lb == std::string::npos || rb == std::string::npos || rb <= lb) return result;
    const std::string inner = line.substr(lb + 1, rb - lb - 1);
    std::size_t pos = 0;
    while (pos < inner.size()) {
      const std::size_t q1 = inner.find('"', pos);
      if (q1 == std::string::npos) break;
      const std::size_t q2 = inner.find('"', q1 + 1);
      if (q2 == std::string::npos) break;
      result.push_back(inner.substr(q1 + 1, q2 - q1 - 1));
      pos = q2 + 1;
    }
    return result;
  };

  // ── Main scan loop ────────────────────────────────────────────────────────
  // We read line-by-line.  Only '#'-prefixed lines are part of the header;
  // the first non-# line means we have left the header section entirely.
  std::string line;
  while (!headerDone && std::getline(fin, line)) {
    if (line.empty() || line[0] != '#') {
      // This is the first data row (or a blank non-comment line).
      // We cannot "put back" a line with std::istream, so we simply mark
      // the header as done.  The caller opens a fresh stream for the data pass.
      headerDone = true;
      break;
    }

    // Strip the leading '#' and trim whitespace.
    const std::string stripped = Trim(line.substr(1));

    // Detect the JSON block boundaries.
    if (stripped == "{") { inJsonBlock = true; continue; }
    if (!inJsonBlock) continue;  // outside JSON block — skip source-file comments etc.
    if (stripped.find("} End JSON") != std::string::npos ||
        stripped.find("}") == 0) {
      // End of JSON block.  Commit the last pending variable and stop.
      commitCurrent();
      headerDone = true;
      break;
    }

    // Inside the JSON block: inspect this line for the three keys we care about.
    // We uppercase the line for the key search so the match is case-insensitive,
    // but pass the original 'stripped' to the extraction helpers so quoted values
    // are not accidentally uppercased.
    const std::string ustripped = ToUpper(stripped);

    if (ustripped.find("\"NAME\"") != std::string::npos) {
      // A new variable is starting.  Commit whatever was being accumulated for
      // the previous variable before resetting the accumulator.
      commitCurrent();
      currentName = extractQuotedValue(stripped);
    }
    else if (ustripped.find("\"START_COLUMN\"") != std::string::npos) {
      // The 0-based column index of this variable's first (or only) element.
      currentStartCol = extractIntValue(stripped);
    }
    else if (ustripped.find("\"ELEMENT_NAMES\"") != std::string::npos) {
      // Optional: names for each component of a vector variable.
      // If present, commitCurrent will create one map entry per element.
      currentElements = extractStringArray(stripped);
    }
    // All other JSON keys (TITLE, LABEL, UNITS, DIMENSION, DESCRIPTION,
    // VALUES, ...) are ignored.
  }

  // Commit any variable whose closing brace may have been the end-of-file.
  commitCurrent();

  if (colMap.empty() && !fileName.empty()) {
    std::cerr << "[TsDriverTable] WARNING: no column-map entries found in header of '"
              << fileName << "'. Proceeding without header-driven column detection.\n";
  }
  return colMap;
}

//======================================================================================
// RequiredColumnsForModel
//======================================================================================
// PURPOSE
//   Returns the set of TsColMap variable names (uppercase) that must be present
//   in the driver file to drive the named Tsyganenko model.  Called by
//   ValidateTsDriverColumns to produce a clear, model-specific error message when
//   a required column is absent.
//
// ARGUMENT
//   model — canonical model name (output of NormalizeTsModelName, already uppercase).
//
// VARIABLE NAME CONVENTIONS
//   The names in the returned vector must match the keys in the TsColMap produced
//   by ParseTsDriverHeader, which are the uppercase "NAME" values from the JSON
//   metadata block.  Standard Qin-Denton / ViRBO mappings:
//
//     ViRBO NAME   TsColMap key   Physical meaning
//     ----------   ------------   -----------------------------------------
//     ByIMF        BYIMF          IMF Y-component [nT]
//     BzIMF        BZIMF          IMF Z-component [nT]
//     Vsw          VSW            Solar wind speed magnitude [km/s]
//     Den_P        DEN_P          Solar wind proton density [cm^-3]
//     Pdyn         PDYN           Solar wind dynamic pressure [nPa]
//     Dst          DST            Dst index [nT]
//     G1..G3       G1..G3         T01 storm-activity G-indices
//     W1..W6       W1..W6         T05/TA15/TA16 storm-time history integrals
//     Bz1..Bz6     BZ1..BZ6       TA15 |Bz_south| running averages [nT]
//
// HOW TO ADD A NEW MODEL
//   1. Add an `if (m == "MODELNAME")` block.
//   2. Build the required list from swBase and any additional variable names.
//   3. Optionally register aliases in NormalizeTsModelName.
//   4. Add the Fortran call in cFieldEvaluator::GetB_T.
//   5. Update PARMOD loading in cFieldEvaluator constructor and ReinitGeopack.
static std::vector<std::string> RequiredColumnsForModel(const std::string& model) {
  // Shared by all external Tsyganenko models: IMF, solar wind, Dst.
  // These are the minimum inputs needed to call any T96/T01/T05/TA15/TA16
  // Fortran routine.
  static const std::vector<std::string> swBase = {
    "BYIMF",  // IMF By [nT]
    "BZIMF",  // IMF Bz [nT]
    "VSW",    // Solar wind speed [km/s] (positive magnitude in file)
    "DEN_P",  // Solar wind proton density [cm^-3]
    "PDYN",   // Solar wind dynamic pressure [nPa]
    "DST"     // Dst index [nT]
  };

  // W1..W6: storm-time history integrals introduced by T05.
  // They encode the time-integrated energy injection into each internal current
  // system (ring current, partial ring current, field-aligned, Chapman-Ferraro,
  // tail) and are also required by TA15 and TA16.
  static const std::vector<std::string> wIntegrals = {
    "W1", "W2", "W3", "W4", "W5", "W6"
  };

  const std::string m = ToUpper(Trim(model));

  // DIPOLE: analytic internal field only; no external Fortran call, no inputs.
  if (m == "DIPOLE") return {};

  // T96: driven by [Pdyn, Dst, By, Bz] only.
  // PARMOD: [0]=Pdyn [1]=Dst [2]=By [3]=Bz [4..10]=0
  if (m == "T96") return swBase;

  // T01: same as T96 plus the G-indices G1, G2, G3 that capture the recent
  // history of strong southward Bz intervals and their effect on the ring current.
  // PARMOD: [0]=Pdyn [1]=Dst [2]=By [3]=Bz [4]=G1 [5]=G2 [6]=G3 [7..10]=0
  if (m == "T01") {
    std::vector<std::string> req = swBase;
    req.insert(req.end(), {"G1", "G2", "G3"});
    return req;
  }

  // T05: adds the six W storm-time integrals to the base set.
  // PARMOD: [0]=Pdyn [1]=Dst [2]=By [3]=Bz [4]=W1 .. [9]=W6 [10]=0
  if (m == "T05") {
    std::vector<std::string> req = swBase;
    req.insert(req.end(), wIntegrals.begin(), wIntegrals.end());
    return req;
  }

  // TA15 (Tsyganenko-Andreeva 2015): two sub-variants (northward / southward IMF).
  // In addition to the W-integrals it needs BZ1..BZ6 — the mean |Bz_south| over
  // 1, 2, 3, 4, 5, 6 hours before the observation epoch.
  if (m == "TA15N" || m == "TA15B") {
    std::vector<std::string> req = swBase;
    req.insert(req.end(), wIntegrals.begin(), wIntegrals.end());
    req.insert(req.end(), {"BZ1", "BZ2", "BZ3", "BZ4", "BZ5", "BZ6"});
    return req;
  }

  // TA16 (Tsyganenko-Andreeva 2016): Fortran interface not yet linked.
  // We register the expected driver columns now so the parser gives a useful
  // validation error rather than silently using wrong values.
  if (m == "TA16") {
    std::vector<std::string> req = swBase;
    req.insert(req.end(), wIntegrals.begin(), wIntegrals.end());
    return req;
  }

  // Unknown model — skip validation and let the solver report the problem.
  std::cerr << "[TsDriverTable] WARNING: no required-column definition for model '"
            << model << "'. Column validation will be skipped.\n";
  return {};
}
//======================================================================================
// ValidateTsDriverColumns
//======================================================================================
// PURPOSE
//   Verify that every column required by `model` is present in `colMap`.
//   Called once after ParseTsDriverHeader, before any data rows are read, so
//   the user gets a clear error message at load time rather than a silent wrong
//   value during the trajectory run.
//
// ERROR MESSAGE
//   On failure the function lists:
//   • each missing column by name, and
//   • all available column names (from the file header) so the user can see
//     whether the file uses a different spelling or simply lacks the variable.
//
// DESIGN NOTE
//   We exit rather than throw because a missing required column is always a
//   configuration error that cannot be recovered from at runtime.  Continuing
//   with a default of zero for, say, Dst or W1 would produce physically wrong
//   cutoff rigidities with no diagnostic.
static void ValidateTsDriverColumns(const TsColMap& colMap,
                                    const std::string& model,
                                    const std::string& fileName) {
  const std::vector<std::string> required = RequiredColumnsForModel(model);

  // Collect all missing names in one pass so the error message is complete.
  std::vector<std::string> missing;
  for (const auto& name : required) {
    if (colMap.find(name) == colMap.end()) missing.push_back(name);
  }
  if (missing.empty()) return;  // all required columns present — OK

  std::ostringstream _m;
  _m << "TS_INPUT_FILE '" << fileName << "' is missing columns required by "
     << "FIELD_MODEL=" << model << ":\n";
  for (const auto& n : missing) _m << "  " << n << "\n";
  _m << "Available columns detected in file header:";
  for (const auto& kv : colMap) _m << " " << kv.first;
  exit(__LINE__,__FILE__,_m.str().c_str());
}
//======================================================================================
// ExtractColumn
//======================================================================================
// PURPOSE
//   Safely retrieve one double value from a pre-tokenised data row using the
//   variable name as a lookup key in the TsColMap.
//
// ARGUMENTS
//   toks       — whitespace-split tokens of the current data line (0-based).
//   colMap     — column map produced by ParseTsDriverHeader.
//   varName    — variable name to look up (case-insensitive; ToUpper applied).
//   defaultVal — value returned when the variable is absent from the map or the
//                token at the mapped column cannot be converted to double.
//
// WHEN DEFAULTVAL IS RETURNED
//   1. varName not in colMap — column does not exist in this file.
//   2. Mapped column index is out of range for this row — row is shorter than
//      expected (can happen for the last line if the file is truncated).
//   3. std::stod throws — the token is not a valid floating-point number.
//
//   Cases (1) and (2) are expected only for optional variables (e.g. G1..G3 in
//   a T96 driver file that was also run through a T01 validation check).  For
//   required columns, ValidateTsDriverColumns has already verified their presence,
//   so case (1) should never occur for those variables at this point.
//
// THREAD SAFETY
//   colMap is read-only; toks is a local copy.  Safe to call concurrently.
static double ExtractColumn(const std::vector<std::string>& toks,
                            const TsColMap& colMap,
                            const std::string& varName,
                            double defaultVal = 0.0) {
  const auto it = colMap.find(ToUpper(varName));
  if (it == colMap.end()) return defaultVal;           // variable not in this file
  const int col = it->second;
  if (col < 0 || col >= static_cast<int>(toks.size())) return defaultVal; // row too short
  try { return std::stod(toks[static_cast<std::size_t>(col)]); }
  catch (...) { return defaultVal; }                   // non-numeric token
}
//======================================================================================
// LoadTsDriverFile
//======================================================================================
// PURPOSE
//   Parse a Qin-Denton / ViRBO formatted driver file and return a populated
//   TsDriverTable suitable for per-trajectory-point interpolation of Tsyganenko
//   model driving parameters.
//
// TWO-PASS ALGORITHM
//
//   Pass 1 — Header scan
//     Open the file once and call ParseTsDriverHeader to build a TsColMap.
//     This maps every variable name (uppercase) to its 0-based column index in
//     the data rows.  Close the stream immediately after the header scan.
//
//   Pass 2 — Data rows
//     Re-open the file and iterate over all lines.  For each non-comment,
//     non-empty line:
//       • Strip '#' and '!' comment characters.
//       • Split into whitespace-separated tokens.
//       • Extract the ISO-8601 timestamp from the column identified in Pass 1.
//       • Reject lines whose timestamp token contains no 'T' (ISO-8601 separator),
//         which filters out residual header lines that survived comment stripping.
//       • Convert the timestamp to SPICE ET (exits if SPICE is unavailable).
//       • Use ExtractColumn to read each physics parameter by name from the TsColMap.
//       • Construct a TsDriverRecord and append it to the table.
//
//   After Pass 2, ValidateTsDriverColumns is called between the passes to ensure
//   that every column required by `modelName` is present before any data is read.
//
// ARGUMENTS
//   fileName  — path to the driver file (from TS_INPUT_FILE in #TEMPORAL).
//   modelName — canonical model name (from NormalizeTsModelName); used only to
//               select the correct required-column set for validation.
//
// COLUMN EXTRACTION DETAILS
//   BYIMF / BZIMF     — IMF components [nT]; stored directly.
//   VSW               — solar wind speed [km/s]; the file stores a positive
//                       magnitude, but the solver convention for swVx_kms is
//                       negative (sunward = -X_GSM direction).  We apply
//                       -abs() so the sign is always correct regardless of
//                       whether a particular file happens to include a sign.
//   DEN_P             — proton density [cm^-3].
//   PDYN              — dynamic pressure [nPa].
//   DST               — Dst index [nT].
//   W1..W6            — T05 storm-time integrals; default to 0.0 if not present
//                       in the file (relevant for T96 driver files used with
//                       a T05 run — caught by ValidateTsDriverColumns, but
//                       ExtractColumn's default provides a safe fallback).
//
// TIMESTAMP DETECTION
//   The TsColMap lookup for the timestamp column tries the following variable
//   names in order: ISODATETIME, DATETIME, TIME, DATE_TIME.  If none are found
//   in the header, column 0 is used as a fallback.
//
// THREAD SAFETY
//   This function is called once from ParseAmpsParamFile (single-threaded) and
//   the resulting TsDriverTable is then used read-only by MPI workers.
static EarthUtil::TsDriverTable LoadTsDriverFile(const std::string& fileName,
                                                  const std::string& modelName) {
  // -----------------------------------------------------------------------
  // Pass 1: Build the column map from the file header.
  // -----------------------------------------------------------------------
  TsColMap colMap;
  {
    std::ifstream hdr(fileName);
    if (!hdr.is_open()) {
      std::ostringstream _m; _m << "Cannot open TS_INPUT_FILE: " << fileName;
      exit(__LINE__,__FILE__,_m.str().c_str());
    }
    colMap = ParseTsDriverHeader(hdr, fileName);
    // hdr is closed here (RAII).
  }

  // Validate that every column required by this model is present before
  // reading a single data row.  Fails with a descriptive error listing the
  // missing columns and all available columns from the header.
  ValidateTsDriverColumns(colMap, modelName, fileName);

  // Resolve the 0-based column index for the ISO-8601 timestamp.
  // Try standard ViRBO name variants in priority order.
  int colTime = -1;
  for (const auto& candidate : {"ISODATETIME", "DATETIME", "TIME", "DATE_TIME"}) {
    const auto it = colMap.find(std::string(candidate));
    if (it != colMap.end()) { colTime = it->second; break; }
  }
  if (colTime < 0) {
    // No timestamp column found in header — fall back to column 0 (the standard
    // first column in all known ViRBO files).
    colTime = 0;
    std::cerr << "[TsDriverTable] WARNING: no timestamp column identified in '"
              << fileName << "' header; defaulting to column 0.\n";
  }

  // -----------------------------------------------------------------------
  // Pass 2: Read data rows.
  // -----------------------------------------------------------------------
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    std::ostringstream _m; _m << "Cannot open TS_INPUT_FILE: " << fileName;
    exit(__LINE__,__FILE__,_m.str().c_str());
  }

  EarthUtil::TsDriverTable table;
  std::string line;
  int lineNo = 0, loaded = 0, skipped = 0;

  while (std::getline(fin, line)) {
    ++lineNo;

    // Strip comment text.  In trajectory / driver files both '#' and '!'
    // introduce comments to end-of-line.  We take the earlier delimiter.
    {
      const std::size_t ph = line.find('#');
      const std::size_t pe = line.find('!');
      const std::size_t pcut = std::min(ph, pe);
      if (pcut != std::string::npos) line = line.substr(0, pcut);
    }
    line = Trim(line);
    if (line.empty()) continue;  // blank or comment-only line

    // Tokenise the remaining text on whitespace.
    std::vector<std::string> toks;
    {
      std::istringstream iss(line);
      std::string t;
      while (iss >> t) toks.push_back(t);
    }
    if (toks.empty()) continue;

    // Extract the timestamp from the column identified above.
    const std::string timeUTC =
        (colTime < static_cast<int>(toks.size()))
          ? toks[static_cast<std::size_t>(colTime)]
          : "";

    // Reject lines that do not look like ISO-8601 timestamps (e.g. column-name
    // rows that survived comment stripping because they had no comment prefix).
    // ISO-8601 requires a 'T' separator between date and time.
    if (timeUTC.find('T') == std::string::npos) { ++skipped; continue; }

    // Convert ISO-8601 to SPICE ET (seconds past J2000).
    // We need ET for two purposes:
    //   (a) storing a monotone numeric time for binary-search in TsDriverTable::Lookup,
    //   (b) potential future use in coordinate transforms.
    double et = 0.0;
    if (!ParseIsoUtcToEt(timeUTC, et)) {
#ifdef _NO_SPICE_CALLS_
      // SPICE is required for timestamp conversion.  Without it the table
      // would have et=0 for every record, making Lookup useless.
      exit(__LINE__,__FILE__,
           "LoadTsDriverFile requires SPICE for timestamp conversion "
           "but the build has _NO_SPICE_CALLS_ defined");
#else
      // SPICE is present but this specific string failed.  Skip the row with
      // a warning; do not abort the entire load.
      std::cerr << "[TsDriverTable] WARNING: cannot parse timestamp '"
                << timeUTC << "' at line " << lineNo
                << " of " << fileName << "; skipping.\n";
      ++skipped; continue;
#endif
    }

    // Extract each physics parameter by variable name.
    // ExtractColumn returns defaultVal (0.0) for any column absent from the
    // TsColMap.  For required columns this should never occur here because
    // ValidateTsDriverColumns has already guaranteed their presence.
    EarthUtil::TsDriverRecord rec;
    rec.et       = et;
    rec.timeUTC  = timeUTC;
    rec.imfBy_nT = ExtractColumn(toks, colMap, "BYIMF");
    rec.imfBz_nT = ExtractColumn(toks, colMap, "BZIMF");
    // Vsw in the file is always a positive speed magnitude (e.g. 587.5 km/s).
    // The solver stores swVx_kms as a signed velocity in GSM X: negative
    // because the solar wind flows anti-sunward (-X_GSM).
    rec.swVx_kms = -std::abs(ExtractColumn(toks, colMap, "VSW"));
    rec.swN_cm3  = ExtractColumn(toks, colMap, "DEN_P");
    rec.pdyn_nPa = ExtractColumn(toks, colMap, "PDYN");
    rec.dst_nT   = ExtractColumn(toks, colMap, "DST");
    // W1..W6 are stored in the array rec.w[0..5].
    for (int i = 0; i < 6; ++i) {
      const std::string wKey = "W" + std::to_string(i + 1);
      rec.w[i] = ExtractColumn(toks, colMap, wKey);
    }

    table.push_back(rec);
    ++loaded;
  }

  if (table.empty()) {
    std::ostringstream _m;
    _m << "LoadTsDriverFile: no valid records in '" << fileName
       << "' (lines=" << lineNo << ", skipped=" << skipped << ")";
    exit(__LINE__,__FILE__,_m.str().c_str());
  }

  std::cout << "[TsDriverTable] Loaded " << loaded << " records from '"
            << fileName << "'"
            << (skipped > 0 ? " (" + std::to_string(skipped) + " lines skipped)" : "")
            << "\n";
  return table;
}

AmpsParam ParseAmpsParamFile(const std::string& fileName) {
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    { std::ostringstream _exit_msg; _exit_msg << "Cannot open input file: "+fileName; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
  }

  AmpsParam p;

  std::string section;
  bool inPointsBlock=false;
  bool inChannelsBlock=false;

  std::string line;
  int lineNo=0;
  while (std::getline(fin,line)) {
    ++lineNo;

    // Preserve the original line because unit hints may be encoded in the
    // trailing comment after '!'. Example: DOMAIN_X_MAX 20 ! Re
    const std::string rawLine=line;
    std::string commentText;
    SplitCodeAndComment(rawLine,line,commentText);
    line=Trim(line);
    commentText=Trim(commentText);
    if (line.empty()) continue;

    // Section header
    if (line.size()>1 && line[0]=='#') {
      section=ToUpper(Trim(line));
      continue;
    }

    // Block delimiters for POINTS
    if (ToUpper(line)=="POINTS_BEGIN") { inPointsBlock=true; continue; }
    if (ToUpper(line)=="POINTS_END") { inPointsBlock=false; continue; }

    // Block delimiters for ENERGY_CHANNELS
    if (ToUpper(line)=="CH_BEGIN") { inChannelsBlock=true; continue; }
    if (ToUpper(line)=="CH_END")   { inChannelsBlock=false; continue; }

    if (inPointsBlock) {
      // Accept both the original bare format "x y z" and the newer "POINT x y z"
      // form (which was added to allow inline units on each token, e.g. POINT 1Re 0 0).
      // The presence of the "POINT" keyword is detected by trying to read the first
      // token: if it converts to a double it is the x-coordinate; otherwise it must
      // be the literal string "POINT" (anything else is an error).
      std::istringstream iss(line);
      std::string firstTok;
      iss >> firstTok;

      std::string xTok, yTok, zTok;
      if (ToUpper(firstTok) == "POINT") {
        // New format: "POINT x y z" — first token is the keyword.
        // Coordinates may carry inline units (e.g. POINT 1Re 0 0) or inherit
        // units from the trailing comment (e.g. POINT 1 0 0 ! Re).
        if (!(iss >> xTok >> yTok >> zTok)) {
          { std::ostringstream _exit_msg; _exit_msg << "Cannot parse POINT coordinates at line "+std::to_string(lineNo)+": "+line; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
        }
      } else {
        // Original format: "x y z" — first token is already the x-coordinate.
        xTok = firstTok;
        if (!(iss >> yTok >> zTok)) {
          { std::ostringstream _exit_msg; _exit_msg << "Cannot parse point coordinates at line "+std::to_string(lineNo)+": "+line; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
        }
      }

      Vec3 v;
      {
        std::vector<double> xyz = ParseLengthListToKm(xTok+" "+yTok+" "+zTok, commentText);
        if (xyz.size()!=3) { std::ostringstream _exit_msg; _exit_msg << "Internal error parsing POINT at line "+std::to_string(lineNo); exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
        v.x=xyz[0]; v.y=xyz[1]; v.z=xyz[2];
      }
      p.output.points.push_back(v);
      continue;
    }

    // ----- ENERGY_CHANNELS block: each line is  NAME  E1_MeV  E2_MeV  -----
    // PHYSICS:
    //   Each channel defines an energy sub-range for integral-flux output.
    //   The solver will integrate  F_ch = 4pi * integral[E1,E2] T(E)*Jb(E) dE.
    // VALIDATION:
    //   E1 and E2 are checked here; cross-validation against DS_EMIN/DS_EMAX
    //   is left to the solver (channels that don't overlap the grid give F=0).
    if (inChannelsBlock) {
      std::istringstream iss(line);
      EnergyChannel ch;
      if (!(iss >> ch.name >> ch.E1_MeV >> ch.E2_MeV)) {
        { std::ostringstream _exit_msg; _exit_msg << "Malformed CH line at line "
            +std::to_string(lineNo)+": expected NAME E1_MeV E2_MeV, got: "+line;
          exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
      }
      if (ch.E1_MeV <= 0.0) {
        { std::ostringstream _exit_msg; _exit_msg << "Energy channel '"+ch.name
            +"' at line "+std::to_string(lineNo)+": E1_MeV must be > 0";
          exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
      }
      if (ch.E2_MeV <= ch.E1_MeV) {
        { std::ostringstream _exit_msg; _exit_msg << "Energy channel '"+ch.name
            +"' at line "+std::to_string(lineNo)+": E2_MeV must be > E1_MeV";
          exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
      }
      p.fluxChannels.push_back(ch);
      continue;
    }

    // Generic KEY VALUE line
    std::string key,val;
    SplitKV(line,key,val);
    if (key.empty()) continue;
    std::string uKey=ToUpper(key);

    auto rememberUnknown=[&](){
      p.unknown[section+":"+uKey]=val;
    };

    if (section=="#RUN_INFO") {
      if (uKey=="RUN_ID") p.runId=val;
      else rememberUnknown();
    }
    else if (section=="#CALCULATION_MODE") {
      if (uKey=="CALC_TARGET") p.calc.target=ToUpper(val);
      else if (uKey=="FIELD_EVAL_METHOD") p.calc.fieldEvalMethod=ToUpper(val);
      else rememberUnknown();
    }
    else if (section=="#CUTOFF_RIGIDITY") {
      if (uKey=="CUTOFF_EMIN") p.cutoff.eMin_MeV=std::stod(val);
      else if (uKey=="CUTOFF_EMAX") p.cutoff.eMax_MeV=std::stod(val);
      else if (uKey=="CUTOFF_NENERGY") p.cutoff.nEnergy=std::stoi(val);
      else if (uKey=="CUTOFF_MAX_PARTICLES") p.cutoff.maxParticlesPerPoint=std::stoi(val);

      // Optional per-trajectory integration cap for the cutoff search.
      // If omitted, the solver will fall back to #NUMERICAL MAX_TRACE_TIME.
      else if (uKey=="CUTOFF_MAX_TRAJ_TIME") p.cutoff.maxTrajTime_s=std::stod(val);

      // Cutoff sampling strategy: VERTICAL or ISOTROPIC.
      // We store the string in uppercase for robust comparisons later.
      else if (uKey=="CUTOFF_SAMPLING") p.cutoff.sampling=ToUpper(val);

      // Directional Rc sky-map controls.
      // Notes:
      //   - These are used only in the gridless cutoff solver.
      //   - A directional map is a diagnostic product; it can be expensive.
      else if (uKey=="DIRECTIONAL_MAP") p.cutoff.directionalMap=ToBool(val);
      else if (uKey=="DIRMAP_LON_RES") p.cutoff.dirMapLonRes_deg=std::stod(val);
      else if (uKey=="DIRMAP_LAT_RES") p.cutoff.dirMapLatRes_deg=std::stod(val);
      else rememberUnknown();
    }
    else if (section=="#PARTICLE_SPECIES") {
      if (uKey=="SPECIES") p.species.name=ToUpper(val);
      else if (uKey=="CHARGE") p.species.charge_e=std::stoi(val);
      else if (uKey=="MASS_AMU") p.species.mass_amu=std::stod(val);
      else rememberUnknown();
    }
    else if (section=="#BACKGROUND_FIELD") {
      if (uKey=="FIELD_MODEL") {
        // Normalize alternate spellings to the canonical model name
        // used throughout the codebase (e.g. TS05 -> T05).
        p.field.model = NormalizeTsModelName(val);
      }
      else if (uKey=="DIPOLE_MOMENT") p.field.dipoleMoment_Me=std::stod(val);
      else if (uKey=="DIPOLE_TILT" || uKey=="DIPOLE_TILT_DEG") p.field.dipoleTilt_deg=std::stod(val);
      else if (uKey=="DST") p.field.dst_nT=std::stod(val);
      else if (uKey=="PDYN") p.field.pdyn_nPa=std::stod(val);
      else if (uKey=="IMF_BY") p.field.imfBy_nT=std::stod(val);
      else if (uKey=="IMF_BZ") p.field.imfBz_nT=std::stod(val);
      else if (uKey=="IMF_BX") p.field.imfBx_nT=std::stod(val);
      else if (uKey=="SW_VX") p.field.swVx_kms=std::stod(val);
      else if (uKey=="SW_N") p.field.swN_cm3=std::stod(val);
      else if (uKey=="W1") p.field.w[0]=std::stod(val);
      else if (uKey=="W2") p.field.w[1]=std::stod(val);
      else if (uKey=="W3") p.field.w[2]=std::stod(val);
      else if (uKey=="W4") p.field.w[3]=std::stod(val);
      else if (uKey=="W5") p.field.w[4]=std::stod(val);
      else if (uKey=="W6") p.field.w[5]=std::stod(val);
      else if (uKey=="EPOCH") p.field.epoch=val;
      else {
        p.field.raw[uKey]=val;
        rememberUnknown();
      }
    }
    else if (section=="#ELECTRIC_FIELD") {
      if (uKey=="EFIELD_MODEL") p.efield.model=ToUpper(val);
      else if (uKey=="COROTATION_SCALE") p.efield.corotationScale=std::stod(val);
      else if (uKey=="VS_POTENTIAL_KV") p.efield.vsPotential_kV=std::stod(val);
      else if (uKey=="VS_GAMMA") p.efield.vsGamma=std::stod(val);
      else if (uKey=="VS_REFERENCE_L") p.efield.vsReferenceL=std::stod(val);
      else if (uKey=="VS_SCALE") p.efield.vsScale=std::stod(val);
      else if (uKey=="EFIELD_RMIN" || uKey=="EFIELD_RMIN_KM") p.efield.rMin_km=ParseLengthToKm(val,commentText);
      else if (uKey=="EFIELD_LMIN") p.efield.lMin=std::stod(val);
      else {
        p.efield.raw[uKey]=val;
        rememberUnknown();
      }
    }
    else if (section=="#DOMAIN_BOUNDARY") {
      if (uKey=="DOMAIN_X_MAX") p.domain.xMax=ParseLengthToKm(val,commentText);
      else if (uKey=="DOMAIN_X_MIN") p.domain.xMin=ParseLengthToKm(val,commentText);
      else if (uKey=="DOMAIN_Y_MAX") p.domain.yMax=ParseLengthToKm(val,commentText);
      else if (uKey=="DOMAIN_Y_MIN") p.domain.yMin=ParseLengthToKm(val,commentText);
      else if (uKey=="DOMAIN_Z_MAX") p.domain.zMax=ParseLengthToKm(val,commentText);
      else if (uKey=="DOMAIN_Z_MIN") p.domain.zMin=ParseLengthToKm(val,commentText);
      else if (uKey=="R_INNER") p.domain.rInner=ParseLengthToKm(val,commentText);
      else rememberUnknown();
    }
    else if (section=="#OUTPUT_DOMAIN") {
      if (uKey=="OUTPUT_MODE") p.output.mode=ToUpper(val);
      else if (uKey=="OUTPUT_COORDS" || uKey=="COORDS") p.output.coords=ToUpper(val);
      else if (uKey=="TRAJ_FRAME") p.output.trajFrame=ToUpper(val);
      else if (uKey=="TRAJ_FILE") p.output.trajFile=Trim(val);
      else if (uKey=="FLUX_DT") p.output.fluxDt_min=std::stod(val);
      else if (uKey=="SHELL_ALTS_KM" || uKey=="SHELL_ALTITUDES") {
        // Although the key name contains _KM (historical), accept explicit Re
        // units inline or in the trailing comment. If no units are provided, km
        // remains the default for backward compatibility.
	auto t=ParseLengthListToKm(val,commentText);
        p.output.shellAlt_km.insert(p.output.shellAlt_km.end(),t.begin(),t.end());  
      }
      else if (uKey=="SHELL_RES_DEG") p.output.shellRes_deg=std::stod(val);
      else {
        p.output.raw[uKey]=val;
        rememberUnknown();
      }
    }
    else if (section=="#NUMERICAL") {
      if (uKey=="DT_TRACE") p.numerics.dtTrace_s=std::stod(val);
      else if (uKey=="MAX_STEPS") p.numerics.maxSteps=std::stoi(val);
      else if (uKey=="MAX_TRACE_TIME") p.numerics.maxTraceTime_s=std::stod(val);
      else if (uKey=="MAX_TRACE_DISTANCE") p.numerics.maxTraceDistance_Re=std::stod(val);
      else rememberUnknown();
    }
    else if (section=="#DENSITY_SPECTRUM") {
      // Density/spectrum workflow controls.
      // All energies are read in MeV (commonly MeV/n in CCMC inputs).
      if (uKey=="DS_EMIN") p.densitySpectrum.Emin_MeV=std::stod(val);
      else if (uKey=="DS_EMAX") p.densitySpectrum.Emax_MeV=std::stod(val);
      else if (uKey=="DS_NINTERVALS") p.densitySpectrum.nIntervals=std::stoi(val);
      else if (uKey=="DS_MAX_PARTICLES") p.densitySpectrum.maxParticlesPerPoint=std::stoi(val);
      else if (uKey=="DS_MAX_TRAJ_TIME") p.densitySpectrum.maxTrajTime_s=std::stod(val);
      else if (uKey=="DS_ENERGY_SPACING") p.densitySpectrum.spacing=ParseEnergySpacingToken(val);
      else if (uKey=="DS_BOUNDARY_MODE") p.densitySpectrum.boundaryMode=ToUpper(val);
      else rememberUnknown();
    }
    else if (section=="#SPECTRUM") {
      p.spectrum[uKey]=val;
    }
    else if (section=="#BOUNDARY_ANISOTROPY") {
      if      (uKey=="BA_PAD_MODEL")        p.anisotropy.padModel        = ToUpper(val);
      else if (uKey=="BA_PAD_EXPONENT")     p.anisotropy.padExponent     = std::stod(val);
      else if (uKey=="BA_SPATIAL_MODEL")    p.anisotropy.spatialModel    = ToUpper(val);
      else if (uKey=="BA_DAYSIDE_FACTOR")   p.anisotropy.daysideFactor   = std::stod(val);
      else if (uKey=="BA_NIGHTSIDE_FACTOR") p.anisotropy.nightsideFactor = std::stod(val);
      else rememberUnknown();
    }
    else if (section=="#OUTPUT_OPTIONS") {
      p.outputOptions[uKey]=val;
    }
    else if (section=="#TEMPORAL") {
      if      (uKey=="TEMPORAL_MODE")    p.temporal.mode=ToUpper(val);
      else if (uKey=="EVENT_START")      p.temporal.eventStart=Trim(val);
      else if (uKey=="EVENT_END")        p.temporal.eventEnd=Trim(val);
      else if (uKey=="FIELD_UPDATE_DT")  p.temporal.fieldUpdateDt_min=std::stod(val);
      else if (uKey=="INJECT_DT")        p.temporal.injectDt_min=std::stod(val);
      else if (uKey=="TS_INPUT_MODE")    p.temporal.tsInputMode=ToUpper(val);
      else if (uKey=="TS_INPUT_FILE")    p.temporal.tsInputFile=Trim(val);
      else rememberUnknown();
    }
    else {
      rememberUnknown();
    }
  }

  

// Validate dipole-specific background field settings if requested.
if (ToUpper(p.field.model)=="DIPOLE") {
  // Moment scaling must be positive.
  if (!(p.field.dipoleMoment_Me>0.0)) {
    exit(__LINE__,__FILE__,"DIPOLE_MOMENT must be > 0 (multiples of M_E)");
  }
  // Keep tilt range conservative; values outside +/-90 deg usually indicate
  // a coordinate/sign convention mistake in the input.
  if (p.field.dipoleTilt_deg < -90.0 || p.field.dipoleTilt_deg > 90.0) {
    exit(__LINE__,__FILE__,"DIPOLE_TILT must be in [-90, 90] degrees");
  }
}

  // Finalize the point-like initialization path.  This is the key step that makes
  // POINTS and TRAJECTORY reuse the same downstream solver kernels:
  //   POINTS     -> one synthetic one-sample trajectory per point
  //   TRAJECTORY -> real multi-sample trajectory loaded from file
  // In both cases the result is a flattened output.points array plus the unified
  // output.trajectories container used for per-sample timestamps.
  InitializePointLikeOutput(p);

  // Load time-varying Tsyganenko driver file when TS_INPUT_MODE = FILE.
  // This must happen after the #TEMPORAL section has been fully parsed and
  // before the solvers start querying the table per trajectory point.
  if (ToUpper(p.temporal.tsInputMode) == "FILE") {
    if (Trim(p.temporal.tsInputFile).empty()) {
      exit(__LINE__,__FILE__,
           "TS_INPUT_MODE=FILE requires TS_INPUT_FILE to be set in #TEMPORAL");
    }
    p.temporal.driverTable = LoadTsDriverFile(
        p.temporal.tsInputFile,
        p.field.model);  // model name drives per-model column validation
  }

  // Validate density/spectrum controls if requested.
  if (ToUpper(p.calc.target)=="DENSITY_SPECTRUM") {
    if (!(p.densitySpectrum.Emin_MeV>0.0)) {
      exit(__LINE__,__FILE__,"DS_EMIN must be > 0 (MeV)");
    }
    if (!(p.densitySpectrum.Emax_MeV>p.densitySpectrum.Emin_MeV)) {
      exit(__LINE__,__FILE__,"DS_EMAX must be > DS_EMIN (MeV)");
    }
    if (!(p.densitySpectrum.nIntervals>=1)) {
      exit(__LINE__,__FILE__,"DS_NINTERVALS must be >= 1");
    }
    if (p.densitySpectrum.maxParticlesPerPoint < 0) {
      exit(__LINE__,__FILE__,"DS_MAX_PARTICLES must be >= 0 (0 means: no cap)");
    }
    if (p.densitySpectrum.maxTrajTime_s < 0.0) {
      exit(__LINE__,__FILE__,"DS_MAX_TRAJ_TIME must be >= 0 (0 means: use MAX_TRACE_TIME)");
    }
    if (p.numerics.maxTraceDistance_Re < 0.0) {
      exit(__LINE__,__FILE__,"MAX_TRACE_DISTANCE must be >= 0 (0 means: disabled)");
    }
    if (ToUpper(p.calc.fieldEvalMethod)!="GRIDLESS") {
      exit(__LINE__,__FILE__,"DENSITY_SPECTRUM currently requires FIELD_EVAL_METHOD = GRIDLESS");
    }
    // Validate DS_BOUNDARY_MODE token.
    const std::string bm = ToUpper(p.densitySpectrum.boundaryMode);
    if (bm!="ISOTROPIC" && bm!="ANISOTROPIC") {
      { std::ostringstream _exit_msg; _exit_msg << "DS_BOUNDARY_MODE must be ISOTROPIC or ANISOTROPIC (got '"+p.densitySpectrum.boundaryMode+"')"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
    }
    // When ANISOTROPIC is requested, validate the anisotropy sub-parameters.
    if (bm=="ANISOTROPIC") {
      const std::string pm = ToUpper(p.anisotropy.padModel);
      if (pm!="ISOTROPIC" && pm!="SINALPHA_N" && pm!="COSALPHA_N" && pm!="BIDIRECTIONAL") {
        { std::ostringstream _exit_msg; _exit_msg << "BA_PAD_MODEL must be ISOTROPIC|SINALPHA_N|COSALPHA_N|BIDIRECTIONAL (got '"+p.anisotropy.padModel+"')"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
      }
      if (p.anisotropy.padExponent < 0.0) {
        exit(__LINE__,__FILE__,"BA_PAD_EXPONENT must be >= 0");
      }
      const std::string sm = ToUpper(p.anisotropy.spatialModel);
      if (sm!="UNIFORM" && sm!="DAYSIDE_NIGHTSIDE") {
        { std::ostringstream _exit_msg; _exit_msg << "BA_SPATIAL_MODEL must be UNIFORM|DAYSIDE_NIGHTSIDE (got '"+p.anisotropy.spatialModel+"')"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }
      }
      if (p.anisotropy.daysideFactor < 0.0 || p.anisotropy.nightsideFactor < 0.0) {
        exit(__LINE__,__FILE__,"BA_DAYSIDE_FACTOR and BA_NIGHTSIDE_FACTOR must be >= 0");
      }
    }
  }

  // Parse the particle spectrum into a typed metadata object first.  The legacy
  // boundary/spectrum initializer still consumes the raw key/value table, but by
  // storing the typed view in AmpsParam we make the selected spectrum explicit in
  // the model initialization sequence.
  p.particleSpectrum = ParseParticleSpectrum(p.spectrum);

  // TABLE spectra support two on-disk formats:
  //   (1) a legacy 2-column table (Energy_MeV, differential flux)
  //   (2) a time-dependent table with one UTC timestamp per row and many energy columns
  // The parser determines the file type here and passes the reference epoch through
  // #SPECTRUM so cSpectrum can call the appropriate loader directly.
  ResolveTableSpectrumFileForInitialization(p);

  // Build/validate the global spectrum object from the prepared #SPECTRUM table because
  // the downstream density/cutoff code still evaluates spectra through ::gSpectrum.
  InitGlobalSpectrumFromKeyValueMap(p.spectrum);

  // Diagnostic output: record the parsed spectrum in a standalone Tecplot file.
  // This is intentionally written immediately after spectrum initialization so the
  // output reflects exactly what injection will use.
  //
  // File: spectrum_input.dat
  //  - Energy is written in MeV
  //  - TABLE spectra preserve the user table points
  //  - Analytic spectra are sampled with log-spaced energies (200 points by default)
  ::WriteSpectrumInputTecplot("spectrum_input.dat", ::gSpectrum);

  return p;
}

}

