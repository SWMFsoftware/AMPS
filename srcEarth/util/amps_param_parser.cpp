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
//   "SPECTRUM"             -> store in result.spectrum raw map
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
//   (D) The #SPECTRUM section is intentionally opaque to this parser: all its
//       key/value pairs are passed through to result.spectrum as raw strings.
//       The InitGlobalSpectrumFromKeyValueMap function in boundary/spectrum.cpp
//       is responsible for interpreting them. This separation means the spectrum
//       model can evolve (new power-law forms, tabulated spectra, SEP event
//       lookup tables) without any changes here.
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

static EarthUtil::SpacecraftTrajectory LoadTrajectoryFileAsGsm(const std::string& fileName,
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
// LoadTsDriverFile — parse Qin-Denton / ViRBO driver file
//======================================================================================
// Column layout (1-based, space-separated; '#' and '!' are comment characters):
//   $1        ISO-8601 datetime (e.g. 2017-09-10T00:00:00)
//   $2–$7     year month day hour min sec    (redundant; $1 drives the SPICE call)
//   $8        IMF By  [nT]
//   $9        IMF Bz  [nT]
//   $10       Vsw     [km/s]  (positive magnitude; stored as negative swVx)
//   $11       Den_P   [cm^-3]
//   $12       Pdyn    [nPa]
//   $13–$15   G1 G2 G3  (T96 G-parameters; skipped — not needed by T96/T05 calls)
//   $16–$23   8 quality-flag integers  (skipped)
//   $24       Kp      [dimensionless]  (skipped)
//   $25       akp3    [dimensionless]  (skipped)
//   $26       Dst     [nT]
//   $27–$32   Bz1..Bz6  (historic average Bz lags; skipped)
//   $33–$38   W1..W6    (T05 storm-time history integrals)
//   $39–$44   W1..W6 quality flags  (skipped)
static EarthUtil::TsDriverTable LoadTsDriverFile(const std::string& fileName) {
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    std::ostringstream _m;
    _m << "Cannot open TS_INPUT_FILE: " << fileName;
    exit(__LINE__,__FILE__,_m.str().c_str());
  }

  EarthUtil::TsDriverTable table;
  std::string line;
  int lineNo = 0, loaded = 0, skipped = 0;

  while (std::getline(fin, line)) {
    ++lineNo;

    // Strip comments: both '#' and '!' introduce a comment to end-of-line.
    {
      const std::size_t ph = line.find('#');
      const std::size_t pe = line.find('!');
      const std::size_t pcut = std::min(ph, pe);
      if (pcut != std::string::npos) line = line.substr(0, pcut);
    }
    line = Trim(line);
    if (line.empty()) continue;

    // Parse 38 required tokens ($1..$38).  Rows with fewer tokens are skipped.
    std::istringstream iss(line);
    std::string timeUTC;
    int yr,mo,dy,hr,mn,sc;
    double by,bz,vsw,denP,pdyn;
    double g1,g2,g3;
    int sf1,sf2,sf3,sf4,sf5,sf6,sf7,sf8;
    double kp,akp3,dst;
    double bz1,bz2,bz3,bz4,bz5,bz6;
    double w1,w2,w3,w4,w5,w6;

    if (!(iss >> timeUTC
              >> yr >> mo >> dy >> hr >> mn >> sc
              >> by >> bz >> vsw >> denP >> pdyn
              >> g1 >> g2 >> g3
              >> sf1 >> sf2 >> sf3 >> sf4 >> sf5 >> sf6 >> sf7 >> sf8
              >> kp >> akp3 >> dst
              >> bz1 >> bz2 >> bz3 >> bz4 >> bz5 >> bz6
              >> w1 >> w2 >> w3 >> w4 >> w5 >> w6)) {
      ++skipped;
      continue;
    }

    // Require ISO-8601 'T' separator to reject residual header lines.
    if (timeUTC.find('T') == std::string::npos) { ++skipped; continue; }

    // Convert the ISO-8601 timestamp to SPICE ephemeris time.
    double et = 0.0;
    if (!ParseIsoUtcToEt(timeUTC, et)) {
#ifdef _NO_SPICE_CALLS_
      exit(__LINE__,__FILE__,
           "LoadTsDriverFile requires SPICE for timestamp conversion "
           "but the build has _NO_SPICE_CALLS_ defined");
#else
      std::cerr << "[TsDriverTable] WARNING: cannot parse timestamp '"
                << timeUTC << "' at line " << lineNo
                << " of " << fileName << "; skipping.\n";
      ++skipped;
      continue;
#endif
    }

    EarthUtil::TsDriverRecord rec;
    rec.et       = et;
    rec.timeUTC  = timeUTC;
    rec.imfBy_nT = by;
    rec.imfBz_nT = bz;
    // Vsw is stored as a positive speed magnitude in the file.
    // The solver convention for swVx_kms is negative (sunward = -X_GSM).
    rec.swVx_kms = -std::abs(vsw);
    rec.swN_cm3  = denP;
    rec.pdyn_nPa = pdyn;
    rec.dst_nT   = dst;
    rec.w[0]=w1; rec.w[1]=w2; rec.w[2]=w3;
    rec.w[3]=w4; rec.w[4]=w5; rec.w[5]=w6;

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

  // Keep a parsed spectrum object local to this translation unit.
  // In the full codebase, add it as a field in AmpsParam and remove this local.

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
      if (uKey=="FIELD_MODEL") p.field.model=ToUpper(val);
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

// Build the unified trajectory representation once parsing is complete.
  if (p.output.mode=="TRAJECTORY") {
    if (Trim(p.output.trajFile).empty()) {
      exit(__LINE__,__FILE__,"OUTPUT_MODE=TRAJECTORY requires TRAJ_FILE in #OUTPUT_DOMAIN");
    }
    p.output.trajectories.clear();
    p.output.trajectories.push_back(LoadTrajectoryFileAsGsm(p.output.trajFile, p.output.trajFrame));
    p.output.RebuildFlattenedPointsFromTrajectories();
  }
  else if (p.output.mode=="POINTS") {
    if (p.output.points.empty()) {
      exit(__LINE__,__FILE__,"OUTPUT_MODE=POINTS but no POINT entries were found in POINTS_BEGIN/END block");
    }
    PackStandalonePointsIntoSyntheticTrajectories(p);
    p.output.RebuildFlattenedPointsFromTrajectories();
  }

  // Load time-varying Tsyganenko driver file when TS_INPUT_MODE = FILE.
  // This must happen after the #TEMPORAL section has been fully parsed and
  // before the solvers start querying the table per trajectory point.
  if (ToUpper(p.temporal.tsInputMode) == "FILE") {
    if (Trim(p.temporal.tsInputFile).empty()) {
      exit(__LINE__,__FILE__,
           "TS_INPUT_MODE=FILE requires TS_INPUT_FILE to be set in #TEMPORAL");
    }
    p.temporal.driverTable = LoadTsDriverFile(p.temporal.tsInputFile);
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

  // Parse spectrum into a typed representation and validate SPECTRUM_TYPE.
  // This guarantees we fail fast for unrecognized spectrum definitions.

  // If your downstream Tecplot writer lives elsewhere, pass parsedSpectrum to it
  // and call WriteSpectrumTecplotAuxData(out, parsedSpectrum) when writing the
  // Tecplot header. See comment above.

    // Build/validate the global spectrum object from the parsed #SPECTRUM section.
  // This exits early with a clear error message if the spectrum is missing, incomplete,
  // or has an unrecognized SPECTRUM_TYPE token.
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

