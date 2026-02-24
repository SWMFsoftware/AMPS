//======================================================================================
// amps_param_parser.cpp
//======================================================================================
// See amps_param_parser.h for a detailed description.
//
// DESIGN CHOICES
//   1) Minimal dependencies: uses only the C++ standard library.
//   2) Permissive parsing: unknown keys are retained rather than rejected.
//   3) Robust comment handling: text after '!' is usually ignored, except that
//      recognized unit hints (km/Re) may be consumed by length parsers.
//   4) Section-sensitive interpretation: the same key in different sections
//      can mean different things; we map only the keys we need.
//
// LIMITATIONS (CURRENT)
//   - GEO/GSE coordinate conversion is NOT performed here.
//     The gridless cutoff solver currently assumes positions are in GSM.
//     We keep the coordinate tag for reporting and future upgrades.
//   - POINT lines are parsed as 3 floating point numbers.
//======================================================================================

#include "amps_param_parser.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <vector>
#include <cstdlib>
#include <cerrno>

namespace EarthUtil {

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
  throw std::runtime_error("Cannot parse boolean token: '"+sIn+"'");
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
  if (t.empty()) throw std::runtime_error("Empty length token");

  const char* begin=t.c_str();
  char* end=nullptr;
  errno=0;
  double value=std::strtod(begin,&end);
  if (begin==end) {
    throw std::runtime_error("Cannot parse numeric length token: '"+t+"'");
  }
  if (errno==ERANGE) {
    throw std::runtime_error("Length token is out of range: '"+t+"'");
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
    throw std::runtime_error("Unsupported inline length unit '"+suffix+"' in token '"+t+"'");
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

  std::vector<std::string> unitHints;
  {
    std::istringstream iss(comment);
    std::string t;
    while (iss>>t) {
      std::string u=ToLower(Trim(t));
      if (u=="km" || u=="re") unitHints.push_back(u);
    }
  }

  std::vector<double> out;
  out.reserve(valueToks.size());

  const bool oneGlobalHint = (unitHints.size()==1);
  const bool perValueHints = (unitHints.size()==valueToks.size());

  for (size_t i=0;i<valueToks.size();++i) {
    std::string hint;
    if (perValueHints) hint=unitHints[i];
    else if (oneGlobalHint) hint=unitHints[0];
    out.push_back(ParseLengthToKm(valueToks[i],hint));
  }
  return out;
}

AmpsParam ParseAmpsParamFile(const std::string& fileName) {
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    throw std::runtime_error("Cannot open input file: "+fileName);
  }

  AmpsParam p;

  std::string section;
  bool inPointsBlock=false;

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

    if (inPointsBlock) {
      std::istringstream iss(line);
      std::string tok; iss >> tok;
      if (ToUpper(tok)!="POINT") {
        throw std::runtime_error("Malformed POINTS block at line "+std::to_string(lineNo)+": "+line);
      }
      // POINT coordinates are length-like values. They may have inline units
      // (e.g. POINT 1Re 0 0) or units specified in the comment after '!':
      //   POINT 1 0 0 ! Re
      //   POINT 6371 0 0 ! km km km
      std::string xTok,yTok,zTok;
      if (!(iss >> xTok >> yTok >> zTok)) {
        throw std::runtime_error("Cannot parse POINT coordinates at line "+std::to_string(lineNo)+": "+line);
      }
      Vec3 v;
      {
        std::vector<double> xyz = ParseLengthListToKm(xTok+" "+yTok+" "+zTok, commentText);
        if (xyz.size()!=3) throw std::runtime_error("Internal error parsing POINT at line "+std::to_string(lineNo));
        v.x=xyz[0]; v.y=xyz[1]; v.z=xyz[2];
      }
      p.output.points.push_back(v);
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
      else if (uKey=="OUTPUT_COORDS") p.output.coords=ToUpper(val);
      else if (uKey=="SHELL_ALTS_KM") {
        // Although the key name contains _KM (historical), accept explicit Re
        // units inline or in the trailing comment. If no units are provided, km
        // remains the default for backward compatibility.
        p.output.shellAlt_km = ParseLengthListToKm(val,commentText);
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
      else rememberUnknown();
    }
    else if (section=="#SPECTRUM") {
      p.spectrum[uKey]=val;
    }
    else if (section=="#OUTPUT_OPTIONS") {
      p.outputOptions[uKey]=val;
    }
    else {
      rememberUnknown();
    }
  }

  if (p.output.mode=="POINTS" && p.output.points.empty()) {
    throw std::runtime_error("OUTPUT_MODE=POINTS but no POINT entries were found in POINTS_BEGIN/END block");
  }

  return p;
}

}

