//======================================================================================
// amps_param_parser.cpp
//======================================================================================
// See amps_param_parser.h for a detailed description.
//
// DESIGN CHOICES
//   1) Minimal dependencies: uses only the C++ standard library.
//   2) Permissive parsing: unknown keys are retained rather than rejected.
//   3) Robust comment stripping: anything after '!' is ignored.
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
    line=StripComment(line);
    line=Trim(line);
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
      Vec3 v;
      if (!(iss >> v.x >> v.y >> v.z)) {
        throw std::runtime_error("Cannot parse POINT coordinates at line "+std::to_string(lineNo)+": "+line);
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
      if (uKey=="DOMAIN_X_MAX") p.domain.xMax=std::stod(val);
      else if (uKey=="DOMAIN_X_MIN") p.domain.xMin=std::stod(val);
      else if (uKey=="DOMAIN_Y_MAX") p.domain.yMax=std::stod(val);
      else if (uKey=="DOMAIN_Y_MIN") p.domain.yMin=std::stod(val);
      else if (uKey=="DOMAIN_Z_MAX") p.domain.zMax=std::stod(val);
      else if (uKey=="DOMAIN_Z_MIN") p.domain.zMin=std::stod(val);
      else if (uKey=="R_INNER") p.domain.rInner=std::stod(val);
      else rememberUnknown();
    }
    else if (section=="#OUTPUT_DOMAIN") {
      if (uKey=="OUTPUT_MODE") p.output.mode=ToUpper(val);
      else if (uKey=="OUTPUT_COORDS") p.output.coords=ToUpper(val);
      else if (uKey=="SHELL_ALTS_KM") {
        std::istringstream iss(val);
        double a;
        p.output.shellAlt_km.clear();
        while (iss >> a) p.output.shellAlt_km.push_back(a);
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
