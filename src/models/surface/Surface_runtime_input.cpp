//$Id$

/*
 * Surface_runtime_input.cpp
 *
 * The runtime input block 'physics: gas_surface_interaction' of the gas-surface interaction models (the vocabulary is
 * described in Surface.h, the grammar in src/INPUT_FILE.md).  The block is owned by this library: an application that
 * reads a runtime input file (AMPS::Input::ParseFile) hands the resolved block to Surface::RuntimeInput::Apply().
 * Only the tree types and helpers of amps_input.h are used here (no file access).
 */

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>

#include "Surface.h"
#include "src/general/input/amps_input.h"

namespace {
  namespace In=AMPS::Input;

  //the per-species keys of each model: one set, or alternative sets (the Nocilla model: the direct parameters or the
  //accommodation coefficients, Surface.h); a species gives exactly one complete set
  std::vector<std::vector<std::string> > ModelKeySets(int model) {
    switch (model) {
    case Surface::Model::Maxwell: return {{"accommodation"}};
    case Surface::Model::CLL: case Surface::Model::QuasiSpecularCLL: return {{"normal_energy_accommodation","tangential_momentum_accommodation"}};
    case Surface::Model::QuasiSpecularNocilla: return {{"drift_speed_m_s","drift_angle_deg","temperature_k"},
        {"normal_momentum_accommodation","tangential_momentum_accommodation","energy_accommodation"}};
    case Surface::Model::QuasiSpecularSchamberg: return {{"energy_accommodation","nu","beam_half_angle_deg"}};
    }

    return {{}};
  }

  std::string KeyList(const std::vector<std::string> &keys) {
    std::string list;

    for (const std::string &p : keys) list+=(list.empty() ? "" : ", ")+p;
    return list;
  }

  bool UsesWallTemperature(int model) {
    switch (model) {
    case Surface::Model::Diffuse: case Surface::Model::Maxwell: case Surface::Model::CLL: case Surface::Model::QuasiSpecularCLL:
    case Surface::Model::QuasiSpecularSchamberg:
      return true;
    }

    return false;
  }

  bool Fail(const In::Assignment &a,const std::string &reason,std::string *error) {
    if (error!=NULL) *error=In::ErrorAt(a.where,reason,a.text).message;
    return false;
  }

  bool Fail(const In::Block &b,const std::string &reason,std::string *error) {
    if (error!=NULL) *error=In::ErrorAt(b.begin,reason).message;
    return false;
  }

  //a number: the whole value must be one floating-point literal ('inf' is accepted only when allowInf)
  bool ParseNumber(const In::Assignment &a,double *x,bool allowInf,std::string *error) {
    const char *s=a.value.c_str();
    char *end=NULL;

    errno=0;
    *x=strtod(s,&end);
    while ((end!=NULL)&&((*end==' ')||(*end=='\t'))) end++;

    if ((end==s)||(end==NULL)||(*end!='\0')||(errno==ERANGE)||std::isnan(*x)) {
      return Fail(a,"key '"+a.key+"': '"+a.value+"' is not a number",error);
    }

    if (std::isinf(*x)&&(!allowInf)) return Fail(a,"key '"+a.key+"': the value must be finite",error);
    return true;
  }

  int SpeciesIndex(const std::string &name) {
    for (int s=0;s<PIC::nTotalSpecies;s++) if (name==PIC::MolecularData::GetChemSymbol(s)) return s;
    return -1;
  }

  std::string SpeciesList() {
    std::string res;

    for (int s=0;s<PIC::nTotalSpecies;s++) res+=std::string((s==0) ? "" : ", ")+PIC::MolecularData::GetChemSymbol(s);
    return res;
  }

  //the parameter tables of all models, saved and restored as one state (Apply() installs all or nothing)
  struct cState {
    int Selected[PIC::nTotalSpecies],NocillaMode[PIC::nTotalSpecies];
    double Temp,p[14][PIC::nTotalSpecies];

    static double* Table(int i) {
      static double *t[14]={Surface::MaxwellReflection::AccommodationCoefficient,
        Surface::CLL::NormalPartKineticEnergyAccommodationCoefficient,Surface::CLL::TangentialMomentumAccommodationCoefficient,
        Surface::QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient,
        Surface::QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient,
        Surface::QuasiSpecularReflection::Nocilla::DriftSpeed,Surface::QuasiSpecularReflection::Nocilla::DriftAngle,
        Surface::QuasiSpecularReflection::Nocilla::Temperature,
        Surface::QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient,
        Surface::QuasiSpecularReflection::Schamberg::Nu,Surface::QuasiSpecularReflection::Schamberg::BeamHalfAngle,
        Surface::QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation,
        Surface::QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation,
        Surface::QuasiSpecularReflection::Nocilla::EnergyAccommodation};
      return t[i];
    }

    void Save() {
      for (int s=0;s<PIC::nTotalSpecies;s++) {
        Selected[s]=Surface::Model::Selected[s],NocillaMode[s]=Surface::QuasiSpecularReflection::Nocilla::Mode[s];
        for (int i=0;i<14;i++) p[i][s]=Table(i)[s];
      }

      Temp=Surface::Temperature::Isothremal::Temp;
    }

    void Restore() const {
      for (int s=0;s<PIC::nTotalSpecies;s++) {
        Surface::Model::Selected[s]=Selected[s],Surface::QuasiSpecularReflection::Nocilla::Mode[s]=NocillaMode[s];
        for (int i=0;i<14;i++) Table(i)[s]=p[i][s];
      }

      Surface::Temperature::Isothremal::Temp=Temp;
    }
  };
}

bool Surface::RuntimeInput::Apply(const AMPS::Input::Block &block,std::string *error) {
  const std::string owner="the gas-surface interaction models (src/models/surface)";
  In::Status status;

  //child blocks: 'species: NAME' of the compiled species only
  status=In::RequireOnlyBlocks(block,{"species"},owner);
  if (!status.ok()) {if (error!=NULL) *error=status.message; return false;}

  for (const In::Block &child : block.blocks) {
    if (child.name.empty()) return Fail(child,"a species block requires the species name, e.g. '#begin species: H2O'",error);
    if (SpeciesIndex(child.name)<0) return Fail(child,"species '"+child.name+"' is not a species of this executable ("+SpeciesList()+")",error);
  }

  //all assignments with qualified keys; the dotted form 'species.H2O.KEY' and the block form are the same key
  std::vector<In::Assignment> items=In::Flatten(block);
  std::map<std::string,const In::Assignment*> keys;

  for (const In::Assignment &a : items) {
    if (keys.count(a.key)!=0) {
      const In::Assignment *first=keys[a.key];
      return Fail(a,"key '"+a.key+"' is given twice (also at "+first->where.file+":"+std::to_string(first->where.line)+")",error);
    }

    keys[a.key]=&a;
  }

  //the model
  if (keys.count("model")==0) return Fail(block,"block '"+block.Label()+"': the key 'model' is required (one of: diffuse, specular, maxwell, cll, quasi_specular_cll, nocilla, schamberg, adsorption)",error);

  const In::Assignment &am=*keys["model"];
  const int model=Model::Code(am.value);

  if (model<0) return Fail(am,"unknown gas-surface interaction model '"+am.value+"' (one of: diffuse, specular, maxwell, cll, quasi_specular_cll, nocilla, schamberg, adsorption)",error);

  //the wall temperature
  double Tw=Surface::Temperature::Isothremal::Temp;

  if (keys.count("wall_temperature_k")!=0) {
    const In::Assignment &a=*keys["wall_temperature_k"];

    if (!ParseNumber(a,&Tw,false,error)) return false;
    if (!(Tw>0.0)) return Fail(a,"the wall temperature must be > 0 K",error);
  }
  else if (UsesWallTemperature(model)) {
    return Fail(am,"model '"+am.value+"' requires the wall temperature 'wall_temperature_k'",error);
  }

  //the per-species parameters
  const std::vector<std::vector<std::string> > sets=ModelKeySets(model);
  std::vector<std::string> allKeys;
  std::map<std::string,double> value[PIC::nTotalSpecies];

  for (const std::vector<std::string> &set : sets) for (const std::string &p : set) allKeys.push_back(p);

  for (auto &k : keys) {
    const In::Assignment &a=*k.second;

    if ((a.key=="model")||(a.key=="wall_temperature_k")) continue;

    const std::string prefix="species.";
    if (a.key.compare(0,prefix.size(),prefix)!=0) return Fail(a,"key '"+a.key+"' is not a key of "+owner,error);

    std::string rest=a.key.substr(prefix.size());
    size_t dot=rest.find('.');
    if (dot==std::string::npos) return Fail(a,"key '"+a.key+"': expected 'species.SPECIES.PARAMETER'",error);

    const std::string sname=rest.substr(0,dot),param=rest.substr(dot+1);
    const int spec=SpeciesIndex(sname);

    if (spec<0) return Fail(a,"species '"+sname+"' is not a species of this executable ("+SpeciesList()+")",error);

    bool known=false;
    for (const std::string &p : allKeys) if (p==param) known=true;

    if (!known) {
      std::string list=KeyList(allKeys);
      return Fail(a,"parameter '"+param+"' is not a parameter of model '"+am.value+"' ("+(list.empty() ? std::string("it has none") : list)+")",error);
    }

    double x;
    if (!ParseNumber(a,&x,param=="nu",error)) return false;
    value[spec][param]=x;
  }

  //the parameter set of each species: with one set every key is required; with alternative sets exactly one complete
  //set and no key of another set
  int setOfSpecies[PIC::nTotalSpecies];

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    const char *sp=PIC::MolecularData::GetChemSymbol(s);

    if (sets.size()==1) {
      for (const std::string &p : sets[0]) if (value[s].count(p)==0) {
        return Fail(block,"model '"+am.value+"' requires 'species."+sp+"."+p+"' (or '"+p+"' in '#begin species: "+sp+"')",error);
      }

      setOfSpecies[s]=0;
      continue;
    }

    int nComplete=0,complete=-1,nUsed=0;

    for (size_t i=0;i<sets.size();i++) {
      int n=0;

      for (const std::string &p : sets[i]) if (value[s].count(p)!=0) n++;
      if (n>0) nUsed++;
      if (n==(int)sets[i].size()) nComplete++,complete=(int)i;
    }

    if ((nComplete!=1)||(nUsed!=1)) {
      std::string alternatives;

      for (size_t i=0;i<sets.size();i++) alternatives+=std::string((i==0) ? "" : ", or ")+"("+KeyList(sets[i])+")";
      return Fail(block,"model '"+am.value+"' requires for species "+sp+" exactly one complete set of parameters: "+alternatives,error);
    }

    setOfSpecies[s]=complete;
  }

  //the wall temperature of the Nocilla model with accommodation coefficients
  if ((model==Model::QuasiSpecularNocilla)&&(keys.count("wall_temperature_k")==0)) {
    for (int s=0;s<PIC::nTotalSpecies;s++) if (setOfSpecies[s]==1) {
      return Fail(am,"model '"+am.value+"' with accommodation coefficients requires the wall temperature 'wall_temperature_k'",error);
    }
  }

  //install everything, check it, and restore the previous state if any value is not valid
  cState saved;
  saved.Save();

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    std::map<std::string,double> &v=value[s];

    switch (model) {
    case Model::Maxwell:
      MaxwellReflection::AccommodationCoefficient[s]=v["accommodation"];
      break;
    case Model::CLL:
      CLL::NormalPartKineticEnergyAccommodationCoefficient[s]=v["normal_energy_accommodation"];
      CLL::TangentialMomentumAccommodationCoefficient[s]=v["tangential_momentum_accommodation"];
      break;
    case Model::QuasiSpecularCLL:
      QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[s]=v["normal_energy_accommodation"];
      QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[s]=v["tangential_momentum_accommodation"];
      break;
    case Model::QuasiSpecularNocilla:
      if (setOfSpecies[s]==0) {
        QuasiSpecularReflection::Nocilla::Mode[s]=QuasiSpecularReflection::Nocilla::Direct;
        QuasiSpecularReflection::Nocilla::DriftSpeed[s]=v["drift_speed_m_s"];
        QuasiSpecularReflection::Nocilla::DriftAngle[s]=v["drift_angle_deg"]*Pi/180.0;
        QuasiSpecularReflection::Nocilla::Temperature[s]=v["temperature_k"];
      }
      else {
        QuasiSpecularReflection::Nocilla::Mode[s]=QuasiSpecularReflection::Nocilla::AccommodationCoefficients;
        QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[s]=v["normal_momentum_accommodation"];
        QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[s]=v["tangential_momentum_accommodation"];
        QuasiSpecularReflection::Nocilla::EnergyAccommodation[s]=v["energy_accommodation"];
      }
      break;
    case Model::QuasiSpecularSchamberg:
      QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[s]=v["energy_accommodation"];
      QuasiSpecularReflection::Schamberg::Nu[s]=v["nu"];
      QuasiSpecularReflection::Schamberg::BeamHalfAngle[s]=v["beam_half_angle_deg"]*Pi/180.0;
      break;
    }

    std::string reason;

    if (CheckParameters(s,model,&reason)==false) {
      saved.Restore();

      //point at the species block, or at the model key when the parameters were given with dotted keys
      const In::Block *sb=block.FindBlock("species",PIC::MolecularData::GetChemSymbol(s));
      return (sb!=NULL) ? Fail(*sb,reason,error) : Fail(am,reason,error);
    }
  }

  Model::Set(-1,model);
  Surface::Temperature::Isothremal::Temp=Tw;
  return true;
}

std::string Surface::RuntimeInput::Describe() {
  std::string res;
  char buf[400];

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    const int model=Model::Selected[s];

    snprintf(buf,sizeof(buf),"species %s: model %s",PIC::MolecularData::GetChemSymbol(s),Model::Name(model));
    res+=buf;

    switch (model) {
    case Model::Maxwell:
      snprintf(buf,sizeof(buf),", accommodation=%.9g",MaxwellReflection::AccommodationCoefficient[s]);
      break;
    case Model::CLL:
      snprintf(buf,sizeof(buf),", alpha_n=%.9g, sigma_t=%.9g",CLL::NormalPartKineticEnergyAccommodationCoefficient[s],CLL::TangentialMomentumAccommodationCoefficient[s]);
      break;
    case Model::QuasiSpecularCLL:
      snprintf(buf,sizeof(buf),", alpha_n=%.9g, sigma_t=%.9g",QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[s],
          QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[s]);
      break;
    case Model::QuasiSpecularNocilla:
      if (QuasiSpecularReflection::Nocilla::Mode[s]==QuasiSpecularReflection::Nocilla::AccommodationCoefficients) {
        snprintf(buf,sizeof(buf),", accommodation coefficients: sigma'(normal momentum)=%.9g, sigma(tangential momentum)=%.9g, epsilon(energy)=%.9g",
            QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[s],QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[s],
            QuasiSpecularReflection::Nocilla::EnergyAccommodation[s]);
      }
      else {
        snprintf(buf,sizeof(buf),", U_r=%.9g m/s, theta_r=%.9g deg, T_r=%.9g K",QuasiSpecularReflection::Nocilla::DriftSpeed[s],
            QuasiSpecularReflection::Nocilla::DriftAngle[s]*180.0/Pi,QuasiSpecularReflection::Nocilla::Temperature[s]);
      }
      break;
    case Model::QuasiSpecularSchamberg:
      snprintf(buf,sizeof(buf),", alpha=%.9g, nu=%.9g, phi0=%.9g deg",QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[s],
          QuasiSpecularReflection::Schamberg::Nu[s],QuasiSpecularReflection::Schamberg::BeamHalfAngle[s]*180.0/Pi);
      break;
    default:
      buf[0]='\0';
    }

    res+=buf;
    snprintf(buf,sizeof(buf),", wall temperature=%.9g K\n",Surface::Temperature::Isothremal::Temp);
    res+=buf;
  }

  return res;
}

std::string Surface::RuntimeInput::ParameterLines() {
  std::string res;
  char buf[400];

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    const char *sp=PIC::MolecularData::GetChemSymbol(s);
    std::vector<std::pair<std::string,double> > p;

    switch (Model::Selected[s]) {
    case Model::Maxwell:
      p={{"accommodation",MaxwellReflection::AccommodationCoefficient[s]}};
      break;
    case Model::CLL:
      p={{"normal_energy_accommodation",CLL::NormalPartKineticEnergyAccommodationCoefficient[s]},
         {"tangential_momentum_accommodation",CLL::TangentialMomentumAccommodationCoefficient[s]}};
      break;
    case Model::QuasiSpecularCLL:
      p={{"normal_energy_accommodation",QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[s]},
         {"tangential_momentum_accommodation",QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[s]}};
      break;
    case Model::QuasiSpecularNocilla:
      if (QuasiSpecularReflection::Nocilla::Mode[s]==QuasiSpecularReflection::Nocilla::AccommodationCoefficients) {
        p={{"normal_momentum_accommodation",QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[s]},
           {"tangential_momentum_accommodation",QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[s]},
           {"energy_accommodation",QuasiSpecularReflection::Nocilla::EnergyAccommodation[s]}};
      }
      else {
        p={{"drift_speed_m_s",QuasiSpecularReflection::Nocilla::DriftSpeed[s]},
           {"drift_angle_deg",QuasiSpecularReflection::Nocilla::DriftAngle[s]*180.0/Pi},
           {"temperature_k",QuasiSpecularReflection::Nocilla::Temperature[s]}};
      }
      break;
    case Model::QuasiSpecularSchamberg:
      p={{"energy_accommodation",QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[s]},
         {"nu",QuasiSpecularReflection::Schamberg::Nu[s]},
         {"beam_half_angle_deg",QuasiSpecularReflection::Schamberg::BeamHalfAngle[s]*180.0/Pi}};
      break;
    }

    for (auto &x : p) {
      snprintf(buf,sizeof(buf),"species.%s.%s=%.9e\n",sp,x.first.c_str(),x.second);
      res+=buf;
    }
  }

  return res;
}
