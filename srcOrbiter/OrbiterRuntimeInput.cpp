//$Id$
//The runtime input of the Orbiter application (OrbiterRuntimeInput.h)

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

#include "pic.h"
#include "Surface.h"
#include "OrbiterRuntimeInput.h"
#include "src/general/input/amps_input.h"

bool Orbiter::RuntimeInput::Loaded=false;
std::string Orbiter::RuntimeInput::LoadedPath,Orbiter::RuntimeInput::SelectedCase;
int Orbiter::RuntimeInput::nTotalIterations=0;

namespace {
  namespace In=AMPS::Input;

  struct cValues {
    double Velocity[3],Temperature,NumberDensity[PIC::nTotalSpecies],CharacteristicSpeed;
    long int nIterations;
  };

  bool Fail(const In::Assignment &a,const std::string &reason,std::string *error) {
    if (error!=NULL) *error=In::ErrorAt(a.where,reason,a.text).message;
    return false;
  }

  bool ParseNumber(const In::Assignment &a,const std::string &text,double *x,std::string *error) {
    const char *s=text.c_str();
    char *end=NULL;

    errno=0;
    *x=strtod(s,&end);
    while ((end!=NULL)&&((*end==' ')||(*end=='\t'))) end++;

    if ((end==s)||(end==NULL)||(*end!='\0')||(errno==ERANGE)||!std::isfinite(*x)) {
      return Fail(a,"key '"+a.key+"': '"+text+"' is not a finite number",error);
    }

    return true;
  }

  int SpeciesIndex(const std::string &name) {
    for (int s=0;s<PIC::nTotalSpecies;s++) if (name==PIC::MolecularData::GetChemSymbol(s)) return s;
    return -1;
  }

  //the keys of the application block (without the gas-surface interaction block)
  bool ParseApplication(const In::Block &app,cValues *v,const In::Block **gsi,std::string *error) {
    In::Status status=In::RequireOnlyBlocks(app,{"upstream","time_step","run","physics"},"the Orbiter application (srcOrbiter)");
    if (!status.ok()) {if (error!=NULL) *error=status.message; return false;}

    *gsi=NULL;
    In::Block rest=app;
    rest.blocks.clear();

    for (const In::Block &child : app.blocks) {
      if (child.kind=="physics") {
        if (child.name!="gas_surface_interaction") {
          if (error!=NULL) *error=In::ErrorAt(child.begin,"block '"+child.Label()+"' is not handled by the Orbiter application (allowed: physics: gas_surface_interaction)").message;
          return false;
        }

        *gsi=&child;
      }
      else rest.blocks.push_back(child);
    }

    if (*gsi==NULL) {
      if (error!=NULL) *error=In::ErrorAt(app.begin,"block '"+app.Label()+"' requires the block 'physics: gas_surface_interaction'").message;
      return false;
    }

    //the keys, with the dotted names of the nested blocks
    std::vector<In::Assignment> items=In::Flatten(rest);
    std::map<std::string,const In::Assignment*> keys;

    for (const In::Assignment &a : items) {
      if (keys.count(a.key)!=0) return Fail(a,"key '"+a.key+"' is given twice",error);
      keys[a.key]=&a;
    }

    std::vector<std::string> required={"schema","upstream.velocity_m_s","upstream.temperature_k","time_step.characteristic_speed_m_s","run.iterations"};
    for (int s=0;s<PIC::nTotalSpecies;s++) required.push_back(std::string("upstream.species.")+PIC::MolecularData::GetChemSymbol(s)+".number_density_m-3");

    for (const std::string &k : required) if (keys.count(k)==0) {
      if (error!=NULL) *error=In::ErrorAt(app.begin,"block '"+app.Label()+"' requires the key '"+k+"'").message;
      return false;
    }

    for (auto &k : keys) {
      bool known=false;

      for (const std::string &r : required) if (r==k.first) known=true;
      if (known) continue;

      const std::string prefix="upstream.species.";
      if (k.first.compare(0,prefix.size(),prefix)==0) {
        std::string sname=k.first.substr(prefix.size());
        sname=sname.substr(0,sname.find('.'));
        if (SpeciesIndex(sname)<0) return Fail(*k.second,"species '"+sname+"' is not a species of this executable",error);
      }

      return Fail(*k.second,"key '"+k.first+"' is not a key of the Orbiter application",error);
    }

    //the values
    double x;
    const In::Assignment *a;

    a=keys["schema"];
    if ((ParseNumber(*a,a->value,&x,error)==false)||(x!=1.0)) return Fail(*a,"only schema = 1 is supported",error);

    a=keys["upstream.velocity_m_s"];
    {
      std::string text=a->value;
      for (char &c : text) if (c==',') c=' ';

      const char *s=text.c_str();
      char *end=NULL;
      int n=0;

      for (n=0;n<3;n++) {
        errno=0;
        v->Velocity[n]=strtod(s,&end);
        if ((end==s)||(errno==ERANGE)||!std::isfinite(v->Velocity[n])) break;
        s=end;
      }

      while ((*s==' ')||(*s=='\t')) s++;
      if ((n!=3)||(*s!='\0')) return Fail(*a,"key '"+a->key+"': three components (m/s) are required, e.g. '-3000, 0, 0'",error);
      if (v->Velocity[0]*v->Velocity[0]+v->Velocity[1]*v->Velocity[1]+v->Velocity[2]*v->Velocity[2]<=0.0) return Fail(*a,"the upstream velocity must not be zero",error);
    }

    a=keys["upstream.temperature_k"];
    if (ParseNumber(*a,a->value,&v->Temperature,error)==false) return false;
    if (!(v->Temperature>0.0)) return Fail(*a,"the upstream temperature must be > 0 K",error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      a=keys[std::string("upstream.species.")+PIC::MolecularData::GetChemSymbol(s)+".number_density_m-3"];
      if (ParseNumber(*a,a->value,v->NumberDensity+s,error)==false) return false;
      if (!(v->NumberDensity[s]>0.0)) return Fail(*a,"the number density must be > 0",error);
    }

    a=keys["time_step.characteristic_speed_m_s"];
    if (ParseNumber(*a,a->value,&v->CharacteristicSpeed,error)==false) return false;
    if (!(v->CharacteristicSpeed>0.0)) return Fail(*a,"the characteristic speed must be > 0",error);

    a=keys["run.iterations"];
    if (ParseNumber(*a,a->value,&x,error)==false) return false;
    if ((x<1.0)||(x>1.0E9)||(x!=floor(x))) return Fail(*a,"the number of iterations must be an integer >= 1",error);
    v->nIterations=(long int)x;

    return true;
  }

  std::uint64_t HashCanonical(const std::string &text) {
    std::uint64_t hash=UINT64_C(1469598103934665603);

    for (unsigned char byte : text) {
      hash^=byte;
      hash*=UINT64_C(1099511628211);
    }

    return hash;
  }
}

bool Orbiter::RuntimeInput::LoadAndApply(const std::string &path,const std::string &caseName,std::string *error) {
  In::Document document;
  In::Block resolved;
  const In::Block *gsi=NULL;
  cValues values;
  std::string localError;

  //parse, select the case, and check the application keys (every process)
  In::Status status=In::ParseFile(path,&document);

  if (status.ok()) status=In::RequireNoAssignments(document.root);
  if (status.ok()) status=In::RequireOnlyBlocks(document.root,{"application"},"the AMPS runtime-input root");

  if (status.ok()) {
    if ((document.root.blocks.size()!=1)||(document.root.blocks.front().name!="orbiter")) {
      status=In::Status::Error(path+": the Orbiter runtime input requires exactly one block '#begin application: orbiter'");
    }
  }

  if (status.ok()) status=In::ResolveCase(document.root.blocks.front(),caseName,&resolved);

  bool ok=status.ok();
  if (ok==false) localError=status.message;
  else ok=ParseApplication(resolved,&values,&gsi,&localError);

  int localOK=ok ? 1 : 0,globalOK=0;
  MPI_Allreduce(&localOK,&globalOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);

  if (globalOK==0) {
    if (error!=NULL) *error=(localOK==1) ? std::string("another MPI process could not read the runtime input") : localError;
    return false;
  }

  //every process must have resolved the same input
  unsigned long long localHash=HashCanonical(In::ToCanonicalText(resolved)),minHash=0,maxHash=0;

  MPI_Allreduce(&localHash,&minHash,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localHash,&maxHash,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  if (minHash!=maxHash) {
    if (error!=NULL) *error="the MPI processes resolved different runtime inputs";
    return false;
  }

  //the gas-surface interaction models (installed only when the whole block is valid)
  ok=Surface::RuntimeInput::Apply(*gsi,&localError);
  localOK=ok ? 1 : 0;
  MPI_Allreduce(&localOK,&globalOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);

  if (globalOK==0) {
    if (error!=NULL) *error=(localOK==1) ? std::string("another MPI process rejected the runtime input") : localError;
    return false;
  }

  //the adsorption model removes the particles: the Orbiter must be compiled with an adsorption model
  if (_ORBITER_ADSORPTION_MODEL__MODE_ == _ORBITER_ADSORPTION_MODEL__MODE__OFF_) {
    for (int s=0;s<PIC::nTotalSpecies;s++) if (Surface::Model::Selected[s]==Surface::Model::Adsorption) {
      if (error!=NULL) *error="the adsorption model is selected, but the executable is compiled without adsorption (_ORBITER_ADSORPTION_MODEL__MODE_)";
      return false;
    }
  }

  //install the application values
  for (int idim=0;idim<3;idim++) Orbiter::UpstreamBC::Velocity[idim]=values.Velocity[idim];
  Orbiter::UpstreamBC::Temperature=values.Temperature;
  for (int s=0;s<PIC::nTotalSpecies;s++) Orbiter::UpstreamBC::NumberDensity[s]=values.NumberDensity[s];
  Orbiter::CharacteristicSpeed=values.CharacteristicSpeed;
  nTotalIterations=(int)values.nIterations;

  Loaded=true;
  LoadedPath=path;
  SelectedCase=caseName;
  return true;
}

std::string Orbiter::RuntimeInput::Describe() {
  std::string res;
  char buf[400];

  snprintf(buf,sizeof(buf),"input=%s\ncase=%s\n",LoadedPath.c_str(),SelectedCase.c_str());
  res+=buf;
  snprintf(buf,sizeof(buf),"upstream_velocity_m_s=%.9e %.9e %.9e\nupstream_temperature_k=%.9e\n",Orbiter::UpstreamBC::Velocity[0],
      Orbiter::UpstreamBC::Velocity[1],Orbiter::UpstreamBC::Velocity[2],Orbiter::UpstreamBC::Temperature);
  res+=buf;

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    snprintf(buf,sizeof(buf),"upstream_number_density_m-3[%s]=%.9e\nmass_kg[%s]=%.9e\n",PIC::MolecularData::GetChemSymbol(s),
        Orbiter::UpstreamBC::NumberDensity[s],PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetMass(s));
    res+=buf;
  }

  snprintf(buf,sizeof(buf),"characteristic_speed_m_s=%.9e\niterations=%i\nboltzmann_constant=%.9e\n",Orbiter::CharacteristicSpeed,nTotalIterations,Kbol);
  res+=buf;
  return res;
}
