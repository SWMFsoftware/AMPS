// Thin adapter for the SAME production state-profile loader/resolver.
// It supplies explicit caller-owned CSV/source identities, never test masses
// or a duplicate state-energy/chemistry formula.
#include "srcMoon/species/PhysicalStateProfile.h"
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace P=Moon::Species::PhysicalStateProfile;
namespace K=Moon::Species::ReducedKinematics;
int Checks=0,Failures=0;
void Check(const std::string& name,bool pass) {
  ++Checks;if(!pass)++Failures;
  std::cout<<"CHECK\t"<<name<<'\t'<<(pass?"PASS":"FAIL")<<'\n';
}
int main(int argc,char** argv) {
  try {
    if(argc<7||(argc-7)%4!=0)throw std::runtime_error("expected three path/hash pairs and source id/path/hash/role groups");
    P::Files files;files.isotopes={argv[1],argv[2]};files.constants={argv[3],argv[4]};files.states={argv[5],argv[6]};
    for(int i=7;i<argc;i+=4)files.sources.push_back({argv[i],argv[i+1],argv[i+2],argv[i+3]});
    P::Registry profile;std::string error;
    if(!profile.Load(files,&error)){std::cout<<"LOAD-REFUSED\t"<<error<<'\n';return 2;}
    std::cout<<std::setprecision(std::numeric_limits<K::Real>::max_digits10);
    std::cout<<"EPSILON\t"<<std::numeric_limits<K::Real>::epsilon()<<'\n';
    std::cout<<"PROFILE\t"<<profile.identity()<<'\n';
    K::ResolvedState first;bool haveFirst=false;
    for(const auto& metadata:profile.states()) {
      K::ResolvedState state=first;
      if(metadata.available) {
        const bool ok=profile.GetState(metadata.id,&state,&error);
        Check("resolve-"+metadata.id,ok&&state.valid());
        if(!ok)continue;
        if(!haveFirst){first=state;haveFirst=true;}
        std::cout<<"STATE\t"<<metadata.id<<'\t'<<state.referenceMassKg()<<'\t'
            <<state.definition().elementaryCharge<<'\t'<<state.definition().formation.joules
            <<'\t'<<state.definition().excitation.joules<<'\n';
        K::MassValue physical;
        const bool mass=profile.GetPhysicalRestMass(metadata.id,&physical,&error);
        Check("physical-mass-"+metadata.id,mass==metadata.physicalRestMassAvailable);
        if(mass)std::cout<<"PHYSICAL\t"<<metadata.id<<'\t'<<physical.kg<<'\n';
      } else {
        const auto old=state;
        const bool ok=profile.GetState(metadata.id,&state,&error);
        Check("unavailable-"+metadata.id,!ok&&!error.empty()&&
            state.definition().identity==old.definition().identity&&
            state.referenceMassKg()==old.referenceMassKg());
        std::cout<<"UNAVAILABLE\t"<<metadata.id<<'\t'<<metadata.missingReason<<'\n';
        K::MassValue physical{-7,"SENTINEL"};
        Check("no-reference-mass-substitution-"+metadata.id,
            !profile.GetPhysicalRestMass(metadata.id,&physical,&error)&&physical.kg==-7&&physical.sourceIdentity=="SENTINEL");
      }
    }
    for(const auto& request:std::vector<std::pair<std::string,std::string>>{
        {"atomic mass constant","kg"},{"electron mass","kg"},
        {"elementary charge","C"},{"speed of light in vacuum","m s^-1"}}) {
      P::Constant value;const bool ok=profile.GetConstant(request.first,request.second,&value,&error);
      Check("constant-"+request.first,ok);
      if(ok)std::cout<<"CONSTANT\t"<<request.first<<'\t'<<value.value<<'\t'<<value.unit<<'\n';
    }
    P::Constant unchanged;unchanged.value=-7;unchanged.unit="SENTINEL";
    Check("wrong-unit-preserves-value",!profile.GetConstant("electron mass","eV",&unchanged,&error)&&
        unchanged.value==-7&&unchanged.unit=="SENTINEL");
    const auto identity=profile.identity();auto bad=files;bad.states.sha256=std::string(64,'0');
    Check("bad-profile-load-preserves-registry",!profile.Load(bad,&error)&&profile.identity()==identity&&profile.loaded());
    bad=files;if(!bad.sources.empty())bad.sources.front().sha256=std::string(64,'0');
    Check("bad-native-source-preserves-registry",!profile.Load(bad,&error)&&profile.identity()==identity&&profile.loaded());
    K::ResolvedState untouched=first;
    Check("unknown-state-preserves-output",!profile.GetState("unknown:state",&untouched,&error)&&
        untouched.definition().identity==first.definition().identity&&untouched.referenceMassKg()==first.referenceMassKg());
    std::cout<<"SUMMARY\t"<<Checks<<'\t'<<Failures<<'\n';return Failures?1:0;
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 2;}
}
