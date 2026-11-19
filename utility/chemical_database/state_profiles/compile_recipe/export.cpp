// Thin source-backed export, not a second mass/charge/state resolver.
#include "srcMoon/species/PhysicalStateProfile.h"
#include "src/general/constants.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
namespace P=Moon::Species::PhysicalStateProfile;
int main(int argc,char**argv){
  if(argc<9)return 2;
  P::Files files;files.isotopes={argv[1],argv[2]};files.constants={argv[3],argv[4]};files.states={argv[5],argv[6]};
  int i=7;for(;i<argc&&std::string(argv[i])!="--states";i+=4){if(i+3>=argc)return 2;files.sources.push_back({argv[i],argv[i+1],argv[i+2],argv[i+3]});}
  if(i>=argc||++i>=argc)return 2;
  P::Registry profile;std::string error;if(!profile.Load(files,&error)){std::cerr<<error<<'\n';return 2;}
  P::Constant charge;if(!profile.GetConstant("elementary charge","C",&charge,&error)){std::cerr<<error<<'\n';return 2;}
  std::cout<<std::setprecision(std::numeric_limits<double>::max_digits10);
  const double unit=static_cast<double>(charge.value);
  if(!(unit>0)||!std::isfinite(unit))return 2;
  std::cout<<"CONSTANT\t"<<unit<<'\t'<<ElectronCharge<<'\t'<<charge.sourceIdentity<<'\t'<<charge.sourceRecord<<'\n';
  for(;i<argc;++i){Moon::Species::ReducedKinematics::ResolvedState state;Moon::Species::ReducedKinematics::MassValue physical;
    if(!profile.GetState(argv[i],&state,&error)||!profile.GetPhysicalRestMass(argv[i],&physical,&error)){std::cerr<<error<<'\n';return 77;}
    const double mass=static_cast<double>(state.referenceMassKg()),rest=static_cast<double>(physical.kg);
    const auto z=state.definition().elementaryCharge;
    if(!(mass>0)||!std::isfinite(mass)||!(rest>0)||!std::isfinite(rest))return 2;
    const double expected=static_cast<double>(z)*unit,token=expected/ElectronCharge;
    if(!std::isfinite(expected)||!std::isfinite(token)||token*ElectronCharge!=expected){std::cerr<<"source charge is not exactly emitted by native token expression\n";return 2;}
    std::cout<<"STATE\t"<<argv[i]<<'\t'<<mass<<'\t'<<rest<<'\t'<<z<<'\t'<<token<<'\t'<<expected<<'\n';
  }
}
