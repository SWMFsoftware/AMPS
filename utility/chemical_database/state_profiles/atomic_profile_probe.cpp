// Extend the maintained profile adapter with the actual Na ground-ion use.
// Original parser/unavailable/rollback checks execute first and are retained.
#define main OriginalProfileProbeMain
#include "probe.cpp"
#undef main

int main(int argc,char** argv) {
  const int original=OriginalProfileProbeMain(argc,argv);
  if(original!=0)return original;
  P::Files files;files.isotopes={argv[1],argv[2]};files.constants={argv[3],argv[4]};files.states={argv[5],argv[6]};
  for(int i=7;i<argc;i+=4)files.sources.push_back({argv[i],argv[i+1],argv[i+2],argv[i+3]});
  P::Registry registry;std::string error;
  if(!registry.Load(files,&error)){std::cerr<<error<<'\n';return 2;}
  K::ResolvedState neutral,ion,electron;
  const bool states=registry.GetState("23Na:ground",&neutral,&error)&&
    registry.GetState("23Na+:ground",&ion,&error)&&registry.GetState("electron:free",&electron,&error);
  Check("atomic-ground-states-resolve",states);
  if(!states)return 1;
  const auto& n=neutral.definition();const auto& i=ion.definition();const auto& e=electron.definition();
  Check("atomic-explicit-isotope-nuclei",n.nuclei.size()==1&&i.nuclei.size()==1&&e.nuclei.empty()&&
    n.nuclei[0].identity=="23Na"&&i.nuclei[0].identity=="23Na"&&n.nuclei[0].count==1&&i.nuclei[0].count==1);
  Check("atomic-integer-charges",n.elementaryCharge==0&&i.elementaryCharge==1&&e.elementaryCharge==-1);
  Check("atomic-ground-excitation-zeros",n.excitation.joules==0&&i.excitation.joules==0&&e.excitation.joules==0);
  K::FormationBalance balance;
  const bool q=K::FormationQ({neutral},{ion,electron},&balance,&error);
  Check("atomic-production-formation-Q",q&&balance.qJ==-i.formation.joules&&i.formation.joules>0);
  K::MassValue physical;
  const bool mass=registry.GetPhysicalRestMass("23Na+:ground",&physical,&error);
  Check("atomic-physical-mass-distinct-reference",mass&&physical.kg>ion.referenceMassKg());
  std::cout<<"ATOMIC\tformation_source\t"<<i.formation.sourceIdentity<<'\n'
    <<"ATOMIC\tformation_J\t"<<i.formation.joules<<'\n'
    <<"ATOMIC\tQ_J\t"<<balance.qJ<<'\n'
    <<"ATOMIC\treference_mass_kg\t"<<ion.referenceMassKg()<<'\n'
    <<"ATOMIC\tphysical_mass_kg\t"<<physical.kg<<'\n'
    <<"ATOMIC-SUMMARY\t"<<Checks<<'\t'<<Failures<<'\n';
  return Failures?1:0;
}
