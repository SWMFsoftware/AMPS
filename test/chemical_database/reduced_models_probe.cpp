// Thin independent-of-PIC adapter: every evaluated model is the shipped shared
// header. All numbers read below are declared analytical fixtures, never data.
#include "../../src/models/chemical_database/ReducedModels.h"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace C=AMPS::ChemicalDatabase;
namespace M=C::ReducedModels;
int main(int argc,char** argv) {
  if(argc!=3) return 2;
  std::ofstream checks(argv[1]),metrics(argv[2]);
  if(!checks||!metrics) return 2;
  checks<<"name,passed\n"; metrics<<"case,key,value\n";
  metrics<<std::setprecision(std::numeric_limits<long double>::max_digits10);
  int failed=0,count=0;
  auto check=[&](const char* name,bool pass){checks<<name<<','<<(pass?1:0)<<'\n';++count;if(!pass)++failed;};
  auto emit=[&](const std::string& name,const std::string& key,long double value){metrics<<name<<','<<key<<','<<value<<'\n';};
  const long double nan=std::numeric_limits<long double>::quiet_NaN();
  const long double inf=std::numeric_limits<long double>::infinity();
  M::IonMass ion{7,8,9};
  check("ion_null_refused",!M::IonPhysicalMass(4,1,1,2,nullptr).ok());
  check("ion_negative_energy_transactional",!M::IonPhysicalMass(4,1,-1,2,&ion).ok()&&ion.referenceMassKg==7&&ion.bindingMassKg==8&&ion.physicalMassKg==9);
  check("ion_no_positive_reference_mass",!M::IonPhysicalMass(1,1,1,2,&ion).ok());
  check("ion_nan_and_infinite_refused",!M::IonPhysicalMass(nan,1,1,2,&ion).ok()&&!M::IonPhysicalMass(4,1,1,inf,&ion).ok());
  M::RadiativeSet radiative{7,8,{9}};
  check("radiative_null_refused",!M::RadiativeBranches({1},nullptr).ok());
  check("radiative_zero_total_outside_domain",M::RadiativeBranches({0,0},&radiative).code==C::Code::OutsideDomain&&radiative.totalRateS==7&&radiative.branches==std::vector<long double>{9});
  check("radiative_negative_transactional",!M::RadiativeBranches({1,-1},&radiative).ok()&&radiative.lifetimeS==8);
  check("radiative_empty_nan_inf_refused",!M::RadiativeBranches({},&radiative).ok()&&!M::RadiativeBranches({nan},&radiative).ok()&&!M::RadiativeBranches({inf},&radiative).ok());
  check("radiative_weak_positive_underflow_not_deleted",M::RadiativeBranches({std::numeric_limits<long double>::denorm_min(),std::numeric_limits<long double>::max()},&radiative).code==C::Code::NumericalFailure);
  long double secondary=7;
  check("opal_null_refused",!M::OpalSecondaryEnergy(10,2,1,.5L,nullptr).ok());
  check("opal_threshold_transactional",!M::OpalSecondaryEnergy(2,2,1,.5L,&secondary).ok()&&secondary==7);
  check("opal_width_and_quantile_refused",!M::OpalSecondaryEnergy(10,2,0,.5L,&secondary).ok()&&!M::OpalSecondaryEnergy(10,2,1,1.1L,&secondary).ok());
  check("opal_nan_refused",!M::OpalSecondaryEnergy(10,2,1,nan,&secondary).ok());
  check("opal_zero_endpoint",M::OpalSecondaryEnergy(10,2,1,0,&secondary).ok()&&secondary==0);
  check("opal_one_endpoint",M::OpalSecondaryEnergy(10,2,1,1,&secondary).ok()&&secondary==4);
  const M::ConservationBounds b{0,1e-12L,0,1e-12L,1e-12L};
  M::TwoBodyState two{};two.totalKineticEnergyJ=91;
  check("two_body_null_refused",!M::TwoBodyProducts({0,0,0},1,0,1,1,{1,0,0},b,nullptr).ok());
  check("two_body_zero_direction_transactional",!M::TwoBodyProducts({0,0,0},1,0,1,1,{0,0,0},b,&two).ok()&&two.totalKineticEnergyJ==91);
  check("two_body_negative_COM_no_clamp",M::TwoBodyProducts({0,0,0},0,-1e-30L,1,1,{1,0,0},b,&two).code==C::Code::OutsideDomain&&two.totalKineticEnergyJ==91);
  check("two_body_invalid_masses",!M::TwoBodyProducts({0,0,0},1,0,0,1,{1,0,0},b,&two).ok());
  M::ConservationBounds invalid=b;invalid.unitVectorTolerance=1;
  check("two_body_invalid_bounds",!M::TwoBodyProducts({0,0,0},1,0,1,1,{1,0,0},invalid,&two).ok());
  check("two_body_nonfinite_refused",!M::TwoBodyProducts({nan,0,0},1,0,1,1,{1,0,0},b,&two).ok());
  const M::ConservationBounds exact{0,0,0,0,0};
  check("two_body_exact_zero_COM",M::TwoBodyProducts({1,0,0},.25L,0,1,1,{1,0,0},exact,&two).ok()&&two.centreOfMassKineticEnergyJ==0&&two.momentumKgMS[0][0]==.5L);
  const long double hugeMass=std::ldexp(1.0L,std::numeric_limits<long double>::max_exponent/2);
  const long double tinyMass=1/hugeMass;
  two.totalKineticEnergyJ=91;
  check("two_body_lost_positive_mass_refused",M::TwoBodyProducts({hugeMass,0,0},hugeMass/2,0,tinyMass,hugeMass,{0,1,0},exact,&two).code==C::Code::NumericalFailure&&two.totalKineticEnergyJ==91);
  const long double smallMomentum=std::ldexp(1.0L,(std::numeric_limits<long double>::min_exponent-std::numeric_limits<long double>::digits)/2-10);
  check("two_body_underflowed_bulk_not_certified_zero",M::TwoBodyProducts({smallMomentum,0,0},0,0,1,1,{1,0,0},exact,&two).code==C::Code::NumericalFailure&&two.totalKineticEnergyJ==91);
  M::InheritedState inherited{};inherited.signedReservoirEnergyJ=99;inherited.signedReservoirChargeC=98;
  check("inherit_null_refused",!M::InheritParentVelocity(2,1,0,1,0,{1,0,0},{0,0,0},0,0,nullptr).ok());
  check("inherit_invalid_mass_transactional",!M::InheritParentVelocity(2,0,0,1,0,{1,0,0},{0,0,0},0,0,&inherited).ok()&&inherited.signedReservoirEnergyJ==99&&inherited.signedReservoirChargeC==98);
  check("inherit_nonfinite_charge_refused",!M::InheritParentVelocity(2,1,nan,1,0,{1,0,0},{0,0,0},0,0,&inherited).ok());
  check("inherit_signed_negative_reservoir",M::InheritParentVelocity(1,2,0,1,-1,{1,0,0},{0,0,0},0,0,&inherited).ok()&&inherited.signedReservoirEnergyJ==-.5L&&inherited.signedReservoirChargeC==-2);
  check("inherit_small_Q_survives_bulk_cancellation",M::InheritParentVelocity(5,5,0,0,0,{1e10L,0,0},{0,0,0},0,1,&inherited).ok()&&inherited.signedReservoirEnergyJ==1);
  M::ThermalState thermal{};thermal.remainingPopulation=99;
  check("thermal_null_refused",!M::ThermalRelease(1,1,1,1,1,1,nullptr).ok());
  check("thermal_zero_temperature_transactional",!M::ThermalRelease(1,1,1,0,1,1,&thermal).ok()&&thermal.remainingPopulation==99);
  check("thermal_negative_stock_refused",!M::ThermalRelease(1,1,1,1,-1,1,&thermal).ok());
  check("thermal_nan_inf_refused",!M::ThermalRelease(nan,1,1,1,1,1,&thermal).ok()&&!M::ThermalRelease(1,1,1,1,1,inf,&thermal).ok());
  check("thermal_zero_attempt_explicit",M::ThermalRelease(0,1,1,1,10,1,&thermal).ok()&&thermal.remainingPopulation==10&&thermal.releasedPopulation==0);
  check("thermal_zero_duration",M::ThermalRelease(1,0,1,1,10,0,&thermal).ok()&&thermal.remainingPopulation==10&&thermal.rateS==1);
  const long double minimum=std::numeric_limits<long double>::min();
  check("thermal_rate_underflow_exposure_preserved",M::ThermalRelease(minimum,1000,1,1,std::numeric_limits<long double>::max()/4,1/minimum,&thermal).ok()&&thermal.rateUnderflowed&&thermal.releasedPopulation>0);
  const long double barrier=-std::log(std::numeric_limits<long double>::denorm_min())+10;
  check("thermal_prefactor_rescues_rate",M::ThermalRelease(std::numeric_limits<long double>::max()/4,barrier,1,1,1,1,&thermal).ok()&&thermal.rateS>0&&!thermal.rateUnderflowed);
  check("thermal_large_exposure_explicit_overflow_diagnostic",M::ThermalRelease(std::numeric_limits<long double>::max(),0,1,1,2,std::numeric_limits<long double>::max(),&thermal).ok()&&!thermal.integratedHazard&&thermal.remainingPopulation==0&&thermal.releasedPopulation==2);

  std::string operation,name;
  while(std::cin>>operation>>name) {
    C::Status status;
    if(operation=="ION") {
      long double n,e,q,c;std::cin>>n>>e>>q>>c;M::IonMass v{};
      status=M::IonPhysicalMass(n,e,q,c,&v);
      if(status.ok()) {emit(name,"reference",v.referenceMassKg);emit(name,"binding",v.bindingMassKg);emit(name,"physical",v.physicalMassKg);}
    } else if(operation=="RAD") {
      std::size_t n;std::cin>>n;std::vector<long double> rates(n);
      for(auto& a:rates) { std::cin>>a; }
      M::RadiativeSet v{};status=M::RadiativeBranches(rates,&v);
      if(status.ok()) {emit(name,"total",v.totalRateS);emit(name,"lifetime",v.lifetimeS);for(std::size_t i=0;i<n;++i)emit(name,"branch"+std::to_string(i),v.branches[i]);}
    } else if(operation=="OPAL") {
      long double e,i,w,u,v;std::cin>>e>>i>>w>>u;status=M::OpalSecondaryEnergy(e,i,w,u,&v);
      if(status.ok()) { emit(name,"energy",v);emit(name,"support_native",(e-i)/2); }
    } else if(operation=="TWO") {
      long double m1,m2,k,q;M::Vector3 p{},d{};std::cin>>m1>>m2>>k>>q;for(auto& x:p)std::cin>>x;for(auto& x:d)std::cin>>x;
      M::TwoBodyState v{};status=M::TwoBodyProducts(p,k,q,m1,m2,d,b,&v);
      if(status.ok()) {emit(name,"com",v.centreOfMassKineticEnergyJ);emit(name,"total",v.totalKineticEnergyJ);
        for(int a=0;a<2;++a){emit(name,"k"+std::to_string(a),v.kineticEnergyJ[a]);for(int j=0;j<3;++j){emit(name,"p"+std::to_string(a)+std::to_string(j),v.momentumKgMS[a][j]);emit(name,"v"+std::to_string(a)+std::to_string(j),v.velocityMS[a][j]);}}}
    } else if(operation=="INHERIT") {
      long double m1,m2,c1,c2,cb,e,q;M::Vector3 v{},p{};std::cin>>m1>>m2>>c1>>c2>>cb;for(auto& x:v)std::cin>>x;for(auto& x:p)std::cin>>x;std::cin>>e>>q;
      M::InheritedState a{};status=M::InheritParentVelocity(m1,m2,c1,c2,cb,v,p,e,q,&a);
      if(status.ok()){emit(name,"heavy_k",a.heavyKineticEnergyJ);emit(name,"bath_e",a.signedReservoirEnergyJ);emit(name,"bath_c",a.signedReservoirChargeC);for(int j=0;j<3;++j){emit(name,"p"+std::to_string(j),a.heavyMomentumKgMS[j]);emit(name,"v"+std::to_string(j),a.heavyVelocityMS[j]);emit(name,"bath_p"+std::to_string(j),a.signedReservoirMomentumKgMS[j]);}}
    } else if(operation=="THERM") {
      long double nu,e,k,t,n,dt;std::cin>>nu>>e>>k>>t>>n>>dt;M::ThermalState v{};status=M::ThermalRelease(nu,e,k,t,n,dt,&v);
      if(status.ok()){emit(name,"rate",v.rateS);emit(name,"remaining",v.remainingPopulation);emit(name,"released",v.releasedPopulation);if(v.integratedHazard)emit(name,"hazard",*v.integratedHazard);}
    } else {std::cerr<<"unknown fixture operation\n";return 2;}
    if(!std::cin){std::cerr<<"malformed analytical fixture input\n";return 2;}
    emit(name,"status",static_cast<int>(status.code));
    if(!status.ok())std::cerr<<name<<": "<<status.message<<'\n';
  }
  std::cout<<"SUMMARY checks="<<count<<" failed="<<failed<<'\n';
  return failed?1:0;
}
