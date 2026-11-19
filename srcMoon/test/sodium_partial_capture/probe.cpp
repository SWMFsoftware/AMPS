// P07 production-kernel controls. Native data and NIST values arrive from the
// maintained Python evidence reader; no test-only capture/recoil law is used.
#include "srcMoon/species/SodiumPartialCapture.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
namespace S=Moon::Species::SodiumPartialCapture;
namespace K=Moon::Species::ReducedKinematics;
namespace {int checks=0,failures=0;
void Check(bool value,const char* name){++checks;if(!value)++failures;std::cout<<(value?"PASS ":"FAIL ")<<name<<'\n';}
void Require(bool value,const std::string& error){if(!value)throw std::runtime_error(error);}
long double Number(const char* text){std::istringstream s(text);long double x=0;if(!(s>>x))throw std::runtime_error("bad numeric argument");return x;}
K::ResolvedState State(const char* identity,const char* state,const char* nucleus,int charge,long double formation,
    const std::vector<K::IsotopeMass>& inventory,const K::MassValue& electron){
  K::StateDefinition d;d.identity=identity;d.state=state;d.nuclei.push_back({nucleus,1});d.elementaryCharge=charge;
  d.formation={formation,"native-NIST-ground-energies-plus-separated-atom-reference-zero"};
  d.excitation={0,"explicit-ground-state-reference-definition"};
  K::ResolvedState resolved;std::string error;Require(K::ResolveState(d,inventory,electron,&resolved,&error),error);return resolved;
}
}
int main(int argc,char** argv){
  if(argc!=10)return 2;
  try{
    const long double hMass=Number(argv[2]),naMass=Number(argv[3]),eMass=Number(argv[4]),eV=Number(argv[5]),u=Number(argv[6]);
    const long double hIon=Number(argv[7])*eV,naIon=Number(argv[8])*eV;
    std::vector<K::IsotopeMass> masses={{{"1H"},{hMass,"native-NIST-isotope-CODATA"}},{{"23Na"},{naMass,"native-NIST-isotope-CODATA"}}};
    K::MassValue electron{eMass,"native-CODATA2022-electron"};
    const auto proton=State("H+","bare-proton","1H",1,hIon,masses,electron);
    const auto sodium=State("Na","3s","23Na",0,0,masses,electron);
    const auto hydrogen=State("H","1s","1H",0,0,masses,electron);
    const auto ion=State("Na+","2p6","23Na",1,naIon,masses,electron);
    S::Constants constants{u,eV,"native-CODATA2022"};
    S::Registry registry;std::string error;
    Check(registry.LoadNativeNodes(argv[1],&error),"load_exact_published_partial_source");
    Require(registry.loaded(),error);Check(registry.nodes().size()==12,"all_twelve_native_nodes_retained");
    Check(!registry.physicalTargetStateQualified(),"initial_target_state_remains_unqualified");
    std::ofstream output(argv[9]);output<<std::setprecision(24)
      <<"energy_eV_per_u,cross_section_m2,relative_speed_m_s,rate_s-1,Q_J,px_in,px_out,kinetic_in_J,kinetic_out_J\n";
    for(const auto& node:registry.nodes()){
      S::Rate rate;Check(registry.RateAtNativeNode(node.energyEvPerU,2.5e6L,constants,&rate,&error),"native_node_rate");
      S::ForwardInput f;f.incomingProton=proton;f.incomingSodium=sodium;f.outgoingHydrogen=hydrogen;f.outgoingSodiumIon=ion;
      f.relativeIncidentDirection={{1,0,0}};f.sourceRate=rate;f.constants=constants;
      f.context={"source-target-rest-inertial","IAEA-D119990-forward-closure-approved-P07"};
      const long double eps=std::numeric_limits<long double>::epsilon();
      f.controls.maximumClassicalSpeedMS=3.0e6L; // mathematical native-source fixture domain, not a lunar validity claim
      f.controls.unitVectorTolerance=512*eps;
      f.controls.referenceMassKg={0,512*eps};f.controls.momentumKgMS={1e-40L,2048*eps};
      f.controls.energyJ={1e-32L,4096*eps};f.controls.thresholdJ={1e-32L,4096*eps};
      K::FinalState result;Check(S::ForwardProducts(registry,f,&result,&error),"source_bound_forward_products");
      Require(result.products.size()==2,error);
      const long double beforeP=proton.referenceMassKg()*rate.relativeSpeedMS();
      const long double beforeK=.5L*proton.referenceMassKg()*rate.relativeSpeedMS()*rate.relativeSpeedMS();
      long double afterP=0,afterK=0;
      for(const auto& product:result.products){afterP+=product.momentumKgMS[0];for(int k=0;k<3;++k)afterK+=.5L*product.state.referenceMassKg()*product.velocityMS[k]*product.velocityMS[k];}
      Check(result.modelIdentity=="EXPERIMENTAL-D119990-ground-partial-forward-COM","partial_model_label_preserved");
      output<<node.energyEvPerU<<','<<node.crossSectionM2<<','<<rate.relativeSpeedMS()<<','<<rate.rateS()<<','<<result.stateQJ<<','<<beforeP<<','<<afterP<<','<<beforeK<<','<<afterK<<'\n';
      K::FinalState preserved=result;
      if(node.energyEvPerU==4000) {
        const auto original=f;
        f.context.frameIdentity="LSO";
        Check(!S::ForwardProducts(registry,f,&result,&error)&&result.modelIdentity==preserved.modelIdentity,"untransformed_frame_label_refused");f=original;
        f.constants.electronVoltJ*=2;
        Check(!S::ForwardProducts(registry,f,&result,&error),"changed_constant_value_refused");f=original;
        f.constants.sourceIdentity="another-constant-source";
        Check(!S::ForwardProducts(registry,f,&result,&error),"changed_constant_identity_refused");f=original;
        f.outgoingSodiumIon=State("Na+","unresolved-state","23Na",1,naIon,masses,electron);
        Check(!S::ForwardProducts(registry,f,&result,&error),"wrong_outgoing_state_refused");f=original;
        auto otherIsotopes=masses;otherIsotopes.push_back({"2H",{2*hMass,"analytical-negative-isotope-fixture-not-data"}});
        f.incomingProton=State("H+","bare-proton","2H",1,hIon,otherIsotopes,electron);
        Check(!S::ForwardProducts(registry,f,&result,&error),"wrong_projectile_isotope_refused");f=original;
        f.relativeIncidentDirection={{0,0,0}};
        Check(!S::ForwardProducts(registry,f,&result,&error),"undefined_forward_direction_refused");f=original;
        Check(!S::ForwardProducts(registry,f,nullptr,&error),"null_output_refused");
      }
      f.sourceRate=S::Rate{};
      Check(!S::ForwardProducts(registry,f,&result,&error)&&result.modelIdentity==preserved.modelIdentity,"forged_rate_source_refused_transactionally");
    }
    S::Rate rate;
    Check(!registry.RateAtNativeNode(5000,1,constants,&rate,&error)&&!rate.valid(),"nonnode_gap_refused_unchanged");
    Check(!registry.RateAtNativeNode(3999,1,constants,&rate,&error),"low_energy_extrapolation_refused");
    Check(!registry.RateAtNativeNode(25001,1,constants,&rate,&error),"high_energy_extrapolation_refused");
    Check(!registry.RateAtNativeNode(4000,-1,constants,&rate,&error),"negative_density_refused");
    Check(registry.RateAtNativeNode(4000,0,constants,&rate,&error)&&rate.rateS()==0,"zero_density_exact_zero_rate");
    constants.sourceIdentity.clear();Check(!registry.RateAtNativeNode(4000,1,constants,&rate,&error),"missing_constant_source_refused");
    Check(!registry.LoadNativeNodes("does-not-exist",&error)&&registry.nodes().size()==12,"failed_reload_preserves_authenticated_registry");
    output.flush();Require(output.good(),"failed comparison output");
    std::cout<<"SUMMARY checks="<<checks<<" failed="<<failures<<'\n';return failures?1:0;
  }catch(const std::exception& e){std::cerr<<"ERROR "<<e.what()<<'\n';return 2;}
}
