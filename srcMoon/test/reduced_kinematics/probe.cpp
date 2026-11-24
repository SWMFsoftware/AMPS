// Thin production adapter. Analytical kg/J fixtures below are mathematical
// controls, not a source of isotope masses or lunar product energies.
#include "srcMoon/species/ReducedKinematics.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
namespace K=Moon::Species::ReducedKinematics;
using Real=K::Real;using Vec=K::Vector3;
int Checks=0,Failed=0;
void Check(const std::string& name,bool pass){++Checks;if(!pass)++Failed;std::cout<<(pass?"PASS ":"FAIL ")<<name<<'\n';}
Real Norm(const Vec& v){return std::hypot(std::hypot(v[0],v[1]),v[2]);}
std::vector<K::IsotopeMass> Inventory() {
  return {{"X",{1,"ANALYTICAL-MASS-X"}},{"H",{1,"ANALYTICAL-MASS-H"}},
      {"O",{16,"ANALYTICAL-MASS-O"}},{"B",{2,"ANALYTICAL-MASS-B"}},
      {"heavy",{1e100L,"ANALYTICAL-EXTREME-HEAVY"}},{"light",{1e-100L,"ANALYTICAL-EXTREME-LIGHT"}}};
}
K::ResolvedState Resolve(const std::string& name,
    std::vector<Moon::Species::ReactionExecution::NucleusCount> nuclei,
    std::int64_t charge=0,Real formation=0,Real excitation=0) {
  K::StateDefinition d;d.identity=name;d.state="explicit-analytical-state";d.nuclei=std::move(nuclei);
  d.elementaryCharge=charge;d.formation={formation,"ANALYTICAL-FORMATION-"+name};
  d.excitation={excitation,"ANALYTICAL-EXCITATION-"+name};
  K::ResolvedState s;std::string error;
  if(!K::ResolveState(d,Inventory(),{.001L,"ANALYTICAL-ELECTRON-MASS"},&s,&error))throw std::runtime_error(error);
  return s;
}
K::NumericControls Controls() {
  const Real eps=std::numeric_limits<Real>::epsilon();K::NumericControls c;
  c.maximumClassicalSpeedMS=1000000;c.unitVectorTolerance=64*eps;
  c.referenceMassKg={0,128*eps};c.momentumKgMS={0,128*eps};
  c.energyJ={0,128*eps};c.thresholdJ={0,128*eps};return c;
}
void Emit(const std::string& name,const K::FinalState& state) {
  std::cout<<"EVENT\t"<<name<<'\t'<<state.availableCOMKineticJ<<'\t'<<state.actualCOMKineticJ
      <<'\t'<<state.totalKineticJ<<'\t'<<state.photonEnergyJ<<'\t'<<state.stateQJ
      <<'\t'<<state.channelThresholdJ<<'\t'<<state.sourceMeanExcessJ<<'\t'<<state.stateAvailableKineticJ<<'\n';
  for(std::size_t i=0;i<state.products.size();++i) {
    const auto& p=state.products[i];std::cout<<"PRODUCT\t"<<name<<'\t'<<i<<'\t'<<p.state.referenceMassKg();
    for(auto x:p.momentumKgMS)std::cout<<'\t'<<x;
    for(auto x:p.velocityMS)std::cout<<'\t'<<x;
    std::cout<<'\n';
  }
  Check(name+"/source-provenance",!state.sourceIdentities.empty()&&!state.frameIdentity.empty());
  for(std::size_t i=0;i<state.products.size();++i)
    Check(name+"/occurrence-"+std::to_string(i),state.products[i].occurrence==i&&state.products[i].state.valid());
}
template<class Input,class Function>
void Refusal(const std::string& name,const Input& input,Function function) {
  K::FinalState output;output.modelIdentity="UNCHANGED";output.frameIdentity="SENTINEL";
  output.sourceIdentities={"UNCHANGED-SOURCE"};output.energyResidualJ=-99;
  std::string error;const bool ok=function(input,&output,&error);
  Check(name,!ok&&!error.empty()&&output.modelIdentity=="UNCHANGED"&&output.frameIdentity=="SENTINEL"&&
      output.energyResidualJ==-99&&output.products.empty()&&output.sourceIdentities==std::vector<std::string>{"UNCHANGED-SOURCE"});
}
int main() {
  try {
    std::cout<<std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout<<"PLATFORM\t"<<std::numeric_limits<Real>::epsilon()<<'\t'<<std::numeric_limits<Real>::digits<<'\t'<<std::numeric_limits<Real>::max_exponent<<'\n';
    auto x=Resolve("X",{{"X",1}}),x2=Resolve("X2",{{"X",2}},0,-3);
    auto ion=Resolve("X+",{{"X",1}},1,5),electron=Resolve("electron",{},-1);
    auto h=Resolve("H",{{"H",1}}),o=Resolve("O",{{"O",1}}),b=Resolve("B",{{"B",1}});
    Check("reference-neutral-mass",x.referenceMassKg()==1);
    Check("reference-ion-electron-additive",std::fabs(ion.referenceMassKg()+electron.referenceMassKg()-x.referenceMassKg())<=4*std::numeric_limits<Real>::epsilon());
    K::ResolvedState kept=x;auto bad=x.definition();bad.nuclei={{"unknown",1}};std::string error;
    Check("unknown-isotope-preserves-state",!K::ResolveState(bad,Inventory(),{.001L,"E"},&kept,&error)&&kept.referenceMassKg()==1&&kept.definition().identity=="X");
    bad=x.definition();bad.excitation.sourceIdentity.clear();
    Check("missing-energy-source-refused",!K::ResolveState(bad,Inventory(),{.001L,"E"},&kept,&error));
    bad=x.definition();bad.nuclei={{"X",1},{"X",1}};
    Check("duplicate-nucleus-refused",!K::ResolveState(bad,Inventory(),{.001L,"E"},&kept,&error));
    auto inventory=Inventory();inventory.push_back(inventory.front());
    Check("duplicate-isotope-refused",!K::ResolveState(x.definition(),inventory,{.001L,"E"},&kept,&error));
    K::FormationBalance q;
    Check("formation-bookkeeping",K::FormationQ({x2},{x,x},&q,&error)&&q.qJ==-3&&q.outgoingStateEnergyJ==0&&q.incomingStateEnergyJ==-3);
    K::FormationBalance untouched;untouched.qJ=-99;
    Check("formation-failure-preserves-output",!K::FormationQ({}, {x},&untouched,&error)&&untouched.qJ==-99);
    auto excited=Resolve("Xexc",{{"X",1}},0,-2,3);
    Check("formation-excitation-counted-once",K::FormationQ({excited},{x},&q,&error)&&q.qJ==1);
    auto large=Resolve("large-formation",{{"X",1}},0,1e20L);
    auto largeExcited=Resolve("large-plus-small",{{"X",1}},0,1e20L,1);
    Check("formation-cancellation-retains-explicit-excitation",K::FormationQ({largeExcited},{large},&q,&error)&&q.qJ==1&&q.terms.size()==2);
    auto huge=Resolve("huge-formation",{{"X",1}},0,1e100L),hugeExcited=Resolve("huge-plus-one",{{"X",1}},0,1e100L,1);
    auto middle=Resolve("middle-formation",{{"X",1}},0,1e50L);
    Check("formation-multiscale-cancellation",K::FormationQ({hugeExcited,middle},{huge,middle},&q,&error)&&q.qJ==1);
    K::ModelErrorReport model;
    Check("independent-model-error-report",K::EvaluateModelErrors(x,{1.001L,"ANALYTICAL-PHYSICAL-MASS"},K::SpeedOfLightMS/10,{.01L,.01L,.01L},&model,&error));
    std::cout<<"ERRORS\t"<<model.relativeInertialMassError<<'\t'<<model.relativeAccelerationError<<'\t'<<model.relativeRelativisticKineticError<<'\n';
    std::cout<<"MODELINPUT\t"<<std::hexfloat<<1.001L<<'\t'<<K::SpeedOfLightMS/10<<std::defaultfloat<<'\n';
    model.speedMS=-7;
    Check("model-error-limit-no-widening",!K::EvaluateModelErrors(x,{1.001L,"ANALYTICAL-PHYSICAL-MASS"},K::SpeedOfLightMS/10,{.0001L,.01L,.01L},&model,&error)&&model.speedMS==-7);
    Check("nonrelativistic-domain-refusal",!K::EvaluateModelErrors(x,{1,"M"},K::SpeedOfLightMS,{.01L,.01L,.01L},&model,&error));
    Check("positive-model-error-underflow-not-zero",!K::EvaluateModelErrors(x,{1,"M"},std::numeric_limits<Real>::denorm_min(),{.01L,.01L,.01L},&model,&error));
    K::TwoBodyInput two;two.products={{x,b}};two.totalMomentumKgMS={{1,-2,3}};
    two.availableCOMKineticJ=4;two.relativeDirection={{0,0,1}};
    two.context={"ANALYTICAL-INERTIAL","ANALYTICAL-TWO-BODY"};two.controls=Controls();
    K::FinalState event;
    if(K::GeneralTwoBodyFinalState(two,&event,&error)){Check("general-two-body",true);Emit("general-two-body",event);}else{Check("general-two-body",false);std::cerr<<error<<'\n';}
    auto near=two;near.relativeDirection[2]+=8*std::numeric_limits<Real>::epsilon();K::FinalState normalized;
    bool sameDirection=K::GeneralTwoBodyFinalState(near,&normalized,&error)&&event.products.size()==2;
    if(sameDirection)for(int d=0;d<3;++d)sameDirection=sameDirection&&normalized.products[0].momentumKgMS[d]==event.products[0].momentumKgMS[d];
    Check("accepted-direction-roundoff-does-not-change-energy",sameDirection);
    auto invalidTwo=two;invalidTwo.availableCOMKineticJ=-1;Refusal("negative-COM-energy",invalidTwo,K::GeneralTwoBodyFinalState);
    invalidTwo=two;invalidTwo.relativeDirection={{0,0,2}};Refusal("nonunit-relative-direction",invalidTwo,K::GeneralTwoBodyFinalState);
    invalidTwo=two;invalidTwo.controls.maximumClassicalSpeedMS=1e-10L;Refusal("explicit-speed-guard",invalidTwo,K::GeneralTwoBodyFinalState);
    invalidTwo=two;invalidTwo.controls=K::NumericControls();Refusal("omitted-bounds-refused",invalidTwo,K::GeneralTwoBodyFinalState);
    invalidTwo=two;invalidTwo.totalMomentumKgMS={{1e100L,0,0}};Refusal("unrepresentable-small-recoil-on-large-drift",invalidTwo,K::GeneralTwoBodyFinalState);
    invalidTwo=two;invalidTwo.relativeDirection={{0,0,0}};invalidTwo.controls.unitVectorTolerance=2;
    Refusal("zero-relative-direction-even-with-wide-bound",invalidTwo,K::GeneralTwoBodyFinalState);
    invalidTwo=two;invalidTwo.controls.unitVectorTolerance=1;
    Refusal("unit-normalization-bound-must-be-less-than-one",invalidTwo,K::GeneralTwoBodyFinalState);
    K::PhotonTwoBodyInput photon;photon.parent=x2;photon.products={{x,x}};
    photon.channelThreshold={3,"ANALYTICAL-THRESHOLD"};photon.meanExcess={2,"ANALYTICAL-MEAN"};
    photon.photonDirection={{1,0,0}};photon.polarOpenUniform=.25L;photon.azimuthOpenUniform=.125L;
    photon.context={"ANALYTICAL-PARENT-REST-INERTIAL","ANALYTICAL-PHOTON"};photon.controls=Controls();
    if(K::PhotonTwoBody(photon,&event,&error)){Check("neutral-photon",true);Emit("neutral-photon",event);Check("photon-momentum-included",event.photonIncluded&&event.photonMomentumKgMS[0]>0);}else{Check("neutral-photon",false);std::cerr<<error<<'\n';}
    auto atomic=photon;atomic.parent=x;atomic.products={{ion,electron}};atomic.channelThreshold={5,"ANALYTICAL-IONIZATION-THRESHOLD"};
    if(K::PhotonTwoBody(atomic,&event,&error)){Check("atomic-ion-electron-photon",true);Emit("atomic-ion-electron-photon",event);}else{Check("atomic-ion-electron-photon",false);std::cerr<<error<<'\n';}
    auto tinyPolar=photon;tinyPolar.polarOpenUniform=1e-25L;
    if(K::PhotonTwoBody(tinyPolar,&event,&error)){Check("tiny-open-polar",true);Emit("tiny-open-polar",event);}else{Check("tiny-open-polar",false);std::cerr<<error<<'\n';}
    auto nearPhoton=photon;nearPhoton.photonDirection[0]+=8*std::numeric_limits<Real>::epsilon();
    Check("accepted-photon-direction-preserves-dispersion",K::PhotonTwoBody(nearPhoton,&normalized,&error)&&normalized.photonMomentumKgMS[0]==normalized.photonEnergyJ/K::SpeedOfLightMS);
    auto tiny=photon;tiny.parent=Resolve("X2-tiny-energy",{{"X",2}},0,-1e-100L);
    tiny.channelThreshold={1e-100L,"ANALYTICAL-TINY-THRESHOLD"};tiny.meanExcess={1e-200L,"ANALYTICAL-TINY-MEAN"};
    if(K::PhotonTwoBody(tiny,&event,&error)){Check("separate-small-excess-preserved",true);Emit("separate-small-excess-preserved",event);}else{Check("separate-small-excess-preserved",false);std::cerr<<error<<'\n';}
    auto badPhoton=photon;badPhoton.parentMomentumKgMS={{1e-100L,0,0}};Refusal("nonzero-parent-momentum-is-not-rest",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.channelThreshold.joules=2;Refusal("threshold-state-crosswalk-mismatch",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.channelThreshold.sourceIdentity.clear();Refusal("missing-threshold-source",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.meanExcess.joules=0;Refusal("negative-photon-recoil-budget-not-clipped",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.polarOpenUniform=0;Refusal("closed-polar-uniform",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.azimuthOpenUniform=1;Refusal("closed-azimuth-uniform",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.photonDirection={{0,0,0}};Refusal("invalid-photon-direction",badPhoton,K::PhotonTwoBody);
    badPhoton.controls.unitVectorTolerance=2;Refusal("zero-photon-direction-even-with-wide-bound",badPhoton,K::PhotonTwoBody);
    badPhoton=photon;badPhoton.parent=x;Refusal("incomplete-material-products",badPhoton,K::PhotonTwoBody);
    // Prescribed analytical fractional energies and directions exercise the
    // charged three-product map. These controls make no claim about the
    // spectrum or states populated by photodissociation of lunar water.
    K::ChargedBreakupInput charged;
    charged.electron=electron;charged.ion=ion;charged.neutral=b;
    charged.availableCOMKineticJ=4;charged.electronRelativeEnergyFraction=.25L;
    charged.electronRelativeToPair={{0,1,0}};charged.ionRelativeToNeutral={{.6L,0,.8L}};
    charged.context={"ANALYTICAL-INERTIAL","ANALYTICAL-COM-ENERGY"};charged.controls=Controls();
    charged.partitionSourceIdentity="ANALYTICAL-DECLARED-PARTITION";
    auto chargedCase=[&](const std::string& name,const K::ChargedBreakupInput& input) {
      const bool ok=K::ChargedBreakupFinalState(input,&event,&error);
      Check(name,ok);
      if(ok)Emit(name,event);else std::cerr<<name<<": "<<error<<'\n';
    };
    chargedCase("charged-rest",charged);
    charged.totalMomentumKgMS={{3,-2,1}};
    chargedCase("charged-drift",charged);
    auto endpoint=charged;endpoint.electronRelativeEnergyFraction=0;endpoint.electronRelativeToPair={};
    chargedCase("charged-heavy-only",endpoint);
    endpoint=charged;endpoint.electronRelativeEnergyFraction=1;endpoint.ionRelativeToNeutral={};
    chargedCase("charged-electron-only",endpoint);
    endpoint.availableCOMKineticJ=0;endpoint.electronRelativeToPair={};
    chargedCase("charged-at-threshold",endpoint);
    auto badCharged=charged;badCharged.partitionSourceIdentity.clear();
    Refusal("unsourced-charged-partition",badCharged,K::ChargedBreakupFinalState);
    badCharged=charged;badCharged.electronRelativeEnergyFraction=1.01L;
    Refusal("nonphysical-charged-partition",badCharged,K::ChargedBreakupFinalState);
    badCharged=charged;std::swap(badCharged.ion,badCharged.neutral);
    Refusal("wrong-charged-product-class",badCharged,K::ChargedBreakupFinalState);
    badCharged=charged;badCharged.electronRelativeToPair={};
    Refusal("missing-active-mode-direction",badCharged,K::ChargedBreakupFinalState);
    badCharged=charged;badCharged.electronRelativeEnergyFraction=1e-100L;
    Refusal("charged-small-mode-lost-in-boost",badCharged,K::ChargedBreakupFinalState);
    K::PhotonChargedBreakupInput fragment;
    fragment.parent=Resolve("XB",{{"X",1},{"B",1}},0,-3);
    fragment.products={{electron,ion,b}};
    fragment.partition=charged;
    fragment.channelThreshold={8,"ANALYTICAL-XB-THRESHOLD"};fragment.meanExcess={4,"ANALYTICAL-XB-EXCESS"};
    fragment.photonDirection={{1,0,0}};fragment.context=charged.context;fragment.controls=charged.controls;
    if(K::PhotonChargedBreakup(fragment,&event,&error)) {
      Check("charged-photon",true);Emit("charged-photon",event);
      Check("charged-complete-state-terms",event.photonIncluded&&event.formationBalance.terms.size()==4&&
            event.products.size()==3&&event.products[0].state.definition().elementaryCharge==-1);
    }
    else {Check("charged-photon",false);std::cerr<<error<<'\n';}
    auto badFragment=fragment;badFragment.products[2]=o;
    Refusal("charged-photon-isotope-mismatch",badFragment,K::PhotonChargedBreakup);
    badFragment=fragment;badFragment.parentMomentumKgMS[0]=1;
    Refusal("charged-photon-nonrest-parent",badFragment,K::PhotonChargedBreakup);
    badFragment=fragment;badFragment.channelThreshold.joules-=1;
    Refusal("charged-photon-incompatible-threshold",badFragment,K::PhotonChargedBreakup);
    K::ThreeBodyInput three;three.products={{o,h,h}};three.availableCOMKineticJ=1;
    three.standardNormals={{{{1,2,-1}},{{2,-1,.5L}},{{-.2L,.7L,-2}}}};
    three.context={"ANALYTICAL-INERTIAL","ANALYTICAL-THREE-BODY"};three.controls=Controls();
    if(K::ThreeNeutralMicrocanonical(three,&event,&error)){Check("microcanonical-deterministic",true);Emit("microcanonical-deterministic",event);Check("identical-H-retained",event.products.size()==3&&event.products[1].state.definition().identity=="H"&&event.products[2].state.definition().identity=="H");}else{Check("microcanonical-deterministic",false);std::cerr<<error<<'\n';}
    auto swap=three;std::swap(swap.standardNormals[1],swap.standardNormals[2]);K::FinalState swapped;
    bool symmetric=K::ThreeNeutralMicrocanonical(swap,&swapped,&error)&&event.products.size()==3;
    if(symmetric)for(int d=0;d<3;++d)symmetric=symmetric&&std::fabs(swapped.products[1].momentumKgMS[d]-event.products[2].momentumKgMS[d])<256*std::numeric_limits<Real>::epsilon();
    Check("identical-H-exchange-symmetry",symmetric);
    auto zero=three;zero.availableCOMKineticJ=0;zero.standardNormals={};zero.totalMomentumKgMS={{18,0,0}};
    if(K::ThreeNeutralMicrocanonical(zero,&event,&error)){Check("zero-energy-deterministic",true);Emit("zero-energy-deterministic",event);}else Check("zero-energy-deterministic",false);
    auto badThree=three;badThree.standardNormals={};Refusal("degenerate-positive-energy-draw",badThree,K::ThreeNeutralMicrocanonical);
    badThree=three;badThree.standardNormals[0][0]=std::numeric_limits<Real>::quiet_NaN();Refusal("nonfinite-normal",badThree,K::ThreeNeutralMicrocanonical);
    badThree=three;badThree.products[0]=ion;Refusal("charged-three-body-not-authorized",badThree,K::ThreeNeutralMicrocanonical);
    auto extreme=three;extreme.products={{Resolve("heavy",{{"heavy",1}}),Resolve("light",{{"light",1}}),Resolve("light",{{"light",1}})}};
    extreme.availableCOMKineticJ=1e-200L;
    if(std::numeric_limits<Real>::max_exponent>=1400) {
      if(K::ThreeNeutralMicrocanonical(extreme,&event,&error)){Check("extreme-mass-energy-scales",true);Emit("extreme-mass-energy-scales",event);}else{Check("extreme-mass-energy-scales",false);std::cerr<<error<<'\n';}
    }
    else Refusal("extreme-individual-energy-range-refused",extreme,K::ThreeNeutralMicrocanonical);
    const int trials=40000;std::mt19937_64 generator(731993);std::normal_distribution<Real> normal;
    Real sumX=0,sumX2=0,energyShare[3]={},cosine2[3]={};int cdf[9]={},completed=0;
    for(int sample=0;sample<trials;++sample) {
      for(auto& g:three.standardNormals)for(auto& component:g)component=normal(generator);
      if(!K::ThreeNeutralMicrocanonical(three,&event,&error)){std::cerr<<error<<'\n';break;}
      Vec difference;for(int d=0;d<3;++d)difference[d]=event.products[1].velocityMS[d]-event.products[2].velocityMS[d];
      const Real jacobi=.25L*Norm(difference)*Norm(difference); // mu(H,H)=1/2kg, K=1J
      sumX+=jacobi;sumX2+=jacobi*jacobi;
      for(int cut=1;cut<10;++cut)if(jacobi<=static_cast<Real>(cut)/10)++cdf[cut-1];
      for(int i=0;i<3;++i) {
        const auto& p=event.products[i];energyShare[i]+=.5L*p.state.referenceMassKg()*Norm(p.velocityMS)*Norm(p.velocityMS);
      }
      const Real length=Norm(difference);for(int d=0;d<3;++d)cosine2[d]+=difference[d]*difference[d]/(length*length);
      ++completed;
    }
    Check("microcanonical-ensemble-all-events-complete",completed==trials);
    std::cout<<"ENSEMBLE\t"<<completed<<'\t'<<sumX/completed<<'\t'<<sumX2/completed;
    for(auto value:energyShare)std::cout<<'\t'<<value/completed;
    for(auto value:cosine2)std::cout<<'\t'<<value/completed;
    for(auto value:cdf)std::cout<<'\t'<<value;
    std::cout<<"\nSUMMARY checks="<<Checks<<" failed="<<Failed<<" seed=731993 trials="<<trials<<'\n';
    return Failed?1:0;
  }
  catch(const std::exception& exception){std::cerr<<"PROBE ERROR: "<<exception.what()<<'\n';return 2;}
}
