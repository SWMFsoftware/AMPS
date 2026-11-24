// Thin adapter for the stationary-helium production coefficient, moment
// and guarded event composer. All masses and thresholds are loaded from
// the authenticated profile rather than declared here.
#include "srcMoon/species/HeliumElectronImpact.h"
#include <iomanip>
#include <iostream>
#include <limits>

namespace H=Moon::Species::HeliumElectronImpact;
namespace P=Moon::Species::PhysicalStateProfile;
namespace R=Moon::Species::ReactionExecution;

int main(int argc,char** argv){
  if(argc<9||(argc-9)%4!=0)return 2;
  P::Files files;
  files.isotopes={argv[3],argv[4]};files.constants={argv[5],argv[6]};files.states={argv[7],argv[8]};
  for(int i=9;i<argc;i+=4)files.sources.push_back({argv[i],argv[i+1],argv[i+2],argv[i+3]});
  P::Registry profile;
  std::string error;
  if(!profile.Load(files,&error)){
    std::cerr<<error<<'\n';
    return 2;
  }
  H::Provider provider;
  auto status=H::Provider::Prepare(profile,argv[1],argv[2],0,1,&provider);
  if(!status.ok()){
    std::cerr<<status.message<<'\n';
    return status.code==H::Code::Unavailable?77:2;
  }
  H::Bath bath;
  bath.thermalEnergyEv=10;
  bath.frameIdentity="prescribed-inertial-target-rest";
  H::Evaluation value;
  std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10);
  for(double density:{0.,1e6,1e7}){
    bath.densityM3=density;
    status=provider.Evaluate(bath,&value);
    if(!status.ok()){
      std::cerr<<status.message<<'\n';
      return 2;
    }
    std::cout<<"DATA\t"<<density<<'\t'<<value.coefficientM3S<<'\t'<<value.rateS
      <<'\t'<<value.logSlope<<'\t'<<value.meanIncomingEnergyJ<<'\t'
      <<value.fitThresholdEv<<'\t'<<value.physicalThresholdJ<<'\n';
  }
  int failed=0;
  auto check=[&](const char* name,bool pass){
    std::cout<<"CHECK\t"<<name<<'\t'<<(pass?"PASS":"FAIL")<<'\n';
    if(!pass)++failed;
  };
  auto unsupported=bath;
  unsupported.bulkVelocityMS[0]=1; // m/s, deliberately outside zero-drift scope
  check("drift-refusal",provider.Evaluate(unsupported,&value).code==H::Code::OutsideDomain);
  unsupported=bath;
  unsupported.targetVelocityMS[1]=1;
  check("moving-target-refusal",provider.Evaluate(unsupported,&value).code==H::Code::OutsideDomain);
  unsupported=bath;
  unsupported.densityM3=-1;
  const auto previousRate=value.rateS;
  check("invalid-density-preserves-output",
      provider.Evaluate(unsupported,&value).code==H::Code::InvalidInput&&value.rateS==previousRate);
  unsupported=bath;
  unsupported.thermalEnergyEv=11;
  check("unreviewed-temperature-refusal",provider.Evaluate(unsupported,&value).code==H::Code::OutsideDomain);

  // A valid accepted-time/configuration request without source-backed
  // model-error admission must produce no plan. The guard is checked before
  // any product can reach the real controller and allocator.
  R::ParentState parent;
  parent.referenceWeight=2;
  parent.individualCorrection=.75;
  R::TraceInterval interval;
  R::OpticalDepthClock ionClock;
  if(!R::TraceInterval::Create(100,2,1,&interval,&error)||
     !R::OpticalDepthClock::FromOpenUniform(.5,&ionClock,&error))return 2;
  H::Composition event;
  event.policy.identity="guard-sentinel";
  check("missing-evidence-preserves-composition",
      provider.Compose(bath,{},parent,interval,101,2,ionClock,&event).code==H::Code::Unavailable&&
      event.policy.identity=="guard-sentinel");
  const auto& reaction=provider.reaction();
  check("one-incoming-two-outgoing-bath-electrons",
      reaction.reactants.size()==2&&reaction.products.size()==2&&
      reaction.reactants[1].population==R::Population::Bath&&reaction.reactants[1].multiplicity==1&&
      reaction.products[1].population==R::Population::Bath&&reaction.products[1].multiplicity==2);
  return failed?1:0;
}
