#include "MoonSodiumTotalCaptureCohort.h"
#include "MoonNativeReactionRegistration.h"
#include "Moon.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace Moon { namespace Chemistry { namespace SodiumTotalCaptureCohort {
namespace {
std::shared_ptr<Prepared> active;
bool Reject(std::string* error,const std::string& message){if(error)*error=message;return false;}
}
bool Prepared::RegisterForNormalStartup(const Request& request,std::shared_ptr<Prepared>* output,std::string* error){
  if(!error||!output||active||request.provider.catalog().size()!=4||request.maximumRootRefinements==0||
     request.provider.frame()!=Moon::Frames::SolarOrbital||request.elementaryChargeC!=1.602176634e-19||
     request.elementaryChargeIdentity.empty())return Reject(error,"require a prepared LSO Na total-capture provider, semantic cohort and explicit solver controls");
  S::Evaluation evaluation;auto status=request.provider.Evaluate(request.beam,{{0,0,0}},&evaluation);
  if(!status.ok())return Reject(error,status.message);
  std::shared_ptr<Prepared> candidate(new Prepared);candidate->request_=request;
  if(!NativeRegistration::SetBeforeParticleBuffer(Install,candidate.get(),S::ModelId,error))return false;
  active=candidate;*output=std::move(candidate);return true;
}
bool Prepared::Install(void* user,std::string* error){
  auto& p=*static_cast<Prepared*>(user);
  if(!SemanticRandom::RegisterStorage(p.request_.random,error))return false;
  Controller::Configuration c;
  // Bind the frozen beam, draws and numerical solver to the migrated history
  // fingerprint. A receipt from another density or velocity cannot resume it.
  std::ostringstream identity;identity.imbue(std::locale::classic());
  identity<<S::ModelId<<";beam-SI="<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<p.request_.beam.protonDensityM3;
  for(auto v:p.request_.beam.protonVelocityMS)identity<<','<<v;
  identity<<";seed="<<p.request_.random.seed<<";cohort="<<p.request_.random.cohort;
  c.identity=identity.str();c.eventFrame=Moon::Frames::SolarOrbital;c.catalog=p.request_.provider.catalog();
  c.elementaryChargeC=p.request_.elementaryChargeC;c.elementaryChargeSourceIdentity=p.request_.elementaryChargeIdentity;
  c.ledger=&p.ledger_;c.userData=&p;c.resolveInitialState=Resolve;c.openUniform=SemanticRandom::DrawForParticle;
  c.initializeProduct=SemanticRandom::InitializeStagedProduct;c.productInitializerIdentity="full-genealogy-HMAC-SHA256-child-key-v1";
  c.exposureHistory=History;c.exposureProviderIdentity=c.identity+";monoenergetic-target-rest-history";
  c.exposureControls.maximumRootRefinements=p.request_.maximumRootRefinements;
  Controller::Channel channel;channel.parentCatalogSpecies=0;channel.reaction=p.request_.provider.reaction();
  channel.scope=Controller::HazardScope::TrajectoryDependent;channel.kinematics=Kinematics;
  channel.hazardEvidenceIdentity=c.exposureProviderIdentity;channel.kinematicProviderIdentity=S::ModelId;c.channels={channel};
  return Controller::Install(c,error);
}
bool Prepared::Resolve(long int particle,Controller::Node*,void* user,std::size_t* species,double* internal,std::string* error){
  const auto& p=*static_cast<Prepared*>(user);SemanticRandom::State key;
  if(!SemanticRandom::Read(particle,&key,error))return false;
  if(PIC::ParticleBuffer::GetI(particle)!=p.request_.provider.catalog()[0].compiledSpecies)
    return Reject(error,"only ground-state neutral sodium may enter as an original capture-cohort particle");
  *species=0;*internal=0;return true; // selected separated neutral atom reference
}
bool Prepared::History(const Controller::ExposureRequest& request,void* user,Controller::ExposureHistory* output,std::string* error){
  const auto& p=*static_cast<Prepared*>(user);
  // The beam rate is constant over this segment only if the Na target stays
  // stationary and unaccelerated in the beam frame (the approved scope).
  // Read the normal force-tracer preview; a varying relative energy cannot be
  // treated as a constant-rate hazard, even when no event is realized.
  if(request.parentCatalogSpecies!=0||request.arc.velocityMS!=R::Vector3{}||
     request.arc.accelerationMS2!=R::Vector3{}||!std::isfinite(request.maximumS)||request.maximumS<0)
    return Reject(error,"I24 total capture requires stationary force-free Na over the entire represented segment");
  S::Evaluation value;const auto status=p.request_.provider.Evaluate(p.request_.beam,request.arc.velocityMS,&value);
  if(!status.ok())return Reject(error,status.message);
  Controller::ExposureHistory result;
  if(request.maximumS>0)result.segments.push_back({0,request.maximumS,{{0,value.rateS,value.rateS}}});
  result.absoluteDepthErrorBound=0; // exact constant rate on this admitted numerical segment
  *output=std::move(result);return true;
}
bool Prepared::Kinematics(const Controller::EventContext& event,R::KinematicPolicy* policy,std::string* error){
  auto& p=*static_cast<Prepared*>(event.userData);
  if(event.channelIndex!=0||!event.node||!event.node->block||event.parent.velocityMS!=R::Vector3{})
    return Reject(error,"unexpected channel, node, or moving parent at a sodium total-capture event");
  // One keyed Kinematics draw selects the final H(n) shell; the angle is fixed.
  double shell=0;
  if(!SemanticRandom::DrawForParticle(Controller::DrawPurpose::Kinematics,event.particle,0,nullptr,&shell,error))return false;
  std::array<R::OpticalDepthClock,2> clocks;std::array<double,2> weights;
  for(std::size_t i=0;i<2;++i){
    double draw=0;
    if(!SemanticRandom::ProspectiveProductClockForParticle(event.particle,i,&draw,error)||
       !R::OpticalDepthClock::FromOpenUniform(draw,&clocks[i],error))return false;
    weights[i]=event.node->block->GetLocalParticleWeight(p.request_.provider.catalog()[i+2].compiledSpecies);
  }
  Receipt receipt;auto status=p.request_.provider.Compose(p.request_.beam,event.parent,event.interval,
      shell,weights,clocks,&receipt.composition);
  if(!status.ok())return Reject(error,status.message);
  if(receipt.composition.evaluation.rateS!=event.selectedInstantaneousRateS)
    return Reject(error,"accepted capture clock and actual-event beam rate disagree");
  if(!SemanticRandom::Read(event.particle,&receipt.key,error))return false;
  receipt.parent=event.parent;receipt.interval=event.interval;receipt.epochEt=event.physicalEpochS;
  receipt.shellUniform=shell;
  *policy=receipt.composition.policy;p.receipts_.push_back(std::move(receipt));
  return true;
}
}}}
