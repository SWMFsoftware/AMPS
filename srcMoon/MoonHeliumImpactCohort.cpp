#include "MoonHeliumImpactCohort.h"
#include "MoonNativeReactionRegistration.h"
#include "Moon.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace Moon { namespace Chemistry { namespace HeliumImpactCohort {
const char OwnerIdentity[]="I23-prescribed-electron-bath-4He-ground-ion-reservoir-v1";
namespace {
std::shared_ptr<Prepared> active;
bool Reject(std::string* error,const std::string& message){if(error)*error=message;return false;}
std::string Status(const H::Status& status){return status.message.empty()?"unspecified helium provider refusal":status.message;}
}

bool Prepared::RegisterForNormalStartup(const Request& request,std::shared_ptr<Prepared>* output,std::string* error){
  if(!error||!output)return false;
  if(active)return Reject(error,"an I23 helium electron-bath owner is already registered");
  const auto& catalog=request.provider.catalog();
  // Catalog order is fixed by HeliumElectronImpact::Prepare: 4He, 4He+, e.
  if(catalog.size()!=3||catalog[0].compiledSpecies<0||catalog[1].compiledSpecies<0||catalog[2].compiledSpecies!=-1)
    return Reject(error,"require a prepared 4He/4He+/bath-electron provider with compiled He and He+ species");
  if(request.bath.frameIdentity!=Moon::Frames::SolarOrbital||request.admission.planner.eventFrame!=Moon::Frames::SolarOrbital||
     request.admission.evidence.frameIdentity!=Moon::Frames::SolarOrbital)
    return Reject(error,"the prescribed bath, admission and event frame must all be the frozen LSO frame");
  if(request.admission.planner.identity.empty())
    return Reject(error,"the reviewed admission must name its kinematic provider identity");
  if(request.maximumRootRefinements==0||request.elementaryChargeC!=1.602176634e-19||request.elementaryChargeIdentity.empty())
    return Reject(error,"explicit exposure root-refinement effort and the exact SI elementary charge are required");
  // Evaluate validates density/temperature/drift/target motion. Invalid
  // plasma is refused here, before the particle buffer exists.
  H::Evaluation evaluation;
  auto status=request.provider.Evaluate(request.bath,&evaluation);
  if(!status.ok())return Reject(error,"prescribed electron bath refused: "+Status(status));
  if(evaluation.rateS>0){
    // Admission preflight on a stationary reference parent: the composer
    // must accept the reviewed evidence before particles move. The result is
    // discarded; it consumes no random draw and creates no particle.
    R::ParentState parent;parent.species=0;parent.referenceWeight=1;parent.individualCorrection=1;
    R::TraceInterval interval;R::OpticalDepthClock clock;std::string local;
    if(!R::TraceInterval::Create(0,1,.5,&interval,&local)||!R::OpticalDepthClock::FromOpenUniform(.5,&clock,&local))
      return Reject(error,"internal admission preflight setup failed: "+local);
    H::Composition probe;
    status=request.provider.Compose(request.bath,request.admission,parent,interval,.5,1,clock,&probe);
    if(!status.ok())return Reject(error,"reviewed helium admission refused before startup: "+Status(status));
  }
  std::shared_ptr<Prepared> candidate(new Prepared);candidate->request_=request;
  if(!NativeRegistration::SetBeforeParticleBuffer(Install,candidate.get(),OwnerIdentity,error))return false;
  active=candidate;*output=std::move(candidate);return true;
}

bool Prepared::Install(void* user,std::string* error){
  auto& p=*static_cast<Prepared*>(user);
  if(!SemanticRandom::RegisterStorage(p.request_.random,error))return false;
  Controller::Configuration c;
  // Bind the frozen bath, admission, draws and numerical solver to the
  // migrated history fingerprint. A receipt with another density,
  // temperature or evidence identity cannot resume this history.
  std::ostringstream identity;identity.imbue(std::locale::classic());
  identity<<OwnerIdentity<<";bath-kTe-eV="<<std::setprecision(std::numeric_limits<double>::max_digits10)
          <<p.request_.bath.thermalEnergyEv<<";ne-m3="<<p.request_.bath.densityM3
          <<";frame="<<p.request_.bath.frameIdentity
          <<";crosswalk="<<p.request_.admission.evidence.stateCrosswalkSha256
          <<";model-error="<<p.request_.admission.evidence.modelErrorEvidenceSha256
          <<";seed="<<p.request_.random.seed<<";cohort="<<p.request_.random.cohort;
  c.identity=identity.str();c.eventFrame=Moon::Frames::SolarOrbital;c.catalog=p.request_.provider.catalog();
  c.elementaryChargeC=p.request_.elementaryChargeC;c.elementaryChargeSourceIdentity=p.request_.elementaryChargeIdentity;
  c.ledger=&p.ledger_;c.userData=&p;c.resolveInitialState=Resolve;c.openUniform=SemanticRandom::DrawForParticle;
  c.initializeProduct=SemanticRandom::InitializeStagedProduct;c.productInitializerIdentity="full-genealogy-HMAC-SHA256-child-key-v1";
  c.exposureHistory=History;c.exposureProviderIdentity=c.identity+";stationary-target-uniform-constant-bath-history";
  c.exposureControls.maximumRootRefinements=p.request_.maximumRootRefinements;
  Controller::Channel channel;channel.parentCatalogSpecies=0;channel.reaction=p.request_.provider.reaction();
  // TrajectoryDependent scope makes the controller ask History for EVERY
  // accepted segment, so the stationary/force-free domain is enforced on the
  // actual flight rather than assumed from the initial state.
  channel.scope=Controller::HazardScope::TrajectoryDependent;channel.kinematics=Kinematics;
  channel.hazardEvidenceIdentity=c.exposureProviderIdentity;
  channel.kinematicProviderIdentity=p.request_.admission.planner.identity;c.channels={channel};
  return Controller::Install(c,error);
}

bool Prepared::Resolve(long int particle,Controller::Node*,void* user,std::size_t* species,double* internal,std::string* error){
  const auto& p=*static_cast<Prepared*>(user);SemanticRandom::State key;
  if(!SemanticRandom::Read(particle,&key,error))return false;
  if(PIC::ParticleBuffer::GetI(particle)!=p.request_.provider.catalog()[0].compiledSpecies)
    return Reject(error,"only ground-state neutral 4He may enter as an original I23 electron-bath parent");
  *species=0;*internal=0;return true; // separated neutral ground atom; formation-energy zero
}

bool Prepared::History(const Controller::ExposureRequest& request,void* user,Controller::ExposureHistory* output,std::string* error){
  const auto& p=*static_cast<Prepared*>(user);
  // The Voronov Maxwellian coefficient applies to a target at rest in the
  // bath frame. Read the normal force-tracer preview: a moving or accelerated
  // He parent is outside the admitted domain even if no event is realized.
  if(request.parentCatalogSpecies!=0||request.arc.velocityMS!=R::Vector3{}||
     request.arc.accelerationMS2!=R::Vector3{}||!std::isfinite(request.maximumS)||request.maximumS<0)
    return Reject(error,"I23 electron-bath hazard requires a stationary, force-free 4He parent over the entire segment");
  H::Evaluation value;const auto status=p.request_.provider.Evaluate(p.request_.bath,&value);
  if(!status.ok())return Reject(error,Status(status));
  Controller::ExposureHistory result;
  // Spatially uniform, time-constant bath: one exact constant segment.
  // A zero-density bath yields a zero hazard and no event.
  if(request.maximumS>0)result.segments.push_back({0,request.maximumS,{{0,value.rateS,value.rateS}}});
  result.absoluteDepthErrorBound=0;
  *output=std::move(result);return true;
}

bool Prepared::Kinematics(const Controller::EventContext& event,R::KinematicPolicy* policy,std::string* error){
  auto& p=*static_cast<Prepared*>(event.userData);
  if(event.channelIndex!=0||!event.node||!event.node->block||event.parent.velocityMS!=R::Vector3{})
    return Reject(error,"unexpected channel, node, or moving parent at an I23 electron-impact event");
  // The He+ product receives one fresh, prospectively keyed optical clock.
  // The reduced reservoir composer itself draws no random number.
  double draw=0;R::OpticalDepthClock clock;
  if(!SemanticRandom::ProspectiveProductClockForParticle(event.particle,0,&draw,error)||
     !R::OpticalDepthClock::FromOpenUniform(draw,&clock,error))return false;
  const double ionWeight=event.node->block->GetLocalParticleWeight(p.request_.provider.catalog()[1].compiledSpecies);
  Receipt receipt;
  H::Evaluation value;auto status=p.request_.provider.Evaluate(p.request_.bath,&value);
  if(!status.ok())return Reject(error,Status(status));
  if(static_cast<long double>(value.rateS)!=event.selectedInstantaneousRateS)
    return Reject(error,"accepted electron-bath clock and actual-event n_e*k(T_e) disagree");
  status=p.request_.provider.Compose(p.request_.bath,p.request_.admission,event.parent,event.interval,
      event.physicalEpochS,ionWeight,clock,&receipt.composition);
  if(!status.ok())return Reject(error,Status(status));
  if(!SemanticRandom::Read(event.particle,&receipt.key,error))return false;
  receipt.parent=event.parent;receipt.interval=event.interval;receipt.epochEt=event.physicalEpochS;receipt.rateS=value.rateS;
  *policy=receipt.composition.policy;p.receipts_.push_back(std::move(receipt));
  return true;
}
}}}
