#include "MoonSodiumCohort.h"
#include "MoonNativeReactionRegistration.h"
#include "MoonSurfaceEnvironment.h"
#include "Moon.h"
#include "MoonTransportClock.h"
#include "src/general/input/amps_input.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
namespace Moon { namespace Chemistry { namespace SodiumCohort {
namespace D=Species::ReactionData;
namespace T=Species::TrajectoryExposure;
namespace {
// Lifetime owner for the one explicitly installed normal startup binding;
// process-local pointers never enter particle bytes or semantic keys.
std::shared_ptr<Prepared> installedOwner;
E::Status Error(E::Code c,const std::string&s){return {c,s};}
bool Positive(double x){return std::isfinite(x)&&x>0;}
bool NativeProfileManifest(const Request& r,std::string* error){
  if(!D::VerifyFileSha256(r.profileNativeManifestPath,r.profileNativeManifestSha256,error))return false;
  AMPS::Input::Document doc;const auto parsed=AMPS::Input::ParseFile(r.profileNativeManifestPath,&doc);
  if(!parsed.ok()){*error=parsed.message;return false;}
  if(doc.root.blocks.size()!=1||doc.root.blocks[0].kind!="physical_state_profile"||!doc.root.assignments.empty()){
    *error="native profile requires one physical_state_profile block";return false;}
  const auto& block=doc.root.blocks[0];std::set<std::string> used;
  auto text=[&](const std::string& key,std::string* value){const auto*a=block.Find(key);if(!a){*error="native profile missing "+key;return false;}used.insert(key);*value=a->value;return true;};
  auto path=[&](const AMPS::Input::Assignment&a){auto s=a.value;if(s.size()>=2&&s.front()=='"'&&s.back()=='"')s=s.substr(1,s.size()-2);std::filesystem::path p(s);if(p.is_relative())p=std::filesystem::path(a.where.file).parent_path()/p;return std::filesystem::absolute(p).lexically_normal().string();};
  auto sameFile=[&](const std::string& key,const std::string& expectedPath,const std::string& expectedHash){std::string hash;if(!text(key+"_sha256",&hash))return false;const auto*a=block.Find(key+"_file");if(!a){*error="native profile missing "+key+"_file";return false;}used.insert(key+"_file");
    if(hash!=expectedHash||path(*a)!=std::filesystem::absolute(expectedPath).lexically_normal().string()){*error="native profile projection does not bind exact request "+key;return false;}return true;};
  std::string schema;
  if(!text("schema",&schema))return false;
  if(schema!="srcMoon-M3-state-profile-native-v1"){*error="unsupported native state-profile projection schema";return false;}
  if(!sameFile("profile_json",r.profileMetadataPath,r.profileMetadataSha256)||
     !sameFile("isotopes",r.profile.isotopes.path,r.profile.isotopes.sha256)||
     !sameFile("constants",r.profile.constants.path,r.profile.constants.sha256)||
     !sameFile("states",r.profile.states.path,r.profile.states.sha256))return false;
  for(const auto&a:block.assignments)if(!used.count(a.key)){*error="unknown native profile key "+a.key;return false;}
  if(block.blocks.size()!=r.profile.sources.size()){*error="native profile source count/order differs from explicit request";return false;}
  for(std::size_t i=0;i<block.blocks.size();++i){const auto&b=block.blocks[i];const auto&s=r.profile.sources[i];
    if(b.kind!="source"||b.name!=s.id||!b.blocks.empty()||b.assignments.size()!=3){*error="native profile ordered source identity/shape mismatch";return false;}
    const auto*p=b.Find("path"),*h=b.Find("sha256"),*role=b.Find("role");
    if(!p||!h||!role||path(*p)!=std::filesystem::absolute(s.path).lexically_normal().string()||h->value!=s.sha256||role->value!=s.role){*error="native profile source path/hash/role differs from exact request";return false;}
  }
  error->clear();return true;
}
bool Crosswalk(const Request& r,const Species::Photochemistry::Registry& source,
    const E::Definition& d,std::string*e){
  std::string bytes;if(!D::VerifyFileSha256(r.crosswalkPath,r.crosswalkSha256,e)||
      !D::ReadFileBytes(r.crosswalkPath,0,&bytes,e))return false;
  D::Csv rows;const D::Row header={"channel_id","parent_state_id","ion_state_id","electron_state_id","energy_reference_identity"};
  const D::Row expected={d.channelId,r.parentStateId,r.ionStateId,r.electronStateId,d.energyReferenceIdentity};
  if(d.photonEnergyPolicy==Species::ReducedKinematics::PhotonEnergyPolicy::SourceEqualsState){
  if(!D::ParseCsv(bytes,header,&rows,e))return false;
  if(rows.size()!=1||rows[0]!=expected){
    *e="reviewed crosswalk must bind exactly this complete source channel, three declared state IDs and common energy reference; no inferred branch/state mapping";return false;}
  }else if(d.photonEnergyPolicy==Species::ReducedKinematics::PhotonEnergyPolicy::SourceEffectivePhoton){
    D::Row extended=header;
    extended.insert(extended.end(),{"photon_energy_policy","source_threshold_eV","physical_threshold_source_identity"});
    if(!D::ParseCsv(bytes,extended,&rows,e))return false;
    if(rows.size()!=1||!std::equal(expected.begin(),expected.end(),rows[0].begin())||
        rows[0][5]!="source-effective-photon"||rows[0][7]!=d.physicalThresholdSourceIdentity){
      *e="source-effective crosswalk must bind exact channel/states/reference, named energy convention and physical formation-energy source";return false;
    }
    // Parse the NATIVE source number only to compare with the same archived
    // registry representation. It is never a fitted or user-tunable delta.
    std::istringstream number(rows[0][6]);number.imbue(std::locale::classic());
    double threshold=0;number>>std::noskipws>>threshold;
    if(!number||number.peek()!=std::char_traits<char>::eof()||!std::isfinite(threshold)||threshold<=0){
      *e="explicit crosswalk source threshold must be one finite complete eV token";return false;
    }
    Species::Photochemistry::Channel channel;
    if(!source.GetChannel(d.channelId,&channel,e))return false;
    if(!channel.hasThreshold||threshold!=channel.thresholdEnergyEv){
      *e="source-effective crosswalk must retain the exact frozen source threshold; no replacement or arbitrary correction";return false;
    }
  }else{*e="unsupported source crosswalk energy policy";return false;}
  return true;
}
bool EventFrame(const Request& request,const T::Field& field,const Controller::EventContext& event,
    E::EventInput* input,std::string* error){
#if defined(_NO_SPICE_CALLS_)
  *error="actual Na event requires normal furnished SPICE state transforms";return false;
#else
  SpiceInt frameCode=0,centre=0,kind=0,classId=0;SpiceBoolean found=SPICEFALSE;
  namfrm_c(request.inertialFrame.c_str(),&frameCode);frinfo_c(frameCode,&centre,&kind,&classId,&found);
  if(!found||kind!=1){*error="declared SUN-origin kinematic axes must be a real SPICE inertial-class frame";return false;}
  SpiceDouble transform[6][6],moon[6],light=0;
  sxform_c(request.inertialFrame.c_str(),Frames::SolarOrbital,field.frozenGeometryEpochS,transform);
  spkezr_c("MOON",field.frozenGeometryEpochS,request.inertialFrame.c_str(),"NONE","SUN",moon,&light);
  if(failed_c()){*error="normal furnished SPICE cannot supply declared SUN-origin event state transform";reset_c();return false;}
  auto& frame=input->frame;frame.inertialIdentity="SUN-origin:"+request.inertialFrame;frame.eventIdentity=Frames::SolarOrbital;
  frame.transformSourceIdentity="normal furnished CSPICE sxform+spkezr MOON/SUN NONE";
  frame.epochConventionIdentity="frozen LSO coefficient basis; advancing actual material transport ET";
  frame.geometryBasisEpochEtS=field.frozenGeometryEpochS;frame.physicalEventEpochEtS=event.physicalEpochS;
  frame.velocityOffsetSupplied=true;
  E::Vector3 relativeInertial{},heliocentric{};
  for(unsigned i=0;i<3;++i){
    for(unsigned j=0;j<3;++j){frame.inertialToEventRotation[3*i+j]=transform[i][j];relativeInertial[i]+=static_cast<long double>(transform[j][i])*event.parent.positionM[j];}
    heliocentric[i]=1000*static_cast<long double>(moon[i])+relativeInertial[i];
  }
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
    frame.velocityOffsetEventMS[i]+=static_cast<long double>(transform[i+3][j])*relativeInertial[j]-transform[i][j]*(1000*static_cast<long double>(moon[j+3]));
  const auto norm=std::hypot(heliocentric[0],heliocentric[1],heliocentric[2]);
  if(!std::isfinite(norm)||norm<=0){*error="nonpositive actual heliocentric photon range";return false;}
  for(unsigned i=0;i<3;++i)input->incomingPhotonDirectionInertial[i]=heliocentric[i]/norm;
  return true;
#endif
}
}
E::Status Prepared::InspectProfile(const Request& r,P::Registry* output){
  if(!output||r.parentStateId.empty()||r.ionStateId.empty()||r.electronStateId.empty())
    return Error(E::Code::InvalidInput,"profile output and explicit3physical state IDs required");
  try {
    if(!std::filesystem::is_regular_file(r.profileMetadataPath)||!std::filesystem::is_regular_file(r.profileNativeManifestPath))return Error(E::Code::Unavailable,"physical profile metadata/native manifest file unavailable");
    for(const auto* file:{&r.profile.isotopes.path,&r.profile.constants.path,&r.profile.states.path})
      if(!std::filesystem::is_regular_file(*file))return Error(E::Code::Unavailable,"required physical profile CSV unavailable: "+*file);
    for(const auto& file:r.profile.sources)if(!std::filesystem::is_regular_file(file.path))return Error(E::Code::Unavailable,"required source evidence unavailable: "+file.path);
    std::string error;if(!D::VerifyFileSha256(r.profileMetadataPath,r.profileMetadataSha256,&error)||!NativeProfileManifest(r,&error))return Error(E::Code::SourceMismatch,error);
    P::Registry profile;if(!profile.Load(r.profile,&error))return Error(E::Code::SourceMismatch,error);
    for(const auto& id:{r.parentStateId,r.ionStateId,r.electronStateId}){
      Species::ReducedKinematics::ResolvedState state;Species::ReducedKinematics::MassValue rest;
      if(!profile.GetState(id,&state,&error)||!profile.GetPhysicalRestMass(id,&rest,&error))
        return Error(E::Code::Unavailable,"source-bound required state/rest mass unavailable: "+id+"; "+error);
    }
    *output=std::move(profile);return {};
  }catch(const std::bad_alloc&){return Error(E::Code::ResourceLimit,"source profile inspection exhausted memory; provide more memory or reduce concurrently retained profiles; no particle was created");}
  catch(const std::exception&e){return Error(E::Code::InvalidInput,std::string("source profile inspection exception: ")+e.what());}
}
E::Status Prepared::Create(const Species::Photochemistry::Registry& source,const Request& r,std::shared_ptr<Prepared>*out){
  if(!out)return Error(E::Code::InvalidInput,"prepared physical cohort output is required");
  try {
    if(r.parentStateId.empty()||r.ionStateId.empty()||r.electronStateId.empty()||r.parentCompiledSpecies<0||r.ionCompiledSpecies<0||
        r.parentCompiledSpecies==r.ionCompiledSpecies||r.inertialFrame.empty()||r.inertialOrigin!="SUN"||
        !Positive(r.sourceReferenceDistanceM)||r.sourceReferenceDistanceIdentity.empty()||
        r.maximumRootRefinements==0||!Positive(r.exposure.rootTimeToleranceS)||
        !std::isfinite(r.exposure.absoluteDepthTolerance)||r.exposure.absoluteDepthTolerance<0||
        !std::isfinite(r.exposure.relativeDepthTolerance)||r.exposure.relativeDepthTolerance<0||
        (r.exposure.absoluteDepthTolerance==0&&r.exposure.relativeDepthTolerance==0)||
        r.exposure.maximumSegments==0||r.exposure.maximumGeometryIntervals==0)
      return Error(E::Code::InvalidInput,"explicit state IDs/compiled mapping, SUN inertial origin/frame, source reference distance and numerical controls are required");
    if(r.prescribedSun&&(r.earthOccultation||r.prescribedSunIdentity.empty()||
        !std::isfinite(r.prescribedSunPositionM[0])||!std::isfinite(r.prescribedSunPositionM[1])||!std::isfinite(r.prescribedSunPositionM[2])))
      return Error(E::Code::InvalidInput,"a prescribed static Sun needs a finite LSO position, a fixture identity and Earth occultation off");
    auto p=std::shared_ptr<Prepared>(new Prepared);p->request_=r;
    const auto inspected=InspectProfile(r,&p->profile_);if(!inspected.ok())return inspected;
    std::string error;
    auto d=r.controls;E::StateBinding* bindings[]={&d.parent,&d.ion,&d.electron};
    const std::string ids[]={r.parentStateId,r.ionStateId,r.electronStateId};
    const int compiled[]={r.parentCompiledSpecies,r.ionCompiledSpecies,-1};d.catalog.resize(3);
    for(unsigned i=0;i<3;++i){
      if(!p->profile_.GetState(ids[i],&bindings[i]->state,&error)||!p->profile_.GetPhysicalRestMass(ids[i],&bindings[i]->physicalRestMass,&error))
        return Error(E::Code::Unavailable,"source-bound required state/rest mass unavailable: "+ids[i]+"; "+error);
      bindings[i]->catalogIndex=i;double internal=0;
      const auto status=E::ReferenceCatalogEntry(bindings[i]->state,compiled[i],&d.catalog[i],&internal);
      if(!status.ok())return status;
      if(i==0)p->parentInternalEnergyJ_=internal;
    }
    if(!std::filesystem::is_regular_file(r.crosswalkPath))return Error(E::Code::Unavailable,"reviewed Na channel/state/threshold crosswalk is unavailable");
    if(d.photonEnergyPolicy==Species::ReducedKinematics::PhotonEnergyPolicy::SourceEffectivePhoton)
      d.physicalThresholdSourceIdentity=d.ion.state.definition().formation.sourceIdentity;
    if(!Crosswalk(r,source,d,&error))return Error(E::Code::SourceMismatch,error);
    P::Constant charge;
    if(!p->profile_.GetConstant("elementary charge","C",&charge,&error))return Error(E::Code::Unavailable,error);
    p->elementaryChargeC_=static_cast<double>(charge.value);p->elementaryChargeIdentity_=charge.sourceIdentity+":"+charge.sourceRecord;
    // The SI eV conversion and compiled charge use exactly the SAME source.
    if(charge.value!=1.602176634e-19L||!Positive(p->elementaryChargeC_))return Error(E::Code::SourceMismatch,"physical source charge does not equal exact SI elementary charge");
    d.electronVoltJ=charge.value;d.electronVoltSourceIdentity=p->elementaryChargeIdentity_;
    d.sourceProfileIdentity=p->profile_.identity();d.sourceProfileSha256=r.profileMetadataSha256;
    d.stateCrosswalkIdentity=r.crosswalkPath;d.stateCrosswalkSha256=r.crosswalkSha256;d.sourceStateCrosswalkAvailable=true;
    const auto prepared=E::SodiumPhotoionization::Create(source,d,&p->provider_);if(!prepared.ok())return prepared;
    p->definition_=std::move(d);*out=std::move(p);return {};
  }catch(const std::bad_alloc&){return Error(E::Code::ResourceLimit,"physical Na profile preparation exhausted memory; reduce source inventory/cohort metadata or provide more memory; no particle was created");}
  catch(const std::exception&e){return Error(E::Code::InvalidInput,std::string("physical Na preparation exception: ")+e.what());}
}
E::Status Prepared::RegisterForNormalStartup(){
  std::string error;
  if(!provider_.prepared()||installedOwner)
    return Error(E::Code::InvalidInput,"unprepared or already installed normal physical Na binding");
  // Retain the exact context for the process lifetime. Neither the normal
  // callback nor future event hooks may borrow a temporary caller object.
  auto owner=shared_from_this();
  if(!NativeRegistration::SetBeforeParticleBuffer(Install,this,definition_.identity,&error))
    return Error(E::Code::InvalidInput,error);
  installedOwner=std::move(owner);return {};
}
bool Prepared::Install(void* user,std::string* error){
  auto& p=*static_cast<Prepared*>(user);
  if(p.request_.prescribedSun){
    const Species::TrajectoryExposure::Vector sun{{p.request_.prescribedSunPositionM[0],
        p.request_.prescribedSunPositionM[1],p.request_.prescribedSunPositionM[2]}};
    if(!Exposure::FrozenProductionField::CreatePrescribedSunAtEpoch(p.request_.exposure,sun,
        Moon::OrbitalMotion::et,&p.geometry_,error))return false;
    p.exposureIdentity_=p.definition_.identity+";prescribed-static-point-Sun-geometry;"+p.request_.prescribedSunIdentity;
  }else{
    if(!Exposure::FrozenProductionField::CreateAtEpoch(p.request_.exposure,p.request_.earthOccultation,
        Moon::OrbitalMotion::et,&p.geometry_,error))return false;
    p.exposureIdentity_=p.definition_.identity+";normal-frozen-point-Sun-geometry";
  }
  p.field_=p.geometry_->geometry();p.field_.referenceDistanceM=p.request_.sourceReferenceDistanceM;
  // Check actual runtime endpoint/source installation, never assume a table
  // loaded during a different preflight became this process's active family.
  const std::string endpoint=p.definition_.solar==Species::Photochemistry::SolarCase::Quiet?"quiet":"active";
  if(Runtime::GetConfiguration().photochemistrySolarCase!=endpoint){
    *error="normal runtime solar endpoint differs from the prepared source-bound cohort";return false;}
  const auto* source=Runtime::GetPhotochemistryRegistry();
  if(!source||!source->loaded()){*error="normal runtime PHIDRATES registry was not installed";return false;}
  E::SodiumPhotoionization native;auto status=E::SodiumPhotoionization::Create(*source,p.definition_,&native);
  if(!status.ok()){*error=status.message;return false;}p.provider_=native;
  if(!SemanticRandom::RegisterStorage(p.request_.random,error))return false;
  Controller::Configuration c;c.identity=p.definition_.identity;c.eventFrame=Frames::SolarOrbital;c.catalog=p.definition_.catalog;
  c.clockPrecision=p.request_.random.clockPrecision;
  c.resolveInitialState=Resolve;c.openUniform=SemanticRandom::DrawForParticle;c.userData=&p;c.ledger=&p.ledger_;
  c.elementaryChargeC=p.elementaryChargeC_;c.elementaryChargeSourceIdentity=p.elementaryChargeIdentity_;
  c.initializeProduct=SemanticRandom::InitializeStagedProduct;c.productInitializerIdentity="full-genealogy-HMAC-SHA256-child-key-v1";
  c.exposureHistory=History;c.exposureProviderIdentity=p.exposureIdentity_;
  c.exposureControls.maximumRootRefinements=p.request_.maximumRootRefinements;
  Controller::Channel channel;channel.parentCatalogSpecies=p.definition_.parent.catalogIndex;channel.reaction=*p.provider_.reaction();
  channel.scope=Controller::HazardScope::TrajectoryDependent;channel.kinematics=Kinematics;
  channel.hazardEvidenceIdentity=p.definition_.channelId;channel.kinematicProviderIdentity=p.definition_.plannerAuthorization.identity;
  c.channels={channel};
  if(!Controller::Install(c,error))return false;
  if(Runtime::GetConfiguration().ionSurfaceBoundary=="ideal-collector") {
    Species::IonCollector::StateDefinition ion;status=p.provider_.IonCollectorState(&ion);
    if(!status.ok()){*error=status.message;return false;}
    if(!SurfaceEnvironment::InstallProductionIonCollector({ion},CollectorState,&p,&p.collectorLedger_,error,
        p.elementaryChargeC_,p.elementaryChargeIdentity_))return false;
  }
  return true;
}
bool Prepared::CollectorState(long int ptr,void*,void* user,std::size_t* state,std::string* error){
  auto& p=*static_cast<Prepared*>(user);Controller::History history;
  if(!Controller::ReadHistory(ptr,&history,error))return false;
  if(history.catalogSpecies!=p.definition_.ion.catalogIndex){*error="collector state is not the actual staged ion history";return false;}
  *state=0;error->clear();return true;
}
bool Prepared::Resolve(long int ptr,Controller::Node*,void* user,std::size_t* state,double*internal,std::string* error){
  auto& p=*static_cast<Prepared*>(user);SemanticRandom::State key;
  if(!SemanticRandom::Read(ptr,&key,error))return false;
  if(PIC::ParticleBuffer::GetI(ptr)!=p.definition_.catalog[p.definition_.parent.catalogIndex].compiledSpecies){
    *error="only explicit initial neutral Na sources may open a new controller history; products must retain staged state";return false;}
  *state=p.definition_.parent.catalogIndex;*internal=p.parentInternalEnergyJ_;error->clear();return true;
}
bool Prepared::History(const Controller::ExposureRequest& request,void* user,Controller::ExposureHistory*out,std::string*error){
  auto& p=*static_cast<Prepared*>(user);double rate=0;auto status=p.provider_.RateS(1,1,&rate);
  if(!status.ok()){*error=status.message;return false;}T::History history;
  if(!T::BuildPhotonHistory(request.arc,request.maximumS,p.field_,{{0,rate}},{},p.request_.exposure,&history,error))return false;
  Controller::ExposureHistory result;result.segments=std::move(history.segments);
  const long double bound=history.interpolationDepthErrorBound+history.shadowBoundaryDepthErrorBound;
  result.absoluteDepthErrorBound=bound==0?0:std::nextafter(bound,std::numeric_limits<long double>::infinity());
  *out=std::move(result);error->clear();return true;
}
bool Prepared::Kinematics(const Controller::EventContext& event,Species::ReactionExecution::KinematicPolicy*out,std::string* error){
  auto& p=*static_cast<Prepared*>(event.userData);E::EventInput input;input.parentAtEvent=event.parent;input.interval=event.interval;
  double polar=0,azimuth=0,u=0;
  if(!SemanticRandom::DrawForParticle(Controller::DrawPurpose::Kinematics,event.particle,event.channelIndex,nullptr,&polar,error)||
      !SemanticRandom::DrawForParticle(Controller::DrawPurpose::Kinematics,event.particle,event.channelIndex,nullptr,&azimuth,error)||
      !SemanticRandom::ProspectiveProductClockForParticle(event.particle,0,&u,error)||
      !Species::ReactionExecution::OpticalDepthClock::FromOpenUniform(u,&input.freshIonClock,error,p.request_.random.clockPrecision))return false;
  input.polarOpenUniform=polar;input.azimuthOpenUniform=azimuth;
  input.ionReferenceWeight=event.node->block->GetLocalParticleWeight(p.definition_.catalog[p.definition_.ion.catalogIndex].compiledSpecies);
  if(!EventFrame(p.request_,p.field_,event,&input,error))return false;
  // The accepted channel is already the positive LEFT-limit photon event.
  // A static re-query at an exact light->dark boundary could see the RIGHT
  // value and reject that SAME event. Preserve the accepted source convention
  // and diagnose its affine rate separately from exact distance-only RateS.
  if(event.selectedInstantaneousRateS<=0||!std::isfinite(event.selectedInstantaneousRateS)||
      event.channelIndex!=0||event.acceptedExposureIdentity!=p.exposureIdentity_||event.rateConventionIdentity!="declared-piecewise-affine-first-hit-LEFT-limit"){
    *error="physical photon composition requires the preserved positive accepted LEFT-limit source decision";return false;}
  input.pointSunVisible=1;
  const double dx=p.field_.sunPositionM[0]-event.parent.positionM[0],dy=p.field_.sunPositionM[1]-event.parent.positionM[1],dz=p.field_.sunPositionM[2]-event.parent.positionM[2];
  input.distanceAu=std::hypot(dx,dy,dz)/p.field_.referenceDistanceM;
  input.photonGeometryIdentity=(p.request_.prescribedSun?"prescribed static point-Sun ("+p.request_.prescribedSunIdentity+")+actual sphere/LOLA BVH; ":
      std::string("normal SPICE+actual sphere/LOLA BVH; point-Sun; "))+p.request_.sourceReferenceDistanceIdentity;
  E::Composition composition;auto status=p.provider_.Compose(input,&composition);if(!status.ok()){*error=status.message;return false;}
  // Reserve receipt space before returning a prospective policy. This is not
  // yet event publication; callers may only publish entries paired to actual
  // committed reaction ledger counts/identities after successful completion.
  ComposedEvent entry;if(!SemanticRandom::Read(event.particle,&entry.key,error))return false;
  entry.parentAtEvent=event.parent;entry.ionAtEvent=composition.plan.products().at(0);
  entry.receipt=std::move(composition.receipt);
  entry.acceptedInstantaneousRateS=event.selectedInstantaneousRateS;
  entry.acceptedExposureIdentity=event.acceptedExposureIdentity;
  entry.rateConventionIdentity=event.rateConventionIdentity;
  const auto* clock=GetCommonTransportClock();
  if(!clock||!clock->step()){*error="physical qualification event requires an active normal common-time step";return false;}
  entry.relativeEventTimeS=clock->step()->relativeBeginS+static_cast<long double>(event.interval.elapsedS());
  p.composed_.push_back(std::move(entry));*out=std::move(composition.policy);error->clear();return true;
}
bool Prepared::InitializeSourceKey(long int particle,std::uint64_t id,std::string* error){
  return SemanticRandom::InitializeSource(particle,id,error);
}
}}}
