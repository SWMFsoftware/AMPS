#include "MoonReactionController.h"
#include "Moon.h"
#include "MoonInput.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <memory>
#include <set>
#include <type_traits>

namespace Moon { namespace Chemistry { namespace Controller {
namespace {
constexpr std::uint64_t HistoryMagic=UINT64_C(0x4d4f4f4e52584e31);
struct Record {
  std::uint64_t magic=0,fingerprint=0;
  History history;
  // Immutable origin of the CURRENT core trace (source, daughter birth or
  // subsequent outer interval). MPI copies these values with the history.
  // Reconstructing exposure from this origin avoids summing rounded adaptive
  // dt values while the mover closes its interval by subtracting remaining dt.
  R::OpticalDepthClock originClock;
  double originRemainingS=0,originElapsedS=0,originEpochS=0;
};
static_assert(std::is_trivially_copyable<Record>::value,
    "controller MPI particle records must contain no owning pointers");
struct Installed {
  Configuration configuration;
  std::vector<std::vector<R::Hazard>> hazards;
  std::uint64_t fingerprint=0;
  std::vector<bool> trajectoryStates;
};
std::unique_ptr<const Installed> Active;
long int Offset=-1;
// Disposable numerical scratch, not a physical history or RNG cache. The
// MPI-only controller accepts one serial core flight at a time. A retry may
// replace this history; only Advance can consume its accepted prefix. Particle
// bytes/clock/ledger never change in Query, and scratch never migrates.
struct ExposureScratch {
  bool valid=false;long int particle=-1;
  double epochS=0,maximumS=0;
  R::OpticalDepthClock clock;
  long double absoluteDepthErrorBound=0;
  std::vector<Species::PiecewiseHazard::Segment> history;
};
thread_local ExposureScratch ExposurePreview;

bool Reject(std::string* error,const std::string& message) {
  if(error)*error=message;
  return false;
}
[[noreturn]] void Fatal(const std::string& message) {
  const std::string text="Moon conservative reaction controller ERROR: "+message;
  exit(__LINE__,__FILE__,text.c_str());std::abort();
}
[[noreturn]] void ResourceFatal() {
  exit(__LINE__,__FILE__,"Moon reaction controller exhausted C++ staging memory before a complete physical event could be published. The run is terminated rather than dropping products, changing weights/rates or redrawing an event. Increase available memory or reduce initial model-particle count / particle-weight resolution, then restart from a retained qualified state.");
  std::abort();
}
void HashBytes(std::uint64_t& hash,const void* bytes,std::size_t count) {
  const auto* data=static_cast<const unsigned char*>(bytes);
  for(std::size_t i=0;i<count;++i){hash^=data[i];hash*=UINT64_C(1099511628211);}
}
void HashString(std::uint64_t& hash,const std::string& value) {
  HashBytes(hash,value.data(),value.size());const unsigned char stop=0;HashBytes(hash,&stop,1);
}
template<class T>void HashValue(std::uint64_t& hash,const T& value){HashBytes(hash,&value,sizeof(value));}
std::uint64_t Fingerprint(const Configuration& c) {
  // Integrity of migrated state, not a replacement for source/data SHA-256
  // provenance. Function/user-data addresses intentionally are NOT serialized.
  std::uint64_t hash=UINT64_C(1469598103934665603);HashString(hash,c.identity);HashString(hash,c.eventFrame);
  HashValue(hash,c.elementaryChargeC);HashString(hash,c.elementaryChargeSourceIdentity);
  HashString(hash,c.productInitializerIdentity);
  const bool hasProductInitializer=c.initializeProduct!=nullptr;HashValue(hash,hasProductInitializer);
  HashString(hash,c.exposureProviderIdentity);HashValue(hash,c.exposureControls.maximumRootRefinements);
  HashValue(hash,c.catalog.size());HashValue(hash,c.channels.size());
  for(const auto& species:c.catalog) {
    HashString(hash,species.identity);HashString(hash,species.state);HashValue(hash,species.compiledSpecies);
    HashValue(hash,species.elementaryCharge);HashValue(hash,species.modelMassKg);
    HashValue(hash,species.nuclei.size());
    for(const auto& nucleus:species.nuclei){HashString(hash,nucleus.identity);HashValue(hash,nucleus.count);}
  }
  for(const auto& channel:c.channels) {
    HashValue(hash,channel.parentCatalogSpecies);HashValue(hash,channel.constantRateS);
    HashValue(hash,channel.scope);
    HashString(hash,channel.reaction.identity);HashString(hash,channel.hazardEvidenceIdentity);HashString(hash,channel.kinematicProviderIdentity);
    for(const auto* terms:{&channel.reaction.reactants,&channel.reaction.products}) {
      const std::size_t size=terms->size();HashValue(hash,size);
      for(const auto& term:*terms){HashValue(hash,term.species);HashValue(hash,term.multiplicity);HashValue(hash,term.population);}
    }
  }return hash;
}
Record Read(long int particle) {
  Record record;std::memcpy(&record,PIC::ParticleBuffer::GetParticleDataPointer(particle)+Offset,sizeof(record));return record;
}
void Write(long int particle,const Record& record) {
  std::memcpy(PIC::ParticleBuffer::GetParticleDataPointer(particle)+Offset,&record,sizeof(record));
}
bool Valid(long int particle,const Record& record,std::string* error) {
  if(!Active || Offset<0 || record.magic!=HistoryMagic || record.fingerprint!=Active->fingerprint ||
      record.history.configurationFingerprint!=Active->fingerprint ||
      record.history.catalogSpecies>=Active->configuration.catalog.size() ||
      !record.history.interval.valid() || !record.history.clock.valid() ||
      !record.originClock.valid() || (record.originClock.consumed() && !record.history.due) ||
      !std::isfinite(record.originRemainingS) || record.originRemainingS<0 ||
      !std::isfinite(record.originElapsedS) || record.originElapsedS<0 ||
      record.originElapsedS>record.history.interval.durationS() || !std::isfinite(record.originEpochS) ||
      record.history.clock.consumed()!=record.history.due || !std::isfinite(record.history.internalEnergyJ))
    return Reject(error,"invalid/corrupt controller history, interval, optical depth or immutable registry fingerprint");
  if(!std::isfinite(record.history.exposureApproximationErrorBound)||record.history.exposureApproximationErrorBound<0||
      !std::isfinite(record.history.clockStorageErrorBound)||record.history.clockStorageErrorBound<0)
    return Reject(error,"invalid accumulated trajectory/clock numerical error diagnostics");
  if(static_cast<int>(PIC::ParticleBuffer::GetI(particle))!=Active->configuration.catalog[record.history.catalogSpecies].compiledSpecies)
    return Reject(error,"particle compiled species changed without an explicit controller state transition");
  if((record.history.due&&(!std::isfinite(record.history.selectedInstantaneousRateS)||record.history.selectedInstantaneousRateS<=0))||
      (!record.history.due&&record.history.selectedInstantaneousRateS!=0))
    return Reject(error,"accepted event lacks its preserved positive instantaneous source rate");
  if(record.history.due && (record.history.selectedChannel>=Active->configuration.channels.size() ||
      Active->configuration.channels[record.history.selectedChannel].parentCatalogSpecies!=record.history.catalogSpecies))
    return Reject(error,"pending selected channel does not belong to the preserved parent catalog state");
  return true;
}
double ClockTolerance(double a,double b) {
  // The application declares the physical time tolerance. The extra term is
  // binary64 epoch representation, not a chemistry/acceptance adjustment.
  return std::max(PIC::Mover::ForceCutSurface::Settings.TimeToleranceS,
      32*std::numeric_limits<double>::epsilon()*std::max({1.0,std::abs(a),std::abs(b)}));
}
void ValidateEpoch(long int particle,const Record& record) {
  const double actual=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
  const double preserved=record.history.interval.startS()+record.history.interval.elapsedS();
  if(!std::isfinite(actual) || std::abs(actual-preserved)>ClockTolerance(actual,preserved))
    Fatal("controller and mover physical epochs disagree; do not reset a continued clock or substitute a species/global timestep");
}
double Draw(DrawPurpose purpose,long int particle,std::size_t channel) {
  double uniform=0;std::string error;
  if(!Active->configuration.openUniform(purpose,particle,channel,Active->configuration.userData,&uniform,&error) ||
      !std::isfinite(uniform) || !(uniform>0 && uniform<1))
    Fatal("explicit RNG provider failed or returned a closed/nonfinite uniform: "+error);
  return uniform;
}
Record Checked(long int particle) {
  if(!Active || Offset<0)Fatal("trajectory controller callback invoked before explicit installation");
  Record record=Read(particle);std::string error;
  if(!Valid(particle,record,&error))Fatal(error);
  return record;
}
long double Sum3(long double a,long double b,long double c) {
  // Neumaier compensation is essential here: a tiny fresh clock must survive
  // e.g. delay - 1 + 1, and must never become a premature zero-time event.
  long double sum=0,correction=0;
  for(long double value:{a,b,c}) {
    const long double next=sum+value;
    correction+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;
    sum=next;
  }
  return sum+correction;
}
void SetOrigin(Record& record,double actualRemaining,double actualEpoch) {
  record.originClock=record.history.clock;record.originRemainingS=actualRemaining;
  record.originElapsedS=record.history.interval.elapsedS();record.originEpochS=actualEpoch;
}
double MoverRemaining(long int particle,const Record& record) {
  const auto time=PIC::Mover::ForceCutSurface::GetTimeContext(particle);
  // A source may be inspected before its first call into Trace; the actual
  // mover has not initialized its interval yet. This read-only case uses the
  // explicit prepared origin, not a local-species/global timestep.
  if(time.RemainingTimeS==0 &&
      (time.RequestedTimeS==0 || time.InitialEpochS!=record.originEpochS))
    return record.originRemainingS;
  if(time.RequestedTimeS!=record.originRemainingS ||
      !std::isfinite(time.RemainingTimeS) || time.RemainingTimeS<0 ||
      time.RemainingTimeS>record.originRemainingS || time.InitialEpochS!=record.originEpochS)
    Fatal("mover/controller immutable trace origin or remaining coordinate disagrees");
  return time.RemainingTimeS;
}
double AnchoredDelay(const Record& record,double remaining,long double* originDelay=nullptr) {
  long double delay=0;std::string error;
  if(!R::ConstantHazardEventDelay(Active->hazards[record.history.catalogSpecies],record.originClock,&delay,&error))Fatal(error);
  if(originDelay)*originDelay=delay;
  // No tolerance: even the next representable delay beyond the immutable
  // caller endpoint must remain pending for the next physical interval.
  if(delay>record.originRemainingS)return std::numeric_limits<double>::infinity();
  const long double local=delay==record.originRemainingS?static_cast<long double>(remaining):
      Sum3(delay,-static_cast<long double>(record.originRemainingS),remaining);
  if(!(local>0) || !std::isfinite(local))
    Fatal("positive constant-hazard residual time is unrepresentable; refine the physical-time coordinate rather than forcing a zero-time reaction");
  double represented=static_cast<double>(local);
  if(!(represented>0) || !std::isfinite(represented))Fatal("positive event residual does not fit binary64 physical-time storage");
  if(static_cast<long double>(represented)<local)
    represented=std::nextafter(represented,std::numeric_limits<double>::infinity());
  return represented;
}
bool Variable(const Record& record){return Active->trajectoryStates[record.history.catalogSpecies];}
void PreviewExposure(long int particle,Node* node,const Record& record,double maximum) {
  ExposurePreview.valid=false;
  ExposureRequest request;request.particle=particle;request.node=node;
  request.parentCatalogSpecies=record.history.catalogSpecies;
  request.compiledSpecies=PIC::ParticleBuffer::GetI(particle);
  request.physicalEpochS=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
  request.maximumS=maximum;request.time=PIC::Mover::ForceCutSurface::GetTimeContext(particle);
  PIC::ParticleBuffer::GetX(request.arc.positionM.data(),particle);
  PIC::ParticleBuffer::GetV(request.arc.velocityMS.data(),particle);
  double* x=request.arc.positionM.data();double* v=request.arc.velocityMS.data();
  if(PIC::Mover::ForceCutSurface::Acceleration)
    PIC::Mover::ForceCutSurface::Acceleration(request.arc.accelerationMS2.data(),
      request.compiledSpecies,particle,x,v,node,request.time);
  else {
    using PIC::Mover::TotalParticleAcceleration_default;
    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(request.arc.accelerationMS2.data(),
      request.compiledSpecies,particle,x,v,node);
  }
  ExposureHistory supplied;std::string error;
  if(!Active->configuration.exposureHistory(request,Active->configuration.userData,&supplied,&error) ||
      !Species::PiecewiseHazard::ValidateHistory(supplied.segments,&error))Fatal("actual trajectory exposure unavailable: "+error);
  if(!std::isfinite(supplied.absoluteDepthErrorBound)||supplied.absoluteDepthErrorBound<0)
    Fatal("trajectory provider omitted/corrupted its absolute numerical depth error bound");
  auto& history=supplied.segments;
  if(history.empty() || history.front().beginS!=0 || history.back().endS!=maximum)
    Fatal("trajectory provider must cover exactly the complete requested physical prefix");
  std::set<std::size_t> expected;
  for(const auto& h:Active->hazards[record.history.catalogSpecies])expected.insert(h.reaction);
  for(const auto& segment:history){std::set<std::size_t> actual;for(const auto& h:segment.channels)actual.insert(h.reaction);
    if(actual!=expected)Fatal("trajectory exposure omitted/added a registered competing channel");}
  ExposurePreview.particle=particle;ExposurePreview.epochS=request.physicalEpochS;
  ExposurePreview.maximumS=maximum;
  ExposurePreview.clock=record.history.clock;
  ExposurePreview.absoluteDepthErrorBound=supplied.absoluteDepthErrorBound;
  ExposurePreview.history=std::move(history);ExposurePreview.valid=true;
}
void PrepareProduct(long int particle,double actualRemaining) {
  Record record;std::string error;record.magic=HistoryMagic;record.fingerprint=Active->fingerprint;
  auto& state=record.history;
  if(!ReadReactionState(particle,&state.interval,&state.clock,&state.catalogSpecies,&state.internalEnergyJ,&error))
    Fatal("committed product lacks the actual adapter's stored state/clock: "+error);
  if(!state.clock.fresh() || state.clock.consumed())Fatal("a daughter must start with its supplied fresh optical depth, not the parent's consumed clock");
  state.configurationFingerprint=Active->fingerprint;
  SetOrigin(record,actualRemaining,PIC::Mover::ForceCutSurface::GetParticleEpochS(particle));
  if(!Valid(particle,record,&error))Fatal(error);
  ValidateEpoch(particle,record);Write(particle,record);
}
}

bool Enabled(){return static_cast<bool>(Active);}
bool Install(const Configuration& c,std::string* error) {
  if(!error)return false;
  if(Active || PIC::ParticleBuffer::ParticleDataBuffer!=nullptr)
    return Reject(error,"install the controller exactly once, before allocating the particle buffer");
#if _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_
  return Reject(error,"constant-history reaction controller requires MPI-only/OMP1 execution");
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  return Reject(error,"conservative product publication requires individual particle weights ON");
#elif _PIC_ELECTRON_IMPACT_IONIZATION_REACTION_MODE_ != _PIC_MODE_OFF_
  return Reject(error,"generic PIC electron execution must remain OFF; it cannot also remove/create controller particles");
#else
  if(!std::isfinite(c.elementaryChargeC) || c.elementaryChargeC<=0 || c.elementaryChargeSourceIdentity.empty() ||
      ((c.initializeProduct!=nullptr)!=(!c.productInitializerIdentity.empty())))
    return Reject(error,"controller needs explicit positive source-bound charge and matching optional product initializer identity");
  if(!Runtime::HasConfiguration() || Runtime::GetConfiguration().particleTracingMover=="legacy" ||
      !PIC::Mover::ForceCutSurface::IsStorageRegistered())
    return Reject(error,"controller requires an explicitly configured/registered first- or second-order Moon physical-time tracer");
  if(PIC::Mover::QueryTrajectoryEvent || PIC::Mover::AdvanceTrajectoryEvent || PIC::Mover::ProcessTrajectoryEvent)
    return Reject(error,"controller refuses to overwrite an existing or partial trajectory-event hook set");
  if(PIC::ParticleSplitting::Mode!=PIC::ParticleSplitting::_disactivated)
    return Reject(error,"active particle splitting needs a separately qualified physical-history/RNG clone policy; it cannot implicitly duplicate pending clocks");
  if(c.identity.empty() || c.eventFrame!=Moon::Frames::SolarOrbital || !c.resolveInitialState || !c.openUniform || !c.ledger)
    return Reject(error,"controller requires explicit immutable identity, LSO event frame, source-state resolver, open-uniform provider and physical ledger");
  if(!R::ValidateCatalog(c.catalog,error))return false;
  std::unique_ptr<Installed> candidate(new Installed);candidate->configuration=c;
  candidate->hazards.resize(c.catalog.size());std::set<std::string> identities;
  candidate->trajectoryStates.resize(c.catalog.size(),false);
  for(std::size_t index=0;index<c.channels.size();++index) {
    const auto& channel=c.channels[index];
    if(channel.scope!=HazardScope::ConstantHistory && channel.scope!=HazardScope::TrajectoryDependent)
      return Reject(error,"unknown hazard scope");
    if(channel.scope==HazardScope::TrajectoryDependent && (!c.exposureHistory ||
        c.exposureProviderIdentity.empty() || c.exposureControls.maximumRootRefinements==0))
      return Reject(error,"trajectory/source/illumination-dependent hazards require an explicit geometry/exposure service and numerical controls; never freeze a photorate along a moving trajectory");
    if(channel.parentCatalogSpecies>=c.catalog.size() || !std::isfinite(channel.constantRateS) || channel.constantRateS<0 ||
        channel.hazardEvidenceIdentity.empty() || channel.kinematicProviderIdentity.empty() || !channel.kinematics ||
        !identities.insert(channel.reaction.identity).second || !R::ValidateReaction(c.catalog,channel.reaction,error))
      return Reject(error,"invalid/missing/duplicate constant channel, topology, source evidence or explicit kinetic-policy provider");
    for(const auto& term:channel.reaction.reactants)
      if(term.population==R::Population::Tracked && term.species!=channel.parentCatalogSpecies)
        return Reject(error,"registered channel parent does not match its explicit tracked reaction state");
    candidate->hazards[channel.parentCatalogSpecies].push_back({index,channel.constantRateS});
    if(channel.scope==HazardScope::TrajectoryDependent)candidate->trajectoryStates[channel.parentCatalogSpecies]=true;
  }
  for(const auto& hazards:candidate->hazards) {
    double probability=0;if(!R::EventProbability(hazards,0,&probability,error))return false;
  }
  for(const auto& species:c.catalog)if(species.compiledSpecies>=0) {
    const int index=species.compiledSpecies;
    if(index>=PIC::nTotalSpecies || PIC::MolecularData::GetMass(index)!=species.modelMassKg ||
        PIC::MolecularData::GetElectricCharge(index)!=species.elementaryCharge*c.elementaryChargeC)
      return Reject(error,"controller catalog mass/charge/index does not match actual generated PIC species data");
  }
  candidate->fingerprint=Fingerprint(c);
  RegisterReactionStorage();
  PIC::ParticleBuffer::RequestDataStorage(Offset,sizeof(Record),alignof(Record));
  Active=std::move(candidate);
  PIC::Mover::QueryTrajectoryEvent=Query;PIC::Mover::AdvanceTrajectoryEvent=Advance;PIC::Mover::ProcessTrajectoryEvent=Process;
  error->clear();return true;
#endif
}
void ResetForNewParticle(PIC::ParticleBuffer::byte* data) {
  if(Offset>=0 && data)std::memset(data+Offset,0,sizeof(Record));
}
bool ReadHistory(long int particle,History* output,std::string* error) {
  if(!output || !error || !Active || particle<0 || particle>=PIC::ParticleBuffer::GetMaxNPart() ||
      !PIC::ParticleBuffer::IsParticleAllocated(particle))return Reject(error,"controller history needs an active controller and allocated particle");
  const Record record=Read(particle);if(!Valid(particle,record,error))return false;
  *output=record.history;error->clear();return true;
}
void PrepareParticle(long int particle,double requested,Node* node) {
  if(!Active)return;
  try {
    if(!node || !node->block || node->Thread!=PIC::ThisThread || !std::isfinite(requested) || requested<0)
      Fatal("invalid actual owner/block or requested physical interval");
    if(PIC::ParticleSplitting::Mode!=PIC::ParticleSplitting::_disactivated)
      Fatal("particle splitting was activated after controller installation without a qualified history-clone policy");
    const double epoch=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
    if(!std::isfinite(epoch))Fatal("Moon must initialize a finite particle physical epoch before controller preparation");
    Record record=Read(particle);std::string error;
    if(record.magic==0) {
      std::size_t state=0;double internal=0;
      if(!Active->configuration.resolveInitialState(particle,node,Active->configuration.userData,&state,&internal,&error) ||
          state>=Active->configuration.catalog.size() || !std::isfinite(internal) ||
          Active->configuration.catalog[state].compiledSpecies!=static_cast<int>(PIC::ParticleBuffer::GetI(particle)))
        Fatal("explicit source-state resolver failed; do not guess a ground/excited state from a compiled species: "+error);
      record.magic=HistoryMagic;record.fingerprint=Active->fingerprint;
      record.history.catalogSpecies=state;record.history.internalEnergyJ=internal;
      record.history.configurationFingerprint=Active->fingerprint;
      if(!R::TraceInterval::Create(epoch,requested,0,&record.history.interval,&error) ||
          !R::OpticalDepthClock::FromOpenUniform(Draw(DrawPurpose::SourceClock,particle,std::numeric_limits<std::size_t>::max()),&record.history.clock,&error))Fatal(error);
      SetOrigin(record,requested,epoch);
    }
    else {
      if(!Valid(particle,record,&error))Fatal(error);
      ValidateEpoch(particle,record);
      const double remaining=record.history.interval.remainingS();
      if(remaining==0) {
        // A completed caller interval is not a new physical species history.
        // Retain pending optical depth/selected endpoint event while opening
        // precisely the next interval supplied by the actual Moon walker.
        if(!R::TraceInterval::Create(epoch,requested,0,&record.history.interval,&error))Fatal(error);
        SetOrigin(record,requested,epoch);
      }
      else if(std::abs(remaining-requested)>ClockTolerance(remaining,requested))
        Fatal("a continued controller history must retain its remaining physical seconds; it cannot restart with a product/species timestep");
    }
    Write(particle,record);
  }catch(const std::bad_alloc&){ResourceFatal();}
  catch(const std::exception& error){Fatal(error.what());}
}
int Query(long int particle,int species,Node* node,double,double maximum,double& eventTime) {
  try {
    const Record record=Checked(particle);ValidateEpoch(particle,record);
    if(species!=static_cast<int>(PIC::ParticleBuffer::GetI(particle)) || !std::isfinite(maximum) || maximum<0)
      Fatal("invalid pure trajectory-event preview inputs");
    // No Write(), callback draw or kinetic-policy call occurs in this path.
    // An accepted crossing cached on a block-face tie survives a zero-time
    // ownership hand-off; the receiving owner processes that SAME decision.
    if(record.history.due){eventTime=0;return PIC::Mover::TrajectoryEventPending;}
    if(Variable(record)) {
      if(maximum==0)return PIC::Mover::TrajectoryEventNone;
      PreviewExposure(particle,node,record,maximum);
      Species::ReactionExposure::Preview preview;std::string error;
      const auto status=Species::ReactionExposure::Locate(ExposurePreview.history,0,maximum,
        record.history.clock,Active->configuration.exposureControls,&preview,&error);
      if(status!=Species::PiecewiseHazard::Status::Complete && status!=Species::PiecewiseHazard::Status::Event)
        Fatal("trajectory optical-clock preview refused: "+error);
      if(preview.event){const long double exact=preview.location.eventTimeS;
        eventTime=static_cast<double>(exact);if(eventTime<exact)eventTime=std::nextafter(eventTime,HUGE_VAL);
        if(!(eventTime>0)||eventTime>maximum)Fatal("trajectory event time cannot be represented within its physical prefix");
        return PIC::Mover::TrajectoryEventPending;}
      return PIC::Mover::TrajectoryEventNone;
    }
    const double delay=AnchoredDelay(record,MoverRemaining(particle,record));
    if(delay<=maximum){eventTime=delay;return PIC::Mover::TrajectoryEventPending;}
    return PIC::Mover::TrajectoryEventNone;
  }catch(const std::bad_alloc&){ResourceFatal();}
  catch(const std::exception& error){Fatal(error.what());}
}
void Advance(long int particle,int species,Node*,double,double elapsed) {
  try {
    Record record=Checked(particle);ValidateEpoch(particle,record);
    if(species!=static_cast<int>(PIC::ParticleBuffer::GetI(particle)) || !std::isfinite(elapsed) || !(elapsed>0) || record.history.due)
      Fatal("accepted flight cannot advance an already due/corrupt history");
    const auto& hazards=Active->hazards[record.history.catalogSpecies];std::string error;
    const double moverRemaining=MoverRemaining(particle,record);long double originDelay=0;
    const bool variable=Variable(record);
    const double delay=variable?HUGE_VAL:AnchoredDelay(record,moverRemaining,&originDelay);
    if(elapsed>moverRemaining)Fatal("accepted flight exceeds the actual mover remaining seconds");
    // At the exact caller endpoint the physical elapsed coordinate is the
    // immutable duration, byte-for-byte. Otherwise use compensated addition;
    // subtracting (remaining-dt) first would erase tiny initial flights.
    long double originElapsed=elapsed==moverRemaining?record.originRemainingS:
      Sum3(record.originRemainingS,-static_cast<long double>(moverRemaining),elapsed);
    // The core invokes Advance only after accepting this actual flight.
    // Its own physical epoch has not yet been incremented, so ValidateEpoch
    // above compares two PRE-segment clocks. Store the new exposure/channel
    // before the core can copy particle bytes to an MPI owner/continuation.
    bool crossing=delay<=elapsed;
    R::HazardClockAdvance outcome;long double selectedRate=0;
    if(variable) {
      if(!ExposurePreview.valid || ExposurePreview.particle!=particle ||
          ExposurePreview.epochS!=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle) ||
          !ExposurePreview.clock.SamePendingDepth(record.history.clock) || elapsed>ExposurePreview.maximumS)
        Fatal("accepted trajectory flight has no matching pure exposure preview; never recompute a different hazard after accepting a flight");
      Species::ReactionExposure::Preview preview;
      auto status=Species::ReactionExposure::Locate(ExposurePreview.history,0,elapsed,
        record.history.clock,Active->configuration.exposureControls,&preview,&error);
      if(status!=Species::PiecewiseHazard::Status::Complete&&status!=Species::PiecewiseHazard::Status::Event)Fatal(error);
      crossing=preview.event;
      const double uniform=crossing?Draw(DrawPurpose::Channel,particle,std::numeric_limits<std::size_t>::max()):NAN;
      Species::ReactionExposure::Decision decision;
      status=Species::ReactionExposure::Evaluate(ExposurePreview.history,0,elapsed,record.history.clock,
        Active->configuration.exposureControls,uniform,&decision,&error);
      if(status!=Species::PiecewiseHazard::Status::Complete&&status!=Species::PiecewiseHazard::Status::Event)Fatal(error);
      outcome.clock=decision.advance.clock;outcome.event=decision.advance.event;outcome.reaction=decision.advance.reaction;
      if(outcome.event) {
        for(const auto& rate:decision.preview.location.instantaneousRates)
          if(rate.reaction==outcome.reaction)selectedRate=rate.rateS;
        if(!std::isfinite(selectedRate)||selectedRate<=0)
          Fatal("selected trajectory channel has no positive accepted left-limit rate");
      }
      auto accumulate=[](long double before,long double increment){if(increment==0)return before;
        const long double bound=std::nextafter(before+increment,std::numeric_limits<long double>::infinity());
        if(!std::isfinite(bound)||bound<0)Fatal("accumulated exposure-error bound overflowed");
        return bound;};
      record.history.exposureApproximationErrorBound=accumulate(record.history.exposureApproximationErrorBound,ExposurePreview.absoluteDepthErrorBound);
      record.history.clockStorageErrorBound=accumulate(record.history.clockStorageErrorBound,decision.advance.residualDepthRoundingErrorBound);
      ExposurePreview.valid=false;
    }
    // A ceil-rounded local event delay can put its reconstructed total offset
    // on a long-double multiplication boundary. Align to the SAME analytically
    // located root only after the accepted local flight has proven crossing.
    if(crossing && originElapsed<originDelay)originElapsed=originDelay;
    if(!variable) {
      const double uniform=crossing?Draw(DrawPurpose::Channel,particle,std::numeric_limits<std::size_t>::max()):NAN;
      if(!R::AdvanceConstantHazardClock(hazards,originElapsed,record.originClock,uniform,&outcome,&error))Fatal(error);
      if(outcome.event)for(const auto& hazard:hazards)if(hazard.reaction==outcome.reaction)selectedRate=hazard.rateS;
    }
    if(outcome.event!=crossing)Fatal("pure compensated delay and accepted combined exposure disagree at representable event time");
    const double remaining=record.history.interval.remainingS();
    if(elapsed>remaining+ClockTolerance(elapsed,remaining))Fatal("accepted physical flight overshoots the preserved controller interval");
    // Only representational endpoint alignment is allowed here. Exposure is
    // still debited for the full accepted elapsed seconds; no rate/time is
    // shortened to make a reaction or model residual agree.
    R::TraceInterval advanced;
    const double totalElapsed=elapsed==moverRemaining?record.history.interval.durationS():
      std::min(record.history.interval.durationS(),static_cast<double>(static_cast<long double>(record.originElapsedS)+originElapsed));
    if(!R::TraceInterval::Create(record.history.interval.startS(),record.history.interval.durationS(),totalElapsed,&advanced,&error))Fatal(error);
    record.history.interval=advanced;record.history.clock=outcome.clock;
    record.history.due=outcome.event;record.history.selectedInstantaneousRateS=outcome.event?selectedRate:0;
    if(outcome.event)record.history.selectedChannel=outcome.reaction;
    Write(particle,record);
  }catch(const std::bad_alloc&){ResourceFatal();}
  catch(const std::exception& error){Fatal(error.what());}
}
int Process(long int particle,double* position,double* velocity,int& species,Node* node,double,double& remaining) {
  try {
    const Record record=Checked(particle);ValidateEpoch(particle,record);
    if(!record.history.due || !node || !node->block || node->Thread!=PIC::ThisThread ||
        !std::isfinite(remaining) || remaining<0 ||
        std::abs(remaining-record.history.interval.remainingS())>ClockTolerance(remaining,record.history.interval.remainingS()))
      Fatal("event callback has no preserved due channel or consistent actual physical interval/owner");
    // Position/velocity and time now belong to the accepted physical event,
    // not a preview/midpoint or a species-dependent nominal step. The cached
    // channel is reused; this callback never redraws which source occurred.
    const auto& channel=Active->configuration.channels[record.history.selectedChannel];
    EventContext context;context.particle=particle;context.node=node;context.channelIndex=record.history.selectedChannel;
    context.catalog=&Active->configuration.catalog;context.reaction=&channel.reaction;context.interval=record.history.interval;
    context.physicalEpochS=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
    context.selectedInstantaneousRateS=record.history.selectedInstantaneousRateS;
    context.acceptedExposureIdentity=Variable(record)?Active->configuration.exposureProviderIdentity:channel.hazardEvidenceIdentity;
    context.rateConventionIdentity=Variable(record)?"declared-piecewise-affine-first-hit-LEFT-limit":"declared-constant-history";
    context.exposureApproximationErrorBound=record.history.exposureApproximationErrorBound;
    context.clockStorageErrorBound=record.history.clockStorageErrorBound;
    context.openUniform=Active->configuration.openUniform;context.userData=Active->configuration.userData;
    context.parent.species=record.history.catalogSpecies;context.parent.referenceWeight=node->block->GetLocalParticleWeight(species);
    context.parent.individualCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
    context.parent.internalEnergyJ=record.history.internalEnergyJ;
    for(int k=0;k<3;++k){context.parent.positionM[k]=position[k];context.parent.velocityMS[k]=velocity[k];}
    R::KinematicPolicy policy;std::string error;
    if(!channel.kinematics(context,&policy,&error) || policy.identity!=channel.kinematicProviderIdentity ||
        policy.eventFrame!=Active->configuration.eventFrame)
      Fatal("selected channel's actual-event kinetic policy is unavailable/mismatched: "+error);
    // The controller owns registered topology, actual parent and interval.
    // Provider kinematics alone cannot substitute an easier reaction/channel.
    R::EventPlan plan;
    if(!R::MakeEventPlan(Active->configuration.catalog,channel.reaction,context.parent,context.interval,
        R::PlanScope::Transport,&policy,&plan,&error))Fatal("actual-event conservative plan refused: "+error);
    const auto time=PIC::Mover::ForceCutSurface::GetTimeContext(particle);
    if(time.IntegrationOrder!=1 && time.IntegrationOrder!=2)Fatal("event has no preserved first/second physical-time order");
    long int unusedHead=-1;
    ParticleTransaction backend(particle,unusedHead,node,Active->configuration.catalog,Active->configuration.ledger,
        ProductPlacement::DetachedContinuation,true,Active->configuration.elementaryChargeC,
        Active->configuration.initializeProduct,Active->configuration.productInitializerUser);
    const auto ledgerBefore=Active->configuration.ledger->events.size();
    R::PendingEvent pending(plan);const auto accepted=pending.TryCommit(&backend,&error);
    if(accepted!=R::AttemptStatus::Committed) {
      if(accepted==R::AttemptStatus::Deferred)
        Fatal("complete reaction packet cannot fit available PIC particle slots (required="+
          std::to_string(plan.products().size())+", available="+
          std::to_string(PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum())+"); "+error+
          ". Increase particle-buffer capacity or reduce initial model-particle count / weight resolution. The run terminates without discarding products, accepting a partial event, or redrawing its channel.");
      Fatal("conservative native event publication failed before continuation: "+error);
    }
    if(Active->configuration.ledger->events.size()!=ledgerBefore+1)Fatal("accepted reaction did not append exactly one event ledger");
    auto& committed=Active->configuration.ledger->events.back();
    committed.physicalEpochS=context.physicalEpochS;
    committed.exposureApproximationErrorBound=context.exposureApproximationErrorBound;
    committed.clockStorageErrorBound=context.clockStorageErrorBound;
    // All products already exist atomically with the accepted event ledger.
    // Opaque adapter state is copied without RNG. Core queues preserve their
    // event epoch and advance ONLY the actual unused physical interval.
    for(long int product:backend.committedProducts()) {
      PrepareProduct(product,remaining);
      PIC::Mover::ForceCutSurface::QueueDetachedParticle(product,remaining,node,time.IntegrationOrder);
    }
    return PIC::Mover::TrajectoryEventDelete; // core deletes retained parent once
  }catch(const std::bad_alloc&){ResourceFatal();}
  catch(const std::exception& error){Fatal(error.what());}
}
} } }
