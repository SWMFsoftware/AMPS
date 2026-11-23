#include "reduced_shock_particles.h"
#include "reduced_shock_geometry.h"
#include "util/sep_run_configuration.h"
#include "util/sep_source_lifecycle.h"
#include <sep_coronal_cme/particle_source.h>
#include <sep_coronal_cme/constants.h>
#include <algorithm>
#include <cstring>
#include <map>
#include <mutex>
#include <numeric>
#include <tuple>
#include <iomanip>

namespace SEP { namespace ReducedParticles { namespace Native {
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace T=Transport;
namespace {
// Measures are physical represented protons and TOTAL represented energies
// [J], not macroparticle counts or energies per numerical representative.
// Birth energy is invariant under transport; event energy is evaluated at the
// resolved return/exit/retirement epoch. Keep both to distinguish cooling from
// disappearance of an injected population.
struct Measure {double number=0,birthJ=0,eventJ=0;};
struct Entry {std::uint64_t id=0;int term=0;Measure measure;double backgroundJ=0,waveJ=0;};
std::mutex ledgerMutex;
std::vector<Entry> pending;
Measure proposed,excluded,born,lost[3];double backgroundWorkJ=0,waveWorkJ=0;
std::map<std::uint64_t,Entry> beforeGeometry;
T::Status Error(const std::string& why){return T::Status::Error(T::StatusCode::InvalidParticleState,why);}
T::Status Agree(T::Status local) {
  // Source preparation is collective even on ranks owning no candidate. The
  // lowest failing rank supplies one diagnostic, so successful ranks cannot
  // publish particles/counters while a peer has already rejected its state.
  int failure=local.ok()?std::numeric_limits<int>::max():PIC::ThisThread,first=0;
  MPI_Allreduce(&failure,&first,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(first==std::numeric_limits<int>::max())return T::Status::Ok();
  char why[2048]={};if(first==PIC::ThisThread)std::snprintf(why,sizeof(why),"reduced particle transaction rank %d: %s",first,local.message.c_str());
  MPI_Bcast(why,sizeof(why),MPI_CHAR,first,MPI_GLOBAL_COMMUNICATOR);return Error(why);
}
double Energy(double p,double mass) {
  // Rationalized relativistic kinetic energy [J], retaining low-energy
  // accuracy without subtracting two nearly equal rest/total energies.
  const double pc=p*SpeedOfLight,rest=mass*SpeedOfLight*SpeedOfLight;
  return pc*(pc/(std::hypot(pc,rest)+rest));
}
SourceRecord Record(const unsigned char* data) {
  SourceRecord r;std::memcpy(&r,data+Offset::ReducedSourceRecord,sizeof(r));
  if(r.schema!=1||!std::isfinite(r.birthKineticEnergyJ)||r.birthKineticEnergyJ<=0)
    throw std::runtime_error("reduced particle has no valid source-cohort birth energy");
  return r;
}
std::map<std::uint64_t,Entry> Inventory() {
  std::map<std::uint64_t,Entry> result;
  for(int line=0;line<FL::nFieldLine;++line)for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext())
    if(segment->Thread==PIC::ThisThread)for(long p=segment->FirstParticleIndex;p!=-1;p=PB::GetNext(p)) {
      const auto* data=PB::GetParticleDataPointer(p);const auto source=Record(data);
      Entry entry;std::memcpy(&entry.id,data+Offset::StableParticleId,sizeof(entry.id));
      const int species=PB::GetI(p);const double mass=PIC::MolecularData::GetMass(species);
      const double speed=std::hypot(PB::GetVParallel(p),PB::GetVNormal(p));
      const auto momentum=T::MomentumFromSpeed(speed,mass,SpeedOfLight);
      if(!momentum.status.ok())throw std::runtime_error(momentum.status.message);
      entry.measure.number=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]*PB::GetIndividualStatWeightCorrection(p);
      entry.measure.birthJ=entry.measure.number*source.birthKineticEnergyJ;
      entry.measure.eventJ=entry.measure.number*Energy(momentum.value,mass);
      if(!result.emplace(entry.id,entry).second)throw std::runtime_error("duplicate reduced particle identity in joined inventory");
    }
  return result;
}
void Add(Measure& a,const Measure& b){a.number+=b.number;a.birthJ+=b.birthJ;a.eventJ+=b.eventJ;}
}
bool Enabled(){return ReducedShock::Enabled()&&Run::Active().get().reducedParticles.enabled;}
T::Status ValidateLedger(const LedgerState& state) {
  for(std::size_t i=0;i<state.size();++i)
    if(!std::isfinite(state[i])||(i<18&&state[i]<0))
      return Error("checkpoint source ledger has an invalid represented measure or signed SI work");
  return T::Status::Ok();
}
T::Status CaptureLedger(LedgerState* state) {
  if(!state)return Error("missing reduced ledger checkpoint destination");
  std::lock_guard<std::mutex> lock(ledgerMutex);
  if(!pending.empty()||!beforeGeometry.empty())return Error("source ledger checkpoint requires joined transport/geometry accounting");
  LedgerState staged{};std::size_t at=0;
  for(const auto* measure:{&proposed,&excluded,&born,&lost[0],&lost[1],&lost[2]}) {
    staged[at++]=measure->number;staged[at++]=measure->birthJ;staged[at++]=measure->eventJ;
  }
  staged[18]=backgroundWorkJ;staged[19]=waveWorkJ;
  const auto valid=ValidateLedger(staged);if(valid.ok())*state=staged;return valid;
}
void RestoreLedger(const LedgerState& state) noexcept {
  // Validation and worker joins belong to the enclosing collective restore.
  // Scalar assignment cannot allocate or silently zero historical losses.
  std::size_t at=0;
  for(auto* measure:{&proposed,&excluded,&born,&lost[0],&lost[1],&lost[2]}) {
    measure->number=state[at++];measure->birthJ=state[at++];measure->eventJ=state[at++];
  }
  backgroundWorkJ=state[18];waveWorkJ=state[19];
}
long int InjectForNativeDriver() {
  if(!Enabled())return 0;
  T::Birth::DriverWindow window;
  auto status=T::PICAdapter::GetNativeDriverWindow(&window);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  const long before=PB::GetAllPartNum();
  status=Inject(window);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  return PB::GetAllPartNum()-before;
}
T::Status Inject(const T::Birth::DriverWindow& window) {
  if(!Enabled())return T::Status::Ok();
  const auto& config=Run::Active().get().reducedParticles;
  auto status=Validate(config);if(!status.ok())return Agree(status);
  const auto canonical=Canonical(config)+" mover="+std::to_string(static_cast<int>(Mover::CurrentProductionMover()));
  char reference[1024]={};
  if(PIC::ThisThread==0&&canonical.size()<sizeof(reference))std::memcpy(reference,canonical.data(),canonical.size());
  MPI_Bcast(reference,sizeof(reference),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  status=Agree(canonical.size()<sizeof(reference)&&canonical==reference?T::Status::Ok():
      Error("MPI ranks selected different reduced source/return configurations"));
  if(!status.ok())return status;
  const auto material=ReducedShock::Native::CaptureMaterialView();
  if(!material.ok())return Agree(Error(material.status.message));
  const auto& mesh=*material.value;
  struct Candidate {std::uint64_t id=0;double coordinate=0,speed=0,mu=0,momentum=0,magneticT=0,number=0,energyJ=0;T::Birth::Record birth;};
  std::vector<Candidate> candidates;
  Measure stagedProposed,stagedExcluded;
  SourceLifecycle::Scheduler stagedScheduler;
  SourceLifecycle::Event event;
  try {
    // This first consumer binds only explicitly compiled protons to the core
    // proton density. An ion mixture needs its own species abundance contract;
    // ne, rho/mp and an arbitrary mass slot are not substitutes for n_p.
    if(PIC::nTotalSpecies!=1||PIC::MolecularData::GetElectricCharge(0)!=ElectronCharge||
        (PIC::MolecularData::GetMass(0)!=CoronalCME::Constants::kProtonMassKg&&PIC::MolecularData::GetMass(0)!=ProtonMass))
      throw std::runtime_error("first reduced net source requires exactly one compiled proton species with a declared proton mass; legacy H_plus=1 amu is not the proton mass");
    if(mesh.interval.DepartureS()!=window.startEpochS||mesh.interval.ArrivalS()<window.endEpochS)
      throw std::runtime_error("source window differs from the retained material realization");
    status=SourceLifecycle::PrepareNativeRestore(SourceLifecycle::CaptureNativeState(),&stagedScheduler);
    if(!status.ok())throw std::runtime_error(status.message);
    SourceLifecycle::Channel channel;channel.kind=SourceLifecycle::Kind::Shock;channel.line=0;channel.species=0;
    status=stagedScheduler.Begin(channel,&event);if(!status.ok())throw std::runtime_error(status.message);
    const double mass=PIC::MolecularData::GetMass(0),rest=mass*SpeedOfLight*SpeedOfLight;
    const double pmin=std::sqrt(config.minimumEnergyJ*(config.minimumEnergyJ+2*rest))/SpeedOfLight;
    const double pmax=std::sqrt(config.maximumEnergyJ*(config.maximumEnergyJ+2*rest))/SpeedOfLight;
    const double flux=FieldLine::FluxTubeGeometry::GetMagneticFluxWb(0);
    for(std::uint64_t j=0;j<config.samplesPerStep;++j) {
      // Stratified temporal Monte Carlo quadrature: one independent draw per
      // time stratum, with physical weight R(t_j)*dt/N. The sum estimates the
      // actual time-dependent NET source integral without assuming a constant
      // rate, redrawing disconnected births, or borrowing PIC's current clock.
      // Purpose-separated keys keep time, momentum and pitch independent and
      // preserve identities when ownership or material segment indices change.
      auto time=SourceLifecycle::Random(event,j,UINT64_C(0x525354494d45));
      const double phase=(j+time.UniformOpen01())/config.samplesPerStep;
      T::Birth::Record birth;status=T::Birth::UniformBirth(window,phase,&birth);
      if(!status.ok())throw std::runtime_error(status.message);
      const double epoch=window.endEpochS-birth.firstMoveDurationS;
      std::vector<ReducedShock::Position> points;
      const auto controls=ReducedShock::Native::Controls(window.durationS,mesh.arc[1]);
      for(const auto& point:mesh.vertices) {
        const auto moved=ReducedShock::MoveVertex(mesh.interval,point,window.startEpochS,epoch,controls);
        if(!moved.ok())throw std::runtime_error(moved.status.message);
        points.push_back(moved.value);
      }
      const auto roots=ReducedShock::ConnectLine(mesh.interval,points,epoch,controls.positionToleranceM);
      if(!roots.ok())throw std::runtime_error(roots.status.message);
      if(std::count_if(roots.value.begin(),roots.value.end(),[](const ReducedShock::Connection& c){return c.front.acceptedShock;})>1)
        throw std::runtime_error("multiple simultaneous shock barriers require a joint first-passage operator; source allocation is rejected before any live mutation");
      const auto releases=ReferenceIntersections(mesh.interval,points,epoch,config.releaseDistanceM,controls.positionToleranceM);
      if(!releases.ok())throw std::runtime_error(releases.status.message);
      for(std::size_t root=0;root<roots.value.size();++root) {
        const auto& connection=roots.value[root];if(!connection.front.acceptedShock)continue;
        if(connection.tangent)throw std::runtime_error("tangent shock source has no qualified finite tube patch");
        const auto rate=NetProtonRate(connection.front,flux,config.releaseFraction);
        if(!rate.ok())throw std::runtime_error(rate.status.message);
        if(rate.value==0)continue;
        // A line may meet both sides of a sphere. Select the first supported
        // upstream offset root in the outward-normal direction of this exact
        // connected shock, rather than borrowing another branch's release.
        const auto& a=points[connection.segment];const auto& b=points[connection.segment+1];
        double projection=0;for(int d=0;d<3;++d)projection+=(b[d]-a[d])*connection.front.outwardNormal[d];
        if(projection==0)throw std::runtime_error("source material chord tangent to shock");
        const double shockS=connection.segment+connection.fraction;
        const ReleasePoint* selected=nullptr;double closest=std::numeric_limits<double>::infinity();
        for(const auto& release:releases.value) {
          const double delta=release.segment+release.fraction-shockS;
          if(delta*projection>0&&std::abs(delta)<closest){selected=&release;closest=std::abs(delta);}
        }
        if(!selected)throw std::runtime_error("configured L_ref has no upstream intersection on the same supported material line");
        Core::Result<CoronalCME::MomentumSpectrum> spectrum;
        if(config.verificationMonoenergeticEnergyJ==0) {
          spectrum=CoronalCME::MomentumSpectrum::LocalCompressionDsa(pmin,pmax,connection.front.compressionRatio);
          if(!spectrum.ok())throw std::runtime_error(spectrum.status.message);
        }
        if(root>std::numeric_limits<std::uint64_t>::max()/config.samplesPerStep)
          throw std::runtime_error("source root/sample identity overflow");
        const std::uint64_t macro=root*config.samplesPerStep+j;
        auto spectral=SourceLifecycle::Random(event,macro,UINT64_C(0x525350));
        auto angular=SourceLifecycle::Random(event,macro,UINT64_C(0x52534d55));
        const double mono=config.verificationMonoenergeticEnergyJ;
        const double p=mono>0?std::sqrt(mono*(mono+2*rest))/SpeedOfLight:
            spectrum.value.InverseCdf(spectral.UniformOpen01());
        Measure proposal;proposal.number=rate.value*window.durationS/config.samplesPerStep;
        proposal.birthJ=proposal.eventJ=proposal.number*Energy(p,mass);
        if(!std::isfinite(proposal.number)||proposal.number<=0||!std::isfinite(proposal.birthJ)||proposal.birthJ<=0)
          throw std::runtime_error("candidate source number/energy overflows finite SI storage");
        Add(stagedProposed,proposal);
        const auto speed=T::SpeedFromMomentum(p,mass,SpeedOfLight);if(!speed.status.ok())throw std::runtime_error(speed.status.message);
        const auto ambient=mesh.interval.PrimitiveAt(selected->xM,epoch);if(!ambient.ok())throw std::runtime_error(ambient.status.message);
        const auto& chord0=points[selected->segment];const auto& chord1=points[selected->segment+1];
        double bn=0,bmag=0,orient=0,advection=-selected->normalFoot.normalSpeedMPerS;
        for(int d=0;d<3;++d){bn+=ambient.value.magneticFieldT[d]*selected->normal[d];bmag=std::hypot(bmag,ambient.value.magneticFieldT[d]);
          orient+=ambient.value.magneticFieldT[d]*(chord1[d]-chord0[d]);advection+=ambient.value.velocityMPerS[d]*selected->normal[d];}
        if(!(bmag>0)||orient==0)throw std::runtime_error("reference source lacks magnetic pitch orientation");
        // This is a prescribed NET escaping spectrum, not a gross emission
        // spectrum to be filtered by an unknown accelerator-to-reference law.
        // Every momentum in its finite support must admit positive outward
        // normal velocity here. That condition is monotone in particle speed,
        // so testing the minimum momentum (or the explicit verification
        // monoenergy) suffices for the entire configured band. Dropping slow
        // sampled proposals would reduce the prescribed net number and change
        // its DSA marginal; redrawing or renormalizing would change it again.
        // Reject the incompatible contract before any allocation/counter
        // commit instead of inventing that missing physical escape law.
        const auto minimumSpeed=T::SpeedFromMomentum(mono>0?p:pmin,mass,SpeedOfLight);
        if(!minimumSpeed.status.ok())throw std::runtime_error(minimumSpeed.status.message);
        const auto minimumEscape=CoronalCME::BuildFocusedEscapeDistribution(
            minimumSpeed.value,orient>0?1:-1,bn/bmag,advection,64,0);
        if(!minimumEscape.ok())throw std::runtime_error(minimumEscape.status.message);
        if(minimumEscape.value.status==CoronalCME::FocusedEscapeStatus::NoFocusedEscape)
          throw std::runtime_error("prescribed net DSA source includes momenta with no outward reference-surface escape; specify a compatible energy band, not a filtered or renormalized net spectrum");
        const auto angularLaw=CoronalCME::BuildFocusedEscapeDistribution(speed.value,orient>0?1:-1,bn/bmag,advection,64,0);
        if(!angularLaw.ok())throw std::runtime_error(angularLaw.status.message);
        if(angularLaw.value.status==CoronalCME::FocusedEscapeStatus::NoFocusedEscape) {
          throw std::runtime_error("realized net source momentum has no outward escape despite minimum-energy preflight");
        }
        const double pitchDraw=angular.UniformOpen01();
        // Focused equations retain the positive-normal-flux angular source.
        // Parker evolves its ANGULAR MARGINAL: its stored diagnostic direction
        // must represent the isotropic projection rather than preserve an
        // injection anisotropy forever in a solver with no pitch-angle state.
        // This projection leaves p, source number, energy and release position
        // identical; it does not introduce a fourth transport prescription.
        const auto mu=Mover::CurrentCapabilities().usesPitchAngleState?
            CoronalCME::SampleFocusedMu(angularLaw.value,pitchDraw):
            Core::Result<double>::Success(2*pitchDraw-1);
        if(!mu.ok())throw std::runtime_error(mu.status.message);
        Candidate candidate;candidate.id=SourceLifecycle::ParticleIdentity(event,macro);
        candidate.coordinate=selected->segment+selected->fraction;candidate.speed=speed.value;candidate.mu=mu.value;
        candidate.momentum=p;candidate.magneticT=bmag;
        candidate.number=rate.value*window.durationS/config.samplesPerStep;candidate.energyJ=Energy(p,mass);candidate.birth=birth;
        if(!std::isfinite(candidate.number)||candidate.number<=0||!std::isfinite(candidate.energyJ)||candidate.energyJ<=0)
          throw std::runtime_error("source weight/energy is not finite positive SI data");
        candidates.push_back(candidate);
      }
    }
  } catch(const std::exception& e){status=Error(e.what());}
  status=Agree(status);if(!status.ok())return status;
  const std::uint64_t eventKey[2]={event.campaign,event.index};std::uint64_t low[2],high[2];
  MPI_Allreduce(eventKey,low,2,MPI_UINT64_T,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(eventKey,high,2,MPI_UINT64_T,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(low[0]!=high[0]||low[1]!=high[1])return Error("source history differs across MPI owners before allocation");
  // Allocate privately, then collectively validate before attaching any live
  // list or committing scheduler/ledger state. A failed preparation does not
  // consume event zero or leave a half-injected source on successful ranks.
  std::vector<std::pair<long,Candidate>> allocations;
  try {
    const double base=PIC::ParticleWeightTimeStep::GlobalParticleWeight[0];
    if(!std::isfinite(base)||base<=0)throw std::runtime_error("source base physical weight is invalid");
    for(const auto& candidate:candidates) {
      const int segment=static_cast<int>(candidate.coordinate);
      if(!FL::FieldLinesAll[0].GetSegment(segment)||candidate.coordinate>=mesh.vertices.size()-1)
        throw std::runtime_error("reference source is at or outside the open native endpoint");
      if(FL::FieldLinesAll[0].GetSegment(segment)->Thread!=PIC::ThisThread)continue;
      const double correction=candidate.number/base;
      if(!std::isfinite(correction)||correction<=0)throw std::runtime_error("source weight correction is unrepresentable");
      const long pointer=PB::GetNewParticle();allocations.emplace_back(pointer,candidate);
      auto* data=PB::GetParticleDataPointer(pointer);std::memset(data,0,PB::ParticleDataLength);
      PB::SetParticleAllocated(data);PB::SetNext(-1,data);PB::SetPrev(-1,data);PB::SetI(0,data);
      PB::SetFieldLineId(0,data);PB::SetFieldLineCoord(candidate.coordinate,data);
      PB::SetIndividualStatWeightCorrection(correction,data);
      PB::SetVParallel(candidate.speed*candidate.mu,data);
      PB::SetVNormal(candidate.speed*std::sqrt((1-candidate.mu)*(1+candidate.mu)),data);
      // Keep optional diagnostic momentum/magnetic-moment slots consistent
      // from birth; the scalar plasma-frame speed/pitch remain mover authority.
      const double perpendicular=candidate.momentum*std::sqrt((1-candidate.mu)*(1+candidate.mu));
      if(PB::_PIC_PARTICLE_DATA__MOMENTUM_PARALLEL_>=0&&PB::_PIC_PARTICLE_DATA__MOMENTUM_NORMAL_>=0) {
        PB::SetMomentumParallel(candidate.momentum*candidate.mu,data);
        PB::SetMomentumNormal(perpendicular,data);
      }
      // Field-line-only configurations often omit the Cartesian guiding-centre
      // moment slot (offset -1). Never write a diagnostic into absent storage.
      if(_PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_>=0)
        PB::SetMagneticMoment(perpendicular*perpendicular/(2*PIC::MolecularData::GetMass(0)*candidate.magneticT),data);
      double x[3];FL::FieldLinesAll[0].GetCartesian(x,candidate.coordinate);PB::SetX(x,data);
      status=T::PICAdapter::InitializeParticleTransportStateWithIdentity(pointer,candidate.id,candidate.birth);
      if(!status.ok())throw std::runtime_error(status.message);
      SourceRecord record;record.schema=1;record.birthKineticEnergyJ=candidate.energyJ;
      static_assert(sizeof(record)==16,"source history must contain two eight-byte scalar slots");
      std::memcpy(data+Offset::ReducedSourceRecord,&record,sizeof(record));
    }
  }catch(const std::exception& e){status=Error(e.what());}
  status=Agree(status);
  if(!status.ok()){for(const auto& a:allocations)PB::DeleteParticle(a.first);return status;}
  for(const auto& a:allocations) {
    auto* segment=FL::FieldLinesAll[0].GetSegment(static_cast<int>(a.second.coordinate));
    PB::SetNext(segment->FirstParticleIndex,a.first);PB::SetPrev(-1,a.first);
    if(segment->FirstParticleIndex!=-1)PB::SetPrev(a.first,segment->FirstParticleIndex);
    segment->FirstParticleIndex=a.first;
    Measure measure;measure.number=a.second.number;measure.birthJ=measure.eventJ=measure.number*a.second.energyJ;Add(born,measure);
  }
  if(PIC::ThisThread==0) {
    Add(proposed,stagedProposed);Add(excluded,stagedExcluded);
    std::cout<<std::setprecision(17)<<"REDUCED_SOURCE event="<<event.index
        <<" candidate="<<stagedProposed.number<<" no_focused_escape="<<stagedExcluded.number
        <<" committed="<<stagedProposed.number-stagedExcluded.number
        <<" reference_distance_m="<<config.releaseDistanceM<<" spectrum="
        <<(config.verificationMonoenergeticEnergyJ>0?"manufactured-monoenergetic-verification":"core-local-compression-dsa")<<std::endl;
  }
  SourceLifecycle::CommitNativeRestore(&stagedScheduler);return T::Status::Ok();
}
void TransportWork(const T::PICAdapter::ParticleContext& context,double before,double after,double waveJ) {
  if(!Enabled())return;
  Entry row;row.id=context.stableParticleId;row.term=-1;
  row.waveJ=context.statisticalWeight*waveJ;
  row.backgroundJ=context.statisticalWeight*(Energy(after,context.state.massKg)-Energy(before,context.state.massKg))-row.waveJ;
  if(!std::isfinite(row.waveJ)||!std::isfinite(row.backgroundJ))
    exit(__LINE__,__FILE__,"nonfinite reduced particle work cannot enter the physical ledger");
  std::lock_guard<std::mutex> lock(ledgerMutex);pending.push_back(row);
}
void RecordTerminal(const T::PICAdapter::ParticleContext& context,Terminal terminal,double momentum) {
  if(!Enabled())return;
  SourceRecord source;
  try{source=Record(context.data);}catch(const std::exception& e){exit(__LINE__,__FILE__,e.what());}
  Entry row;row.id=context.stableParticleId;row.term=static_cast<int>(terminal);
  if(row.term<0||row.term>2)exit(__LINE__,__FILE__,"unknown reduced terminal-loss channel");
  row.measure.number=context.statisticalWeight;row.measure.birthJ=row.measure.number*source.birthKineticEnergyJ;
  row.measure.eventJ=row.measure.number*Energy(momentum,context.state.massKg);
  std::lock_guard<std::mutex> lock(ledgerMutex);pending.push_back(row);
}
T::Status BeforeGeometry(){
  if(!Enabled())return T::Status::Ok();
  std::map<std::uint64_t,Entry> staged;T::Status status=T::Status::Ok();
  try{staged=Inventory();}catch(const std::exception& e){status=Error(e.what());}
  status=Agree(status);if(status.ok())beforeGeometry.swap(staged);return status;
}
T::Status AfterGeometry(double epoch) {
  if(!Enabled())return T::Status::Ok();
  std::map<std::uint64_t,Entry> inventory;T::Status status=T::Status::Ok();
  try{inventory=Inventory();}catch(const std::exception& e){status=Error(e.what());}
  status=Agree(status);if(!status.ok())return status;
  try {
    Measure retained;
    for(const auto& pair:beforeGeometry)if(!inventory.count(pair.first))Add(lost[2],pair.second.measure);
    beforeGeometry.clear();
    for(const auto& pair:inventory)Add(retained,pair.second.measure);
    std::stable_sort(pending.begin(),pending.end(),[](const Entry& a,const Entry& b){return std::tie(a.id,a.term)<std::tie(b.id,b.term);});
    for(const auto& row:pending)if(row.term<0){backgroundWorkJ+=row.backgroundJ;waveWorkJ+=row.waveJ;}else Add(lost[row.term],row.measure);
    pending.clear();
    double local[14]={born.number,born.birthJ,retained.number,retained.birthJ,retained.eventJ,
      lost[0].number,lost[0].birthJ,lost[0].eventJ,lost[1].number,lost[1].birthJ,lost[1].eventJ,
      lost[2].number,lost[2].birthJ,lost[2].eventJ};double total[14];
    MPI_Allreduce(local,total,14,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    const double sourceMeasures[4]={proposed.number,proposed.birthJ,excluded.number,excluded.birthJ};double sourceTotal[4];
    MPI_Allreduce(sourceMeasures,sourceTotal,4,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    const double work[2]={backgroundWorkJ,waveWorkJ};double globalWork[2];
    MPI_Allreduce(work,globalWork,2,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    for(double value:total)if(!std::isfinite(value))return Error("reduced number/energy inventory overflowed finite SI storage");
    for(double value:sourceTotal)if(!std::isfinite(value))return Error("reduced source/exclusion ledger is nonfinite");
    for(double value:globalWork)if(!std::isfinite(value))return Error("reduced prescribed-background work ledger is nonfinite");
    // This is an open transport system: the already accelerated source and
    // prescribed scatterers supply explicit energy. Number and birth-energy
    // budgets close without work terms; the current kinetic-energy budget
    // includes signed deterministic and wave-frame work. No loss is repaired
    // by scaling the survivors or redefining the requested source strength.
    const double numberResidual=total[0]-total[2]-total[5]-total[8]-total[11];
    const double birthResidual=total[1]-total[3]-total[6]-total[9]-total[12];
    const double energyResidual=total[1]+globalWork[0]+globalWork[1]-total[4]-total[7]-total[10]-total[13];
    const double energyScale=std::abs(total[1])+std::abs(globalWork[0])+std::abs(globalWork[1]);
    if(std::abs(sourceTotal[0]-sourceTotal[2]-total[0])>1e-10*sourceTotal[0]||
        std::abs(sourceTotal[1]-sourceTotal[3]-total[1])>1e-10*sourceTotal[1]||
        std::abs(numberResidual)>1e-10*total[0]||std::abs(birthResidual)>1e-10*total[1]||std::abs(energyResidual)>1e-10*energyScale)
      return Error("reduced source/return/domain/geometry represented number or signed energy ledger does not close");
    if(PIC::ThisThread==0)std::cout<<std::setprecision(17)<<"REDUCED_SOURCE_LEDGER epoch_s="<<epoch
      <<" born="<<total[0]<<" retained="<<total[2]<<" returned="<<total[5]<<" domain_exit="<<total[8]<<" geometry_retired="<<total[11]
      // Birth energy and current kinetic energy differ after transport. Both
      // independently accumulated measures are needed by an external audit.
      <<" born_kinetic_j="<<total[1]<<" retained_birth_j="<<total[3]<<" retained_kinetic_j="<<total[4]
      <<" background_work_j="<<globalWork[0]<<" wave_work_j="<<globalWork[1]
      <<" returned_birth_j="<<total[6]<<" returned_event_j="<<total[7]
      <<" domain_birth_j="<<total[9]<<" domain_event_j="<<total[10]
      <<" retired_birth_j="<<total[12]<<" retired_event_j="<<total[13]
      <<" number_residual="<<numberResidual<<" birth_energy_residual_j="<<birthResidual<<" energy_residual_j="<<energyResidual<<std::endl;
  }catch(const std::exception& e){return Error(e.what());}
  return T::Status::Ok();
}
} } }
