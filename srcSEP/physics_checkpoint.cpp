#include "physics_checkpoint.h"
#include "physics_checkpoint_state.h"
#include "sep.h"
#include "transport_common.h"
#include "turbulence_production_adapter.h"
#include "polarized_wave_native_adapter.h"
#include "flux_tube_checkpoint_state.h"
#include "util/sep_background_runtime.h"
#include "util/sep_initialization.h"
#include "util/sep_run_configuration.h"
#include "util/sep_source_lifecycle.h"
#include "turbulence_wave_number_resolved/turbulence_wave_number_resolved.h"
#include "pic_persistent_format.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

// Application-owned corrected physics checkpoint, not an AMR mesh dump. The
// first supported realization is fixed open analytic field lines on identical
// CPU-MPI rank ownership, with explicitly typed coefficients/global numerics.
// SI physical buffers and complete stochastic/source histories are encoded;
// native objects, pointers, atomics and allocation/list handles are not.
// Decode/layout/physics/global-identity checks and replacement allocation all
// precede the collective commit barrier. The old population remains live until
// every rank can commit through safe PIC allocation/link APIs. Output/sample
// continuation is explicitly excluded. Active births and live four-spectrum
// authority are rejected rather than approximated by this retained codec.
namespace SEP { namespace PhysicsCheckpoint { namespace {
namespace T=Transport;
namespace A=Transport::PICAdapter;
namespace FL=PIC::FieldLine;
namespace PB=PIC::ParticleBuffer;
namespace W=SEP::AlfvenTurbulence_Kolmogorov;
namespace ShockOwner=SEP::ParticleSource::ShockWave::Tenishev2005;
namespace TurbOwner=SEP::Turbulence::PICAdapter;
namespace FluxOwner=SEP::FieldLine::FluxTubeGeometry;
const std::uint64_t kSchema=4;
const std::uint64_t kMagic=UINT64_C(0x5345505048595331);

T::Status Error(const std::string& message) {
  return T::Status::Error(T::StatusCode::InvalidParticleState,message);
}
struct Writer {
  std::string bytes;
  void U(std::uint64_t value) {
    for(unsigned i=0;i<8;++i)bytes.push_back(static_cast<char>(value>>(8*i)));
  }
  void D(double value) {std::uint64_t bits;std::memcpy(&bits,&value,8);U(bits);}
  void S(const std::string& value) {U(value.size());bytes.append(value);}
  void Raw(const void* data,std::size_t size) {
    if(size==0){S(std::string());return;}
    if(!data)throw std::runtime_error("missing checkpoint native buffer");
    S(std::string(static_cast<const char*>(data),size));
  }
  void Doubles(const std::vector<double>& values) {
    U(values.size());for(double value:values)D(value);
  }
};
struct Reader {
  const std::string& bytes;std::size_t at=0;
  explicit Reader(const std::string& input):bytes(input){}
  void Need(std::uint64_t size) {
    if(size>bytes.size()-at)throw std::runtime_error("truncated physics checkpoint section");
  }
  std::uint64_t U() {
    Need(8);std::uint64_t value=0;
    for(unsigned i=0;i<8;++i)value|=static_cast<std::uint64_t>(
        static_cast<unsigned char>(bytes[at++]))<<(8*i);
    return value;
  }
  double D() {const auto bits=U();double value;std::memcpy(&value,&bits,8);return value;}
  std::string S() {const auto size=U();Need(size);const auto begin=at;at+=size;return bytes.substr(begin,size);}
  std::size_t Count(std::uint64_t minimumBytes) {
    const auto count=U();
    if(minimumBytes&&count>(bytes.size()-at)/minimumBytes)
      throw std::runtime_error("physics checkpoint table exceeds payload");
    return static_cast<std::size_t>(count);
  }
  std::vector<double> Doubles() {
    const auto count=Count(8);std::vector<double> values(count);
    for(double& value:values)value=D();
    return values;
  }
  void End() {if(at!=bytes.size())throw std::runtime_error("trailing physics checkpoint data");}
};
void PutSource(Writer& out,const SourceLifecycle::State& state) {
  out.U(state.schema);out.U(state.campaign);out.S(state.fixedLineLayout);
  out.U(state.counters.size());
  for(const auto& row:state.counters) {
    out.U(static_cast<std::uint64_t>(row.channel.kind));out.U(row.channel.line);
    out.U(row.channel.species);out.U(row.next);
  }
}
SourceLifecycle::State GetSource(Reader& in) {
  SourceLifecycle::State state;state.schema=in.U();state.campaign=in.U();state.fixedLineLayout=in.S();
  const auto count=in.Count(32);state.counters.resize(count);
  for(auto& row:state.counters) {
    row.channel.kind=static_cast<SourceLifecycle::Kind>(in.U());
    row.channel.line=in.U();row.channel.species=in.U();row.next=in.U();
  }
  return state;
}
void PutBoundary(Writer& out,const BoundaryCache& cache) {
  out.U(static_cast<std::uint64_t>(cache.lineCount+1));out.U(cache.components);
  out.Doubles(cache.values);out.U(cache.spectra.size());
  for(const auto& spectrum:cache.spectra)out.Doubles(spectrum);
  out.Raw(cache.isSet.data(),cache.isSet.size());
}
BoundaryCache GetBoundary(Reader& in) {
  BoundaryCache cache;const auto lines=in.U();
  if(lines>static_cast<std::uint64_t>(INT_MAX))throw std::runtime_error("invalid checkpoint boundary line count");
  cache.lineCount=static_cast<int>(lines)-1;cache.components=in.U();cache.values=in.Doubles();
  const auto count=in.Count(8);cache.spectra.resize(count);
  for(auto& spectrum:cache.spectra)spectrum=in.Doubles();
  const auto flags=in.S();cache.isSet.assign(flags.begin(),flags.end());return cache;
}
void PutTurbulence(Writer& out,const TurbOwner::CheckpointCache& cache) {
  const auto& a=cache.ledger;
  out.D(a.shockSourceJ);out.D(a.particleExchangeJ);out.D(a.boundaryExchangeJ);
  out.D(a.physicalDissipationJ);out.D(a.limiterCorrectionJ);out.D(a.closureResidualJ);
  out.U(a.coreSubcycles);out.U(a.limiterActivations);out.U(a.rejectedUpdates);
  out.U(a.advectionFaceUpdates);out.U(a.reflectionCellUpdates);out.U(a.cascadeCellUpdates);out.U(a.spectralBinUpdates);
  const auto& d=cache.shockDiagnostics;
  out.U(d.accepted);out.U(d.noIntersection);out.U(d.invalidGeometry);out.U(d.invalidPhysics);
  out.U(cache.shockContributions.size());
  for(const auto& source:cache.shockContributions) {
    out.U(source.fieldLine);out.U(source.segment);out.D(source.plusJ);out.D(source.minusJ);out.S(source.provenance);
  }
}
TurbOwner::CheckpointCache GetTurbulence(Reader& in) {
  TurbOwner::CheckpointCache cache;auto& a=cache.ledger;
  a.shockSourceJ=in.D();a.particleExchangeJ=in.D();a.boundaryExchangeJ=in.D();
  a.physicalDissipationJ=in.D();a.limiterCorrectionJ=in.D();a.closureResidualJ=in.D();
  a.coreSubcycles=in.U();a.limiterActivations=in.U();a.rejectedUpdates=in.U();
  a.advectionFaceUpdates=in.U();a.reflectionCellUpdates=in.U();a.cascadeCellUpdates=in.U();a.spectralBinUpdates=in.U();
  auto& d=cache.shockDiagnostics;
  d.accepted=in.U();d.noIntersection=in.U();d.invalidGeometry=in.U();d.invalidPhysics=in.U();
  const auto count=in.Count(40);cache.shockContributions.resize(count);
  for(auto& source:cache.shockContributions) {
    const auto line=in.U(),segment=in.U();
    if(line>INT_MAX||segment>INT_MAX)throw std::runtime_error("invalid checkpoint source index");
    source.fieldLine=line;source.segment=segment;source.plusJ=in.D();source.minusJ=in.D();source.provenance=in.S();
  }
  return cache;
}

struct SegmentImage {
  std::uint64_t owner=0;
  std::string data;
  std::vector<double> weights;
  std::vector<std::string> particles;
  long int stagedHead=-1; // private allocation list; never encoded
};
struct LineImage {
  // Vertex positions live outside PIC's associated-data buffer. Persist them
  // fieldwise so decoded geometry, rather than the mutable receiving process,
  // is the oracle for every derived particle Cartesian position. This first
  // fixed-geometry capability still requires exact live-position agreement.
  std::vector<std::array<double,3>> positions;
  std::vector<std::string> vertices;
  std::vector<SegmentImage> segments;
};
struct Image {
  std::uint64_t completed=0,rank=0,ranks=0,mover=0,purpose=0;
  double clock=0.0;
  std::string runConfiguration,activeConfiguration,initialization,coefficientFingerprint,backgroundFingerprint;
  SourceLifecycle::State source;
  std::shared_ptr<const Background::BackgroundSnapshot> snapshot;
  std::vector<double> steps,weights;
  bool stepsInitialized=false;
  unsigned long seed=0;
  std::mt19937 generator;
  std::uniform_real_distribution<double> distribution;
  bool hasNormal=false;
  double normal=0.0;
  BoundaryCache integrated,spectral;
  AnalyticShockCache shock;
  HostStepCache host;
  TurbOwner::CheckpointCache turbulence;
  FieldLine::FluxTubeState::State flux;
  std::vector<LineImage> lines;
  SourceLifecycle::Scheduler stagedSource; // allocated before physical commit
  FluxOwner::CheckpointRestorePlan stagedFlux;
};

std::string ActualConfiguration() {
  // The native adapters consume these installed typed configurations. Record
  // their actual values as well as the frozen startup manifest, so a controlled
  // fixture or host cannot change a coefficient/operator behind the manifest.
  auto actual=Run::Active().get();actual.mover.value=SEP::Mover::CurrentProductionMover();
  actual.coefficients=T::Coefficient::ActiveConfiguration();actual.turbulence=Turbulence::ActiveConfiguration();
  Run::FrozenConfiguration frozen;auto status=Run::FrozenConfiguration::Create(actual,&frozen);
  if(!status.ok())throw std::runtime_error(status.message);
  std::string typed;status=Run::Serialize(frozen,&typed);if(!status.ok())throw std::runtime_error(status.message);
  Writer out;out.S(typed);out.U(A::CampaignRandomSeed);out.U(SEP::AccountAdiabaticCoolingFlag);
  out.D(PIC::CPLR::SWMF::MeanPlasmaAtomicMass);
  out.U(W::ActiveFlag);out.U(W::ParticleCouplingMode);out.U(W::WaveNumberResolved::IsActive());
  namespace I=SEP::FieldLine::InjectionParameters;
  out.U(I::nParticlesPerIteration);out.U(I::InjectLocation);out.U(I::InjectionMomentumModel);out.U(I::UseAnalyticShockModel);
  out.D(I::PowerIndex);out.D(I::emin);out.D(I::emax);out.D(I::InjectionEfficiency);
  out.D(I::ConstEnergyInjectionValue);out.D(I::ConstSpeedInjectionValue);out.D(I::ConstMuInjectionValue);
  for(int species=0;species<PIC::nTotalSpecies;++species) {
    out.D(PIC::MolecularData::GetMass(species));out.D(PIC::MolecularData::GetElectricCharge(species));
  }
  return out.bytes;
}

std::uint64_t Layout() {
  Writer layout;layout.U(kMagic);layout.U(kSchema);layout.U(A::ParticleTransportSchema);
  layout.U(sizeof(T::OperatorRandomState));layout.U(PB::ParticleDataLength);
  layout.U(SEP::Offset::TransportSchema);layout.U(SEP::Offset::StableParticleId);
  layout.U(SEP::Offset::OperatorRandomStates);layout.U(SEP::Offset::MfpOpticalDepth);
  layout.U(SEP::Offset::MfpEventIndex);layout.U(SEP::Offset::SplitEventCounter);
  layout.U(SEP::Offset::BirthRecord);layout.U(sizeof(T::Birth::Record));
  layout.U(PB::_PIC_PARTICLE_DATA__MOMENTUM_NORMAL_);
  layout.U(PB::_PIC_PARTICLE_DATA__MOMENTUM_PARALLEL_);
  layout.U(PIC::nTotalSpecies);
  const auto add=[&layout](long long value){layout.U(static_cast<std::uint64_t>(value));};
  add(DIM);add(sizeof(bool));add(sizeof(unsigned long));
  add(sizeof(void*));add(sizeof(int));add(sizeof(long));add(sizeof(long long));
  add(sizeof(std::size_t));add(sizeof(float));add(sizeof(double));
  add(PIC::PersistentFormat::HostLittleEndian()?1:2);
  add(_COMPILATION_MODE_);add(_SIMULATION_TIME_STEP_MODE_);add(_SIMULATION_PARTICLE_WEIGHT_MODE_);
  add(_PIC_GLOBAL_TIME_COUNTER_MODE_);
  add(_STATE_VECTOR_MODE_);add(_INDIVIDUAL_PARTICLE_WEIGHT_MODE_);add(_PIC_FIELD_LINE_MODE_);
  add(_PIC_PARTICLE_LIST_ATTACHING_);add(_PIC_PARTICLE_TRACKER_MODE_);add(_RND_MODE_);
  add(_PIC_PARTICLE_DATA__NEXT_OFFSET_);add(_PIC_PARTICLE_DATA__PREV_OFFSET_);
  add(_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_);add(_PIC_PARTICLE_DATA__VELOCITY_OFFSET_);
  add(_PIC_PARTICLE_DATA__POSITION_OFFSET_);add(_PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_);
  add(_PIC_PARTICLE_DATA__BASIC_DATA_LENGTH_);add(_PIC_PARTICLE_DATA__FIELD_LINE_ID_OFFSET_);
  add(_PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_);add(_PIC_PARTICLE_DATA__V_PARALLEL_OFFSET_);
  add(_PIC_PARTICLE_DATA__V_NORMAL_OFFSET_);add(_PIC_PARTICLE_DATA__FULL_DATA_LENGTH_);
  add(_PIC_PARTICLE_DATA__DUST_GRAIN_MASS_OFFSET_);add(_PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_OFFSET_);
  add(_PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_);add(_PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_);
  add(_PIC_PARTICLE_DATA__INJECTION_FACE_OFFSET_);add(_PIC_PARTICLE_DATA__WEIGHT_OVER_TIME_STEP_OFFSET_);
  add(_PIC_PARTICLE_MOVER__RELATIVITY_MODE_);add(_PIC_MOVER_INTEGRATOR_MODE_);
  add(PB::OptionalParticleFieldAllocationManager.MomentumParallelNormal);
  auto datums=[&layout](const std::vector<PIC::Datum::cDatumStored*>& rows) {
    layout.U(rows.size());for(const auto* row:rows) {
      layout.S(row->name);layout.U(row->offset);layout.U(row->length);layout.U(row->type);
    }
  };
  datums(FL::DataStoredAtVertex);datums(FL::cFieldLineSegment::DataStoredAtSegment);
  for(int l=0;l<FL::nFieldLine;++l) {
    auto& line=FL::FieldLinesAll[l];layout.U(line.GetTotalSegmentNumber());
    layout.U(line.GetFirstVertex()->AssociatedDataLength());
    layout.U(line.GetFirstSegment()->AssociatedDataLength());
  }
  std::uint64_t hash=UINT64_C(1469598103934665603);
  for(unsigned char byte:layout.bytes){hash^=byte;hash*=UINT64_C(1099511628211);}
  return hash;
}
T::Status ValidateLiveLists() {
  std::set<long int> seen;
  for(int line=0;line<FL::nFieldLine;++line)
    for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext()) {
      if(segment->tempFirstParticleIndex.load()!=-1)return Error("checkpoint requires completed native moving lists");
      long int previous=-1;
      for(long int particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        if(particle<0||particle>=PB::GetMaxNPart()||segment->Thread!=PIC::ThisThread||
            !seen.insert(particle).second||!PB::IsParticleAllocated(particle)||PB::GetPrev(particle)!=previous)
          return Error("checkpoint live native lists are invalid, cyclic, or owned by another rank");
        previous=particle;
      }
    }
  if(seen.size()!=static_cast<std::size_t>(PB::GetAllPartNum()))
    return Error("checkpoint refuses unlisted or Cartesian native particles");
  return T::Status::Ok();
}
T::Status Capability() {
  // The retained native codec stores the two-branch state and its reservoirs.
  // Once four independently evolving polarizations own the model, aggregate
  // E+/E- diagnostics cannot reconstruct that authority or its provenance.
  // Reject before reading/writing any checkpoint until its complete codec and
  // enclosing collective restore transaction have been implemented/qualified.
  if(Turbulence::Polarized::PICAdapter::HasNativePolarizedState())
    return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
        "physics checkpoint does not yet serialize the live four-polarization wave authority; aggregate two-branch diagnostics cannot represent its complete state");
  if(A::NativeDriverWindowActive())return Error("physics checkpoint requires a completed host step, not an active source birth window");
  if(FluxOwner::HasOpaqueAreaProfile())return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
      "physics checkpoint cannot capture an opaque flux-tube area callback");
  const auto closure=A::ProductionClosurePreflight();if(!closure.ok())return closure;
#if _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_ || _CUDA_MODE_ == _ON_
  return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
      "physics checkpoint currently supports CPU MPI-only execution");
#endif
#if _PIC_FL_MOVER_MULTITHREAD_ == _PIC_MODE_ON_
  return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
      "physics checkpoint first boundary excludes native std::thread field-line workers");
#endif
#if (_SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_ && _SIMULATION_TIME_STEP_MODE_ != _SINGLE_GLOBAL_TIME_STEP_) || _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
      "physics checkpoint requires a declared global timestep mode and species-global particle weight");
#endif
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
      "physics checkpoint does not yet serialize external particle trajectory buffers");
#endif
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||
      SEP::ShockModelType!=SEP::cShockModelType::Analytic1D)
    return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
        "physics checkpoint first boundary supports the standalone analytic provider only");
  const auto& coefficient=T::Coefficient::ActiveConfiguration();
  if(coefficient.pitchAngle==T::Coefficient::PitchAngleKind::Configured||
      coefficient.spatial==T::Coefficient::SpatialKind::ParallelDiffusionLibrary)
    return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
        "physics checkpoint needs typed coefficient models; opaque legacy callbacks/library state are not captured");
  if(!FL::FieldLinesAll||FL::nFieldLine<=0)
    return Error("physics checkpoint requires initialized fixed native field lines");
  for(int line=0;line<FL::nFieldLine;++line)
    if(!FL::FieldLinesAll[line].IsInitialized()||FL::FieldLinesAll[line].is_loop())
      return Error("physics checkpoint first boundary requires initialized open fixed field lines");
  if(!A::PendingWaveContributionsEmpty())return Error("checkpoint requires an empty joined particle-wave queue");
  const auto snapshot=Background::SnapshotStore::Instance().Current();
  if(!snapshot||snapshot->provider()!=Background::Provider::Analytic||
      snapshot->ownership()!=Background::Ownership::ModelOwned)
    return Error("physics checkpoint requires the frozen analytic background realization");
  try {Background::SnapshotStore::Instance().AssertProviderMayWrite(Background::Provider::Analytic);}
  catch(const std::exception& e){return Error(e.what());}
  return ValidateLiveLists();
}

std::string Encode(const Image& image,bool commonOnly=false) {
  Writer out;out.U(kMagic);out.U(kSchema);out.U(0); // output_continuation=false
  out.U(image.purpose); // run_purpose: Production=0, ControlledDiagnostic=1
  out.U(image.completed);out.U(image.ranks);out.U(image.mover);out.D(image.clock);
  out.S(image.runConfiguration);out.S(image.activeConfiguration);out.S(image.initialization);out.S(image.coefficientFingerprint);
  out.S(image.backgroundFingerprint);PutSource(out,image.source);
  const auto& s=*image.snapshot;
  out.U(static_cast<std::uint64_t>(s.provider()));out.U(static_cast<std::uint64_t>(s.ownership()));
  out.D(s.epoch_seconds());out.D(s.valid_from_seconds());out.D(s.valid_until_seconds());
  out.U(s.field_line_generation());out.S(s.configuration_fingerprint());out.S(s.provenance());
  out.D(s.previous_physical_epoch_seconds());out.D(s.current_physical_epoch_seconds());
  out.Doubles(image.steps);out.Doubles(image.weights);out.U(image.stepsInitialized);
  std::string flux;const auto fluxStatus=FieldLine::FluxTubeState::Serialize(image.flux,&flux);
  if(!fluxStatus.ok())throw std::runtime_error(fluxStatus.message);
  out.S(flux);
  if(commonOnly)return out.bytes;
  out.U(image.rank);out.U(image.seed);
  std::ostringstream random;
  random<<std::setprecision(std::numeric_limits<double>::max_digits10)
        <<image.generator<<'\n'<<image.distribution;
  out.S(random.str());out.U(image.hasNormal);out.D(image.normal);
  PutBoundary(out,image.integrated);PutBoundary(out,image.spectral);
  out.U(image.shock.initialized);out.D(image.shock.radiusM);out.D(image.shock.minimumLineRadiusM);
  std::string trajectory;
  if(image.shock.initialized) {
    const auto status=SEP::Shock::SerializeTrajectory(image.shock.trajectory,&trajectory);
    if(!status.ok())throw std::runtime_error(status.message);
  }
  out.S(trajectory);out.U(image.host.adiabaticInitialized);out.U(image.host.shockHistoryInitialized);out.D(image.host.previousShockRadiusM);
  PutTurbulence(out,image.turbulence);out.U(image.lines.size());
  for(const auto& line:image.lines) {
    out.U(line.vertices.size());
    for(std::size_t vertex=0;vertex<line.vertices.size();++vertex) {
      for(double component:line.positions.at(vertex))out.D(component);
      out.S(line.vertices[vertex]);
    }
    out.U(line.segments.size());for(const auto& segment:line.segments) {
      out.U(segment.owner);out.S(segment.data);out.Doubles(segment.weights);
      out.U(segment.particles.size());for(const auto& particle:segment.particles)out.S(particle);
    }
  }
  return out.bytes;
}
Image Decode(const std::string& text) {
  Reader in(text);Image image;
  if(in.U()!=kMagic||in.U()!=kSchema||in.U()!=0)throw std::runtime_error("incompatible corrected physics checkpoint schema");
  image.purpose=in.U();
  image.completed=in.U();image.ranks=in.U();image.mover=in.U();image.clock=in.D();
  image.runConfiguration=in.S();image.activeConfiguration=in.S();image.initialization=in.S();image.coefficientFingerprint=in.S();
  image.backgroundFingerprint=in.S();image.source=GetSource(in);
  const auto provider=in.U(),ownership=in.U();
  const double epoch=in.D(),from=in.D(),until=in.D();const auto generation=in.U();
  const auto fingerprint=in.S(),provenance=in.S();
  const double previous=in.D(),current=in.D();
  image.snapshot.reset(new Background::BackgroundSnapshot(
      static_cast<Background::Provider>(provider),static_cast<Background::Ownership>(ownership),
      epoch,from,until,generation,fingerprint,provenance,previous,current));
  image.steps=in.Doubles();image.weights=in.Doubles();const auto initialized=in.U();
  if(initialized>1)throw std::runtime_error("invalid checkpoint timestep initialized flag");
  image.stepsInitialized=initialized;
  const auto fluxStatus=FieldLine::FluxTubeState::Deserialize(in.S(),&image.flux);
  if(!fluxStatus.ok())throw std::runtime_error(fluxStatus.message);
  image.rank=in.U();const auto seed=in.U();
  if(seed>std::numeric_limits<unsigned long>::max())throw std::runtime_error("checkpoint RNG seed exceeds native width");
  image.seed=seed;std::istringstream random(in.S());random>>image.generator>>image.distribution;
  if(!random)throw std::runtime_error("invalid complete native RNG checkpoint state");
  random>>std::ws;if(!random.eof())throw std::runtime_error("trailing native RNG checkpoint state");
  const auto normal=in.U();if(normal>1)throw std::runtime_error("invalid native normal-cache flag");
  image.hasNormal=normal;image.normal=in.D();image.integrated=GetBoundary(in);image.spectral=GetBoundary(in);
  const auto shock=in.U();if(shock>1)throw std::runtime_error("invalid shock cache initialized flag");
  image.shock.initialized=shock;image.shock.radiusM=in.D();image.shock.minimumLineRadiusM=in.D();
  const auto trajectory=in.S();
  if(shock) {
    const auto status=SEP::Shock::DeserializeTrajectory(trajectory,&image.shock.trajectory);
    if(!status.ok())throw std::runtime_error(status.message);
  } else if(!trajectory.empty())throw std::runtime_error("uninitialized shock cache has hidden trajectory state");
  const auto adiabatic=in.U(),history=in.U();
  if(adiabatic>1||history>1)throw std::runtime_error("invalid host-step checkpoint flag");
  image.host.adiabaticInitialized=adiabatic;image.host.shockHistoryInitialized=history;image.host.previousShockRadiusM=in.D();
  image.turbulence=GetTurbulence(in);
  const auto lines=in.Count(16);image.lines.resize(lines);
  for(auto& line:image.lines) {
    const auto vertices=in.Count(32);line.vertices.resize(vertices);line.positions.resize(vertices);
    for(std::size_t vertex=0;vertex<vertices;++vertex) {
      for(auto& component:line.positions[vertex])component=in.D();
      line.vertices[vertex]=in.S();
    }
    const auto segments=in.Count(32);line.segments.resize(segments);
    for(auto& segment:line.segments) {
      segment.owner=in.U();segment.data=in.S();segment.weights=in.Doubles();
      const auto particles=in.Count(8);segment.particles.resize(particles);
      for(auto& particle:segment.particles)particle=in.S();
    }
  }
  in.End();return image;
}

template<class Value> Value Word(const std::string& particle,int offset) {
  Value value;std::memcpy(&value,particle.data()+offset,sizeof(value));return value;
}
T::Status ReadVertexPrimitive(const std::string& buffer,
    const PIC::Datum::cDatumStored& datum,std::size_t components,double* output) {
  // Associated buffers carry doubles, not native objects. Decode through
  // memcpy so a byte string is never dereferenced as an aligned double array.
  // The complete layout fingerprint has already pinned datum meaning/units.
  if(datum.offset<0||datum.length<0||static_cast<std::size_t>(datum.length)!=components||
      static_cast<std::size_t>(datum.offset)>buffer.size()||
      components>(buffer.size()-static_cast<std::size_t>(datum.offset))/sizeof(double))
    return Error("checkpoint current native vertex primitive layout is incompatible");
  std::memcpy(output,buffer.data()+datum.offset,components*sizeof(double));
  return T::Status::Ok();
}
T::Status ValidateCurrentVertexPrimitives(const std::string& buffer) {
  double magneticT[3],velocityMPerS[3],numberDensityPerM3=0.0;
  auto status=ReadVertexPrimitive(buffer,FL::DatumAtVertexMagneticField,3,magneticT);
  if(!status.ok())return status;
  status=ReadVertexPrimitive(buffer,FL::DatumAtVertexPlasmaVelocity,3,velocityMPerS);
  if(!status.ok())return status;
  status=ReadVertexPrimitive(buffer,FL::DatumAtVertexPlasmaDensity,1,&numberDensityPerM3);
  if(!status.ok())return status;
  if(!std::isfinite(numberDensityPerM3)||numberDensityPerM3<=0.0)
    return Error("checkpoint current native vertex density must be finite and positive [m^-3]");
  const double magneticMagnitudeT=std::hypot(std::hypot(magneticT[0],magneticT[1]),magneticT[2]);
  if(!std::isfinite(magneticT[0])||!std::isfinite(magneticT[1])||!std::isfinite(magneticT[2])||
      !std::isfinite(magneticMagnitudeT)||magneticMagnitudeT<=0.0)
    return Error("checkpoint current native vertex magnetic field must be finite and nonzero [T]");
  const double plasmaSpeedMPerS=std::hypot(std::hypot(velocityMPerS[0],velocityMPerS[1]),velocityMPerS[2]);
  if(!std::isfinite(velocityMPerS[0])||!std::isfinite(velocityMPerS[1])||!std::isfinite(velocityMPerS[2])||
      !std::isfinite(plasmaSpeedMPerS)||plasmaSpeedMPerS>=SpeedOfLight)
    return Error("checkpoint current native vertex plasma velocity must be finite and subluminal [m/s]");
  // This codec is restricted to the analytic authority. Its frozen mean
  // plasma atomic mass is part of ActualConfiguration(), so use the same
  // declared conversion and Alfven domain as the native background adapter.
  // A positive byte value whose mass density underflows is not usable state.
  const double massDensityKgM3=numberDensityPerM3*PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
  const double magneticDenominator=VacuumPermeability*massDensityKgM3;
  if(!std::isfinite(massDensityKgM3)||massDensityKgM3<=0.0||
      !std::isfinite(magneticDenominator)||magneticDenominator<=0.0)
    return Error("checkpoint current native vertex density/composition has no representable positive Alfven mass density [kg/m^3]");
  const double alfvenMPerS=magneticMagnitudeT/std::sqrt(magneticDenominator);
  if(!std::isfinite(alfvenMPerS)||alfvenMPerS<=0.0||alfvenMPerS>=SpeedOfLight)
    return Error("checkpoint current native vertex Alfven speed is outside the finite subluminal MHD domain [m/s]");
  return T::Status::Ok();
}
T::Status ValidateLineFluxGeometry(const Image& image,int lineIndex) {
  const auto& saved=image.lines[lineIndex];double first[3];
  auto status=ReadVertexPrimitive(saved.vertices[0],FL::DatumAtVertexMagneticField,3,first);
  if(!status.ok())return status;
  const auto magnitude=[](const double* b){return std::hypot(std::hypot(b[0],b[1]),b[2]);};
  double flux=0.0;
  if(static_cast<std::size_t>(lineIndex)<image.flux.records.size()&&
      !FieldLine::FluxTubeState::Uninitialized(image.flux.records[lineIndex]))
    flux=image.flux.records[lineIndex].magneticFluxWb;
  else flux=image.flux.referenceAreaM2*magnitude(first);
  if(!std::isfinite(flux)||flux<=0.0)return Error("checkpoint flux cannot represent positive tube normalization [Wb]");
  for(std::size_t s=0;s<saved.segments.size();++s) {
    double begin[3],end[3],middle[3];
    status=ReadVertexPrimitive(saved.vertices[s],FL::DatumAtVertexMagneticField,3,begin);if(!status.ok())return status;
    status=ReadVertexPrimitive(saved.vertices[s+1],FL::DatumAtVertexMagneticField,3,end);if(!status.ok())return status;
    for(int d=0;d<3;++d)middle[d]=0.5*begin[d]+0.5*end[d];
    const double a=flux/magnitude(begin),b=flux/magnitude(middle),c=flux/magnitude(end);
    if(!std::isfinite(a)||a<=0.0||!std::isfinite(b)||b<=0.0||!std::isfinite(c)||c<=0.0)
      return Error("checkpoint conserved flux/background cannot represent finite positive tube areas [m^2]");
    // Validate the production Simpson-volume representation against decoded B,
    // not whichever lazy flux table happens to be warm in this process. This
    // is observational: no missing record is materialized while preflighting.
    const double volume=FieldLine::FluxTubeGeometryCore::IntegrateSimpsonArea(
        Units::AreaM2(a),Units::AreaM2(b),Units::AreaM2(c),
        Units::LengthM(FL::FieldLinesAll[lineIndex].GetSegment(static_cast<int>(s))->GetLength())).Value();
    if(!std::isfinite(volume)||volume<=0.0)return Error("checkpoint tube volume is not finite and positive [m^3]");
  }
  return T::Status::Ok();
}
T::Status ValidateRecord(const std::string& record,int line,int segment,
    const std::array<double,3>& begin,const std::array<double,3>& end,
    std::uint64_t campaign,std::set<std::uint64_t>* ids) {
  if(record.size()!=static_cast<std::size_t>(PB::ParticleDataLength))return Error("incompatible checkpoint particle record length");
  auto* data=reinterpret_cast<PB::byte*>(const_cast<char*>(record.data()));
  const int species=PB::GetI(data);const double coordinate=PB::GetFieldLineCoord(data);
  const double parallel=PB::GetVParallel(data),normal=PB::GetVNormal(data);
  if(!PB::IsParticleAllocated(data)||species<0||species>=PIC::nTotalSpecies||
      PB::GetFieldLineId(data)!=line||!std::isfinite(coordinate)||coordinate<0.0||
      coordinate>=FL::FieldLinesAll[line].GetTotalSegmentNumber()||
      FL::FieldLinesAll[line].GetSegment(coordinate)!=FL::FieldLinesAll[line].GetSegment(segment)||
      !std::isfinite(parallel)||!std::isfinite(normal)||normal<0.0||
      !(std::hypot(parallel,normal)<SpeedOfLight)||
      !(PB::GetIndividualStatWeightCorrection(data)>0.0)||
      !std::isfinite(PB::GetIndividualStatWeightCorrection(data)))
    return Error("checkpoint particle violates native field-line phase-space contract");
  double position[3];PB::GetX(position,data);
  const double fraction=coordinate-static_cast<int>(coordinate);
  for(unsigned component=0;component<3;++component) {
    // This is the pinned native piecewise-linear GetCartesian map, evaluated
    // on DECODED endpoints. A finite but stale X is as invalid as NaN: generic
    // PIC consumers may read it before the first post-restart mover call.
    const double expected=(1.0-fraction)*begin[component]+fraction*end[component];
    if(!std::isfinite(expected)||!std::isfinite(position[component])||position[component]!=expected)
      return Error("checkpoint particle Cartesian view is nonfinite or inconsistent with decoded field-line geometry");
  }
  const auto id=Word<std::uint64_t>(record,SEP::Offset::StableParticleId);
  if(Word<std::uint64_t>(record,SEP::Offset::TransportSchema)!=A::ParticleTransportSchema||
      id==0||!ids->insert(id).second)return Error("checkpoint particle has incompatible schema or duplicate identity");
  T::Birth::Record birth;std::memcpy(&birth,record.data()+SEP::Offset::BirthRecord,sizeof(birth));
  const auto birthStatus=T::Birth::ValidateRecord(birth);if(!birthStatus.ok())return birthStatus;
  if(birth.pendingFirstMove)return Error("completed physics checkpoint cannot contain an unconsumed scheduled first move");
  T::OperatorRandomState states[3];std::memcpy(states,record.data()+SEP::Offset::OperatorRandomStates,sizeof(states));
  for(int operation=0;operation<3;++operation) {
    T::KeyedRandomStream expected(campaign,id,7+operation,0);
    const auto status=expected.RestoreState(states[operation].stream);if(!status.ok())return status;
  }
  const double hazard=Word<double>(record,SEP::Offset::MfpOpticalDepth);
  const auto event=Word<std::uint64_t>(record,SEP::Offset::MfpEventIndex);
  if((std::isnan(hazard)&&(event!=0||states[2].acceptedIntervals!=0))||
      (!std::isnan(hazard)&&(!std::isfinite(hazard)||hazard<=0.0)))
    return Error("checkpoint particle has lost residual optical-depth history");
  return T::Status::Ok();
}

T::Status Capture(std::uint64_t completed,Image* image) {
  const auto allowed=Capability();if(!allowed.ok())return allowed;
  image->completed=completed;image->rank=PIC::ThisThread;image->ranks=PIC::nTotalThreads;
  image->purpose=static_cast<std::uint64_t>(A::CurrentRunPurpose());
  image->mover=static_cast<std::uint64_t>(SEP::Mover::CurrentProductionMover());
  image->clock=Background::SimulationTimeSeconds();
  auto status=Run::Serialize(Run::Active(),&image->runConfiguration);if(!status.ok())return status;
  image->activeConfiguration=ActualConfiguration();
  image->initialization=Initialization::HasActive()?Initialization::Fingerprint(Initialization::Active()):"legacy-fixed-native-initialization";
  image->coefficientFingerprint=T::Coefficient::ConfigurationFingerprint(T::Coefficient::ActiveConfiguration());
  image->backgroundFingerprint=Background::CurrentConfigurationFingerprint();
  image->source=SourceLifecycle::CaptureNativeState();
  status=SourceLifecycle::ValidateNativeState(image->source);if(!status.ok())return status;
  image->snapshot=Background::SnapshotStore::Instance().Current();
  image->steps.assign(PIC::ParticleWeightTimeStep::GlobalTimeStep,PIC::ParticleWeightTimeStep::GlobalTimeStep+PIC::nTotalSpecies);
  image->weights.assign(PIC::ParticleWeightTimeStep::GlobalParticleWeight,PIC::ParticleWeightTimeStep::GlobalParticleWeight+PIC::nTotalSpecies);
  image->stepsInitialized=PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized;
  image->seed=RandomNumberGenerator::rndLastSeed;image->generator=RandomNumberGenerator::gen;
  image->distribution=RandomNumberGenerator::dist;image->hasNormal=RandomNumberGenerator::NormalSpareAvailable;
  image->normal=RandomNumberGenerator::NormalSpare;
  image->integrated=W::CaptureIntegratedBoundaryCache();image->spectral=W::WaveNumberResolved::CaptureSpectralBoundaryCache();
  status=ShockOwner::CaptureAnalyticShockCache(&image->shock);if(!status.ok())return status;
  image->host=CaptureHostStepCache();image->turbulence=TurbOwner::CaptureCheckpointCache();
  status=FluxOwner::CaptureCheckpointState(&image->flux);if(!status.ok())return status;
  image->lines.resize(FL::nFieldLine);
  std::set<long int> pointers;std::set<std::uint64_t> ids;
  for(int l=0;l<FL::nFieldLine;++l) {
    auto& line=FL::FieldLinesAll[l];auto& saved=image->lines[l];const int count=line.GetTotalSegmentNumber();
    saved.vertices.resize(count+1);saved.positions.resize(count+1);saved.segments.resize(count);
    for(int v=0;v<=count;++v) {
      auto* vertex=line.GetVertex(v);const int bytes=vertex->AssociatedDataLength();
      if(bytes<0||(bytes&&!vertex->GetAssociatedDataBufferPointer()))return Error("checkpoint vertex associated state is absent");
      saved.vertices[v].assign(vertex->GetAssociatedDataBufferPointer(),bytes);
      vertex->GetX(saved.positions[v].data());
      for(double component:saved.positions[v])if(!std::isfinite(component))
        return Error("checkpoint native vertex position must be finite [m]");
    }
    for(int s=0;s<count;++s) {
      auto* segment=line.GetSegment(s);auto& row=saved.segments[s];row.owner=segment->Thread;
      const int bytes=segment->AssociatedDataLength();
      if(bytes<0||(bytes&&!segment->GetAssociatedDataBufferPointer()))return Error("checkpoint segment associated state is absent");
      row.data.assign(segment->GetAssociatedDataBufferPointer(),bytes);
      for(int species=0;species<PIC::nTotalSpecies;++species)row.weights.push_back(segment->GetWeight(species));
      if(segment->tempFirstParticleIndex.load()!=-1)return Error("checkpoint encountered unfinished native moving lists");
      for(long int p=segment->FirstParticleIndex;p!=-1;p=PB::GetNext(p)) {
        if(segment->Thread!=PIC::ThisThread||!pointers.insert(p).second)return Error("checkpoint encountered non-owner or cyclic native particle list");
        std::string record(reinterpret_cast<const char*>(PB::GetParticleDataPointer(p)),PB::ParticleDataLength);
        status=ValidateRecord(record,l,s,saved.positions[s],saved.positions[s+1],A::CampaignRandomSeed,&ids);if(!status.ok())return status;
        // Links are allocation bookkeeping. Preserve semantic list order and
        // rebuild both directions through GetNewParticle/CloneParticle on read.
        auto* raw=reinterpret_cast<PB::byte*>(&record[0]);PB::SetNext(-1,raw);PB::SetPrev(-1,raw);
        row.particles.push_back(std::move(record));
      }
    }
  }
  if(pointers.size()!=static_cast<std::size_t>(PB::GetAllPartNum()))return Error("checkpoint refuses unlisted or Cartesian native particles");
  return T::Status::Ok();
}

T::Status Validate(const Image& image,bool replacing=false) {
  const auto allowed=Capability();if(!allowed.ok())return allowed;
  if(image.purpose!=static_cast<std::uint64_t>(A::CurrentRunPurpose()))
    return Error("checkpoint run purpose differs; controlled diagnostics cannot become production realizations");
  if(image.rank!=static_cast<std::uint64_t>(PIC::ThisThread)||image.ranks!=static_cast<std::uint64_t>(PIC::nTotalThreads)||
      image.mover!=static_cast<std::uint64_t>(SEP::Mover::CurrentProductionMover()))
    return Error("physics checkpoint requires identical MPI rank count, ownership, and selected mover");
  std::string run;auto status=Run::Serialize(Run::Active(),&run);if(!status.ok())return status;
  if(image.runConfiguration!=run||image.activeConfiguration!=ActualConfiguration()||image.coefficientFingerprint!=T::Coefficient::ConfigurationFingerprint(T::Coefficient::ActiveConfiguration())||
      image.backgroundFingerprint!=Background::CurrentConfigurationFingerprint()||
      image.initialization!=(Initialization::HasActive()?Initialization::Fingerprint(Initialization::Active()):"legacy-fixed-native-initialization"))
    return Error("physics checkpoint frozen configuration differs from requested run");
  status=SourceLifecycle::ValidateNativeState(image.source);if(!status.ok())return status;
  if(!std::isfinite(image.clock)||!image.snapshot||!image.snapshot->Covers(image.clock)||
      image.snapshot->provider()!=Background::Provider::Analytic||image.snapshot->ownership()!=Background::Ownership::ModelOwned||
      image.snapshot->configuration_fingerprint()!=image.backgroundFingerprint||!std::isfinite(image.normal))
    return Error("physics checkpoint clock/background/RNG history is invalid");
  try {Background::SnapshotStore::Instance().ValidateCheckpointRestore(*image.snapshot);}
  catch(const std::exception& e){return Error(e.what());}
  if(image.steps.size()!=static_cast<std::size_t>(PIC::nTotalSpecies)||image.weights.size()!=image.steps.size()||!image.stepsInitialized)
    return Error("physics checkpoint global timestep/weight layout is invalid");
  for(std::size_t i=0;i<image.steps.size();++i) {
    // A single-global build consumes table[0] for every species. Other table
    // entries may retain core sentinels and are preserved byte-for-byte rather
    // than promoted to active timesteps. Species-global consumes every entry.
    const bool activeStep=_SIMULATION_TIME_STEP_MODE_!=_SINGLE_GLOBAL_TIME_STEP_||i==0;
    if((activeStep&&(!std::isfinite(image.steps[i])||image.steps[i]<=0.0))||
        !std::isfinite(image.weights[i])||image.weights[i]<=0.0)
      return Error("physics checkpoint global timestep/weight is nonphysical");
  }
  status=W::ValidateIntegratedBoundaryCache(image.integrated);if(!status.ok())return status;
  status=W::WaveNumberResolved::ValidateSpectralBoundaryCache(image.spectral);if(!status.ok())return status;
  status=ShockOwner::ValidateAnalyticShockCache(image.shock);if(!status.ok())return status;
  status=ValidateHostStepCache(image.host);if(!status.ok())return status;
  status=TurbOwner::ValidateCheckpointCache(image.turbulence);if(!status.ok())return status;
  status=FluxOwner::ValidateCheckpointState(image.flux);if(!status.ok())return status;
  if(image.lines.size()!=static_cast<std::size_t>(FL::nFieldLine))return Error("checkpoint line topology differs from native geometry");
  std::size_t particles=0;std::set<std::uint64_t> ids;
  for(int l=0;l<FL::nFieldLine;++l) {
    auto& line=FL::FieldLinesAll[l];const auto& saved=image.lines[l];const int count=line.GetTotalSegmentNumber();
    if(saved.vertices.size()!=static_cast<std::size_t>(count+1)||saved.positions.size()!=saved.vertices.size()||saved.segments.size()!=static_cast<std::size_t>(count))
      return Error("checkpoint vertex/segment topology differs from native geometry");
    for(int v=0;v<=count;++v) {
      const double* current=line.GetVertex(v)->GetX();
      for(unsigned component=0;component<3;++component)
        if(!std::isfinite(saved.positions[v][component])||saved.positions[v][component]!=current[component])
          return Error("checkpoint decoded vertex position differs from the supported fixed physical geometry");
      if(saved.vertices[v].size()!=static_cast<std::size_t>(line.GetVertex(v)->AssociatedDataLength()))return Error("checkpoint vertex data layout differs");
      // Preflight decoded primitives before PrepareReplacement allocates any
      // PIC record and before Commit can overwrite native vertex storage.
      // Thermal/previous-epoch slots are preserved; this first analytic codec
      // does not invent a missing thermal closure or use old density as div U.
      status=ValidateCurrentVertexPrimitives(saved.vertices[v]);if(!status.ok())return status;
    }
    status=ValidateLineFluxGeometry(image,l);if(!status.ok())return status;
    for(int s=0;s<count;++s) {
      const auto& row=saved.segments[s];auto* segment=line.GetSegment(s);
      if(row.owner!=static_cast<std::uint64_t>(segment->Thread)||row.data.size()!=static_cast<std::size_t>(segment->AssociatedDataLength())||
          row.weights.size()!=static_cast<std::size_t>(PIC::nTotalSpecies)||
          (!row.particles.empty()&&segment->Thread!=PIC::ThisThread))
        return Error("checkpoint segment ownership or native data layout differs");
      for(double weight:row.weights)if(!std::isfinite(weight)||weight<0.0)return Error("invalid checkpoint segment source weight");
      for(const auto& record:row.particles) {
        status=ValidateRecord(record,l,s,saved.positions[s],saved.positions[s+1],A::CampaignRandomSeed,&ids);if(!status.ok())return status;++particles;
      }
    }
  }
  const auto occupied=PB::GetAllPartNum(),capacity=PB::GetMaxNPart();
  if(occupied<0||capacity<occupied||particles>static_cast<std::size_t>(capacity)||
      (replacing&&particles>static_cast<std::size_t>(capacity-occupied)))
    return Error("checkpoint restore needs spare native capacity for a fully staged replacement population");
  return T::Status::Ok();
}

T::Status Agreed(const T::Status& local) {
  const int candidate=local.ok()?INT_MAX:PIC::ThisThread;int first=0;
  MPI_Allreduce(&candidate,&first,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(first==INT_MAX)return T::Status::Ok();
  std::string message=PIC::ThisThread==first?local.message:std::string();
  int length=static_cast<int>(message.size());MPI_Bcast(&length,1,MPI_INT,first,MPI_GLOBAL_COMMUNICATOR);
  int code=PIC::ThisThread==first?static_cast<int>(local.code):0;
  MPI_Bcast(&code,1,MPI_INT,first,MPI_GLOBAL_COMMUNICATOR);
  message.resize(length);MPI_Bcast(length?&message[0]:nullptr,length,MPI_CHAR,first,MPI_GLOBAL_COMMUNICATOR);
  return T::Status::Error(static_cast<T::StatusCode>(code),
      "collective checkpoint rejected on rank "+std::to_string(first)+": "+message);
}

T::Status AgreeText(const std::string& local,const char* failure) {
  auto status=Agreed(local.size()>INT_MAX?Error("checkpoint agreement text exceeds MPI limit"):T::Status::Ok());
  if(!status.ok())return status;
  int size=PIC::ThisThread==0?static_cast<int>(local.size()):0;
  MPI_Bcast(&size,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  std::string expected;T::Status allocated;
  try {expected=PIC::ThisThread==0?local:std::string(size,'\0');}
  catch(const std::exception& e){allocated=Error(e.what());}
  status=Agreed(allocated);if(!status.ok())return status;
  MPI_Bcast(size?&expected[0]:nullptr,size,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  return Agreed(local==expected?T::Status::Ok():Error(failure));
}
T::Status AgreeLayout() {
  std::uint64_t local=0;T::Status valid;
  try {local=Layout();}catch(const std::exception& e){valid=Error(e.what());}
  auto status=Agreed(valid);if(!status.ok())return status;
  std::uint64_t expected=PIC::ThisThread==0?local:0;
  MPI_Bcast(&expected,1,MPI_UINT64_T,0,MPI_GLOBAL_COMMUNICATOR);
  return Agreed(local==expected?T::Status::Ok():Error("MPI owners have different native checkpoint layouts"));
}
T::Status UniqueIdentities(const Image& image) {
  std::vector<std::uint64_t> local;T::Status allocated;
  try {
    for(const auto& line:image.lines)for(const auto& segment:line.segments)
      for(const auto& record:segment.particles)local.push_back(Word<std::uint64_t>(record,SEP::Offset::StableParticleId));
    if(local.size()>INT_MAX)allocated=Error("checkpoint identities exceed MPI count limit");
  }catch(const std::exception& e){allocated=Error(e.what());}
  auto status=Agreed(allocated);
  if(!status.ok())return status;
  const int count=static_cast<int>(local.size());std::vector<int> counts,offsets;
  try {counts.resize(PIC::nTotalThreads);offsets.resize(PIC::nTotalThreads);}
  catch(const std::exception& e){allocated=Error(e.what());}
  status=Agreed(allocated);if(!status.ok())return status;
  MPI_Allgather(&count,1,MPI_INT,counts.data(),1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
  std::size_t total=0;bool overflow=false;
  for(int rank=0;rank<PIC::nTotalThreads;++rank) {
    if(counts[rank]<0||total>static_cast<std::size_t>(INT_MAX-counts[rank])){overflow=true;break;}
    offsets[rank]=static_cast<int>(total);total+=counts[rank];
  }
  status=Agreed(overflow?Error("global checkpoint identities exceed MPI count limit"):T::Status::Ok());if(!status.ok())return status;
  std::vector<std::uint64_t> all;
  try {all.resize(total);}catch(const std::exception& e){allocated=Error(e.what());}
  status=Agreed(allocated);if(!status.ok())return status;
  MPI_Allgatherv(local.data(),count,MPI_UINT64_T,all.data(),counts.data(),offsets.data(),MPI_UINT64_T,MPI_GLOBAL_COMMUNICATOR);
  std::sort(all.begin(),all.end());
  return Agreed(std::adjacent_find(all.begin(),all.end())==all.end()?T::Status::Ok():Error("global checkpoint particle identity collision"));
}

T::Status PrepareReplacement(Image* image) {
  auto status=SourceLifecycle::PrepareNativeRestore(image->source,&image->stagedSource);
  if(!status.ok())return status;
  status=FluxOwner::StageCheckpointRestore(image->flux,&image->stagedFlux);if(!status.ok())return status;
  for(auto& line:image->lines)for(auto& segment:line.segments)
    for(auto record=segment.particles.rbegin();record!=segment.particles.rend();++record) {
      const auto particle=PB::GetNewParticle(segment.stagedHead);
      PB::CloneParticle(PB::GetParticleDataPointer(particle),reinterpret_cast<PB::byte*>(&(*record)[0]));
      PB::SetParticleAllocated(particle);
    }
  return T::Status::Ok();
}
void DiscardReplacement(Image* image) {
  for(auto& line:image->lines)for(auto& segment:line.segments) {
    auto particle=segment.stagedHead;
    while(particle!=-1){const auto next=PB::GetNext(particle);PB::DeleteParticle(particle);particle=next;}
    segment.stagedHead=-1;
  }
}
void Commit(Image* image) {
  // Every malformed/incompatible-state check precedes the all-rank agreement.
  // Owner caches were allocated while staging. Allocation links are rebuilt by
  // safe PIC APIs; raw pointer/atomic/list-head objects are never deserialized.
  if(!FluxOwner::CommitCheckpointRestore(&image->stagedFlux))
    throw std::runtime_error("flux owner became stale after the collective checkpoint barrier");
  for(int l=0;l<FL::nFieldLine;++l)for(auto* segment=FL::FieldLinesAll[l].GetFirstSegment();segment;segment=segment->GetNext()) {
    long int particle=segment->FirstParticleIndex;
    while(particle!=-1){const auto next=PB::GetNext(particle);PB::DeleteParticle(particle);particle=next;}
    segment->FirstParticleIndex=-1;segment->tempFirstParticleIndex=-1;
  }
  for(int l=0;l<FL::nFieldLine;++l) {
    auto& line=FL::FieldLinesAll[l];auto& saved=image->lines[l];
    for(std::size_t v=0;v<saved.vertices.size();++v)
      std::memcpy(line.GetVertex(static_cast<int>(v))->GetAssociatedDataBufferPointer(),saved.vertices[v].data(),saved.vertices[v].size());
    for(std::size_t s=0;s<saved.segments.size();++s) {
      auto* segment=line.GetSegment(static_cast<int>(s));auto& row=saved.segments[s];
      std::memcpy(segment->GetAssociatedDataBufferPointer(),row.data.data(),row.data.size());
      for(int species=0;species<PIC::nTotalSpecies;++species)segment->SetWeight(row.weights[species],species);
      segment->FirstParticleIndex=row.stagedHead;row.stagedHead=-1;
    }
  }
  std::copy(image->steps.begin(),image->steps.end(),PIC::ParticleWeightTimeStep::GlobalTimeStep);
  std::copy(image->weights.begin(),image->weights.end(),PIC::ParticleWeightTimeStep::GlobalParticleWeight);
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=image->stepsInitialized;
  RandomNumberGenerator::rndLastSeed=image->seed;RandomNumberGenerator::gen=image->generator;
  RandomNumberGenerator::dist=image->distribution;RandomNumberGenerator::NormalSpareAvailable=image->hasNormal;RandomNumberGenerator::NormalSpare=image->normal;
  W::RestoreIntegratedBoundaryCache(&image->integrated);W::WaveNumberResolved::RestoreSpectralBoundaryCache(&image->spectral);
  ShockOwner::RestoreAnalyticShockCache(&image->shock);RestoreHostStepCache(image->host);TurbOwner::RestoreCheckpointCache(&image->turbulence);
  SourceLifecycle::CommitNativeRestore(&image->stagedSource);
  PIC::SimulationTime::SetInitialValue(image->clock);
  Background::SnapshotStore::Instance().RestoreCheckpoint(image->snapshot);
}

}  // namespace

T::Status LocalPhysicsImage(std::uint64_t completed,std::string* text) {
  if(!text)return Error("null physics image destination");
  try {Image image;auto status=Capture(completed,&image);if(!status.ok())return status;
    status=Validate(image);if(!status.ok())return status;*text=Encode(image);return T::Status::Ok();}
  catch(const std::exception& e){return Error(e.what());}
}

T::Status SaveCollective(const std::string& path,std::uint64_t completed) {
  T::Status local;std::string bytes,common;Image image;
  try {
    if(path.empty())local=Error("empty physics checkpoint destination");
    else {local=Capture(completed,&image);if(local.ok())local=Validate(image);
      if(local.ok()){bytes=Encode(image);common=Encode(image,true);}}
    if(bytes.size()>INT_MAX||common.size()>INT_MAX)local=Error("checkpoint rank payload exceeds first MPI codec size limit");
  } catch(const std::exception& e){local=Error(e.what());}
  auto status=Agreed(local);if(!status.ok())return status;
  status=AgreeText(path,"ranks disagree on checkpoint destination");if(!status.ok())return status;
  status=AgreeLayout();if(!status.ok())return status;
  status=AgreeText(common,"ranks disagree on physics epoch/configuration/source scheduler");if(!status.ok())return status;
  status=UniqueIdentities(image);if(!status.ok())return status;
  const int count=static_cast<int>(bytes.size());std::vector<int> counts,offsets;
  try {counts.resize(PIC::nTotalThreads);offsets.resize(PIC::nTotalThreads);}
  catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);if(!status.ok())return status;
  MPI_Allgather(&count,1,MPI_INT,counts.data(),1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
  std::uint64_t total=0;bool overflow=false;
  for(int rank=0;rank<PIC::nTotalThreads;++rank) {
    if(counts[rank]<0||total>static_cast<std::uint64_t>(INT_MAX-counts[rank])){overflow=true;break;}
    offsets[rank]=static_cast<int>(total);total+=counts[rank];
  }
  const std::uint64_t framing=24+8*static_cast<std::uint64_t>(PIC::nTotalThreads);
  status=Agreed(overflow||framing>INT_MAX||total>INT_MAX-framing?
      Error("checkpoint total payload exceeds first MPI codec size limit"):T::Status::Ok());if(!status.ok())return status;
  std::string gathered;
  try {if(PIC::ThisThread==0)gathered.resize(static_cast<std::size_t>(total),'\0');}
  catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);if(!status.ok())return status;
  MPI_Gatherv(count?&bytes[0]:nullptr,count,MPI_CHAR,PIC::ThisThread==0?&gathered[0]:nullptr,
      counts.data(),offsets.data(),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread==0)try {
    Writer payload;payload.U(kMagic);payload.U(kSchema);payload.U(PIC::nTotalThreads);
    for(int rank=0;rank<PIC::nTotalThreads;++rank)payload.S(gathered.substr(offsets[rank],counts[rank]));
    std::vector<unsigned char> encoded(payload.bytes.begin(),payload.bytes.end());
    const auto written=PIC::PersistentFormat::WriteAtomic(path,PIC::PersistentFormat::Kind::ParticleRestart,
        PIC::PersistentFormat::NativeLayout,Layout(),encoded);
    local=written==PIC::PersistentFormat::Status::Success?T::Status::Ok():Error("atomic corrected physics checkpoint publication failed");
  } catch(const std::exception& e){local=Error(e.what());}
  return Agreed(local);
}

T::Status LoadCollective(const std::string& path,std::uint64_t* completed) {
  T::Status local;
  try {local=completed?Capability():Error("null restored step destination");}
  catch(const std::exception& e){local=Error(e.what());}
  auto status=Agreed(local);if(!status.ok())return status;
  status=AgreeText(path,"ranks disagree on checkpoint source");if(!status.ok())return status;
  status=AgreeLayout();if(!status.ok())return status;
  Image staged;
  std::string loadedBytes;
  if(PIC::ThisThread==0)try {
    FILE* file=nullptr;std::uint64_t size=0;
    const auto opened=PIC::PersistentFormat::OpenValidated(path,PIC::PersistentFormat::Kind::ParticleRestart,Layout(),false,file,size);
    if(opened!=PIC::PersistentFormat::Status::Success)local=Error("incompatible, corrupt, missing, or legacy physics checkpoint");
    else if(size>INT_MAX){std::fclose(file);local=Error("physics checkpoint exceeds first collective codec size limit");}
    else {
      loadedBytes.resize(static_cast<std::size_t>(size),'\0');
      const bool read=size==0||std::fread(&loadedBytes[0],1,loadedBytes.size(),file)==loadedBytes.size();std::fclose(file);
      if(!read)local=Error("cannot read validated physics checkpoint payload");
    }
  } catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);if(!status.ok())return status;
  // Only rank zero touches the checkpoint filesystem. Broadcast the one
  // CRC-validated payload, so a multi-node run does not assume node-local /tmp
  // or a checkpoint path is mounted identically on every MPI owner.
  int size=static_cast<int>(loadedBytes.size());MPI_Bcast(&size,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  try {if(PIC::ThisThread!=0)loadedBytes.resize(size);}
  catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);if(!status.ok())return status;
  MPI_Bcast(size?&loadedBytes[0]:nullptr,size,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  try {
    Reader input(loadedBytes);
    if(input.U()!=kMagic||input.U()!=kSchema||input.U()!=static_cast<std::uint64_t>(PIC::nTotalThreads))
      throw std::runtime_error("incompatible corrected physics schema or MPI rank count");
    std::string selected;
    for(int rank=0;rank<PIC::nTotalThreads;++rank) {
      auto section=input.S();if(rank==PIC::ThisThread)selected=std::move(section);
    }
    input.End();staged=Decode(selected);local=Validate(staged,true);
  }catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);if(!status.ok())return status;
  std::string common;
  try {common=Encode(staged,true);}catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);if(!status.ok())return status;
  status=AgreeText(common,"ranks disagree on restored replicated physics/source/flux state");if(!status.ok())return status;
  status=UniqueIdentities(staged);if(!status.ok())return status;
  try {local=PrepareReplacement(&staged);
    if(local.ok())local=FluxOwner::ValidateCheckpointReady(staged.stagedFlux);
  }catch(const std::exception& e){local=Error(e.what());}
  status=Agreed(local);
  if(!status.ok()){DiscardReplacement(&staged);return status;}
  // All ranks have fully decoded and validated before any old particle or
  // provider state is touched. Commit is quiescent; failures after this point
  // indicate a violated host ownership/resource invariant and are fatal.
  try {Commit(&staged);}catch(const std::exception& e) {
    exit(__LINE__,__FILE__,("physics checkpoint commit invariant failed: "+std::string(e.what())).c_str());
  }
  *completed=staged.completed;MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  return T::Status::Ok();
}

} }  // namespace SEP::PhysicsCheckpoint
