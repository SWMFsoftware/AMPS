#include "sep_maxwellian_native_validation.h"
#include "../transport_common.h"
#include "sep_source_lifecycle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <set>
#include <vector>

namespace SEP { namespace Testing { namespace {
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace T=Transport;
namespace A=T::PICAdapter;

// The legacy multiplicative scalar generator has modulus 2^31. Even seeds
// inhabit shorter subcycles (2^30 is a fixed point), so the test-only bridge
// selects its positive odd residue class explicitly. This finite projection
// is not a claim that the legacy sampler has the production keyed-stream law.
unsigned long LegacyScalarSeed(std::uint64_t bits) {
  return static_cast<unsigned long>((bits&UINT64_C(0x7fffffff))|UINT64_C(1));
}

// Only this fixture's actual native records are deleted. The empty-population
// prerequisite avoids detaching unrelated records or confusing physical source
// state with a convenience duplicate. Geometry and core RNG are restored even
// if the alias regression fires. No mover or production source is certified.
struct Fixture {
  struct Point {FL::cFieldLineVertex* vertex;std::array<double,3> x;};
  std::vector<Point> points;
  std::vector<double> steps,weights;
  bool initialized=PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized;
  unsigned long nativeSeed=RandomNumberGenerator::rndLastSeed;
  std::mt19937 generator=RandomNumberGenerator::gen;
  std::uniform_real_distribution<double> distribution=RandomNumberGenerator::dist;
  bool normalFlag=RandomNumberGenerator::NormalSpareAvailable;
  double normal=RandomNumberGenerator::NormalSpare;
  FL::cFieldLineSegment* first=nullptr;
  int originalOwner=-1;
  Fixture() {
    steps.assign(PIC::ParticleWeightTimeStep::GlobalTimeStep,
        PIC::ParticleWeightTimeStep::GlobalTimeStep+PIC::nTotalSpecies);
    weights.assign(PIC::ParticleWeightTimeStep::GlobalParticleWeight,
        PIC::ParticleWeightTimeStep::GlobalParticleWeight+PIC::nTotalSpecies);
    for(int line=0;line<FL::nFieldLine;++line)
      for(int vertex=0;vertex<=FL::FieldLinesAll[line].GetTotalSegmentNumber();++vertex) {
        Point point;point.vertex=FL::FieldLinesAll[line].GetVertex(vertex);
        point.vertex->GetX(point.x.data());points.push_back(point);
      }
    first=FL::FieldLinesAll[0].GetFirstSegment();originalOwner=first->Thread;
  }
  ~Fixture() {
    std::set<long int> created;
    for(int line=0;line<FL::nFieldLine;++line)
      for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext()) {
        for(auto particle=segment->FirstParticleIndex;particle!=-1;) {
          if(!created.insert(particle).second)break;
          particle=PB::GetNext(particle);
        }
        for(auto particle=segment->tempFirstParticleIndex.load();particle!=-1;) {
          if(!created.insert(particle).second)break;
          particle=PB::GetNext(particle);
        }
        segment->FirstParticleIndex=-1;segment->tempFirstParticleIndex=-1;
      }
    for(auto particle:created)PB::DeleteParticle(particle);
    if(first)first->Thread=originalOwner;
    for(const auto& point:points) {
      double restored[3]={point.x[0],point.x[1],point.x[2]};point.vertex->SetX(restored);
    }
    std::copy(steps.begin(),steps.end(),PIC::ParticleWeightTimeStep::GlobalTimeStep);
    std::copy(weights.begin(),weights.end(),PIC::ParticleWeightTimeStep::GlobalParticleWeight);
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=initialized;
    RandomNumberGenerator::rndLastSeed=nativeSeed;RandomNumberGenerator::gen=generator;
    RandomNumberGenerator::dist=distribution;RandomNumberGenerator::NormalSpareAvailable=normalFlag;
    RandomNumberGenerator::NormalSpare=normal;
  }
  bool GeometryUnchanged() const {
    for(const auto& point:points) {
      double x[3];point.vertex->GetX(x);
      if(std::memcmp(x,point.x.data(),sizeof(x))!=0)return false;
    }
    return true;
  }
};
struct BirthContext {SourceLifecycle::Event event;std::uint64_t initialized=0;};
void Initialize(long int particle,std::uint64_t macro,void* data) {
  auto& context=*static_cast<BirthContext*>(data);
  A::InitializeParticleTransportStateWithIdentity(particle,
      SourceLifecycle::ParticleIdentity(context.event,macro));
  ++context.initialized;
}
template<class Value> Value Read(long int particle,int offset) {
  Value value;std::memcpy(&value,PB::GetParticleDataPointer(particle)+offset,sizeof(value));return value;
}

Result ActualBirthGeometryAndOwner() {
  Result result;
  if(!FL::FieldLinesAll||FL::nFieldLine<1||!FL::FieldLinesAll[0].IsInitialized()) {
    result.status=Status::Skip;result.message="native Maxwellian diagnostic requires an initialized fixed field line";return result;
  }
  const long int initial=PB::GetAllPartNum();long int totalInitial=0;
  MPI_Allreduce(&initial,&totalInitial,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(totalInitial!=0) {
    result.status=Status::Skip;result.message="run SOURCE02 in a fresh empty native test process";return result;
  }
#if (_SIMULATION_TIME_STEP_MODE_ != _SINGLE_GLOBAL_TIME_STEP_ && _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_) || _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  result.status=Status::Skip;result.message="SOURCE02 isolates global-numerics alias/owner behavior; local-owner rejection is a separate gate";return result;
#endif
  unsigned failures=0;
  // Exercise the actual native scalar recurrence at the old bridge's worst
  // bit pattern. The old mask-only bridge produces the constant sequence 1/2;
  // the odd-residue bridge must produce distinct strict-open draws instead.
  cRndSeedContainer bad;bad.Seed=UINT64_C(1)<<30;
  if(rnd(&bad)!=0.5||rnd(&bad)!=0.5)++failures;
  for(const auto bits:{UINT64_C(0),UINT64_C(1)<<30,UINT64_C(0x7fffffff),UINT64_MAX}) {
    cRndSeedContainer good;good.Seed=LegacyScalarSeed(bits);
    const double first=rnd(&good),second=rnd(&good);
    if(!(first>0.0&&first<1.0&&second>0.0&&second<1.0&&first!=second))++failures;
  }
  std::uint64_t births=0;
  double firstParallel=0.0,firstCoordinate=0.0;
  int forcedOwner=0,amrOwner=-1;
  const auto seed=TestSeed(2);
  {
    Fixture fixture;
    double x[3];fixture.first->GetBegin()->GetX(x);
    auto* node=PIC::Mesh::mesh?PIC::Mesh::mesh->findTreeNode(x,nullptr):nullptr;
    amrOwner=node?node->Thread:-1;
    forcedOwner=PIC::nTotalThreads>1?
        ((amrOwner>=0?amrOwner:fixture.originalOwner)+1)%PIC::nTotalThreads:fixture.originalOwner;
    fixture.first->Thread=forcedOwner;
    // Even a speed arbitrarily close to c remains well inside this segment
    // during the declared diagnostic interval. No clipping of S is performed.
    const double dt=std::min(1.0e-6,fixture.first->GetLength()/(100.0*SpeedOfLight));
    for(int species=0;species<PIC::nTotalSpecies;++species)
      PIC::ParticleWeightTimeStep::GlobalTimeStep[species]=dt;
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=true;
    double direction[3];fixture.first->GetDir(direction);
    double bulk[3]={4.0e5*direction[0],4.0e5*direction[1],4.0e5*direction[2]};
    double normal[3]={-direction[0],-direction[1],-direction[2]};
    const double flux=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(1.0e6,1.0e5,bulk,normal,0);
    const double area=SEP::FieldLine::FluxTubeGeometry::AreaAtVertexM2(fixture.first->GetBegin(),0);
    const double weight=0.5*dt*flux*area;
    const int localReady=std::isfinite(weight)&&weight>0.0&&std::isfinite(dt)&&dt>0.0;
    int ready=0;MPI_Allreduce(&localReady,&ready,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    if(!ready)++failures;
    else {
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[0]=weight;
      // The existing core scalar generator is 31-bit. Its projection is
      // explicit and reported; this diagnostic does not claim purpose-keyed
      // production Maxwellian sampling or decomposition-independent noise.
      T::KeyedRandomStream mapping(seed,UINT64_C(0x534f555243453032),UINT64_C(0x434f5245524e47),0);
      const auto bits=mapping.SaveState().integerState;
      RandomNumberGenerator::rndLastSeed=LegacyScalarSeed(bits);
      RandomNumberGenerator::gen.seed(static_cast<std::uint32_t>(bits));
      RandomNumberGenerator::dist=std::uniform_real_distribution<double>(0.0,1.0);
      RandomNumberGenerator::NormalSpareAvailable=false;RandomNumberGenerator::NormalSpare=0.0;
      const auto beforeSeed=RandomNumberGenerator::rndLastSeed;
      const auto beforeGenerator=RandomNumberGenerator::gen;
      const auto beforeNormal=RandomNumberGenerator::NormalSpare;
      const auto preflight=FL::PreflightMaxwellianLineBeginning(0,0);
      using Code=FL::MaxwellianLineBeginningPreflight::Code;
      if(preflight.ownerRank!=forcedOwner||preflight.code!=(PIC::ThisThread==forcedOwner?Code::Ready:Code::NonOwner))++failures;
      const auto invalid=FL::PreflightMaxwellianLineBeginning(-1,0);
      if(invalid.code!=Code::InvalidInput||PB::GetAllPartNum()!=0||
          beforeSeed!=RandomNumberGenerator::rndLastSeed||beforeGenerator!=RandomNumberGenerator::gen||
          std::memcmp(&beforeNormal,&RandomNumberGenerator::NormalSpare,sizeof(beforeNormal))!=0)++failures;
      BirthContext context;context.event.campaign=seed;context.event.index=UINT64_C(0x534f555243453032);
      context.event.channel.kind=SourceLifecycle::Kind::MaxwellianBeginning;
      const int created=FL::InjectMaxwellianLineBeginning(0,1.0e6,1.0e5,bulk,area,0,4,Initialize,&context);
      if(created!=(PIC::ThisThread==forcedOwner?2:0)||context.initialized!=static_cast<std::uint64_t>(created))++failures;
      for(auto particle=fixture.first->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        ++births;
        const double coordinate=PB::GetFieldLineCoord(particle);
        if(PIC::ThisThread!=forcedOwner||PB::GetFieldLineId(particle)!=0||coordinate<0.0||coordinate>=1.0||
            Read<std::uint64_t>(particle,SEP::Offset::TransportSchema)!=A::ParticleTransportSchema||
            Read<std::uint64_t>(particle,SEP::Offset::StableParticleId)==0)++failures;
        // The generic raw template must preserve the allocator's fresh payload
        // contract for extensions that the postFresh RNG hook does not touch.
        if(Read<double>(particle,SEP::Offset::MeanFreePath)!=0.0||
            PB::TestInitFlag(PB::GetParticleDataPointer(particle)))++failures;
        if(firstParallel==0.0){firstParallel=PB::GetVParallel(particle);firstCoordinate=coordinate;}
      }
      if(!fixture.GeometryUnchanged())++failures;
    }
  }
  if(PB::GetAllPartNum()!=initial)++failures;
  unsigned globalFailures=0;std::uint64_t globalBirths=0;
  MPI_Allreduce(&failures,&globalFailures,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&births,&globalBirths,1,MPI_UINT64_T,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  double observedParallel=0.0,observedCoordinate=0.0;
  MPI_Allreduce(&firstParallel,&observedParallel,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&firstCoordinate,&observedCoordinate,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  result.status=globalFailures?Status::Fail:Status::Pass;
  result.message=globalFailures?"native Maxwellian alias/owner diagnostic failed":
      "actual native Maxwellian births preserve physical vertices and follow the field-line segment owner";
  result.hasSeed=true;result.seed=seed;
  result.configuration={"evidence_level=native-generic-PIC-Maxwellian-birth;alias-and-initial-owner-only",
      "run_purpose=controlled-software-diagnostic;source_physics_qualified=false",
      "limitations=lab-frame-legacy-law;short-dt-no-segment-crossing;birth-timing-and-keyed-sampler-unqualified",
      "core_seed_mapping=keyed-state-positive-odd31;bulk_m_per_s=4e5;T_K=1e5;n_m3=1e6;expected_macros=2",
      "forced_segment_owner="+std::to_string(forcedOwner)+";AMR_endpoint_owner="+std::to_string(amrOwner)};
  result.metrics.push_back({"assertion_failures",static_cast<double>(globalFailures),0,"<=","count"});
  result.metrics.push_back({"initialized_births",static_cast<double>(globalBirths),2,"==","count"});
  result.metrics.push_back({"actual_sample_parallel_velocity",observedParallel,0,">","m/s"});
  result.metrics.push_back({"actual_sample_coordinate",observedCoordinate,0,">","segment-coordinate"});
  return result;
}
} // namespace
std::vector<Descriptor> MaxwellianBirthNativeDescriptors() {
  Descriptor descriptor;descriptor.id="SOURCE02";descriptor.name="Native Maxwellian geometry and initial-owner isolation";
  descriptor.group="source-native";descriptor.description=descriptor.name;
  descriptor.initialization=InitializationLevel::FieldLineModel;descriptor.runtime=RuntimeClass::Routine;
  descriptor.supportedBuildModes="native initialized fixed-line model with global timestep/weight";
  descriptor.seedPolicy="requested test seed mapped explicitly into legacy core RNG; original full RNG restored";
  descriptor.seedCapability=SeedCapability::Context;descriptor.repeatSupported=true;
  descriptor.stateIsolation="empty native process; all created records deleted; geometry/numerics/owner/RNG restored";
  descriptor.callback=ActualBirthGeometryAndOwner;return {descriptor};
}
} } // namespace SEP::Testing
