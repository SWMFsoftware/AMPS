#include "sep_population_native_validation.h"
#include "../transport_common.h"
#include "sep_population_control.h"
#include "sep_population_production_adapter.h"
#include "sep_source_lifecycle.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SEP { namespace Testing { namespace {
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace Split=PIC::ParticleSplitting::FledLine;
namespace T=SEP::Transport;
namespace Adapter=T::PICAdapter;

template<class Value> Value Read(long int particle,int offset) {
  Value value;
  std::memcpy(&value,PB::GetParticleDataPointer(particle)+offset,sizeof(value));
  return value;
}
template<class Value> void Write(long int particle,int offset,const Value& value) {
  std::memcpy(PB::GetParticleDataPointer(particle)+offset,&value,sizeof(value));
}
bool Same(const T::KeyedRandomStreamState& a,const T::KeyedRandomStreamState& b) {
  return a.schema==b.schema&&a.keyIdentity==b.keyIdentity&&
      a.integerState==b.integerState&&a.hasSpareNormal==b.hasSpareNormal&&
      a.spareNormal==b.spareNormal;
}
struct Checks {
  unsigned failures=0;
  std::vector<std::string> messages;
  void Require(bool value,const char* message) {
    if(value)return;
    ++failures;
    if(messages.size()<8)messages.push_back(message);
  }
};

// Detach original lists without changing any original particle byte. All test
// allocations belong to this scope; even an exception deletes only these records
// and restores original heads, temp heads, and optional callbacks. The next
// production population pass refreshes its private identity reservation cache.
class Fixture {
 public:
  struct SavedList {FL::cFieldLineSegment* segment;long int head,temp;};
  std::vector<SavedList> saved;
  FL::cFieldLineSegment* segment=nullptr;
  int line=-1,segmentIndex=-1;
  std::vector<long int> unlisted;
  Split::fParticleIdentity identity=Split::ParticleIdentity;
  Split::fPrepareExactSplit prepare=Split::PrepareExactSplit;
  Split::fPostExactSplit post=Split::PostExactSplit;
  Split::fMergeController merger=Split::MergeController;
  explicit Fixture(int minimumSegmentIndex=0) {
    if(!FL::FieldLinesAll)return;
    for(int l=0;l<FL::nFieldLine;++l) {
      if(!FL::FieldLinesAll[l].IsInitialized())continue;
      int index=0;
      for(auto* current=FL::FieldLinesAll[l].GetFirstSegment();current;
          current=current->GetNext(),++index) {
        saved.push_back({current,current->FirstParticleIndex,
                         current->tempFirstParticleIndex.load()});
        if(PIC::ThisThread==0&&current->Thread==0&&!segment&&index>=minimumSegmentIndex) {
          segment=current;line=l;segmentIndex=index;
        }
        current->FirstParticleIndex=-1;
        current->tempFirstParticleIndex=-1;
      }
    }
  }
  ~Fixture() {
    std::set<long int> owned(unlisted.begin(),unlisted.end());
    for(const auto& entry:saved) {
      long int particle=entry.segment->FirstParticleIndex;
      while(particle!=-1) {
        if(!owned.insert(particle).second)break;
        particle=PB::GetNext(particle);
      }
      entry.segment->FirstParticleIndex=-1;
      entry.segment->tempFirstParticleIndex=-1;
    }
    for(long int particle:owned)PB::DeleteParticle(particle);
    for(const auto& entry:saved) {
      entry.segment->FirstParticleIndex=entry.head;
      entry.segment->tempFirstParticleIndex=entry.temp;
    }
    Split::ParticleIdentity=identity;
    Split::PrepareExactSplit=prepare;
    Split::PostExactSplit=post;
    Split::MergeController=merger;
  }
  std::vector<long int> Particles() const {
    std::vector<long int> particles;
    if(segment)for(long int p=segment->FirstParticleIndex;p!=-1;p=PB::GetNext(p))
      particles.push_back(p);
    return particles;
  }
  long int Add(unsigned trial,double weight,bool linked=true) {
    const long int particle=linked?PB::GetNewParticle(segment->FirstParticleIndex):
                                   PB::GetNewParticle();
    if(!linked)unlisted.push_back(particle);
    PB::SetI(0,particle);
    PB::SetFieldLineId(line,particle);
    PB::SetFieldLineCoord(segmentIndex+0.1+0.07*trial,particle);
    PB::SetIndividualStatWeightCorrection(weight,particle);
    const double parallel=1.0e7+1.0e6*trial,normal=2.0e7+1.0e6*trial;
    PB::SetVParallel(parallel,particle);PB::SetVNormal(normal,particle);
    double position[3];
    FL::FieldLinesAll[line].GetCartesian(position,PB::GetFieldLineCoord(particle));
    PB::SetX(position,particle);
    const double mass=PIC::MolecularData::GetMass(0);
    const double speed=std::hypot(parallel,normal);
    const auto momentum=T::MomentumFromSpeed(speed,mass,SpeedOfLight);
    if(!momentum.status.ok())throw std::runtime_error(momentum.status.message);
    const double p[3]={momentum.value*parallel/speed,
                       momentum.value*normal/speed,0.0};
    // Rank zero is only the fixture owner. Physical line/segment/trial keys,
    // rather than a worker/rank, define these source identities.
    const std::uint64_t source=UINT64_C(0x504f504e00000000)|
        static_cast<std::uint64_t>(line);
    const std::uint64_t sequence=(static_cast<std::uint64_t>(segmentIndex)<<32)|trial;
    Adapter::InitializeParticleTransportState(particle,source,sequence,p);
    return particle;
  }
};

bool Ready(const Fixture& fixture,Checks& checks) {
  const int local=(PIC::ThisThread!=0||fixture.segment)?1:0;
  int ready=0;
  MPI_Allreduce(&local,&ready,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  checks.Require(ready!=0,"rank-zero-owned native field-line segment is unavailable");
  return ready!=0;
}
void GiveHistory(long int particle) {
  const auto id=Read<std::uint64_t>(particle,SEP::Offset::StableParticleId);
  T::OperatorRandomState states[3];
  for(int operation=0;operation<3;++operation) {
    T::KeyedRandomStream random(Adapter::CampaignRandomSeed,id,7+operation,0);
    random.UniformOpen01();random.Normal01();
    states[operation].stream=random.SaveState();
    states[operation].acceptedIntervals=17+operation;
  }
  std::memcpy(PB::GetParticleDataPointer(particle)+SEP::Offset::OperatorRandomStates,
              states,sizeof(states));
  Write(particle,SEP::Offset::MfpOpticalDepth,0.37);
  Write(particle,SEP::Offset::MfpEventIndex,UINT64_C(23));
  Write(particle,SEP::Offset::MeanFreePath,42.0);
}
struct Observation {
  long int particle;
  double coordinate,parallel,normal,weight,hazard;
  std::uint64_t id,event,counter;
  T::OperatorRandomState states[3];
};
Observation Observe(long int particle) {
  Observation o;
  o.particle=particle;o.coordinate=PB::GetFieldLineCoord(particle);
  o.parallel=PB::GetVParallel(particle);o.normal=PB::GetVNormal(particle);
  o.weight=PB::GetIndividualStatWeightCorrection(particle);
  o.id=Read<std::uint64_t>(particle,SEP::Offset::StableParticleId);
  o.event=Read<std::uint64_t>(particle,SEP::Offset::MfpEventIndex);
  o.counter=Read<std::uint64_t>(particle,SEP::Offset::SplitEventCounter);
  o.hazard=Read<double>(particle,SEP::Offset::MfpOpticalDepth);
  std::memcpy(o.states,PB::GetParticleDataPointer(particle)+
      SEP::Offset::OperatorRandomStates,sizeof(o.states));
  return o;
}
long double Energy(const std::vector<long int>& particles) {
  long double total=0.0;
  const long double c=SpeedOfLight,mass=PIC::MolecularData::GetMass(0);
  for(long int p:particles) {
    const long double parallel=PB::GetVParallel(p),normal=PB::GetVNormal(p);
    const long double gamma=1.0L/std::sqrt(1.0L-
        (parallel*parallel+normal*normal)/(c*c));
    total+=PB::GetIndividualStatWeightCorrection(p)*mass*c*c*(gamma-1.0L);
  }
  return total;
}
Result Finish(const Checks& checks,const char* success) {
  Result result;
  result.status=checks.failures?Status::Fail:Status::Pass;
  std::ostringstream message;
  if(!checks.failures)message<<success;
  else for(std::size_t i=0;i<checks.messages.size();++i) {
    if(i)message<<"; ";
    message<<checks.messages[i];
  }
  result.message=message.str();result.hasSeed=true;result.seed=Adapter::CampaignRandomSeed;
  result.configuration.push_back("evidence_level=native-PIC-record/controller;MPI-identity-preflight");
  result.configuration.push_back("fixture_owner_rank=0;physical_keys=line,segment,trial");
  result.metrics.push_back({"assertion_failures",static_cast<double>(checks.failures),0,"<=","count"});
  return result;
}

Result SplitLifecycle() {
  Checks checks;
  const long int initialCount=PB::GetAllPartNum();
  double energyError=0.0;
  {
    Fixture fixture;
    if(Ready(fixture,checks)) {
      std::vector<Observation> original;
      long double beforeEnergy=0.0;
      // Setup failure must be reduced before the next collective preflight, so
      // one owner cannot throw while another rank enters MPI_Allgatherv.
      int localSetup=1;
      if(PIC::ThisThread==0)try {
        for(unsigned i=0;i<10;++i) {
          const long int particle=fixture.Add(i,i==0?8.0:1.0);
          GiveHistory(particle);original.push_back(Observe(particle));
        }
        beforeEnergy=Energy(fixture.Particles());
      } catch(const std::exception& exception) {
        checks.Require(false,exception.what());localSetup=0;
      }
      int setup=0;
      MPI_Allreduce(&localSetup,&setup,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
      if(setup) {
        Adapter::RegisterPopulationControlHooks();
        Adapter::BeginPopulationControl();
        if(PIC::ThisThread==0) {
          double vpRange=0,vnRange=0;
          fixture.segment->WeightedParticleSplitting(0,fixture.segment->FirstParticleIndex,
              2,2,2,1.0,vpRange,vnRange,10,20,2);
          const auto particles=fixture.Particles();
          checks.Require(particles.size()==12,"actual first core split did not produce twelve records");
          std::set<std::uint64_t> ids;
          for(long int particle:particles) {
            const auto current=Observe(particle);
            checks.Require(ids.insert(current.id).second,"split produced a duplicate stable identity");
            const auto parent=std::find_if(original.begin(),original.end(),
                [&current](const Observation& o){return o.coordinate==current.coordinate;});
            checks.Require(parent!=original.end(),"split changed a physical coordinate");
            if(parent==original.end())continue;
            checks.Require(current.parallel==parent->parallel&&current.normal==parent->normal,
                           "exact split changed intrinsic velocity");
            checks.Require(Read<double>(particle,SEP::Offset::MeanFreePath)==42.0,
                           "split lost an unrelated requested extension value");
            if(particle==parent->particle) {
              checks.Require(current.id==parent->id&&current.event==parent->event&&
                             current.hazard==parent->hazard,"true parent lost its identity/event/hazard history");
              checks.Require(current.counter<=1&&current.weight==
                  (current.counter?0.5*parent->weight:parent->weight),
                  "true parent weight or split serial changed inconsistently");
              for(int op=0;op<3;++op)checks.Require(Same(current.states[op].stream,parent->states[op].stream)&&
                  current.states[op].acceptedIntervals==parent->states[op].acceptedIntervals,
                  "true parent lost generator/cache/interval history");
            } else {
              checks.Require(current.id==PopulationControl::DeriveChildId(parent->id,1,1,
                  Adapter::ParticleTransportSchema),"child identity differs from persisted lineage key");
              checks.Require(current.weight==0.5*parent->weight&&current.counter==0&&
                  current.event==0&&std::isnan(current.hazard),"child weight or fresh future history is invalid");
              for(int op=0;op<3;++op) {
                T::KeyedRandomStream expected(Adapter::CampaignRandomSeed,current.id,7+op,0);
                checks.Require(Same(current.states[op].stream,expected.SaveState())&&
                    current.states[op].acceptedIntervals==0,"child keyed stream/cache was not initialized independently");
                checks.Require(current.states[op].stream.keyIdentity!=parent->states[op].stream.keyIdentity,
                               "child retained its parent's stochastic key");
              }
            }
          }
        }
        // The second controller sees twelve initial representatives. Both the
        // original heavy parent and its first child have weight four, so both
        // are selected once; this tests parent serial persistence across calls.
        Adapter::BeginPopulationControl();
        if(PIC::ThisThread==0) {
          double vpRange=0,vnRange=0;
          fixture.segment->WeightedParticleSplitting(0,fixture.segment->FirstParticleIndex,
              2,2,2,1.0,vpRange,vnRange,12,20,2);
          const auto particles=fixture.Particles();
          checks.Require(particles.size()==14,"actual second core split did not produce fourteen records");
          const auto heavy=Observe(original.front().particle);
          checks.Require(heavy.id==original.front().id&&heavy.counter==2&&heavy.weight==2.0,
                         "retained parent split serial did not persist across native calls");
          for(int op=0;op<3;++op)checks.Require(Same(heavy.states[op].stream,
              original.front().states[op].stream),"second split changed retained parent RNG history");
          double totalWeight=0.0;
          for(long int particle:particles)totalWeight+=PB::GetIndividualStatWeightCorrection(particle);
          checks.Require(totalWeight==17.0,"exact splitting changed represented particle number");
          energyError=static_cast<double>(std::fabs((Energy(particles)-beforeEnergy)/beforeEnergy));
          checks.Require(energyError<=64*std::numeric_limits<double>::epsilon(),
                         "actual PIC split changed relativistic kinetic energy");
        }
      } else checks.Require(false,"native fixture setup failed on another rank");
    }
  }
  checks.Require(PB::GetAllPartNum()==initialCount,"native split fixture leaked PIC allocations");
  Result result=Finish(checks,"actual PIC exact split retains parent history and initializes independent child lineage");
  result.metrics.push_back({"relativistic_energy_relative_error",energyError,
      64*std::numeric_limits<double>::epsilon(),"<=","dimensionless"});
  result.configuration.push_back("initial_records=10;core_pass_targets=12,14");
  return result;
}

Result ByteCloneHistory() {
  Checks checks;
  const long int initialCount=PB::GetAllPartNum();
  {
    Fixture fixture;
    if(Ready(fixture,checks)&&PIC::ThisThread==0) {
      const long int parent=fixture.Add(0,3.0,false);
      GiveHistory(parent);
      Write(parent,SEP::Offset::SplitEventCounter,UINT64_C(17));
      const long int copy=fixture.Add(1,1.0,false);
      const long int next=PB::GetNext(copy),previous=PB::GetPrev(copy);
      std::vector<PB::byte> expected(PB::ParticleDataLength);
      std::memcpy(expected.data(),PB::GetParticleDataPointer(parent),expected.size());
      PB::SetNext(next,expected.data());PB::SetPrev(previous,expected.data());
      PB::CloneParticle(copy,parent);
      checks.Require(std::memcmp(expected.data(),PB::GetParticleDataPointer(copy),expected.size())==0,
                     "ordinary PIC CloneParticle changed requested history bytes or destination links");
      const auto source=Observe(parent),cloned=Observe(copy);
      checks.Require(source.id==cloned.id&&cloned.counter==17&&source.event==cloned.event&&
                     source.hazard==cloned.hazard,"ordinary clone rekeyed a migration/restart history");
      for(int op=0;op<3;++op) {
        T::KeyedRandomStream a(Adapter::CampaignRandomSeed,source.id,7+op,0);
        T::KeyedRandomStream b(Adapter::CampaignRandomSeed,cloned.id,7+op,0);
        checks.Require(a.RestoreState(source.states[op].stream).ok()&&
            b.RestoreState(cloned.states[op].stream).ok(),"copied native generator state cannot resume");
        for(unsigned draw=0;draw<12;++draw)
          checks.Require(a.Normal01()==b.Normal01(),"byte clone lost the cached-normal continuation");
      }
    }
  }
  checks.Require(PB::GetAllPartNum()==initialCount,"native byte clone fixture leaked PIC allocations");
  return Finish(checks,"ordinary PIC clone preserves full lineage, cached RNG and hazard bytes");
}
class SourceScope {
 public:
  SourceLifecycle::State saved=SourceLifecycle::CaptureNativeState();
  Checks& checks;
  explicit SourceScope(Checks& result):checks(result) {}
  ~SourceScope() {
    if(!saved.fixedLineLayout.empty())
      checks.Require(SourceLifecycle::RestoreNativeState(saved).ok(),
                     "native merger fixture failed to restore source scheduler history");
  }
};
Result ConservativeNativeMerge() {
  Checks checks;
  const long int initialCount=PB::GetAllPartNum();
  double energyError=0;
  {
    SourceScope source(checks);
    Fixture fixture(1); // test local fractions on a nonzero global segment
    if(Ready(fixture,checks)) {
      int setup=source.saved.fixedLineLayout.empty()?0:1;
      std::map<std::uint64_t,std::vector<PB::byte>> original;
      long double beforeEnergy=0;
      if(PIC::ThisThread==0&&setup)try {
        for(unsigned i=0;i<10;++i) {
          const long int particle=fixture.Add(i,1.0+.2*i);
          const double speed=(.1+.07*i)*SpeedOfLight,mu=.9*std::cos(.5+1.3*i);
          PB::SetVParallel(speed*mu,particle);
          PB::SetVNormal(speed*std::sqrt(1-mu*mu),particle);
          GiveHistory(particle);
          const auto id=Read<std::uint64_t>(particle,SEP::Offset::StableParticleId);
          auto& bytes=original[id];bytes.resize(PB::ParticleDataLength);
          std::memcpy(bytes.data(),PB::GetParticleDataPointer(particle),bytes.size());
        }
        beforeEnergy=Energy(fixture.Particles());
      } catch(const std::exception& exception) {
        checks.Require(false,exception.what());setup=0;
      }
      int allSetup=0;
      MPI_Allreduce(&setup,&allSetup,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
      if(allSetup) {
        PopulationControl::RegisterConservativeMerger();
        PopulationControl::BeginConservativePopulationPass();
        if(PIC::ThisThread==0) {
          double vp=0,vn=0;
          fixture.segment->WeightedParticleMerging(0,fixture.segment->FirstParticleIndex,
              1,1,1,1.0,vp,vn,6,10,2);
          const auto particles=fixture.Particles();
          checks.Require(particles.size()<=8&&particles.size()>=6,
                         "native conservative merger did not reduce its feasible existing bin");
          for(long int particle:particles) {
            const auto id=Read<std::uint64_t>(particle,SEP::Offset::StableParticleId);
            const auto old=original.find(id);
            checks.Require(old!=original.end(),"native merger invented a phase-point identity");
            if(old==original.end())continue;
            auto expected=old->second;
            PB::SetIndividualStatWeightCorrection(PB::GetIndividualStatWeightCorrection(particle),expected.data());
            PB::SetNext(PB::GetNext(particle),expected.data());
            PB::SetPrev(PB::GetPrev(particle),expected.data());
            checks.Require(std::memcmp(expected.data(),PB::GetParticleDataPointer(particle),expected.size())==0,
                "merger changed a surviving phase point, ID, RNG cache, hazard or unrelated byte");
          }
          energyError=static_cast<double>(std::fabs((Energy(particles)-beforeEnergy)/beforeEnergy));
          checks.Require(energyError<=1e-12,"native conservative merger lost relativistic kinetic energy");
          const auto& reports=PopulationControl::ConservativeMergeReports();
          checks.Require(!reports.empty()&&reports.back().status==PopulationControl::ConservativeMergeStatus::Reduced&&
              reports.back().maximumRelativeResidual<=1e-12,
              "native conservative merger lacks its six-moment transaction report");
        }
        // A separate pass with deliberately sparse configured bins must retain
        // all weights/records and report no reduction. It must never execute the
        // previous 2->1 velocity-average merger as a fallback.
        PopulationControl::BeginConservativePopulationPass();
        if(PIC::ThisThread==0) {
          std::map<long int,std::vector<PB::byte>> before;
          const auto particles=fixture.Particles();
          for(long int particle:particles) {
            auto& bytes=before[particle];bytes.resize(PB::ParticleDataLength);
            std::memcpy(bytes.data(),PB::GetParticleDataPointer(particle),bytes.size());
          }
          double vp=0,vn=0;
          const int count=static_cast<int>(particles.size());
          fixture.segment->WeightedParticleMerging(0,fixture.segment->FirstParticleIndex,
              100,100,100,1.0,vp,vn,2,count,2);
          checks.Require(fixture.Particles()==particles,"sparse merger changed its original native list");
          for(const auto& entry:before)checks.Require(std::memcmp(entry.second.data(),
              PB::GetParticleDataPointer(entry.first),entry.second.size())==0,
              "sparse no-reduction changed native record bytes");
          const auto& reports=PopulationControl::ConservativeMergeReports();
          checks.Require(!reports.empty()&&reports.back().status==PopulationControl::ConservativeMergeStatus::SparseNoReduction&&
              reports.back().before==reports.back().after&&reports.back().sparseBins>0,
              "sparse no-reduction was not explicit or fell back to legacy merging");
        }
      } else checks.Require(false,"native conservative fixture lacks initialized source scheduler history");
    }
  }
  checks.Require(PB::GetAllPartNum()==initialCount,"native conservative merger fixture leaked PIC allocations");
  Result result=Finish(checks,"actual native merger reweights existing histories conservatively and keeps sparse bins unchanged without fallback");
  result.metrics.push_back({"native_relativistic_energy_relative_error",energyError,1e-12,"<=","dimensionless"});
  return result;
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;
  descriptor.id=id;descriptor.name=name;descriptor.group="population-native";
  descriptor.description=name;descriptor.callback=callback;
  descriptor.initialization=InitializationLevel::FieldLineModel;
  descriptor.runtime=RuntimeClass::Routine;
  descriptor.supportedBuildModes="configured native srcSEP PIC field-line model";
  descriptor.seedCapability=SeedCapability::Context;descriptor.repeatSupported=true;
  descriptor.seedPolicy="native TestSeed scope; active campaign seed; explicit physical fixture source/trial identity";
  descriptor.stateIsolation="native allocation/list/hook scopes restored; no production timestep, MPI transfer or checkpoint claim";
  return descriptor;
}
} // namespace
std::vector<Descriptor> PopulationNativeDescriptors() {
  return {Make("POPN01","Native exact-split stochastic lifecycle",SplitLifecycle),
          Make("POPN02","Native byte-clone stochastic continuity",ByteCloneHistory),
          Make("POPN03","Native conservative existing-point merger",ConservativeNativeMerge)};
}
}} // namespace SEP::Testing
