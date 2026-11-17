#include "sep_standard_batsrus_wave_validation.h"
#include "../standard_batsrus_wave_adapter.h"
#include "../transport_common.h"
#include "../polarized_wave_native_adapter.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>

namespace SEP { namespace Testing { namespace {
namespace T=Transport;
namespace A=T::PICAdapter;
namespace W=Waves::StandardBatsrus;
namespace FL=PIC::FieldLine;
struct Checks {
  unsigned failures=0;
  void Require(bool condition){if(!condition)++failures;}
  void Near(double actual,double expected){Require(std::isfinite(actual)&&std::abs(actual-expected)<=8*std::numeric_limits<double>::epsilon()*std::abs(expected));}
};
bool Same(const W::MagneticBranchView& a,const W::MagneticBranchView& b) {
  return a.plusMagneticJPerM3==b.plusMagneticJPerM3&&a.minusMagneticJPerM3==b.minusMagneticJPerM3&&
      a.plusVarianceT2==b.plusVarianceT2&&a.minusVarianceT2==b.minusVarianceT2&&a.totalVarianceT2==b.totalVarianceT2&&
      a.generation==b.generation&&a.checksum==b.checksum&&a.branchDirectionsKnown==b.branchDirectionsKnown&&
      a.circularPolarizationsKnown==b.circularPolarizationsKnown;
}
class Fixture {
 public:
  Checks& checks;
  std::shared_ptr<const Background::BackgroundSnapshot> original=Background::SnapshotStore::Instance().Current();
  FL::cFieldLineVertex* vertices[2];std::string saved[2];
  unsigned generation=0;bool replaced=false;
  const long int particleCount=PIC::ParticleBuffer::GetAllPartNum();
  const std::mt19937 generator=RandomNumberGenerator::gen;
  const std::uniform_real_distribution<double> distribution=RandomNumberGenerator::dist;
  const unsigned long seed=RandomNumberGenerator::rndLastSeed;
  const bool normalAvailable=RandomNumberGenerator::NormalSpareAvailable;
  const double normal=RandomNumberGenerator::NormalSpare;
  explicit Fixture(Checks& value):checks(value) {
    auto* segment=FL::FieldLinesAll[0].GetFirstSegment();vertices[0]=segment->GetBegin();vertices[1]=segment->GetEnd();
    for(unsigned i=0;i<2;++i)saved[i].assign(vertices[i]->GetAssociatedDataBufferPointer(),vertices[i]->AssociatedDataLength());
  }
  ~Fixture() {
    auto& store=Background::SnapshotStore::Instance();
    try {
      if(replaced)store.ResetForTests();
      for(unsigned i=0;i<2;++i)std::memcpy(vertices[i]->GetAssociatedDataBufferPointer(),saved[i].data(),saved[i].size());
      if(replaced&&original) {
        // Reset/temporary publication belong only to this isolated controlled
        // fixture. Restore the original shared_ptr as well as metadata, so no
        // prepared owner's background lifetime identity is changed by a test.
        store.Publish(*original);store.RestoreCheckpoint(original);
      }
      checks.Require(store.Current()==original);
    }catch(const std::exception&){checks.Require(false);}
    checks.Require(PIC::ParticleBuffer::GetAllPartNum()==particleCount&&RandomNumberGenerator::gen==generator&&
        RandomNumberGenerator::dist==distribution&&RandomNumberGenerator::rndLastSeed==seed&&
        RandomNumberGenerator::NormalSpareAvailable==normalAvailable&&std::memcmp(&normal,&RandomNumberGenerator::NormalSpare,sizeof(normal))==0);
  }
  void Publish(const std::array<double,2>& first,const std::array<double,2>& last,int polarity=1) {
    auto& store=Background::SnapshotStore::Instance();
    if(!replaced){store.ResetForTests();replaced=true;}
    else store.AssertProviderMayWrite(Background::Provider::Swmf);
    // Explicit synthetic importer-owned updates occur only while quiescent,
    // before publishing the next immutable generation. They are not produced
    // by a live BATSRUS run and do not qualify a full SWMF plasma provider.
    const std::array<double,2> values[2]={first,last};
    for(unsigned i=0;i<2;++i) {
      double raw[2]={values[i][0],values[i][1]};vertices[i]->SetDatum(FL::DatumAtVertexPlasmaWaves,raw);
      double magnetic[3]={static_cast<double>(polarity)*3e-8,0.0,0.0};vertices[i]->SetDatum(FL::DatumAtVertexMagneticField,magnetic);
    }
    ++generation;const double epoch=generation;
    store.Publish(Background::BackgroundSnapshot(Background::Provider::Swmf,Background::Ownership::ImportedReadOnly,
        epoch,epoch,epoch+1.0,generation,"controlled-standard-batsrus-wave-fixture",W::ContractIdentity(),epoch-1.0,epoch));
  }
  A::ParticleContext Context(double fraction)const {
    A::ParticleContext context;context.segment=FL::FieldLinesAll[0].GetFirstSegment();context.state.fieldLineId=0;
    context.state.coordinate=fraction;context.snapshotGeneration=generation;return context;
  }
  T::Status Read(const A::ParticleContext& context,W::MagneticBranchView* output) {
    std::string before[2];for(unsigned i=0;i<2;++i)before[i].assign(vertices[i]->GetAssociatedDataBufferPointer(),vertices[i]->AssociatedDataLength());
    T::Status status;
    {auto read=Background::SnapshotStore::Instance().BeginParticleRead(generation);status=A::ReadStandardBatsrusWaveView(context,output);}
    for(unsigned i=0;i<2;++i)checks.Require(std::memcmp(before[i].data(),vertices[i]->GetAssociatedDataBufferPointer(),before[i].size())==0);
    return status;
  }
};
bool Supported(Result* result) {
  int local=FL::FieldLinesAll&&FL::nFieldLine>0&&FL::FieldLinesAll[0].GetFirstSegment()&&
      FL::FieldLinesAll[0].GetFirstSegment()->GetBegin()&&FL::FieldLinesAll[0].GetFirstSegment()->GetEnd()&&
      FL::DatumAtVertexPlasmaWaves.offset>=0&&FL::DatumAtVertexPlasmaWaves.length==2&&
      FL::DatumAtVertexMagneticField.offset>=0&&FL::DatumAtVertexMagneticField.length==3&&
      !A::NativeDriverWindowActive()&&!Turbulence::Polarized::PICAdapter::HasNativePolarizedState()&&
      PIC::ParticleBuffer::GetAllPartNum()==0;
  int all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(!all){result->status=Status::Skip;result->message="standard BATSRUS native fixture needs empty initialized native lines/two-wave vertex storage and no live source/four-wave transaction";return false;}
  return true;
}
Result Finish(const Checks& checks,const char* message) {
  unsigned failures=0;MPI_Allreduce(&checks.failures,&failures,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  Result result;result.status=failures?Status::Fail:Status::Pass;result.message=failures?"native standard BATSRUS wave accessor assertion failed":message;
  result.configuration={"evidence_level=actual-native-vertex-arrays-and-production-read-only-accessor",
      "synthetic-importer=true;standard-total-J/m3-to-magnetic-half;Wave1=along-B;Wave2=opposite-B",
      "immutable-raw-array-and-RNG-and-PB-count=true;original-snapshot-pointer-restored=true",
      "full-SWMF-coupling-gradient-composition-and-production-physics-qualified=false"};
  result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0,"<=","count"});return result;
}
Result NativeView() {
  Result missing;if(!Supported(&missing))return missing;Checks checks;
  {
    Fixture fixture(checks);const std::array<double,2> first{{8e-9,2e-9}},last{{4e-9,10e-9}};
    W::MagneticBranchView forward,flipped,reversed;
    fixture.Publish(first,last,1);checks.Require(fixture.Read(fixture.Context(0.25),&forward).ok());
    checks.Near(forward.plusVarianceT2,1.256637062e-6*7e-9);checks.Near(forward.minusVarianceT2,1.256637062e-6*4e-9);
    fixture.Publish(first,last,-1);checks.Require(fixture.Read(fixture.Context(0.25),&flipped).ok());
    checks.Require(forward.plusVarianceT2==flipped.plusVarianceT2&&forward.minusVarianceT2==flipped.minusVarianceT2);
    fixture.Publish(last,first,-1);checks.Require(fixture.Read(fixture.Context(0.75),&reversed).ok());
    checks.Require(forward.plusVarianceT2==reversed.plusVarianceT2&&forward.minusVarianceT2==reversed.minusVarianceT2&&
        reversed.branchDirectionsKnown&&!reversed.circularPolarizationsKnown&&forward.checksum!=flipped.checksum);
  }
  return Finish(checks,"actual read-only vertex accessor preserves along/opposite-B branches under magnetic-polarity and endpoint/fraction reversal without raw mutation");
}
Result NativeRejection() {
  Result missing;if(!Supported(&missing))return missing;Checks checks;
  {
    Fixture fixture(checks);W::MagneticBranchView sentinel;sentinel.checksum=19;sentinel.plusVarianceT2=7;
    auto context=fixture.Context(0.25);auto output=sentinel;
    checks.Require(!A::ReadStandardBatsrusWaveView(context,&output).ok()&&Same(output,sentinel));
    // Exercise a foreign provider in every configured build, even if startup
    // itself was SWMF. This is explicit test-only metadata replacement; the
    // fixture restores the exact original pointer and all backing bytes.
    auto& store=Background::SnapshotStore::Instance();store.ResetForTests();fixture.replaced=true;
    store.Publish(Background::BackgroundSnapshot(Background::Provider::Analytic,Background::Ownership::ModelOwned,
        1.0,1.0,2.0,1,"controlled-foreign-wave-source","foreign-provider rejection fixture"));
    context.snapshotGeneration=1;
    {auto read=store.BeginParticleRead(1.0);
     checks.Require(!A::ReadStandardBatsrusWaveView(context,&output).ok()&&Same(output,sentinel));}
    store.ResetForTests();
    const std::array<double,2> valid{{8e-9,2e-9}};
    fixture.Publish(valid,valid);context=fixture.Context(0.25);++context.snapshotGeneration;
    checks.Require(!fixture.Read(context,&output).ok()&&Same(output,sentinel));
    context=fixture.Context(0.25);context.state.coordinate=-1;
    checks.Require(!fixture.Read(context,&output).ok()&&Same(output,sentinel));
    for(double invalid:{-1.0,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::denorm_min(),1e-320}) {
      std::array<double,2> bad=valid;bad[0]=invalid;fixture.Publish(valid,bad);
      checks.Require(!fixture.Read(fixture.Context(0.0),&output).ok()&&Same(output,sentinel));
    }
  }
  return Finish(checks,"absent read phase, foreign provider, stale generation, invalid coordinate and invalid unweighted raw endpoints reject without imported/RNG/particle mutation");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.description=name;
  descriptor.group="standard-batsrus-waves-native";descriptor.initialization=InitializationLevel::FieldLineModel;
  descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
  descriptor.seedPolicy="deterministic explicit imported values; legacy/keyed RNG unchanged";
  descriptor.stateIsolation="empty native model; two associated vertex buffers and exact background pointer restored; no particle allocation";
  descriptor.supportedBuildModes="native model with allocated wave datums; controlled synthetic import only";descriptor.callback=callback;return descriptor;
}
} // namespace
std::vector<Descriptor> StandardBatsrusWaveNativeDescriptors() {
  return {Make("SBWN01","Native standard BATSRUS read-only magnetic wave view",NativeView),
      Make("SBWN02","Native imported wave ownership and atomic rejection",NativeRejection)};
}
} }
