#include "sep_birth_native_validation.h"
#include "../transport_common.h"
#include "../coefficient_providers.h"
#include "../physics_checkpoint.h"
#include "sep_background_runtime.h"
#include "sep_source_lifecycle.h"
#include "sep_analytic_radial_background.h"
#include "sep_initialization.h"
#include "../adapters/swcme1d_adapter.h"
#include "../reduced_shock_geometry.h"
#include "../reduced_shock_particles.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <unistd.h>

namespace SEP { namespace Testing { namespace {
namespace T=Transport;
namespace A=T::PICAdapter;
namespace B=T::Birth;
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace W=AlfvenTurbulence_Kolmogorov;

struct Checks {
  unsigned failures=0;std::string message;
  double maximumCartesianViewResidualM=0.0;
  void Require(bool condition,const std::string& reason) {
    if(condition)return;
    ++failures;if(message.size()<1800){if(!message.empty())message+="; ";message+=reason;}
  }
  bool Status(const T::Status& status){Require(status.ok(),status.message);return status.ok();}
  bool Agree(bool local) {
    int value=local?1:0,all=0;MPI_Allreduce(&value,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    if(!all)Require(false,"birth fixture prerequisite failed on an MPI rank");
    return all!=0;
  }
};

// This diagnostic owns a REAL radial material map, while deliberately retaining
// the application's documented departure-snapshot Lie composition.  The native
// particle operator sees one immutable departure chart R0. Only after workers
// and readers join is M_D published through the common field-line transaction.
// Thus P_R0(h) is tested independently of M_D: no dense-time particle equations,
// additional U displacement, or fractional physical sampling are introduced.
class RadialMaterial {
 public:
  using Point=FL::Geometry::Point;
  Checks& check;
  Background::AnalyticRadial::Wind wind=Background::AnalyticRadial::DeclaredNativeWind();
  double epoch=0,duration=0.01;
  std::vector<Point> departure;
  explicit RadialMaterial(Checks& c):check(c){}
  Point Map(const Point& a,double elapsed) const {
    Point x;double r=0;
    for(unsigned d=0;d<3;++d)r=std::hypot(r,a[d]-wind.originM[d]);
    if(!std::isfinite(r)||r<=0||r+wind.speedMPerS*elapsed<=0)
      throw std::runtime_error("radial material map is outside its declared positive-radius domain");
    for(unsigned d=0;d<3;++d)x[d]=wind.originM[d]+(a[d]-wind.originM[d])*(1+wind.speedMPerS*elapsed/r);
    return x;
  }
  Point InitialB(const Point& a) const {
    if(!Initialization::HasActive())throw std::runtime_error("birth material fixture requires explicit schema6 startup");
    const auto& c=Initialization::Active();Point b;
    ParkerSpiral::GetB(b.data(),a.data(),wind.originM.data(),c.innerRadiusM,
        c.solarWindSpeedMPerS,c.solarRotationRateRadPerS,
        SW1DAdapter::GetConfigurationSummary().parker_radial_field_at_one_au_t);
    return b;
  }
  // The declared initial field is frozen into the radial material motion:
  // F=ee^T+(r/r0)(I-ee^T), J=(r/r0)^2, B=F B0/J, n=n0/J.
  // This is an ideal induction/continuity push, not the unrelated operation of
  // evaluating an unchanged Eulerian Parker spiral at a moved vertex. Thermal
  // and wave evolution are not exercised by these timing/operator fixtures.
  Point PushedB(const Point& a,const Point& b,double elapsed) const {
    Point e;double r=0;for(unsigned d=0;d<3;++d)r=std::hypot(r,a[d]-wind.originM[d]);
    const double scale=1+wind.speedMPerS*elapsed/r;
    double parallel=0;for(unsigned d=0;d<3;++d){e[d]=(a[d]-wind.originM[d])/r;parallel+=e[d]*b[d];}
    Point result;for(unsigned d=0;d<3;++d)result[d]=(parallel*e[d]+scale*(b[d]-parallel*e[d]))/(scale*scale);
    return result;
  }
  bool Install(double start) {
    epoch=start;PIC::SimulationTime::SetInitialValue(start);
    auto& line=FL::FieldLinesAll[0];
    for(int i=0;i<=line.GetTotalSegmentNumber();++i) {
      const auto* raw=line.GetVertex(i)->GetX();const Point a{{raw[0],raw[1],raw[2]}};
      departure.push_back(a);Point b;line.GetVertex(i)->GetMagneticField(b.data());
      const auto expected=InitialB(a);
      double scale=0;for(double component:expected)scale=std::hypot(scale,component);
      for(unsigned d=0;d<3;++d)check.Require(std::fabs(b[d]-expected[d])<=2e-12*scale,
          "birth initial native B differs from explicitly configured Parker field");
      Background::AnalyticRadial::Sample sample;
      check.Status(Background::AnalyticRadial::Evaluate(wind,a.data(),&sample));
      check.Status(Background::AnalyticRadial::ValidateVelocity(wind,a.data(),line.GetVertex(i)->GetPlasmaVelocity(),2e-12));
      // Independently differentiated flow map verifies the full Cartesian
      // tensor, including transverse strain which cannot be inferred from a
      // density history or merely from the scalar divergence.
      for(unsigned row=0;row<3;++row)for(unsigned column=0;column<3;++column) {
        const double eRow=(a[row]-wind.originM[row])/sample.radiusM;
        const double eColumn=(a[column]-wind.originM[column])/sample.radiusM;
        const double expectedGradient=wind.speedMPerS/sample.radiusM*((row==column?1.:0.)-eRow*eColumn);
        check.Require(std::fabs(sample.gradientPerS[3*row+column]-expectedGradient)<=
            16*std::numeric_limits<double>::epsilon()*wind.speedMPerS/sample.radiusM,
            "declared radial tensor differs from independent Cartesian flow-map derivative");
      }
    }
    const auto registered=line.RegisterGeometryOwner("native-birth-radial-operator-stage-v1",epoch,0);
    check.Require(registered.ok(),registered.message);if(!registered.ok())return false;
    // This is a NEW explicitly declared diagnostic R0, before any particle or
    // source event. Startup previously initialized a fixed-geometry namespace;
    // retaining that namespace after real owner registration would be false.
    // Do not reset a continuing scheduler: even an empty particle population
    // does not imply unused source history, so require its counter table empty
    // and no active window, collectively, before constructing the fresh owner.
    const auto originalSource=SourceLifecycle::CaptureNativeState();
    const bool fresh=PB::GetAllPartNum()==0&&!A::NativeDriverWindowActive()&&originalSource.counters.empty();
    if(!check.Agree(fresh))return false;
    SourceLifecycle::Scheduler preparedSource;bool sourcePrepared=false;
    try {
    SourceLifecycle::Scheduler candidate;
    const auto& identity=line.GetGeometryOwnerIdentity();
    // This diagnostic deliberately couples to field_line.cpp's PRIVATE
    // one-line material-namespace format. PrepareNativeRestore below recomputes
    // the canonical namespace independently from the live registered geometry.
    // A future format change must visibly reject this fixture; it must NEVER
    // guess another key, copy fixed-grid history or bypass that validation.
    const std::string layout="fixed-native-lines-v1:1;material:0:"+std::to_string(identity.size())+":"+identity;
    if(!check.Status(candidate.Initialize(A::CampaignRandomSeed,layout)))throw std::runtime_error("cannot initialize fresh diagnostic source owner");
    const auto freshSource=candidate.Capture();
    for(bool wrongCampaign:{false,true}) {
      auto incompatible=freshSource;
      if(wrongCampaign)incompatible.campaign^=UINT64_C(1);
      else incompatible.fixedLineLayout+=";wrong-namespace";
      SourceLifecycle::Scheduler rejected;
      const auto invalid=SourceLifecycle::PrepareNativeRestore(incompatible,&rejected);
      const auto unchanged=SourceLifecycle::CaptureNativeState();
      check.Require(!invalid.ok()&&unchanged.campaign==originalSource.campaign&&
          unchanged.fixedLineLayout==originalSource.fixedLineLayout&&unchanged.counters.empty(),
          "wrong fresh material source namespace/campaign was accepted or mutated the original scheduler");
    }
    const auto sourceReady=SourceLifecycle::PrepareNativeRestore(freshSource,&preparedSource);
    check.Require(sourceReady.ok(),sourceReady.message);
    sourcePrepared=sourceReady.ok()&&check.failures==0;
    }catch(const std::exception& error){check.Require(false,error.what());}
    if(!check.Agree(sourcePrepared))return false;
    SourceLifecycle::CommitNativeRestore(&preparedSource);
    if(!check.Status(SourceLifecycle::InitializeNativeFresh()))return false;
    check.Status(SourceLifecycle::ValidateNativeState(SourceLifecycle::CaptureNativeState()));
    Background::PublishModelOwnedSnapshot(Background::Provider::Analytic,epoch,epoch+duration,
        "native birth operator R0; actual radial material owner; frozen departure composition");
    Background::PrepareSnapshotForParticleStep();
    return check.failures==0;
  }
  static FL::Geometry::Status Put(std::vector<char>* bytes,const PIC::Datum::cDatumStored& datum,
      const double* value,int count,bool required=true) {
    if(datum.offset<0)return required?FL::Geometry::Status::Error("birth material datum is unallocated"):FL::Geometry::Status{};
    if(datum.length!=count||datum.offset+count*sizeof(double)>bytes->size())
      return FL::Geometry::Status::Error("birth material datum layout mismatch");
    std::memcpy(bytes->data()+datum.offset,value,count*sizeof(double));return {};
  }
  bool Publish() {
    auto& line=FL::FieldLinesAll[0];auto& store=Background::SnapshotStore::Instance();
    bool writable=true;
    try{store.AssertProviderMayWrite(Background::Provider::Analytic);}
    catch(const std::exception& error){check.Require(false,error.what());writable=false;}
    if(!check.Agree(writable))return false;
    const auto& c=Initialization::Active();FL::Geometry::Context context;
    context.departureS=epoch;context.arrivalS=epoch+duration;context.centerM=wind.originM;
    context.minimumRadiusM=c.innerRadiusM;context.maximumRadiusM=c.outerRadiusM;
    context.toleranceM=std::max(0.001,256*std::numeric_limits<double>::epsilon()*c.outerRadiusM);
    context.maximumInsertions=100000;context.expectedGeneration=line.GeometryGeneration;
    FL::cFieldLine::GeometryPlan plan;
    // Capture authoritative D separately. At epoch 1e12, arrival-departure is
    // about .010009765625; it MUST NOT replace the declared physical .01 s.
    const auto prepared=line.PrepareGeometryStep(context,
      [&](std::size_t,const Point& a,double,double,Point* x){*x=Map(a,duration);return FL::Geometry::Status{};},
      [&](const Point& first,const Point& chord,double distance,double at,Point* x){
        FL::Geometry::TraceControl controls;controls.maximumStepM=distance/8;
        controls.toleranceM=context.toleranceM;controls.maximumAttempts=100000;
        return FL::Geometry::TraceAlongField(first,chord,distance,at,
          [&](const Point& point,double,Point* b){const auto a=Map(point,-duration);*b=PushedB(a,InitialB(a),duration);return FL::Geometry::Status{};},controls,x);
      },[&](int index,const Point& x,double,std::vector<char>* data){
        data->assign(FL::cFieldLineVertex::StoredDataBytes(),0);
        Point a,b0;double n0=0;
        if(index>=0) {
          auto* vertex=line.GetVertex(index);
          std::memcpy(data->data(),vertex->GetAssociatedDataBufferPointer(),data->size());
          a=departure[index];vertex->GetMagneticField(b0.data());vertex->GetPlasmaDensity(n0);
        }else {
          // New computational vertices sample the declared initial continuum
          // on the inverse material map; no neighbouring history or particles
          // are cloned. The existing source rules remain their only births.
          a=Map(x,-duration);b0=InitialB(a);
          Background::AnalyticRadial::Sample before;
          const auto status=Background::AnalyticRadial::Evaluate(wind,a.data(),&before);
          if(!status.ok())return FL::Geometry::Status::Error(status.message);
          n0=5e6*std::pow(1.495978707e11/before.radiusM,2);
        }
        Background::AnalyticRadial::Sample previous,current;
        auto status=Background::AnalyticRadial::Evaluate(wind,a.data(),&previous);
        if(!status.ok())return FL::Geometry::Status::Error(status.message);
        status=Background::AnalyticRadial::Evaluate(wind,x.data(),&current);
        if(!status.ok())return FL::Geometry::Status::Error(status.message);
        const double radialScale=1+wind.speedMPerS*duration/previous.radiusM;
        const double n=n0/(radialScale*radialScale);const auto b=PushedB(a,b0,duration);
        for(const auto result:{Put(data,FL::DatumAtVertexMagneticField,b.data(),3),
            Put(data,FL::DatumAtVertexPlasmaVelocity,current.velocityMPerS.data(),3),Put(data,FL::DatumAtVertexPlasmaDensity,&n,1),
            Put(data,FL::DatumAtVertexPrevious::DatumAtVertexMagneticField,b0.data(),3),
            Put(data,FL::DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity,previous.velocityMPerS.data(),3),
            Put(data,FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&n0,1)})if(!result.ok())return result;
        // E is an allocated diagnostic in this application. When present it
        // must describe the same ideal state, not retain the departure value
        // after B changes. No thermal or wave-evolution closure is inferred.
        for(bool previousState:{false,true}) {
          const auto& velocity=previousState?previous.velocityMPerS:current.velocityMPerS;
          const auto& magnetic=previousState?b0:b;
          const double e[3]={velocity[2]*magnetic[1]-velocity[1]*magnetic[2],
              velocity[0]*magnetic[2]-velocity[2]*magnetic[0],velocity[1]*magnetic[0]-velocity[0]*magnetic[1]};
          const auto result=Put(data,previousState?FL::DatumAtVertexPrevious::DatumAtVertexElectricField:FL::DatumAtVertexElectricField,e,3,false);
          if(!result.ok())return result;
        }
        return FL::Geometry::Status{};
      },&plan,[&](int index,const Point&,const Point&,double,std::vector<char>* data){
        data->assign(FL::cFieldLineSegment::StoredDataBytes(),0);
        if(index>=0)std::memcpy(data->data(),line.GetSegment(index)->GetAssociatedDataBufferPointer(),data->size());
        return FL::Geometry::Status{};
      });
    if(!check.Agree(prepared.ok())){check.Require(prepared.ok(),prepared.message);return false;}
    // Stage background publication before committing geometry so allocation
    // failures cannot leave a new chart carrying old epoch/generation metadata.
    const Background::BackgroundSnapshot snapshot(Background::Provider::Analytic,
        Background::Ownership::ModelOwned,context.arrivalS,context.arrivalS,context.arrivalS+duration,
        line.GeometryGeneration+2,Background::CurrentConfigurationFingerprint(),
        "native birth completed radial material map; frozen departure particle operator",epoch,context.arrivalS);
    Background::SnapshotStore::Publication publication;bool publicationReady=true;
    try{publication=store.PreparePublication(snapshot);}
    catch(const std::exception& error){check.Require(false,error.what());publicationReady=false;}
    // Allocation/metadata preparation is rank-local and may throw. EVERY rank
    // must agree before the first geometry mutation, rather than discovering a
    // peer failure after some ranks have already committed the new chart.
    if(!check.Agree(publicationReady))return false;
    const auto committed=line.CommitGeometryStep(&plan);
    if(!check.Agree(committed.ok())){check.Require(committed.ok(),committed.message);return false;}
    store.CommitPublication(&publication);
    const int prefix=static_cast<int>(plan.Positions().insertedPositions.size());
    for(std::size_t index=0;index<departure.size()-plan.Positions().removedBack;++index) {
      const auto expected=Map(departure[index],duration);auto* vertex=line.GetVertex(static_cast<int>(index)+prefix);
      check.Require(expected!=departure[index],"material fixture failed to move an existing vertex, including an endpoint");
      for(unsigned d=0;d<3;++d)check.Require(vertex->GetX()[d]==expected[d],"native material transaction lost the declared full-duration vertex map");
      double n;vertex->GetPlasmaDensity(n);const auto r=FL::Geometry::Norm(FL::Geometry::Difference(expected,wind.originM));
      const double density=5e6*std::pow(1.495978707e11/r,2);
      check.Require(std::fabs(n-density)<=32*std::numeric_limits<double>::epsilon()*density,"material density violates independent radial continuity");
      // Independent tensor-vector contraction grades the stored field against
      // the Cauchy induction law. It does not call PushedB to generate its own
      // expected value; parallel and transverse deformation are both retained.
      Point previousB,currentB;vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexMagneticField,previousB.data());
      vertex->GetMagneticField(currentB.data());
      const double r0=FL::Geometry::Norm(FL::Geometry::Difference(departure[index],wind.originM));
      const double scale=r/r0,jacobian=scale*scale;
      for(unsigned row=0;row<3;++row) {
        double oracle=0;
        for(unsigned column=0;column<3;++column) {
          const double projector=(departure[index][row]-wind.originM[row])*(departure[index][column]-wind.originM[column])/(r0*r0);
          oracle+=(projector+scale*((row==column?1.:0.)-projector))*previousB[column]/jacobian;
        }
        const double magnitude=FL::Geometry::Norm(currentB);
        check.Require(std::fabs(currentB[row]-oracle)<=64*std::numeric_limits<double>::epsilon()*magnitude,
            "completed material magnetic field violates independent full-tensor ideal induction");
      }
    }
    check.Require(line.GeometryGeneration==1&&line.GeometryEpochS==context.arrivalS,
        "birth completed material publication lost owner epoch/generation");
    return true;
  }
};

// Each role is a freshly initialized, initially empty process. In particular,
// an operator control never changes live GlobalTimeStep and never pretends to
// complete a shorter physical host/sampling interval. Its state is compared
// against the pending-birth operator on the identical frozen parent chart.
class Fixture {
 public:
  Checks& check;
  T::Coefficient::Configuration originalCoefficient=T::Coefficient::ActiveConfiguration(),typed=originalCoefficient;
  Mover::ProductionMover originalMover=Mover::CurrentProductionMover();
  bool originalCoupling=W::ParticleCouplingMode;
  bool saved=false;
  BirthNativeRequest request;
  std::unique_ptr<RadialMaterial> radial;
  std::string path;
  int segmentIndex=3;
  double initialParallel=1.5e5,initialNormal=2.0e5;
  std::uint64_t identity=0;
  SourceLifecycle::Event event;
  explicit Fixture(Checks& checks):check(checks) {
    request=ReadBirthNativeRequest(GetExecutionContext().testId,GetExecutionContext().inputPath);
    check.Require(request.selected,"native birth case requires an explicit reviewed startup/role request");
    if(!request.selected)return;
    path=request.prefix+".original";
    typed.source=T::Coefficient::SourceMode::Prescribed;
    typed.pitchAngle=T::Coefficient::PitchAngleKind::Constant;typed.constantDmumuPerS=0.2;
    typed.spatial=T::Coefficient::SpatialKind::FromPitchAngle;
    typed.meanFreePath=request.material?T::Coefficient::MeanFreePathKind::FromSpatial:T::Coefficient::MeanFreePathKind::Tenishev2005;
    if(!check.Status(T::Coefficient::SetActiveConfiguration(typed)))return;
    // Only the existing qualified reduced owner has a material checkpoint
    // codec. Radial timing roles are isolated processes, not unsupported
    // material Save/Load tests disguised as cleanup mementos.
    saved=request.material?check.Status(PhysicsCheckpoint::SaveCollective(path,0)):true;
    event.campaign=TestSeed(5);event.index=UINT64_C(0x42495254484e);
    event.channel.kind=SourceLifecycle::Kind::FixedBeginning;
    identity=SourceLifecycle::ParticleIdentity(event,17);
    auto random=SourceLifecycle::Random(event,17,3);
    initialParallel+=1.0e5*random.UniformOpen01();
  }
  void DeleteTestParticles() {
    std::set<long int> particles;
    for(int line=0;line<FL::nFieldLine;++line)
      for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext()) {
        for(const auto head:{segment->FirstParticleIndex,segment->tempFirstParticleIndex.load()})
          for(auto particle=head;particle!=-1;) {
            if(!particles.insert(particle).second)break;
            particle=PB::GetNext(particle);
          }
        segment->FirstParticleIndex=-1;segment->tempFirstParticleIndex=-1;
      }
    for(auto particle:particles)PB::DeleteParticle(particle);
  }
  ~Fixture() {
    DeleteTestParticles();
    if(A::NativeDriverWindowActive()) {
      B::DriverWindow window;
      if(check.Status(A::GetNativeDriverWindow(&window))) {
        PIC::SimulationTime::SetInitialValue(window.endEpochS);
        check.Status(A::CloseNativeDriverWindow());
      }
    }
    check.Require(A::PendingWaveContributionsEmpty(),"birth fixture left physical wave work for a second owner");
    Mover::SelectProductionMover(originalMover);W::ParticleCouplingMode=originalCoupling;
    check.Status(T::Coefficient::SetActiveConfiguration(typed));
    if(saved&&request.material) {std::uint64_t completed=0;check.Status(PhysicsCheckpoint::LoadCollective(path,&completed));}
    check.Status(T::Coefficient::SetActiveConfiguration(originalCoefficient));
    if(PIC::ThisThread==0)std::remove(path.c_str());
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  bool Prepare(Mover::ProductionMover selected,double epochS,double hostDtS) {
    if(!saved)return false;
    Mover::SelectProductionMover(selected);W::ParticleCouplingMode=false;
    check.Require(PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized,
        "native birth startup did not initialize fixed global dt");
    for(int species=0;species<PIC::nTotalSpecies;++species)
      check.Require(PIC::ParticleWeightTimeStep::GlobalTimeStep[species]==hostDtS,
          "native birth role dt differs from its initialized scientific startup; live dt changes are forbidden");
    if(check.failures)return false;
    if(request.material) {
      check.Require(Background::ConfiguredProvider()==Background::Provider::ReducedShock&&
          FL::FieldLinesAll[0].ExternallyManagedGeometry&&FL::FieldLinesAll[0].GeometryEpochS==epochS,
          "checkpoint birth role lacks the actual registered reduced owner at its declared epoch");
      Background::PrepareSnapshotForParticleStep();return check.failures==0;
    }
    radial.reset(new RadialMaterial(check));return radial->Install(epochS);
  }
  bool PublishMaterial() {
    if(!radial)return true;
    PIC::SimulationTime::SetInitialValue(radial->epoch+radial->duration);
    const bool result=radial->Publish();if(result)VerifyCartesianViews(false);return result;
  }
  std::string DepartureImage() const {
    std::ostringstream metadata;const auto snapshot=Background::SnapshotStore::Instance().Current();
    metadata<<std::setprecision(17)<<"BIRTH-R0-physical-chart-v2;"<<snapshot->epoch_seconds()<<';'
        <<snapshot->valid_from_seconds()<<';'<<snapshot->valid_until_seconds()<<';'
        <<snapshot->field_line_generation()<<';'<<snapshot->configuration_fingerprint()<<';'
        <<Background::ProviderName(snapshot->provider())<<';'<<Background::OwnershipName(snapshot->ownership())<<';'
        <<T::Coefficient::ConfigurationFingerprint(T::Coefficient::ActiveConfiguration())<<';';
    std::string image=metadata.str();
    for(int line=0;line<FL::nFieldLine;++line) {
      const auto& owner=FL::FieldLinesAll[line];
      image+=owner.GetGeometryOwnerIdentity();image+=';';
      image.append(reinterpret_cast<const char*>(&owner.GeometryEpochS),sizeof(owner.GeometryEpochS));
      image.append(reinterpret_cast<const char*>(&owner.GeometryGeneration),sizeof(owner.GeometryGeneration));
      image.append(reinterpret_cast<const char*>(&owner.FirstMaterialLabel),sizeof(owner.FirstMaterialLabel));
      for(auto* vertex=FL::FieldLinesAll[line].GetFirstVertex();vertex;vertex=vertex->GetNext()) {
        const auto* x=vertex->GetX();image.append(reinterpret_cast<const char*>(x),3*sizeof(double));
        image.append(vertex->GetAssociatedDataBufferPointer(),FL::cFieldLineVertex::StoredDataBytes());
      }
      for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext())
        image.append(segment->GetAssociatedDataBufferPointer(),FL::cFieldLineSegment::StoredDataBytes());
    }
    return image;
  }
  std::string ArtifactPath(const std::string& role,const std::string& phase) const {
    return request.prefix+"."+role+"."+phase+".rank"+std::to_string(PIC::ThisThread)+".bin";
  }
  bool WriteArtifact(const std::string& phase,const std::string& bytes) {
    const auto target=ArtifactPath(request.role,phase),temporary=target+".tmp";
    std::ofstream out(temporary,std::ios::binary|std::ios::trunc);
    out.write(bytes.data(),bytes.size());out.flush();bool valid=out.good();out.close();
    valid=valid&&std::rename(temporary.c_str(),target.c_str())==0;
    check.Require(valid,"cannot transactionally write independent birth operator artifact");return check.Agree(valid);
  }
  std::string ReadArtifact(const std::string& role,const std::string& phase) {
    std::ifstream in(ArtifactPath(role,phase),std::ios::binary);
    check.Require(in.good(),"missing independent birth role artifact: "+ArtifactPath(role,phase));
    return std::string(std::istreambuf_iterator<char>(in),{});
  }
  bool Spawn(const B::Record& birth,std::uint64_t macro=0,int allocatingRank=-1,double weight=1.0) {
    auto* segment=FL::FieldLinesAll[0].GetSegment(segmentIndex);bool valid=true;
    const int allocationOwner=allocatingRank<0?segment->Thread:allocatingRank;
    if(allocationOwner==PIC::ThisThread) {
      const auto particle=PB::GetNewParticle(segment->FirstParticleIndex);
      PB::SetI(0,particle);PB::SetFieldLineId(0,particle);
      PB::SetFieldLineCoord(segmentIndex+0.45,particle);
      // Synthetic input is explicitly local-plasma-frame gyrotropic velocity;
      // this test supplies no Maxwellian/lab-frame distribution or conversion.
      PB::SetVParallel(initialParallel,particle);PB::SetVNormal(initialNormal,particle);
      PB::SetIndividualStatWeightCorrection(weight,particle);
      double x[3];FL::FieldLinesAll[0].GetCartesian(x,PB::GetFieldLineCoord(particle));PB::SetX(x,particle);
      valid=check.Status(A::InitializeParticleTransportStateWithIdentity(particle,
          SourceLifecycle::ParticleIdentity(event,17+macro),birth));
    }
    return check.Agree(valid);
  }
  std::string ParticleImage(bool includeBirthHistory=false)const {
    std::string image;
    for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
      for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        std::string record(reinterpret_cast<const char*>(PB::GetParticleDataPointer(particle)),PB::ParticleDataLength);
        auto* raw=reinterpret_cast<PB::byte*>(&record[0]);PB::SetNext(-1,raw);PB::SetPrev(-1,raw);
        // Compare identical phase-space/RNG/hazard/identity observables while
        // intentionally distinguishing initial-condition vs source history.
        if(!includeBirthHistory)
          std::fill(record.begin()+Offset::BirthRecord,record.begin()+Offset::BirthRecord+sizeof(B::Record),0);
        image+=record;
      }
    return image;
  }
  bool FixedPhysicalDt() {
    bool valid=PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized;
    for(int species=0;species<PIC::nTotalSpecies;++species)
      valid=valid&&PIC::ParticleWeightTimeStep::GlobalTimeStep[species]==0.01;
    check.Require(valid,"operator role changed the initialized .01-second physical AMPS timestep");return check.Agree(valid);
  }
  void DirectOperator(double dtS) {
    // Reuse the maintained native dispatch/join protocol, omitting sampling,
    // host clocks and source windows. The explicit duration is an independent
    // ControlledDiagnostic operator query, not a shortened physical AMPS step.
    // Snapshot handles before dispatch: movers can change segment membership.
    std::vector<long int> particles;
    for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext()) {
      for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle))particles.push_back(particle);
      segment->FirstParticleIndex=-1;
    }
    for(auto particle:particles) {
      PB::SetNext(-1,particle);PB::SetPrev(-1,particle);
      Mover::DispatchProductionMover(particle,dtS,nullptr);
    }
    for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
      segment->FirstParticleIndex=segment->tempFirstParticleIndex.exchange(-1);
    PIC::ParallelFieldLines::ExchangeFieldLineParticles();
  }
  bool Move(double readEpochS) {
    auto read=Background::SnapshotStore::Instance().BeginParticleRead(readEpochS);
    PIC::Mover::FieldLine::MoveParticles();
    return true;
  }
  void VerifyCartesianViews(bool requireNonzeroMotion) {
    int moved=0;
    for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
      for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        const double coordinate=PB::GetFieldLineCoord(particle);
        double derived[3],stored[3];FL::FieldLinesAll[0].GetCartesian(derived,coordinate);PB::GetX(stored,particle);
        for(unsigned component=0;component<3;++component) {
          check.maximumCartesianViewResidualM=std::max(check.maximumCartesianViewResidualM,std::abs(derived[component]-stored[component]));
          check.Require(std::isfinite(stored[component])&&stored[component]==derived[component],
              "native mover committed field-line coordinate with stale derived Cartesian position");
        }
        if(coordinate!=segmentIndex+0.45)moved=1;
      }
    int any=0;MPI_Allreduce(&moved,&any,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    if(requireNonzeroMotion)check.Require(any!=0,"Cartesian-view regression requires actual nonzero coordinate motion");
  }
  void RejectNonfiniteCartesianCommit(double readEpochS,double hostDurationS) {
    std::vector<A::ParticleContext> contexts;
    {
      auto read=Background::SnapshotStore::Instance().BeginParticleRead(readEpochS);
      for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
        for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
          A::ParticleContext context;double effective=0.0;
          if(!check.Status(A::LoadParticle(particle,&context))||!check.Status(A::PrepareParticleStep(&context,hostDurationS,&effective)))continue;
          // A valid privately staged phase/RNG/birth candidate makes partial
          // publication observable. No stochastic history or pending flag may
          // escape if its coordinate's derived Cartesian view is invalid.
          context.state.coordinate+=0.01;
          T::KeyedRandomStream random(A::CampaignRandomSeed,context.stableParticleId,7,0);
          std::uint64_t intervals=0;
          if(!check.Status(A::ResumeRandomStream(context,A::StochasticOperator::Parker,&random,&intervals)))continue;
          random.UniformOpen01();
          if(!check.Status(A::StageRandomStream(&context,A::StochasticOperator::Parker,random,intervals+1)))continue;
          contexts.push_back(context);
        }
    }
    // Fault injection is outside every particle-read phase. Only this fixture
    // owns its vertices/particles, and all bytes are restored before another
    // read, source event or checkpoint. It is not a geometry update policy.
    Background::SnapshotStore::Instance().AssertProviderMayWrite(Background::Provider::Analytic);
    for(auto& context:contexts) {
      auto* segment=context.segment;const auto oldTemp=segment->tempFirstParticleIndex.load();
      check.Require(oldTemp==-1,"Cartesian rejection fixture requires an empty temporary attachment list");
      if(oldTemp!=-1)continue;
      double savedVertex[3];segment->GetBegin()->GetX(savedVertex);
      double badVertex[3]={std::numeric_limits<double>::quiet_NaN(),savedVertex[1],savedVertex[2]};
      const std::string before(reinterpret_cast<const char*>(context.data),PB::ParticleDataLength);
      segment->GetBegin()->SetX(badVertex);
      const auto rejected=A::CommitAndAttach(context);
      const bool unchanged=std::memcmp(context.data,before.data(),PB::ParticleDataLength)==0&&segment->tempFirstParticleIndex.load()==oldTemp;
      // A fail-before build still needs safe teardown/continued assertions:
      // restore the deliberately corrupted geometry and the old implementation's
      // observable partial publication, including its attachment head.
      segment->GetBegin()->SetX(savedVertex);
      std::memcpy(context.data,before.data(),PB::ParticleDataLength);segment->tempFirstParticleIndex=oldTemp;
      check.Require(!rejected.ok()&&rejected.code==T::StatusCode::InvalidParticleState&&unchanged,
          "nonfinite derived Cartesian position did not reject before phase/RNG/birth/list publication");
    }
  }
};

Result Finish(const Checks& check,const char* success) {
  unsigned failures=0;MPI_Allreduce(&check.failures,&failures,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  Result result;result.status=failures?Status::Fail:Status::Pass;
  result.message=failures?(check.message.empty()?"native birth fixture failed on another MPI rank":check.message):success;
  result.hasSeed=true;result.seed=TestSeed(5);
  result.configuration={"evidence_level=native-PIC-three-mover-timing-and-persistence",
      "run_purpose=controlled-software-diagnostic;source_physics_qualified=false",
      "input_velocity_frame=explicit-local-plasma;Maxwellian-sampler=not-exercised",
      "mutable_clock=role-explicit;departure=explicit-frozen-window;finite-snapshot",
      "compiled_timestep_mode="+std::to_string(_SIMULATION_TIME_STEP_MODE_)+
          ";compiled_counter_mode="+std::to_string(_PIC_GLOBAL_TIME_COUNTER_MODE_)+
          ";compiled_coupler_mode="+std::to_string(_PIC_COUPLER_MODE_)};
  result.configuration.push_back("external_campaign_replicate_index="+std::to_string(GetExecutionContext().replicateIndex)+";native_role_invocations=1");
  result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0,"<=","count"});
  double largest=0.0;MPI_Allreduce(&check.maximumCartesianViewResidualM,&largest,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  result.metrics.push_back({"maximum_cartesian_view_residual",largest,0,"==","m"});
  return result;
}
bool Supported(Result* result) {
  if(GetExecutionContext().replicateCount>1) {
    result->status=Status::Error;
    result->message="native birth repetitions require run_birth_operator_stage_tests.py --repeat-count so every operator role has a fresh process and immutable initialized timestep";
    return false;
  }
  const auto clock=A::NativeDriverClockCapability();
  if(!clock.ok()){result->status=Status::Error;result->message=clock.message;return false;}
  const auto request=ReadBirthNativeRequest(GetExecutionContext().testId,GetExecutionContext().inputPath);
  // The outer routine launcher expands repetitions into fresh role processes.
  // Rebind ONLY the existing test-campaign context before its first TestSeed:
  // this reproduces the registry's exact versioned schedule, not seed+index.
  // Registry::RunCampaign restores its saved context at the process boundary.
  auto context=GetExecutionContext();
  if(request.campaignReplicateIndex&&!context.seedActive) {
    result->status=Status::Error;result->message="repeated birth role requires an explicit campaign seed";return false;
  }
  context.caseSeed=DeriveTestCaseSeed(context.campaignSeed,context.testId,request.campaignReplicateIndex);
  context.replicateIndex=request.campaignReplicateIndex;SetExecutionContext(context);
  const auto expected=request.material?Background::Provider::ReducedShock:Background::Provider::Analytic;
  if(!request.selected||Background::ConfiguredProvider()!=expected||!W::ActiveFlag||
      !FL::FieldLinesAll||FL::nFieldLine!=1||FL::FieldLinesAll[0].GetTotalSegmentNumber()<7) {
    result->status=Status::Error;result->message="native birth role requires its reviewed startup, one actual supported line with at least seven segments, and allocated wave storage";return false;
  }
  long int particles=PB::GetAllPartNum(),total=0;
  MPI_Allreduce(&particles,&total,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(total||A::NativeDriverWindowActive()) {
    result->status=Status::Error;result->message="run native birth fixture in a fresh empty test process outside a host window";return false;
  }
  return true;
}
Result RemainingDuration(Mover::ProductionMover selected) {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;double remainingS=0.0025;std::string role;std::vector<std::string> artifacts;
  {
    Fixture fixture(check);role=fixture.request.role;
    const double start=1.0,hostDt=0.01;
    const bool actual=role=="actual";
    if(fixture.Prepare(selected,start,hostDt)&&fixture.FixedPhysicalDt()) {
      const auto frozen=fixture.DepartureImage();
      B::Record birth;
      bool ready=true;
      if(actual) {
        ready=check.Status(A::FreezeNativeDriverWindow());
        if(ready)ready=check.Status(A::MakeUniformBirthRecord(0.75,&birth));
        remainingS=birth.firstMoveDurationS;
        check.Require(remainingS==0.0025,"native scheduled birth did not retain exact .0025-second remainder");
      }else ready=check.Status(B::InitialCondition(start,&birth));
      if(ready&&fixture.Spawn(birth)) {
        // Initial phase space, immutable identity and complete RNG/hazard
        // history are part of R0 parity, in addition to all consumed native
        // vertex/segment state and declared provider/coefficient metadata.
        const auto initial=frozen+fixture.ParticleImage();fixture.WriteArtifact("R0",initial);
        // All processes retain exactly the same initialized .01 physical dt.
        // Controls use the existing explicit-duration diagnostic dispatch with
        // no source window, counter advancement or sampling. Only the actual
        // role executes the ordinary native manager/end-counter ordering.
        {auto read=Background::SnapshotStore::Instance().BeginParticleRead(start);
         if(actual) {
           PIC::SimulationTime::SetInitialValue(start+hostDt);
           PIC::Mover::FieldLine::MoveParticles();
         }else fixture.DirectOperator(role=="direct"?remainingS:hostDt);}
        if(actual)check.Status(A::CloseNativeDriverWindow());
        fixture.VerifyCartesianViews(true);const auto moved=fixture.ParticleImage();
        fixture.WriteArtifact("operator",moved);
        if(actual) {
          check.Require(initial==fixture.ReadArtifact("direct","R0")&&initial==fixture.ReadArtifact("full","R0"),
              "independent birth controls did not use identical frozen parent realization");
          check.Require(moved==fixture.ReadArtifact("direct","operator"),
              "pending native birth differs from independent identical-R0 exact-duration operator reference");
          const int differs=moved!=fixture.ReadArtifact("full","operator")?1:0;int any=0;
          MPI_Allreduce(&differs,&any,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
          check.Require(any!=0,"negative control did not distinguish full-host movement from remaining-duration movement");
          for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
            for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
              B::Record consumed;std::memcpy(&consumed,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(consumed));
              check.Require(consumed.pendingFirstMove==0&&consumed.firstMoveDurationS==remainingS,
                  "native mover lost provenance or failed atomic first-move consumption");
            }
        }
        // Publish the SAME complete parent map M_.01 after every operator
        // experiment. The controls establish an operator identity, not a new
        // physical trajectory starting at the birth epoch. Endpoint canonical
        // X and the original byte comparison must agree on this second chart.
        if(fixture.PublishMaterial()) {
          const auto completed=fixture.ParticleImage();fixture.WriteArtifact("completed",completed);
          if(actual)check.Require(completed==fixture.ReadArtifact("direct","completed"),
              "equal remaining-duration operators differ after identical full material-map publication");
        }
      }
      fixture.FixedPhysicalDt();
      for(const auto* phase:{"R0","operator","completed"})artifacts.push_back(fixture.ArtifactPath(role,phase));
    }
  }
  auto result=Finish(check,"native birth operator role preserves exact remaining duration and unchanged frozen-departure composition with real completed material publication");
  result.configuration.push_back("birth_operator_role="+role+";parent_physical_duration_s=.01;operator_remaining_s=.0025;reference_is_operator_stage=true;no_fractional_sampling=true");
  result.configuration.push_back("native_particle_operator_duration_s="+std::string(role=="full"?".01":".0025")+
      ";mutable_clock_during_operator="+(role=="actual"?"parent-end":"parent-departure")+
      ";completed_geometry_clock=parent-end;initialized_physical_dt_s=.01");
  result.metrics.push_back({"first_move_duration",remainingS,0.0025,"==","s"});result.artifacts=artifacts;return result;
}
Result Parker(){return RemainingDuration(Mover::ProductionMover::Parker);}
Result Dmumu(){return RemainingDuration(Mover::ProductionMover::FocusedTransportDiffusion);}
Result Mfp(){return RemainingDuration(Mover::ProductionMover::FocusedTransportMeanFreePath);}

Result AtomicTimingRejection() {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;double lateDuration=0;
  {
    Fixture fixture(check);const double epoch=1.0e12,dt=0.01;
    if(fixture.Prepare(Mover::ProductionMover::FocusedTransportMeanFreePath,epoch,dt)&&
        check.Status(A::FreezeNativeDriverWindow())) {
      B::Record birth;check.Status(A::MakeUniformBirthRecord(std::nextafter(1.0,0.0),&birth));
      lateDuration=birth.firstMoveDurationS;
      check.Require(lateDuration>0&&epoch+dt-lateDuration==epoch+dt,
          "late-birth fixture did not expose absolute-epoch rounding loss");
      if(fixture.Spawn(birth)) {
        fixture.RejectNonfiniteCartesianCommit(epoch,dt);
        const auto before=fixture.ParticleImage(true);
        {auto read=Background::SnapshotStore::Instance().BeginParticleRead(epoch);
         PIC::SimulationTime::SetInitialValue(epoch+dt);
         const auto rejected=A::CloseNativeDriverWindow();
         check.Require(!rejected.ok()&&rejected.message.find("during particle motion")!=std::string::npos&&
             A::NativeDriverWindowActive()&&before==fixture.ParticleImage(true),
             "closing an active particle read phase changed or consumed its source window/history");
         for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
           for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
             auto* data=PB::GetParticleDataPointer(particle);B::Record stored;
             std::memcpy(&stored,data+Offset::BirthRecord,sizeof(stored));
             auto unknown=stored;unknown.schema=0;
             std::memcpy(data+Offset::BirthRecord,&unknown,sizeof(unknown));
             A::ParticleContext incompatible;
             const auto rejectedHistory=A::LoadParticle(particle,&incompatible);
             std::memcpy(data+Offset::BirthRecord,&stored,sizeof(stored));
             check.Require(!rejectedHistory.ok()&&before==fixture.ParticleImage(true),
                 "native load guessed missing birth history or mutated the restored record");
             A::ParticleContext context;check.Status(A::LoadParticle(particle,&context));
             double effective=123;
             const auto wrong=A::PrepareParticleStep(&context,0.5*dt,&effective);
             check.Require(!wrong.ok()&&effective==123&&before==fixture.ParticleImage(true),
                 "incompatible first-duration request changed source history or output");
           }
         PIC::Mover::FieldLine::MoveParticles();}
        check.Status(A::CloseNativeDriverWindow());
        for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            B::Record consumed;T::OperatorRandomState random[3];
            const auto* data=PB::GetParticleDataPointer(particle);
            std::memcpy(&consumed,data+Offset::BirthRecord,sizeof(consumed));
            std::memcpy(random,data+Offset::OperatorRandomStates,sizeof(random));
            check.Require(!consumed.pendingFirstMove&&random[2].acceptedIntervals>0,
                "strictly positive late birth was discarded below the adaptive floor");
          }
      }
      fixture.PublishMaterial();
    }
  }
  auto result=Finish(check,"read-phase close and incompatible requests reject without mutation; a representable late birth still executes an actual MFP interval below the adaptive refinement floor");
  result.configuration.push_back("fixture_start_epoch_s=1000000000000;new_initial_realization=true;physical_duration_s=.01;relative_remainder_preserved=true");
  result.metrics.push_back({"preserved_late_first_move_duration",lateDuration,0,">","s"});return result;
}

Result CloneSplitMigration() {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;std::uint64_t splitPopulation=0;double declaredEpochS=0;
  {
    Fixture fixture(check);const double epoch=1.0,dt=0.01;
    declaredEpochS=epoch;
    if(fixture.Prepare(Mover::ProductionMover::FocusedTransportMeanFreePath,epoch,dt)&&
        check.Status(A::FreezeNativeDriverWindow())) {
      B::Record birth;check.Status(A::MakeUniformBirthRecord(0.75,&birth));
      auto* segment=FL::FieldLinesAll[0].GetSegment(fixture.segmentIndex);
      const int owner=segment->Thread;
      const int originalRank=PIC::nTotalThreads>1?(owner+1)%PIC::nTotalThreads:owner;
      bool created=true;
      for(unsigned macro=0;macro<5;++macro)
        created=fixture.Spawn(birth,macro,originalRank,macro==0?8.0:1.0)&&created;
      if(created) {
        if(PIC::ThisThread==originalRank) {
          // Generic cloning is a transfer/copy primitive. It must preserve
          // physical birth age and pending time; only the actual statistical
          // split callback may assign a child lineage/future random stream.
          const auto parent=segment->FirstParticleIndex,copy=PB::GetNewParticle();
          PB::CloneParticle(copy,parent);
          const auto* a=PB::GetParticleDataPointer(parent);const auto* b=PB::GetParticleDataPointer(copy);
          check.Require(std::memcmp(a+Offset::BirthRecord,b+Offset::BirthRecord,sizeof(B::Record))==0&&
              std::memcmp(a+Offset::OperatorRandomStates,b+Offset::OperatorRandomStates,3*sizeof(T::OperatorRandomState))==0,
              "generic CloneParticle changed pending physical birth or stochastic history");
          PB::DeleteParticle(copy);
        }
        // Place source records on another rank deliberately, then execute the
        // actual native MPI migration before their first movement. The full
        // particle payload, including all64 birth bytes, must reach the owner.
        PIC::ParallelFieldLines::ExchangeFieldLineParticles();
        int localCount=0;bool migrationValid=true;
        for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
          ++localCount;B::Record migrated;
          std::memcpy(&migrated,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(migrated));
          migrationValid=migrationValid&&PIC::ThisThread==owner&&std::memcmp(&birth,&migrated,sizeof(birth))==0;
        }
        check.Require(localCount==(PIC::ThisThread==owner?5:0),"native exchange changed source population/owner");
        if(check.Agree(migrationValid)) {
          {auto read=Background::SnapshotStore::Instance().BeginParticleRead(epoch);
           PIC::SimulationTime::SetInitialValue(epoch+dt);PIC::Mover::FieldLine::MoveParticles();}
          check.Status(A::CloseNativeDriverWindow());
          A::RegisterPopulationControlHooks();A::BeginPopulationControl();
          if(PIC::ThisThread==owner) {
            double parallelRange=0,normalRange=0;
            // Existing exact half-weight core controller: five representatives
            // become six. This is a statistical copy after the closed physical
            // step, not a new physical source event or a reset birth clock.
            segment->WeightedParticleSplitting(0,segment->FirstParticleIndex,2,2,2,1.0,
                parallelRange,normalRange,5,7,2);
          }
          auto expected=birth;expected.pendingFirstMove=0;
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            ++splitPopulation;B::Record inherited;
            std::memcpy(&inherited,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(inherited));
            check.Require(std::memcmp(&inherited,&expected,sizeof(inherited))==0,
                "actual exact split invented a new birth age or lost consumed source-window provenance");
          }
        }
      }
      fixture.PublishMaterial();
    }
  }
  std::uint64_t total=0;MPI_Allreduce(&splitPopulation,&total,1,MPI_UINT64_T,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  check.Require(total==6,"actual exact split did not retain six physical-history representatives");
  auto result=Finish(check,"generic clone, actual native MPI exchange and actual exact half-weight split preserve physical birth provenance and pending/consumed history");
  result.configuration.push_back("fixture_start_epoch_s="+std::to_string(declaredEpochS)+";isolated_material_role=true;no_live_dt_change=true");
  result.configuration.push_back(std::string("actual_interrank_transfer=")+(PIC::nTotalThreads>1?"true":"false;requires-MPI2-for-transfer-gate"));
  result.metrics.push_back({"post_split_representatives",static_cast<double>(total),6,"==","count"});return result;
}

std::string ReadRootFile(const std::string& path) {
  if(PIC::ThisThread!=0)return std::string();
  std::ifstream input(path.c_str(),std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
}
Result CompletedCheckpointRejection() {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;double declaredEpochS=0;long int survivingSources=0;
  {
    Fixture fixture(check);const auto originalFile=ReadRootFile(fixture.path);
    const double epoch=Background::SimulationTimeSeconds(),dt=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
    declaredEpochS=epoch;
    if(fixture.Prepare(Mover::ProductionMover::FocusedTransportMeanFreePath,epoch,dt)&&
        check.Status(A::FreezeNativeDriverWindow())) {
      B::DriverWindow window;bool ready=check.Status(A::GetNativeDriverWindow(&window));
      // This host has an enabled, qualified reduced source. Therefore records
      // must come from its ACTUAL injection API, including immutable cohort
      // energy and source ledgers. A generic manufactured particle would fail
      // the source-history guard before reaching the pending-birth negative;
      // fabricating SourceRecord bytes would invalidate that checkpoint test.
      if(ready)ready=check.Status(ReducedParticles::Native::Inject(window));
      long int localBorn=PB::GetAllPartNum(),born=0;
      MPI_Allreduce(&localBorn,&born,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      check.Require(born>0,"qualified checkpoint birth fixture injected no actual source records");
      if(ready&&born>0) {
        const auto before=fixture.ParticleImage(true);std::uint64_t completed=123;
        const auto activeSave=PhysicsCheckpoint::SaveCollective(fixture.path,0);
        const auto activeLoad=PhysicsCheckpoint::LoadCollective(fixture.path,&completed);
        check.Require(!activeSave.ok()&&!activeLoad.ok()&&
            activeSave.message.find("active source birth window")!=std::string::npos&&
            activeLoad.message.find("active source birth window")!=std::string::npos&&completed==123&&
            before==fixture.ParticleImage(true)&&A::NativeDriverWindowActive(),
            "active-window checkpoint saved/loaded partial source movement or mutated physics");
        check.Require(ReadRootFile(fixture.path)==originalFile,"rejected active-window save changed the existing checkpoint destination");
        {auto read=Background::SnapshotStore::Instance().BeginParticleRead(epoch);
         PIC::SimulationTime::SetInitialValue(epoch+dt);PIC::Mover::FieldLine::MoveParticles();
         // This fixture calls the actual manager directly, so it must also
         // execute PIC::TimeStep's ordinary joined postmover migration. Genuine
         // source records can traverse MPI-owned segments over a full300 s;
         // omitting this phase leaves foreign-owned lists and invalidates both
         // the inventory ledger and the completed checkpoint. Migration copies
         // the complete real source/birth/RNG payload, without reinitializing it.
         PIC::ParallelFieldLines::ExchangeFieldLineParticles();}
        check.Status(A::CloseNativeDriverWindow());
        // Preserve the real source inventory/loss accounting and complete its
        // material chart before claiming a completed-state checkpoint. The
        // unmodified absorbing boundary may delete genuine returning sources;
        // the diagnostic requires a positive surviving cohort, never disables
        // absorption or substitutes an arbitrary exact survivor count.
        check.Status(ReducedParticles::Native::BeforeGeometry());
        const auto advanced=ReducedShock::Native::Advance(epoch,epoch+dt);
        check.Require(advanced.ok(),advanced.message);
        check.Status(ReducedParticles::Native::AfterGeometry(epoch+dt));
        fixture.VerifyCartesianViews(false);
        const long int localSurvivors=PB::GetAllPartNum();
        MPI_Allreduce(&localSurvivors,&survivingSources,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        check.Require(survivingSources>0,"checkpoint pending-history negative has no surviving real source record");
        // A positive save proves every ordinary completed-state prerequisite
        // before changing just the pending flag. It uses a separate artifact,
        // preserving the original startup memento/destination byte oracle.
        check.Status(PhysicsCheckpoint::SaveCollective(fixture.path+".completed",1));
        std::vector<std::pair<long int,B::Record>> consumedRecords;
        for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            B::Record consumed;auto* data=PB::GetParticleDataPointer(particle);
            std::memcpy(&consumed,data+Offset::BirthRecord,sizeof(consumed));
            check.Require(!consumed.pendingFirstMove,"actual source movement left an unconsumed record before deliberate fault injection");
            consumedRecords.push_back({particle,consumed});
            ReducedParticles::Native::SourceRecord cohort;
            std::memcpy(&cohort,data+Offset::ReducedSourceRecord,sizeof(cohort));
            check.Require(cohort.schema==1&&std::isfinite(cohort.birthKineticEnergyJ)&&cohort.birthKineticEnergyJ>0,
                "actual surviving checkpoint source lacks its valid immutable birth-energy cohort");
            std::uint64_t stableId=0;std::memcpy(&stableId,data+Offset::StableParticleId,sizeof(stableId));
            // Fault injection changes ONLY the pending bit of the actual
            // preserved source history. It cannot invent an energy, origin,
            // birth epoch, source event, duration or a compatible old record.
            auto pending=consumed;pending.pendingFirstMove=1;
            std::memcpy(data+Offset::BirthRecord,&pending,sizeof(pending));
            std::uint64_t afterId=0;std::memcpy(&afterId,data+Offset::StableParticleId,sizeof(afterId));
            check.Require(stableId!=0&&afterId==stableId&&std::memcmp(data+Offset::ReducedSourceRecord,&cohort,sizeof(cohort))==0,
                "pending-bit fault injection changed immutable source identity/energy/cohort bytes");
          }
        const auto pendingImage=fixture.ParticleImage(true);
        const auto pendingSave=PhysicsCheckpoint::SaveCollective(fixture.path,1);
        check.Require(!pendingSave.ok()&&pendingSave.message.find("unconsumed scheduled first move")!=std::string::npos&&
            pendingImage==fixture.ParticleImage(true),
            "completed checkpoint accepted/consumed pending history or rejected for a different prerequisite: "+pendingSave.message);
        check.Require(ReadRootFile(fixture.path)==originalFile,"rejected pending-record save changed the checkpoint destination");
        for(const auto& record:consumedRecords)
          std::memcpy(PB::GetParticleDataPointer(record.first)+Offset::BirthRecord,&record.second,sizeof(record.second));
      }
    }
  }
  auto result=Finish(check,"completed checkpoint API rejects active windows and unconsumed actual source records without changing particles, caller state or an existing destination");
  result.configuration.push_back("fixture_start_epoch_s="+std::to_string(declaredEpochS)+";actual_reduced_source_and_material_owner=true;no_live_dt_change=true");
  result.metrics.push_back({"surviving_actual_source_records",static_cast<double>(survivingSources),0,">","count"});
  return result;
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.description=name;
  descriptor.group="birth-window-native";descriptor.initialization=std::string(id)=="BIRTHN06"?InitializationLevel::FieldLineModel:InitializationLevel::NativeStorage;descriptor.tags.push_back("birth-operator-startup-v1");
  descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Context;
  descriptor.repeatSupported=true;descriptor.tags.push_back("repeat:process-isolated");descriptor.seedPolicy="requested seed drives actual synthetic source identity and local-plasma velocity in an isolated fixed-startup role";
  descriptor.supportedBuildModes="explicit radial material operator role or approved reduced checkpoint role; initialized fixed global dt";
  descriptor.stateIsolation="fresh process per operator role; no live timestep change; actual reduced checkpoint restored only for BIRTHN06";
  descriptor.callback=callback;return descriptor;
}
} // namespace
BirthNativeRequest ReadBirthNativeRequest(const std::string& id,const std::string& path) {
  BirthNativeRequest request;
  if(id.size()!=8||id.compare(0,7,"BIRTHN0")!=0||id[7]<'1'||id[7]>'6'||path.empty())return request;
  request.material=id=="BIRTHN06";
  std::ifstream input(path);std::string magic,roleKey,prefixKey,inputKey,eventKey;unsigned version=0;
  if(!(input>>magic>>version>>roleKey>>std::quoted(request.role)>>prefixKey>>std::quoted(request.prefix)
       >>inputKey>>std::quoted(request.inputPath))||roleKey!="role"||prefixKey!="prefix"||inputKey!="startup_input"||
       request.prefix.empty()||request.prefix[0]!='/'||request.inputPath.empty()||request.inputPath[0]!='/')
    throw std::runtime_error("native birth startup requires explicit role, absolute prefix and startup_input");
  if(request.material) {
    if(magic!="SEP_CHECKPOINT_FRESH_PROCESS"||version!=2||request.role!="routine"||
        !(input>>eventKey>>std::quoted(request.eventPath))||eventKey!="event"||request.eventPath.empty()||request.eventPath[0]!='/')
      throw std::runtime_error("BIRTHN06 requires the existing reduced material checkpoint protocol2 and actual event");
  }else if(magic!="SEP_BIRTH_OPERATOR_STAGE"||version!=1||
      (request.role!="actual"&&request.role!="direct"&&request.role!="full")||
      (id>"BIRTHN03"&&request.role!="actual"))
    throw std::runtime_error("native birth operator protocol1 requires actual role or BIRTHN01-N03 direct/full controls");
  input>>std::ws;
  if(!input.eof()) {
    std::string key;unsigned long long index=0;
    if(!(input>>key>>index)||key!="campaign_replicate_index"||index>std::numeric_limits<std::uint32_t>::max())
      throw std::runtime_error("invalid native birth external campaign replicate index");
    request.campaignReplicateIndex=static_cast<std::uint32_t>(index);
    input>>std::ws;if(!input.eof())throw std::runtime_error("trailing native birth startup protocol fields");
  }
  request.selected=true;return request;
}
std::vector<Descriptor> BirthNativeDescriptors() {
  return {Make("BIRTHN01","Native Parker consumes source remaining duration",Parker),
      Make("BIRTHN02","Native Dmumu consumes source remaining duration",Dmumu),
      Make("BIRTHN03","Native MFP consumes source remaining duration",Mfp),
      Make("BIRTHN04","Native birth timing rejection and late interval",AtomicTimingRejection),
      Make("BIRTHN05","Native clone split and migration retain birth history",CloneSplitMigration),
      Make("BIRTHN06","Completed checkpoint refuses pending birth movement",CompletedCheckpointRejection)};
}
} }
