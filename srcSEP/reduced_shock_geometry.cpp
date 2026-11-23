#include "sep.h"
#include "reduced_shock_geometry.h"
#include "adapters/reduced_material_background_adapter.h"
#include "util/sep_initialization.h"
#include "util/sep_background_runtime.h"
#include "util/sep_run_configuration.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>

namespace SEP {namespace ReducedShock {namespace Native {
namespace FL=PIC::FieldLine;
namespace {
using LineView=MaterialView;
std::shared_ptr<const LineView> lineView;
Core::Status Error(const std::string& why) {return Core::Status::Failure(Core::StatusCode::InvalidState,why);}
Core::Result<std::shared_ptr<const LineView>> PrepareLineView(
    double epoch,std::uint64_t generation,const std::vector<Position>& points,
    const IntervalView* restoredInterval=nullptr) {
  using Return=Core::Result<std::shared_ptr<const LineView>>;
  if(points.size()<2)return Return::Failure(Core::StatusCode::InvalidState,"core material connectivity requires a nonempty line");
  const auto interval=restoredInterval?Core::Result<IntervalView>::Success(*restoredInterval):
      CapturePhysicalStep(epoch,std::min(EventEndS(),epoch+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]),EventIdentity());
  if(!interval.ok())return Return::Failure(interval.status.code,interval.status.message);
  auto candidate=std::make_shared<LineView>();candidate->interval=interval.value;candidate->vertices=points;
  candidate->geometryGeneration=generation;candidate->arc.push_back(0);
  for(std::size_t i=1;i<points.size();++i)candidate->arc.push_back(candidate->arc.back()+
      FL::Geometry::Norm(FL::Geometry::Difference(points[i],points[i-1])));
  // Native metadata has its own generation (cache/geometry join), distinct
  // from both the physical material labels and the core cache cadence.
  const auto material=MaterialAdapter::Capture(interval.value,generation,
      std::max(interval.value.Generation(),generation+1),Background::CurrentConfigurationFingerprint(),
      Controls(interval.value.ArrivalS()-epoch,candidate->arc[1]));
  if(!material.ok())return Return::Failure(material.status.code,material.status.message);
  candidate->material=material.value;
  const auto connected=ConnectLine(interval.value,points,epoch,Controls(0,candidate->arc[1]).positionToleranceM);
  if(!connected.ok())return Return::Failure(connected.status.code,connected.status.message);
  candidate->departureConnections=connected.value;return Return::Success(std::move(candidate));
}
Core::Status Agree(const Core::Status& local) {
  int failed=local.ok()?std::numeric_limits<int>::max():PIC::ThisThread,first=0;
  MPI_Allreduce(&failed,&first,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(first==std::numeric_limits<int>::max())return Core::Status::Success();
  char message[2048]={};
  if(PIC::ThisThread==first)std::snprintf(message,sizeof(message),"reduced material preparation rank %d: %s",first,local.message.c_str());
  MPI_Bcast(message,sizeof(message),MPI_CHAR,first,MPI_GLOBAL_COMMUNICATOR);
  return Error(message);
}
Core::Status AgreeRealization() {
  // The resolved event identity includes the full physical configuration and
  // asset bytes. Matching coordinates alone would miss rank-dependent EOS or
  // temperature inputs. Compare the complete SHA-256 identity, never a lossy
  // numerical hash, before staging any core/native candidate.
  char reference[65]={};const auto& identity=EventIdentity();
  if(PIC::ThisThread==0&&identity.size()==64)std::memcpy(reference,identity.data(),64);
  MPI_Bcast(reference,65,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  auto status=Agree(identity.size()==64&&identity==reference?Core::Status::Success():
      Error("MPI ranks selected different reduced core event/asset identities"));
  if(!status.ok())return status;
  const std::uint64_t ids[2]={Metadata().generation,FL::nFieldLine==1?FL::FieldLinesAll[0].GeometryGeneration:
      std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t minimum[2],maximum[2];
  MPI_Allreduce(ids,minimum,2,MPI_UINT64_T,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(ids,maximum,2,MPI_UINT64_T,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  for(int i=0;i<2;++i)if(minimum[i]!=maximum[i])return Error("MPI core-cache/material generations differ before preparation");
  return Core::Status::Success();
}
Core::Status Put(std::vector<char>& data,const PIC::Datum::cDatumStored& datum,
    const double* values,int length,bool required=true) {
  if(datum.offset<0)return required?Error("required reduced material datum is not allocated"):Core::Status::Success();
  if(datum.length!=length||static_cast<std::size_t>(datum.offset+length*sizeof(double))>data.size())
    return Error("reduced material datum layout mismatch");
  std::memcpy(data.data()+datum.offset,values,length*sizeof(double));return Core::Status::Success();
}
Core::Status PutState(std::vector<char>& data,const AmbientSample& sample,bool previous) {
  if(!sample.parkerExterior)return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
      "material native geometry is qualified only in the core Parker exterior; PFSS evolution requires its own induction contract");
  auto& b=previous?FL::DatumAtVertexPrevious::DatumAtVertexMagneticField:FL::DatumAtVertexMagneticField;
  auto& u=previous?FL::DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity:FL::DatumAtVertexPlasmaVelocity;
  auto& n=previous?FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity:FL::DatumAtVertexPlasmaDensity;
  auto& p=previous?FL::DatumAtVertexPrevious::DatumAtVertexPlasmaPressure:FL::DatumAtVertexPlasmaPressure;
  auto& t=previous?FL::DatumAtVertexPrevious::DatumAtVertexPlasmaTemperature:FL::DatumAtVertexPlasmaTemperature;
  for(auto status:{Put(data,b,sample.magneticFieldT.data(),3),Put(data,u,sample.velocityMPerS.data(),3),
      Put(data,n,&sample.protonNumberDensityM3,1),Put(data,p,&sample.pressurePa,1),Put(data,t,&sample.protonTemperatureK,1)})
    if(!status.ok())return status;
  // The prescribed core ambient is ideal: E=-U cross B. The focused particle
  // equations do not use E directly, but allocated diagnostic E must describe
  // this same state rather than retain a competing Eulerian initialization.
  const auto& v=sample.velocityMPerS;const auto& field=sample.magneticFieldT;
  const double e[3]={v[2]*field[1]-v[1]*field[2],v[0]*field[2]-v[2]*field[0],v[1]*field[0]-v[0]*field[1]};
  return Put(data,previous?FL::DatumAtVertexPrevious::DatumAtVertexElectricField:FL::DatumAtVertexElectricField,e,3,false);
}
}
GeometryControls Controls(double interval,double length) {
  // A coordinate-scale roundoff allowance is an explicit numerical floor,
  // not a magnetic/density floor. Both trajectory and magnetic tracing use
  // tighter accuracy than the discretized line chord. Work limits reject.
  const double scale=Initialization::HasActive()?Initialization::Active().outerRadiusM:std::max(1.,length);
  return {std::max(0.001,256*std::numeric_limits<double>::epsilon()*scale),
      std::max(std::abs(interval)/8,1e-6),length/2,100000};
}
Core::Status BuildInitialLine(std::list<SEP::cFieldLine>* output) {
  if(!output||!Initialization::HasActive()||SEP::Domain_nTotalParkerSpirals!=1)
    return Error("reduced material owner requires one explicitly configured runtime field line");
  const auto& configuration=Initialization::Active();
  if(configuration.parkerOriginM.x!=0||configuration.parkerOriginM.y!=0||configuration.parkerOriginM.z!=0)
    return Error("core reduced material geometry requires the declared zero-centred HCI origin");
  const double epoch=Background::SimulationTimeSeconds();std::string error;
  if(!Prepare(epoch,&error))return Error(error);
  const auto view=CaptureInterval(epoch,epoch,EventIdentity());if(!view.ok())return view.status;
  const Position initial{{configuration.parkerInitialPointM.x,configuration.parkerInitialPointM.y,configuration.parkerInitialPointM.z}};
  const double length=configuration.parkerLengthM/(configuration.parkerPointCount-1);
  const auto points=ReducedShock::BuildInitialLine(view.value,initial,configuration.parkerLengthM,
      configuration.parkerPointCount,epoch,Controls(0,length));
  if(!points.ok())return points.status;
  std::list<SEP::cFieldLine> staged;
  for(const auto& x:points.value) {
    const auto sample=view.value.PrimitiveAt(x,epoch);if(!sample.ok())return sample.status;
    SEP::cFieldLine point;
    std::copy(x.begin(),x.end(),point.x);std::copy(sample.value.magneticFieldT.begin(),sample.value.magneticFieldT.end(),point.B);
    std::copy(sample.value.velocityMPerS.begin(),sample.value.velocityMPerS.end(),point.U);staged.push_back(point);
  }
  output->swap(staged);return Core::Status::Success();
}
Core::Status RegisterLines(double epoch) {
  if(!Enabled())return Core::Status::Success();
  if(FL::nFieldLine!=1)return Error("first reduced material owner supports one configured physical line");
  const auto registered=FL::FieldLinesAll[0].RegisterGeometryOwner(EventIdentity()+":configured-line:0",epoch);
  return registered.ok()?Core::Status::Success():Error(registered.message);
}
Core::Status InstallBackground(double epoch) {
  if(!Enabled())return Core::Status::Success();
  const auto binding=AgreeRealization();if(!binding.ok())return binding;
  auto& store=Background::SnapshotStore::Instance();
  std::vector<std::vector<char>> staged;
  std::shared_ptr<const LineView> preparedView;
  Background::SnapshotStore::Publication publication;
  CadencePublication cadencePublication;
  Core::Status local=Core::Status::Success();bool reuseMetadata=false;
  try {
    if(FL::nFieldLine!=1)throw std::runtime_error("initial reduced publication requires one registered material line");
    const auto cadence=PrepareCadencePublication(epoch);
    if(!cadence.ok())throw std::runtime_error(cadence.status.message);
    cadencePublication=cadence.value;
    auto& line=FL::FieldLinesAll[0];
    // The actual registered owner establishes lineage, even at valid initial
    // generation zero. Reject a merely allocated/unregistered line rather than
    // inventing a positive counter or accepting a default metadata identity.
    if(!line.ExternallyManagedGeometry||line.GetGeometryOwnerIdentity()!=EventIdentity()+":configured-line:0")
      throw std::runtime_error("reduced publication requires the actual registered core material owner");
    if(line.GeometryEpochS!=epoch)throw std::runtime_error("initial background and material geometry epochs differ");
    const auto current=store.Current();
    const double until=std::min(EventEndS(),epoch+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);
    const auto generation=std::max(cadencePublication.Generation(),line.GeometryGeneration+1);
    reuseMetadata=current&&current->provider()==Background::Provider::ReducedShock&&
        current->epoch_seconds()==epoch&&current->valid_until_seconds()==until&&
        current->field_line_generation()==generation;
    store.AssertProviderMayWrite(Background::Provider::ReducedShock);
    // The candidate cadence may be newer than the live cache (especially
    // when background_dt < fixed AMPS dt). Bind to that exact staged owner,
    // whose cache generation is used by the candidate native snapshot above.
    // Reading the live owner here would join two different realizations and
    // can falsely reject a correct publication when the cache counter exceeds
    // the geometry counter. The immutable physical functions and full fixed
    // AMPS interval are unchanged; neither candidate is published yet.
    const auto view=cadencePublication.PhysicalInterval(epoch,until);
    if(!view.ok())throw std::runtime_error(view.status.message);
    std::vector<Position> points;
    for(auto* vertex=line.GetFirstVertex();vertex;vertex=vertex->GetNext()) {
      const double* x=vertex->GetX();points.push_back({{x[0],x[1],x[2]}});
      const auto sample=view.value.AmbientAt(points.back(),epoch);if(!sample.ok())throw std::runtime_error(sample.status.message);
      staged.emplace_back(FL::cFieldLineVertex::StoredDataBytes());
      std::memcpy(staged.back().data(),vertex->GetAssociatedDataBufferPointer(),staged.back().size());
      // Re-preparing cadence metadata at an unchanged physical epoch must not
      // replace the genuine previous material state. At the first epoch both
      // slots describe this initial condition, without inventing prehistory.
      if(!current||current->epoch_seconds()!=epoch) {
        if(current) {
          for(auto pair:{std::make_pair(&FL::DatumAtVertexMagneticField,&FL::DatumAtVertexPrevious::DatumAtVertexMagneticField),
              std::make_pair(&FL::DatumAtVertexPlasmaVelocity,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity),
              std::make_pair(&FL::DatumAtVertexPlasmaDensity,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity),
              std::make_pair(&FL::DatumAtVertexPlasmaPressure,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaPressure),
              std::make_pair(&FL::DatumAtVertexPlasmaTemperature,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaTemperature),
              std::make_pair(&FL::DatumAtVertexElectricField,&FL::DatumAtVertexPrevious::DatumAtVertexElectricField)}) {
            if(pair.first->offset>=0&&pair.second->offset>=0)
              std::memcpy(staged.back().data()+pair.second->offset,staged.back().data()+pair.first->offset,pair.first->length*sizeof(double));
          }
        } else {const auto status=PutState(staged.back(),sample.value,true);if(!status.ok())throw std::runtime_error(status.message);}
      }
      const auto status=PutState(staged.back(),sample.value,false);if(!status.ok())throw std::runtime_error(status.message);
    }
    const auto connected=PrepareLineView(epoch,line.GeometryGeneration,points,&view.value);
    if(!connected.ok())throw std::runtime_error(connected.status.message);
    preparedView=connected.value;
    const Background::BackgroundSnapshot snapshot(
        Background::Provider::ReducedShock,Background::Ownership::ModelOwned,epoch,epoch,until,generation,
        Background::CurrentConfigurationFingerprint(),"core reduced Parker material initial state",
        current?(current->epoch_seconds()==epoch?current->previous_physical_epoch_seconds():current->current_physical_epoch_seconds()):epoch,epoch);
    // Reuse validates the actual retained store object. Validating only the
    // hypothetical replacement above would let a foreign fingerprint or
    // ownership survive when no publication is performed. Reject before any
    // native datum/cache/view changes; no implicit reconfiguration is inferred.
    const auto compatible=preparedView->material->BindSnapshot(reuseMetadata?*current:snapshot);
    if(!compatible.ok())throw std::runtime_error(compatible.message);
    if(!reuseMetadata)publication=store.PreparePublication(snapshot);
  } catch(const std::exception& error){local=Error(error.what());}
  // All callbacks, queries and allocations finish before collective readiness.
  // A rank-local late failure cannot leave some native arrays or metadata at
  // the next epoch. Geometry/particle writes remain forbidden during reads.
  const auto ready=Agree(local);if(!ready.ok())return ready;
  const int count=static_cast<int>(staged.size());int minimum=0,maximum=0;
  MPI_Allreduce(&count,&minimum,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&count,&maximum,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(minimum!=maximum)return Error("reduced initial vertex counts differ across MPI ranks");
  // Compare physical stored primitives rather than sample buffers after
  // publication. Sampling numerators are local and deliberately not compared.
  for(std::size_t i=0;i<staged.size();++i) {
    // Up to 27 SI doubles: position and current/previous B,U,np,p,T,E.
    // Optional E contributes only when allocated. No heap allocation follows
    // readiness, and rank-local sampling buffers are never MPI-summed here.
    double values[27]={},minima[27],maxima[27];int count=3;const auto* bytes=staged[i].data();
    std::copy(preparedView->vertices[i].begin(),preparedView->vertices[i].end(),values);
    for(const auto* datum:{&FL::DatumAtVertexMagneticField,&FL::DatumAtVertexPlasmaVelocity,
        &FL::DatumAtVertexPlasmaDensity,&FL::DatumAtVertexPlasmaPressure,&FL::DatumAtVertexPlasmaTemperature,
        &FL::DatumAtVertexElectricField,&FL::DatumAtVertexPrevious::DatumAtVertexMagneticField,
        &FL::DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,
        &FL::DatumAtVertexPrevious::DatumAtVertexPlasmaPressure,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaTemperature,
        &FL::DatumAtVertexPrevious::DatumAtVertexElectricField})if(datum->offset>=0) {
      std::memcpy(values+count,bytes+datum->offset,datum->length*sizeof(double));count+=datum->length;
    }
    MPI_Allreduce(values,minima,count,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(values,maxima,count,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    for(int j=0;j<count;++j)if(minima[j]!=maxima[j])return Error("reduced initial core states differ across MPI ranks");
  }
  const auto cacheCommitted=CommitCadencePublication(&cadencePublication);
  if(!cacheCommitted.ok())return cacheCommitted;
  auto* vertex=FL::FieldLinesAll[0].GetFirstVertex();
  for(const auto& data:staged){std::memcpy(vertex->GetAssociatedDataBufferPointer(),data.data(),data.size());vertex=vertex->GetNext();}
  if(!reuseMetadata)store.CommitPublication(&publication);
  lineView.swap(preparedView);
  const auto updated=UpdateShockDistances(false);if(!updated.ok())return updated;
  if(PIC::ThisThread==0) {
    const auto front=CurrentFrontSummary();auto* endpoint=FL::FieldLinesAll[0].GetLastVertex();
    double b[3],u[3],n=0,t=0,p=0;endpoint->GetMagneticField(b);endpoint->GetPlasmaVelocity(u);
    endpoint->GetPlasmaDensity(n);endpoint->GetPlasmaTemperature(t);endpoint->GetPlasmaPressure(p);
    std::cout<<"REDUCED_BACKGROUND epoch_s="<<epoch<<" generation="<<Metadata().generation<<" vertices="<<staged.size()
        <<" apex_radius_m="<<front.apexRadiusM<<" accepted_area_m2="<<front.acceptedShockAreaM2
        <<" apex_shock_accepted="<<front.apexShockAccepted<<" endpoint_b_t="<<b[0]<<','<<b[1]<<','<<b[2]
        <<" endpoint_u_m_s="<<u[0]<<','<<u[1]<<','<<u[2]<<" endpoint_n_m3="<<n<<" endpoint_t_k="<<t
        <<" endpoint_p_pa="<<p<<" event="<<EventIdentity()<<std::endl;
  }
  return Core::Status::Success();
}
Core::Status RefreshConnections(double epoch) {
  if(!Enabled())return Core::Status::Success();
  if(FL::nFieldLine!=1)return Error("reduced connectivity requires its registered native line");
  Background::SnapshotStore::Instance().AssertProviderMayWrite(Background::Provider::ReducedShock);
  auto& line=FL::FieldLinesAll[0];
  if(line.GeometryEpochS!=epoch)return Error("connectivity epoch differs from native material geometry");
  std::vector<Position> points;
  for(auto* v=line.GetFirstVertex();v;v=v->GetNext()){const double* x=v->GetX();points.push_back({{x[0],x[1],x[2]}});}
  const auto prepared=PrepareLineView(epoch,line.GeometryGeneration,points);
  if(!prepared.ok())return prepared.status;
  lineView=prepared.value;return UpdateShockDistances(false);
}
Core::Result<std::vector<Connection>> ConnectionsAt(double epoch) {
  using Return=Core::Result<std::vector<Connection>>;
  if(!lineView||FL::nFieldLine!=1||lineView->geometryGeneration!=FL::FieldLinesAll[0].GeometryGeneration)
    return Return::Failure(Core::StatusCode::InvalidState,"reduced material connectivity is stale or unprepared");
  if(epoch==lineView->interval.DepartureS())return Return::Success(lineView->departureConnections);
  return ConnectMaterialLine(lineView->interval,lineView->vertices,epoch,
      Controls(0,lineView->arc[1]).positionToleranceM,
      Controls(lineView->interval.ArrivalS()-lineView->interval.DepartureS(),lineView->arc[1]));
}
Core::Result<std::shared_ptr<const MaterialView>> CaptureMaterialView() {
  using R=Core::Result<std::shared_ptr<const MaterialView>>;
  if(!lineView||FL::nFieldLine!=1||lineView->geometryGeneration!=FL::FieldLinesAll[0].GeometryGeneration)
    return R::Failure(Core::StatusCode::InvalidState,"reduced source/boundary requested stale material geometry");
  return R::Success(lineView);
}
Core::Result<std::shared_ptr<const MaterialView>> PrepareCheckpointView(
    const IntervalView& interval,std::uint64_t generation,const std::vector<Position>& vertices) {
  return PrepareLineView(interval.DepartureS(),generation,vertices,&interval);
}
void CommitCheckpointView(std::shared_ptr<const MaterialView>* prepared) noexcept {
  lineView.swap(*prepared);
}
Core::Status UpdateShockDistances(bool report) {
  if(!lineView||FL::nFieldLine!=1||lineView->geometryGeneration!=FL::FieldLinesAll[0].GeometryGeneration)
    return Error("cannot publish derived shock distances from stale material connectivity");
  auto& line=FL::FieldLinesAll[0];bool accepted=false;
  for(const auto& hit:lineView->departureConnections)accepted=accepted||hit.front.acceptedShock;
  for(std::size_t i=0;i<lineView->vertices.size();++i) {
    double distance=0;bool found=false;
    for(const auto& hit:lineView->departureConnections)if(hit.front.acceptedShock) {
      const double candidate=lineView->arc[i]-hit.arcLengthM;
      if(!found||std::abs(candidate)<std::abs(distance)){found=true;distance=candidate;}
    }
    auto* vertex=line.GetVertex(static_cast<int>(i));
    // Zero distance with an absent flag is NO MEASUREMENT, not an invented
    // shock at the line start. The reduced MFP limiter uses the full root list
    // below, so switching the nearest root cannot manufacture a zero crossing.
    if(FL::DatumAtVertexShockLocation.is_active())vertex->SetDatum(FL::DatumAtVertexShockLocation,accepted?1.:0.);
    if(FL::DatumAtVertexShockLocationDistance.is_active())vertex->SetDatum(FL::DatumAtVertexShockLocationDistance,distance/_AU_);
  }
  if(report&&PIC::ThisThread==0)std::cout<<"REDUCED_CONNECTIONS epoch_s="<<lineView->interval.DepartureS()
      <<" geometry_generation="<<lineView->geometryGeneration<<" supported_intersections="<<lineView->departureConnections.size()<<std::endl;
  return Core::Status::Success();
}
Core::Result<double> CrossingTimeEstimate(double coordinate,double streaming) {
  using Return=Core::Result<double>;
  if(!lineView||FL::nFieldLine!=1||lineView->geometryGeneration!=FL::FieldLinesAll[0].GeometryGeneration||
      !std::isfinite(coordinate)||!std::isfinite(streaming)||coordinate<0||coordinate>=lineView->vertices.size()-1)
    return Return::Failure(Core::StatusCode::InvalidState,"invalid reduced crossing-estimate context");
  const auto segment=static_cast<std::size_t>(coordinate);const double fraction=coordinate-segment;
  const double particleArc=lineView->arc[segment]+fraction*(lineView->arc[segment+1]-lineView->arc[segment]);
  double limit=std::numeric_limits<double>::infinity();
  for(const auto& hit:lineView->departureConnections)if(hit.front.acceptedShock) {
    const auto chord=FL::Geometry::Difference(lineView->vertices[hit.segment+1],lineView->vertices[hit.segment]);
    const double length=FL::Geometry::Norm(chord),projection=FL::Geometry::Dot(chord,hit.front.outwardNormal)/length;
    if(hit.tangent||projection==0)return Return::Failure(Core::StatusCode::UnsupportedCapability,
        "tangent material shock connection requires a resolved space-time crossing, not an obliquity floor");
    const double frontRelative=(hit.front.normalSpeedMPerS-
        FL::Geometry::Dot(hit.front.upstream.velocityMPerS,hit.front.outwardNormal))/projection;
    const double relative=streaming-frontRelative,distance=particleArc-hit.arcLengthM;
    if(std::abs(distance)<=1e-9*length)continue; // same on-front numerical resolution as the maintained MFP shell
    if(distance*relative<0)limit=std::min(limit,std::abs(distance/relative));
  }
  return Return::Success(limit);
}
Core::Status Advance(double departure,double arrival) {
  if(!Enabled())return Core::Status::Success();
  const auto binding=AgreeRealization();if(!binding.ok())return binding;
  FL::cFieldLine::GeometryPlan plan;
  Background::SnapshotStore::Publication publication;
  std::shared_ptr<const LineView> nextLineView;
  Core::Status local=Core::Status::Success();std::string archive,temporary;
  try {
    if(!Initialization::HasActive()||FL::nFieldLine!=1||arrival<=departure)
      local=Error("reduced material advance requires a configured open line and positive full interval");
    if(local.ok()&&Run::Active().get().turbulence.source!=Turbulence::Source::Prescribed)
      local=Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
          "moving local-wave state requires a qualified wave remap; this reduced stage supports prescribed background only");
    const auto view=CapturePhysicalStep(departure,arrival,EventIdentity());
    if(local.ok()&&!view.ok())local=view.status;
    if(local.ok()) {
      auto& line=FL::FieldLinesAll[0];
      if(line.GeometryGeneration>std::numeric_limits<std::uint64_t>::max()-2)
        throw std::overflow_error("reduced snapshot geometry generation exhausted");
      const auto controls=Controls(arrival-departure,line.GetFirstSegment()->GetLength());
      const auto material=MaterialAdapter::Capture(view.value,line.GeometryGeneration,
          std::max(view.value.Generation(),line.GeometryGeneration+1),Background::CurrentConfigurationFingerprint(),controls);
      if(!material.ok())throw std::runtime_error(material.status.message);
      const auto complete=material.value->RequireStep(departure,arrival,
          PIC::ParticleWeightTimeStep::GlobalTimeStep[0],Background::Material::Region::ParkerExterior);
      if(!complete.ok())throw std::runtime_error(complete.message);
      FL::Geometry::Context context;context.departureS=departure;context.arrivalS=arrival;
      context.minimumRadiusM=Initialization::Active().innerRadiusM;context.maximumRadiusM=Initialization::Active().outerRadiusM;
      context.toleranceM=controls.positionToleranceM;context.maximumInsertions=100000;context.expectedGeneration=line.GeometryGeneration;
      const auto prepared=line.PrepareGeometryStep(context,
          [&](std::size_t,const Position& x,double from,double to,Position* target) {
            const auto moved=MoveVertex(view.value,x,from,to,controls);
            if(!moved.ok())return FL::Geometry::Status::Error(moved.status.message);
            *target=moved.value;return FL::Geometry::Status{};
          },[&](const Position& first,const Position& chord,double distance,double epoch,Position* target) {
            const auto traced=TraceMagnetic(view.value,first,chord,distance,epoch,controls);
            if(!traced.ok())return FL::Geometry::Status::Error(traced.status.message);
            *target=traced.value;return FL::Geometry::Status{};
          },[&](int oldIndex,const Position& x,double epoch,std::vector<char>* data) {
            data->assign(FL::cFieldLineVertex::StoredDataBytes(),0);
            Position oldPosition;
            if(oldIndex>=0) {
              auto* vertex=line.GetVertex(oldIndex);
              std::memcpy(data->data(),vertex->GetAssociatedDataBufferPointer(),data->size());
              const double* old=vertex->GetX();oldPosition={{old[0],old[1],old[2]}};
            } else {
              // Newly entering computational material already has an actual
              // plasma history in this analytical provider. Evaluate it on the
              // inverse flow map; never copy a neighbouring vertex's past or
              // pretend that new simulation particles existed before entry.
              const auto inverse=MoveVertex(view.value,x,arrival,departure,controls);
              if(!inverse.ok())return FL::Geometry::Status::Error(inverse.status.message);
              oldPosition=inverse.value;
            }
            const auto current=view.value.AmbientAt(x,epoch),previous=view.value.AmbientAt(oldPosition,departure);
            if(!current.ok()||!previous.ok())return FL::Geometry::Status::Error(!current.ok()?current.status.message:previous.status.message);
            auto status=PutState(*data,current.value,false);if(status.ok())status=PutState(*data,previous.value,true);
            return status.ok()?FL::Geometry::Status{}:FL::Geometry::Status::Error(status.message);
          },&plan,[&](int oldIndex,const Position&,const Position&,double,std::vector<char>* data) {
            // Prescribed/background-only mode evolves no local wave state.
            // Preserve old unowned edge diagnostics; new inactive storage is
            // empty, and is not a declaration of a physical turbulence model.
            data->assign(FL::cFieldLineSegment::StoredDataBytes(),0);
            if(oldIndex>=0)std::memcpy(data->data(),line.GetSegment(oldIndex)->GetAssociatedDataBufferPointer(),data->size());
            return FL::Geometry::Status{};
          });
      if(!prepared.ok())local=Error(prepared.message);
      if(local.ok()) {
        const auto& positions=plan.Positions();
        std::vector<Position> target(positions.insertedPositions.rbegin(),positions.insertedPositions.rend());
        target.insert(target.end(),positions.movedPositions.begin(),positions.movedPositions.end());
        target.resize(target.size()-positions.removedBack);
        const auto connected=PrepareLineView(arrival,line.GeometryGeneration+1,target);
        if(!connected.ok())local=connected.status;else nextLineView=connected.value;
      }
      if(local.ok()) {
        // Native backing state and metadata advance together. Allocate and
        // validate the replacement before collective readiness, so a later
        // bad_alloc cannot leave old metadata attached to new geometry.
        auto& store=Background::SnapshotStore::Instance();
        store.AssertProviderMayWrite(Background::Provider::ReducedShock);
        const Background::BackgroundSnapshot snapshot(
            Background::Provider::ReducedShock,Background::Ownership::ModelOwned,
            arrival,arrival,std::min(EventEndS(),arrival+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]),
            std::max(Metadata().generation,line.GeometryGeneration+2),
            Background::CurrentConfigurationFingerprint(),"core reduced material geometry and state",
            departure,arrival);
        const auto compatible=nextLineView->material->BindSnapshot(snapshot);
        if(!compatible.ok())throw std::runtime_error(compatible.message);
        publication=store.PreparePublication(snapshot);
      }
      // Stage rank-local retirement history before readiness. Histories differ
      // across particle owners; preserve them separately instead of summing
      // exposure as if each MPI rank represented extra elapsed physical time.
      if(local.ok()&&!plan.RetiredSamples().empty()) {
        archive=std::string(PIC::OutputDataFileDirectory)+"/reduced-material-retired-"+
            std::to_string(line.GeometryGeneration+1)+"-rank"+std::to_string(PIC::ThisThread)+".dat";
        temporary=archive+".tmp";std::ofstream out(temporary);
        out<<std::setprecision(17)<<"# reduced-material-retirement-v1 local numerators; exposure is not MPI-additive\n";
        out<<"# event "<<EventIdentity()<<" departure_s "<<departure<<" arrival_s "<<arrival<<'\n';
        for(const auto& record:plan.RetiredSamples()) {
          out<<"segment "<<record.label<<" counts "<<record.exposure.collecting<<' '<<record.exposure.completed
             <<" seconds "<<record.exposure.collectingSeconds<<' '<<record.exposure.completedSeconds<<" fixed_dt_s "<<record.exposure.fixedStepSeconds<<'\n';
          for(const auto& points:{record.departureEndpoints,record.targetEndpoints}){out<<"endpoints";for(const auto& point:points)for(double value:point)out<<' '<<value;out<<'\n';}
          for(const auto& values:{record.collecting,record.completed}){out<<"numerator";for(double value:values)out<<' '<<value;out<<'\n';}
          out<<"opaque_edge_hex ";for(unsigned char value:record.associatedData)out<<std::hex<<std::setw(2)<<std::setfill('0')<<static_cast<unsigned>(value);out<<std::dec<<'\n';
        }
        out.flush();if(!out.good())local=Error("cannot stage reduced material retirement history");
      }
    }
  }catch(const std::exception& error){local=Error(error.what());}
  const auto ready=Agree(local);
  if(!ready.ok()){if(!temporary.empty())std::remove(temporary.c_str());return ready;}
  // Compare the whole prepared coordinate realization, not only its endpoint.
  // All workers/migration/read phases must already have joined at this call.
  const auto& positions=plan.Positions();
  int sizes[3]={static_cast<int>(positions.movedPositions.size()),static_cast<int>(positions.insertedPositions.size()),static_cast<int>(positions.removedBack)},lo[3],hi[3];
  MPI_Allreduce(sizes,lo,3,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);MPI_Allreduce(sizes,hi,3,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  for(int i=0;i<3;++i)if(lo[i]!=hi[i]){if(!temporary.empty())std::remove(temporary.c_str());return Error("MPI reduced material candidate topology differs");}
  for(const auto* table:{&positions.movedPositions,&positions.insertedPositions})for(const auto& point:*table) {
    double minimum[3],maximum[3];MPI_Allreduce(point.data(),minimum,3,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);MPI_Allreduce(point.data(),maximum,3,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    for(int d=0;d<3;++d)if(minimum[d]!=maximum[d]){if(!temporary.empty())std::remove(temporary.c_str());return Error("MPI reduced material candidate positions differ");}
  }
  const auto committed=FL::FieldLinesAll[0].CommitGeometryStep(&plan);
  if(!committed.ok())return Error(committed.message);
  Background::SnapshotStore::Instance().CommitPublication(&publication);
  lineView.swap(nextLineView);
  const auto distances=UpdateShockDistances(false);if(!distances.ok())return distances;
  if(!temporary.empty()&&std::rename(temporary.c_str(),archive.c_str())!=0)return Error("cannot commit reduced material retirement history file");
  for(const auto& loss:plan.Losses())if(loss.macros)
    std::cout<<"REDUCED_MATERIAL_LOSS rank="<<PIC::ThisThread<<" species="<<loss.species<<" macros="<<loss.macros<<" represented="<<loss.represented<<" plasma_frame_kinetic_j="<<loss.kineticJ<<'\n';
  if(PIC::ThisThread==0)std::cout<<"REDUCED_GEOMETRY epoch_s="<<arrival<<" generation="<<FL::FieldLinesAll[0].GeometryGeneration
      <<" inserted="<<positions.insertedPositions.size()<<" removed="<<positions.removedBack<<" segments="<<FL::FieldLinesAll[0].GetTotalSegmentNumber()<<std::endl;
  return Core::Status::Success();
}
}}}
