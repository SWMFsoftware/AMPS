//$Id$

//implementation of methods in namespace FieldLine

#include "pic.h"
#include <cstdio>

//static variables of class cFieldLineVertex
int PIC::FieldLine::cFieldLineVertex::totalAssociatedDataLength=-1;
int PIC::FieldLine::cFieldLineVertex::sampleDataLength=-1;
int PIC::FieldLine::cFieldLineVertex::CollectingSamplingOffset=-1;
int PIC::FieldLine::cFieldLineVertex::CompletedSamplingOffset=-1;

int PIC::FieldLine::cFieldLineSegment::totalAssociatedDataLength=-1;
int PIC::FieldLine::cFieldLineSegment::sampleDataLength=-1;
int PIC::FieldLine::cFieldLineSegment::CollectingSamplingOffset=-1;
int PIC::FieldLine::cFieldLineSegment::CompletedSamplingOffset=-1;
vector<PIC::Datum::cDatumStored*> PIC::FieldLine::cFieldLineSegment::DataStoredAtSegment;
vector<PIC::Datum::cDatumSampled*> PIC::FieldLine::cFieldLineSegment::DataSampledAtSegment;
PIC::Datum::cDatumSampled* PIC::FieldLine::cFieldLineSegment::SamplingWeightDatum=nullptr;

//hooks for calculating the magnertic tube radius and the volume of the field line segment
PIC::FieldLine::fMagneticTubeRadius PIC::FieldLine::MagneticTubeRadius=NULL;
PIC::FieldLine::fSegmentVolume PIC::FieldLine::SegmentVolume=NULL;

//process data associated with a filed line before it is saved in a file
PIC::FieldLine::fUserDefinedfDataProcessingManager PIC::FieldLine::UserDefinedfDataProcessingManager=NULL;

//the following is used to output the distance from the beginning of the 
//field line in units other than SI
//first -> the conversion factor
//second -> the string contains the unit symbol
std::pair<double,string> PIC::FieldLine::cFieldLine::OutputLengthConversionFactor(1.0,"");

//user-defined function that defiens an title that is printed in the Tecplot output file (e.g., simulation time of the file)
PIC::FieldLine::fUserDefinedTecplotFileTitle PIC::FieldLine::UserDefinedTecplotFileTitle=NULL; 

//sample cycle counter 
int PIC::FieldLine::SampleCycleCounter=0;
long int PIC::FieldLine::LastDataOutputFileNumber=-1;

//offset of the field line related data in a particle state vector
int PIC::FieldLine::ParticleDataOffset=-1; 

      void PIC::FieldLine::cFieldLineSegment::DeleteAttachedParticles() {
        if (_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_) {
          long int ptr_next,ptr=FirstParticleIndex;

          while (ptr!=-1) {
            ptr_next=PIC::ParticleBuffer::GetNext(ptr);
            PIC::ParticleBuffer::DeleteParticle(ptr);
            ptr=ptr_next;
          }

          FirstParticleIndex=-1;
        }
      }


namespace PIC {
  namespace FieldLine {

// Preparation owns detached records until commit transfers them into the
// line. Shared ownership makes copied/moved plan handles safe to discard;
// stale or rejected plans return every allocated record without touching live
// topology. The prepared MPI rows are likewise allocation-complete here.
struct cFieldLine::GeometryResources {
  std::vector<cFieldLineVertex*> vertices;
  std::vector<cFieldLineSegment*> segments;
  std::vector<cFieldLineSegment*> index;
  std::function<void()> publishOwnership;
  std::uint64_t ownershipGeneration=0;
  bool adopted=false;
  ~GeometryResources() {
    if(adopted)return;
    for(auto* v:vertices)VerticesAll.deleteElement(v);
    for(auto* edge:segments)SegmentsAll.deleteElement(edge);
  }
};

Geometry::Status cFieldLine::RegisterGeometryOwner(const std::string& identity,
    double epochS,std::int64_t firstLabel) {
#if _CUDA_MODE_ == _ON_
  return Geometry::Status::Error("material geometry requires constructed CPU segment storage; CUDA owner is not implemented");
#endif
  if(identity.empty()||!std::isfinite(epochS)||ExternallyManagedGeometry||nSegment<1||is_loop())
    return Geometry::Status::Error("invalid or duplicate material geometry ownership registration");
  if(firstLabel>std::numeric_limits<std::int64_t>::max()-nSegment)
    return Geometry::Status::Error("material node labels exceed the signed identity range");
  for(int line=0;line<nFieldLine;++line)
    if(&FieldLinesAll[line]!=this&&FieldLinesAll[line].ExternallyManagedGeometry&&
        FieldLinesAll[line].GetGeometryOwnerIdentity()==identity)
      return Geometry::Status::Error("duplicate physical material line identity");
  for(auto* edge=FirstSegment;edge;edge=edge->GetNext())
    if(edge->SamplingExposure.collecting||edge->SamplingExposure.completed)
      return Geometry::Status::Error("cannot reinterpret an existing fixed-grid sample history as material exposure");
  GeometryOwnerIdentity=identity;GeometryEpochS=epochS;FirstMaterialLabel=firstLabel;
  ExternallyManagedGeometry=true;return Geometry::Status{};
}

Geometry::Status cFieldLine::PrepareGeometryStep(const Geometry::Context& context,
    const Geometry::Motion& motion,const Geometry::InnerTrace& trace,
    const VertexStateInitializer& initialize,GeometryPlan* output,
    const SegmentStateInitializer& initializeSegment) {
  using Geometry::Status;
  if(!output||!initialize||!ExternallyManagedGeometry||is_loop()||nSegment<1||
      context.expectedGeneration!=GeometryGeneration||
      (ExternallyManagedGeometry&&context.departureS!=GeometryEpochS))
    return Status::Error("material geometry owner/epoch/generation/open-line precondition failed");
  if(_PIC_PARTICLE_LIST_ATTACHING_!=_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_)
    return Status::Error("material geometry requires explicit segment particle attachment");
  GeometryPlan candidate;candidate.owner=this;
  std::vector<Geometry::Point> old;
  for(int i=0;i<=nSegment;++i) {
    const double* x=GetVertex(i)->GetX();old.push_back({{x[0],x[1],x[2]}});
  }
  const auto prepared=Geometry::Prepare(old,context,motion,trace,&candidate.positions);
  if(!prepared.ok())return prepared;
  const std::size_t prefix=candidate.positions.insertedPositions.size();
  if(prefix>static_cast<std::size_t>(std::numeric_limits<int>::max()-nSegment)||
      prefix>(static_cast<std::uint64_t>(FirstMaterialLabel)-
          static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::min()))||
      GeometryGeneration==std::numeric_limits<std::uint64_t>::max())
    return Status::Error("material topology label/index/generation overflow");
  try {
    // Stage complete physical datum prefixes. New records explicitly request
    // index -1; the provider must supply a valid entry state, not a guessed
    // previous history. Existing sampled buffers are never copied by this API.
    const int stored=cFieldLineVertex::StoredDataBytes();
    if(stored<0)return Status::Error("material vertex datum layout is uninitialized");
    candidate.storedVertexData.resize(old.size()+prefix);
    const std::size_t keptVertices=old.size()+prefix-candidate.positions.removedBack;
    for(std::size_t i=0;i<old.size()+prefix;++i) {
      const bool inserted=i>=old.size();
      const std::size_t finalIndex=inserted?prefix-1-(i-old.size()):prefix+i;
      // Do not query unsupported target-epoch primitives for records which
      // are being retired. Their departure data are retained in the ledger.
      if(finalIndex>=keptVertices)continue;
      const auto& x=inserted?candidate.positions.insertedPositions[i-old.size()]:candidate.positions.movedPositions[i];
      const auto status=initialize(inserted?-1:static_cast<int>(i),x,context.arrivalS,&candidate.storedVertexData[i]);
      if(!status.ok())return status;
      if(candidate.storedVertexData[i].size()!=static_cast<std::size_t>(stored))
        return Status::Error("material vertex initializer returned the wrong physical datum layout");
      // Stored field-line datums are double arrays. Invalid optional data must
      // be rejected by the owner as well; no NaN replacement is allowed here.
      for(int at=0;at+static_cast<int>(sizeof(double))<=stored;at+=sizeof(double)) {
        double value;std::memcpy(&value,candidate.storedVertexData[i].data()+at,sizeof(value));
        if(!std::isfinite(value))return Status::Error("nonfinite staged material vertex datum");
      }
    }
    const int edgeBytes=cFieldLineSegment::StoredDataBytes();
    if(edgeBytes<0)return Status::Error("material edge datum layout is uninitialized");
    if(edgeBytes>0&&!initializeSegment)
      return Status::Error("material geometry requires an explicit initializer for stored segment physics");
    candidate.storedSegmentData.resize(keptVertices-1);
    for(std::size_t i=0;i+1<keptVertices;++i) {
      if(!initializeSegment)continue; // no stored edge physics is allocated
      const auto point=[&](std::size_t index)->const Geometry::Point& {
        return index<prefix?candidate.positions.insertedPositions[prefix-1-index]:
            candidate.positions.movedPositions[index-prefix];
      };
      const int oldIndex=i<prefix?-1:static_cast<int>(i-prefix);
      const auto status=initializeSegment(oldIndex,point(i),point(i+1),context.arrivalS,&candidate.storedSegmentData[i]);
      if(!status.ok())return status;
      if(candidate.storedSegmentData[i].size()!=static_cast<std::size_t>(edgeBytes))
        return Status::Error("material edge initializer returned the wrong physical datum layout");
      for(int at=0;at+static_cast<int>(sizeof(double))<=edgeBytes;at+=sizeof(double)) {
        double value;std::memcpy(&value,candidate.storedSegmentData[i].data()+at,sizeof(value));
        if(!std::isfinite(value))return Status::Error("nonfinite staged material edge datum");
      }
    }
    candidate.losses.resize(PIC::nTotalSpecies);
    for(int spec=0;spec<PIC::nTotalSpecies;++spec)candidate.losses[spec]={spec,0,0,0,0};
    const int retained=nSegment-static_cast<int>(candidate.positions.removedBack);
    for(int segment=0;segment<nSegment;++segment) {
      auto* edge=GetSegment(segment);
      candidate.departureHeads.push_back(edge->FirstParticleIndex);
      candidate.departureExposure.push_back(edge->SamplingExposure);
      if(edge->tempFirstParticleIndex.load()!=-1)
        return Status::Error("material geometry update overlaps pending particle migration");
      if(segment>=retained) {
        GeometryPlan::RetiredSample retired;
        retired.label=FirstMaterialLabel+segment;retired.exposure=edge->SamplingExposure;
        retired.collecting=edge->MaterialCollecting;retired.completed=edge->MaterialCompleted;
        const int bytes=edge->AssociatedDataLength();
        if(bytes>0)retired.associatedData.assign(edge->GetAssociatedDataBufferPointer(),
            edge->GetAssociatedDataBufferPointer()+bytes);
        retired.departureEndpoints={{old[segment],old[segment+1]}};
        retired.targetEndpoints={{candidate.positions.movedPositions[segment],candidate.positions.movedPositions[segment+1]}};
        candidate.retired.push_back(std::move(retired));
      }
      for(long int p=edge->FirstParticleIndex;p!=-1;p=PIC::ParticleBuffer::GetNext(p)) {
        const double s=PIC::ParticleBuffer::GetFieldLineCoord(p);
        if(!std::isfinite(s)||s<segment||s>=segment+1)
          return Status::Error("particle material coordinate disagrees with segment ownership");
        const auto speciesValue=PIC::ParticleBuffer::GetI(p);
        const double individual=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        const double parallel=PIC::ParticleBuffer::GetVParallel(p),normal=PIC::ParticleBuffer::GetVNormal(p);
        if(speciesValue>=static_cast<unsigned int>(PIC::nTotalSpecies)||
            !std::isfinite(individual)||individual<=0||!std::isfinite(parallel)||!std::isfinite(normal)||normal<0)
          return Status::Error("invalid species/weight/velocity in a material particle record");
        candidate.departureParticles.push_back({p,PIC::ParticleBuffer::GetNext(p),s,
            individual,parallel,normal,static_cast<int>(speciesValue)});
        // Head insertion changes the integer part of the stored global line
        // coordinate. Binary64 addition S+prefix can round its fractional part
        // when the index crosses a power of two (e.g. 15.x -> 16.x). X [m] is
        // a DERIVED view of that published coordinate: interpolate using its
        // post-addition fraction, exactly as cFieldLine::GetCartesian does.
        // Using the old fraction left finite one-ULP contradictions in native
        // records, which strict checkpoint/restart correctly rejected.
        // This adds no plasma-relative displacement or solar-wind term: the
        // owner already supplied the target-epoch endpoints once above.
        const double mapped=s+prefix;
        const double mappedSegment=static_cast<double>(segment)+prefix;
        const double fraction=mapped-mappedSegment;
        // Never silently clip a rounded coordinate onto another segment. If
        // the address can no longer represent this material segment, reject
        // the still-private transaction before any geometry/particle mutation.
        if(!std::isfinite(mapped)||fraction<0.0||fraction>=1.0)
          return Status::Error("material coordinate loses segment identity after prefix insertion");
        const auto& a=candidate.positions.movedPositions[segment];
        const auto& b=candidate.positions.movedPositions[segment+1];
        Geometry::Point x;
        for(int d=0;d<3;++d)x[d]=(1-fraction)*a[d]+fraction*b[d];
        if(!Geometry::Finite(x))return Status::Error("invalid staged material particle Cartesian position");
        if(segment<retained) {
          candidate.survivors.push_back({p,mapped,x});
        }else {
          candidate.retiringParticles.push_back({p,mapped,x});
          const int species=static_cast<int>(speciesValue);
          if(species<0||species>=PIC::nTotalSpecies)return Status::Error("invalid species in material boundary loss");
          const double weight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]*
              PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
          auto& loss=candidate.losses[species];
          if(!std::isfinite(weight)||weight<=0||!std::isfinite(loss.represented+weight)||
              loss.macros==std::numeric_limits<std::uint64_t>::max())
            return Status::Error("invalid/overflowing represented material boundary loss");
          ++loss.macros;loss.represented+=weight;
          const double parallel=PIC::ParticleBuffer::GetVParallel(p),normal=PIC::ParticleBuffer::GetVNormal(p);
          const double cLight=299792458.0,mass=PIC::MolecularData::GetMass(species);
          const double beta2=(parallel/cLight)*(parallel/cLight)+(normal/cLight)*(normal/cLight);
          if(!std::isfinite(beta2)||beta2>=1||!std::isfinite(mass)||mass<=0)
            return Status::Error("invalid relativistic particle state in material boundary loss");
          const double gamma=1/std::sqrt(1-beta2);
          // gamma-1 in this form avoids cancellation for low particle speeds.
          loss.kineticJ+=weight*mass*cLight*cLight*beta2*gamma*gamma/(gamma+1);
          loss.parallelMomentumKgMPerS+=weight*gamma*mass*parallel;
          if(!std::isfinite(loss.kineticJ)||!std::isfinite(loss.parallelMomentumKgMPerS))
            return Status::Error("material boundary energy/momentum ledger overflow");
        }
      }
    }
    // Allocate every detached record and ownership row during PREPARATION,
    // before the host's all-rank readiness agreement. Commit below performs no
    // allocation. New records own neither particles nor old sample history.
    candidate.resources=std::make_shared<GeometryResources>();
    auto& resources=*candidate.resources;
    resources.vertices.reserve(prefix);resources.segments.reserve(prefix);
    for(std::size_t i=0;i<prefix;++i) {
      auto* vertex=VerticesAll.newElement();resources.vertices.push_back(vertex);
      vertex->SetPrev(nullptr);vertex->SetNext(nullptr);vertex->cleanDataBuffer();
      auto* edge=SegmentsAll.newElement();resources.segments.push_back(edge);
      edge->SetPrev(nullptr);edge->SetNext(nullptr);edge->FirstParticleIndex=-1;
      edge->tempFirstParticleIndex=-1;edge->ResetMaterialSampling();edge->cleanDataBuffer();
      edge->Thread=FirstSegment->Thread;
      // A recycled injection normalization must not enter a new empty edge.
      for(int species=0;species<PIC::nTotalSpecies;++species)edge->SetWeight(0,species);
    }
    resources.index.assign(resources.segments.rbegin(),resources.segments.rend());
    resources.index.insert(resources.index.end(),SegmentPointerTable.begin(),SegmentPointerTable.end());
    std::vector<cFieldLineSegment*> finalEdges=resources.index;
    finalEdges.resize(finalEdges.size()-candidate.positions.removedBack);
    const auto ownership=PIC::ParallelFieldLines::PrepareThreadSegmentTableReplacement(this,finalEdges,&resources.publishOwnership);
    if(!ownership.ok())throw std::runtime_error(ownership.message);
    resources.ownershipGeneration=PIC::ParallelFieldLines::GetThreadSegmentTableGeneration();
  }catch(const std::exception& e){return Status::Error(std::string("material state preparation failed: ")+e.what());}
  *output=std::move(candidate);return Status{};
}

Geometry::Status cFieldLine::CommitGeometryStep(GeometryPlan* plan) {
  using Geometry::Status;
  if(!plan||plan->owner!=this||plan->consumed||
      plan->positions.context.expectedGeneration!=GeometryGeneration||
      plan->positions.oldPositions.size()!=static_cast<std::size_t>(nSegment+1))
    return Status::Error("stale, foreign or already consumed material geometry plan");
  for(int i=0;i<=nSegment;++i) {
    const double* x=GetVertex(i)->GetX();
    for(int d=0;d<3;++d)if(x[d]!=plan->positions.oldPositions[i][d])
      return Status::Error("material geometry changed after plan preparation");
  }
  // A plan is valid only at the joined host boundary that prepared it. Even
  // unchanged geometry cannot authorize committing stale particle views after
  // injection, migration, weight changes or particle motion in between.
  for(int i=0;i<nSegment;++i) {
    if(GetSegment(i)->FirstParticleIndex!=plan->departureHeads[i]||GetSegment(i)->tempFirstParticleIndex.load()!=-1)
      return Status::Error("material particle ownership changed after preparation");
    const auto& now=GetSegment(i)->SamplingExposure;const auto& old=plan->departureExposure[i];
    if(now.collecting!=old.collecting||now.completed!=old.completed||
        now.collectingSeconds!=old.collectingSeconds||now.completedSeconds!=old.completedSeconds)
      return Status::Error("material sampling changed after geometry preparation");
  }
  for(const auto& p:plan->departureParticles)
    if(PIC::ParticleBuffer::GetNext(p.particle)!=p.next||
        PIC::ParticleBuffer::GetFieldLineCoord(p.particle)!=p.coordinate||
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p.particle)!=p.weight||
        PIC::ParticleBuffer::GetVParallel(p.particle)!=p.parallel||
        PIC::ParticleBuffer::GetVNormal(p.particle)!=p.normal||
        PIC::ParticleBuffer::GetI(p.particle)!=p.species)
      return Status::Error("material particle state changed after preparation");
  const std::size_t prefix=plan->positions.insertedPositions.size();
  if(!plan->resources||plan->resources->adopted||
      plan->resources->ownershipGeneration!=PIC::ParallelFieldLines::GetThreadSegmentTableGeneration())
    return Status::Error("material plan has stale or consumed prepared ownership resources");
  auto& vertices=plan->resources->vertices;auto& segments=plan->resources->segments;
  const int oldVertices=nSegment+1;
  for(int i=0;i<oldVertices;++i) {
    auto* vertex=GetVertex(i);
    vertex->SetX(plan->positions.movedPositions[i].data());
    if(!plan->storedVertexData[i].empty())std::memcpy(vertex->GetAssociatedDataBufferPointer(),
        plan->storedVertexData[i].data(),plan->storedVertexData[i].size());
  }
  for(std::size_t i=0;i<prefix;++i) {
    auto* vertex=vertices[i];auto* segment=segments[i];
    vertex->SetX(plan->positions.insertedPositions[i].data());
    if(!plan->storedVertexData[oldVertices+i].empty())std::memcpy(vertex->GetAssociatedDataBufferPointer(),
        plan->storedVertexData[oldVertices+i].data(),plan->storedVertexData[oldVertices+i].size());
    vertex->SetNext(FirstVertex);FirstVertex->SetPrev(vertex);
    segment->SetVertices(vertex,FirstVertex);segment->SetNext(FirstSegment);
    FirstSegment->SetPrev(segment);FirstVertex=vertex;FirstSegment=segment;
    ++nSegment;
  }
  // Swapping a fully prepared index avoids both allocation during commit and
  // invalidation of the live index's storage during an unsuccessful prepare.
  SegmentPointerTable.swap(plan->resources->index);
  // All new records now belong to the live line. Destruction of any copied
  // plan must not recycle records that have been attached by this commit.
  plan->resources->adopted=true;
  // Publish the moved view even for records about to leave. Trajectory/death
  // diagnostics must see the target-epoch location, not the old Cartesian X.
  // No particle's momentum, statistical weight or stochastic history changes.
  for(auto* views:{&plan->survivors,&plan->retiringParticles})for(auto& p:*views) {
    PIC::ParticleBuffer::SetFieldLineCoord(p.coordinate,p.particle);
    PIC::ParticleBuffer::SetX(p.position.data(),p.particle);
  }
  // Trim the complete candidate. A large displacement may also discard some
  // empty prefix edges; old physical losses/samples remain available in plan.
  CutBack(static_cast<int>(plan->positions.removedBack));
  for(std::size_t i=0;i<plan->storedSegmentData.size();++i)
    if(!plan->storedSegmentData[i].empty())std::memcpy(GetSegment(static_cast<int>(i))->GetAssociatedDataBufferPointer(),
        plan->storedSegmentData[i].data(),plan->storedSegmentData[i].size());
  FirstMaterialLabel-=static_cast<std::int64_t>(prefix);
  UpdateLength();++GeometryGeneration;GeometryEpochS=plan->positions.context.arrivalS;
  plan->resources->publishOwnership();
  ExternallyManagedGeometry=true;plan->consumed=true;
  return Status{};
}

    cDatumStored DatumAtVertexElectricField(3,"\"Ex [V/m]\",\"Ey [V/m]\",\"Ez [V/m]\"",false);
    cDatumStored DatumAtVertexMagneticField(3,"\"Bx [nT]\",\"By [nT]\",\"Bz [nT]\"",true);
    cDatumStored DatumAtVertexPlasmaVelocity(3,"\"Plasma Vx [m/s]\",\"Plasma Vy [m/s]\",\"Plasma Vz [m/s]\"",true);
    cDatumStored DatumAtVertexPlasmaDensity(1,"\"Plasma number density [1/m^3]\"", true);
    cDatumStored DatumAtVertexPlasmaTemperature(1,"\"Plasma Temperature [K]\"", true);
    cDatumStored DatumAtVertexPlasmaPressure(1,"\"Plasma pressure [Pa]\"", true);
    cDatumStored DatumAtVertexMagneticFluxFunction(1,"\"MagneticFluxFunction [nT*m]\"", true);
    cDatumStored DatumAtVertexPlasmaWaves(2,"\"Wave1\",\"Wave2\"",true);
    cDatumStored DatumAtVertexShockLocation(1,"\"Shock Location\"",true);

    cDatumStored DatumAtVertexPrevious::DatumAtVertexElectricField(3,"\"Ex [V/m]\",\"Ey [V/m]\",\"Ez [V/m]\"",false);
    cDatumStored DatumAtVertexPrevious::DatumAtVertexMagneticField(3,"\"Bx [nT]\",\"By [nT]\",\"Bz [nT]\"",false);
    cDatumStored DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity(3,"\"Plasma Vx [m/s]\",\"Plasma Vy [m/s]\",\"Plasma Vz [m/s]\"",false);
    cDatumStored DatumAtVertexPrevious::DatumAtVertexPlasmaDensity(1,"\"Plasma number density [1/m^3]\"", false);
    cDatumStored DatumAtVertexPrevious::DatumAtVertexPlasmaTemperature(1,"\"Plasma Temperature [K]\"", false);
    cDatumStored DatumAtVertexPrevious::DatumAtVertexPlasmaPressure(1,"\"Plasma pressure [Pa]\"", false);
    cDatumStored DatumAtVertexPrevious::DatumAtVertexPlasmaWaves(2,"\"Wave1\",\"Wave2\"",false);

    cDatumStored DatumAtVertexFluence(1,"\"Fluence\"", true); 
    cDatumStored DatumAtVertexShockLocationDistance(1,"\"Distance To Shock [AU]\"",true);

    cDatumTimed DatumAtVertexParticleWeight(1,"\"Particle Weight\"",false);
    cDatumTimed DatumAtVertexParticleNumber(1,"\"Particle Number\"",true);
    cDatumTimed DatumAtVertexNumberDensity(1,"\"Number Density[1/m^3]\"",true);
    cDatumWeighted DatumAtVertexParticleEnergy(1,"\"Kinetic energy [MeV]\"",true);
    cDatumWeighted DatumAtVertexParticleSpeed(1,"\"Speed [m/s]\"",true);

    cDatumTimed DatumAtVertexNumberDensityEnergyBinned("Number Density[1/m^3]",10,0.1,600.0,PIC::Datum::cDatum::DistTypeLog); // units for the energy limits are MeV 

    cDatumWeighted DatumAtVertexParticleCosPitchAngle(1,"\"Cos Pitch Angle\"",true);
    cDatumWeighted DatumAtVertexParticleAbsCosPitchAngle(1,"\"Abs Cos Pitch Angle\"",true);

    cDatumTimed DatumAtVertexNumberDensity_mu_positive(1,"\"Number Density (mu positive) [1/m^3]\"",true); 
    cDatumTimed DatumAtVertexNumberDensity_mu_negative(1,"\"Number Density (mu negative) [1/m^3]\"",true);

    cDatumTimed DatumAtVertexParticleFlux_mu_positive(1,"\"Particle Flux (mu positive) [1/m^3]\"",true);
    cDatumTimed DatumAtVertexParticleFlux_mu_negative(1,"\"Particle Flux (mu negative) [1/m^3]\"",true);

    cDatumWeighted DatumAtGridParticleEnergy(1,"\"Kinetic energy [J]\"",true);

    cVertexAllocationManager VertexAllocationManager;

    vector<cDatumStored*> DataStoredAtVertex;
    vector<cDatumSampled*> DataSampledAtVertex;

    cAssociatedDataAMRstack<cFieldLineVertex>  VerticesAll;
    cAssociatedDataAMRstack<cFieldLineSegment> SegmentsAll;
    cFieldLine *FieldLinesAll = NULL;
    
    long int nFieldLine=0;

    double TimeLastUpdate = -1;

    //=========================================================================
    bool cFieldLine::is_broken() {
      int count;
      cFieldLineSegment* Segment;
      //  cFieldLineVertex*  Vertex; 
      
      //check connectivity forward
      count = 1, Segment = FirstSegment;//, Vertex  = FirstVertex;

      for (int iSegment=0; iSegment<nSegment; iSegment++) {
        if (Segment->GetNext()==NULL) break; //|| Vertex->GetNext()==NULL) break;

        Segment = Segment->GetNext();
        //    Vertex  = Vertex->GetNext();

        count++;
      }
      
      if ((count<nSegment) || (Segment != LastSegment)) return true; // || Vertex != LastVertex)

      //check connectivity backward
      count = 1, Segment = LastSegment;//, Vertex  = LastVertex;

      for (int iSegment=0; iSegment<nSegment; iSegment++) {
	      if (Segment->GetPrev()==NULL) break; // || Vertex->GetPrev()==NULL) break;

	      Segment = Segment->GetPrev();
	      //    Vertex  = Vertex->GetPrev();

	      count++;
      }
      
      if ((count<nSegment) || (Segment != FirstSegment)) return true; // || Vertex != FirstVertex)

      //the line is fully connected
      return false;
    }
    
    //delete the filed line
    void  cFieldLine::Delete() {
      //delete vetexes
      cFieldLineVertex *v,*v_next; 
      
      for (v=FirstVertex;v!=NULL;v=v_next) {
        v_next=v->GetNext();
        VerticesAll.deleteElement(v);
      } 

      FirstVertex=NULL,LastVertex=NULL;

      //delete segments 
      cFieldLineSegment *s,*s_next;

      for (s=FirstSegment;s!=NULL;s=s_next) {
        s_next=s->GetNext();

        //delete particles that attached to the segment 
        s->DeleteAttachedParticles();

        //delete the segment
        SegmentsAll.deleteElement(s);
      }  

      nSegment=0; 
      FirstSegment=NULL,LastSegment=NULL;
      TotalLength=0.0,IsSet=0;
      SegmentPointerTable.clear();
      // Deleting/replacing a line invalidates outstanding geometry plans even
      // if an identical coordinate list is subsequently allocated in its slot.
      if(GeometryGeneration==std::numeric_limits<std::uint64_t>::max())
        exit(__LINE__,__FILE__,"field-line generation overflow during retirement");
      ++GeometryGeneration;ExternallyManagedGeometry=false;GeometryEpochS=0;FirstMaterialLabel=0;
      GeometryOwnerIdentity.clear();
    }

    //=========================================================================
    cFieldLineVertex* cFieldLine::AddFront(double *xIn) {
      if (nSegment==0) {
        return Add(xIn); 
      }
      else {
        cFieldLineVertex* newVertex;
        cFieldLineSegment* newSegment;

        //allocate the new vertex
        newVertex = VerticesAll.newElement();
        newVertex->SetX(xIn);

        //link the list of vertexes
        newVertex->SetNext(FirstVertex);
        newVertex->SetPrev(NULL);
        FirstVertex->SetPrev(newVertex);

        //allocate a new segement  
        newSegment=SegmentsAll.newElement();
        newSegment->ResetMaterialSampling();
        newSegment->SetVertices(newVertex,FirstVertex);

        newSegment->SetPrev(NULL);
        newSegment->SetNext(FirstSegment);

        FirstSegment->SetPrev(newSegment);

        FirstVertex=newVertex;
        FirstSegment=newSegment; 

        SegmentPointerTable.insert(SegmentPointerTable.begin(),FirstSegment);        

        nSegment++; 
      }

      return FirstVertex;
    } 

    cFieldLineVertex* cFieldLine::AddBack(double *xIn) {
       return Add(xIn);
    }

    // VerticesAll and SegmentsAll are stacks that hand out released elements as
    // they are: a vertex or segment released by Delete(), CutBack() or CutFront()
    // keeps its old prev/next links.  Add() therefore sets every link of the
    // element it allocates.  It used to leave 'next' of the new last vertex and
    // segment (and both links of the first ones) as they were, so a field line
    // rebuilt from released elements -- e.g. the import in amps_get_bline_c_()
    // after DeleteAll() -- ended in a stale link: the last segment pointed back
    // to the previous one, and every loop 'for (s=FirstSegment;s!=NULL;
    // s=s->GetNext())' (UpdateLength(), Delete(), ...) never ended.
    cFieldLineVertex* cFieldLine::Add(double *xIn) {

      if (PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal==false) {
        exit(__LINE__,__FILE__,"Error: when field line module is used, PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal has to be set to TRUE");
      }

      // check if field lineis unset
      if (IsSet == 0) {
        if (FirstVertex==NULL) {
          //allocate the first vertex
          LastVertex = (FirstVertex = VerticesAll.newElement());
          FirstVertex->SetX(xIn);
          FirstVertex->SetPrev(NULL);
          FirstVertex->SetNext(NULL);

          //the first segment can't be allocated yet
          nSegment = 0; TotalLength = 0.0; IsSet = 0;
        }
        else{
         //allocate the second vertex
         LastVertex = VerticesAll.newElement();
         LastVertex->SetX(xIn);
         LastVertex->SetPrev(FirstVertex);
         LastVertex->SetNext(NULL);
         FirstVertex->SetNext(LastVertex);

         //allocate the first segment
         LastSegment = (FirstSegment = SegmentsAll.newElement());
         LastSegment->ResetMaterialSampling();
         LastSegment->SetVertices(LastVertex->GetPrev(), LastVertex);
         LastSegment->SetPrev(NULL);
         LastSegment->SetNext(NULL);
         nSegment++; TotalLength+= LastSegment->GetLength(); IsSet = 1;

         SegmentPointerTable.push_back(FirstSegment);
       }

       return LastVertex;
     }
      
     //allocate vertex
     cFieldLineVertex* newVertex;
     newVertex = VerticesAll.newElement();
     newVertex->SetX(xIn);
      
     //connect it to the last vertex in the field line
     LastVertex->SetNext(newVertex);
     newVertex->SetPrev(LastVertex);
     newVertex->SetNext(NULL);
      
     //now the new vertex is the last one
     LastVertex = newVertex;
      
     //allocate segment
     cFieldLineSegment* newSegment;
     newSegment = SegmentsAll.newElement();
     newSegment->ResetMaterialSampling();
     newSegment->SetVertices(LastVertex->GetPrev(), LastVertex);
      
     //connect it to the last segment in the line
     LastSegment->SetNext(newSegment);
     newSegment->SetPrev(LastSegment);
     newSegment->SetNext(NULL);
      
     //now the new segment is the last one
     LastSegment = newSegment;
      
     //update house-keeping data
     nSegment++;
     TotalLength += LastSegment->GetLength();

     SegmentPointerTable.push_back(LastSegment);

     return LastVertex;
   }

   void cFieldLine::CutFront(int nDeletedSegments) {

     auto RemoveSingleSegment = [&] (cFieldLineSegment* s) {
       if ((s==FirstSegment)&&(s==LastSegment)) {
         //delete particles that attached to the segment 
         s->DeleteAttachedParticles();

         //remove the segment 
         VerticesAll.deleteElement(FirstVertex);
         VerticesAll.deleteElement(LastVertex);
         SegmentsAll.deleteElement(s);

         FirstVertex=NULL,LastVertex=NULL;
         FirstSegment=NULL,LastSegment=NULL;
 
         nSegment=0,TotalLength=0.0,IsSet=0;
         SegmentPointerTable.clear();
       }
       else {
         //remove segment at the beginning of the line
         TotalLength-=s->GetLength();
         nSegment--;

         //delete particles that attached to the segment
         s->DeleteAttachedParticles();

         auto next_vertex=FirstVertex->GetNext();
         VerticesAll.deleteElement(FirstVertex);
         FirstVertex=next_vertex;
         FirstVertex->SetPrev(NULL);

         FirstSegment=FirstSegment->GetNext();
         SegmentsAll.deleteElement(s);

         FirstSegment->SetPrev(NULL);
         SegmentPointerTable.erase(SegmentPointerTable.begin());
       }
     };

     for (int i=0;i<nDeletedSegments;i++) if (nSegment!=0) RemoveSingleSegment(FirstSegment);
   }


   void cFieldLine::CutBack(int nDeletedSegments) {
     auto RemoveSingleSegment = [&] (cFieldLineSegment* s) {
       if ((s==FirstSegment)&&(s==LastSegment)) {
         //delete particles that attached to the segment
         s->DeleteAttachedParticles();
         
         //the line consists only of one segment -> remove all segments of the line
         VerticesAll.deleteElement(FirstVertex);
         VerticesAll.deleteElement(LastVertex);
         SegmentsAll.deleteElement(s);

         FirstVertex=NULL,LastVertex=NULL;
         FirstSegment=NULL,LastSegment=NULL;

         nSegment=0,TotalLength=0.0,IsSet=0;
         SegmentPointerTable.clear();
       }
       else {
         //delete particles that attached to the segment
         s->DeleteAttachedParticles();
         
         //remove segment at the end of the line
         TotalLength-=s->GetLength();
         nSegment--;

         auto prev_vertex=LastVertex->GetPrev();
         VerticesAll.deleteElement(LastVertex);
         LastVertex=prev_vertex;
         LastVertex->SetNext(NULL);

         LastSegment=LastSegment->GetPrev();
         SegmentsAll.deleteElement(s);

         LastSegment->SetNext(NULL);
         SegmentPointerTable.pop_back();
       }
     };

     for (int i=0;i<nDeletedSegments;i++) if (nSegment!=0) RemoveSingleSegment(LastSegment);
   }




    //=========================================================================
    void cFieldLine::ResetSegmentWeights() {
      cFieldLineSegment *Segment = FirstSegment;
      double w[PIC::nTotalSpecies];

      //compute weights and normalize them
      for(int spec=0; spec < PIC::nTotalSpecies; spec++) TotalWeight[spec] = 0;

      for(int iSegment=0; iSegment<nSegment; iSegment++) {
        //compute weights
        _FIELDLINE_SEGMENT_WEIGHT_(w, Segment);

        //set segment's weight
        Segment->SetWeight(w);
        for(int spec=0; spec < PIC::nTotalSpecies; spec++) TotalWeight[spec] += w[spec];

        //get the next one
        Segment = Segment->GetNext();
      }
    }

    //=========================================================================
    // choose a random segment on field line and return pointer SegmentOut;
    // since statistical weights may vary, compute a correction factor
    // for a particle to be injected on this segment
    //-----------------------------------------------------------------------
    // choose uniformly a segment, i.e. weight_uniform = 1.0 / nSegment
    void cFieldLine::GetSegmentRandom(int& iSegment,double& WeightCorrectionFactor,int spec) {
      int iSegmentChoice = (int)(nSegment * rnd());

      WeightCorrectionFactor = 1.0;
      iSegment = iSegmentChoice;
      // cycle through segment until get the chosen one
      //      SegmentOut = FirstSegment;
      //      for(int iSegment=0; iSegment < iSegmentChoice; iSegment++)
      //	SegmentOut = SegmentOut->GetNext();

      // SegmentOut now is a pointer to a uniformly chosen segment;
      // find the correction factor: weight_segment / weight_uniform
      //-----------------------------------------------------------------------
      //      WeightCorrectionFactor = 
      //      	SegmentOut->GetWeight(spec) / TotalWeight[spec] * nSegment;
    }
    
    //=========================================================================
    long int InjectParticle(int spec) {
      // this is a wrapper that can call either the default injection procedure
      // or a user-defined procedure
      //      return InjectParticle_default(spec);
      return _PIC_INJECT_PARTICLE_ONTO_FIELD_LINE_(spec);
    }
    
    //=========================================================================
 long int InjectParticle_default(int spec) {
   double WeightCorrection,p[3],v[3];
   int iFieldLine,iSegment;

   iFieldLine = (int)(nFieldLine * rnd());
   FieldLinesAll[iFieldLine].GetSegmentRandom(iSegment,WeightCorrection, spec);

   int q=3;
   double vmin=1e6, vmax=1e6;
   double pvmin = pow(vmin, 1-q), pvmax = pow(vmax, 1-q);
   double r= rnd();
   double absv = pow( (1-r)*pvmin + r*pvmax, 1.0/(1-q));

   double m0=PIC::MolecularData::GetMass(spec);

   switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
   case _PIC_MODE_OFF_:
     Vector3D::Distribution::Uniform(p,absv*m0);
     break;

   case _PIC_MODE_ON_:
     Vector3D::Distribution::Uniform(v,absv);
     Relativistic::Vel2Momentum(p,v,PIC::MolecularData::GetMass(spec));
     break;
   }

   return InjectParticle_default(spec,p,1.0,iFieldLine,iSegment);
}



    long int InjectParticle_default(int spec,double *p,double ParticleWeightCorrectionFactor,int iFieldLine,int iSegment,double sIn) {
      //namespace alias
      namespace PB = PIC::ParticleBuffer;

      // pointer particle to the particle to be injected
      PB::byte ptrData[PB::ParticleDataLength];
      
      //default settings
      PB::SetIndividualStatWeightCorrection(ParticleWeightCorrectionFactor,ptrData);
      PB::SetI(spec, ptrData);
      
      // pick a random field line
      PB::SetFieldLineId(iFieldLine, ptrData);

      // inject particle onto this field line
      double WeightCorrection;
      //    FieldLinesAll[iFieldLine].GetSegmentRandom(iSegment,
      //					       WeightCorrection, spec);
      
      //Inject at the beginning of the field line FOR PARKER SPIRAL
      cFieldLineSegment* Segment=FieldLinesAll[iFieldLine].GetSegment(iSegment);
      double S = (sIn>=0) ? sIn : iSegment + rnd();


      PB::SetFieldLineCoord(S, ptrData);
      double x[3], v[3];

      FieldLinesAll[iFieldLine].GetCartesian(x, S);
      PB::SetX(x, ptrData);

      //verify that the 'x' belongs to the subdomain of the current MPI process
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->findTreeNode(x,NULL);       

      if (_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_NODE_) {
        if (node==NULL) return -1;
        if (node->block==NULL) return -1;
        if (node->Thread!=PIC::ThisThread) return -1;  
      }	
      
      // fix kintic energy and pitch angle far now
      //    double cosPhi=pow(2,-0.5);
      //    double vpar = 1e6 * cosPhi;
      //    double KinEnergyPerp=5E-17*(1-cosPhi*cosPhi);
      //direction of velocity: INJECT ALONG HE FIELD LINE FOR PARKER SPIRAL
      //    vpar*=(rnd()<0.5)?1:-1;
      
      // inject power law, 
      // see Numerical Studies of the Solar Energetic Particle 
      // Transport and Acceleration, Tenishev et al.
      //int q=3;
      //double vmin=1e6, vmax=1e6;
      //double pvmin = pow(vmin, 1-q), pvmax = pow(vmax, 1-q);
      //double r= rnd();
      //double absv = pow( (1-r)*pvmin + r*pvmax, 1.0/(1-q));

      //direction of velocity: INJECT ALONG HE FIELD LINE FOR PARKER SPIRAL
      //double cosPhi =  1 - 2*rnd();

      //double vpar = absv*cosPhi;
      //double KinEnergyPerp = 9.1E-31 * 0.5 * (absv*absv - vpar*vpar);


      //velocity is paral to the field line
      //Segment->GetDir(v);
      //for (int i=0; i<3; i++) v[i]*=vpar;
      //PB::SetV(v, ptrData);
     

      int idim;
      double vpar[3],vperp[3],l[3],c0=0.0,m0=PIC::MolecularData::GetMass(spec);
      double pPerpAbs,pParAbs,pPerpAbs2=0.0,pParAbs2=0.0;

      Segment->GetDir(l);

      switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
      case _PIC_MODE_OFF_:
        for (c0=0.0,idim=0;idim<3;idim++) {
          v[idim]=p[idim]/m0;
          c0+=p[idim]*l[idim];
        }

        for (pPerpAbs2=0.0,pParAbs2=0.0,idim=0;idim<3;idim++) {
          double t;

          t=p[idim]-c0*l[idim];
          pPerpAbs2+=t*t;

          t=c0*l[idim];
          pParAbs2+=t*t; 
        }

        pPerpAbs=sqrt(pPerpAbs2);
        pParAbs=sqrt(pParAbs2);
        break;

      case _PIC_MODE_ON_:
        ::Relativistic::Momentum2Vel(v,p,m0);
        c0=Vector3D::DotProduct(p,l);

        for (pPerpAbs2=0.0,pParAbs2=0.0,idim=0;idim<3;idim++) {
          double t;

          t=p[idim]-c0*l[idim];
          pPerpAbs2+=t*t;

          t=c0*l[idim];
          pParAbs2+=t*t;
        } 
        
        pPerpAbs=sqrt(pPerpAbs2);
        pParAbs=sqrt(pParAbs2);
        break;
      }


      //PB::SetV(v,ptrData);

      double vParallel,vNormal;

      Vector3D::GetComponents(vParallel,vNormal,v,l);
      PB::SetVParallel(vParallel,ptrData);
      PB::SetVNormal(vNormal,ptrData);

      PB::SetMomentumParallel(pParAbs,ptrData);
      PB::SetMomentumNormal(pPerpAbs,ptrData);
 
      //magnetic field
      double B[3], AbsB;
      Segment->GetMagneticField(S-(int)S, B);
      AbsB = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2], 0.5)+1E-15;
      
      //magnetic moment
      double mu=pPerpAbs2/(2.0*m0*AbsB);
      PB::SetMagneticMoment(mu, ptrData);

//      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
//      node=PIC::Mesh::mesh->findTreeNode(x);
      
      //generate a new particles
      long int NewParticle;
      PB::byte* NewParticleData;

      switch (_PIC_PARTICLE_LIST_ATTACHING_) {
      case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
        NewParticle=PB::GetNewParticle(Segment->FirstParticleIndex);
        NewParticleData=PB::GetParticleDataPointer(NewParticle); 

        PB::CloneParticle(NewParticleData,ptrData);
        PB::SetParticleAllocated(NewParticleData);
        break;
      case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
        NewParticle=PB::InitiateParticle(x,NULL,NULL,NULL,ptrData,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }  

      return NewParticle;
    }

    //=========================================================================
    void cFieldLine::SetMagneticField(double *BIn, int iVertex){
      cFieldLineVertex *Vertex;
      int nVertex = (is_loop()) ? nSegment : nSegment+1;

      if (iVertex == -1) Vertex = LastVertex;
      else if ((iVertex > 0.5*nSegment) && (iVertex <= nVertex)) {
        Vertex = LastVertex;

        for (int i=nSegment; i>iVertex; i--) Vertex = Vertex->GetPrev();
      }
      else if (iVertex >= 0) {
        Vertex = FirstVertex;

        for (int i=0; i<iVertex; i++) Vertex = Vertex->GetNext();
      }
      else exit(__LINE__, __FILE__, "ERROR: invalid index of vertex");

      Vertex->SetMagneticField(BIn);
    }
    //=========================================================================
    void cFieldLine::GetMagneticField(double* BOut, double S) {
      // check correctness
      if ((S < 0.0) || (S > nSegment))exit(__LINE__,__FILE__,"ERROR: trying to get magnetic field at an invalid location");

      // interpolate the magnetic field at the location S:
      //  floor(S) is the number of the segment,
      //  S - floor(S) is the location along segment (between 0 & 1)
      //-----------------------------------------------------------------------

      // number of the begin vertex
      int iSegment = (int) S;

      cFieldLineSegment* Segment = GetSegment(iSegment);
      Segment->GetMagneticField(S - iSegment, BOut);
    }    

    //=========================================================================
    void cFieldLine::GetPlasmaVelocity(double* VelOut, double S) {
      // check correctness
      if ((S < 0.0) || (S > nSegment))exit(__LINE__,__FILE__,"ERROR: trying to get magnetic field at an invalid location");
    
      // interpolate the magnetic field at the location S:
      //  floor(S) is the number of the segment,
      //  S - floor(S) is the location along segment (between 0 & 1)
      //-----------------------------------------------------------------------
   
      // number of the begin vertex
      int iSegment = (int) S;

      cFieldLineSegment* Segment = GetSegment(iSegment);
      Segment->GetPlasmaVelocity(S - iSegment, VelOut);
    }

    double cFieldLine::GetPlasmaDensity(double S) {
      // check correctness
      if ((S < 0.0) || (S > nSegment))exit(__LINE__,__FILE__,"ERROR: trying to get magnetic field at an invalid location");

      // interpolate the magnetic field at the location S:
      //  floor(S) is the number of the segment,
      //  S - floor(S) is the location along segment (between 0 & 1)
      //-----------------------------------------------------------------------

      // number of the begin vertex
      int iSegment = (int) S;
      double VelOut; 

      cFieldLineSegment* Segment = GetSegment(iSegment);
      Segment->GetPlasmaDensity(S - iSegment, VelOut);

      return VelOut;
    }


    //=========================================================================
    void cFieldLine::Output(FILE* fout, bool OutputGeometryOnly,bool collective) {
      namespace PB = PIC::ParticleBuffer;

      cFieldLineVertex *Vertex=FirstVertex;
      int nVertex = (is_loop()) ? nSegment : nSegment+1;
      int iVertex,iSegment;
      cFieldLineSegment* Segment;

      double *DistanceTable=new double [nVertex]; 

      DistanceTable[0]=0.0;
      
      double LengthConversionFactor=PIC::FieldLine::cFieldLine::OutputLengthConversionFactor.first;

      //calculate length of vertexes from the beginning of the field line 
      for (iSegment=0,Segment=FirstSegment; iSegment<nSegment; iSegment++,Segment=Segment->GetNext()) {
        DistanceTable[iSegment+1]=DistanceTable[iSegment]+Segment->GetLength();
      } 

      // Create mapping from vertices to adjacent segments for edge data averaging
      std::vector<std::vector<cFieldLineSegment*>> vertexToSegmentsMap(nVertex);

      if (!OutputGeometryOnly) {
      // Build the vertex-to-segments mapping
        for (iSegment=0, Segment=FirstSegment; iSegment<nSegment; iSegment++, Segment=Segment->GetNext()) {
          // Each segment connects vertex iSegment to vertex iSegment+1
          vertexToSegmentsMap[iSegment].push_back(Segment);

          if (iSegment+1 < nVertex) {
            vertexToSegmentsMap[iSegment+1].push_back(Segment);
          }

          // Special handling for closed loops
          if (is_loop() && iSegment == nSegment-1) {
            // Last segment connects back to first vertex
            vertexToSegmentsMap[0].push_back(Segment);
          }
        }
      }


      if (OutputGeometryOnly==true) { 
        for (int iVertex=0; iVertex<nVertex; iVertex++) {
          //print coordinates
          double x[DIM];

          Vertex->GetX(x);

          for (int idim=0; idim<DIM; idim++) fprintf(fout, "%e ", LengthConversionFactor*x[idim]);

          fprintf(fout,"%e \n",LengthConversionFactor*DistanceTable[iVertex]);
          Vertex = Vertex->GetNext();
        }

        delete [] DistanceTable;
        return;
      }
      
      // the longest datum printed at a vertex: sizes the per-vertex value buffer below
      int MaxVertexDatumLength=3;

      for (auto d : DataStoredAtVertex) if (d->length>MaxVertexDatumLength) MaxVertexDatumLength=d->length;
      for (auto d : DataSampledAtVertex) if (d->length>MaxVertexDatumLength) MaxVertexDatumLength=d->length;

      // Counts and speed SUMS belong to the actual particle owner. Reduce
      // them before computing a mean; averaging owner-local means would give
      // empty ranks equal statistical weight. Use wide counters for MPI sums.
      std::vector<long long> nSegmentParticles(nSegment),nSegmentParticles_mu_negative(nSegment),nSegmentParticles_mu_positive(nSegment);
      std::vector<double> ModelParticleSpeedTable(nSegment);double v[3];

      for (iSegment=0;iSegment<nSegment; iSegment++) {
        nSegmentParticles[iSegment]=0,nSegmentParticles_mu_negative[iSegment]=0,nSegmentParticles_mu_positive[iSegment]=0;
        ModelParticleSpeedTable[iSegment]=0.0; 
      }
   

      if (_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_) {
        for (iSegment=0,Segment=FirstSegment; iSegment<nSegment; iSegment++,Segment=Segment->GetNext()) {
          nSegmentParticles[iSegment]=0; 
          ModelParticleSpeedTable[iSegment]=0.0;
     
          long int ptr=Segment->FirstParticleIndex;

          while (ptr!=-1) {
            nSegmentParticles[iSegment]++;

            //v=PB::GetV(ptr);

            v[0]=PB::GetVParallel(ptr);
            v[1]=PB::GetVNormal(ptr); 
            v[2]=0.0;


            ModelParticleSpeedTable[iSegment]+=Vector3D::Length(v);

            if (v[0]>0.0) {
              nSegmentParticles_mu_positive[iSegment]++;
            }
            else {
              nSegmentParticles_mu_negative[iSegment]++;
            }

            ptr=PB::GetNext(ptr);
          }

        }
      }

      std::vector<double> fluence;
      if(collective) {
        for(auto* counts:{&nSegmentParticles,&nSegmentParticles_mu_negative,&nSegmentParticles_mu_positive})
          MPI_Allreduce(MPI_IN_PLACE,counts->data(),nSegment,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(MPI_IN_PLACE,ModelParticleSpeedTable.data(),nSegment,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        if(DatumAtVertexFluence.offset>=0) {
          fluence.resize(nVertex);auto* vertex=FirstVertex;
          for(int i=0;i<nVertex;++i,vertex=vertex->GetNext())vertex->GetDatum(DatumAtVertexFluence,fluence[i]);
          MPI_Allreduce(MPI_IN_PLACE,fluence.data(),nVertex,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          // This is a temporary output view. Writing a global cumulative sum
          // back to every local primitive would count the same past fluence
          // again at the next output. The checkpoint retains local integrals.
          for(double value:fluence)if(!std::isfinite(value))
            exit(__LINE__,__FILE__,"nonfinite collective field-line fluence");
        }
      }
      for(int i=0;i<nSegment;++i)if(nSegmentParticles[i]!=0)
        ModelParticleSpeedTable[i]/=nSegmentParticles[i];


      for (iVertex=0; iVertex<nVertex; iVertex++) {
        //print coordinates
        double x[DIM];
        // A vertex datum can be longer than a 3-vector (e.g., the energy-binned
        // number density DatumAtVertexNumberDensityEnergyBinned has 10 bins);
        // a fixed 'double Value[3]' was overrun by GetDatum*() (stack smashing)
        std::vector<double> ValueBuffer(MaxVertexDatumLength);
        double *Value=ValueBuffer.data();
        vector<cDatumStored*>::iterator itrDatumStored;

        Vertex->GetX(x);

        for (int idim=0; idim<DIM; idim++) fprintf(fout, "%e ", LengthConversionFactor*x[idim]);

        fprintf(fout, "%e ", LengthConversionFactor*DistanceTable[iVertex]);

        if (_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_) {
          if (iVertex==0) {
            fprintf(fout, "%e  %e  %e  %e  %d  ",(double)nSegmentParticles[0],
              (double)nSegmentParticles_mu_positive[0],(double)nSegmentParticles_mu_negative[0],ModelParticleSpeedTable[0],iVertex);
          }
          else if (iVertex==nVertex-1) {
            fprintf(fout, "%e  %e  %e  %e %d  ",(double)nSegmentParticles[iVertex-1],
              (double)nSegmentParticles_mu_positive[iVertex-1],(double)nSegmentParticles_mu_negative[iVertex-1],ModelParticleSpeedTable[iVertex-1],iVertex); 
          }
          else {
            fprintf(fout, "%e  %e  %e  %e  %d  ",0.5*(nSegmentParticles[iVertex-1]+nSegmentParticles[iVertex]), 
              0.5*(nSegmentParticles_mu_positive[iVertex-1]+nSegmentParticles_mu_positive[iVertex]),
              0.5*(nSegmentParticles_mu_negative[iVertex-1]+nSegmentParticles_mu_negative[iVertex]),
              0.5*(ModelParticleSpeedTable[iVertex-1]+ModelParticleSpeedTable[iVertex]),iVertex);
          }
        }

        for (itrDatumStored = DataStoredAtVertex.begin();itrDatumStored!= DataStoredAtVertex.end();  itrDatumStored++) {
          if ((*itrDatumStored)->doPrint) {
            Vertex->GetDatum(*(*itrDatumStored), Value);
            // Vertex B is authoritative SI [T] for coupling and transport.
            // This established output interface explicitly advertises [nT];
            // convert only the serialized value, never the stored primitive.
            // The previous-state vector has the same units. Other datums keep
            // their existing independent output contract. PrintName supplies
            // separators, so vector labels above must not add a trailing comma.
            if(collective&&(*itrDatumStored)==&DatumAtVertexFluence)Value[0]=fluence.at(iVertex);
            const bool magnetic=(*itrDatumStored)==&DatumAtVertexMagneticField||
                (*itrDatumStored)==&DatumAtVertexPrevious::DatumAtVertexMagneticField;
            for (int i=0; i<(*itrDatumStored)->length; i++) fprintf(fout, "%e ", magnetic?1e9*Value[i]:Value[i]);
          }
        }

        vector<cDatumSampled*>::iterator itrDatum;

        for(itrDatum = DataSampledAtVertex.begin();itrDatum!= DataSampledAtVertex.end();itrDatum++) if ((*itrDatum)->doPrint) {
          cDatumTimed*    ptrDatumTimed;
          cDatumWeighted* ptrDatumWeighted;

          if(ExternallyManagedGeometry) {
            // Keep the original endpoint shape-function projection, but divide
            // each material history by its own exposure BEFORE adding it to a
            // neighboring history. For weighted moments both numerator and
            // sampled weight receive that division. Equal lifetimes reduce to
            // the original static projection exactly, without age dilution.
            const auto& datum=**itrDatum;
            const bool weighted=datum.type==PIC::Datum::cDatumSampled::Weighted_;
            const std::size_t width=cFieldLineVertex::SampleDataBytes()/sizeof(double);
            for(int component=0;component<datum.length;++component) {
              double numerator=0,denominator=0;bool measured=false;
              for(auto* edge:vertexToSegmentsMap[iVertex]) {
                if(!edge->SamplingExposure.HasCompletedMeasurement())continue;
                measured=true;
                const std::size_t endpoint=edge->GetBegin()==Vertex?0:1;
                const double count=static_cast<double>(edge->SamplingExposure.completed);
                numerator+=edge->MaterialCompleted.at(endpoint*width+datum.offset/sizeof(double)+component)/count;
                if(weighted)denominator+=edge->MaterialCompleted.at(endpoint*width+
                    DatumAtVertexParticleWeight.offset/sizeof(double))/count;
              }
              Value[component]=!measured?std::numeric_limits<double>::quiet_NaN():
                  (weighted?(denominator>0?numerator/denominator:0):numerator);
            }
          }else if ((*itrDatum)->type == PIC::Datum::cDatumSampled::Timed_) {
            ptrDatumTimed = static_cast<cDatumTimed*> ((*itrDatum));
            Vertex->GetDatumAverage(*ptrDatumTimed, Value, 0,SampleCycleCounter);
          }
          else {
            ptrDatumWeighted = static_cast<cDatumWeighted*> ((*itrDatum));
            Vertex->GetDatumAverage(*ptrDatumWeighted, Value, 0);
          }

          for(int i=0; i<(*itrDatum)->length; i++) {
            fprintf(fout, "%e ", Value[i]);
          }
        }

    for (auto it : PIC::FieldLine::cFieldLineSegment::DataStoredAtSegment) {
      if (it->doPrint) {
        // Initialize averaging arrays
        std::vector<double> avgValue(it->length, 0.0);
        double totalWeight = 0.0;
        
        // Average over all segments connected to this vertex
        for (auto segment : vertexToSegmentsMap[iVertex]) {
          double segmentValue[it->length];
          segment->GetDatum(*it, segmentValue);
          
          // Weight by segment length for better averaging
          double weight = segment->GetLength();
          totalWeight += weight;
          
          for (int i = 0; i < it->length; i++) {
            avgValue[i] += weight * segmentValue[i];
          }
        }
        
        // Normalize by total weight
        if (totalWeight > 0.0) {
          for (int i = 0; i < it->length; i++) {
            avgValue[i] /= totalWeight;
            fprintf(fout, "%e ", avgValue[i]);
          }
        } else {
          // No adjacent segments - output zeros
          for (int i = 0; i < it->length; i++) {
            fprintf(fout, "%e ", 0.0);
          }
        }
      }
    }

    // Output segment-sampled data (averaged over adjacent segments)
    for (auto it : PIC::FieldLine::cFieldLineSegment::DataSampledAtSegment) {
      if (it->doPrint) {
        // Initialize averaging arrays
        std::vector<double> avgValue(it->length, 0.0);
        double totalWeight = 0.0;
        
        // Average over all segments connected to this vertex
        for (auto segment : vertexToSegmentsMap[iVertex]) {
          double segmentValue[it->length];
          
          // Handle different sampling types
          if (it->type == PIC::Datum::cDatumSampled::Timed_) {
            cDatumTimed* ptrDatumTimed = static_cast<cDatumTimed*>(it);
            // Segment-local sampled datums have their own actual lifetime,
            // even when their vertex projection uses geometric length weights.
            const auto count=segment->SamplingExposure.completed;
            if(count>static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
              exit(__LINE__,__FILE__,"segment sampling count exceeds legacy datum API");
            segment->GetDatumAverage(*ptrDatumTimed, segmentValue, 0, static_cast<int>(count));
          } else if (it->type == PIC::Datum::cDatumSampled::Weighted_) {
            cDatumWeighted* ptrDatumWeighted = static_cast<cDatumWeighted*>(it);
            segment->GetDatumAverage(*ptrDatumWeighted, segmentValue, 0);
          } else {
            exit(__LINE__,__FILE__,"Error: not implemented");
          }
          
          // Weight by segment length
          double weight = segment->GetLength();
          totalWeight += weight;
          
          for (int i = 0; i < it->length; i++) {
            avgValue[i] += weight * segmentValue[i];
          }
        }
        
        // Normalize by total weight
        if (totalWeight > 0.0) {
          for (int i = 0; i < it->length; i++) {
            avgValue[i] /= totalWeight;
            fprintf(fout, "%e ", avgValue[i]);
          }
        } else {
          // No adjacent segments - output zeros
          for (int i = 0; i < it->length; i++) {
            fprintf(fout, "%e ", 0.0);
          }
        }
      }
    }


        fprintf(fout,"\n");

        //flush completed sampling buffer
        Vertex->flushCollectingSamplingBuffer();
        Vertex = Vertex->GetNext();
      }

      delete [] DistanceTable;
    }
    
    //=========================================================================
    void Init() {

      static bool init_flag=false;

      if (init_flag==true) return;
      init_flag=true;
      
    //  if(PIC::nTotalThreads > 1) exit(__LINE__, __FILE__,"Not implemented for multiple processors");
      
      // allocate container for field lines
      FieldLinesAll = new cFieldLine [nFieldLineMax];

      //activate segment data storage 
      cFieldLineSegment::InitDatum();

      // activate vertex data storage
      long int Offset = 0;

      // activate data that are stored but NOT sampled
      if (VertexAllocationManager.MagneticField==true)        DatumAtVertexMagneticField.       activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.ElectricField==true)        DatumAtVertexElectricField.       activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PlasmaVelocity==true)       DatumAtVertexPlasmaVelocity.      activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PlasmaDensity==true)        DatumAtVertexPlasmaDensity.       activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PlasmaTemperature==true)    DatumAtVertexPlasmaTemperature.   activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PlasmaPressure==true)       DatumAtVertexPlasmaPressure.      activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.MagneticFluxFunction==true) DatumAtVertexMagneticFluxFunction.activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PlasmaWaves==true)          DatumAtVertexPlasmaWaves.         activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.ShockLocation==true)        DatumAtVertexShockLocation.       activate(Offset, &DataStoredAtVertex);

      if (VertexAllocationManager.PreviousVertexData.MagneticField==true)        DatumAtVertexPrevious::DatumAtVertexMagneticField.       activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PreviousVertexData.ElectricField==true)        DatumAtVertexPrevious::DatumAtVertexElectricField.       activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PreviousVertexData.PlasmaVelocity==true)       DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity.      activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PreviousVertexData.PlasmaDensity==true)        DatumAtVertexPrevious::DatumAtVertexPlasmaDensity.       activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PreviousVertexData.PlasmaTemperature==true)    DatumAtVertexPrevious::DatumAtVertexPlasmaTemperature.   activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PreviousVertexData.PlasmaPressure==true)       DatumAtVertexPrevious::DatumAtVertexPlasmaPressure.      activate(Offset, &DataStoredAtVertex);
      if (VertexAllocationManager.PreviousVertexData.PlasmaWaves==true)          DatumAtVertexPrevious::DatumAtVertexPlasmaWaves.         activate(Offset, &DataStoredAtVertex);

      if (VertexAllocationManager.Fluence==true) DatumAtVertexFluence.activate(Offset, &DataStoredAtVertex);

      if (VertexAllocationManager.DistanceToShockLocation==true) DatumAtVertexShockLocationDistance.activate(Offset, &DataStoredAtVertex);

      // activate data that is sampled
      long int SamplingOffset = Offset;
      Offset = 0;

      DatumAtVertexParticleWeight.  activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleNumber.  activate(Offset, &DataSampledAtVertex);
      DatumAtVertexNumberDensity.   activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleEnergy.  activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleSpeed.   activate(Offset, &DataSampledAtVertex); 

      DatumAtVertexNumberDensityEnergyBinned.   activate(Offset, &DataSampledAtVertex);

      DatumAtVertexParticleCosPitchAngle.   activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleAbsCosPitchAngle.activate(Offset, &DataSampledAtVertex);

      DatumAtVertexNumberDensity_mu_positive.activate(Offset, &DataSampledAtVertex);
      DatumAtVertexNumberDensity_mu_negative.activate(Offset, &DataSampledAtVertex);

      DatumAtVertexParticleFlux_mu_positive.activate(Offset, &DataSampledAtVertex);
      DatumAtVertexParticleFlux_mu_negative.activate(Offset, &DataSampledAtVertex);


      // assign offsets and data length
      cFieldLineVertex::SetDataOffsets(SamplingOffset, Offset);
      
      PIC::IndividualModelSampling::DataSampledList.push_back(&DatumAtGridParticleEnergy);

      //request data in the particle state vector to keep the field line-related data
      long int offset;
      int DataLength;

      DataLength=sizeof(cParticleFieldLineData);
      PIC::ParticleBuffer::RequestDataStorage(offset,DataLength,
                                              alignof(cParticleFieldLineData));

      ParticleDataOffset=offset;
    }
    
    //=========================================================================
    void Output(const char* fname, bool OutputGeometryOnly) {

      if ((OutputGeometryOnly==false)&&(LastDataOutputFileNumber!=DataOutputFileNumber)) {
        // Validate replicated topology before entering per-segment collectives;
        // different loop lengths would otherwise hang or mix physical edges.
        long lineCount=nFieldLine,lineLow=0,lineHigh=0;
        MPI_Allreduce(&lineCount,&lineLow,1,MPI_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(&lineCount,&lineHigh,1,MPI_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
        if(lineLow!=lineHigh)exit(__LINE__,__FILE__,"field-line sampling topology disagrees over MPI");
        // Complete each segment's actual lifetime exposure together with its
        // own buffers. Replicated geometry counts must agree over MPI; summing
        // them would multiply physical exposure by the process count.
        cFieldLineSegment::swapSamplingBuffers();
        for(int line=0;line<nFieldLine;++line) {
          auto& fieldLine=FieldLinesAll[line];
          int shape[2]={fieldLine.GetTotalSegmentNumber(),fieldLine.ExternallyManagedGeometry?1:0},shapeLow[2],shapeHigh[2];
          MPI_Allreduce(shape,shapeLow,2,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
          MPI_Allreduce(shape,shapeHigh,2,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
          if(shapeLow[0]!=shapeHigh[0]||shapeLow[1]!=shapeHigh[1])
            exit(__LINE__,__FILE__,"material sampling owner/segment counts disagree over MPI");
          for(auto* edge=fieldLine.GetFirstSegment();edge;edge=edge->GetNext()) {
            edge->SamplingExposure.Complete();
            edge->flushCollectingSamplingBuffer();
            if(!fieldLine.ExternallyManagedGeometry)continue;
            unsigned long long count=edge->SamplingExposure.completed,low=0,high=0;
            MPI_Allreduce(&count,&low,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
            MPI_Allreduce(&count,&high,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
            if(low!=high)exit(__LINE__,__FILE__,"material segment sampling exposures disagree over MPI");
            double duration[2]={edge->SamplingExposure.completedSeconds,edge->SamplingExposure.fixedStepSeconds};
            double durationLow[2],durationHigh[2];
            MPI_Allreduce(duration,durationLow,2,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
            MPI_Allreduce(duration,durationHigh,2,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
            if(durationLow[0]!=durationHigh[0]||durationLow[1]!=durationHigh[1])
              exit(__LINE__,__FILE__,"material sampling durations/timesteps disagree over MPI");
            edge->MaterialCompleted.swap(edge->MaterialCollecting);
            std::fill(edge->MaterialCollecting.begin(),edge->MaterialCollecting.end(),0.0);
            if(!edge->MaterialCompleted.empty()) {
              if(edge->MaterialCompleted.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
                exit(__LINE__,__FILE__,"material sample MPI count overflow");
              MPI_Allreduce(MPI_IN_PLACE,edge->MaterialCompleted.data(),
                  static_cast<int>(edge->MaterialCompleted.size()),MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
              for(double value:edge->MaterialCompleted)if(!std::isfinite(value))
                exit(__LINE__,__FILE__,"material sample reduction produced a nonfinite numerator");
            }
          }
        }
        // Legacy vertex buffers retain their existing layout. Material output
        // below uses separately normalized segment-endpoint histories instead.
        cFieldLineVertex::swapSamplingBuffers();
	LastDataOutputFileNumber=DataOutputFileNumber;

        PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(&PIC::FieldLine::DatumAtVertexParticleWeight);
      
        for(auto itrDatum = DataSampledAtVertex.begin();itrDatum!= DataSampledAtVertex.end();itrDatum++) if ((*itrDatum)->doPrint) {
          cDatumTimed*    ptrDatumTimed;
          cDatumWeighted* ptrDatumWeighted;
 
          if ((*itrDatum)->type == PIC::Datum::cDatumSampled::Timed_) {
            ptrDatumTimed = static_cast<cDatumTimed*> ((*itrDatum));
	    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(ptrDatumTimed);
          }
          else {
            ptrDatumWeighted = static_cast<cDatumWeighted*> ((*itrDatum));
	    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(ptrDatumWeighted);
          }
	}
  
        //execute a user-defined function to processes sampled data before output
        if (UserDefinedfDataProcessingManager!=NULL) UserDefinedfDataProcessingManager();
      }


      // Reductions and per-line traversal remain collective, but the shared
      // destination has exactly one writer. Other ranks serialize into their
      // private /dev/null stream so existing local formatting callbacks remain
      // valid without concurrent fopen/truncation of the requested data file.
      const std::string stagedOutput=std::string(fname)+".tmp";
      FILE* fout=std::fopen(PIC::ThisThread==0?stagedOutput.c_str():"/dev/null","w");
      if(!fout)exit(__LINE__,__FILE__,"cannot stage field-line output");
      
      fprintf(fout, "TITLE=\"Field line geometry\"\n");

      //if no units are defined in OutputLengthConversionFactor, then set the default to [m] 
      if (PIC::FieldLine::cFieldLine::OutputLengthConversionFactor.second=="") {
        PIC::FieldLine::cFieldLine::OutputLengthConversionFactor.second="m";
      }

      auto *units=PIC::FieldLine::cFieldLine::OutputLengthConversionFactor.second.c_str();
      
switch (DIM) {
  case 3:
    {
      fprintf(fout,"VARIABLES=\"x[%s]\",\"y[%s]\",\"z[%s]\",\"Distance from the beginning[%s]\"",units,units,units,units);
      
      vector<cDatumStored*>::iterator itrDatumStored;
      vector<cDatumSampled*>::iterator itrDatum;

      if (OutputGeometryOnly==false) {
        if (_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_) {
          fprintf(fout,",\"nSegmentParticles\",\"nSegmentParticles_mu_positive\", \"nSegmentParticles_mu_negative\",\"Model Particle Speed\",\"iSegment\"");
        }

        for (itrDatumStored = DataStoredAtVertex.begin();itrDatumStored!= DataStoredAtVertex.end(); itrDatumStored++) {
          if ((*itrDatumStored)->doPrint) (*itrDatumStored)->PrintName(fout);
        }

        for (itrDatum = DataSampledAtVertex.begin(); itrDatum!= DataSampledAtVertex.end(); itrDatum++) {
          if ((*itrDatum)->doPrint) (*itrDatum)->PrintName(fout);
        }

	for (auto it : PIC::FieldLine::cFieldLineSegment::DataStoredAtSegment) {
          if (it->doPrint==true) it->PrintName(fout);  
	}

	for (auto it : PIC::FieldLine::cFieldLineSegment::DataSampledAtSegment) {
          if (it->doPrint==true) it->PrintName(fout);
        }
      }

      fprintf(fout,"\n");

      //output title of the output file if such is defined 
      if (UserDefinedTecplotFileTitle!=NULL) {
        char title[_MAX_STRING_LENGTH_PIC_];

        UserDefinedTecplotFileTitle(title);
        fprintf(fout, "TITLE=\"%s\"\n",title);
      } 
    }
      break;
  default:
      exit(__LINE__,__FILE__,"not implemented");
      break;
}
      
      for (int iFieldLine=0; iFieldLine<nFieldLine; iFieldLine++) {
        fprintf(fout,"ZONE T=\"Field-line %i\" F=POINT\n",iFieldLine);
        FieldLinesAll[iFieldLine].Output(fout, OutputGeometryOnly,!OutputGeometryOnly);
      }
      
      //reset the sanple counpter 
      SampleCycleCounter=0;

      bool outputFailed=std::ferror(fout)!=0;
      outputFailed=std::fclose(fout)!=0||outputFailed;
      if(PIC::ThisThread==0&&!outputFailed)outputFailed=std::rename(stagedOutput.c_str(),fname)!=0;
      if(outputFailed)exit(__LINE__,__FILE__,"cannot publish complete field-line output");
      if(!OutputGeometryOnly&&PIC::ThisThread==0) {
        bool material=false;
        for(int line=0;line<nFieldLine;++line)material=material||FieldLinesAll[line].ExternallyManagedGeometry;
        if(material) {
          // Export actual lifetimes independently of Tecplot's vertex rows.
          // A completed count of zero is explicitly "no measurement", not a
          // vacuum density. Atomic replacement keeps this metadata file
          // from being mistaken for a complete sample after a failed write.
          const std::string destination=std::string(fname)+".sampling";
          const std::string temporary=destination+".tmp";
          FILE* metadata=std::fopen(temporary.c_str(),"w");
          if(!metadata)exit(__LINE__,__FILE__,"cannot open material sampling metadata");
          bool failed=std::fprintf(metadata,
              "# material-segment exposure v1; SI seconds; count=0 means no measurement\n"
              "# line material_label geometry_epoch_s generation completed_samples completed_seconds fixed_dt_s valid\n")<0;
          for(int line=0;line<nFieldLine;++line) {
            const auto& geometry=FieldLinesAll[line];if(!geometry.ExternallyManagedGeometry)continue;
            int index=0;
            for(auto* edge=FieldLinesAll[line].GetFirstSegment();edge;edge=edge->GetNext(),++index) {
              const auto& e=edge->SamplingExposure;
              failed=std::fprintf(metadata,"%d %lld %.17g %llu %llu %.17g %.17g %d\n",line,
                  static_cast<long long>(geometry.FirstMaterialLabel+index),geometry.GeometryEpochS,
                  static_cast<unsigned long long>(geometry.GeometryGeneration),
                  static_cast<unsigned long long>(e.completed),e.completedSeconds,e.fixedStepSeconds,
                  e.HasCompletedMeasurement()?1:0)<0||failed;
            }
          }
          failed=std::fclose(metadata)!=0||failed;
          if(failed||std::rename(temporary.c_str(),destination.c_str())!=0)
            exit(__LINE__,__FILE__,"cannot publish complete material sampling metadata");
        }
      }
    }

    //=========================================================================
    void Sampling() {
      SampleCycleCounter++;

      int iFieldLine,spec;
      PIC::FieldLine::cFieldLineSegment *Segment;
      long int ptr;
      double LocalParticleWeight;

      if (_SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_) {
        exit(__LINE__,__FILE__,"Error: the function is implemented only for the case _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_. Changed input file.");
      }

      /*
      if (DatumAtVertexFluence.is_active()==true) {
        if (_SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_) {
           exit(__LINE__,__FILE__,"Error: the function is implemented only for the case _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_. Change input file."); 
        }
      }
      */


      for (iFieldLine=0;iFieldLine<nFieldLine;iFieldLine++) {
        for (Segment=FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
          // One exposure tick per eligible segment, including a vacuum.
          // Never count particles, species or mover substeps as extra samples.
          Segment->SamplingExposure.Sample(PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);
          if(FieldLinesAll[iFieldLine].ExternallyManagedGeometry) {
            const std::size_t width=2*cFieldLineVertex::SampleDataBytes()/sizeof(double);
            if(Segment->MaterialCollecting.size()!=width)
              Segment->MaterialCollecting.assign(width,0.0);
          }
          // Joined mover lists have a scheduling-dependent insertion order.
          // Material restart retains accumulated binary64 samples, so even a
          // one-ULP change in addition order breaks exact continuation while
          // particle trajectories remain identical. Sample a READ-ONLY index
          // in a deterministic order of the physical inputs used below:
          // species, line coordinate, parallel/normal velocity and individual
          // weight. Equal keys produce identical contributions. No particle,
          // list, sampling formula, exposure count or physical dt is changed.
          // Compare raw scalar bits to define a total order even for invalid
          // values; the ordinary sampler's validity checks still own rejection.
          // Keep this key synchronized with any future per-particle sampling
          // input. Nonmaterial legacy lines retain their historical traversal.
          std::vector<std::pair<std::array<std::uint64_t,5>,long int>> materialOrder;
          if(FieldLinesAll[iFieldLine].ExternallyManagedGeometry) {
            for(ptr=Segment->FirstParticleIndex;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
              const double values[4]={PIC::ParticleBuffer::GetFieldLineCoord(ptr),
                  PIC::ParticleBuffer::GetVParallel(ptr),PIC::ParticleBuffer::GetVNormal(ptr),
                  PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr)};
              std::array<std::uint64_t,5> key{{static_cast<std::uint64_t>(PIC::ParticleBuffer::GetI(ptr)),0,0,0,0}};
              static_assert(sizeof(double)==sizeof(std::uint64_t),"material sampling requires binary64 keys");
              for(unsigned i=0;i<4;++i)std::memcpy(&key[i+1],&values[i],sizeof(double));
              materialOrder.push_back({key,ptr});
            }
            std::sort(materialOrder.begin(),materialOrder.end(),[](const auto& a,const auto& b){return a.first<b.first;});
          }
          std::size_t materialAt=0;
          ptr=FieldLinesAll[iFieldLine].ExternallyManagedGeometry?
              (materialOrder.empty()?-1:materialOrder.front().second):Segment->FirstParticleIndex;

          while (ptr!=-1) {
            spec=PIC::ParticleBuffer::GetI(ptr);

            LocalParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
            LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

            PIC::FieldLine::Sampling(ptr,LocalParticleWeight,NULL);
            ptr=FieldLinesAll[iFieldLine].ExternallyManagedGeometry?
                (++materialAt<materialOrder.size()?materialOrder[materialAt].second:-1):PIC::ParticleBuffer::GetNext(ptr);
          }
        }
      }
    }

    namespace {
      // ---------------------------------------------------------------------
      // Defensive helpers for field-line vertex sampling.
      //
      // The field-line output samples particle quantities to vertices and later
      // reduces those sampled quantities over MPI.  A single NaN in a sampled
      // diagnostic such as DatumAtVertexParticleEnergy is sticky: subsequent
      // additions keep the value NaN and the failure may only be observed much
      // later in PIC::FieldLine::Parallel::PackAllFieldLinesVertexData().
      //
      // Sampling is diagnostic and must not repair or alter the particle state.
      // Therefore the policy here is: validate every quantity used to form a
      // sampled value; skip the particle for this diagnostic sample if the state
      // is outside the safe numerical range; and cap only the diagnostic speed
      // used in the relativistic energy conversion so that Speed2E() is never
      // evaluated at v >= c.  The mover remains responsible for preventing bad
      // particles from being produced in the first place.
      // ---------------------------------------------------------------------
      inline bool IsSafeFieldLineSamplingValue(double v) {
        const double MaxAbsSamplingValue = 1.0e300;
        return std::isfinite(v) && (std::fabs(v) < MaxAbsSamplingValue);
      }

      inline double ClampFieldLineSamplingPitchAngleCosine(double mu) {
        if (mu >  1.0) return  1.0;
        if (mu < -1.0) return -1.0;
        return mu;
      }

      inline bool ShouldPrintFieldLineSamplingWarning() {
        // Keep diagnostics useful without flooding a parallel run.  The skipped
        // samples are diagnostic-only events; a small number of warnings is
        // enough to identify the failure mode while keeping stdout/stderr usable.
        static int nWarningsPrinted = 0;
        if (nWarningsPrinted < 25) {
          nWarningsPrinted++;
          return true;
        }
        return false;
      }

      inline void PrintFieldLineSamplingWarning(const char* reason,long int ptr,int spec,
          int iFieldLine,double S,double Weight,double vParallel,double vNormal,
          double Speed,double extra0=0.0,double extra1=0.0) {
        if ((_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) && ShouldPrintFieldLineSamplingWarning()) {
          std::cerr << "Warning: skipping field-line vertex particle sample: "
                    << reason
                    << "; rank=" << PIC::ThisThread
                    << "; ptr=" << ptr
                    << "; spec=" << spec
                    << "; field_line=" << iFieldLine
                    << "; s=" << S
                    << "; weight=" << Weight
                    << "; vParallel=" << vParallel
                    << "; vNormal=" << vNormal
                    << "; speed=" << Speed
                    << "; extra0=" << extra0
                    << "; extra1=" << extra1
                    << std::endl;
        }
      }
    }

    struct MaterialSamplingSink {
      cFieldLineVertex* vertex;
      cFieldLineSegment* segment;
      bool managed;
      int endpoint;
      void SampleDatum(cDatumSampled& datum,double value,int species,double shape) {
        vertex->SampleDatum(datum,value,species,shape);
        if(managed)segment->SampleMaterialDatum(datum,value,species,endpoint,shape);
      }
      void SampleDatum(cDatumSampled& datum,double value,int species,int component,double shape) {
        vertex->SampleDatum(datum,value,species,component,shape);
        if(managed)segment->SampleMaterialDatum(datum,value,species,endpoint,shape,component);
      }
    };

    void Sampling(long int ptr, double Weight, char* CellSamplingBuffer){
      // namespace alias
      namespace PB = PIC::ParticleBuffer;

      // find the field line and location along it
      int iFieldLine = PB::GetFieldLineId(ptr);
      double TimeStep,S = PB::GetFieldLineCoord(ptr);

      // Validate the field-line coordinate before it is used to access the
      // field-line container.  Invalid field-line coordinates can occur for
      // particles that have just left a field line or for particles whose mover
      // state was corrupted earlier.  Sampling is diagnostic, so the safe action
      // is to skip this particle instead of letting it contaminate the sampled
      // vertex arrays.
      if ((iFieldLine<0)||(iFieldLine>=nFieldLine)||
          (IsSafeFieldLineSamplingValue(S)==false)) {
        PrintFieldLineSamplingWarning("invalid field-line coordinate",ptr,-1,
            iFieldLine,S,Weight,0.0,0.0,0.0);
        return;
      }

      cFieldLineSegment* CurrentSegment=FieldLinesAll[iFieldLine].GetSegment(S);
      if (CurrentSegment==NULL) {
        PrintFieldLineSamplingWarning("null segment for field-line coordinate",ptr,-1,
            iFieldLine,S,Weight,0.0,0.0,0.0);
        return;
      }

      // weights of vertices.  The interpolation weight is a diagnostic weight
      // between the begin and end vertices of the containing segment.  Roundoff
      // near a segment boundary can put the fractional part marginally outside
      // [0,1]; clamp that small numerical overshoot so a valid particle does not
      // generate a negative sampled contribution.
      double w = S - floor(S);
      if (w<0.0) w=0.0;
      if (w>1.0) w=1.0;

      //magnetic field
      double B[3],l[3];

      FieldLinesAll[iFieldLine].GetMagneticField(B,S);
      CurrentSegment->GetDir(l);

      double AbsB = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2], 0.5);
      if ((IsSafeFieldLineSamplingValue(B[0])==false)||
          (IsSafeFieldLineSamplingValue(B[1])==false)||
          (IsSafeFieldLineSamplingValue(B[2])==false)||
          (IsSafeFieldLineSamplingValue(l[0])==false)||
          (IsSafeFieldLineSamplingValue(l[1])==false)||
          (IsSafeFieldLineSamplingValue(l[2])==false)||
          (IsSafeFieldLineSamplingValue(AbsB)==false)) {
        PrintFieldLineSamplingWarning("invalid magnetic field or field-line direction",ptr,-1,
            iFieldLine,S,Weight,0.0,0.0,0.0,AbsB);
        return;
      }

      //magnetic moment, mass, velocity and energy
      double mu = PB::GetMagneticMoment(ptr);
      int spec = PB::GetI(ptr);
      double v[3],Speed;
      double x[3];

      if ((spec<0)||(spec>=PIC::nTotalSpecies)) {
        PrintFieldLineSamplingWarning("invalid species index",ptr,spec,
            iFieldLine,S,Weight,0.0,0.0,0.0);
        return;
      }

      if ((IsSafeFieldLineSamplingValue(Weight)==false)||(Weight<=0.0)) {
        PrintFieldLineSamplingWarning("invalid particle statistical weight",ptr,spec,
            iFieldLine,S,Weight,0.0,0.0,0.0);
        return;
      }

      //PB::GetV(v, ptr);

      //time step
      switch (_SIMULATION_TIME_STEP_MODE_) {
      case _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_:
        TimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
	break;
      case _SINGLE_GLOBAL_TIME_STEP_:
        TimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
        break;
      default:
	exit(__LINE__,__FILE__,"Error: not implemented");
      }

      if (IsSafeFieldLineSamplingValue(TimeStep)==false) {
        PrintFieldLineSamplingWarning("invalid sampling time step",ptr,spec,
            iFieldLine,S,Weight,0.0,0.0,0.0,TimeStep);
        return;
      }

      v[0]=PB::GetVParallel(ptr);
      v[1]=PB::GetVNormal(ptr);
      v[2]=0.0;

      if ((IsSafeFieldLineSamplingValue(v[0])==false)||
          (IsSafeFieldLineSamplingValue(v[1])==false)) {
        PrintFieldLineSamplingWarning("invalid particle velocity components",ptr,spec,
            iFieldLine,S,Weight,v[0],v[1],0.0);
        return;
      }

      Speed=Vector3D::Length(v);

      // Keep the sampler out of the singular point of relativistic energy
      // conversion.  This cap is used only for diagnostic conversion to energy;
      // it does not change the stored particle velocity.  The mover has its own
      // responsibility to keep particles subluminal.
      const double MinFieldLineSamplingSpeed=1.0;
      const double MaxFieldLineSamplingSpeed=(1.0-1.0e-12)*SpeedOfLight;
      if ((IsSafeFieldLineSamplingValue(Speed)==false)||(Speed<MinFieldLineSamplingSpeed)) {
        PrintFieldLineSamplingWarning("invalid or too small particle speed",ptr,spec,
            iFieldLine,S,Weight,v[0],v[1],Speed);
        return;
      }

      double DiagnosticSpeed=Speed;
      if (DiagnosticSpeed>MaxFieldLineSamplingSpeed) DiagnosticSpeed=MaxFieldLineSamplingSpeed;

      double CosPitchAngleSign=(v[0]*Vector3D::DotProduct(B,l)>0.0) ? 1.0 : -1.0;

      PB::GetX(x, ptr);

      double m0= PIC::MolecularData::GetMass(spec);
      if ((IsSafeFieldLineSamplingValue(m0)==false)||(m0<=0.0)) {
        PrintFieldLineSamplingWarning("invalid species mass",ptr,spec,
            iFieldLine,S,Weight,v[0],v[1],Speed,m0);
        return;
      }

      double E=Relativistic::Speed2E(DiagnosticSpeed,m0);
      if ((IsSafeFieldLineSamplingValue(E)==false)||(E<0.0)) {
        PrintFieldLineSamplingWarning("invalid diagnostic kinetic energy",ptr,spec,
            iFieldLine,S,Weight,v[0],v[1],Speed,E);
        return;
      }

      double CosPitchAngle=ClampFieldLineSamplingPitchAngleCosine(v[0]/Speed);

      // volume
      double volume = 0.0;

      if (SegmentVolume!=NULL) {
        volume=SegmentVolume(CurrentSegment,iFieldLine);
      }
      else {
        exit(__LINE__,__FILE__,"Error: hook SegmentVolume is not initialized; provide a function for calculating the volumne of a segment of a field line (e.g., SEP::FieldLine::GetSegmentVolume in srcSEP/shock_injection.cpp");
      }

      if ((IsSafeFieldLineSamplingValue(volume)==false)||(volume<=0.0)) {
        PrintFieldLineSamplingWarning("invalid field-line segment volume",ptr,spec,
            iFieldLine,S,Weight,v[0],v[1],Speed,volume);
        return;
      }


      //sample to vertices
      cFieldLineVertex* V=CurrentSegment->GetBegin();
      cFieldLineVertex* Vnext=(V!=NULL) ? V->GetNext() : NULL;
      if ((V==NULL)||(Vnext==NULL)) {
        PrintFieldLineSamplingWarning("invalid vertex connectivity for sampling",ptr,spec,
            iFieldLine,S,Weight,v[0],v[1],Speed);
        return;
      }

      const double SafeFluxSpeed=Speed;
      MaterialSamplingSink sample{V,CurrentSegment,
          FieldLinesAll[iFieldLine].ExternallyManagedGeometry,0};

      // The actual calls below use the guarded cFieldLineVertex::SampleDatum()
      // implementation.  The per-contribution guard in cFieldLineVertex protects
      // every sampled datum from NaN/Inf, while the checks above make sure the
      // quantities used here are physically and numerically well-defined.
      sample.SampleDatum(DatumAtVertexNumberDensity,Weight/volume, spec, (1-w));
      sample.SampleDatum(DatumAtVertexParticleWeight,Weight,spec, (1-w));
      sample.SampleDatum(DatumAtVertexParticleNumber,1.0,spec, (1-w));
      sample.SampleDatum(DatumAtVertexParticleEnergy,Weight*E*J2MeV,spec, (1-w));
      sample.SampleDatum(DatumAtVertexParticleSpeed,Weight*Speed,spec, (1-w));

      //sample binned number density
      int ibin=DatumAtVertexNumberDensityEnergyBinned.GetBinIndex(E*J2MeV);
      if (ibin>=0) sample.SampleDatum(DatumAtVertexNumberDensityEnergyBinned,Weight/volume, spec, ibin,(1-w));

      sample.SampleDatum(DatumAtVertexParticleCosPitchAngle,Weight*CosPitchAngle,spec, (1-w));
      sample.SampleDatum(DatumAtVertexParticleAbsCosPitchAngle,Weight*fabs(CosPitchAngle),spec, (1-w));

      if (CosPitchAngleSign>0.0) {
        sample.SampleDatum(DatumAtVertexNumberDensity_mu_positive,Weight/volume, spec, (1-w));
        sample.SampleDatum(DatumAtVertexParticleFlux_mu_positive,Weight*SafeFluxSpeed/volume, spec, (1-w));
      }
      else {
        sample.SampleDatum(DatumAtVertexNumberDensity_mu_negative,Weight/volume, spec, (1-w));
        sample.SampleDatum(DatumAtVertexParticleFlux_mu_negative,Weight*SafeFluxSpeed/volume, spec, (1-w));
      }

      if (DatumAtVertexFluence.is_active()==true) {
        double t;

        V->GetDatum(DatumAtVertexFluence,t);
        if (IsSafeFieldLineSamplingValue(t)==false) t=0.0;
        t+=Weight/volume*SafeFluxSpeed*(1-w)*TimeStep;
        if (IsSafeFieldLineSamplingValue(t)==true) V->SetDatum(DatumAtVertexFluence,t);
      }


      V = Vnext;sample.vertex=Vnext;sample.endpoint=1;
      sample.SampleDatum(DatumAtVertexNumberDensity,Weight/volume,spec,  (w));
      sample.SampleDatum(DatumAtVertexParticleWeight,Weight,spec, (w));
      sample.SampleDatum(DatumAtVertexParticleNumber,1.0,spec, (w));
      sample.SampleDatum(DatumAtVertexParticleEnergy,Weight*E*J2MeV,spec, w);
      sample.SampleDatum(DatumAtVertexParticleSpeed,Weight*Speed,spec, w);

      if (ibin>=0) sample.SampleDatum(DatumAtVertexNumberDensityEnergyBinned,Weight/volume, spec, ibin,w);

      sample.SampleDatum(DatumAtVertexParticleCosPitchAngle,Weight*CosPitchAngle,spec, w);
      sample.SampleDatum(DatumAtVertexParticleAbsCosPitchAngle,Weight*fabs(CosPitchAngle),spec, w);

      if (CosPitchAngleSign>0.0) {
        sample.SampleDatum(DatumAtVertexNumberDensity_mu_positive,Weight/volume, spec, (w));
        sample.SampleDatum(DatumAtVertexParticleFlux_mu_positive,Weight*SafeFluxSpeed/volume, spec, (w));
      }
      else {
        sample.SampleDatum(DatumAtVertexNumberDensity_mu_negative,Weight/volume, spec, (w));
        sample.SampleDatum(DatumAtVertexParticleFlux_mu_negative,Weight*SafeFluxSpeed/volume, spec, (w));
      }

      if (DatumAtVertexFluence.is_active()==true) {
        double t;

        V->GetDatum(DatumAtVertexFluence,t);
        if (IsSafeFieldLineSamplingValue(t)==false) t=0.0;
        t+=Weight/volume*SafeFluxSpeed*(w)*TimeStep;
        if (IsSafeFieldLineSamplingValue(t)==true) V->SetDatum(DatumAtVertexFluence,t);
      }


      //.........................
      if ((CellSamplingBuffer!=NULL)&&(DatumAtGridParticleEnergy.offset>=0)) {
        double* grid_energy=(double*)(CellSamplingBuffer+DatumAtGridParticleEnergy.offset);
        if (IsSafeFieldLineSamplingValue(*grid_energy)==false) *grid_energy=0.0;
        const double dEGrid=Weight*E;
        if (IsSafeFieldLineSamplingValue(dEGrid)==true) {
          const double updated=*grid_energy+dEGrid;
          if (IsSafeFieldLineSamplingValue(updated)==true) *grid_energy=updated;
        }
      }
    }

    //=========================================================================
    void FieldLineWeight_Uniform(double* Weight, cFieldLineSegment* Segment) {
      for(int spec=0; spec<PIC::nTotalSpecies; spec++) Weight[spec] = 1.0;
    }

    //=========================================================================
    void InitLoop2D(double *xStart,  //start loop here
		    double DArc,     //increment in the angle of arc 
		    double DMin,     //min allowed length of the arc
		    double DMax      //max allowed length of the arc
		    ) {
      // in 2D case (e.g. cylindrical symmetry) generate a magnetic field line;
      // magnetic flux function, psi, is utilized:
      // loop is a line of constant value
      //-----------------------------------------------------------------------
      // check correctness
      if (_PIC_SYMMETRY_MODE_ != _PIC_SYMMETRY_MODE__AXIAL_) {
      exit(__LINE__, __FILE__,"ERROR: implemented only for axial symmetry!");
      }

      // check for limit of number of field lines allowed
      if (nFieldLine == nFieldLineMax) exit(__LINE__,__FILE__,"ERROR: reached limit for field line number");

      // mark time of extraction
      TimeLastUpdate = PIC::SimulationTime::Get();

      // list of variables utilized
      //.......................................................................
      // current location
      double x[3] = {xStart[0], xStart[1], xStart[2]};

      // candidate location
      //      double xNew[3] = {0.0, 0.0, 0.0};
      // the first location; used to close the loop
      double xFirst[3] = {xStart[0], xStart[1], xStart[2]};

      // min value resolved for magnetic field (squared)
      const double epsB2 = 1E-30;
      const double epsB  = 1E-15;

      // plasma velocity, stored at vertex
      double V[3] = {0.0, 0.0, 0.0};

      //plasma velocity
      double T = 0;

      // magnetic field vector, unit vector and magnitude
      // NOTE: y-component is ignored!!!
      double B[3] = {0.0, 0.0, 0.0};
      double b[3] = {0.0, 0.0, 0.0};
      double absB =  0.0;

      // flux function
      double Psi0 = 0.0, Psi = 0.0;

      //housekeeping for controlling loop's generation
      // Arc = curvature * Length
      double Arc = 0.0;

      // direction of new segment
      double Dir[3] = {0.0, 0.0, 0.0};

      // angle swiped by the loop so far
      double Angle = 0.0;

      // estimate for the next segment's length
      double Length = 0.5*(DMin + DMax);

      // position and magnetic field and next candidate vertex
      double xNew[3] = {0.0,0.0,0.0};
      double BNew[3] = {0.0,0.0,0.0};
      double bNew[3] = {0.0,0.0,0.0};
      double absBNew;

      // new vertex
      cFieldLineVertex *Vertex;
      //.......................................................................

      //increase counter of field lines
      nFieldLine++;

      //rotate to the y=0 plane
      x[0] = pow(x[0]*x[0] + x[1]*x[1], 0.5) * ( (x[0]>0) ? 1 : -1);
      x[1] = 0.0;

      // get magnetic field at the current location
      CPLR::InitInterpolationStencil(x);
      CPLR::GetBackgroundMagneticField(B);
      absB = pow(B[0]*B[0]+B[2]*B[2],0.5);
      b[0] = B[0]/absB; b[1] = 0; b[2] = B[2]/absB;

      // target value of the flux function
      Psi0 = CPLR::GetBackgroundMagneticFluxFunction();

      // check correctness: need to have non-zero components in x-z plane
      if (B[0]*B[0] + B[2]*B[2] < epsB2) exit(__LINE__,__FILE__, "ERROR: magnetic field magnitude is below min resolved value");
      
      //add the initial vertex
      Vertex = FieldLinesAll[nFieldLine-1].Add(x);

      //      FieldLinesAll[nFieldLine-1].SetMagneticField(B,0);
      Vertex -> SetMagneticField(B);
      Vertex -> SetDatum(DatumAtVertexMagneticFluxFunction, Psi0);
      CPLR::GetBackgroundPlasmaVelocity(V);
      Vertex -> SetDatum(DatumAtVertexPlasmaVelocity, V);
      T = CPLR::GetBackgroundPlasmaTemperature();
      Vertex -> SetDatum(DatumAtVertexPlasmaTemperature, T);	

      //generate the loop
      //......................................................................
      // new vertices are added until:
      // - 2*Pi angle is swiped AND 
      //   candidate vertex is spatially close to the first one => SUCCESS
      // - 4*Pi angle is swiped => FAILURE
      // - exited the domain => FAILURE
      //......................................................................
      // new point is generated in 2 steps
      // 1) make the first guess, x_0, in the direction of B
      // 2) iterate towards point with same value of FLUX via gradient descent
      //......................................................................
      
      while (fabs(Angle) < 4*Pi) {
        // get the first guess
        //------------------------------------------------------------------
        // control inner while loop
        int count;
        const  int countMax = 100;
        count = 0;

        //predictor
        bool DoBreak = false;

        while(true) {
          count++;

          if(count > countMax) exit(__LINE__, __FILE__,"ERROR: can't generate field-line loop ");

          // get new location
          xNew[0] = x[0] + Length * b[0];
          xNew[2] = x[2] + Length * b[2];

          // get magnetic field at this location
          CPLR::InitInterpolationStencil(xNew);
          CPLR::GetBackgroundMagneticField(BNew);
          absBNew = pow(BNew[0]*BNew[0]+BNew[2]*BNew[2],0.5);

          // need to have non-zero field
          if (absBNew < epsB) exit(__LINE__,__FILE__, "ERROR: magnetic field magnitude is below min resolved");

          bNew[0] = BNew[0]/absBNew; bNew[1] = 0; bNew[2] = BNew[2]/absBNew;

          // find angle of arc to the new location
          Arc = fabs(b[0]*bNew[2] - b[2]*bNew[0]);
          if (DoBreak) break;

          //check condition to exit loop
          if ((Arc < 0.5 * DArc) || (Arc > DArc)) {
            // too small or too large step; factor -> 1 as count grows
            // to avoid reaching stationary points and infinite loop
            Length *= 1. + (0.75 * DArc / max(Arc, 0.075*DArc) - 1.) / count;
          }
          else break;

          if (Length < DMin) {Length = DMin; DoBreak = true;}
          if (Length > DMax) {Length = DMax; DoBreak = true;}
        }

        // corrector
        {
          Dir[0] = 0.5*(BNew[0]+B[0]);
          Dir[1] = 0;
          Dir[2] = 0.5*(BNew[2]+B[2]);

          double misc = pow(Dir[0]*Dir[0] + Dir[2]*Dir[2], 0.5);
          Dir[0]/= misc; Dir[2]/= misc;

          //new vertex location
          xNew[0] = x[0] +  Length * Dir[0];
          xNew[2] = x[2] +  Length * Dir[2];
        }
        //------------------------------------------------------------------

        // iterate towards location with needed value of flux function
        //------------------------------------------------------------------
        //get magnetic field at the candidate location
        CPLR::InitInterpolationStencil(xNew);
        CPLR::GetBackgroundMagneticField(BNew);

        // need to have non-zero field
        if(BNew[0]*BNew[0]+BNew[2]*BNew[2] < epsB2) exit(__LINE__,__FILE__,"ERROR: magnetic field magnitude is below min resolved");

        // flux function
        Psi = CPLR::GetBackgroundMagneticFluxFunction();

        // correct the initial guess until flux function is close enough
        // to the original value
        //------------------------------------------------------------------
        while(fabs(2*(Psi-Psi0)/(Psi+Psi0)) > 1E-5) {
          double misc = BNew[0]*BNew[0] + BNew[2]*BNew[2];

          xNew[0]+=-0.001 * (Psi-Psi0) * BNew[2] / misc;
          xNew[2]+= 0.001 * (Psi-Psi0) * BNew[0] / misc;
          CPLR::InitInterpolationStencil(xNew);
          CPLR::GetBackgroundMagneticField(BNew);
          Psi = CPLR::GetBackgroundMagneticFluxFunction();
        }
        //------------------------------------------------------------------

        // next vertex is found;
        // add it to the field line and update information
        absBNew = pow(BNew[0]*BNew[0]+BNew[2]*BNew[2],0.5);
        bNew[0] = BNew[0]/absBNew; bNew[1] = 0; bNew[2] = BNew[2]/absBNew;

        //housekeeping
        Angle += asin(b[0]*bNew[2] - b[2]*bNew[0]);

        //add this vertex
        Vertex = FieldLinesAll[nFieldLine-1].Add(xNew);

        //	FieldLinesAll[nFieldLine-1].SetMagneticField(BNew);
        Vertex -> SetMagneticField(BNew);
        Vertex -> SetDatum(DatumAtVertexMagneticFluxFunction, Psi0);
        CPLR::GetBackgroundPlasmaVelocity(V);
        Vertex -> SetDatum(DatumAtVertexPlasmaVelocity, V);
        T = CPLR::GetBackgroundPlasmaTemperature();
        Vertex -> SetDatum(DatumAtVertexPlasmaTemperature, T);

        //check conditions for completing the loop
        cFieldLineSegment* First=FieldLinesAll[nFieldLine-1].GetFirstSegment();
        cFieldLineSegment* Last =FieldLinesAll[nFieldLine-1].GetLastSegment();

        double Dist = 0.0; //distance between last and first verticies
        double DirFirst[3], DirLast[3];

        First->GetDir(DirFirst); Last->GetDir(DirLast);

        for (int idim=0; idim<DIM; idim++) {
          Dist  += pow(xNew[idim]-xFirst[idim], 2);
        }

        Dist   = pow(Dist, 0.5);

        if ((fabs( 0.5*fabs(Angle)/Pi - 1) < 0.1) && (Dist < Length)) {
          FieldLinesAll[nFieldLine-1].close_loop();
          return;
        }

        //save current location and magnetic field
        x[0] = xNew[0]; x[1] = xNew[1]; x[2] = xNew[2];
        b[0] = bNew[0]; b[1] = bNew[1]; b[2] = bNew[2];
        B[0] = BNew[0]; B[1] = BNew[1]; B[2] = BNew[2];
        absB = absBNew;
      }

      exit(__LINE__, __FILE__,"ERROR: can't generate field-line loop ");
    }

    //=========================================================================
    // update field lines 
    void Update(){
      cFieldLineVertex*  Vertex;
      cFieldLineSegment* Segment;
      double dt = PIC::SimulationTime::Get() - TimeLastUpdate;

      TimeLastUpdate += dt;

      for (int iFieldLine=0; iFieldLine<nFieldLine; iFieldLine++) {
        if(FieldLinesAll[iFieldLine].ExternallyManagedGeometry)continue;
        int cnt = 0;
        bool done = false;

        for (Vertex = FieldLinesAll[iFieldLine].GetFirstVertex();true;Vertex = Vertex->GetNext(), cnt++) {
          double X[3],B[3],V[3];

          done = ((Vertex==FieldLinesAll[iFieldLine].GetLastVertex()) && (cnt!=0)) ? true : false;
          if (done) break;

          Vertex->GetX(X);
          Vertex->GetPlasmaVelocity(V);

          for (int i=0; i<DIM; i++) X[i] += V[i] * dt;

          if (_PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_) {
          // rotate to the y=0 plane
          X[0] = pow(X[0]*X[0]+X[1]*X[1], 0.5);
          X[1] = 0.0;
          } //if (_PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_)

          if (_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_) {
           // flux function
          CPLR::InitInterpolationStencil(X);
          CPLR::GetBackgroundMagneticField(B);

          double Psi = CPLR::GetBackgroundMagneticFluxFunction();
          double Psi0;
          Vertex->GetDatum(DatumAtVertexMagneticFluxFunction, Psi0);

          // correct the initial guess until flux function is close enough
          // to the original value
          //------------------------------------------------------------------
          while(fabs(2*(Psi-Psi0)/(Psi+Psi0)) > 1E-5) {
            double misc = B[0]*B[0] + B[2]*B[2];

            X[0]+=-0.001 * (Psi-Psi0) * B[2] / misc;
            X[2]+= 0.001 * (Psi-Psi0) * B[0] / misc;

            CPLR::InitInterpolationStencil(X);
            CPLR::GetBackgroundMagneticField(B);
            Psi = CPLR::GetBackgroundMagneticFluxFunction();
          }
	        //------------------------------------------------------------------
         }//_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_

          Vertex->SetX(X);

          // update background data
          PIC::CPLR::InitInterpolationStencil(X);
          PIC::CPLR::GetBackgroundMagneticField(B);
          PIC::CPLR::GetBackgroundPlasmaVelocity(V);
          Vertex->SetDatum(DatumAtVertexMagneticField, B);
          Vertex->SetDatum(DatumAtVertexPlasmaVelocity,V);

          double T = CPLR::GetBackgroundPlasmaTemperature();
          Vertex -> SetDatum(DatumAtVertexPlasmaTemperature, T);
	     }

	      cnt = 0; done = false;

	      for (Segment = FieldLinesAll[iFieldLine].GetFirstSegment();!done; Segment = Segment->GetNext(), cnt++) {
	        done = Segment==FieldLinesAll[iFieldLine].GetLastSegment()&&cnt!=0;
	        Segment->SetVertices(Segment->GetBegin(),Segment->GetEnd());
        }
      }
    }

    //=========================================================================
    void InitSimpleParkerSpiral(double *xStart) {
      
if (DIM != 3) {
      exit(__LINE__,__FILE__,"Implemetned only for 3D case");
}

      if(nFieldLine == nFieldLineMax) exit(__LINE__,__FILE__,"ERROR: reached limit for field line number");

      //increase counter of field lines
      nFieldLine++;
      
      //spatial step and # of segments
      const double Ds=1E9;
      int nSegment = 1000;
           
      //Parker spiral starts @ 10 Sun radii
      const double R0 = 10 *_SUN__RADIUS_;

      //magnetic field at this location
      const double B0 = 1.83E-6;

      //velocity of solar wind
      const double SpeedSolarWind = 4.0E5;

      //rate of solar rotation at equator (rad/sec)
      const double Omega = 2.97211E-6;      

      //position of the Sun
      // Field-line geometry is embedded in physical 3-space independently of
      // the AMR mesh dimension.  These vectors are explicitly indexed through
      // component 2 below, so sizing them by DIM corrupts the 1-D/2-D compile
      // contract even when this optional initializer is never selected.
      double xSun[3]={0.0,0.0,0.0};

      //relative position to Sun
      double x2Sun[3]={0.0,0.0,0.0};

      //distance to Sun
      double R2Sun=0;

      //distance to sun in equatorial plane
      double r2Sun=0;

      //find original latitude
      for(int idim=0; idim<DIM; idim++) x2Sun[idim] = xStart[idim] - xSun[idim];
      R2Sun = sqrt(x2Sun[0]*x2Sun[0]+x2Sun[1]*x2Sun[1]+x2Sun[2]*x2Sun[2]);

      //original latitude 
      double costheta = x2Sun[2] / R2Sun;
      double sintheta = sqrt(1 - costheta*costheta); // >=0 always     

      //start position
      double x[3]={xStart[0],xStart[1],xStart[2]};


      //if too close to the Sun => move to distance R0
      if (R2Sun < R0) for (int idim=0; idim<DIM; idim++) x[idim] = xSun[idim] + x2Sun[idim] * R0/R2Sun;

      for (int iSegment=0;iSegment<nSegment; iSegment++) {
        //distances to sun
        for(int idim=0; idim<DIM; idim++) x2Sun[idim]=x[idim]-xSun[idim];

        R2Sun = sqrt(x2Sun[0]*x2Sun[0]+x2Sun[1]*x2Sun[1]+x2Sun[2]*x2Sun[2]);
        r2Sun = sqrt(x2Sun[0]*x2Sun[0]+x2Sun[1]*x2Sun[1]);

        //angle in the equatorial plane
        double cosphi = x2Sun[0]/r2Sun, sinphi=x2Sun[1]/r2Sun;

        double misc1 = B0 * pow(R0/R2Sun, 2);
        double misc2 = (R2Sun - R0) * Omega/SpeedSolarWind * sintheta;
      
        //magnetic field at location x
        double B[3];

        B[0] = misc1 * (sintheta*cosphi + misc2 * sinphi);
        B[1] = misc1 * (sintheta*sinphi - misc2 * cosphi);
        B[2] = misc1 *  costheta;
      
        //add vertex to the spiral
        FieldLinesAll[nFieldLine-1].Add(x);
        FieldLinesAll[nFieldLine-1].SetMagneticField(B,iSegment);
      
        //magnitude of the magnetic field
        double AbsB = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);

        //direction of the magnetic field
        double b[3] = {B[0]/AbsB, B[1]/AbsB, B[2]/AbsB};

        // next vertex of the spiral
        x[0]+=b[0]*Ds; x[1]+=b[1]*Ds; x[2]+=b[2]*Ds;
      }
    }
  }

  //---------------------------------------------------------------------------

namespace Mover{
namespace FieldLine{
  // A segment-owned particle belongs to the FIELD-LINE decomposition, which
  // need not match the Cartesian AMR block containing its derived position.
  // Its complete physical interval [s] is therefore the configured global
  // species interval, exactly as in the serial segment sweep. Looking up a
  // local Cartesian block can fail on another MPI owner, or select an unrelated
  // spatial timestep. Share this selection across all segment worker paths;
  // worker count must never change the sampling interval or transport physics.
  // Spatially local timesteps have no defined segment-owned contract here and
  // remain an explicit unsupported configuration, never a global fallback.
  static double SegmentParticleTimeStep(int species) {
    if(species<0||species>=PIC::nTotalSpecies)
      exit(__LINE__,__FILE__,"invalid species in field-line timestep selection");
    double dt=0.0;
    #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    dt=PIC::ParticleWeightTimeStep::GlobalTimeStep[species];
    #elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
    dt=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
    #else
    exit(__LINE__,__FILE__,"segment-owned field-line motion requires a global physical timestep");
    #endif
    if(!std::isfinite(dt)||dt<=0.0)
      exit(__LINE__,__FILE__,"invalid global physical timestep for field-line motion");
    return dt;
  }

  //manager executing the particle moving procedure when the partice lists are attached to the field line segments
  void MoveParticlesMultiThread(int this_thread_id,int thread_id_table_size) {
    namespace FL = PIC::FieldLine;
    namespace PB = PIC::ParticleBuffer;

    long int SegmentCounter=0;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL;
    long int ptr,ptr_next;
    double LocalTimeStep;
    int spec;
    PB::byte* ParticleData;
    FL::cFieldLineSegment* Segment; 

    for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
      for (int i=0;i<FL::FieldLinesAll[iFieldLine].GetTotalSegmentNumber();i++) if (++SegmentCounter%thread_id_table_size==this_thread_id) {
        Segment=FL::FieldLinesAll[iFieldLine].SegmentPointerTable[i];

        ptr=Segment->FirstParticleIndex;

        while (ptr!=-1) {
          ptr_next=PB::GetNext(ptr);
 
          ParticleData=PB::GetParticleDataPointer(ptr);

          spec=PB::GetI(ParticleData);

          // This manager is entered only for segment-attached particles.
          // Pass no Eulerian node, just as the serial field-line manager does.
          // The mover consumes the immutable line snapshot, not a mesh block.
          LocalTimeStep=SegmentParticleTimeStep(spec);
          PIC::Mover::DispatchParticleMover(ptr,LocalTimeStep,node);

          ptr=ptr_next;
        }
      }
    }
  }


  void MoveParticles () {
    namespace FL = PIC::FieldLine;
    namespace PB = PIC::ParticleBuffer; 
  
    if ((_PIC_PARTICLE_LIST_ATTACHING_!=_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_)||(_PIC_FIELD_LINE_MODE_!=_PIC_MODE_ON_)) {
      exit(__LINE__,__FILE__,"Error: wrong options");
    } 

    if (_PIC_FL_MOVER_MULTITHREAD_== _PIC_MODE_ON_) {
      std::thread tTable[_PIC_NUMBER_STD_THREADS_];

      for (int i=1;i<_PIC_NUMBER_STD_THREADS_;i++) {
        tTable[i]=std::thread(MoveParticlesMultiThread,i,_PIC_NUMBER_STD_THREADS_);
      }

      MoveParticlesMultiThread(0,_PIC_NUMBER_STD_THREADS_);

      for (int i=1;i<_PIC_NUMBER_STD_THREADS_;i++)  tTable[i].join();

      //update particle list indexes 
      for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
        for (auto Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
          Segment->FirstParticleIndex=Segment->tempFirstParticleIndex;
          Segment->tempFirstParticleIndex=-1;
        }
      }

      return;
    }


const int _process_by_particles=0;
const int _process_by_segments=1;

int _process_mode=_process_by_segments;


  auto CountParticle = [&] (FL::cFieldLineSegment* Segment) {
    long int ptr;
    int cnt=0;

    ptr=Segment->FirstParticleIndex;

    while (ptr!=-1) {
      cnt++;
      ptr=PB::GetNext(ptr);
    }

    return cnt;
  };

  auto PopulateParticleTable = [&] (long int *ParticleTable,FL::cFieldLineSegment* Segment) {
    long int ptr;
    int cnt=0;

    ptr=Segment->FirstParticleIndex;

    while (ptr!=-1) {
      ParticleTable[cnt]=ptr;

      cnt++;
      ptr=PB::GetNext(ptr);
    }

    return cnt;
  };
      

  auto ProcessSegment = [&] (FL::cFieldLineSegment* Segment) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL;
    long int ptr,ptr_next;
    double *x,LocalTimeStep;
    int spec;
    PB::byte* ParticleData;

    ptr=Segment->FirstParticleIndex;


    long int *ParticleTable=NULL;
    int ParticleTableLength=0;

    int SegmentParticleNumber=CountParticle(Segment);

    if (ParticleTableLength<SegmentParticleNumber) {
      if (ParticleTable!=NULL) delete [] ParticleTable;

      ParticleTable=new long int [SegmentParticleNumber]; 
      ParticleTableLength=SegmentParticleNumber;
    }

    PopulateParticleTable(ParticleTable,Segment);

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    #pragma omp parallel for default (none) private (ptr,x,spec, ParticleData,LocalTimeStep) firstprivate(node,SegmentParticleNumber,PIC::Mesh::mesh,ParticleTable) 
    #endif
    for (int ii=0;ii<SegmentParticleNumber;ii++) {
      ptr=ParticleTable[ii];
      ParticleData=PB::GetParticleDataPointer(ptr);

      x=PB::GetX(ParticleData);
      spec=PB::GetI(ParticleData);

      switch (_PIC_PARTICLE_LIST_ATTACHING_) {
      case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
        node=NULL;
        LocalTimeStep=SegmentParticleTimeStep(spec);
        break;
      case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
        node=PIC::Mesh::mesh->findTreeNode(x,node);
        if (node==NULL) exit(__LINE__,__FILE__,"Error: the point is not found");
        if (node->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not allocated");
	LocalTimeStep=node->block->GetLocalTimeStep(spec);
	break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }

      PIC::Mover::DispatchParticleMover(ptr,LocalTimeStep,node);
    }


    if (ParticleTable!=NULL) delete [] ParticleTable;
  };

    auto ProcessEntireSegment = [&] (FL::cFieldLineSegment* Segment) {
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL;
      long int ptr,ptr_next;
      double *x,LocalTimeStep;
      int spec;
      PB::byte* ParticleData;

      ptr=Segment->FirstParticleIndex;

      while (ptr!=-1) {
        ptr_next=PB::GetNext(ptr);

       ParticleData=PB::GetParticleDataPointer(ptr);

       x=PB::GetX(ParticleData);
       spec=PB::GetI(ParticleData);

      switch (_PIC_PARTICLE_LIST_ATTACHING_) {
      case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
        node=NULL;

        LocalTimeStep=SegmentParticleTimeStep(spec);

        break;
      case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
        node=PIC::Mesh::mesh->findTreeNode(x,node);
        if (node==NULL) exit(__LINE__,__FILE__,"Error: the point is not found");
        if (node->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not allocated");
        LocalTimeStep=node->block->GetLocalTimeStep(spec);

        break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }

       PIC::Mover::DispatchParticleMover(ptr,LocalTimeStep,node);

       ptr=ptr_next;
    }
   };


    //loop through the field lines 
    int iFieldLine;
    FL::cFieldLineSegment *Segment; 
    int thread=0,SegmentCounter=0;

    
    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    thread=omp_get_thread_num();
    #endif

    for (iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
      switch (_process_mode) {
      case _process_by_segments:
        #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
        #pragma omp parallel for schedule(dynamic) private (Segment) firstprivate(FL::FieldLinesAll,iFieldLine) 
        #endif
        for (int i=0;i<FL::FieldLinesAll[iFieldLine].GetTotalSegmentNumber();i++) {
          Segment=FL::FieldLinesAll[iFieldLine].SegmentPointerTable[i]; 
          ProcessEntireSegment(Segment);
        }

        break;
      case _process_by_particles:
        for (Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
          ProcessSegment(Segment);
          SegmentCounter++;
        }

        break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }
    } 
 

    //update particle list indexes 
    for (iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
      for (Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
        Segment->FirstParticleIndex=Segment->tempFirstParticleIndex;
        Segment->tempFirstParticleIndex=-1;                
      }
    }
  }

  // procedure that returns parameters of the guiding center motion
  void GuidingCenterMotion(double& ForceParal, double& AbsB,int spec,long int ptr,int   iFieldLine,double FieldLineCoord) {
    /* function returns guiding center velocity in direction perpendicular
     * to the magnetic field and the force parallel to it
     * for Lorentz force (considered here)
     *
     * v_{guide_perp} = 
     *  E\times B / B^2 + 
     *  \mu/(q\gamma) B\times\nabla B / B^2 +
     *  (p_{\parallel}^2)/(q\gamma m_0) B\times(B\cdot\nabla)B/B^4
     *
     * dp_{\parallel}/dt=
     *  q E_{\parallel} - \mu/\gamma \nabla_{\parallel}B
     *
     * \mu = p_{\perp}^2/(2m_0B)
     *
     * \gamma = 1/\sqrt{1-v^2/c^2}
     *********************************************************************/

    // namespace alias
    namespace FL = PIC::FieldLine;

if (_PIC__IDEAL_MHD_MODE_ == _PIC_MODE_OFF_) {
    exit(__LINE__,__FILE__,"Error: calculation of the electric field component of the Lorentz force is not implemented");
}

    double B[3],B0[3],B1[3], AbsBDeriv;
    double mu     = PIC::ParticleBuffer::GetMagneticMoment(ptr);
    double q      = PIC::MolecularData::GetElectricCharge(spec);

    FL::FieldLinesAll[iFieldLine].GetMagneticField(B0, (int)FieldLineCoord);
    FL::FieldLinesAll[iFieldLine].GetMagneticField(B,       FieldLineCoord);
    FL::FieldLinesAll[iFieldLine].GetMagneticField(B1, (int)FieldLineCoord+1-1E-7);
    AbsB   = pow(B[0]*B[0] + B[1]*B[1] + B[2]*B[2], 0.5);
    
    AbsBDeriv = (pow(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2], 0.5) -
        pow(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2], 0.5)) /  FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

    //parallel force
switch(_PIC__IDEAL_MHD_MODE_) {
    case _PIC_MODE_ON_:
        // in this case E = - V \cross B => E_{\paral} = E*b = 0
        ForceParal = - mu * AbsBDeriv;
        break;
    default:
        exit(__LINE__, __FILE__, "not implemented");
}
  }
  
  // mover itself
  int Mover_SecondOrder(long int ptr, double dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    //aliases
    namespace FL = PIC::FieldLine;
    namespace PB = PIC::ParticleBuffer;

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
    double dtTemp;
    PIC::ParticleBuffer::byte *ParticleData;
    double AbsBInit=0.0;
//    double vInit[3]={0.0,0.0,0.0},xInit[3]={0.0,0.0,0.0};
    double AbsBMiddle=0.0;
    double vFinal[3]={0.0,0.0,0.0},xFinal[3]={0.0,0.0,0.0};

    int   iFieldLine = -1;
    double FieldLineCoordInit   = -1.0;
    double FieldLineCoordMiddle = -1.0;
    double FieldLineCoordFinal  = -1.0;
    
    double ForceParalInit=0.0, ForceParalMiddle=0.0;
    double dirInit[3]={0.0,0.0,0.0};
    double vParInit=0.0, vParMiddle=0.0, vParFinal=0.0;
    double pParInit,pPerpInit,pParMiddle,pParFinal;
    int i,j,k,spec;
	
    ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    pPerpInit=PIC::ParticleBuffer::GetMomentumNormal(ParticleData);
    pParInit=PIC::ParticleBuffer::GetMomentumParallel(ParticleData);


    spec=PIC::ParticleBuffer::GetI(ParticleData);

    double m0 = PIC::MolecularData::GetMass(spec);
    double mu = PIC::ParticleBuffer::GetMagneticMoment(ptr);
    double v[3];
    double p[3];

    iFieldLine = PIC::ParticleBuffer::GetFieldLineId(ptr);
    FieldLineCoordInit = PIC::ParticleBuffer::GetFieldLineCoord(ptr);

    static long int nCall=0;
    nCall++;
    
    // predictor step
if (_PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_) {
    GuidingCenterMotion(ForceParalInit,AbsBInit,spec,ptr,iFieldLine,FieldLineCoordInit);
}
    
    FL::FieldLinesAll[iFieldLine].GetSegmentDirection(dirInit, FieldLineCoordInit);

    switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
    case _PIC_MODE_OFF_:
      vParInit=pParInit/m0;
      break;
    case _PIC_MODE_ON_:
      p[0]=pParInit,p[1]=pPerpInit,p[2]=0.0; 

      ::Relativistic::Momentum2Vel(v,p,m0);
      vParInit=v[0];
    }
    
    dtTemp=dtTotal/2.0;

    // advance coordinates half-step
    FieldLineCoordMiddle = FL::FieldLinesAll[iFieldLine].move(FieldLineCoordInit, dtTemp*vParInit);

    // advance momentum half-step
    pParMiddle=pParInit+dtTemp*ForceParalInit;
    
    if (FL::FieldLinesAll[iFieldLine].is_loop()) FL::FieldLinesAll[iFieldLine].fix_coord(FieldLineCoordMiddle);
    
    // check if a particle has left the domain
    if (FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoordMiddle)==NULL) { 
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;
      
      //call the function to processes particles that left the domain
      switch(code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }
    
    
    // corrector step
if (_PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_) {
    GuidingCenterMotion(ForceParalMiddle,AbsBMiddle,spec,ptr,iFieldLine,FieldLineCoordMiddle);
}
    

    // advance coordinates full-step

    switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
    case _PIC_MODE_OFF_:
      vParMiddle=pParMiddle/m0;
      break;
    case _PIC_MODE_ON_:
      p[0]=pParMiddle,p[1]=pPerpInit,p[2]=0.0;
      ::Relativistic::Momentum2Vel(v,p,m0);
      vParMiddle=v[0];
    }

    FieldLineCoordFinal = FL::FieldLinesAll[iFieldLine].move(FieldLineCoordInit, dtTotal*vParMiddle);

    // advance momentum full-step
    pParFinal=pParInit+dtTotal*ForceParalMiddle;

    if (FL::FieldLinesAll[iFieldLine].is_loop()) FL::FieldLinesAll[iFieldLine].fix_coord(FieldLineCoordFinal);

    FL::cFieldLineSegment *Segment;
   
    //advance the particle's position and velocity
    //interaction with the faces of the block and internal surfaces
    if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoordFinal))==NULL) {
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;
      
      //call the function that process particles that leaved the coputational domain
      switch (code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }

    FL::FieldLinesAll[iFieldLine].GetCartesian(xFinal, FieldLineCoordFinal);
    newNode=PIC::Mesh::Search::FindBlock(xFinal);


    ////PIC::ParticleBuffer::GetMomentumNormal(pPerpFinal,ParticleData) --> pPerpFinal was not updated 
    PIC::ParticleBuffer::SetMomentumParallel(pParFinal,ParticleData);


    double l[3],e0[3],e1[3],SpeedPar,SpeedPerp;

    Segment->GetDir(l);
    Vector3D::GetNormFrame(e0,e1,l);

    double vParallel,vNormal;

    switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
    case _PIC_MODE_OFF_:
//      vFinal[0]=pParFinal/m0,vFinal[1]=pPerpInit/m0,vFinal[2]=0.0;

      vParallel=pParFinal/m0;
      vNormal=pPerpInit/m0;

      for (int idim=0;idim<3;idim++) vFinal[idim]=SpeedPar*l[idim]+SpeedPerp*e0[idim]; 
      break;
    case _PIC_MODE_ON_:
      p[0]=pParFinal,p[1]=pPerpInit,p[2]=0.0;

      //for (int idim=0;idim<3;idim++) p[idim]=pParFinal*l[idim]+pPerpInit*e0[idim]; 
      ::Relativistic::Momentum2Vel(vFinal,p,m0);

      vParallel=vFinal[0];
      vNormal=vFinal[1];
    }

    //PIC::ParticleBuffer::SetV(vFinal,ParticleData);
    PIC::ParticleBuffer::SetVParallel(vParallel,ParticleData);
    PIC::ParticleBuffer::SetVNormal(vNormal,ParticleData);


    PIC::ParticleBuffer::SetX(xFinal,ParticleData);
    PB::SetFieldLineCoord(FieldLineCoordFinal, ParticleData);
  
    //adjust the value of 'startNode'
    startNode=newNode;

  //save the trajectory point
if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
  PIC::ParticleTracker::RecordTrajectoryPoint(xFinal,vFinal,spec,ParticleData,(void*)startNode);
}


if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
  if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_) {
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)startNode);
  }
}


    if (PIC::Mesh::mesh->FindCellIndex(xFinal,i,j,k,startNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell where the particle is located");
    
    switch (_PIC_PARTICLE_LIST_ATTACHING_) {
    case  _PIC_PARTICLE_LIST_ATTACHING_NODE_:   
      #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_ 
        { 
        PIC::Mesh::cDataBlockAMR *block=startNode->block;
        long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    
        PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
        PIC::ParticleBuffer::SetPrev(-1,ParticleData);
    
        if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
        block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;
        }
      #else  
        exit(__LINE__,__FILE__,"Error: the block need to be generalized to the case when _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_");
      #endif

      break;
    case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
       PIC::ParticleBuffer::SetNext(Segment->tempFirstParticleIndex,ParticleData);
       PIC::ParticleBuffer::SetPrev(-1,ParticleData);

       if (Segment->tempFirstParticleIndex!=-1) PIC::ParticleBuffer::SetPrev(ptr,Segment->tempFirstParticleIndex);
       Segment->tempFirstParticleIndex=ptr;
       break;
     default:
       exit(__LINE__,__FILE__,"Error: the option is unknown");
     }
 
    return _PARTICLE_MOTION_FINISHED_;
  }

}
}
  //namespace Mover -----------------------------------------------------------

}


//check consistency of the particles lists in case the lists are attached to the filed line segments  
void PIC::FieldLine::CheckParticleList() {
  namespace FL = PIC::FieldLine;
  namespace PB = PIC::ParticleBuffer;

  if ((_PIC_PARTICLE_LIST_ATTACHING_!=_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_)||(_PIC_FIELD_LINE_MODE_!=_PIC_MODE_ON_)) {
    exit(__LINE__,__FILE__,"Error: wrong options");
  }

  auto ProcessSegment = [&] (FL::cFieldLineSegment* Segment) {
    long int ptr;
    int cnt=0;

    ptr=Segment->FirstParticleIndex;

    while (ptr!=-1) {
      cnt++;
      ptr=PB::GetNext(ptr);
    }

    return cnt;
  };
 
  //loop through the field lines 
  int iFieldLine;
  FL::cFieldLineSegment *Segment;
  long int nTotalParticles=0;

  for (iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    for (Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
      nTotalParticles+=ProcessSegment(Segment);
    }
  }

  if (nTotalParticles!=PB::GetAllPartNum()) exit(__LINE__,__FILE__,"Error: the total number of particles stored in the lists is different from that stored in the particle buffer");
} 

void PIC::FieldLine::PopulateSegment(int spec,double NumberDensity,double Temperature,double* BulkVelocity,double Volume,int iSegment,int iFieldLine,int nMaxInjectedParticles) {
  namespace FL = PIC::FieldLine;
  namespace PB = PIC::ParticleBuffer;

  double l[3],anpart,ParticleWeight=0,x[3][3],WeightCorrectionFactor=1.0;
  int npart,cnt=0;
  cFieldLineSegment* Segment=FieldLinesAll[iFieldLine].GetSegment(iSegment);

  //determine the injected particle number 
  Segment->GetBegin()->GetX(x[0]);  
  Segment->GetEnd()->GetX(x[1]); 

  for (int idim=0;idim<3;idim++) {
    x[2][idim]=0.5*(x[0][idim]+x[1][idim]); 
    l[idim]=x[1][idim]-x[0][idim];
  }    

  for (int ip=0;ip<3;ip++) {
    //determine the local particle weight at the point x[ip]
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->findTreeNode(x[ip],NULL);    

    if (node!=NULL) {
      cnt++;
      ParticleWeight+=node->block->GetLocalParticleWeight(spec);
    }
  }

  if (cnt==0) exit(__LINE__,__FILE__,"Error: no segment point is within the domain");

  ParticleWeight/=cnt;      
  anpart=NumberDensity*Volume/ParticleWeight; 

  if ((nMaxInjectedParticles!=-1)&&(anpart>nMaxInjectedParticles)) {
    //limit the number of injected particles 
    WeightCorrectionFactor=anpart/nMaxInjectedParticles;
    anpart=nMaxInjectedParticles;  

    if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_) exit(__LINE__,__FILE__,"Error: _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ must be set _INDIVIDUAL_PARTICLE_WEIGHT_ON_"); 
  }

  npart=(int)anpart;
  if (anpart-npart>rnd()) npart++;


  //populate new particles 
  for (;npart>0;npart--) {
    //get particle location 
    double xp[3],vp[3],vParallel,vNormal,S=rnd();
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
    
    for (int idim=0;idim<3;idim++) xp[idim]=x[0][idim]+S*l[idim]; 
 
    node=PIC::Mesh::mesh->findTreeNode(xp,NULL);     

    if (node==NULL) continue; 
    else if (node->Thread!=PIC::ThisThread) continue; 

    //get porticle's microscopic properties 
    PIC::Distribution::MaxwellianVelocityDistribution(vp,BulkVelocity,Temperature,spec);
    Vector3D::GetComponents(vParallel,vNormal,vp,l);
    
    //add the particle to the particle list   
    // pointer particle to the particle to be injected
    PB::byte ptrData[PB::ParticleDataLength];

    // A fresh source template follows the same zero-payload contract as
    // GetNewParticle. Cloning uninitialized extension bytes would otherwise
    // undo the allocator's initialization before an application birth hook.
    memset(ptrData,0,PB::ParticleDataLength);

    //default settings
    PB::SetIndividualStatWeightCorrection(WeightCorrectionFactor,ptrData);
    PB::SetI(spec, ptrData);

    //set field line
    PB::SetFieldLineId(iFieldLine, ptrData);

    //set the localtion and velocity 
    S+=iSegment;
    PB::SetFieldLineCoord(S, ptrData);
    PB::SetX(xp, ptrData);

    PB::SetVParallel(vParallel,ptrData);
    PB::SetVNormal(vNormal,ptrData);


    long int NewParticle;
    PB::byte* NewParticleData;

    switch (_PIC_PARTICLE_LIST_ATTACHING_) {
    case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
      NewParticle=PB::GetNewParticle(Segment->FirstParticleIndex);
      NewParticleData=PB::GetParticleDataPointer(NewParticle);

      PB::CloneParticle(NewParticleData,ptrData);
      PB::SetParticleAllocated(NewParticleData);
      break;
    case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
      NewParticle=PB::InitiateParticle(xp,NULL,NULL,NULL,ptrData,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
  } 
}

PIC::FieldLine::MaxwellianLineBeginningPreflight
PIC::FieldLine::PreflightMaxwellianLineBeginning(int spec,int iFieldLine) {
  MaxwellianLineBeginningPreflight result;
  using Code=MaxwellianLineBeginningPreflight::Code;
  if(spec<0||spec>=PIC::nTotalSpecies||!FieldLinesAll||iFieldLine<0||iFieldLine>=nFieldLine||
      !FieldLinesAll[iFieldLine].GetFirstSegment()||!FieldLinesAll[iFieldLine].GetFirstVertex()) {
    result.message="Maxwellian line source has invalid species or geometry";return result;
  }
  auto* segment=FieldLinesAll[iFieldLine].GetFirstSegment();
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=nullptr;
  if(_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_) {
    result.ownerRank=segment->Thread;
    if(result.ownerRank<0||result.ownerRank>=PIC::nTotalThreads) {
      result.message="Maxwellian line source has an invalid segment owner";return result;
    }
    if(result.ownerRank!=PIC::ThisThread){result.code=Code::NonOwner;return result;}
  }
  const bool localTime=_SIMULATION_TIME_STEP_MODE_==_SPECIES_DEPENDENT_LOCAL_TIME_STEP_;
  const bool localWeight=_SIMULATION_PARTICLE_WEIGHT_MODE_==_SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_;
  if(localTime||localWeight||_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_NODE_) {
    if(!PIC::Mesh::mesh) {
      result.code=Code::UnsupportedConfiguration;result.message="Maxwellian line source needs an initialized local AMR owner";return result;
    }
    double x[3];FieldLinesAll[iFieldLine].GetFirstVertex()->GetX(x);
    node=PIC::Mesh::mesh->findTreeNode(x,nullptr);
    if(!node){result.message="Maxwellian line source endpoint is outside the AMR domain";return result;}
    if(_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_NODE_) {
      result.ownerRank=node->Thread;
      if(result.ownerRank!=PIC::ThisThread){result.code=Code::NonOwner;return result;}
    }
    if(node->Thread!=PIC::ThisThread||!node->block) {
      result.code=Code::UnsupportedConfiguration;
      result.message="Maxwellian field-line source local timestep/weight block is absent on the segment owner; remote lookup is unspecified";
      return result;
    }
  }
  if(_PIC_PARTICLE_LIST_ATTACHING_!=_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_&&
      _PIC_PARTICLE_LIST_ATTACHING_!=_PIC_PARTICLE_LIST_ATTACHING_NODE_) {
    result.code=Code::UnsupportedConfiguration;result.message="unsupported Maxwellian source particle attachment";return result;
  }
  switch(_SIMULATION_TIME_STEP_MODE_) {
    case _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_:result.timeStepS=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];break;
    case _SINGLE_GLOBAL_TIME_STEP_:result.timeStepS=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];break;
    case _SPECIES_DEPENDENT_LOCAL_TIME_STEP_:result.timeStepS=node->block->GetLocalTimeStep(spec);break;
    default:result.code=Code::UnsupportedConfiguration;result.message="unsupported Maxwellian source timestep mode";return result;
  }
  switch(_SIMULATION_PARTICLE_WEIGHT_MODE_) {
    case _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_:result.particleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];break;
    case _SINGLE_GLOBAL_PARTICLE_WEIGHT_:result.particleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[0];break;
    case _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_:result.particleWeight=node->block->GetLocalParticleWeight(spec);break;
    default:result.code=Code::UnsupportedConfiguration;result.message="unsupported Maxwellian source weight mode";return result;
  }
  if(!std::isfinite(result.timeStepS)||result.timeStepS<=0.0||
      !std::isfinite(result.particleWeight)||result.particleWeight<=0.0) {
    result.message="Maxwellian line source requires finite positive timestep and particle weight";return result;
  }
  result.code=Code::Ready;return result;
}

//inject particles at the beginning of a segment
int PIC::FieldLine::InjectMaxwellianLineBeginning(int spec,double NumberDensity,double Temperature,double* BulkVelocity,double InjectionArea,int iFieldLine,int nMaxInjectedParticles,
    void (*postFreshParticle)(long int,std::uint64_t,void*),void* sourceContext) {
  namespace FL = PIC::FieldLine;
  namespace PB = PIC::ParticleBuffer; 
  const auto preflight=PreflightMaxwellianLineBeginning(spec,iFieldLine);
  if(preflight.code==MaxwellianLineBeginningPreflight::Code::NonOwner)return 0;
  if(preflight.code!=MaxwellianLineBeginningPreflight::Code::Ready)
    exit(__LINE__,__FILE__,preflight.message.c_str());
  
  int nInjectedParticles=0;
  int npart;
  double dtTotal=preflight.timeStepS,SourceRate,n_sw,t_sw,v_v[3],l[3],anpart,
      ParticleWeight=preflight.particleWeight,x[3],GlobalWeightCorrectionFactor=1.0;
  auto Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();
  FL::cFieldLineVertex* FirstVertex=Segment->GetBegin();
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;
  
  //get the direction of the first segment of the magnetic field line 
  Segment->GetDir(l);

  if (Vector3D::DotProduct(l,BulkVelocity)>0.0) { 
    for (int idim=0;idim<3;idim++) l[idim]*=-1.0;
  }

  //determine the corresponding source rate 
  SourceRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(NumberDensity,Temperature,BulkVelocity,l,spec);

  //determine the number of the injected model particles 
  // Cartesian output belongs to the newborn, never to the provider vertex.
  // The former pointer alias let GetCartesian(x,S) silently move the first
  // physical vertex at every injection without rebuilding its segment metric.
  FirstVertex->GetX(x);
  if(_PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_NODE_)
    node=PIC::Mesh::mesh->findTreeNode(x,nullptr);
  anpart=dtTotal*SourceRate/ParticleWeight*InjectionArea;

  if ((nMaxInjectedParticles!=-1)&&(anpart>nMaxInjectedParticles)) {
    //limit the number of injected particles 
    GlobalWeightCorrectionFactor=anpart/nMaxInjectedParticles;
    anpart=nMaxInjectedParticles;

    if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_) exit(__LINE__,__FILE__,"Error: _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ must be set _INDIVIDUAL_PARTICLE_WEIGHT_ON_");
  }
  
  npart=(int)anpart;
  if (anpart-npart>rnd()) npart++;
  
  nInjectedParticles=npart;

  //populate new particles 
  std::uint64_t sourceMacro=0;
  for (;npart>0;npart--,++sourceMacro) {
    //get particle location 
    double vp[3],vParallel,vNormal,WeightCorrectionFactor,S=0.0;

    //get porticle's microscopic properties         
    WeightCorrectionFactor=PIC::Distribution::InjectMaxwellianDistribution(vp,BulkVelocity,Temperature,l,spec,_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_);
    Vector3D::GetComponents(vParallel,vNormal,vp,l);
    vParallel*=-1.0;
    
    //add the particle to the particle list   
    // pointer particle to the particle to be injected
    PB::byte ptrData[PB::ParticleDataLength];

    // Preserve the complete fresh-record initialization contract, including
    // application extensions not interpreted by this generic source helper.
    memset(ptrData,0,PB::ParticleDataLength);

    //default settings
    PB::SetIndividualStatWeightCorrection(WeightCorrectionFactor*GlobalWeightCorrectionFactor,ptrData);
    PB::SetI(spec, ptrData);

    //set field line
    PB::SetFieldLineId(iFieldLine, ptrData);
    
    //move particle along the field line for a fraction of a time step
    S=FL::FieldLinesAll[iFieldLine].move(S,dtTotal*rnd()*vParallel);
    FL::FieldLinesAll[iFieldLine].GetCartesian(x,S);

    //set the localtion and velocity 
    PB::SetFieldLineCoord(S,ptrData);
    PB::SetX(x, ptrData);

    PB::SetVParallel(vParallel,ptrData);
    PB::SetVNormal(vNormal,ptrData);

    long int NewParticle;
    PB::byte* NewParticleData;

    switch (_PIC_PARTICLE_LIST_ATTACHING_) {
    case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
      NewParticle=PB::GetNewParticle(Segment->FirstParticleIndex);
      NewParticleData=PB::GetParticleDataPointer(NewParticle);

      PB::CloneParticle(NewParticleData,ptrData);
      PB::SetParticleAllocated(NewParticleData);
      break;
    case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
      NewParticle=PB::InitiateParticle(x,NULL,NULL,NULL,ptrData,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    if(postFreshParticle&&NewParticle!=-1)postFreshParticle(NewParticle,sourceMacro,sourceContext);
  } 
  
  return nInjectedParticles;
}
  
long int  PIC::FieldLine::TraverseAllFieldLines(void (*processParticle)(long int)) {
    namespace FL = PIC::FieldLine;
    long int cnt=0;

    if (FL::FieldLinesAll == nullptr) {
        if (PIC::ThisThread == 0) printf("TraverseAllFieldLines: ERROR -> FieldLinesAll is NULL\n");
        return 0;
    }

    // Loop through all field lines
    for (int iFieldLine = 0; iFieldLine < FL::nFieldLine; iFieldLine++) {
        FL::cFieldLine& currentFieldLine = FL::FieldLinesAll[iFieldLine];
        FL::cFieldLineSegment* currentSegment = currentFieldLine.GetFirstSegment();

        // Loop through all segments in current field line
        while (currentSegment != nullptr) {
            // Only process segments assigned to this thread
            if (currentSegment->Thread == PIC::ThisThread) {
                // Get particles in current segment
                long int particlePtr = currentSegment->FirstParticleIndex;

                // Loop through all particles in current segment
                while (particlePtr != -1) {
                    // Process the particle using the provided function
                    processParticle(particlePtr);
                    cnt++;

                    // Move to next particle
                    particlePtr = PIC::ParticleBuffer::GetNext(particlePtr);
                }
            }

            currentSegment = currentSegment->GetNext();
        }
    }

    return cnt;
}
  
