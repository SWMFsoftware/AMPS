#include "pic.h"
#include "util/sep_shock_source_core.h"
#include <iostream>

#if _PIC_FIELD_LINE_MODE_ != _PIC_MODE_ON_ || _PIC_PARTICLE_LIST_ATTACHING_ != _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_
#error "Native material fixture requires a configured field-line/segment-attached build"
#endif
#if _CUDA_MODE_ == _ON_
#error "Native material geometry fixture currently qualifies constructed CPU storage only"
#endif

// A linked native-container test, deliberately independent of a solar-wind
// model. The owner callbacks prescribe exact affine motion/entry data; the
// real cFieldLine, pool allocator and particle buffer perform the transaction.
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  std::cout<<"NATIVE_CONTAINER_STARTED\n"<<std::flush;
  namespace FL=PIC::FieldLine;namespace PB=PIC::ParticleBuffer;
  namespace G=FL::Geometry;
  MPI_GLOBAL_COMMUNICATOR=MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD,&PIC::ThisThread);
  MPI_Comm_size(MPI_COMM_WORLD,&PIC::nTotalThreads);
  PB::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;
  FL::VertexAllocationManager.MagneticField=true;
  FL::VertexAllocationManager.PlasmaVelocity=true;
  FL::VertexAllocationManager.Fluence=true;
  FL::Init();PB::Init(128);
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[0]=2;
  PIC::MolecularData::SetMass(1.67262192369e-27,0);
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0]=0.75;
  // An explicitly prescribed 2 m^3 test volume makes the fluence oracle
  // distinguish density*speed*dt from the former erroneous extra 1/volume.
  FL::SegmentVolume=[](FL::cFieldLineSegment*,int){return 2.0;};
  auto& line=FL::FieldLinesAll[0];FL::nFieldLine=1;
  for(double r:{10.0,14.0,18.0,22.0}) {
    double x[3]={r,0,0},b[3]={1e-9,0,0},u[3]={5/0.75,0,0};
    auto* vertex=line.Add(x);vertex->SetMagneticField(b);vertex->SetPlasmaVelocity(u);
  }
  if(!line.RegisterGeometryOwner("controlled-material-line-0",0).ok())return 1;
  auto particle=[&](int segment,double s) {
    const long int p=PB::GetNewParticle();PB::SetI(0,p);
    PB::SetFieldLineCoord(s,p);PB::SetFieldLineId(0,p);
    PB::SetIndividualStatWeightCorrection(3,p);PB::SetNext(-1,p);
    PB::SetVParallel(2,p);PB::SetVNormal(3,p);
    line.GetSegment(segment)->FirstParticleIndex=p;return p;
  };
  const long int kept=particle(0,0.25);particle(2,2.5);
  FL::Sampling();
  double fluence=0;line.GetFirstVertex()->GetDatum(FL::DatumAtVertexFluence,fluence);
  bool sampled=line.GetSegment(1)->SamplingExposure.collecting==1&&
      std::fabs(fluence-(6.0/2)*std::sqrt(13.0)*0.75*0.75)<1e-12;
  FL::cFieldLine::GeometryPlan plan;
  G::Context c;c.arrivalS=0.75;c.minimumRadiusM=10;c.maximumRadiusM=20;
  c.toleranceM=1e-10;c.maximumInsertions=10;
  const G::Motion move=[](std::size_t,const G::Point& x,double departure,double arrival,G::Point* out) {
    *out=x;(*out)[0]+=(5/0.75)*(arrival-departure);return G::Status{};
  };
  const G::InnerTrace trace=[](const G::Point& x,const G::Point&,double length,double,G::Point* out) {
    *out=x;(*out)[0]-=length;return G::Status{};
  };
  const FL::cFieldLine::VertexStateInitializer state=[&](int index,const G::Point&,double,std::vector<char>* out) {
    out->assign(FL::cFieldLineVertex::StoredDataBytes(),0);
    // Existing integral diagnostics belong to the surviving material vertex.
    // Only a genuinely new entry record starts with a zero integral history.
    if(index>=0)std::memcpy(out->data(),line.GetVertex(index)->GetAssociatedDataBufferPointer(),out->size());
    const double b[3]={1e-9,0,0};
    std::memcpy(out->data()+FL::DatumAtVertexMagneticField.offset,b,sizeof(b));
    const double u[3]={5/0.75,0,0};
    std::memcpy(out->data()+FL::DatumAtVertexPlasmaVelocity.offset,u,sizeof(u));
    return G::Status{};
  };
  // No wave/edge physics is excited in this fixture. Any allocated edge datum
  // is an explicitly prescribed zero state, not a production fallback.
  const FL::cFieldLine::SegmentStateInitializer edgeState=[](int,const G::Point&,const G::Point&,double,std::vector<char>* out) {
    out->assign(FL::cFieldLineSegment::StoredDataBytes(),0);return G::Status{};
  };
  FL::cFieldLine::GeometryPlan rejectedPlan;
  const FL::cFieldLine::VertexStateInitializer rejectedState=[&](int index,const G::Point& x,double epoch,std::vector<char>* out) {
    return index==2?G::Status::Error("controlled late retained vertex rejection"):state(index,x,epoch,out);
  };
  sampled=sampled&&!line.PrepareGeometryStep(c,move,trace,rejectedState,&rejectedPlan,edgeState).ok()&&
      line.GetFirstVertex()->GetX()[0]==10&&line.GetTotalSegmentNumber()==3&&
      PB::GetFieldLineCoord(kept)==0.25&&line.GeometryGeneration==0;
  auto prepared=line.PrepareGeometryStep(c,move,trace,state,&plan,edgeState);
  bool pass=prepared.ok()&&sampled;
  if(!prepared.ok())std::cerr<<prepared.message<<'\n';
  pass=pass&&line.GetFirstVertex()->GetX()[0]==10&&line.GetTotalSegmentNumber()==3;
  if(!pass) {
    std::cerr<<"FAIL native material preparation or independent fluence oracle\n";
    line.Delete();FL::nFieldLine=0;MPI_Finalize();return 1;
  }
  auto committed=line.CommitGeometryStep(&plan);
  if(!committed.ok())std::cerr<<committed.message<<'\n';
  if(!committed.ok()){line.Delete();FL::nFieldLine=0;MPI_Finalize();return 1;}
  double x[3];PB::GetX(x,kept);
  double retainedFluence=0;line.GetVertex(2)->GetDatum(FL::DatumAtVertexFluence,retainedFluence);
  pass=pass&&committed.ok()&&line.GetFirstVertex()->GetX()[0]==7&&
      line.GetLastVertex()->GetX()[0]==23&&line.GetTotalSegmentNumber()==4&&
      line.GetSegment(2)->FirstParticleIndex==kept&&PB::GetFieldLineCoord(kept)==2.25&&
      x[0]==16&&retainedFluence==fluence&&plan.Losses()[0].macros==1&&plan.Losses()[0].represented==6&&
      plan.Losses()[0].kineticJ>0&&plan.Losses()[0].parallelMomentumKgMPerS>0&&
      plan.RetiredSamples().size()==1&&plan.RetiredSamples()[0].exposure.collecting==1&&
      !plan.RetiredSamples()[0].collecting.empty()&&
      line.GetSegment(0)->SamplingExposure.collecting==0&&
      line.GetSegment(2)->SamplingExposure.collecting==1&&
      PB::GetVParallel(kept)==2&&PB::GetVNormal(kept)==3&&
      PB::GetIndividualStatWeightCorrection(kept)==3&&
      line.FirstMaterialLabel==-2&&line.GeometryGeneration==1&&
      !line.CommitGeometryStep(&plan).ok();
  auto* table=PIC::ParallelFieldLines::ThreadSegmentTable;
  pass=pass&&table&&table[0].TableLength[0]==4;
  if(table)for(int i=0;i<4;++i)pass=pass&&table[0].Table[0][i]==line.GetSegment(i);
  // The actual moved native geometry must intersect the independently timed
  // heliocentric shock. The material-relative travel is (7-5)=2 m; after two
  // prefix insertions this is S=2+2/4=2.5, at Cartesian radius 17 m.
  SEP::Shock::TrajectoryConfiguration trajectory;trajectory.launchRadiusM=10;
  SEP::Shock::SpeedKnot knot;knot.radiusM=0;knot.speedMPerS=7/0.75;
  trajectory.knots.push_back(knot);SEP::Shock::TrajectoryState shock;
  pass=pass&&SEP::Shock::StateAtEpoch(trajectory,c.arrivalS,&shock).ok();
  std::vector<SEP::Shock::Vector3> shockLine;
  for(int i=0;i<=line.GetTotalSegmentNumber();++i) {
    const double* point=line.GetVertex(i)->GetX();shockLine.push_back({point[0],point[1],point[2]});
  }
  const auto hit=SEP::Shock::IntersectSphere(shockLine,shock.radiusM,1e-10,SEP::Shock::IntersectionPolicy::FirstOutward);
  pass=pass&&hit.status==SEP::Shock::IntersectionStatus::Ok;
  if(hit.status==SEP::Shock::IntersectionStatus::Ok)
    pass=pass&&std::fabs(hit.intersections[0].segment+hit.intersections[0].fraction-2.5)<1e-12&&
        std::fabs(hit.intersections[0].positionM.x-17)<1e-12;
  // Refresh is a supported operation after topology changes; the previous
  // one-shot table constructor rejected precisely this required lifecycle.
  const auto tableGeneration=PIC::ParallelFieldLines::GetThreadSegmentTableGeneration();
  PIC::ParallelFieldLines::GenerateThreadSegmentTable();
  pass=pass&&PIC::ParallelFieldLines::GetThreadSegmentTableGeneration()==tableGeneration+1;
  // Explicit manufactured population after proving that insertion itself
  // created no particles. This gives old/new segments identical density but
  // different exposure and exercises the actual native output normalization.
  particle(0,0.25);
  PIC::SimulationTime::SetInitialValue(c.arrivalS);
  FL::Sampling();
  pass=pass&&line.GetSegment(0)->SamplingExposure.collecting==1&&
      line.GetSegment(2)->SamplingExposure.collecting==2;
  FL::Output("native-material-sampling.dat",false);
  const std::size_t width=FL::cFieldLineVertex::SampleDataBytes()/sizeof(double);
  const std::size_t densityOffset=FL::DatumAtVertexNumberDensity.offset/sizeof(double);
  const auto* fresh=line.GetSegment(0);const auto* old=line.GetSegment(2);
  pass=pass&&fresh->SamplingExposure.completed==1&&old->SamplingExposure.completed==2&&
      std::fabs(fresh->SamplingExposure.Average(fresh->MaterialCompleted.at(densityOffset))-2.25)<1e-12&&
      std::fabs(old->SamplingExposure.Average(old->MaterialCompleted.at(densityOffset))-2.25)<1e-12&&
      fresh->MaterialCompleted.size()==2*width;
  std::cout<<(pass?"PASS":"FAIL")<<" native material insertion/removal/particle-coordinate transaction\n";
  line.Delete();FL::nFieldLine=0;PIC::ParallelFieldLines::GenerateThreadSegmentTable();
  MPI_Finalize();return pass?0:1;
}
