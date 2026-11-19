// Supporting production force/mover integration at a DECLARED ANALYTICAL
// frozen LSO snapshot. RK4 calls the actual force and therefore verifies the
// integrator/dispatch, not independent physical force truth or observations.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonSurfaceEnvironment.h"
#include "species/ReactionExecution.h"

#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
#error Moon force integration fixture requires actual orbital force terms ON
#endif
namespace T=PIC::Mover::ForceCutSurface;
namespace R=Moon::Species::ReactionExecution;
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
using Vec=std::array<long double,3>;
using State=std::array<long double,6>;
namespace {
int Order=0,Checks=0,Failures=0;
long int ClockOffset=-1;
std::ofstream CheckFile,Rows,Terms,ReferenceProof;
const double X0[3]={2.11e6,123456,234567},V0[3]={730,-410,95};
const double Electric[3]={1e-6,-2e-6,3e-6},Magnetic[3]={3e-9,-4e-9,5e-9};
constexpr double Duration=4,Epoch=1000;
void Require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
void Check(bool condition,const char* name) {
  const int local=condition?1:0;int global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  ++Checks;if(!global)++Failures;
  if(PIC::ThisThread==0){CheckFile<<name<<','<<(global?"PASS":"FAIL")<<'\n';CheckFile.flush();}
}
long double Norm(const Vec& a){return std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);}
Vec Cross(const Vec& a,const Vec& b){return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};}
Vec Difference(const State& a,const State& b,int offset){return {{a[offset]-b[offset],a[offset+1]-b[offset+1],a[offset+2]-b[offset+2]}};}
double Resolution(double*){return 1e6;}
double UnitLoad(Node*){return 1;}
void Collect(Node* node,std::vector<Node*>& leaves) {
  if(node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_){leaves.push_back(node);return;}
  for(int k=0;k<8;++k)if(node->downNode[k])Collect(node->downNode[k],leaves);
}
void Setup() {
  PIC::Init_BeforeParser();
  Moon::Runtime::Configuration c;
  c.particleTracingMover=Order==1?"force-cut-surface-first-order":"force-cut-surface-second-order";
  c.simulationEnvironmentMode="frozen";c.tracingFieldFrame="lso";
  c.tracingRelativeTolerance=1e-9;c.tracingAbsolutePositionToleranceM=1e-9;
  c.tracingAbsoluteVelocityToleranceMS=1e-10;c.tracingGeometryToleranceM=1e-8;
  c.tracingTimeToleranceS=1e-12;c.tracingMinimumSubstepS=1e-12;c.tracingMaximumSubstepS=1;
  c.tracingMaximumSubsteps=10000;c.tracingMaximumRefinements=64;c.tracingMaximumInteractions=1000;
  c.tracingMaximumSamePointInteractions=16;c.tracingAdaptive=false;
  for(int k=0;k<3;++k){c.tracingElectricFieldVM[k]=Electric[k];c.tracingMagneticFieldT[k]=Magnetic[k];}
  std::string error;Require(Moon::Runtime::InstallConfiguration(c,&error),error);
  Moon::RegisterParticleTracing();Moon::ConfigureParticleTracingFields();
  PIC::ParticleBuffer::RequestDataStorage(ClockOffset,sizeof(R::OpticalDepthClock));
  double lower[3]={1.8e6,-5e5,-5e5},upper[3]={2.8e6,5e5,5e5};
  PIC::Mesh::mesh->AllowBlockAllocation=false;PIC::Mesh::mesh->init(lower,upper,Resolution);
  Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree),"root split failed");
  for(int k=0;k<8;++k)Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree->downNode[k]),"child split failed");
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::initCellSamplingDataBuffer();PIC::Mesh::mesh->AllowBlockAllocation=true;PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();PIC::Init_AfterParser();PIC::DomainBlockDecomposition::UpdateBlockTable();
  std::vector<Node*> leaves;Collect(PIC::Mesh::mesh->rootTree,leaves);
  for(Node* node:leaves)if(node->block)for(int species=0;species<PIC::nTotalSpecies;++species) {
    node->block->SetLocalParticleWeight(1+species,species);
    node->block->SetLocalTimeStep(.125/(1+species),species);
  }
  PIC::ParticleBuffer::Init(64);
  // Explicit analytical snapshot, not furnished ephemerides. All mutable
  // quantities have SI components in the same frozen Moon-centred LSO frame.
  Exosphere::xObjectRadial=_AU_;
  const double sun[3]={_AU_,0,0},earth[3]={-3.8e8,2e7,5e6},omega[3]={1e-7,-3e-7,2e-6},objectV[3]={0,30000,0};
  for(int k=0;k<3;++k){Exosphere::xSun_SO[k]=sun[k];Exosphere::xEarth_SO[k]=earth[k];
    Exosphere::RotationVector_SO_FROZEN[k]=omega[k];Exosphere::vObject_SO_FROZEN[k]=objectV[k];}
  T::Acceleration=nullptr;
  PIC::Mover::QueryTrajectoryEvent=nullptr;PIC::Mover::AdvanceTrajectoryEvent=nullptr;PIC::Mover::ProcessTrajectoryEvent=nullptr;
}
State Initial(){return {{X0[0],X0[1],X0[2],V0[0],V0[1],V0[2]}};}
State Derivative(const State& s,int species,Node* node) {
  double x[3],v[3],a[3];for(int k=0;k<3;++k){x[k]=static_cast<double>(s[k]);v[k]=static_cast<double>(s[k+3]);}
  Moon::TotalParticleAcceleration(a,species,-1,x,v,node);
  return {{s[3],s[4],s[5],a[0],a[1],a[2]}};
}
State Reference(int species,Node* node,int steps) {
  // Independently coded classical RK4; NEVER installed as the particle mover.
  // Long-double accumulation keeps the lunar-distance offset from erasing a
  // small convergence error while force evaluation remains actual production.
  State s=Initial();const long double h=static_cast<long double>(Duration)/steps;
  for(int step=0;step<steps;++step) {
    const auto a=Derivative(s,species,node);State t;
    for(int k=0;k<6;++k)t[k]=s[k]+h*a[k]/2;
    const auto b=Derivative(t,species,node);
    for(int k=0;k<6;++k)t[k]=s[k]+h*b[k]/2;
    const auto c=Derivative(t,species,node);
    for(int k=0;k<6;++k)t[k]=s[k]+h*c[k];
    const auto d=Derivative(t,species,node);
    for(int k=0;k<6;++k)s[k]+=h*(a[k]+2*b[k]+2*c[k]+d[k])/6;
  }return s;
}
void ForceWiring(int species,Node* node) {
  const Vec x{{X0[0],X0[1],X0[2]}},v{{V0[0],V0[1],V0[2]}};
  const Vec omega{{Exosphere::RotationVector_SO_FROZEN[0],Exosphere::RotationVector_SO_FROZEN[1],Exosphere::RotationVector_SO_FROZEN[2]}};
  const Vec sun{{_AU_,0,0}},earth{{Exosphere::xEarth_SO[0],Exosphere::xEarth_SO[1],Exosphere::xEarth_SO[2]}};
  const Vec omegaX=Cross(omega,x),omega2X=Cross(omega,omegaX),omegaV=Cross(omega,v);
  Vec lunar{},solar{},terrestrial{},centrifugal{},coriolis{},radiation{},lorentz{},sum{};
  const auto differential=[&](const Vec& attractor,double mass) {
    Vec displacement;for(int k=0;k<3;++k)displacement[k]=attractor[k]-x[k];
    const long double r=Norm(attractor),d=Norm(displacement);Vec out;
    for(int k=0;k<3;++k)out[k]=static_cast<long double>(GravityConstant)*mass*(displacement[k]/(d*d*d)-attractor[k]/(r*r*r));
    return out;
  };
  solar=differential(sun,_MASS_(_SUN_));terrestrial=differential(earth,_MASS_(_EARTH_));
  const long double radius=Norm(x);
  for(int k=0;k<3;++k){lunar[k]=-static_cast<long double>(GravityConstant)*_MASS_(_MOON_)*x[k]/(radius*radius*radius);
    centrifugal[k]=-omega2X[k];coriolis[k]=-2*omegaV[k];}
  if(species==_NA_SPEC_) {
    Require(!Moon::EarthShadowCheck(const_cast<double*>(X0)),"declared snapshot unexpectedly Earth-shadowed");
    Vec relativeV,displacement;
    for(int k=0;k<3;++k){relativeV[k]=v[k]+omegaX[k]+Exosphere::vObject_SO_FROZEN[k];displacement[k]=x[k]-sun[k];}
    const long double distance=Norm(displacement);long double radial=0;
    for(int k=0;k<3;++k)radial+=relativeV[k]*displacement[k]/distance;
    const double scalar=SodiumRadiationPressureAcceleration__Combi_1997_icarus(static_cast<double>(radial),static_cast<double>(distance));
    for(int k=0;k<3;++k)radiation[k]=scalar*displacement[k]/distance;
  }
  else {
    const Vec e{{Electric[0],Electric[1],Electric[2]}},b{{Magnetic[0],Magnetic[1],Magnetic[2]}},vxB=Cross(v,b);
    const long double qm=PIC::MolecularData::GetElectricCharge(species)/PIC::MolecularData::GetMass(species);
    for(int k=0;k<3;++k)lorentz[k]=qm*(e[k]+vxB[k]);
  }
  const Vec vectors[]={lunar,solar,terrestrial,centrifugal,coriolis,radiation,lorentz};
  const char* names[]={"lunar","solar-differential","earth-differential","centrifugal","coriolis","radiation","lorentz"};
  bool nonzero=true;
  for(int j=0;j<7;++j) {
    if(j<5 || (j==5&&species==_NA_SPEC_) || (j==6&&species==_NA_PLUS_SPEC_))nonzero=nonzero && Norm(vectors[j])>0;
    if(PIC::ThisThread==0)Terms<<species<<','<<names[j]<<','<<static_cast<double>(vectors[j][0])<<','
      <<static_cast<double>(vectors[j][1])<<','<<static_cast<double>(vectors[j][2])<<'\n';
    for(int k=0;k<3;++k)sum[k]+=vectors[j][k];
  }
  double localProduction[3]={},production[3];
  if(node->Thread==PIC::ThisThread) {
    const auto actual=Derivative(Initial(),species,node);
    for(int k=0;k<3;++k)localProduction[k]=static_cast<double>(actual[k+3]);
  }
  MPI_Allreduce(localProduction,production,3,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  bool agree=true;
  for(int k=0;k<3;++k)agree=agree && std::abs(production[k]-sum[k])<=1e-13L+5e-12L*std::abs(sum[k]);
  Check(nonzero,"all_applicable_force_terms_nonzero");Check(agree,"actual_total_force_initial_wiring");
}
State Advance(int species,double h,Node* node) {
  T::Settings.MaximumSubstepS=h;const bool owner=node->Thread==PIC::ThisThread;
  R::OpticalDepthClock incoming;std::string error;Require(R::OpticalDepthClock::FromOpenUniform(.37,&incoming,&error),error);
  long int particle=-1;double local[6]={};int count=0,invalid=0;
  if(owner) {
    particle=PIC::ParticleBuffer::GetNewParticle();auto* data=PIC::ParticleBuffer::GetParticleDataPointer(particle);
    PIC::ParticleBuffer::SetI(species,data);PIC::ParticleBuffer::SetX(const_cast<double*>(X0),data);PIC::ParticleBuffer::SetV(const_cast<double*>(V0),data);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1/node->block->GetLocalParticleWeight(species),data);
    PIC::ParticleBuffer::SetNext(-1,data);PIC::ParticleBuffer::SetPrev(-1,data);
    std::memcpy(data+ClockOffset,&incoming,sizeof(incoming));T::ResetParticleState(data);T::SetParticleEpochS(particle,Epoch);
    const int code=Moon::MoveParticle(particle,Duration,node);
    if(code!=_PARTICLE_MOTION_FINISHED_)invalid=1;
    const auto time=T::GetTimeContext(particle);
    if(PIC::ParticleBuffer::GetI(particle)!=species || std::memcmp(data+ClockOffset,&incoming,sizeof(incoming)) ||
        std::abs(node->block->GetLocalParticleWeight(species)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle)-1)>1e-14 ||
        T::GetParticleEpochS(particle)!=Epoch+Duration || time.InitialEpochS!=Epoch || time.RequestedTimeS!=Duration ||
        time.ElapsedTimeS!=Duration || time.RemainingTimeS!=0 || time.IntegrationOrder!=Order)invalid=1;
    PIC::ParticleBuffer::GetX(local,particle);PIC::ParticleBuffer::GetV(local+3,particle);count=1;
  }
  // This geometry intentionally keeps the full path in its initial owned
  // leaf. Collective layout/cohort invariance is tested; migration is not
  // inferred from a no-crossing fixture. Core mover migration has other gates.
  PIC::Mover::CommitTempParticleMovingListsToFirstCellParticleTable();
  double global[6];int total=0,bad=0;
  MPI_Allreduce(local,global,6,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&count,&total,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&invalid,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  Check(total==1 && !bad,"global_species_weight_clock_and_physical_interval");
  Check(T::Acceleration==nullptr,"normal_total_acceleration_macro_dispatch");
  PIC::ParticleBuffer::DeleteAllParticles();
  State state;for(int k=0;k<6;++k)state[k]=global[k];return state;
}
void Cases() {
  Node* node=PIC::Mesh::mesh->findTreeNode(const_cast<double*>(X0));Require(node!=nullptr,"initial point outside mesh");
  // Only the actual owner evaluates production force/reference trajectories;
  // other ranks never pass an unallocated replica/halo as a physical node.
  // Broadcast the independently integrated state without reducing precision.
  for(int species : {_NA_SPEC_,_NA_PLUS_SPEC_}) {
    ForceWiring(species,node);State coarse{},reference{};
    if(node->Thread==PIC::ThisThread){coarse=Reference(species,node,4096);reference=Reference(species,node,8192);}
    MPI_Bcast(coarse.data(),6,MPI_LONG_DOUBLE,node->Thread,MPI_GLOBAL_COMMUNICATOR);
    MPI_Bcast(reference.data(),6,MPI_LONG_DOUBLE,node->Thread,MPI_GLOBAL_COMMUNICATOR);
    double position[3],velocity[3];const double steps[3]={1,.5,.25};
    for(int stage=0;stage<3;++stage) {
      const auto actual=Advance(species,steps[stage],node);
      position[stage]=static_cast<double>(Norm(Difference(actual,reference,0)));
      velocity[stage]=static_cast<double>(Norm(Difference(actual,reference,3)));
      if(PIC::ThisThread==0){Rows<<species<<','<<Order<<','<<steps[stage]<<','<<position[stage]<<','<<velocity[stage];
        for(int k=0;k<6;++k)Rows<<','<<static_cast<double>(actual[k]);
        for(int k=0;k<6;++k)Rows<<','<<static_cast<double>(reference[k]);Rows<<'\n';}
    }
    const double referenceX=static_cast<double>(Norm(Difference(coarse,reference,0))),referenceV=static_cast<double>(Norm(Difference(coarse,reference,3)));
    if(PIC::ThisThread==0) {
      ReferenceProof<<species<<",4096,8192,"<<referenceX<<','<<referenceV<<','<<position[2]<<','<<velocity[2]<<','
        <<referenceX/position[2]<<','<<referenceV/velocity[2]<<'\n';
      ReferenceProof.flush(); // preserve the failed oracle metric too
    }
    Require(referenceX<=.001*position[2] && referenceV<=.001*velocity[2],"RK4 reference insufficiently resolved; preserve outputs, do not relax gate");
    Check(position[0]>position[1] && position[1]>position[2] && velocity[0]>velocity[1] && velocity[1]>velocity[2],
        "position_and_velocity_errors_decrease");
    bool orders=true;
    for(int i=0;i<2;++i){const double px=std::log2(position[i]/position[i+1]),pv=std::log2(velocity[i]/velocity[i+1]);
      orders=orders && px>=Order-.2 && px<=Order+.2 && pv>=Order-.2 && pv<=Order+.2;}
    Check(orders,"declared_first_or_second_order_convergence");
  }
  // Independent source-time invariant: a stationary source born within the
  // species step must end at that declared step's endpoint. Start with a
  // genuinely uninitialized/recycled tracer record and invoke EACH actual
  // configured boundary macro, rather than setting its birth epoch ourselves.
  // These binary fractions make the epoch identities exactly representable.
  Moon::SurfaceEnvironment::BeginProductionParticleStep(Epoch);
  for(int withFace=0;withFace<2;++withFace) {
    bool valid=true;
    if(node->Thread==PIC::ThisThread) {
      const int species=_NA_SPEC_;
      const double duration=node->block->GetLocalTimeStep(species),remaining=duration/4;
      const long int particle=PIC::ParticleBuffer::GetNewParticle();
      auto* data=PIC::ParticleBuffer::GetParticleDataPointer(particle);
      PIC::ParticleBuffer::SetI(species,data);PIC::ParticleBuffer::SetX(const_cast<double*>(X0),data);
      PIC::ParticleBuffer::SetV(const_cast<double*>(V0),data);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1/node->block->GetLocalParticleWeight(species),data);
      PIC::ParticleBuffer::SetNext(-1,data);PIC::ParticleBuffer::SetPrev(-1,data);
      valid=!std::isfinite(T::GetParticleEpochS(particle));
      if(withFace) {
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_FACE_(particle,remaining,node,true,nullptr);
      }
      else {
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(particle,remaining,node,true);
      }
      const auto time=T::GetTimeContext(particle);
      valid=valid && T::GetParticleEpochS(particle)==Epoch+duration &&
          time.InitialEpochS==Epoch+duration-remaining && time.RequestedTimeS==remaining &&
          time.RemainingTimeS==0 && time.IntegrationOrder==Order;
    }
    PIC::Mover::CommitTempParticleMovingListsToFirstCellParticleTable();
    Check(valid,withFace?"configured_boundary_face_birth_clock_and_step_endpoint":
        "configured_boundary_birth_clock_and_step_endpoint");
    PIC::ParticleBuffer::DeleteAllParticles();
  }
}
}
int main(int argc,char** argv) {
  PIC::InitMPI();SetExitErrorCode(__LINE__,1);if(argc!=3){MPI_Finalize();return 2;}
  Order=std::stoi(argv[1]);if(Order!=1&&Order!=2){MPI_Finalize();return 2;}
  if(PIC::ThisThread==0) {
    CheckFile.open(std::string(argv[2])+".checks.csv");Rows.open(std::string(argv[2])+".states.csv");Terms.open(std::string(argv[2])+".terms.csv");
    ReferenceProof.open(std::string(argv[2])+".reference.csv");
    CheckFile<<"check,status\n";Rows<<std::setprecision(17)<<"species,order,h_s,position_error_m,velocity_error_m_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,ref_x_m,ref_y_m,ref_z_m,ref_vx_m_s,ref_vy_m_s,ref_vz_m_s\n";
    Terms<<std::setprecision(17)<<"species,term,ax_m_s2,ay_m_s2,az_m_s2\n";
    ReferenceProof<<std::setprecision(17)<<"species,coarse_steps,fine_steps,position_difference_m,velocity_difference_m_s,finest_position_error_m,finest_velocity_error_m_s,position_error_fraction,velocity_error_fraction\n";
  }
  int error=0;try {Setup();Cases();}catch(const std::exception& exception){error=1;std::cerr<<"ERROR Moon-force tracer: "<<exception.what()<<'\n';}
  int bad=0,failed=0;MPI_Allreduce(&error,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&Failures,&failed,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread==0){CheckFile.flush();Rows.flush();Terms.flush();std::cout<<"M3 actual Moon force tracer order="<<Order<<" checks="<<Checks<<" status="<<(bad?"ERROR":failed?"FAIL":"PASS")<<'\n';}
  MPI_Finalize();return bad?2:failed?1:0;
}
