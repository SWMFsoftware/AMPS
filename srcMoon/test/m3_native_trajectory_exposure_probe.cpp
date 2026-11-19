// Reuse the actual native fixture's PIC setup, source, transaction and packet
// kinematics. Only the analytical hazard changes. This is NOT a physical Na
// product/source profile and never supplies an observation or formal M3 PASS.
#include "MoonReactionController.h"
#include "MoonInput.h"
#include "MoonTrajectoryExposure.h"
#include "MoonSurfaceEnvironment.h"
namespace {void ConfigureTrajectory(Moon::Chemistry::Controller::Configuration&);void ConfigureRuntime(Moon::Runtime::Configuration&);}
#define MOON_CONTROLLER_CONFIGURATION_CALLBACK ConfigureTrajectory
#define MOON_CONTROLLER_RUNTIME_CONFIGURATION_CALLBACK ConfigureRuntime
#define main OriginalConstantControllerMain
#include "m3_native_reaction_controller_probe.cpp"
#undef main
#undef MOON_CONTROLLER_CONFIGURATION_CALLBACK
#undef MOON_CONTROLLER_RUNTIME_CONFIGURATION_CALLBACK
namespace {
int PreviewCalls=0;
void ConfigureRuntime(Moon::Runtime::Configuration&c){c.surfaceGeometry=Moon::Runtime::SurfaceGeometry::Lola;}
bool AffineExposure(const C::ExposureRequest& request,void*,
    C::ExposureHistory* output,std::string*error){
  ++PreviewCalls;
  // The explicitly zero-force native flight has x=x0+730*t, so this spatial
  // source has lambda(x)=1+(x-x0)/730 s^-1. Its independent cumulative hazard
  // is t+t²/2, NOT a constant-rate clock or a duplicated trajectory solver.
  if(request.arc.accelerationMS2[0]!=0 || request.arc.accelerationMS2[1]!=0 || request.arc.accelerationMS2[2]!=0){*error="analytical affine fixture requires the declared zero-force arc";return false;}
  const long double begin=1+(static_cast<long double>(request.arc.positionM[0])-Start[0])/730;
  const long double end=begin+request.arc.velocityMS[0]*request.maximumS/730;
  if(begin<0||end<0){*error="fixture left its declared nonnegative spatial domain";return false;}
  Moon::Species::PiecewiseHazard::Segment segment;segment.endS=request.maximumS;
  const std::size_t first=2*request.parentCatalogSpecies;
  if(request.parentCatalogSpecies==1)
    segment.channels={{first,.25L*7e-6,.25L*7e-6},{first+1,.75L*7e-6,.75L*7e-6}};
  else if(request.parentCatalogSpecies==2){
    const long double lightEnd=1001-static_cast<long double>(request.physicalEpochS);
    segment.channels={{first,lightEnd>0?.25L:0,lightEnd>0?.25L:0},{first+1,lightEnd>0?.75L:0,lightEnd>0?.75L:0}};
    if(lightEnd>0&&lightEnd<request.maximumS){auto dark=segment;segment.endS=lightEnd;dark.beginS=lightEnd;for(auto&rate:dark.channels)rate.beginRateS=rate.endRateS=0;output->segments={segment,dark};}
    else output->segments={segment};
    output->absoluteDepthErrorBound=0;return true;
  }
  else segment.channels={{first,.25L*begin,.25L*end},{first+1,.75L*begin,.75L*end}};
  output->segments={segment};output->absoluteDepthErrorBound=0;return true;
}
bool EndpointSeen=false,EndpointAcceptedContext=false;
bool AcceptedKinematics(const C::EventContext& event,R::KinematicPolicy*policy,std::string*error){
  if(Tag(event.particle)==2){EndpointSeen=true;
    EndpointAcceptedContext=event.selectedInstantaneousRateS==.25&&
      event.rateConventionIdentity=="declared-piecewise-affine-first-hit-LEFT-limit"&&
      event.acceptedExposureIdentity=="explicit-linear-position-hazard-v1";
  }
  return Kinematics(event,policy,error);
}
void ConfigureTrajectory(C::Configuration& c){
  c.identity="ANALYTICAL-SPATIAL-HAZARD-NOT-PHYSICAL-NA-v1";
  c.exposureProviderIdentity="explicit-linear-position-hazard-v1";
  c.exposureHistory=AffineExposure;c.exposureControls.maximumRootRefinements=128;
  for(auto& channel:c.channels){channel.scope=C::HazardScope::TrajectoryDependent;channel.kinematics=AcceptedKinematics;}
}
void TrajectoryCase(){
  T::Acceleration=DeclaredZeroAcceleration;Node* initial=PIC::Mesh::mesh->findTreeNode(Start.data());
  std::string error;bool pure=true;const int owner=initial->Thread;
  if(owner==PIC::ThisThread){
    const long int particle=Source(0,initial);C::PrepareParticle(particle,1,initial);
    C::History before,after;Require(C::ReadHistory(particle,&before,&error),error);
    const auto draws=Counts;
    for(double horizon:{1.0,.5,.25,1.0}){double event=-1;C::Query(particle,_NA_SPEC_,initial,1,horizon,event);}
    Require(C::ReadHistory(particle,&after,&error),error);
    pure=before.clock.remainingExtended()==after.clock.remainingExtended()&&before.interval.elapsedS()==after.interval.elapsedS()&&
      Counts.source==draws.source&&Counts.channel==draws.channel&&Counts.product==draws.product;
    int i,j,k;PIC::Mesh::mesh->FindCellIndex(Start.data(),i,j,k,initial,false);
    long int& head=initial->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    Moon::Chemistry::ProcessParticle(particle,head,initial);
  }
  Check(pure,"variable_preview_changes_no_physical_clock_or_random_draw");
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  const double expected=1000+std::sqrt(1+2*(-std::log(std::exp(-.4))))-1;
  bool candidate=true;for(std::size_t i=0;i<CandidateEpochs.size();++i)
    candidate=candidate&&std::abs(CandidateEpochs[i]-expected)<1e-10&&CandidateChannels[i]==0;
  Check(candidate&&Sum(CandidateEpochs.size())==1,"actual_accepted_spatial_hazard_event_matches_independent_integral");
  double children=0;bool completed=true;
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){Node* n=PIC::DomainBlockDecomposition::BlockTable[b];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int p=n->block->FirstCellParticleTable[cell];p!=-1;p=PIC::ParticleBuffer::GetNext(p)){
        ++children;C::History history;Require(C::ReadHistory(p,&history,&error),error);
        completed=completed&&!history.due&&history.interval.remainingS()==0&&T::GetParticleEpochS(p)==1001;
      }}
  Check(completed&&Sum(children)==2&&Sum(Ledger.events.size())==1,"variable_event_complete_atomic_packet_and_unused_daughter_intervals");
  Check(Sum(Counts.source)==1&&Sum(Counts.channel)==1&&Sum(Counts.product)==2,"variable_event_exact_source_channel_product_draw_ownership");
  Check(Sum(PreviewCalls)>0&&T::CountPending()==0,"variable_provider_called_by_normal_mover_and_all_continuations_drained");
}
void WeakClockCase(){
  PIC::ParticleBuffer::DeleteAllParticles();Ledger.events.clear();Counts={};
  CandidateEpochs.clear();CandidateRemaining.clear();CandidateTags.clear();CandidateChannels.clear();CandidateOwners.clear();
  Node* initial=PIC::Mesh::mesh->findTreeNode(Start.data());
  if(initial->Thread==PIC::ThisThread){const long int ptr=Source(1,initial);C::PrepareParticle(ptr,1,initial);
    int i,j,k;double x[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::Mesh::mesh->FindCellIndex(x,i,j,k,initial,false);
    long int& head=initial->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];Moon::Chemistry::ProcessParticle(ptr,head,initial);}
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  bool state=true;double count=0;std::string error;
  const long double expected=-std::log(std::exp(-.6))-2*static_cast<long double>(7e-6);
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){Node*n=PIC::DomainBlockDecomposition::BlockTable[b];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=n->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)){
        ++count;C::History h;Require(C::ReadHistory(ptr,&h,&error),error);
        state=state&&h.catalogSpecies==1&&!h.due&&!h.clock.consumed()&&h.interval.remainingS()==0&&
          T::GetParticleEpochS(ptr)==1002&&std::abs(h.clock.remainingExtended()-expected)<4*std::numeric_limits<double>::epsilon();
      }}
  Check(state&&Sum(count)==1&&Sum(Ledger.events.size())==0,"weak_realistic_hazard_preserves_compensated_clock_over_many_substeps_and_outer_intervals");
  Check(Sum(Counts.source)==1&&Sum(Counts.resolver)==1&&Sum(Counts.channel)==0&&Sum(Counts.product)==0,
      "weak_hazard_continuation_and_mpi_do_not_redraw_or_drop_low_component");
}
void LeftLimitCase(){
  // The complete accepted exposure reaches depth1 EXACTLY as illumination
  // turns off. The declared stochastic convention uses the LEFT-limit CDF;
  // a physical provider must not re-gate the already selected photon at t+.
  PIC::ParticleBuffer::DeleteAllParticles();Ledger.events.clear();Counts={};
  CandidateEpochs.clear();CandidateRemaining.clear();CandidateTags.clear();CandidateChannels.clear();CandidateOwners.clear();
  EndpointSeen=false;EndpointAcceptedContext=false;Phase="light-to-dark-endpoint";
  Node* initial=PIC::Mesh::mesh->findTreeNode(Start.data());
  if(initial->Thread==PIC::ThisThread){const long int ptr=Source(2,initial);double velocity[3]={0,0,0};PIC::ParticleBuffer::SetV(velocity,ptr);
    C::PrepareParticle(ptr,1,initial);int i,j,k;double x[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::Mesh::mesh->FindCellIndex(x,i,j,k,initial,false);
    long int&head=initial->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];Moon::Chemistry::ProcessParticle(ptr,head,initial);}
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  bool result=!EndpointSeen||EndpointAcceptedContext;
  for(double epoch:CandidateEpochs)result=result&&epoch==1001;
  Check(result&&Sum(EndpointSeen?1:0)==1&&Sum(CandidateEpochs.size())==1&&Sum(Ledger.events.size())==1&&
      Sum(Counts.source)==1&&Sum(Counts.channel)==1&&Sum(Counts.product)==2&&T::CountPending()==0,
      "accepted_photon_keeps_left_limit_rate_identity_at_exact_light_to_dark_endpoint");
}
void ProductionPhotonCase(){
  // A declared screen on identity frozen axes isolates the ACTUAL production
  // terrain/BVH-to-exposure bridge. It is not a LOLA data qualification or an
  // actual SPICE orientation. Native shadow geometry remains independent of
  // the algebraic tangency/screen integral used as the oracle below.
  for(int r=0;r<6;++r)for(int c=0;c<6;++c){Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[r][c]=r==c?1:0;Moon::OrbitalMotion::SO_to_IAU_TransformationMartix[r][c]=r==c?1:0;}
  Moon::Surface::Triangulation mesh;
  const double plane=Start[0]+1000,centre=Start[1];
  mesh.verticesM={{{plane,centre-10,Start[2]-1000}},{{plane,centre+10,Start[2]-1000}},
      {{plane,centre+10,Start[2]+1000}},{{plane,centre-10,Start[2]+1000}}};
  mesh.faces={{{0,1,2}},{{0,2,3}}};std::string error;
  Require(Moon::SurfaceEnvironment::InitializeProductionTerrain(mesh,Moon::Runtime::GetConfiguration(),&error),error);
  Moon::Species::TrajectoryExposure::Controls controls;controls.rootTimeToleranceS=1e-12;
  controls.absoluteDepthTolerance=1e-10;controls.relativeDepthTolerance=1e-7;controls.maximumSegments=10000;controls.maximumGeometryIntervals=1000000;
  std::shared_ptr<const Moon::Chemistry::Exposure::FrozenProductionField> field;
  const bool ready=Moon::Chemistry::Exposure::FrozenProductionField::Create(controls,false,&field,&error);
  Check(ready,"actual_frozen_production_mesh_and_bvh_installation");Require(ready,error);
  C::ExposureRequest request;request.node=PIC::Mesh::mesh->rootTree;request.physicalEpochS=1000;request.maximumS=1;
  // Find one allocated local block on every rank. The source-space operator
  // uses an immutable declared analytical arc, not a transported population.
  request.node=PIC::DomainBlockDecomposition::BlockTable[0];
  request.arc.positionM={{Start[0],centre-40,Start[2]}};request.arc.velocityMS={{0,80,0}};
  Moon::Species::TrajectoryExposure::History history;
  Require(field->Build(request,{{0,2}},{},&history,&error),error);
  Moon::Species::PiecewiseHazard::Exposure exposure;
  Require(Moon::Species::PiecewiseHazard::Integrate(history.segments,0,1,&exposure,&error)==Moon::Species::PiecewiseHazard::Status::Complete,error);
  const long double dx=Moon::xSun_SO[0]-Start[0],z=Start[2],scale=(Moon::xSun_SO[0]-plane)/dx;
  const long double entry=((centre-10)/scale-(centre-40))/80,leave=((centre+10)/scale-(centre-40))/80;
  const long double b=std::sqrt(dx*dx+z*z),ref=static_cast<long double>(_AU_)*_AU_;
  auto primitive=[&](long double t){return 2*ref/(80*b)*std::atan((centre-40+80*t)/b);};
  const long double oracle=primitive(entry)-primitive(0)+primitive(1)-primitive(leave);
  const long double bound=history.interpolationDepthErrorBound+history.shadowBoundaryDepthErrorBound+1e-10L;
  Check(entry>0&&leave<1&&std::abs(exposure.integratedHazard-oracle)<=bound,
      "actual_production_bvh_shadow_exposure_matches_independent_screen_integral");
  Check(field->geometry().frame==Moon::Frames::SolarOrbital&&history.silhouetteRoots>=2,
      "frozen_external_geometry_separate_from_actual_transport_epoch");
}
}
int main(int argc,char**argv){
  PIC::InitMPI();int code=2;
  try{
    if(argc!=3)throw std::runtime_error("usage: native-trajectory-probe ORDER OUTPUT_PREFIX");
    Order=std::stoi(argv[1]);if(Order!=1&&Order!=2)throw std::runtime_error("order must be1or2");
    const std::string prefix=argv[2];
    if(PIC::ThisThread==0){CheckFile.open(prefix+".checks.csv");CheckFile<<"check,status\n";}
    EventFile.open(prefix+".rank"+std::to_string(PIC::ThisThread)+".events.csv");EventFile<<std::setprecision(17)<<"phase,tag,channel,event_epoch_s,remaining_s,owner\n";
    Setup();TrajectoryCase();WeakClockCase();LeftLimitCase();ProductionPhotonCase();code=Failures?1:0;
    if(PIC::ThisThread==0)std::cout<<"M3 trajectory exposure order="<<Order<<" checks="<<Checks<<" status="<<(Failures?"FAIL":"PASS")<<std::endl;
  }catch(const std::exception&e){std::cerr<<"ERROR native trajectory exposure: "<<e.what()<<std::endl;}
  MPI_Finalize();return code;
}
