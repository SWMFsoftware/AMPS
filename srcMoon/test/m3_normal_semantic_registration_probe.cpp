// Actual normal startup, furnished SPICE and semantic storage admission. This
// is an analytical infrastructure control; no physical channel is installed.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonSemanticRandom.h"
#include "MoonNativeReactionRegistration.h"
#include "MoonTrajectoryExposure.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
void amps_init();
namespace {
namespace S=Moon::Chemistry::SemanticRandom;
namespace N=Moon::Chemistry::NativeRegistration;
int Calls=0,Checks=0,Failures=0;bool Before=false,Geometry=false;
std::shared_ptr<const Moon::Chemistry::Exposure::FrozenProductionField> Field;
std::ofstream Output;
S::Configuration Config(){S::Configuration c;c.seed=23;
  c.cohort="DECLARED-SEMANTIC-PIC-COMPOSITION-NOT-LUNAR-PHYSICS";c.genealogyCapacity=4;return c;}
void Require(bool ok,const std::string& error){if(!ok)throw std::runtime_error(error);}
void Check(bool ok,const char* name){int local=ok,all=0;
  MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);++Checks;if(!all)++Failures;
  if(PIC::ThisThread==0){Output<<name<<','<<(all?"PASS":"FAIL")<<'\n';Output.flush();}}
bool Install(void*,std::string* error){
  ++Calls;Before=!PIC::ParticleBuffer::ParticleDataBuffer && PIC::nTotalSpecies>0;
  if(!Before || !S::RegisterStorage(Config(),error))return false;
#ifdef _NO_SPICE_CALLS_
  *error="normal semantic geometry control needs actual furnished SPICE";return false;
#else
  Moon::Species::TrajectoryExposure::Controls controls;controls.rootTimeToleranceS=1e-10;
  controls.absoluteDepthTolerance=1e-9L;controls.relativeDepthTolerance=1e-9L;
  controls.maximumSegments=10000;controls.maximumGeometryIntervals=10000;
  if(!Moon::Chemistry::Exposure::FrozenProductionField::CreateAtEpoch(controls,true,Moon::OrbitalMotion::et,&Field,error))return false;
  const auto& field=Field->geometry();double sun[6],earth[6],lightTime=0,expectedEpoch=0;
  str2et_c(Moon::Runtime::GetConfiguration().simulationStartUtc.c_str(),&expectedEpoch);
  spkezr_c("SUN",Moon::OrbitalMotion::et,Moon::Frames::SolarOrbital,"NONE","MOON",sun,&lightTime);
  spkezr_c("EARTH",Moon::OrbitalMotion::et,Moon::Frames::SolarOrbital,"NONE","MOON",earth,&lightTime);
  // CreateAtEpoch's documented reference is the SIX-component SPICE STATE.
  // A separate direct state request verifies epoch/frame/observer and km->m
  // wiring exactly. spkpos follows another floating-point path for dynamic
  // LSO frames: its near-zero Sun transverse component can differ by tens of
  // micrometres at1AU. Bitwise equality BETWEEN those APIs is not their
  // contract. We neither introduce a physical tolerance nor alter geometry.
  Geometry=field.frozenGeometryEpochS==Moon::OrbitalMotion::et && Moon::OrbitalMotion::et==expectedEpoch && !field.opaqueSpheres.empty() &&
    std::hypot(sun[0],sun[1],sun[2])>0 && std::hypot(earth[0],earth[1],earth[2])>0;
  for(int k=0;k<3;++k)Geometry=Geometry && field.sunPositionM[k]==1000*sun[k] &&
    field.opaqueSpheres[0].centreM[k]==1000*earth[k];
  if(!Geometry){*error="pre-buffer provider differs from independently furnished SPICE positions/epoch";return false;}
  return true;
#endif
}
}
int main(int argc,char** argv){
  if(argc!=4 || (std::string(argv[3])!="1" && std::string(argv[3])!="2"))return 2;
  const std::string expectedMover=std::string(argv[3])=="1"?"force-cut-surface-first-order":"force-cut-surface-second-order";
  bool mpi=false;int infrastructure=0;
  try{
    Moon::Runtime::Configuration input;std::string error;
    Require(Moon::Runtime::ParseApplicationInput(argv[1],&input,&error),error);
    Require(Moon::Runtime::InstallConfiguration(input,&error),error);
    Require(N::SetBeforeParticleBuffer(Install,nullptr,"ANALYTICAL-normal-semantic-registration",&error),error);
    amps_init();mpi=true;SetExitErrorCode(__LINE__,1);
    if(PIC::ThisThread==0){Output.open(std::string(argv[2])+".checks.csv");Output<<"check,status\n";}
    Check(Calls==1 && Before && S::Registered() && PIC::ParticleBuffer::ParticleDataBuffer &&
        Moon::Runtime::GetConfiguration().particleTracingMover==expectedMover,
      "normal_amps_init_executes_semantic_registration_before_real_buffer");
    Check(Geometry && Field && Field->geometry().frozenGeometryEpochS==Moon::OrbitalMotion::et,
      "normal_prebuffer_photon_field_uses_actual_furnished_spice_epoch_and_positions");
    Check(!N::RunBeforeParticleBuffer(&error) && Calls==1,"normal_prebuffer_callback_cannot_run_twice");
    Check(!N::SetBeforeParticleBuffer(Install,nullptr,"late",&error),"normal_prebuffer_callback_cannot_be_replaced_after_allocation");
    const long int ptr=PIC::ParticleBuffer::GetNewParticle();Require(S::InitializeSource(ptr,73,&error),error);
    S::State state;Require(S::Read(ptr,&state,&error),error);
    Check(state.originalParent==73 && state.genealogy.empty() && state.acceptedDraws==std::array<std::uint64_t,4>{{0,0,0,0}},
      "normal_initialized_pic_buffer_accepts_explicit_source_semantic_state");
    PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
#ifndef _NO_SPICE_CALLS_
    if(PIC::ThisThread==0){std::ofstream dependencies(std::string(argv[2])+".spice_files.txt");SpiceInt count=0;ktotal_c("ALL",&count);
      for(SpiceInt i=0;i<count;++i){SpiceChar path[4096],kind[128],source[4096];SpiceInt handle;SpiceBoolean found;
        kdata_c(i,"ALL",sizeof(path),sizeof(kind),sizeof(source),path,kind,source,&handle,&found);
        Require(found,"missing furnished SPICE dependency");dependencies<<path<<'\n';}}
#endif
  }catch(const std::exception& error){infrastructure=1;std::cerr<<"ERROR normal semantic registration: "<<error.what()<<'\n';}
  if(!mpi)return 2;int bad=0,failed=0;MPI_Allreduce(&infrastructure,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&Failures,&failed,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread==0){Output.flush();std::cout<<"M3 normal semantic registration checks="<<Checks<<" status="<<(bad?"ERROR":failed?"FAIL":"PASS")<<'\n';}
  MPI_Finalize();return bad?2:failed?1:0;
}
