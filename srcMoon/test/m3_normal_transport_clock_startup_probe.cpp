// Thin NORMAL application startup driver: same parser, amps_init and exactly
// two real amps_time_step calls. No duplicate initializer, force/source law,
// synthetic observation or physical reaction provider is installed here.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonTransportClock.h"
#include "MoonSurfaceEnvironment.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
void amps_init();
void amps_time_step();
namespace {
int Checks=0,Failures=0,ExpectedStep=0;long int InitialSamples=0;
long long Emissions=0;bool SourceOrder=true;std::ofstream CheckFile,StepFile;
PIC::ParticleEvents::ObserverFunction PreviousObserver=nullptr;
void Require(bool valid,const std::string& message){if(!valid)throw std::runtime_error(message);}
void Check(bool valid,const char* name) {
  int local=valid?1:0,all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  ++Checks;if(!all)++Failures;
  if(PIC::ThisThread==0){CheckFile<<name<<','<<(all?"PASS":"FAIL")<<'\n';CheckFile.flush();}
}
void Observe(PIC::ParticleEvents::Event event,int species,double weight,int detail) {
  // Preserve any production inventory/budget callback exactly. Our passive
  // counters must not disable the application's ordinary accepted ledger.
  if(PreviousObserver)PreviousObserver(event,species,weight,detail);
  if(event!=PIC::ParticleEvents::Emission || species!=_NA_SPEC_)return;
  // Passive observation only: no RNG, callbacks, resampling, clock mutation
  // or MPI inside the source callback. Real collective checks follow the step.
  ++Emissions;const auto* clock=Moon::GetCommonTransportClock();
  SourceOrder=SourceOrder && clock && clock->active() &&
    clock->step()->sequence==static_cast<std::uint64_t>(ExpectedStep) &&
    PIC::CollectingSampleCounter==InitialSamples+ExpectedStep;
}
bool Endpoints(double expected,long long* count) {
  bool valid=true;long long local=0,all=0;
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  for(unsigned block=0;block<PIC::DomainBlockDecomposition::nLocalBlocks;++block){auto node=PIC::DomainBlockDecomposition::BlockTable[block];
    if(!node || !node->block || node->Thread!=PIC::ThisThread)continue;
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        ++local;valid=valid && Moon::TransportTime::EpochProjectionAgrees(PIC::Mover::ForceCutSurface::GetParticleEpochS(ptr),expected) &&
          PIC::Mover::ForceCutSurface::GetTimeContext(ptr).RemainingTimeS==0;
      }}
  MPI_Allreduce(&local,&all,1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);*count=all;return valid;
}
}
int main(int argc,char** argv) {
  if(argc!=3)return 2;
  Moon::Runtime::Configuration input;std::string error;
  if(!Moon::Runtime::ParseApplicationInput(argv[1],&input,&error) || !Moon::Runtime::InstallConfiguration(input,&error)) {
    std::cerr<<"ERROR normal common startup input: "<<error<<'\n';return 2;
  }
  amps_init(); // Calls the NORMAL MPI/SPICE/mesh/source/weight/buffer startup.
  SetExitErrorCode(__LINE__,1);const std::string prefix(argv[2]);
  if(PIC::ThisThread==0){CheckFile.open(prefix+".checks.csv");CheckFile<<"check,status\n";
    StepFile.open(prefix+".steps.csv");StepFile<<std::setprecision(21)<<"step,relative_elapsed_s,expected_relative_elapsed_s,projected_endpoint_et_s,external_et_s,sampling_counter,emission_count,gas_count\n";}
  int errorFlag=0;
  try {
    Check(PIC::nTotalSpecies==1 && _NA_SPEC_==0 && Moon::CommonTransportTimeSelected() && !input.m2EnvironmentEnabled,
        "canonical_single_na_common_mode_after_normal_initializer");
#ifdef _NO_SPICE_CALLS_
    throw std::runtime_error("normal startup control requires actual configured SPICE kernels");
#else
    SpiceDouble expectedOrigin=0;str2et_c(input.simulationStartUtc.c_str(),&expectedOrigin);
    const auto* clock=Moon::GetCommonTransportClock();Require(clock!=nullptr,"normal startup did not initialize common clock");
    Check(clock->originEtS()==expectedOrigin && Moon::OrbitalMotion::et==expectedOrigin,
        "normal_spice_utc_origin_matches_explicit_transport_origin");
    Check(clock->completedSteps()==0 && clock->elapsedS()==0 && !clock->active(),
        "normal_startup_has_no_accepted_or_provisional_transport_step");
    const double dt=PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
    Check(std::isfinite(dt) && dt>0 && PIC::RequiredSampleLength>2,
        "normal_common_dt_is_positive_and_two_steps_do_not_cross_output_window");
    PreviousObserver=PIC::ParticleEvents::Observer;
    Require(PreviousObserver!=Observe,"startup observer is already installed");
    PIC::ParticleEvents::Observer=Observe;InitialSamples=PIC::CollectingSampleCounter;
    long long totalSources=0;
    for(ExpectedStep=1;ExpectedStep<=2;++ExpectedStep) {
      Emissions=0;SourceOrder=true;Moon::nIterationCounter=ExpectedStep-1;
      amps_time_step(); // Includes real Begin, Sampling, sources, mover, Commit.
      long long globalEmissions=0;MPI_Allreduce(&Emissions,&globalEmissions,1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      totalSources+=globalEmissions;
      Check(globalEmissions>0,ExpectedStep==1?"first_normal_source_cohort_recorded":"second_normal_source_cohort_recorded");
      const long double expectedElapsed=static_cast<long double>(ExpectedStep)*dt;
      double endpoint=0;Require(clock->Project(expectedElapsed,&endpoint,&error),error);
      Check(clock->completedSteps()==static_cast<std::uint64_t>(ExpectedStep) && clock->elapsedS()==expectedElapsed && !clock->active(),
          ExpectedStep==1?"first_actual_time_step_commits_full_common_duration":"second_actual_time_step_commits_full_common_duration");
      Check(SourceOrder && PIC::CollectingSampleCounter==InitialSamples+ExpectedStep,
          ExpectedStep==1?"first_normal_sampling_precedes_actual_source_callbacks":"second_normal_sampling_precedes_actual_source_callbacks");
      long long gas=0;const bool endpoints=Endpoints(endpoint,&gas);
      Check(endpoints && PIC::Mover::ForceCutSurface::CountPending()==0,
          ExpectedStep==1?"first_actual_population_finishes_common_endpoint":"second_actual_population_finishes_common_endpoint");
      Check(Moon::OrbitalMotion::et==expectedOrigin && Moon::SurfaceEnvironment::ProductionSourceEpochS()==expectedOrigin,
          ExpectedStep==1?"first_external_geometry_coefficient_epoch_remains_frozen":"second_external_geometry_coefficient_epoch_remains_frozen");
      if(PIC::ThisThread==0)StepFile<<ExpectedStep<<','<<clock->elapsedS()<<','<<expectedElapsed<<','<<endpoint<<','
        <<Moon::OrbitalMotion::et<<','<<PIC::CollectingSampleCounter<<','<<globalEmissions<<','<<gas<<'\n';
    }
    PIC::ParticleEvents::Observer=PreviousObserver;
    // This is execution coverage of the fixed-seed canonical source callback,
    // not an observational rate/agreement or statistical-precision criterion.
    Check(totalSources>0,"actual_normal_source_callbacks_were_exercised");
    SpiceInt count=0;ktotal_c("ALL",&count);bool kernels=count>0;
    if(PIC::ThisThread==0) {
      std::ofstream dependencies(prefix+".spice_files.txt");
      for(SpiceInt index=0;index<count;++index) {
        SpiceChar file[4096],kind[128],source[4096];SpiceInt handle;SpiceBoolean found;
        kdata_c(index,"ALL",sizeof(file),sizeof(kind),sizeof(source),file,kind,source,&handle,&found);
        if(!found){kernels=false;continue;}dependencies<<file<<'\n';
      }
      std::cout<<"SPICE_TOOLKIT "<<tkvrsn_c("TOOLKIT")<<'\n';
    }
    Check(kernels,"actual_furnished_spice_dependencies_recorded");
#endif
  }catch(const std::exception& exception){errorFlag=1;std::cerr<<"ERROR normal common startup: "<<exception.what()<<'\n';}
  if(PIC::ParticleEvents::Observer==Observe)PIC::ParticleEvents::Observer=PreviousObserver;
  int bad=0,failed=0;MPI_Allreduce(&errorFlag,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);MPI_Allreduce(&Failures,&failed,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread==0){CheckFile.flush();StepFile.flush();std::cout<<"M3 normal common startup checks="<<Checks<<" status="<<(bad?"ERROR":failed?"FAIL":"PASS")<<'\n';}
  MPI_Finalize();return bad?2:failed?1:0;
}
