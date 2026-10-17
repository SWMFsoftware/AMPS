//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf



//the particle class
#include "pic.h"
#include "constants.h"
#include "MoonInput.h"
#include "MoonIntegrationTests.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <new>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


void amps_init();
void amps_time_step();

namespace {
// std::new_handler is the last in-process opportunity to emit a controlled
// AMPS diagnostic before C++ allocation failure unwinds through legacy code.
// There is deliberately no generic campaign memory or wall-time limit. I11
// applies its reviewed 5000-particle-step bound inside the linked diagnostic;
// I17 periodic thermal spin-up still runs until its convergence contract is
// met. This handler is reached only when
// the allocator reports actual exhaustion.  Keep the text static because
// constructing a std::string at that point could itself require the
// unavailable memory.
void MoonAllocationFailure() {
  ::exit(__LINE__,__FILE__,
      "srcMoon terminated because a C++ memory allocation failed: the host "
      "or job memory available to this process was exhausted while building "
      "or advancing the requested AMR/thermal model. No M1 or I17 acceptance "
      "result may be reported from this incomplete run. Preserve stdout/stderr and the "
      "resource_metrics.csv artifact. To reduce required resources, increase "
      "the requested AMR cell sizes, narrow the near-surface/far-field or "
      "angular/polar refinement regions, or increase the surface-mesh target "
      "length. For I11, changing the approved mesh sequence creates a different "
      "campaign and must be documented as SKIPPED until the roadmap contract "
      "is revised; alternatively request a job with more memory and rerun the "
      "unchanged receipt. For I17, do not "
      "loosen the spin-up temperature or enthalpy tolerances after seeing "
      "this error; request more memory or define and qualify a lower-cost "
      "thermal grid in a separately reviewed input/convergence campaign.");
}
}

int main(int argc,char **argv) {
  std::set_new_handler(MoonAllocationFailure);
  // Keep the historical no-argument invocation for the regression test.  A
  // production input file is supplied explicitly so that srcMoon does not
  // accidentally interpret an AMPS configuration file intended for another
  // application.  Parse before PIC::InitMPI(): syntax/path failures then
  // produce one deterministic diagnostic without partially initialized ranks.
  std::string applicationInput,integrationTest,integrationOutput;
  for (int argument=1;argument<argc;argument++) {
    const std::string option(argv[argument]);
    if ((option=="-input") && argument+1<argc) {
      applicationInput=argv[++argument];
    }
    else if ((option=="--moon-integration-test") && argument+1<argc) {
      integrationTest=argv[++argument];
    }
    else if ((option=="--integration-output") && argument+1<argc) {
      integrationOutput=argv[++argument];
    }
    else {
      std::cerr << "Usage: " << argv[0]
                << " [-input application.in]"
                << " [--moon-integration-test U01|U01_INSIDE|U06|I02|I03|I05|I06|I07|I08|I09|I10|I11|I12|I13|I14|I16|I17|I18"
                << " --integration-output receipt.txt]\n";
      return EXIT_FAILURE;
    }
  }

  if (integrationTest.empty()!=integrationOutput.empty()) {
    std::cerr << "--moon-integration-test and --integration-output must be used together\n";
    return EXIT_FAILURE;
  }

  if (!applicationInput.empty()) {
    Moon::Runtime::Configuration configuration;
    std::string error;
    if (!Moon::Runtime::ParseApplicationInput(applicationInput,
                                               &configuration,&error)) {
      std::cerr << error << '\n';
      return EXIT_FAILURE;
    }
    // Installation transfers the validated receipt to amps_init().  It is a
    // one-shot operation so geometry cannot change after MPI/AMR startup.
    if (!Moon::Runtime::InstallConfiguration(configuration,&error)) {
      std::cerr << error << '\n';
      return EXIT_FAILURE;
    }
  }
  
  clock_t runtime =-clock();

PIC::Debugger::cGenericTimer t;
  
t.Start("main",__LINE__);
  amps_init();

  // M0 diagnostics execute only after the complete production initializer.
  // Consequently a PASS can establish linked wiring rather than source-level
  // presence.  The normal no-argument science/regression path is unchanged.
  if (!integrationTest.empty()) {
    const int status=Moon::IntegrationTests::Run(integrationTest,
                                                 integrationOutput);
    MPI_Finalize();
    return status;
  }

t.SwitchTimeSegment(__LINE__,"first switch"); 

//  PIC::Debugger::logger.InitLogger(PIC::ThisThread); 

  if (_PIC_LOGGER_MODE_==_PIC_MODE_ON_) {
    PIC::Debugger::LoggerData.erase();
    PIC::Debugger::logger.func_enter(__LINE__,"main()",&PIC::Debugger::LoggerData,0,5);
  }

  int nTotalIterations=100000001;
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=100;

  //time step
  static int LastDataOutputFileNumber=0;

  

  for (Moon::nIterationCounter=0;Moon::nIterationCounter<nTotalIterations;Moon::nIterationCounter++) {

    if (_PIC_LOGGER_MODE_==_PIC_MODE_ON_) {
      PIC::Debugger::LoggerData.erase();
      snprintf(PIC::Debugger::LoggerData.msg,sizeof(PIC::Debugger::LoggerData.msg),
          "line=%d,iter=%d",__LINE__,Moon::nIterationCounter);
      PIC::Debugger::logger.add_data_point(__LINE__,&PIC::Debugger::LoggerData);
    }

t.SwitchTimeSegment(__LINE__);
    amps_time_step();

t.SwitchTimeSegment(__LINE__,"another switch");

    if (PIC::Mesh::mesh->ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf(": (%i/%i %i:%i:%i), Iteration: %i  (currect sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,Moon::nIterationCounter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>40000) PIC::RequiredSampleLength=40000;


      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
    }

t.SwitchTimeSegment(__LINE__);

  }

t.SwitchTimeSegment(__LINE__);

  // OutputDataFileDirectory can be longer than the historical 400-byte local
  // buffer, so retain the complete path instead of truncating or overflowing it.
  const std::string fname=std::string(PIC::OutputDataFileDirectory)+"/test_Moon.dat";
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname.c_str());

t.SwitchTimeSegment(__LINE__);
t.Stop(__LINE__);


if (PIC::ThisThread==0) t.PrintSampledData(); 

t.PrintSampledDataMPI();

t.Start("part2",__LINE__);
t.Stop(__LINE__);
if (PIC::ThisThread==0) t.PrintSampledData();


t.clear();
t.Start("part3",__LINE__,__FILE__);
t.Stop(__LINE__);
if (PIC::ThisThread==0) t.PrintSampledData();


//  cout << "End of the run:" << PIC::nTotalSpecies << endl;
  MPI_Finalize();

  runtime+=clock();
  
  if(PIC::Mesh::mesh->ThisThread==0)
    cout << "Total AMPS runtime is "
	 << (double)runtime / CLOCKS_PER_SEC << endl;

  return EXIT_SUCCESS;
}
