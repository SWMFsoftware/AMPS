#include "pic.h"
#include "CometChemistry.h"
#include "CometModel.h"
#include "CometPhysics.h"
#include "CometRadiativeCooling.h"
#include "LocalTimeStepValidation.h"
#include "RadialMesh.h"
#include "RadialSimulation.h"

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <string>

namespace {
void PrintUsage(const char *program) {
  std::printf("usage: %s [--validate-local-timestep | --validate-chemistry | "
      "--validate-weight-adaptation] "
      "[--validation-particles N] [--validation-mover first|second|both] "
      "[--validation-skip-multiface] "
      "[--chemistry-validation-events N]\n",program);
}
}

int main(int argc,char **argv) {
#if DIM != 1
  std::fprintf(stderr,"Comet RadialMain requires DIM=1\n");
  return EXIT_FAILURE;
#else
  bool validateLocalTimeStep=false;
  bool validateChemistry=false;
  bool validateWeightAdaptation=false;
  long int chemistryValidationEvents=200000;
  Comet::LocalTimeStepValidation::cOptions validationOptions;
  for (int argument=1;argument<argc;argument++) {
    const std::string option=argv[argument];
    if (option=="--validate-local-timestep") {
      validateLocalTimeStep=true;
    }
    else if (option=="--validate-chemistry") {
      validateChemistry=true;
    }
    else if (option=="--validate-weight-adaptation") {
      validateWeightAdaptation=true;
    }
    else if (option=="--validation-particles") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      char *end=NULL;
      errno=0;
      const long value=std::strtol(argv[argument],&end,10);
      if (errno!=0 || end==argv[argument] || *end!='\0' || value<=0) {
        std::fprintf(stderr,"invalid --validation-particles value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
      validationOptions.ParticlesPerCase=value;
    }
    else if (option=="--validation-mover") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      const std::string mover=argv[argument];
      validationOptions.TestFirstOrder=mover=="first" || mover=="both";
      validationOptions.TestSecondOrder=mover=="second" || mover=="both";
      if (!validationOptions.TestFirstOrder &&
          !validationOptions.TestSecondOrder) {
        std::fprintf(stderr,"invalid validation mover: %s\n",argv[argument]);
        return EXIT_FAILURE;
      }
    }
    else if (option=="--validation-skip-multiface") {
      validationOptions.TestMultiFace=false;
    }
    else if (option=="--chemistry-validation-events") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      char *end=NULL;
      errno=0;
      const long value=std::strtol(argv[argument],&end,10);
      if (errno!=0 || end==argv[argument] || *end!='\0' || value<=0) {
        std::fprintf(stderr,
            "invalid --chemistry-validation-events value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
      chemistryValidationEvents=value;
    }
    else if (option=="--help" || option=="-h") {
      PrintUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else {
      std::fprintf(stderr,"unknown option: %s\n",argv[argument]);
      PrintUsage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  const int validationModeCount=(validateLocalTimeStep ? 1 : 0)+
      (validateChemistry ? 1 : 0)+(validateWeightAdaptation ? 1 : 0);
  if (validationModeCount>1) {
    std::fprintf(stderr,
        "select only one focused validation mode per invocation\n");
    return EXIT_FAILURE;
  }

  // AMPS input preprocessing has already compiled case values into the source
  // tree.  This call initializes MPI, molecular data, and core services that
  // the validators use; no runtime text input is parsed here.
  PIC::Init_BeforeParser();

  std::string validationError;
  if (!Comet::Model::Validate(&validationError,
      Comet::RadialSimulation::Enabled)) {
    std::fprintf(stderr,"Invalid comet model: %s\n",validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::Chemistry::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid comet chemistry: %s\n",validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::Physics::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid comet physics database: %s\n",
        validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::RadiativeCooling::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid comet radiative cooling: %s\n",
        validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::RadialSimulation::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid radial simulation: %s\n",
        validationError.c_str());
    return EXIT_FAILURE;
  }

  // The computational coordinate is cometocentric radius.  DIM1 spherical
  // symmetry in the AMR core supplies shell volumes and radial face areas.
  Comet::RadialModel::InitMeshProfile();

  double xmin[1]={Comet::RadialModel::NucleusRadius};
  double xmax[1]={Comet::RadialModel::GetDomainOuterRadius()};

  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,Comet::RadialModel::LocalResolution);

  // One rank constructs and serializes the AMR tree.  Other ranks read the
  // same tree before AMPS creates its load-distribution lists.
  const char *meshFile="comet.radial.mesh.msh";
  if (PIC::ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile(meshFile);
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread!=0) PIC::Mesh::mesh->readMeshFile(meshFile);

  PIC::Mesh::mesh->SetParallelLoadMeasure([](
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) -> double {return 1.0;});
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  Comet::RadialModel::PrintMeshDiagnostics();

  if (Comet::RadialSimulation::Enabled) {
    Comet::RadialSimulation::Initialize();
    if (validateLocalTimeStep) {
      const bool passed=Comet::LocalTimeStepValidation::Run(validationOptions);
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateChemistry) {
      const bool passed=
          Comet::Chemistry::ValidateEnsemble(chemistryValidationEvents);
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateWeightAdaptation) {
      const bool passed=
          Comet::RadialSimulation::ValidateParticleWeightAdaptation();
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    Comet::RadialSimulation::Run();
  }
  else if (validationModeCount!=0) {
    std::fprintf(stderr,
        "focused validation requires RadialSimulationMode=on\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
#endif
}
