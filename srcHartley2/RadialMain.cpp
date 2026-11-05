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
#include <cstdint>
#include <string>

namespace {
void PrintUsage(const char *program) {
  std::printf("usage: %s [--validate-local-timestep | --validate-chemistry | "
      "--validate-weight-adaptation | --validate-correction-windows | "
      "--validate-systematic-residual | --validate-daughter-packets | "
      "--validate-channel-importance | --validate-total-reaction-forcing | "
      "--validate-unsupported-simplified-ntc] "
      "[--validation-particles N] [--validation-mover first|second|both] "
      "[--validation-skip-multiface] [--validation-unit-corrections] "
      "[--chemistry-validation-events N] [--random-seed N] "
      "[--chemistry-residual-sampling bernoulli|systematic]\n",program);
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
  bool validateCorrectionWindows=false;
  bool validateSystematicResidual=false;
  bool validateDaughterPackets=false;
  bool validateChannelImportance=false;
  bool validateTotalReactionForcing=false;
  bool validateUnsupportedSimplifiedNTC=false;
  bool randomSeedOverride=false;
  bool residualSamplingOverride=false;
  long int chemistryValidationEvents=200000;
  std::uint64_t randomSeed=Comet::Model::RandomSeed;
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
    else if (option=="--validate-correction-windows") {
      validateCorrectionWindows=true;
    }
    else if (option=="--validate-systematic-residual") {
      validateSystematicResidual=true;
    }
    else if (option=="--validate-daughter-packets") {
      validateDaughterPackets=true;
    }
    else if (option=="--validate-channel-importance") {
      validateChannelImportance=true;
    }
    else if (option=="--validate-total-reaction-forcing") {
      validateTotalReactionForcing=true;
    }
    else if (option=="--validate-unsupported-simplified-ntc") {
      validateUnsupportedSimplifiedNTC=true;
    }
    else if (option=="--random-seed") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      char *end=NULL;
      errno=0;
      const unsigned long long value=std::strtoull(argv[argument],&end,10);
      if (errno!=0 || end==argv[argument] || *end!='\0' || value<1 ||
          value>2147483000ULL) {
        std::fprintf(stderr,"invalid --random-seed value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
      randomSeed=static_cast<std::uint64_t>(value);
      randomSeedOverride=true;
    }
    else if (option=="--chemistry-residual-sampling") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      const std::string mode=argv[argument];
      if (mode=="bernoulli" || mode=="BERNOULLI")
        Comet::Chemistry::ResidualSamplingMode=
            Comet::Chemistry::ResidualSampling::IndependentBernoulli;
      else if (mode=="systematic" || mode=="SYSTEMATIC")
        Comet::Chemistry::ResidualSamplingMode=
            Comet::Chemistry::ResidualSampling::RandomPhaseSystematic;
      else {
        std::fprintf(stderr,
            "invalid --chemistry-residual-sampling value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
      residualSamplingOverride=true;
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
    else if (option=="--validation-unit-corrections") {
      validationOptions.UseUnitCorrections=true;
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
      (validateChemistry ? 1 : 0)+(validateWeightAdaptation ? 1 : 0)+
      (validateCorrectionWindows ? 1 : 0)+
      (validateSystematicResidual ? 1 : 0)+
      (validateDaughterPackets ? 1 : 0)+
      (validateChannelImportance ? 1 : 0)+
      (validateTotalReactionForcing ? 1 : 0)+
      (validateUnsupportedSimplifiedNTC ? 1 : 0);
  if (validationModeCount>1) {
    std::fprintf(stderr,
        "select only one focused validation mode per invocation\n");
    return EXIT_FAILURE;
  }
  if (validationOptions.UseUnitCorrections && !validateLocalTimeStep) {
    std::fprintf(stderr,
        "--validation-unit-corrections requires --validate-local-timestep\n");
    return EXIT_FAILURE;
  }

  // AMPS input preprocessing has already compiled case values into the source
  // tree.  This call initializes MPI, molecular data, and core services that
  // the validators use; no runtime text input is parsed here.
  // The science-gate residual override is deliberately applied before this
  // point and, more importantly, before RadialSimulation::Initialize() calls
  // Chemistry::RequestBlockData().  Therefore one frozen executable can run
  // both arms: Bernoulli reserves no residual bytes, while SYSTEMATIC reserves
  // and initializes the normal owner-block state.  Physics and chemistry
  // tables remain byte-identical between arms.
  PIC::Init_BeforeParser();
  if (residualSamplingOverride &&
      Comet::Chemistry::ResidualSamplingMode==
          Comet::Chemistry::ResidualSampling::RandomPhaseSystematic &&
      !Comet::Model::RandomSeedConfigured && !randomSeedOverride) {
    std::fprintf(stderr,
        "systematic command-line residual sampling requires an explicit seed\n");
    return EXIT_FAILURE;
  }
  if (Comet::Model::RandomSeedConfigured || randomSeedOverride) {
    Comet::Model::RandomSeed=randomSeed;
    // A command-line override is just as explicit as an input record and must
    // be recorded as valid provenance in Stage-3 center-cell tables.
    Comet::Model::RandomSeedConfigured=true;
    Comet::Chemistry::ResidualPhaseSeed=randomSeed;
    rnd_seed(static_cast<int>(randomSeed)+PIC::ThisThread);
    if (PIC::ThisThread==0) {
      std::printf("$PREFIX: comet random seed=%llu rank-stream=seed+rank residual-phase-stream=keyed\n",
          static_cast<unsigned long long>(randomSeed));
      std::fflush(stdout);
    }
  }

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
    if (validateCorrectionWindows) {
      const bool passed=
          Comet::RadialSimulation::ValidateIndividualCorrectionWindows();
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateSystematicResidual) {
      const bool passed=Comet::Chemistry::ValidateSystematicResidualState();
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateDaughterPackets) {
      const bool passed=Comet::Chemistry::ValidateDaughterPacketSplitting(
          chemistryValidationEvents);
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateChannelImportance) {
      const bool passed=Comet::Chemistry::
          ValidateChannelImportanceSampling(chemistryValidationEvents);
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateTotalReactionForcing) {
      const bool passed=Comet::Chemistry::ValidateTotalReactionForcing(
          chemistryValidationEvents);
      MPI_Finalize();
      return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (validateUnsupportedSimplifiedNTC) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      // The function must terminate through its correction-on setup guard
      // before dereferencing the deliberately null node.  Returning would
      // mean that an unsupported collision executor could silently run.
      PIC::MolecularCollisions::ParticleCollisionModel::
          ModelCellCollisions_ntc_simplified(0,0,0,NULL);
      std::fprintf(stderr,
          "unsupported simplified NTC unexpectedly accepted corrections\n");
#else
      std::fprintf(stderr,
          "unsupported simplified NTC guard test requires corrections on\n");
#endif
      MPI_Finalize();
      return EXIT_FAILURE;
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
