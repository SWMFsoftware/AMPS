#include "pic.h"
#include "CometChemistry.h"
#include "CometBrightness.h"
#include "CometModel.h"
#include "CometPhysics.h"
#include "CometRadiativeCooling.h"
#include "LocalTimeStepValidation.h"
#include "RadialMesh.h"
#include "RadialSimulation.h"

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace {
void PrintUsage(const char *program) {
  std::printf("usage: %s [--validate-local-timestep | --validate-chemistry | "
      "--validate-weight-adaptation | --validate-correction-windows | "
      "--validate-systematic-residual | --validate-daughter-packets | "
      "--validate-channel-importance | --validate-total-reaction-forcing | "
      "--validate-transported-daughter-collision | "
      "--validate-unsupported-simplified-ntc | "
      "--validate-surface-injection] "
      "[--validation-particles N] [--validation-mover first|second|both] "
      "[--validation-skip-multiface] [--validation-unit-corrections] "
      "[--chemistry-validation-events N] [--random-seed N] "
      "[--chemistry-residual-sampling bernoulli|systematic] "
      "[--chemistry-estimator-mode configured|analog|enhanced]\n",program);
}

// Return the compiled species number without depending on application-specific
// _H2O_SPEC_, _CO_SPEC_, or _CO2_SPEC_ macros.  The I3 input deliberately
// defines all three parents, but retaining a runtime lookup gives a clear test
// failure if somebody tries to run the focused validation with an incompatible
// preprocessed executable.
int FindSpecies(const char *requestedSymbol) {
  char symbol[_MAX_STRING_LENGTH_PIC_];
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    PIC::MolecularData::GetChemSymbol(symbol,species);
    if (std::string(symbol)==requestedSymbol) return species;
  }
  return -1;
}

double EmpiricalKolmogorovDistance(std::vector<double> &samples,
    double (*cdf)(double)) {
  std::sort(samples.begin(),samples.end());
  const double inverseCount=1.0/static_cast<double>(samples.size());
  double distance=0.0;
  for (std::size_t sample=0;sample<samples.size();sample++) {
    const double expected=cdf(samples[sample]);
    const double lower=static_cast<double>(sample)*inverseCount;
    const double upper=static_cast<double>(sample+1)*inverseCount;
    distance=std::max(distance,std::fabs(expected-lower));
    distance=std::max(distance,std::fabs(upper-expected));
  }
  return distance;
}

double FluxEnergyCdf(double q) {
  // With q=m*v^2/(2*k*T), a surface-crossing Maxwellian has the Gamma(2,1)
  // density q*exp(-q), rather than the volume-Maxwellian energy law.
  return q<=0.0 ? 0.0 : 1.0-(1.0+q)*std::exp(-q);
}

double CosineLawCdf(double mu) {
  // Maxwellian flux through a plane has p(mu)=2*mu on [0,1].
  if (mu<=0.0) return 0.0;
  if (mu>=1.0) return 1.0;
  return mu*mu;
}

bool ValidateSurfaceInjection(long int particlesPerSpecies) {
  const double temperature=175.0;
  const double bulkVelocity[3]={0.0,0.0,0.0};
  double externalNormal[3]={-1.0,0.0,0.0};
  const char *parentSymbols[3]={"H2O","CO","CO2"};
  bool localPass=true;

  if (particlesPerSpecies<1000) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "surface-injection validation requires at least 1000 particles per species\n");
    return false;
  }

  // A fixed per-rank stream makes this a reproducible regression while still
  // detecting rank-local initialization errors in an MPI invocation.
  rnd_seed(872341+PIC::ThisThread);
  for (const char *parentSymbol:parentSymbols) {
    const int species=FindSpecies(parentSymbol);
    if (species<0) {
      if (PIC::ThisThread==0)
        std::fprintf(stderr,
            "surface-injection validation requires compiled species %s\n",
            parentSymbol);
      localPass=false;
      continue;
    }

    const double mass=PIC::MolecularData::GetMass(species);
    const double beta=std::sqrt(mass/(2.0*Kbol*temperature));
    const double tangentialVariance=1.0/(2.0*beta*beta);
    const double expectedNormalMean=std::sqrt(Pi)/(2.0*beta);
    const double expectedNormalSecondMoment=1.0/(beta*beta);
    std::vector<double> fluxEnergies;
    std::vector<double> directionCosines;
    fluxEnergies.reserve(static_cast<std::size_t>(particlesPerSpecies));
    directionCosines.reserve(static_cast<std::size_t>(particlesPerSpecies));
    double sumNormal=0.0,sumNormalSquared=0.0;
    double sumTangential1=0.0,sumTangential2=0.0;
    double sumTangential1Squared=0.0,sumTangential2Squared=0.0;
    long int outwardFailures=0;

    for (long int particle=0;particle<particlesPerSpecies;particle++) {
      double velocity[3]={0.0,0.0,0.0};
      const double correction=PIC::Distribution::InjectMaxwellianDistribution(
          velocity,bulkVelocity,temperature,externalNormal,species);
      const double speedSquared=velocity[0]*velocity[0]+
          velocity[1]*velocity[1]+velocity[2]*velocity[2];
      const double speed=std::sqrt(speedSquared);
      const double normalVelocity=-
          (velocity[0]*externalNormal[0]+velocity[1]*externalNormal[1]+
           velocity[2]*externalNormal[2]);
      if (!(normalVelocity>0.0) || !(correction==1.0) || !(speed>0.0))
        outwardFailures++;
      fluxEnergies.push_back(beta*beta*speedSquared);
      directionCosines.push_back(normalVelocity/speed);
      sumNormal+=normalVelocity;
      sumNormalSquared+=normalVelocity*normalVelocity;
      sumTangential1+=velocity[1];
      sumTangential2+=velocity[2];
      sumTangential1Squared+=velocity[1]*velocity[1];
      sumTangential2Squared+=velocity[2]*velocity[2];
    }

    const double inverseCount=1.0/static_cast<double>(particlesPerSpecies);
    const double normalMean=sumNormal*inverseCount;
    const double normalSecondMoment=sumNormalSquared*inverseCount;
    const double tangentialMean1=sumTangential1*inverseCount;
    const double tangentialMean2=sumTangential2*inverseCount;
    const double tangentialVariance1=sumTangential1Squared*inverseCount-
        tangentialMean1*tangentialMean1;
    const double tangentialVariance2=sumTangential2Squared*inverseCount-
        tangentialMean2*tangentialMean2;
    const double energyKs=EmpiricalKolmogorovDistance(fluxEnergies,
        FluxEnergyCdf);
    const double cosineKs=EmpiricalKolmogorovDistance(directionCosines,
        CosineLawCdf);

    // Eight asymptotic standard errors are intentionally conservative.  The
    // small fixed floors keep the test meaningful at its documented minimum
    // sample count without making normal random variation a flaky failure.
    const double inverseRootCount=1.0/std::sqrt(
        static_cast<double>(particlesPerSpecies));
    const double ksTolerance=std::max(8.0*inverseRootCount,0.012);
    const double momentTolerance=std::max(8.0*inverseRootCount,0.02);
    const double tangentialMeanTolerance=
        8.0*std::sqrt(tangentialVariance)*inverseRootCount;
    const double tangentialVarianceTolerance=
        std::max(8.0*std::sqrt(2.0)*inverseRootCount,0.025);
    const bool speciesPass=outwardFailures==0 &&
        energyKs<=ksTolerance && cosineKs<=ksTolerance &&
        std::fabs(normalMean/expectedNormalMean-1.0)<=momentTolerance &&
        std::fabs(normalSecondMoment/expectedNormalSecondMoment-1.0)<=
            momentTolerance &&
        std::fabs(tangentialMean1)<=tangentialMeanTolerance &&
        std::fabs(tangentialMean2)<=tangentialMeanTolerance &&
        std::fabs(tangentialVariance1/tangentialVariance-1.0)<=
            tangentialVarianceTolerance &&
        std::fabs(tangentialVariance2/tangentialVariance-1.0)<=
            tangentialVarianceTolerance;
    localPass=localPass && speciesPass;
    if (PIC::ThisThread==0) {
      std::printf("%s surface-injection species=%s T=175K samples=%ld "
          "energy-KS=%.6e cosine-KS=%.6e normal-mean-ratio=%.6e "
          "normal-second-ratio=%.6e tangent-means=(%.6e,%.6e) "
          "tangent-variance-ratios=(%.6e,%.6e) outward-failures=%ld\n",
          speciesPass ? "PASS" : "FAIL",parentSymbol,particlesPerSpecies,
          energyKs,cosineKs,normalMean/expectedNormalMean,
          normalSecondMoment/expectedNormalSecondMoment,tangentialMean1,
          tangentialMean2,tangentialVariance1/tangentialVariance,
          tangentialVariance2/tangentialVariance,outwardFailures);
    }
  }

  int localResult=localPass ? 1 : 0,globalResult=0;
  MPI_Allreduce(&localResult,&globalResult,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    std::printf("%s I3 surface-injection validation mpi-ranks=%d "
        "species=H2O,CO,CO2 temperature-K=175 particles-per-species=%ld\n",
        globalResult!=0 ? "PASS" : "FAIL",PIC::nTotalThreads,
        particlesPerSpecies);
    std::fflush(stdout);
  }
  return globalResult!=0;
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
  bool validateTransportedDaughterCollision=false;
  bool validateUnsupportedSimplifiedNTC=false;
  bool validateSurfaceInjection=false;
  bool randomSeedOverride=false;
  bool residualSamplingOverride=false;
  bool estimatorModeOverride=false;
  int estimatorMode=Comet::Chemistry::EstimatorMode::Configured;
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
    else if (option=="--validate-transported-daughter-collision") {
      validateTransportedDaughterCollision=true;
    }
    else if (option=="--validate-unsupported-simplified-ntc") {
      validateUnsupportedSimplifiedNTC=true;
    }
    else if (option=="--validate-surface-injection") {
      validateSurfaceInjection=true;
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
    else if (option=="--chemistry-estimator-mode") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      const std::string mode=argv[argument];
      if (mode=="configured" || mode=="CONFIGURED")
        estimatorMode=Comet::Chemistry::EstimatorMode::Configured;
      else if (mode=="analog" || mode=="ANALOG")
        estimatorMode=Comet::Chemistry::EstimatorMode::Analog;
      else if (mode=="enhanced" || mode=="ENHANCED")
        estimatorMode=Comet::Chemistry::EstimatorMode::Enhanced;
      else {
        std::fprintf(stderr,
            "invalid --chemistry-estimator-mode value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
      estimatorModeOverride=true;
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
      (validateTransportedDaughterCollision ? 1 : 0)+
      (validateUnsupportedSimplifiedNTC ? 1 : 0)+
      (validateSurfaceInjection ? 1 : 0);
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
  if (estimatorModeOverride) {
    // The generated input has now initialized the compiled Stage-6/7/8
    // policy flags.  Select only the runtime estimator arm; never rewrite
    // those flags, because both arms must reserve the same particle bytes,
    // install the same mover callbacks, and pass the same generated-policy
    // validation before the first block or particle is allocated.
    std::string estimatorError;
    if (!Comet::Chemistry::ConfigureEstimatorMode(estimatorMode,
        &estimatorError)) {
      std::fprintf(stderr,"invalid chemistry estimator mode: %s\n",
          estimatorError.c_str());
      return EXIT_FAILURE;
    }
    if (PIC::ThisThread==0) {
      std::printf("$PREFIX: comet chemistry estimator mode=%s compiled-packets=%d compiled-channel-importance=%d compiled-total-forcing=%d identical-layout=on\n",
          Comet::Chemistry::GetEstimatorModeName(),
          Comet::Chemistry::DaughterPacketSplittingEnabled ? 1 : 0,
          Comet::Chemistry::ChannelImportanceSamplingEnabled ? 1 : 0,
          Comet::Chemistry::TotalReactionForcingEnabled ? 1 : 0);
      std::fflush(stdout);
    }
  }
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

  // The sampler needs only initialized MPI/random/molecular-data services.
  // Run before mesh construction so this focused distribution regression is
  // fast and cannot be confounded by radial transport, chemistry, or weights.
  if (validateSurfaceInjection) {
    const bool passed=ValidateSurfaceInjection(
        validationOptions.ParticlesPerCase);
    MPI_Finalize();
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
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
  if (!Comet::Brightness::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid comet brightness model: %s\n",
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
    if (validateTransportedDaughterCollision) {
      // Stage 9 starts with the complete Stage-8 production-mover fixture,
      // then feeds its exact transported daughter records to the real general
      // NTC cell kernel.  It is deliberately a separate mode so the accepted
      // Stage-8 parser and archive remain unchanged.
      const bool passed=Comet::Chemistry::ValidateTotalReactionForcing(
          chemistryValidationEvents,true);
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
