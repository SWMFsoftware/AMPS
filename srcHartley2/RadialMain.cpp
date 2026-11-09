#include "pic.h"
#include "CometChemistry.h"
#include "CometBrightness.h"
#include "CometBrightnessMath.h"
#include "CometModel.h"
#include "CometPhysics.h"
#include "CometRadiativeCooling.h"
#include "CometRuntimeInput.h"
#include "LocalTimeStepValidation.h"
#include "RadialMesh.h"
#include "RadialSimulation.h"
#include "resonance_fluorescence.h"

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
      "--validate-surface-injection | --validate-brightness "
      "[--brightness-out-of-domain exclude|abort]] "
      "[--validation-particles N] [--validation-mover first|second|both] "
      "[--validation-skip-multiface] [--validation-unit-corrections] "
      "[--chemistry-validation-events N] [--random-seed N] "
      "[--chemistry-residual-sampling bernoulli|systematic] "
      "[--chemistry-estimator-mode configured|analog|enhanced] "
      "[--validate-runtime-input] [-i FILE [-case CASE]]\n",program);
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

// ---------------------------------------------------------------------------
// Particle-resolved brightness sampler fixture (--validate-brightness).
//
// The production sampler is driven by PIC's external-sampling lifecycle and
// only sees whatever particles happen to be in owner cells.  This fixture
// replaces the transport-dependent population with a deterministic global
// particle set whose membership is defined by a global index alone, so the
// union over MPI ranks is identical for every decomposition:
//
//   * out=0 is written after one sampling call with no particles at all and
//     must be identically zero (empty-window regression for the free-list
//     aliasing defect fixed on 2026-11-07);
//   * out=1 samples the seeded set.  The Python checker compares out=1 across
//     1/4/12-rank runs (decomposition equivalence) and against the
//     independent column densities written here with the PIC-free
//     Math::VisibleShellChordLength() helper;
//   * the first and the last physical radial cell are always populated;
//   * all seeded velocities stay inside every selected core model's
//     published radial-velocity interval, so the production ABORT policy is
//     not triggered by the regular set;
//   * with --brightness-out-of-domain, one extra H particle with a 500 km/s
//     local radial velocity is placed at r=1e7 m.  "exclude" selects
//     EXCLUDE_AND_REPORT and must finish with positive excluded counts;
//     "abort" selects ABORT and must terminate through MPI_Abort.
//
// The fixture never advances particles and never writes ordinary pic.* cell
// output.  It must not be used inside a science run.
namespace BrightnessFixture {
enum OutOfDomainMode {None=0,Exclude=1,Abort=2};

constexpr long int ParticlesPerSpecies=4000;
constexpr double OutOfDomainRadiusM=1.0e7;
constexpr double OutOfDomainLocalSpeedMPS=5.0e5;
constexpr int SamplingRepeats=3;

double FixtureImpactParameter(long int point) {
  // Same expression as the production grid in CometBrightness.cpp.  The
  // Python checker additionally compares this column with the Tecplot one.
  using namespace Comet::Brightness;
  if (ImpactParameterCount==1) return MinimumImpactParameter;
  const double fraction=static_cast<double>(point)/
      static_cast<double>(ImpactParameterCount-1);
  if (ImpactParameterGridMode==GridMode::Logarithmic)
    return std::exp(std::log(MinimumImpactParameter)+fraction*
        (std::log(MaximumImpactParameter)-std::log(MinimumImpactParameter)));
  return MinimumImpactParameter+fraction*
      (MaximumImpactParameter-MinimumImpactParameter);
}

struct cSeed {
  double Radius;
  int Species;
  double LocalVelocity[3];
};

struct cInserted {
  double InnerRadius,OuterRadius,Measure,Weight;
  int Species;
};

// Insert one global seed if this rank owns the block that contains it.
// Returns true when inserted locally.
bool InsertIfOwned(const cSeed &seed,std::vector<cInserted> &inserted) {
  double position[3]={seed.Radius,0.0,0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
      PIC::Mesh::mesh->findTreeNode(position);
  if (node==NULL || node->block==NULL || node->Thread!=PIC::ThisThread)
    return false;
  const double cellWidth=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  long int cell=static_cast<long int>(
      std::floor((seed.Radius-node->xmin[0])/cellWidth));
  cell=std::max(0L,std::min(cell,static_cast<long int>(_BLOCK_CELLS_X_-1)));
  // Dense physical-cell index for the particle list; ghost-aware index for
  // the center node.  See the matching comment in CometBrightness.cpp.
  PIC::Mesh::cDataCenterNode *centerNode=node->block->GetCenterNode(
      _getCenterNodeLocalNumber(cell,0,0));
  if (centerNode==NULL || !(centerNode->Measure>0.0)) return false;

  long int &first=node->block->FirstCellParticleTable[cell];
  const long int ptr=PIC::ParticleBuffer::GetNewParticle(first);
  PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::SetX(position,data);
  PIC::ParticleBuffer::SetV(const_cast<double*>(seed.LocalVelocity),data);
  PIC::ParticleBuffer::SetI(seed.Species,data);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
#endif
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
  PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
#endif
  cInserted record;
  record.InnerRadius=node->xmin[0]+cell*cellWidth;
  record.OuterRadius=record.InnerRadius+cellWidth;
  record.Measure=centerNode->Measure;
  record.Weight=node->block->GetLocalParticleWeight(seed.Species);
  record.Species=seed.Species;
  inserted.push_back(record);
  return true;
}

bool Run(int outOfDomainMode) {
  using namespace Comet::Brightness;
  bool localPass=true;
  auto fail=[&](const char *message) {
    if (PIC::ThisThread==0) std::fprintf(stderr,"brightness fixture: %s\n",message);
    localPass=false;
  };
  if (!Enabled || EmissionLineCount<=0) {
    fail("brightness is not enabled in the runtime input");
    return false;
  }

  // Emitting species, in first-appearance order of the configured lines.
  std::vector<int> species;
  for (int line=0;line<EmissionLineCount;line++)
    if (std::find(species.begin(),species.end(),EmissionSpecies[line])==
        species.end()) species.push_back(EmissionSpecies[line]);

  // Velocity amplitude: keep |e_r.(V_comet+R v)| inside every selected core
  // model's published interval.  |e_r.(R v)|<=|v|, so a bound on |v| is
  // sufficient for every quadrature direction.
  double headroom=15.0e3;
  for (int line=0;line<EmissionLineCount;line++) {
    const AMPS::ResonanceFluorescence::ModelMetadata *metadata=NULL;
    if (GFactorModel[line]==GFactor::HydrogenLymanAlphaIKL2020)
      metadata=AMPS::ResonanceFluorescence::GetMetadata(
          AMPS::ResonanceFluorescence::Model::HydrogenLymanAlphaIKL2020);
    else if (GFactorModel[line]==
        GFactor::HydroxylAX00SchleicherAhearn1988Unquenched)
      metadata=AMPS::ResonanceFluorescence::GetMetadata(AMPS::
          ResonanceFluorescence::Model::HydroxylAX00SchleicherAhearn1988Unquenched);
    if (metadata==NULL) continue;
    headroom=std::min(headroom,std::min(
        metadata->maximumRadialVelocityMPS-BulkHeliocentricRadialVelocity,
        BulkHeliocentricRadialVelocity-metadata->minimumRadialVelocityMPS));
  }
  if (!(headroom>1.0e3)) {
    fail("bulk comet velocity leaves no in-domain headroom for the fixture");
    return false;
  }
  const double amplitude=0.5*headroom;

  // ---- empty window: out=0 must be identically zero ----
  int localEmpty=PIC::ParticleBuffer::GetAllPartNum()==0 ? 1 : 0,globalEmpty=0;
  MPI_Allreduce(&localEmpty,&globalEmpty,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalEmpty==0) fail("particle buffer is not empty before seeding");
  SampleParticleResolvedBrightness();
  PIC::LastSampleLength=1;
  OutputSampledBrightness(0);

  // ---- deterministic global seed list ----
  const double rn=Comet::RadialModel::NucleusRadius;
  const double rout=Comet::RadialModel::GetDomainOuterRadius();
  std::vector<cSeed> seeds;
  for (std::size_t s=0;s<species.size();s++) {
    for (long int k=0;k<ParticlesPerSpecies+2;k++) {
      cSeed seed;
      seed.Species=species[s];
      if (k==ParticlesPerSpecies) seed.Radius=rn+1.0e-3;          // first cell
      else if (k==ParticlesPerSpecies+1) seed.Radius=rout-1.0e-3; // last cell
      else {
        const double u=(k+0.5)/static_cast<double>(ParticlesPerSpecies);
        seed.Radius=rn*std::exp(u*std::log(rout/rn));
      }
      const double phase=static_cast<double>(k)+0.37*static_cast<double>(s);
      seed.LocalVelocity[0]=0.75*amplitude*std::sin(1.7*phase+0.3);
      seed.LocalVelocity[1]=0.50*amplitude*std::cos(2.3*phase+0.1);
      seed.LocalVelocity[2]=0.35*amplitude*std::sin(0.7*phase+1.1);
      seeds.push_back(seed);
    }
  }
  if (outOfDomainMode!=None) {
    int hydrogen=-1;
    for (int line=0;line<EmissionLineCount;line++)
      if (GFactorModel[line]==GFactor::HydrogenLymanAlphaIKL2020)
        hydrogen=EmissionSpecies[line];
    if (hydrogen<0) {
      fail("out-of-domain fixture requires the H_LYA_IKL2020 line");
      return false;
    }
    cSeed seed;
    seed.Species=hydrogen;
    seed.Radius=OutOfDomainRadiusM;
    seed.LocalVelocity[0]=OutOfDomainLocalSpeedMPS;
    seed.LocalVelocity[1]=0.0;
    seed.LocalVelocity[2]=0.0;
    seeds.push_back(seed);
    OutOfRangeVelocityPolicy=outOfDomainMode==Exclude ?
        VelocityDomainPolicy::ExcludeAndReport : VelocityDomainPolicy::Abort;
  }

  std::vector<cInserted> inserted;
  long int localFirst=0,localLast=0;
  for (std::size_t i=0;i<seeds.size();i++) {
    if (!InsertIfOwned(seeds[i],inserted)) continue;
    if (seeds[i].Radius==rn+1.0e-3) localFirst++;
    if (seeds[i].Radius==rout-1.0e-3) localLast++;
  }
  long int localCounts[3]={static_cast<long int>(inserted.size()),
      localFirst,localLast},globalCounts[3]={0,0,0};
  MPI_Allreduce(localCounts,globalCounts,3,MPI_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  const long int expectedEdge=static_cast<long int>(species.size());
  if (globalCounts[0]!=static_cast<long int>(seeds.size()))
    fail("not every seeded particle found exactly one owner block");
  if (globalCounts[1]!=expectedEdge || globalCounts[2]!=expectedEdge)
    fail("first/last physical cell coverage is incomplete");

  // ---- sample the seeded set; timing is the max over ranks ----
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  const double t0=MPI_Wtime();
  for (int repeat=0;repeat<SamplingRepeats;repeat++)
    SampleParticleResolvedBrightness();
  double localSeconds=(MPI_Wtime()-t0)/SamplingRepeats,maxSeconds=0.0;
  MPI_Reduce(&localSeconds,&maxSeconds,1,MPI_DOUBLE,MPI_MAX,0,
      MPI_GLOBAL_COMMUNICATOR);
  PIC::LastSampleLength=SamplingRepeats;
  OutputSampledBrightness(1);

  // ---- independent column density: sum_p W/V * chord(b) ----
  const std::size_t rays=static_cast<std::size_t>(ImpactParameterCount);
  std::vector<double> localColumn(rays*species.size(),0.0),
      globalColumn(rays*species.size(),0.0);
  for (const cInserted &record:inserted) {
    const std::size_t s=static_cast<std::size_t>(std::find(species.begin(),
        species.end(),record.Species)-species.begin());
    for (std::size_t ray=0;ray<rays;ray++) {
      const double chord=Math::VisibleShellChordLength(record.InnerRadius,
          record.OuterRadius,FixtureImpactParameter(ray),rn);
      localColumn[ray*species.size()+s]+=record.Weight/record.Measure*chord;
    }
  }
  MPI_Reduce(localColumn.data(),globalColumn.data(),
      static_cast<int>(localColumn.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    const std::string path=std::string(PIC::OutputDataFileDirectory)+
        "/brightness_validation_expected_columns.csv";
    FILE *csv=std::fopen(path.c_str(),"w");
    if (csv==NULL) fail("cannot write expected-column CSV");
    else {
      std::fprintf(csv,"impact_index,impact_parameter_m,species_index,expected_column_m2\n");
      for (std::size_t ray=0;ray<rays;ray++)
        for (std::size_t s=0;s<species.size();s++)
          std::fprintf(csv,"%zu,%.17e,%d,%.17e\n",ray,
              FixtureImpactParameter(ray),species[s],
              globalColumn[ray*species.size()+s]);
      std::fclose(csv);
    }
  }

  int localResult=localPass ? 1 : 0,globalResult=0;
  MPI_Allreduce(&localResult,&globalResult,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    std::printf("%s I3 particle-resolved brightness fixture mpi-ranks=%d "
        "seeded=%ld first-cell=%ld last-cell=%ld velocity-amplitude-mps=%.6e "
        "out-of-domain=%s rays=%ld lines=%d sampler-seconds-per-call=%.6e\n",
        globalResult!=0 ? "PASS" : "FAIL",PIC::nTotalThreads,globalCounts[0],
        globalCounts[1],globalCounts[2],amplitude,
        outOfDomainMode==None ? "none" :
            (outOfDomainMode==Exclude ? "exclude" : "abort"),
        ImpactParameterCount,EmissionLineCount,maxSeconds);
    std::fflush(stdout);
  }
  return globalResult!=0;
}
} // namespace BrightnessFixture
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
  bool validateBrightness=false;
  int brightnessOutOfDomainMode=BrightnessFixture::None;
  bool validateRuntimeInput=false;
  bool randomSeedOverride=false;
  bool residualSamplingOverride=false;
  bool estimatorModeOverride=false;
  int estimatorMode=Comet::Chemistry::EstimatorMode::Configured;
  int residualSamplingMode=Comet::Chemistry::ResidualSamplingMode;
  long int chemistryValidationEvents=200000;
  std::uint64_t randomSeed=Comet::Model::RandomSeed;
  std::string runtimeInputPath,runtimeCaseName;
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
    else if (option=="--validate-brightness") {
      validateBrightness=true;
    }
    else if (option=="--brightness-out-of-domain") {
      if (++argument>=argc) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
      const std::string mode=argv[argument];
      if (mode=="exclude") brightnessOutOfDomainMode=BrightnessFixture::Exclude;
      else if (mode=="abort") brightnessOutOfDomainMode=BrightnessFixture::Abort;
      else {
        std::fprintf(stderr,"invalid --brightness-out-of-domain value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
    }
    else if (option=="--validate-runtime-input") {
      // This mode deliberately performs the same strict parse, selected-case
      // overlay, MPI hash agreement, and complete model validation as a
      // production run, then stops before AMR/particle allocation. It gives
      // every case in one runtime file a fast executable-level regression.
      validateRuntimeInput=true;
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
        residualSamplingMode=
            Comet::Chemistry::ResidualSampling::IndependentBernoulli;
      else if (mode=="systematic" || mode=="SYSTEMATIC")
        residualSamplingMode=
            Comet::Chemistry::ResidualSampling::RandomPhaseSystematic;
      else {
        std::fprintf(stderr,
            "invalid --chemistry-residual-sampling value: %s\n",
            argv[argument]);
        return EXIT_FAILURE;
      }
      residualSamplingOverride=true;
    }
    else if (option=="-i" || option=="--input") {
      if (++argument>=argc || !runtimeInputPath.empty()) {
        std::fprintf(stderr,"-i/--input requires exactly one file\n");
        return EXIT_FAILURE;
      }
      runtimeInputPath=argv[argument];
    }
    else if (option=="-case" || option=="--case") {
      if (++argument>=argc || !runtimeCaseName.empty()) {
        std::fprintf(stderr,"-case/--case requires exactly one case name\n");
        return EXIT_FAILURE;
      }
      runtimeCaseName=argv[argument];
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
      (validateSurfaceInjection ? 1 : 0)+
      (validateBrightness ? 1 : 0)+
      (validateRuntimeInput ? 1 : 0);
  if (validationModeCount>1) {
    std::fprintf(stderr,
        "select only one focused validation mode per invocation\n");
    return EXIT_FAILURE;
  }
  if (brightnessOutOfDomainMode!=BrightnessFixture::None &&
      !validateBrightness) {
    std::fprintf(stderr,
        "--brightness-out-of-domain requires --validate-brightness\n");
    return EXIT_FAILURE;
  }
  if (validationOptions.UseUnitCorrections && !validateLocalTimeStep) {
    std::fprintf(stderr,
        "--validation-unit-corrections requires --validate-local-timestep\n");
    return EXIT_FAILURE;
  }
  if (Comet::RuntimeInput::Required && runtimeInputPath.empty()) {
    std::fprintf(stderr,
        "this executable requires runtime model input: use -i FILE -case CASE\n");
    return EXIT_FAILURE;
  }
  if (runtimeInputPath.empty() && !runtimeCaseName.empty()) {
    std::fprintf(stderr,"-case requires -i FILE\n");
    return EXIT_FAILURE;
  }

  // Compile-time preprocessing fixes species, dimension, enabled core models,
  // block layout, and maximum allocation sizes.  A modern two-layer comet
  // executable then loads the physical run receipt below. Legacy inputs keep
  // RuntimeInput::Required=false and retain their generated values.
  // The science-gate residual override is deliberately applied before this
  // point and, more importantly, before RadialSimulation::Initialize() calls
  // Chemistry::RequestBlockData().  Therefore one frozen executable can run
  // both arms: Bernoulli reserves no residual bytes, while SYSTEMATIC reserves
  // and initializes the normal owner-block state.  Physics and chemistry
  // tables remain byte-identical between arms.
  PIC::Init_BeforeParser();
  if (!runtimeInputPath.empty()) {
    std::string inputError;
    if (!Comet::RuntimeInput::LoadAndApply(runtimeInputPath,runtimeCaseName,
        &inputError)) {
      if (PIC::ThisThread==0)
        std::fprintf(stderr,"invalid comet runtime input: %s\n",
            inputError.c_str());
      return EXIT_FAILURE;
    }
    if (PIC::ThisThread==0) {
      std::printf("$PREFIX: comet runtime input=%s case=%s\n",
          runtimeInputPath.c_str(),runtimeCaseName.c_str());
      std::fflush(stdout);
    }
  }
  if (residualSamplingOverride)
    Comet::Chemistry::ResidualSamplingMode=residualSamplingMode;
  if (!randomSeedOverride) randomSeed=Comet::Model::RandomSeed;
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

  if (validateRuntimeInput) {
    // Stop only after all model owners have validated the installed state.
    // In particular this constructs/checks the runtime chemistry network and
    // validates the generated collision database and runtime brightness
    // metadata, but it creates no mesh, block, particle, PT, or plot output.
    if (PIC::ThisThread==0) {
      std::printf("PASS comet runtime input validation file=%s case=%s\n",
          runtimeInputPath.empty() ? "<compiled>" : runtimeInputPath.c_str(),
          runtimeCaseName.empty() ? "<none>" : runtimeCaseName.c_str());
      std::fflush(stdout);
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
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
    if (validateBrightness) {
      const bool passed=BrightnessFixture::Run(brightnessOutOfDomainMode);
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
