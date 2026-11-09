#include "CometBrightness.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "CometBrightnessMath.h"
#include "CometModel.h"
#include "RadialMesh.h"
#include "pic.h"
#include "resonance_fluorescence.h"

bool Comet::Brightness::Enabled=false;
double Comet::Brightness::MinimumImpactParameter=-1.0;
double Comet::Brightness::MaximumImpactParameter=-1.0;
long int Comet::Brightness::ImpactParameterCount=0;
int Comet::Brightness::ImpactParameterGridMode=Comet::Brightness::GridMode::Linear;
int Comet::Brightness::EmissionLineCount=0;
int Comet::Brightness::EmissionSpecies[Comet::Brightness::MaxEmissionLines]={0};
double Comet::Brightness::PhotonGFactor[Comet::Brightness::MaxEmissionLines]={0.0};
int Comet::Brightness::GFactorModel[Comet::Brightness::MaxEmissionLines]={0};
double Comet::Brightness::EffectivePhotonWavelength[Comet::Brightness::MaxEmissionLines]={0.0};
double Comet::Brightness::CometHeliocentricVelocity[3]={0.0,0.0,0.0};
double Comet::Brightness::HeliocentricRadialDirection[3]={1.0,0.0,0.0};
double Comet::Brightness::LineOfSightDirection[3]={0.0,0.0,1.0};
double Comet::Brightness::ImpactParameterDirection[3]={1.0,0.0,0.0};
bool Comet::Brightness::ObservationGeometryConfigured=false;
double Comet::Brightness::BulkHeliocentricRadialVelocity=0.0;
double Comet::Brightness::SolarLymanAlphaFlux1AU=0.0;
// Keep this generated setting on one physical line.  hartley2.pl replaces the
// complete declaration, so splitting the initializer would leave an orphaned
// continuation in the configured source tree.
int Comet::Brightness::OutOfRangeVelocityPolicy=Comet::Brightness::VelocityDomainPolicy::Abort;
int Comet::Brightness::TauDiagnosticMode[Comet::Brightness::MaxEmissionLines]={0};
double Comet::Brightness::LineCenterCrossSection[Comet::Brightness::MaxEmissionLines]={0.0};
double Comet::Brightness::TauFlagThreshold[Comet::Brightness::MaxEmissionLines]={0.0};

namespace {
char TransitionName[Comet::Brightness::MaxEmissionLines]
    [Comet::Brightness::MetadataStringLength]={{'\0'}};
char GFactorProvenance[Comet::Brightness::MaxEmissionLines]
    [Comet::Brightness::MetadataStringLength]={{'\0'}};
char TauProvenance[Comet::Brightness::MaxEmissionLines]
    [Comet::Brightness::MetadataStringLength]={{'\0'}};

bool MetadataLoaded=false;
bool SamplingInitialized=false;
bool SamplingRegistered=false;

// Everything in this record is fixed for the life of a configured run.  It
// is calculated once during InitializeSampling(), rather than rebuilding the
// sky-plane offset and b^2 for every occupied shell in every sampling step.
// The shell intersections still depend on the AMR cell boundaries and are
// therefore evaluated where the owner cell is visited below.
struct cRayGeometry {
  double ImpactParameter;
  double ImpactParameterSquared;
  double SkyPlaneOffset[3];
};
std::vector<cRayGeometry> RayGeometry;
std::vector<double> SampledColumnDensity;
std::vector<double> SampledBrightnessRayleigh;
std::vector<double> ExcludedColumnDensity;
std::vector<unsigned long long> ExcludedEvaluationCount;

constexpr double AstronomicalUnitM=1.495978707e11;
constexpr double InverseSqrtThree=0.577350269189625764509148780502;

// Load-balance accounting for the sampler.  Radial comet builds that select
// ForceRepeatableSimulationPath use particle-number balancing: the mover adds
// one unit to a block's ParallelLoadMeasure per particle move.  The sampler
// cost is not proportional to particle number.  It is
//   (emitting particles in a cell) x (rays crossing that shell)
//   x (2 or 4 quadrature points) g-factor evaluations,
// and it is concentrated in outer shells, which almost every ray crosses.
// In the 2026-11-08 I3 pilots, the one rank that owned the outer H/OH front
// did almost all of the 200-ray work while 63 others waited in the next
// collective (128 ranks were no faster than 64). One evaluation (rotation
// plus core model) measured about 60 ns, versus about 1.5 us for one
// particle's complete iteration, hence 0.04 particle-move equivalents per
// evaluation.  The count is deterministic, so decomposition stays
// repeatable.  It changes MPI block ownership only, never a sampled value.
constexpr double LoadPerGFactorEvaluation=0.04;

std::size_t SampleIndex(long int impact,int line) {
  return static_cast<std::size_t>(impact)*
      Comet::Brightness::EmissionLineCount+line;
}

bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}

bool NormalizeDirection(const double input[3],double output[3],
    const char *name,std::string *error) {
  if (input==NULL || output==NULL)
    return Fail(std::string(name)+" direction pointer is null",error);
  double norm2=0.0;
  for (int i=0;i<3;i++) {
    if (!std::isfinite(input[i]))
      return Fail(std::string(name)+" direction is not finite",error);
    norm2+=input[i]*input[i];
  }
  if (!(norm2>0.0) || !std::isfinite(norm2))
    return Fail(std::string(name)+" direction is zero",error);
  const double inverseNorm=1.0/std::sqrt(norm2);
  for (int i=0;i<3;i++) output[i]=input[i]*inverseNorm;
  return true;
}

const AMPS::ResonanceFluorescence::ModelMetadata *CoreMetadata(int mode) {
  using AMPS::ResonanceFluorescence::Model;
  using AMPS::ResonanceFluorescence::GetMetadata;
  switch (mode) {
  case Comet::Brightness::GFactor::HydrogenLymanAlphaIKL2020:
    return GetMetadata(Model::HydrogenLymanAlphaIKL2020);
  case Comet::Brightness::GFactor::HydroxylAX00SchleicherAhearn1988Unquenched:
    return GetMetadata(Model::HydroxylAX00SchleicherAhearn1988Unquenched);
  default: return NULL;
  }
}

bool ModelForLine(int line,AMPS::ResonanceFluorescence::Model &model) {
  using AMPS::ResonanceFluorescence::Model;
  switch (Comet::Brightness::GFactorModel[line]) {
  case Comet::Brightness::GFactor::HydrogenLymanAlphaIKL2020:
    model=Model::HydrogenLymanAlphaIKL2020;
    return true;
  case Comet::Brightness::GFactor::HydroxylAX00SchleicherAhearn1988Unquenched:
    model=Model::HydroxylAX00SchleicherAhearn1988Unquenched;
    return true;
  default: return false;
  }
}

AMPS::ResonanceFluorescence::Status EvaluatePhotonGFactor(int line,
    double radialVelocity,double &photonRate) {
  using namespace AMPS::ResonanceFluorescence;
  if (Comet::Brightness::GFactorModel[line]==
      Comet::Brightness::GFactor::Fixed) {
    photonRate=Comet::Brightness::PhotonGFactor[line];
    return Status::Ok;
  }

  Model model;
  if (!ModelForLine(line,model)) return Status::UnknownModel;
  const Conditions conditions={radialVelocity,
      Comet::Model::HeliocentricDistanceAU*AstronomicalUnitM,
      Comet::Brightness::SolarLymanAlphaFlux1AU};
  const Result result=Evaluate(model,conditions);
  if (result.status!=Status::Ok) return result.status;
  if (model==Model::HydrogenLymanAlphaIKL2020) {
    photonRate=result.photonRateS1;
  }
  else {
    photonRate=Comet::Brightness::Math::EnergyRateToPhotonRate(
        result.emittedPowerW,
        Comet::Brightness::EffectivePhotonWavelength[line]);
  }
  return (!std::isfinite(photonRate) || !(photonRate>0.0)) ?
      Status::InvalidModelResult : Status::Ok;
}

void SetGeneratedMetadata(int line,const char *transition,
    const char *provenance,const char *tauProvenance) {
  if (line<0 || line>=Comet::Brightness::MaxEmissionLines ||
      transition==NULL || provenance==NULL || tauProvenance==NULL ||
      std::strlen(transition)>=Comet::Brightness::MetadataStringLength ||
      std::strlen(provenance)>=Comet::Brightness::MetadataStringLength ||
      std::strlen(tauProvenance)>=Comet::Brightness::MetadataStringLength) {
    std::fprintf(stderr,"invalid generated comet-brightness metadata\n");
    std::abort();
  }
  std::strcpy(TransitionName[line],transition);
  std::strcpy(GFactorProvenance[line],provenance);
  std::strcpy(TauProvenance[line],tauProvenance);
}

void LoadGeneratedMetadata() {
  if (MetadataLoaded) return;
  MetadataLoaded=true;
  // COMET_BRIGHTNESS_INPUT_METADATA_MARKER -- generated calls follow.
}

double ImpactParameter(long int point) {
  if (Comet::Brightness::ImpactParameterCount==1)
    return Comet::Brightness::MinimumImpactParameter;
  const double fraction=static_cast<double>(point)/
      static_cast<double>(Comet::Brightness::ImpactParameterCount-1);
  if (Comet::Brightness::ImpactParameterGridMode==
      Comet::Brightness::GridMode::Logarithmic) {
    return std::exp(std::log(Comet::Brightness::MinimumImpactParameter)+
        fraction*(std::log(Comet::Brightness::MaximumImpactParameter)-
        std::log(Comet::Brightness::MinimumImpactParameter)));
  }
  return Comet::Brightness::MinimumImpactParameter+
      fraction*(Comet::Brightness::MaximumImpactParameter-
      Comet::Brightness::MinimumImpactParameter);
}

struct cParticleSample {
  int Species;
  // Number of physical particles represented by this model particle.  This
  // includes both the block-local base weight and, when enabled, the
  // individual correction.  Dividing by the exact spherical cell measure
  // converts it to the piecewise-constant DIM1 density estimator.
  double RepresentedPopulation;
  // AMPS DIM1 retains all three velocity components in the local radial
  // basis even though it stores position as (r,0,0).  The brightness sampler
  // must rotate this vector to the actual LOS quadrature point before taking
  // a heliocentric projection.
  double LocalVelocity[3];
};

void SampleQuadraturePoint(const std::vector<cParticleSample> &particles,
    long int impact,const cRayGeometry &ray,double lineCoordinate,
    double pathWeight,double inverseCellMeasure) {
  double position[3],radialDirection[3];
  double radius2=0.0;
  for (int i=0;i<3;i++) {
    // LOS convention: positive lineCoordinate points from the comet toward
    // the observer.  ray.SkyPlaneOffset is b*IMPACT and is perpendicular to
    // LOS by configuration validation.
    position[i]=ray.SkyPlaneOffset[i]+lineCoordinate*
        Comet::Brightness::LineOfSightDirection[i];
    radius2+=position[i]*position[i];
  }
  if (!(radius2>0.0) || !std::isfinite(radius2)) {
    std::fprintf(stderr,"invalid comet-brightness quadrature position\n");
    MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
    std::abort();
  }
  const double inverseRadius=1.0/std::sqrt(radius2);
  for (int i=0;i<3;i++) radialDirection[i]=position[i]*inverseRadius;

  for (std::size_t particle=0;particle<particles.size();particle++) {
    const cParticleSample &sample=particles[particle];
    double cartesianVelocity[3];
    if (!Comet::Brightness::Math::RotateLocalVelocityToCartesian(
            sample.LocalVelocity,radialDirection,cartesianVelocity)) {
      std::fprintf(stderr,"invalid velocity in comet-brightness sampler\n");
      MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
      std::abort();
    }
    double heliocentricVelocity[3];
    for (int i=0;i<3;i++) heliocentricVelocity[i]=
        Comet::Brightness::CometHeliocentricVelocity[i]+
        cartesianVelocity[i];
    const double radialVelocity=Comet::Brightness::Math::Dot3(
        heliocentricVelocity,
        Comet::Brightness::HeliocentricRadialDirection);

    for (int line=0;line<Comet::Brightness::EmissionLineCount;line++) {
      if (Comet::Brightness::EmissionSpecies[line]!=sample.Species ||
          Comet::Brightness::GFactorModel[line]==
              Comet::Brightness::GFactor::Fixed) continue;
      const std::size_t index=SampleIndex(impact,line);
      double photonRate=0.0;
      const AMPS::ResonanceFluorescence::Status status=
          EvaluatePhotonGFactor(line,radialVelocity,photonRate);
      const double columnContribution=sample.RepresentedPopulation*
          inverseCellMeasure*pathWeight;
      if (status!=AMPS::ResonanceFluorescence::Status::Ok) {
        ExcludedColumnDensity[index]+=columnContribution;
        ExcludedEvaluationCount[index]++;
        if (Comet::Brightness::OutOfRangeVelocityPolicy==
            Comet::Brightness::VelocityDomainPolicy::Abort) {
          const AMPS::ResonanceFluorescence::ModelMetadata *metadata=
              CoreMetadata(Comet::Brightness::GFactorModel[line]);
          std::fprintf(stderr,
              "brightness g-factor evaluation failed on rank %d: species=%d transition=%s model=%s velocity=%.17e m/s status=%s\n",
              PIC::ThisThread,sample.Species,TransitionName[line],
              metadata==NULL ? "UNKNOWN" : metadata->id,radialVelocity,
              AMPS::ResonanceFluorescence::StatusMessage(status));
          MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
          std::abort();
        }
        continue;
      }
      SampledBrightnessRayleigh[index]+=
          columnContribution*photonRate/1.0e10;
    }
  }
}

void SampleSegment(const std::vector<cParticleSample> &particles,
    long int impact,const cRayGeometry &ray,double begin,double end,
    double inverseCellMeasure) {
  if (!(end>begin)) return;
  // The two points +/-1/sqrt(3) integrate a quadratic exactly.  More
  // importantly here, they evaluate the Doppler-sensitive g factor at two
  // physical radial directions instead of assigning a whole curved shell
  // segment the midpoint direction.  halfLength is also each quadrature
  // point's path weight, so the two weights sum to the segment length.
  const double midpoint=0.5*(begin+end);
  const double halfLength=0.5*(end-begin);
  SampleQuadraturePoint(particles,impact,ray,
      midpoint-halfLength*InverseSqrtThree,halfLength,inverseCellMeasure);
  SampleQuadraturePoint(particles,impact,ray,
      midpoint+halfLength*InverseSqrtThree,halfLength,inverseCellMeasure);
}

std::string TecplotQuoted(const char *value) {
  std::string result="\"";
  if (value!=NULL) for (const char *c=value;*c!='\0';c++) {
    if (*c=='\"') result+='\"';
    result+=*c;
  }
  result+='\"';
  return result;
}

void ClearSamplingBuffers() {
  std::fill(SampledColumnDensity.begin(),SampledColumnDensity.end(),0.0);
  std::fill(SampledBrightnessRayleigh.begin(),
      SampledBrightnessRayleigh.end(),0.0);
  std::fill(ExcludedColumnDensity.begin(),ExcludedColumnDensity.end(),0.0);
  std::fill(ExcludedEvaluationCount.begin(),
      ExcludedEvaluationCount.end(),0ULL);
}
} // namespace

bool Comet::Brightness::SetCometHeliocentricVelocity(
    const double velocity[3],std::string *error) {
  if (velocity==NULL) return Fail("comet heliocentric velocity is null",error);
  for (int i=0;i<3;i++) {
    if (!std::isfinite(velocity[i]))
      return Fail("comet heliocentric velocity is not finite",error);
    CometHeliocentricVelocity[i]=velocity[i];
  }
  BulkHeliocentricRadialVelocity=Math::Dot3(CometHeliocentricVelocity,
      HeliocentricRadialDirection);
  return true;
}

bool Comet::Brightness::SetObservationGeometry(
    const double radialDirection[3],const double losDirection[3],
    const double impactDirection[3],std::string *error) {
  double radial[3],los[3],impact[3];
  if (!NormalizeDirection(radialDirection,radial,"heliocentric radial",error) ||
      !NormalizeDirection(losDirection,los,"line of sight",error) ||
      !NormalizeDirection(impactDirection,impact,"impact parameter",error))
    return false;
  if (std::fabs(Math::Dot3(los,impact))>1.0e-12)
    return Fail("brightness LOS and impact directions must be orthogonal",error);
  for (int i=0;i<3;i++) {
    HeliocentricRadialDirection[i]=radial[i];
    LineOfSightDirection[i]=los[i];
    ImpactParameterDirection[i]=impact[i];
  }
  ObservationGeometryConfigured=true;
  BulkHeliocentricRadialVelocity=Math::Dot3(CometHeliocentricVelocity,
      HeliocentricRadialDirection);
  return true;
}

bool Comet::Brightness::ResetRuntimeConfiguration(std::string *error) {
  if (SamplingInitialized || SamplingRegistered)
    return Fail("brightness configuration cannot change after sampler initialization",error);
  Enabled=false;
  MinimumImpactParameter=-1.0;
  MaximumImpactParameter=-1.0;
  ImpactParameterCount=0;
  ImpactParameterGridMode=GridMode::Linear;
  EmissionLineCount=0;
  SolarLymanAlphaFlux1AU=0.0;
  OutOfRangeVelocityPolicy=VelocityDomainPolicy::Abort;
  ObservationGeometryConfigured=false;
  BulkHeliocentricRadialVelocity=0.0;
  for (int component=0;component<3;component++) {
    CometHeliocentricVelocity[component]=0.0;
    HeliocentricRadialDirection[component]=component==0 ? 1.0 : 0.0;
    LineOfSightDirection[component]=component==2 ? 1.0 : 0.0;
    ImpactParameterDirection[component]=component==0 ? 1.0 : 0.0;
  }
  for (int line=0;line<MaxEmissionLines;line++) {
    EmissionSpecies[line]=0;
    PhotonGFactor[line]=0.0;
    GFactorModel[line]=GFactor::Fixed;
    EffectivePhotonWavelength[line]=0.0;
    TauDiagnosticMode[line]=OpticalDepthDiagnostic::None;
    LineCenterCrossSection[line]=0.0;
    TauFlagThreshold[line]=0.0;
    TransitionName[line][0]='\0';
    GFactorProvenance[line][0]='\0';
    TauProvenance[line][0]='\0';
  }
  MetadataLoaded=true;  // Runtime metadata replaces generated metadata.
  return true;
}

bool Comet::Brightness::SetEmissionMetadata(int line,const char *transition,
    const char *gFactorProvenance,const char *tauProvenance,
    std::string *error) {
  if (line<0 || line>=MaxEmissionLines || transition==NULL ||
      gFactorProvenance==NULL || tauProvenance==NULL)
    return Fail("invalid brightness emission metadata arguments",error);
  if (std::strlen(transition)>=MetadataStringLength ||
      std::strlen(gFactorProvenance)>=MetadataStringLength ||
      std::strlen(tauProvenance)>=MetadataStringLength)
    return Fail("brightness emission metadata exceeds compiled storage",error);
  std::strcpy(TransitionName[line],transition);
  std::strcpy(GFactorProvenance[line],gFactorProvenance);
  std::strcpy(TauProvenance[line],tauProvenance);
  MetadataLoaded=true;
  return true;
}

bool Comet::Brightness::Validate(std::string *error) {
  LoadGeneratedMetadata();
  if (!Enabled) return true;
#if DIM != 1
  return Fail("comet brightness is currently implemented only for DIM1",error);
#else
  if (!ObservationGeometryConfigured)
    return Fail("particle-resolved brightness requires explicit observation geometry",error);
  double radial[3],los[3],impact[3];
  if (!NormalizeDirection(HeliocentricRadialDirection,radial,
          "heliocentric radial",error) ||
      !NormalizeDirection(LineOfSightDirection,los,"line of sight",error) ||
      !NormalizeDirection(ImpactParameterDirection,impact,
          "impact parameter",error)) return false;
  for (int i=0;i<3;i++) {
    if (!std::isfinite(CometHeliocentricVelocity[i]))
      return Fail("comet heliocentric velocity is not finite",error);
    HeliocentricRadialDirection[i]=radial[i];
    LineOfSightDirection[i]=los[i];
    ImpactParameterDirection[i]=impact[i];
  }
  if (std::fabs(Math::Dot3(LineOfSightDirection,
          ImpactParameterDirection))>1.0e-12)
    return Fail("brightness LOS and impact directions must be orthogonal",error);
  BulkHeliocentricRadialVelocity=Math::Dot3(CometHeliocentricVelocity,
      HeliocentricRadialDirection);
  if (OutOfRangeVelocityPolicy!=VelocityDomainPolicy::Abort &&
      OutOfRangeVelocityPolicy!=VelocityDomainPolicy::ExcludeAndReport)
    return Fail("unknown brightness velocity-domain policy",error);
  if (EmissionLineCount<=0 || EmissionLineCount>MaxEmissionLines)
    return Fail("brightness requires 1..16 configured emission lines",error);
  if (ImpactParameterCount<=0)
    return Fail("brightness impact-parameter count must be positive",error);
  if (!std::isfinite(MinimumImpactParameter) ||
      !std::isfinite(MaximumImpactParameter) ||
      MinimumImpactParameter<0.0 ||
      MaximumImpactParameter<MinimumImpactParameter ||
      (ImpactParameterCount>1 &&
       !(MaximumImpactParameter>MinimumImpactParameter)))
    return Fail("invalid brightness impact-parameter range",error);
  if (MaximumImpactParameter>Comet::RadialModel::GetDomainOuterRadius())
    return Fail("brightness impact grid extends beyond the radial domain",error);
  if (ImpactParameterGridMode!=GridMode::Linear &&
      ImpactParameterGridMode!=GridMode::Logarithmic)
    return Fail("unknown brightness impact-parameter grid",error);
  if (ImpactParameterGridMode==GridMode::Logarithmic &&
      !(MinimumImpactParameter>0.0))
    return Fail("a logarithmic brightness grid requires positive minimum impact parameter",error);
  bool usesCoreModel=false;
  for (int line=0;line<EmissionLineCount;line++)
    usesCoreModel=usesCoreModel || GFactorModel[line]!=GFactor::Fixed;
  if (usesCoreModel) {
    const char *referenceError=NULL;
    if (!AMPS::ResonanceFluorescence::ValidateReferenceData(&referenceError))
      return Fail(std::string("invalid resonance-fluorescence reference data: ")+
          (referenceError==NULL ? "unknown error" : referenceError),error);
  }
  for (int line=0;line<EmissionLineCount;line++) {
    if (GFactorModel[line]==GFactor::Fixed) {
      if (!std::isfinite(PhotonGFactor[line]) || !(PhotonGFactor[line]>0.0) ||
          EffectivePhotonWavelength[line]!=0.0)
        return Fail("invalid fixed brightness photon g factor",error);
    }
    else {
      double nominalPhotonRate=0.0;
      const AMPS::ResonanceFluorescence::Status status=
          EvaluatePhotonGFactor(line,BulkHeliocentricRadialVelocity,
              nominalPhotonRate);
      if (status!=AMPS::ResonanceFluorescence::Status::Ok)
        return Fail(std::string("brightness model bulk-state validation failed: ")+
            AMPS::ResonanceFluorescence::StatusMessage(status),error);
      PhotonGFactor[line]=nominalPhotonRate;
    }
    if (EmissionSpecies[line]<0 || EmissionSpecies[line]>=PIC::nTotalSpecies ||
        TransitionName[line][0]=='\0' || GFactorProvenance[line][0]=='\0' ||
        TauProvenance[line][0]=='\0')
      return Fail("invalid generated brightness emission line",error);
    if (TauDiagnosticMode[line]==OpticalDepthDiagnostic::None) {
      if (LineCenterCrossSection[line]!=0.0 || TauFlagThreshold[line]!=0.0 ||
          std::strcmp(TauProvenance[line],"NONE")!=0)
        return Fail("brightness tau NONE requires zero values and provenance NONE",error);
    }
    else if (TauDiagnosticMode[line]==OpticalDepthDiagnostic::LineCenter) {
      if (!std::isfinite(LineCenterCrossSection[line]) ||
          !(LineCenterCrossSection[line]>0.0) ||
          !std::isfinite(TauFlagThreshold[line]) ||
          !(TauFlagThreshold[line]>0.0) ||
          std::strcmp(TauProvenance[line],"NONE")==0)
        return Fail("line-center tau diagnostics require positive values and provenance",error);
    }
    else return Fail("unknown brightness optical-depth diagnostic",error);
  }
  return true;
#endif
}

void Comet::Brightness::InitializeSampling() {
  if (!Enabled) return;
#if DIM == 1
  if (SamplingInitialized) return;
  const std::size_t sampleCount=static_cast<std::size_t>(ImpactParameterCount)*
      EmissionLineCount;
  if (sampleCount==0 || sampleCount>static_cast<std::size_t>(INT_MAX)) {
    std::fprintf(stderr,"comet-brightness sampling buffer exceeds MPI limits\n");
    MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
    std::abort();
  }
  RayGeometry.resize(ImpactParameterCount);
  for (long int point=0;point<ImpactParameterCount;point++) {
    cRayGeometry &ray=RayGeometry[point];
    ray.ImpactParameter=ImpactParameter(point);
    ray.ImpactParameterSquared=ray.ImpactParameter*ray.ImpactParameter;
    for (int i=0;i<3;i++) ray.SkyPlaneOffset[i]=ray.ImpactParameter*
        ImpactParameterDirection[i];
  }
  SampledColumnDensity.assign(sampleCount,0.0);
  SampledBrightnessRayleigh.assign(sampleCount,0.0);
  ExcludedColumnDensity.assign(sampleCount,0.0);
  ExcludedEvaluationCount.assign(sampleCount,0ULL);
  SamplingInitialized=true;
  if (!SamplingRegistered) {
    PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(
        SampleParticleResolvedBrightness,OutputSampledBrightness);
    SamplingRegistered=true;
  }
#endif
}

void Comet::Brightness::SampleParticleResolvedBrightness() {
  if (!Enabled) return;
#if DIM == 1
  if (!SamplingInitialized) {
    std::fprintf(stderr,"comet-brightness sampler was not initialized\n");
    MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
    std::abort();
  }

  const double nucleusRadius=Comet::RadialModel::NucleusRadius;
  for (long int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL || node->Thread!=PIC::ThisThread)
      continue;
    const double cellWidth=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
    double blockEvaluations=0.0;
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    const double blockStartTime=MPI_Wtime();
#endif
    for (int localCell=0;localCell<_BLOCK_CELLS_X_;localCell++) {
      // These two AMPS arrays deliberately use different index spaces.
      // GetCenterNode() addresses the center-node array, which includes
      // ghost cells; in this DIM1 build its index is localCell+
      // _GHOST_CELLS_X_.  FirstCellParticleTable contains only the physical
      // cells and is therefore indexed densely from 0 to
      // _BLOCK_CELLS_X_-1.  Reusing the ghost-aware center-node index for
      // the particle table reads past its end.  In an empty simulation that
      // adjacent value can be zero, which is the head of the particle-buffer
      // free list, and the sampler will then mistake the complete free list
      // for active particles.  Keep the names and the bounds check separate
      // so a future multidimensional extension must define its own dense
      // physical-cell mapping explicitly.
      const long int centerNodeNumber=
          _getCenterNodeLocalNumber(localCell,0,0);
      const long int particleCellNumber=localCell;
      if (particleCellNumber<0 || particleCellNumber>=_BLOCK_CELLS_X_) {
        std::fprintf(stderr,
            "invalid physical-cell index in comet-brightness sampler\n");
        MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
        std::abort();
      }
      PIC::Mesh::cDataCenterNode *cell=
          node->block->GetCenterNode(centerNodeNumber);
      if (cell==NULL || !(cell->Measure>0.0) ||
          !std::isfinite(cell->Measure)) {
        std::fprintf(stderr,"invalid shell measure in comet-brightness sampler\n");
        MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
        std::abort();
      }
      const double innerRadius=node->xmin[0]+localCell*cellWidth;
      const double outerRadius=innerRadius+cellWidth;
      const double inverseCellMeasure=1.0/cell->Measure;
      std::vector<cParticleSample> particles;
      double representedPopulation[MaxEmissionLines]={0.0};

      for (long int ptr=
               node->block->FirstCellParticleTable[particleCellNumber];
           ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        const int species=PIC::ParticleBuffer::GetI(data);
        bool emitted=false;
        for (int line=0;line<EmissionLineCount;line++)
          emitted=emitted || EmissionSpecies[line]==species;
        if (!emitted) continue;
        cParticleSample sample;
        sample.Species=species;
        sample.RepresentedPopulation=
            node->block->GetLocalParticleWeight(species);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        sample.RepresentedPopulation*=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
#endif
        PIC::ParticleBuffer::GetV(sample.LocalVelocity,data);
        if (!(sample.RepresentedPopulation>0.0) ||
            !std::isfinite(sample.RepresentedPopulation)) {
          std::fprintf(stderr,"invalid particle weight in brightness sampler\n");
          MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
          std::abort();
        }
        particles.push_back(sample);
        for (int line=0;line<EmissionLineCount;line++)
          if (EmissionSpecies[line]==species)
            representedPopulation[line]+=sample.RepresentedPopulation;
      }
      if (particles.empty()) continue;

      for (long int impact=0;impact<ImpactParameterCount;impact++) {
        const cRayGeometry &ray=RayGeometry[impact];
        const double b=ray.ImpactParameter;
        if (!(b<outerRadius)) continue;
        const double clippedInner=std::max(innerRadius,nucleusRadius);
        if (!(outerRadius>clippedInner)) continue;
        const double innerAlongRay=std::max(clippedInner,b);
        const double zOuter=std::sqrt(std::max(0.0,
            outerRadius*outerRadius-ray.ImpactParameterSquared));
        const double zInner=std::sqrt(std::max(0.0,
            innerAlongRay*innerAlongRay-ray.ImpactParameterSquared));
        if (!(zOuter>zInner)) continue;

        // A ray outside or tangent to the opaque nucleus sees both sides of
        // the shell.  For b<Rn the far-side interval is hidden; the near
        // side (positive LOS coordinate, toward the observer) remains.
        const bool farSideVisible=b>=nucleusRadius;
        const double visibleLength=(farSideVisible ? 2.0 : 1.0)*
            (zOuter-zInner);
        for (int line=0;line<EmissionLineCount;line++) {
          const std::size_t index=SampleIndex(impact,line);
          const double column=representedPopulation[line]*
              inverseCellMeasure*visibleLength;
          SampledColumnDensity[index]+=column;
          if (GFactorModel[line]==GFactor::Fixed)
            SampledBrightnessRayleigh[index]+=
                column*PhotonGFactor[line]/1.0e10;
        }
        SampleSegment(particles,impact,ray,zInner,zOuter,inverseCellMeasure);
        if (farSideVisible)
          SampleSegment(particles,impact,ray,-zOuter,-zInner,
              inverseCellMeasure);
        blockEvaluations+=(farSideVisible ? 4.0 : 2.0)*
            static_cast<double>(particles.size());
      }
    }
    // See LoadPerGFactorEvaluation above.  Only owner blocks reach this
    // point, and the core consumes and resets the measure at every
    // CreateNewParallelDistributionLists() call.
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    node->ParallelLoadMeasure+=LoadPerGFactorEvaluation*blockEvaluations;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    (void)blockEvaluations;
    node->ParallelLoadMeasure+=MPI_Wtime()-blockStartTime;
#else
    (void)blockEvaluations;
#endif
  }
#endif
}

void Comet::Brightness::OutputSampledBrightness(int outputIndex) {
  if (!Enabled) return;
#if DIM == 1
  if (!SamplingInitialized) return;
  const std::size_t count=SampledColumnDensity.size();
  std::vector<double> globalColumn(count,0.0),globalBrightness(count,0.0),
      globalExcludedColumn(count,0.0);
  std::vector<unsigned long long> globalExcludedCount(count,0ULL);
  MPI_Reduce(SampledColumnDensity.data(),globalColumn.data(),
      static_cast<int>(count),MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(SampledBrightnessRayleigh.data(),globalBrightness.data(),
      static_cast<int>(count),MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(ExcludedColumnDensity.data(),globalExcludedColumn.data(),
      static_cast<int>(count),MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(ExcludedEvaluationCount.data(),globalExcludedCount.data(),
      static_cast<int>(count),MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    const double inverseSamples=PIC::LastSampleLength>0 ?
        1.0/PIC::LastSampleLength : 0.0;
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      std::vector<int> lines;
      for (int line=0;line<EmissionLineCount;line++)
        if (EmissionSpecies[line]==species) lines.push_back(line);
      if (lines.empty()) continue;
      char symbol[_MAX_STRING_LENGTH_PIC_];
      PIC::MolecularData::GetChemSymbol(symbol,species);
      const std::string fileName=std::string(PIC::OutputDataFileDirectory)+
          "/pic."+symbol+".brightness.s="+std::to_string(species)+
          ".out="+std::to_string(outputIndex)+".dat";
      FILE *output=std::fopen(fileName.c_str(),"w");
      if (output==NULL) {
        std::fprintf(stderr,"cannot open comet-brightness output %s: %s\n",
            fileName.c_str(),std::strerror(errno));
        MPI_Abort(MPI_GLOBAL_COMMUNICATOR,EXIT_FAILURE);
        std::abort();
      }
      std::fprintf(output,
          "TITLE=\"Particle-resolved comet brightness: %s; output=%d; samples=%ld\"\n",
          symbol,outputIndex,PIC::LastSampleLength);
      std::fprintf(output,"VARIABLES=\"Impact Parameter [m]\"");
      for (std::size_t entry=0;entry<lines.size();entry++) {
        const int line=lines[entry];
        std::fprintf(output,
            ", \"Column Density %s [m^-2]\", \"Effective Photon g Factor %s [s^-1]\", \"Brightness %s [R]\", \"Brightness %s [kR]\", \"Line Center Tau %s\", \"Tau Threshold Exceeded %s\", \"Excluded Column Fraction %s\", \"Excluded g Evaluations %s\"",
            TransitionName[line],TransitionName[line],TransitionName[line],
            TransitionName[line],TransitionName[line],TransitionName[line],
            TransitionName[line],TransitionName[line]);
      }
      std::fprintf(output,"\n");
      std::fprintf(output,
          "DATASETAUXDATA HeliocentricDistanceAU=\"%.17e\"\n",
          Comet::Model::HeliocentricDistanceAU);
      std::fprintf(output,
          "DATASETAUXDATA SamplingEstimator=\"PARTICLE_RESOLVED_DIM1_SHELL_GAUSS2\"\n");
      std::fprintf(output,
          "DATASETAUXDATA VelocityDomainPolicy=\"%s\"\n",
          OutOfRangeVelocityPolicy==VelocityDomainPolicy::Abort ?
              "ABORT" : "EXCLUDE_AND_REPORT");
      std::fprintf(output,
          "DATASETAUXDATA CometHeliocentricVelocityMPS=\"%.17e %.17e %.17e\"\n",
          CometHeliocentricVelocity[0],CometHeliocentricVelocity[1],
          CometHeliocentricVelocity[2]);
      std::fprintf(output,
          "DATASETAUXDATA HeliocentricRadialDirection=\"%.17e %.17e %.17e\"\n",
          HeliocentricRadialDirection[0],HeliocentricRadialDirection[1],
          HeliocentricRadialDirection[2]);
      std::fprintf(output,
          "DATASETAUXDATA LineOfSightDirectionTowardObserver=\"%.17e %.17e %.17e\"\n",
          LineOfSightDirection[0],LineOfSightDirection[1],
          LineOfSightDirection[2]);
      std::fprintf(output,
          "DATASETAUXDATA ImpactParameterDirection=\"%.17e %.17e %.17e\"\n",
          ImpactParameterDirection[0],ImpactParameterDirection[1],
          ImpactParameterDirection[2]);
      for (std::size_t entry=0;entry<lines.size();entry++) {
        const int line=lines[entry];
        const AMPS::ResonanceFluorescence::ModelMetadata *metadata=
            CoreMetadata(GFactorModel[line]);
        std::fprintf(output,
            "DATASETAUXDATA Line%zu=\"transition=%s; model=%s; provenance=%s; reference_sha256=%s; tau_provenance=%s\"\n",
            entry,TransitionName[line],metadata==NULL ? "FIXED" : metadata->id,
            GFactorProvenance[line],
            metadata==NULL ? "" : metadata->referenceDataSha256,
            TauProvenance[line]);
      }
      std::fprintf(output,
          "ZONE T=%s, I=%ld, DATAPACKING=POINT\n",
          TecplotQuoted(symbol).c_str(),ImpactParameterCount);
      for (long int impact=0;impact<ImpactParameterCount;impact++) {
        std::fprintf(output,"%.17e",RayGeometry[impact].ImpactParameter);
        for (std::size_t entry=0;entry<lines.size();entry++) {
          const int line=lines[entry];
          const std::size_t index=SampleIndex(impact,line);
          const double column=globalColumn[index]*inverseSamples;
          const double brightness=globalBrightness[index]*inverseSamples;
          const double excluded=globalExcludedColumn[index]*inverseSamples;
          const double effectiveG=column>0.0 ?
              brightness*1.0e10/column : 0.0;
          const double tau=TauDiagnosticMode[line]==
              OpticalDepthDiagnostic::LineCenter ?
              Math::LineCenterOpticalDepth(column,
                  LineCenterCrossSection[line]) : 0.0;
          const bool tauFlag=TauDiagnosticMode[line]==
              OpticalDepthDiagnostic::LineCenter &&
              Math::OpticalDepthThresholdExceeded(tau,
                  TauFlagThreshold[line]);
          const double excludedFraction=column>0.0 ? excluded/column : 0.0;
          std::fprintf(output,
              " %.17e %.17e %.17e %.17e %.17e %d %.17e %llu",
              column,effectiveG,brightness,brightness/1.0e3,tau,
              tauFlag ? 1 : 0,excludedFraction,
              globalExcludedCount[index]);
        }
        std::fprintf(output,"\n");
      }
      std::fclose(output);
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: wrote particle-resolved comet brightness: %s\n",
          fileName.c_str());
      std::fflush(PIC::DiagnospticMessageStream);
    }
  }

  ClearSamplingBuffers();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
#else
  (void)outputIndex;
#endif
}
