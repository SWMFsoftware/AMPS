#include "resonance_fluorescence.h"

#include <cmath>
#include <cstring>
#include <limits>

namespace AMPS {
namespace ResonanceFluorescence {
namespace {

constexpr double AstronomicalUnitM=1.495978707e11;
constexpr double PlanckJS=6.62607015e-34;
constexpr double LightSpeedMPS=299792458.0;
constexpr double LymanAlphaWavelengthM=121.56701e-9;
constexpr double SolarGMM3S2=1.327124421864553e20;
constexpr double HydrogenMassKG=1.6735328e-27;
constexpr double IKLMeanIrradiancePhotonCM2S=4.25e11;
constexpr double ErgPerSecondToWatt=1.0e-7;

struct IKLCoefficient {
  const char *name;
  double beta;
  double alpha;
};

#include "generated/reference_tables.inc"

constexpr std::size_t HCoefficientCount=
    sizeof(HIKL2020Coefficients)/sizeof(HIKL2020Coefficients[0]);
constexpr std::size_t OHCompiledLength=
    sizeof(OHEnergyGFactor1AUErgPerSecond)/
    sizeof(OHEnergyGFactor1AUErgPerSecond[0]);
constexpr double OHVelocityMaximumMPS=
    OHVelocityMinimumMPS+(OHTableLength-1)*OHVelocityStepMPS;

const ModelMetadata HMetadata={
  Model::HydrogenLymanAlphaIKL2020,
  "H_LYA_IKL2020",
  "H",
  "LYMAN_ALPHA_121.56701_NM",
  "Kowalska-Leszczynska et al. 2020, ApJS 247:62",
  "https://users.cbk.waw.pl/~ikowalska/LyaModel/files/table1v4.txt",
  HReferenceSha256,
  AstronomicalUnitM,-200.0e3,200.0e3,true,true,true,true
};

const ModelMetadata OHMetadata={
  Model::HydroxylAX00SchleicherAhearn1988Unquenched,
  "OH_AX00_SA1988_LOWELL_UNQUENCHED",
  "OH",
  "A2SIGMA_X2PI_00",
  "Schleicher and A'Hearn 1988, ApJ 331:1058",
  "https://asteroid.lowell.edu/api/comet/calc/gfactor",
  OHReferenceSha256,
  AstronomicalUnitM,OHVelocityMinimumMPS,OHVelocityMaximumMPS,
  false,false,true,false
};

double MissingValue() {
  return std::numeric_limits<double>::quiet_NaN();
}

Result Failure(Status status) {
  return Result{status,MissingValue(),MissingValue(),MissingValue()};
}

double DistanceScale(double heliocentricDistanceM) {
  const double ratio=AstronomicalUnitM/heliocentricDistanceM;
  return ratio*ratio;
}

double HParameter(std::size_t index,double irradiancePhotonCM2S) {
  const IKLCoefficient &coefficient=HIKL2020Coefficients[index];
  return coefficient.beta*(1.0+coefficient.alpha*
      irradiancePhotonCM2S/IKLMeanIrradiancePhotonCM2S);
}

Result EvaluateH(const Conditions &conditions) {
  const double irradianceCM2=
      conditions.solarLymanAlphaFlux1AUPhotonM2S/1.0e4;
  if (!std::isfinite(irradianceCM2) || !(irradianceCM2>0.0))
    return Failure(Status::MissingSolarLymanAlphaFlux);

  const double velocityKMS=conditions.radialVelocityMPS/1.0e3;
  const double ak=HParameter(0,irradianceCM2);
  const double muk=HParameter(1,irradianceCM2);
  const double sigmak=HParameter(2,irradianceCM2);
  const double kappa=HParameter(3,irradianceCM2);
  const double ar=HParameter(4,irradianceCM2);
  const double dmu=HParameter(5,irradianceCM2);
  const double sigmar=HParameter(6,irradianceCM2);
  const double bbkg=HParameter(7,irradianceCM2);
  const double abkg=HParameter(8,irradianceCM2);
  const double pi=std::acos(-1.0);

  if (!(sigmak>0.0) || !(sigmar>0.0) || !(kappa>0.0))
    return Failure(Status::InvalidModelResult);
  const double normalized=(velocityKMS-muk)*(velocityKMS-muk)/
      (2.0*sigmak*sigmak*kappa);
  const double kappaComponent=ak*std::pow(1.0+normalized,-kappa-1.0);
  const double reversalArgument=velocityKMS-(muk+dmu);
  const double reversal=ar/(sigmar*std::sqrt(2.0*pi))*
      std::exp(-reversalArgument*reversalArgument/(2.0*sigmar*sigmar));
  const double radiationPressureRatio=
      kappaComponent-reversal+abkg*velocityKMS+bbkg;
  if (!std::isfinite(radiationPressureRatio) ||
      !(radiationPressureRatio>0.0))
    return Failure(Status::InvalidModelResult);

  // A resonant scattering event transfers h*nu/c of radial photon momentum.
  // Equating that force to mu*G*M_sun*m_H/AU^2 converts the published IKL
  // radiation-pressure profile into the 1-AU photon scattering rate.
  const double frequency=LightSpeedMPS/LymanAlphaWavelengthM;
  const double photonMomentum=PlanckJS*frequency/LightSpeedMPS;
  const double solarGravityForce1AU=
      SolarGMM3S2*HydrogenMassKG/(AstronomicalUnitM*AstronomicalUnitM);
  const double photonRate1AU=
      radiationPressureRatio*solarGravityForce1AU/photonMomentum;
  const double photonRate=photonRate1AU*
      DistanceScale(conditions.heliocentricDistanceM);
  const double power=photonRate*PlanckJS*frequency;
  if (!std::isfinite(photonRate) || !(photonRate>0.0) ||
      !std::isfinite(power) || !(power>0.0))
    return Failure(Status::InvalidModelResult);
  return Result{Status::Ok,photonRate,power,radiationPressureRatio};
}

Result EvaluateOH(const Conditions &conditions) {
  const double coordinate=(conditions.radialVelocityMPS-
      OHVelocityMinimumMPS)/OHVelocityStepMPS;
  std::size_t lower=static_cast<std::size_t>(std::floor(coordinate));
  double energyErgPerSecond;
  if (lower>=OHTableLength-1) {
    // The public range check admits the exact upper endpoint only.  Handle it
    // directly so interpolation never addresses one element past the table.
    energyErgPerSecond=OHEnergyGFactor1AUErgPerSecond[OHTableLength-1];
  }
  else {
    const double fraction=coordinate-static_cast<double>(lower);
    energyErgPerSecond=OHEnergyGFactor1AUErgPerSecond[lower]+fraction*
        (OHEnergyGFactor1AUErgPerSecond[lower+1]-
         OHEnergyGFactor1AUErgPerSecond[lower]);
  }
  const double power=energyErgPerSecond*ErgPerSecondToWatt*
      DistanceScale(conditions.heliocentricDistanceM);
  if (!std::isfinite(power) || !(power>0.0))
    return Failure(Status::InvalidModelResult);
  return Result{Status::Ok,MissingValue(),power,MissingValue()};
}

} // namespace

bool ModelFromName(const char *name,Model *model) {
  if (name==NULL || model==NULL) return false;
  if (std::strcmp(name,HMetadata.id)==0) {
    *model=HMetadata.model;
    return true;
  }
  if (std::strcmp(name,OHMetadata.id)==0) {
    *model=OHMetadata.model;
    return true;
  }
  return false;
}

const char *ModelName(Model model) {
  const ModelMetadata *metadata=GetMetadata(model);
  return metadata==NULL ? "UNKNOWN" : metadata->id;
}

const char *StatusMessage(Status status) {
  switch (status) {
  case Status::Ok: return "ok";
  case Status::UnknownModel: return "unknown model";
  case Status::InvalidHeliocentricDistance:
    return "heliocentric distance must be finite and positive";
  case Status::InvalidRadialVelocity:
    return "heliocentric radial velocity must be finite";
  case Status::RadialVelocityOutOfRange:
    return "heliocentric radial velocity is outside the model range";
  case Status::MissingSolarLymanAlphaFlux:
    return "positive finite 1-AU solar Lyman-alpha flux is required";
  case Status::InvalidModelResult: return "model produced an invalid result";
  }
  return "invalid status";
}

const ModelMetadata *GetMetadata(Model model) {
  switch (model) {
  case Model::HydrogenLymanAlphaIKL2020: return &HMetadata;
  case Model::HydroxylAX00SchleicherAhearn1988Unquenched:
    return &OHMetadata;
  }
  return NULL;
}

Result Evaluate(Model model,const Conditions &conditions) {
  const ModelMetadata *metadata=GetMetadata(model);
  if (metadata==NULL) return Failure(Status::UnknownModel);
  if (!std::isfinite(conditions.heliocentricDistanceM) ||
      !(conditions.heliocentricDistanceM>0.0))
    return Failure(Status::InvalidHeliocentricDistance);
  if (!std::isfinite(conditions.radialVelocityMPS))
    return Failure(Status::InvalidRadialVelocity);
  if (conditions.radialVelocityMPS<metadata->minimumRadialVelocityMPS ||
      conditions.radialVelocityMPS>metadata->maximumRadialVelocityMPS)
    return Failure(Status::RadialVelocityOutOfRange);

  switch (model) {
  case Model::HydrogenLymanAlphaIKL2020: return EvaluateH(conditions);
  case Model::HydroxylAX00SchleicherAhearn1988Unquenched:
    return EvaluateOH(conditions);
  }
  return Failure(Status::UnknownModel);
}

bool ValidateReferenceData(const char **error) {
  if (error!=NULL) *error=NULL;
  if (HCoefficientCount!=9) {
    if (error!=NULL) *error="IKL-2020 coefficient table must contain nine rows";
    return false;
  }
  if (OHCompiledLength!=OHTableLength || OHTableLength!=1201 ||
      OHVelocityMinimumMPS!=-60.0e3 || OHVelocityMaximumMPS!=60.0e3 ||
      OHVelocityStepMPS!=100.0) {
    if (error!=NULL) *error="Lowell OH table geometry is invalid";
    return false;
  }
  for (std::size_t i=0;i<OHTableLength;i++) {
    if (!std::isfinite(OHEnergyGFactor1AUErgPerSecond[i]) ||
        !(OHEnergyGFactor1AUErgPerSecond[i]>0.0)) {
      if (error!=NULL) *error="Lowell OH table contains an invalid value";
      return false;
    }
  }
  if (std::strlen(HReferenceSha256)!=64 ||
      std::strlen(OHReferenceSha256)!=64) {
    if (error!=NULL) *error="reference-data SHA-256 metadata is invalid";
    return false;
  }
  return true;
}

} // namespace ResonanceFluorescence
} // namespace AMPS
