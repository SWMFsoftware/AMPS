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
constexpr std::size_t OKVB2022CompiledLength=
    sizeof(OKVB2022DirectPhotonRate1AUS1)/
    sizeof(OKVB2022DirectPhotonRate1AUS1[0]);
constexpr double OKVB2022VelocityMaximumMPS=
    OKVB2022VelocityMinimumMPS+(OKVB2022TableLength-1)*OKVB2022VelocityStepMPS;
constexpr std::size_t OSUMER2001DirectCompiledLength=
    sizeof(OSUMER2001DirectRatePerTripletFlux1AU)/
    sizeof(OSUMER2001DirectRatePerTripletFlux1AU[0]);
constexpr std::size_t OSUMER2001BowenCompiledLength=
    sizeof(OSUMER2001BowenRatePerLymanBetaFlux1AU)/
    sizeof(OSUMER2001BowenRatePerLymanBetaFlux1AU[0]);
constexpr double OSUMER2001VelocityMaximumMPS=
    OSUMER2001VelocityMinimumMPS+(OSUMER2001TableLength-1)*
    OSUMER2001VelocityStepMPS;

const ModelMetadata HMetadata={
  Model::HydrogenLymanAlphaIKL2020,
  "H_LYA_IKL2020",
  "H",
  "LYMAN_ALPHA_121.56701_NM",
  "Kowalska-Leszczynska et al. 2020, ApJS 247:62",
  "https://users.cbk.waw.pl/~ikowalska/LyaModel/files/table1v4.txt",
  HReferenceSha256,
  AstronomicalUnitM,-200.0e3,200.0e3,true,false,false,true,true,true
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
  false,false,false,false,true,false
};

// The three O I models share one transition label: every one of them counts
// photons of the complete 130.4-nm triplet emitted by 2p4 3P2 absorbers.
const ModelMetadata OKVB2022DirectMetadata={
  Model::OxygenI1304DirectKVB2022GroundJ2,
  "O_I_1304_KVB2022_DIRECT_J2",
  "O",
  "OI_130.4_NM_TRIPLET_FROM_3P2_DIRECT_130.217_NM_PUMP",
  "Killen, Vervack and Burger 2022, ApJS 263:37 (Fig. 14 data; "
  "Gladstone 1992 profiles, SORCE SOLSTICE 2009-03-01); NIST ASD",
  "https://content.cld.iop.org/journals/0067-0049/263/2/37/revision1/"
  "apjsac9eabf14.tar.gz",
  OKVB2022ReferenceSha256,
  AstronomicalUnitM,OKVB2022VelocityMinimumMPS,OKVB2022VelocityMaximumMPS,
  false,true,false,true,true,false
};

const ModelMetadata OSUMER2001DirectMetadata={
  Model::OxygenI1304DirectSUMER2001GroundJ2,
  "O_I_1304_SUMER2001_DIRECT_J2",
  "O",
  "OI_130.4_NM_TRIPLET_FROM_3P2_DIRECT_130.217_NM_PUMP",
  "Curdt et al. 2001, A&A 375:591 (SOHO/SUMER disk-centre atlas, "
  "BASS2000); Killen et al. 2022 line fraction; NIST ASD",
  "https://bass2000.obspm.fr/php/getSolarSpectrumDB.php?"
  "WL=1297&DW=14&resol=0.001&fmt=txt",
  OSUMER2001TripletReferenceSha256,
  AstronomicalUnitM,OSUMER2001VelocityMinimumMPS,OSUMER2001VelocityMaximumMPS,
  false,true,false,true,true,false
};

const ModelMetadata OSUMER2001BowenMetadata={
  Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2,
  "O_I_1304_SUMER2001_LYB_BOWEN_J2",
  "O",
  "OI_130.4_NM_TRIPLET_FROM_3P2_LYB_102.576_NM_PUMP_VIA_1128.7_844.6_NM",
  "Curdt et al. 2001, A&A 375:591 (SOHO/SUMER disk-centre atlas, "
  "BASS2000); NIST ASD; Bowen 1947, PASP 59:196",
  "https://bass2000.obspm.fr/php/getSolarSpectrumDB.php?"
  "WL=1019&DW=13&resol=0.001&fmt=txt",
  OSUMER2001LymanBetaReferenceSha256,
  AstronomicalUnitM,OSUMER2001VelocityMinimumMPS,OSUMER2001VelocityMaximumMPS,
  false,false,true,true,true,false
};

const ModelMetadata *const AllMetadata[]={
  &HMetadata,&OHMetadata,&OKVB2022DirectMetadata,&OSUMER2001DirectMetadata,
  &OSUMER2001BowenMetadata
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

bool PositiveFinite(double value) {
  return std::isfinite(value) && value>0.0;
}

// Linear interpolation on a uniform velocity table.  Evaluate() has already
// rejected velocities outside the table, so only the exact upper endpoint
// needs special handling.
double InterpolateUniform(const double *table,std::size_t length,
    double minimumMPS,double stepMPS,double velocityMPS) {
  const double coordinate=(velocityMPS-minimumMPS)/stepMPS;
  const std::size_t lower=static_cast<std::size_t>(std::floor(coordinate));
  if (lower>=length-1) return table[length-1];
  const double fraction=coordinate-static_cast<double>(lower);
  return table[lower]+fraction*(table[lower+1]-table[lower]);
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

// Inverse-square scaling and the emitted triplet power (mean photon energy
// weighted by the 3s 3S1 -> 3P2,1,0 Einstein coefficients).
Result OxygenTripletResult(double photonRate1AU,const Conditions &conditions) {
  const double photonRate=photonRate1AU*
      DistanceScale(conditions.heliocentricDistanceM);
  const double power=photonRate*OTripletMeanPhotonEnergyJ;
  if (!PositiveFinite(photonRate) || !PositiveFinite(power))
    return Failure(Status::InvalidModelResult);
  return Result{Status::Ok,photonRate,power,MissingValue()};
}

Result EvaluateOxygenKVB2022Direct(const Conditions &conditions) {
  const double flux=conditions.solarOxygen1304TripletFlux1AUPhotonM2S;
  if (!PositiveFinite(flux)) return Failure(Status::MissingSolarOxygen1304Flux);
  // The table holds 1-AU rates for the Killen et al. SOLSTICE 2009-03-01
  // irradiance; Gladstone (1992) found the line shapes and the flux shares
  // of the three components constant over the solar cycle.
  const double referenceRate=InterpolateUniform(OKVB2022DirectPhotonRate1AUS1,
      OKVB2022TableLength,OKVB2022VelocityMinimumMPS,OKVB2022VelocityStepMPS,
      conditions.radialVelocityMPS);
  return OxygenTripletResult(
      referenceRate*flux/OKVB2022ReferenceTripletFlux1AUPhotonM2S,conditions);
}

Result EvaluateOxygenSUMER2001Direct(const Conditions &conditions) {
  const double flux=conditions.solarOxygen1304TripletFlux1AUPhotonM2S;
  if (!PositiveFinite(flux)) return Failure(Status::MissingSolarOxygen1304Flux);
  const double ratePerFlux=InterpolateUniform(
      OSUMER2001DirectRatePerTripletFlux1AU,OSUMER2001TableLength,
      OSUMER2001VelocityMinimumMPS,OSUMER2001VelocityStepMPS,
      conditions.radialVelocityMPS);
  return OxygenTripletResult(ratePerFlux*flux,conditions);
}

Result EvaluateOxygenSUMER2001Bowen(const Conditions &conditions) {
  const double flux=conditions.solarLymanBetaFlux1AUPhotonM2S;
  if (!PositiveFinite(flux)) return Failure(Status::MissingSolarLymanBetaFlux);
  const double ratePerFlux=InterpolateUniform(
      OSUMER2001BowenRatePerLymanBetaFlux1AU,OSUMER2001TableLength,
      OSUMER2001VelocityMinimumMPS,OSUMER2001VelocityStepMPS,
      conditions.radialVelocityMPS);
  return OxygenTripletResult(ratePerFlux*flux,conditions);
}

bool ValidTable(const double *table,std::size_t length) {
  for (std::size_t i=0;i<length;i++)
    if (!PositiveFinite(table[i])) return false;
  return true;
}

} // namespace

bool ModelFromName(const char *name,Model *model) {
  if (name==NULL || model==NULL) return false;
  for (const ModelMetadata *metadata:AllMetadata) {
    if (std::strcmp(name,metadata->id)==0) {
      *model=metadata->model;
      return true;
    }
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
  case Status::MissingSolarOxygen1304Flux:
    return "positive finite 1-AU solar O I 130.4-nm triplet flux is required";
  case Status::MissingSolarLymanBetaFlux:
    return "positive finite 1-AU solar Lyman-beta flux is required";
  case Status::InvalidModelResult: return "model produced an invalid result";
  }
  return "invalid status";
}

const ModelMetadata *GetMetadata(Model model) {
  switch (model) {
  case Model::HydrogenLymanAlphaIKL2020: return &HMetadata;
  case Model::HydroxylAX00SchleicherAhearn1988Unquenched:
    return &OHMetadata;
  case Model::OxygenI1304DirectKVB2022GroundJ2:
    return &OKVB2022DirectMetadata;
  case Model::OxygenI1304DirectSUMER2001GroundJ2:
    return &OSUMER2001DirectMetadata;
  case Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2:
    return &OSUMER2001BowenMetadata;
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
  case Model::OxygenI1304DirectKVB2022GroundJ2:
    return EvaluateOxygenKVB2022Direct(conditions);
  case Model::OxygenI1304DirectSUMER2001GroundJ2:
    return EvaluateOxygenSUMER2001Direct(conditions);
  case Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2:
    return EvaluateOxygenSUMER2001Bowen(conditions);
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
  if (OKVB2022CompiledLength!=OKVB2022TableLength ||
      OKVB2022TableLength!=961 || OKVB2022VelocityMinimumMPS!=-48.0e3 ||
      OKVB2022VelocityMaximumMPS!=48.0e3 || OKVB2022VelocityStepMPS!=100.0 ||
      OSUMER2001DirectCompiledLength!=OSUMER2001TableLength ||
      OSUMER2001BowenCompiledLength!=OSUMER2001TableLength ||
      OSUMER2001TableLength!=2001 || OSUMER2001VelocityMinimumMPS!=-100.0e3 ||
      OSUMER2001VelocityMaximumMPS!=100.0e3 ||
      OSUMER2001VelocityStepMPS!=100.0) {
    if (error!=NULL) *error="O I 130.4-nm table geometry is invalid";
    return false;
  }
  if (!ValidTable(OKVB2022DirectPhotonRate1AUS1,OKVB2022TableLength) ||
      !ValidTable(OSUMER2001DirectRatePerTripletFlux1AU,OSUMER2001TableLength) ||
      !ValidTable(OSUMER2001BowenRatePerLymanBetaFlux1AU,OSUMER2001TableLength) ||
      !PositiveFinite(OKVB2022ReferenceTripletFlux1AUPhotonM2S) ||
      !PositiveFinite(OTripletMeanPhotonEnergyJ)) {
    if (error!=NULL) *error="O I 130.4-nm table contains an invalid value";
    return false;
  }
  if (std::strlen(HReferenceSha256)!=64 ||
      std::strlen(OHReferenceSha256)!=64 ||
      std::strlen(OKVB2022ReferenceSha256)!=64 ||
      std::strlen(OSUMER2001TripletReferenceSha256)!=64 ||
      std::strlen(OSUMER2001LymanBetaReferenceSha256)!=64 ||
      std::strlen(ONISTUltravioletReferenceSha256)!=64 ||
      std::strlen(ONISTInfraredReferenceSha256)!=64) {
    if (error!=NULL) *error="reference-data SHA-256 metadata is invalid";
    return false;
  }
  return true;
}

} // namespace ResonanceFluorescence
} // namespace AMPS
