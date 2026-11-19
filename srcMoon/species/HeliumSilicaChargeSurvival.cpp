#include "HeliumSilicaChargeSurvival.h"
#include "Photochemistry.h" // existing source-identity service, no duplicate SHA256

#include <cmath>
#include <limits>
#include <utility>

namespace Moon { namespace Species { namespace HeliumSilicaChargeSurvival {
namespace {
constexpr const char* PublicationSha="bf554c667bb01d03952d151f86e4e348b706f20630c8fe896eedf15f6fa30ed2";
constexpr const char* TableSha="bffa29b2d3a608f05210b01da52bb70fc42e64b0f4399390c40764ce967d8576";
struct Row {
  Target target;Mechanism mechanism;
  double a,da,vc,dvc,r2,ratio,dratio;
  bool r2Available,ratioAvailable;
};
// Tellez etal(2012), TableI. SI Vc and its AUTHOR-REPORTED uncertainty use
// exactly the published scale1e5m/s. Intercepts belong to Eq4 corrected ion
// yield; including exp(a) in Eq3 would not be a probability. Missing r2/ratio
// cells are represented by explicit availability flags, never synthetic data.
const Row Rows[]={
 {Target::Silicon,Mechanism::CollisionInducedNeutralization,15.78,.06,291000,6000,.99404,0,0,true,false},
 {Target::Oxygen16,Mechanism::CollisionInducedNeutralization,16.40,.16,226000,29000,.96737,1.86,.19,true,true},
 {Target::Oxygen18,Mechanism::CollisionInducedNeutralization,16.54,.15,267000,28000,.97880,2.14,.19,true,true},
 {Target::PublishedOxygenMean,Mechanism::CollisionInducedNeutralization,16.47,.16,247000,28000,0,2.00,.19,false,true},
 {Target::Oxygen16,Mechanism::AugerNeutralization,16.40,.16,188000,3000,.99998,0,0,true,false},
 {Target::Oxygen18,Mechanism::AugerNeutralization,16.54,.15,191000,2000,.99999,0,0,true,false},
 {Target::PublishedOxygenMean,Mechanism::AugerNeutralization,16.47,.16,190000,2000,0,0,0,false,false}
};
Status Refuse(Status status,const char* reason,std::string* error) {
  if(error)*error=reason;
  return status;
}
bool Known(Target target) {
  return target==Target::Silicon || target==Target::Oxygen16 ||
      target==Target::Oxygen18 || target==Target::PublishedOxygenMean;
}
bool Known(Mechanism mechanism) {
  return mechanism==Mechanism::AugerNeutralization ||
      mechanism==Mechanism::CollisionInducedNeutralization;
}
Fit Metadata(const Row& row) {
  Fit fit;fit.target=row.target;fit.mechanism=row.mechanism;
  fit.id=std::string(TargetName(row.target))+"_"+MechanismName(row.mechanism);
  fit.minimumIncidentEnergyEv=row.mechanism==Mechanism::AugerNeutralization ? 600:2000;
  fit.minimumIncluded=row.mechanism==Mechanism::AugerNeutralization;
  fit.maximumIncidentEnergyEv=row.mechanism==Mechanism::AugerNeutralization ? 800:7000;
  fit.maximumIncluded=row.mechanism==Mechanism::CollisionInducedNeutralization;
  if(row.target==Target::Silicon) {
    fit.minimumIncidentEnergyEv=600;fit.minimumIncluded=true;
  }
  fit.characteristicVelocityMS=row.vc;fit.reportedVelocityUncertaintyMS=row.dvc;
  fit.correctedYieldIntercept=row.a;fit.reportedInterceptUncertainty=row.da;
  fit.hasRSquared=row.r2Available;fit.rSquared=row.r2;
  fit.hasSurfaceOxygenSiliconRatio=row.ratioAvailable;
  fit.surfaceOxygenSiliconRatio=row.ratio;fit.reportedSurfaceRatioUncertainty=row.dratio;
  return fit;
}
bool Contains(const Fit& fit,double energy) {
  return (energy>fit.minimumIncidentEnergyEv ||
      (fit.minimumIncluded && energy==fit.minimumIncidentEnergyEv)) &&
      (energy<fit.maximumIncidentEnergyEv ||
      (fit.maximumIncluded && energy==fit.maximumIncidentEnergyEv));
}
bool Exponential(long double exponent,double* output) {
  if(!std::isfinite(exponent) || exponent>0)return false;
  const long double value=std::exp(exponent);
  const double represented=static_cast<double>(value);
  // A mathematically positive probability must not silently become zero.
  // Its true underflow is an explicit arithmetic refusal, not a new sink.
  if(!std::isfinite(value) || !std::isfinite(represented) || represented<=0 || represented>1)return false;
  *output=represented;return true;
}
}

const char* TargetName(Target target) {
  switch(target) {
    case Target::Silicon:return "Si";
    case Target::Oxygen16:return "O16";
    case Target::Oxygen18:return "O18";
    case Target::PublishedOxygenMean:return "O_mean";
  }
  return "UNSUPPORTED";
}
const char* MechanismName(Mechanism mechanism) {
  switch(mechanism) {
    case Mechanism::AugerNeutralization:return "AN";
    case Mechanism::CollisionInducedNeutralization:return "CIN";
  }
  return "UNSUPPORTED";
}
SourceEvidence FrozenV1Sources(const std::string& publication,const std::string& table) {
  SourceEvidence source;source.publicationPath=publication;source.processedTablePath=table;
  source.publicationSha256=PublicationSha;source.processedTableSha256=TableSha;
  return source;
}
bool Registry::LoadPublishedTable(const SourceEvidence& source,std::string* error) {
  if(!error)return false;
  const auto frozen=FrozenV1Sources("","");
  if(source.schema!=frozen.schema || source.publicationSha256!=PublicationSha ||
      source.processedTableSha256!=TableSha) {
    *error="Tellez2012 source schema/hash is not the frozen published laboratory table";return false;
  }
  if(!ReactionData::VerifyFileSha256(source.publicationPath,PublicationSha,error) ||
      !ReactionData::VerifyFileSha256(source.processedTablePath,TableSha,error))return false;
  // Commit only after BOTH source identities succeed; an unsuccessfully loaded
  // alternative cannot erase the previous correctly initialized registry.
  SourceEvidence candidate=source;source_=std::move(candidate);loaded_=true;
  error->clear();return true;
}
Status Registry::GetFit(Target target,Mechanism mechanism,Fit* output,std::string* error) const {
  if(!output || !Known(target) || !Known(mechanism))
    return Refuse(Status::InvalidInput,"fit query needs output and valid source target/mechanism identities",error);
  if(!loaded_)return Refuse(Status::Unavailable,"source PDF and processed TableI identities have not been loaded",error);
  for(const auto& row:Rows)if(row.target==target && row.mechanism==mechanism) {
    Fit value=Metadata(row);*output=std::move(value);if(error)error->clear();return Status::Evaluated;
  }
  return Refuse(Status::Unavailable,"TableI supplies no fit for this target/mechanism",error);
}
Status Registry::SelectFit(Target target,double energy,Fit* output,std::string* error) const {
  if(!output || !Known(target) || !std::isfinite(energy) || energy<=0)
    return Refuse(Status::InvalidInput,"source fit selection requires target and positive finite incident energy[eV]",error);
  if(!loaded_)return Refuse(Status::Unavailable,"source PDF and processed TableI identities have not been loaded",error);
  for(const auto& row:Rows)if(row.target==target) {
    Fit value=Metadata(row);
    if(Contains(value,energy)){*output=std::move(value);if(error)error->clear();return Status::Evaluated;}
  }
  return Refuse(Status::Unavailable,"incident energy lies outside TableI/study support or in an oxygen no-fit gap; no crossover/interpolation/extrapolation is defined",error);
}
Status Registry::EvaluatePublishedCoordinate(Target target,double energy,double inverse,
    Result* output,std::string* error) const {
  if(!output || !std::isfinite(inverse) || inverse<=0)
    return Refuse(Status::InvalidInput,"source reciprocal-velocity coordinate must be finite positive[s/m]",error);
  Result value;
  const auto status=SelectFit(target,energy,&value.fit,error);
  if(status!=Status::Evaluated)return status;
  const long double coordinate=inverse,vc=value.fit.characteristicVelocityMS,
      uncertainty=value.fit.reportedVelocityUncertaintyMS;
  const long double logarithm=-vc*coordinate;
  if(!std::isfinite(logarithm) || std::abs(logarithm)>std::numeric_limits<double>::max() ||
      !Exponential(logarithm,&value.chargedFraction) ||
      !Exponential(-(vc+uncertainty)*coordinate,&value.fractionAtVcPlusReportedUncertainty) ||
      !Exponential(-(vc-uncertainty)*coordinate,&value.fractionAtVcMinusReportedUncertainty))
    return Refuse(Status::NumericalError,"published exponential/parameter envelope is not representable; no zero probability or fitted clipping substituted",error);
  value.sourceReciprocalVelocitySumSM=inverse;value.logChargedFraction=static_cast<double>(logarithm);
  *output=std::move(value);if(error)error->clear();return Status::Evaluated;
}
Status Registry::EvaluateFromVelocityVectors(Target target,double energy,
    const std::array<double,3>& incoming,const std::array<double,3>& outgoing,
    const std::array<double,3>& normal,Result* output,std::string* error) const {
  if(!output)return Refuse(Status::InvalidInput,"physical velocity query needs output storage",error);
  for(int k=0;k<3;++k)if(!std::isfinite(incoming[k]) || !std::isfinite(outgoing[k]) || !std::isfinite(normal[k]))
    return Refuse(Status::InvalidInput,"physical velocity/normal components must be finite",error);
  Fit fit;const auto status=SelectFit(target,energy,&fit,error);
  if(status!=Status::Evaluated)return status;
  return Refuse(Status::Unavailable,"TableI fitted reciprocal-velocity convention has not been established from authoritative source metadata; do not guess normal/full-speed projections, elastic outgoing energy, or a lunar ion outcome",error);
}

}}}
