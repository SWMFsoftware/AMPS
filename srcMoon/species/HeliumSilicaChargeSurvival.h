#ifndef SRCMOON_SPECIES_HELIUM_SILICA_CHARGE_SURVIVAL_H
#define SRCMOON_SPECIES_HELIUM_SILICA_CHARGE_SURVIVAL_H

#include <array>
#include <string>

namespace Moon { namespace Species { namespace HeliumSilicaChargeSurvival {

// This source-specific operator describes the inferred He+ fraction among
// projectiles backscattered into the published laboratory class. It includes
// survival/reionization and is NOT an unconditional reflection, implantation,
// sticking, retention, or lunar-regolith outcome probability.
enum class Target { Silicon, Oxygen16, Oxygen18, PublishedOxygenMean };
enum class Mechanism { AugerNeutralization, CollisionInducedNeutralization };
enum class Status { Evaluated, Unavailable, InvalidInput, NumericalError };

struct SourceEvidence {
  std::string schema="srcMoon-Tellez2012-laboratory-coordinate-v1";
  std::string publicationPath,processedTablePath;
  std::string publicationSha256,processedTableSha256;
};
// The CSV is generated from the immutable publisher PDF by the committed
// helium_silica/prepare_source.py. Frozen identities cannot be updated by a
// caller to authorize other coefficients under this source version.
SourceEvidence FrozenV1Sources(const std::string& publicationPath,
                              const std::string& processedTablePath);

struct Fit {
  std::string id;
  Target target=Target::Silicon;
  Mechanism mechanism=Mechanism::CollisionInducedNeutralization;
  // TableI inequalities intersect the actual experimental600..7000eV range.
  // In particular O800eV and2000eV have no admissible fit. No epsilon,
  // interpolation or hand-selected crossover is applied at either boundary.
  double minimumIncidentEnergyEv=0,maximumIncidentEnergyEv=0;
  bool minimumIncluded=false,maximumIncluded=false;
  double characteristicVelocityMS=0;
  double reportedVelocityUncertaintyMS=0;
  double correctedYieldIntercept=0,reportedInterceptUncertainty=0;
  bool hasRSquared=false;
  double rSquared=0;
  bool hasSurfaceOxygenSiliconRatio=false;
  double surfaceOxygenSiliconRatio=0,reportedSurfaceRatioUncertainty=0;
  // These source uncertainties have no locally established confidence level
  // or covariance. The mean-O row is a separate published fit, not an
  // isotope-mixture law inferred from its two neighbouring isotope rows.
};

struct Result {
  Fit fit;
  double sourceReciprocalVelocitySumSM=0; // s/m in the EXACT source coordinate
  double logChargedFraction=0,chargedFraction=0;
  // Monotone propagation of the author's Vc+/-parameter values only. This is
  // neither a statistical confidence interval nor independent model error.
  double fractionAtVcPlusReportedUncertainty=0;
  double fractionAtVcMinusReportedUncertainty=0;
  bool physicalVelocityMappingAvailable=false;
};

class Registry {
 public:
  // Source/hash validation is transactional. Coefficients are frozen native
  // TableI transcriptions, checked against two native text layouts and the
  // original rendered table; no fit or parameter inference occurs here.
  bool LoadPublishedTable(const SourceEvidence& source,std::string* error);
  Status GetFit(Target target,Mechanism mechanism,Fit* output,
                std::string* error) const;
  Status SelectFit(Target target,double incidentEnergyEv,Fit* output,
                   std::string* error) const;
  // Eq3: P+=exp(-Vc*x), x=(1/vi+1/vf). Caller MUST supply the published,
  // mechanism-specific reciprocal-velocity coordinate, in seconds/metre.
  // Its construction from actual velocities is deliberately not guessed.
  Status EvaluatePublishedCoordinate(Target target,double incidentEnergyEv,
      double sourceReciprocalVelocitySumSM,Result* output,std::string* error) const;
  // Explicit unavailable branch: the source notes that velocity definitions
  // depend on mechanism but does not establish the projection convention of
  // each fitted TableI coefficient. Surface-rest vectors/normal may not be
  // silently converted using a normal/full-speed or elastic-scattering model.
  Status EvaluateFromVelocityVectors(Target target,double incidentEnergyEv,
      const std::array<double,3>& incomingMS,const std::array<double,3>& outgoingMS,
      const std::array<double,3>& surfaceNormal,Result* output,std::string* error) const;
  bool loaded() const { return loaded_; }
  const SourceEvidence& sources() const { return source_; }
 private:
  bool loaded_=false;
  SourceEvidence source_;
};

// Names are source identities, not PIC species indices or inferred isotopes.
const char* TargetName(Target target);
const char* MechanismName(Mechanism mechanism);

}}}
#endif
