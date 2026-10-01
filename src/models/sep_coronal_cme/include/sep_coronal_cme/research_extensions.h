#ifndef SEP_CORONAL_CME_RESEARCH_EXTENSIONS_H
#define SEP_CORONAL_CME_RESEARCH_EXTENSIONS_H
#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_status.h"
#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME { namespace Research {

// Research products have their own schema/fingerprint. These types do not
// add a schema-5 input selector or claim a qualified production AMPS adapter.
constexpr int ResearchSchemaMajor = 6;
struct CapabilityIdentity {
  int schemaMajor = ResearchSchemaMajor;
  std::string capabilityId, algorithmVersion, coefficientFingerprint;
  std::uint64_t backgroundGeneration = 0, shockGeneration = 0;
  bool enabled = false;
};
Core::Status ValidateCapability(const CapabilityIdentity& identity);

struct ForeshockProxy {
  CapabilityIdentity identity;
  double referenceDistanceM = 0, supportWidthM = 0, minimumFactor = 1;
};
struct ForeshockEvaluation {
  double factor = 1, meanFreePathM = 0;
  std::uint64_t shockGeneration = 0;
  bool upstream = false;
};
// Signed one-sided distance is measured from the current front. The unresolved
// 0<d<L_ref layer remains owned by eta*g(p), and is rejected, not reprocessed.
Core::Result<ForeshockEvaluation> EvaluateForeshockProxy(
    const ForeshockProxy&, double signedFrontDistanceM,
    std::uint64_t currentShockGeneration, double ambientMeanFreePathM);

struct DynamicGeometryState {
  Vec3 centerM, centerVelocityMPerS, rotationAxis;
  double angleRad = 0, angularRateRadPerS = 0;
  std::array<double,3> axesM{}, axisRatesMPerS{};
  double epochS = 0;
  std::string frameId;
  // An independently reconstructed body-to-inertial attitude may provide the
  // complete matrix and inertial angular velocity instead of a fixed axis.
  // The shape tensor Q=R diag(a^-2) R^T remains a distinct derived quantity.
  bool fullAttitude = false;
  std::array<double,9> bodyToInertial{1,0,0,0,1,0,0,0,1};
  Vec3 inertialAngularVelocityPerS;
};
// The default Rodrigues route preserves the original fixed-axis reduction.
// General histories supply a proper SO(3) matrix and its independently verified
// angular velocity; Build checks orthogonality/orientation before geometry use.
class DynamicEllipsoid {
 public:
  static Core::Result<DynamicEllipsoid> Build(const DynamicGeometryState&);
  Vec3 Point(double theta, double phi) const;
  Vec3 SurfaceVelocity(double theta, double phi) const;
  Core::Result<SurfaceEvaluation> Evaluate(Vec3 positionM) const;
  Vec3 Rotate(Vec3 body) const;
  Vec3 InverseRotate(Vec3 inertial) const;
 private:
  DynamicGeometryState state_;
  Vec3 Omega() const;
};

struct ReferenceStratum {
  std::string speciesId;
  double momentumLowerSi = 0, momentumUpperSi = 0;
  std::uint64_t patchId = 0, geometryGeneration = 0;
  double timeLowerS = 0, timeUpperS = 0;
};
struct ReferenceFamilyMember {
  ReferenceStratum stratum;
  double offsetM = 0, integratedPeclet = 0, jointMeasure = 0;
  bool sourceIsSeparable = false;
};
// Positive u/kappa is integrated on the actual normal-coordinate support.
// Geometry certification receives the resulting offset and must test reach,
// folds, overlaps, masks and clearance against its immutable geometry. No
// scalar switch is allowed to bypass that independent authority.
Core::Result<ReferenceFamilyMember> SolveReferenceFamilyMember(
    const ReferenceStratum&, double targetPeclet, double maximumOffsetM,
    double quadratureTolerance, const std::function<double(double)>& inflowOverDiffusivity,
    const std::function<Core::Result<double>(double)>& geometryAndJointMeasure);
Core::Status ValidateReferenceFamily(const std::vector<ReferenceFamilyMember>&);

enum class RenewalOutcome { Absorbed, Downstream, ReReleased };
struct RenewalBranch {
  RenewalOutcome outcome = RenewalOutcome::Absorbed;
  double probability = 0, residenceTimeS = 0;
  Vec3 outgoingMomentumKgMPerS;
  double declaredShockWorkJ = 0;
  Vec3 declaredShockImpulseKgMPerS;
};
struct RenewalBalance {
  double incomingNumber = 0, absorbedNumber = 0, downstreamNumber = 0, reReleasedNumber = 0;
  double incomingEnergyJ = 0, outgoingEnergyJ = 0, shockWorkJ = 0;
  Vec3 incomingMomentum, outgoingMomentum, shockImpulse;
};
// All energies include rest energy so this is a four-momentum census. An
// absorbed branch retains energy/momentum in the downstream reservoir rather
// than deleting them. The immutable first-passage census is a separate object.
Core::Result<RenewalBalance> EvaluateRenewal(
    double representedNumber, double restMassKg, Vec3 incomingMomentum,
    const std::vector<RenewalBranch>&, bool validatedWorkAuthority,
    double relativeTolerance = 1e-10);

enum class ImpulsiveFrontRelation { NoFrontYet, Upstream, Downstream };
struct ImpulsiveBirth {
  std::uint64_t immutableBirthId = 0, shockGeneration = 0;
  double representedNumber = 0, kineticEnergyJ = 0;
  ImpulsiveFrontRelation frontRelation = ImpulsiveFrontRelation::NoFrontYet;
  std::string sourceOrigin = "ImpulsiveCoronalRelease";
};
class ImpulsiveLedger {
 public:
  Core::Status Birth(const ImpulsiveBirth&);
  Core::Status AbsorbFrontEncounter(std::uint64_t birthId, std::uint64_t currentGeneration);
  double BornNumber() const { return bornNumber_; }
  double FrontEncounterNumber() const { return frontNumber_; }
  double ShockFirstPassageNumber() const { return 0; }
  double BornKineticEnergyJ() const { return bornEnergy_; }
  double FrontEncounterKineticEnergyJ() const { return frontEnergy_; }
 private:
  std::map<std::uint64_t,ImpulsiveBirth> births_;
  std::map<std::uint64_t,bool> removed_;
  double bornNumber_ = 0, frontNumber_ = 0;
  double bornEnergy_ = 0, frontEnergy_ = 0;
};

} } }
#endif
